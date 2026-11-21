/**
 * @file moment_engine.cpp
 * @brief Implementation of MomentEngine::Run (serial and state-blocked execution).
 */
#include "core/chebyshev/moment_engine.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

#include "core/algebra/linalg.hpp"
#include "core/algebra/memory_placement.hpp"
#include "core/algebra/vector_block.hpp"
#include "core/chebyshev/checkpoint.hpp"
#include "core/chebyshev/chebyshev_recursion.hpp"
#include "core/chebyshev/chebyshev_vectors.hpp"

namespace lsquant {

namespace {

using SeedOp = std::variant<ApplyOperator, ApplySeries, ApplyDerivedSeries>;

/**
 * @brief Trace-state checkpoint sink shared by the executors: after every finished state
 *        (or state block) persists the partial tensor, then optionally aborts for the
 *        fault-injection test (`LSQ_CHECKPOINT_ABORT_AFTER`). Inactive when `path` is empty.
 */
struct CheckpointHook {
  std::string path;
  Index statesDone = 0;
  Index abortAfter = -1;

  void operator()(const ChebyshevMoments& mu, Index finished) {
    if (path.empty()) {
      return;
    }
    statesDone += finished;
    WriteCheckpoint(path, mu, statesDone);
    if (abortAfter >= 0 && statesDone >= abortAfter) {
      throw std::runtime_error("LSQ_CHECKPOINT_ABORT_AFTER: simulated kill after " +
                               std::to_string(statesDone) + " states");
    }
  }
};

struct SideProgram {
  std::vector<SeedOp> preOps;
  bool hasExpand = false;
  int numMoments = 0;
  std::vector<int> postOperators;
};

// Everything `Run` needs from `plan`, fully validated.
struct ParsedPlan {
  bool hasTimeAxis = false;
  TimeAxis timeAxisSpec;
  std::vector<SeedOp> seedKetOps;
  std::vector<SeedOp> seedBraOps;
  bool seedMirror = false;
  SideProgram ketProg;
  SideProgram braProg;
  bool seenMirror = false;
  std::vector<Scalar> traceWeights;
  std::vector<int> dims;
  bool scalarWeight = true;
  bool validMomentWeights = false;
};

ParsedPlan ParsePlan(const Plan& plan, int numOperators) {
  if (plan.code.empty()) {
    throw std::invalid_argument("MomentEngine::Run: plan code must not be empty");
  }

  std::size_t timeAxisPos = plan.code.size();
  for (std::size_t pc = 0; pc < plan.code.size(); ++pc) {
    if (std::holds_alternative<TimeAxis>(plan.code[pc])) {
      if (timeAxisPos != plan.code.size()) {
        throw std::invalid_argument("MomentEngine::Run: TimeAxis may appear at most once");
      }
      timeAxisPos = pc;
    }
  }

  ParsedPlan pp;
  pp.hasTimeAxis = (timeAxisPos != plan.code.size());

  if (pp.hasTimeAxis) {
    const auto* ta = std::get_if<TimeAxis>(&plan.code[timeAxisPos]);
    if (ta->numSteps <= 0) {
      throw std::invalid_argument("MomentEngine::Run: TimeAxis numSteps must be positive");
    }
    for (const AxisStep& ev : ta->evolve) {
      if (const auto* ed = std::get_if<EvolveDerived>(&ev)) {
        if (ed->operatorIndex < 0 || ed->operatorIndex >= numOperators) {
          throw std::invalid_argument(
              "MomentEngine::Run: EvolveDerived operatorIndex out of range");
        }
      }
    }
    pp.timeAxisSpec = *ta;

    std::vector<SeedOp>* activeSeedOps = &pp.seedKetOps;
    for (std::size_t pc = 0; pc < timeAxisPos; ++pc) {
      const Instruction& inst = plan.code[pc];
      if (const auto* app = std::get_if<ApplyOperator>(&inst)) {
        if (app->operatorIndex < 0 || app->operatorIndex >= numOperators) {
          throw std::invalid_argument("MomentEngine::Run: operatorIndex out of range");
        }
        activeSeedOps->push_back(*app);
      } else if (const auto* ser = std::get_if<ApplySeries>(&inst)) {
        activeSeedOps->push_back(*ser);
      } else if (const auto* dser = std::get_if<ApplyDerivedSeries>(&inst)) {
        if (dser->operatorIndex < 0 || dser->operatorIndex >= numOperators) {
          throw std::invalid_argument("MomentEngine::Run: operatorIndex out of range");
        }
        activeSeedOps->push_back(*dser);
      } else if (std::holds_alternative<MirrorToBra>(inst)) {
        if (pp.seedMirror) {
          throw std::invalid_argument(
              "MomentEngine::Run: MirrorToBra may appear at most once before TimeAxis");
        }
        pp.seedMirror = true;
        activeSeedOps = &pp.seedBraOps;
      } else {
        throw std::invalid_argument(
            "MomentEngine::Run: Expand and Trace cannot precede TimeAxis");
      }
    }
  }

  SideProgram* activeSide = &pp.ketProg;
  bool seenTrace = false;

  const std::size_t bodyStart = pp.hasTimeAxis ? (timeAxisPos + 1) : 0;
  for (std::size_t pc = bodyStart; pc < plan.code.size(); ++pc) {
    if (seenTrace) {
      throw std::invalid_argument("MomentEngine::Run: Trace must be the final instruction");
    }
    const Instruction& inst = plan.code[pc];
    if (const auto* app = std::get_if<ApplyOperator>(&inst)) {
      if (app->operatorIndex < 0 || app->operatorIndex >= numOperators) {
        throw std::invalid_argument("MomentEngine::Run: operatorIndex out of range");
      }
      if (!activeSide->hasExpand) {
        activeSide->preOps.push_back(*app);
      } else {
        activeSide->postOperators.push_back(app->operatorIndex);
      }
    } else if (const auto* ser = std::get_if<ApplySeries>(&inst)) {
      if (activeSide->hasExpand) {
        throw std::invalid_argument("MomentEngine::Run: ApplySeries must precede Expand");
      }
      activeSide->preOps.push_back(*ser);
    } else if (const auto* dser = std::get_if<ApplyDerivedSeries>(&inst)) {
      if (dser->operatorIndex < 0 || dser->operatorIndex >= numOperators) {
        throw std::invalid_argument("MomentEngine::Run: operatorIndex out of range");
      }
      if (activeSide->hasExpand) {
        throw std::invalid_argument("MomentEngine::Run: ApplyDerivedSeries must precede Expand");
      }
      activeSide->preOps.push_back(*dser);
    } else if (const auto* exp = std::get_if<Expand>(&inst)) {
      if (exp->numMoments <= 0) {
        throw std::invalid_argument("MomentEngine::Run: Expand numMoments must be positive");
      }
      if (activeSide->hasExpand) {
        throw std::runtime_error("Expand: rank > 2 not implemented");
      }
      activeSide->hasExpand = true;
      activeSide->numMoments = exp->numMoments;
    } else if (std::holds_alternative<MirrorToBra>(inst)) {
      if (pp.seenMirror) {
        throw std::invalid_argument("MomentEngine::Run: MirrorToBra may appear at most once");
      }
      pp.seenMirror = true;
      activeSide = &pp.braProg;
    } else if (const auto* tr = std::get_if<Trace>(&inst)) {
      seenTrace = true;
      pp.traceWeights = tr->weights;
    }
  }

  if (!seenTrace) {
    throw std::invalid_argument("MomentEngine::Run: plan must end with a Trace instruction");
  }

  if (pp.braProg.hasExpand) {
    pp.dims.push_back(pp.braProg.numMoments);
  }
  if (pp.ketProg.hasExpand) {
    pp.dims.push_back(pp.ketProg.numMoments);
  }
  if (pp.hasTimeAxis) {
    if (pp.braProg.hasExpand) {
      throw std::invalid_argument(
          "MomentEngine::Run: TimeAxis body does not support bra Expand");
    }
    pp.dims.push_back(pp.timeAxisSpec.numSteps);
  }

  std::size_t tensorSize = 1;
  for (int d : pp.dims) {
    tensorSize *= static_cast<std::size_t>(d);
  }
  pp.validMomentWeights =
      pp.hasTimeAxis && pp.ketProg.hasExpand &&
      (pp.traceWeights.size() == static_cast<std::size_t>(pp.ketProg.numMoments));
  if (pp.traceWeights.size() != 1 && pp.traceWeights.size() != tensorSize &&
      !pp.validMomentWeights) {
    throw std::invalid_argument(
        "MomentEngine::Run: Trace weights size must be 1 or match tensor size");
  }
  pp.scalarWeight = (pp.traceWeights.size() == 1);
  return pp;
}

void ApplySingleVectorOperators(const std::vector<const LinearOperator*>& operators,
                                const std::vector<int>& opIndices, const Vector& init,
                                Vector& out, Vector& scratch) {
  linalg::copy(init, out);
  for (int idx : opIndices) {
    const LinearOperator& op = *operators[static_cast<std::size_t>(idx)];
    if (op.IsIdentity()) {
      continue;
    }
    linalg::copy(out, scratch);
    op.Multiply(Scalar(1.0, 0.0), scratch, Scalar(0.0, 0.0), out);
  }
}

void RunApplySeries(const LinearOperator& rescaledH, const std::vector<Scalar>& coeffs,
                    Scalar scale, Vector& phi, Vector& cheb0, Vector& cheb1) {
  if (coeffs.empty()) {
    linalg::scal(Scalar(0.0, 0.0), phi);
    return;
  }
  linalg::copy(phi, cheb0);
  rescaledH.Multiply(Scalar(1.0, 0.0), cheb0, Scalar(0.0, 0.0), cheb1);
  linalg::scal(Scalar(0.0, 0.0), phi);
  for (std::size_t n = 0; n < coeffs.size(); ++n) {
    linalg::axpy(coeffs[n], cheb0, phi);
    rescaledH.Multiply(Scalar(2.0, 0.0), cheb1, Scalar(-1.0, 0.0), cheb0);
    cheb0.swap(cheb1);
  }
  linalg::scal(scale, phi);
}

void RunApplyDerivedSeries(const LinearOperator& rescaledH, const LinearOperator& commOp,
                           const std::vector<Scalar>& coeffs, Scalar scale, Vector& phi,
                           Vector& cheb0, Vector& cheb1, Vector& chebx0, Vector& chebx1) {
  if (coeffs.empty()) {
    linalg::scal(Scalar(0.0, 0.0), phi);
    return;
  }
  linalg::copy(phi, cheb0);
  linalg::copy(phi, chebx0);
  linalg::scal(Scalar(0.0, 0.0), chebx0);
  rescaledH.Multiply(Scalar(1.0, 0.0), cheb0, Scalar(0.0, 0.0), cheb1);
  commOp.Multiply(Scalar(1.0, 0.0), cheb0, Scalar(0.0, 0.0), chebx1);
  linalg::scal(Scalar(0.0, 0.0), phi);
  for (std::size_t n = 0; n < coeffs.size(); ++n) {
    const Scalar cn = (n == 0) ? Scalar(0.0, 0.0) : coeffs[n];
    linalg::axpy(cn, chebx0, phi);
    rescaledH.Multiply(Scalar(2.0, 0.0), cheb1, Scalar(-1.0, 0.0), cheb0);
    rescaledH.Multiply(Scalar(2.0, 0.0), chebx1, Scalar(-1.0, 0.0), chebx0);
    commOp.Multiply(Scalar(2.0, 0.0), cheb1, Scalar(1.0, 0.0), chebx0);
    cheb0.swap(cheb1);
    chebx0.swap(chebx1);
  }
  linalg::scal(scale, phi);
}

void ApplySeedOps(const LinearOperator& rescaledH,
                  const std::vector<const LinearOperator*>& operators,
                  const std::vector<SeedOp>& ops, const Vector& init, Vector& out,
                  Vector& scratch, Vector& cheb0, Vector& cheb1, Vector& chebx0,
                  Vector& chebx1) {
  linalg::copy(init, out);
  for (const SeedOp& op : ops) {
    if (const auto* app = std::get_if<ApplyOperator>(&op)) {
      const LinearOperator& linOp = *operators[static_cast<std::size_t>(app->operatorIndex)];
      if (linOp.IsIdentity()) {
        continue;
      }
      linalg::copy(out, scratch);
      linOp.Multiply(Scalar(1.0, 0.0), scratch, Scalar(0.0, 0.0), out);
    } else if (const auto* ser = std::get_if<ApplySeries>(&op)) {
      RunApplySeries(rescaledH, ser->coeffs, ser->scale, out, cheb0, cheb1);
    } else if (const auto* dser = std::get_if<ApplyDerivedSeries>(&op)) {
      const LinearOperator& commOp = *operators[static_cast<std::size_t>(dser->operatorIndex)];
      RunApplyDerivedSeries(rescaledH, commOp, dser->coeffs, dser->scale, out, cheb0, cheb1,
                            chebx0, chebx1);
    }
  }
}

void AccumulateWeightedDot(Scalar& target, const Scalar& weight, const Vector& bra,
                           const Vector& ket) {
  const Scalar d = linalg::dot(bra, ket);
  if (weight.imag() == 0.0) {
    target += weight.real() * d;
  } else {
    target += weight * d;
  }
}

void AccumulateWeightedDotPtr(Scalar& target, const Scalar& weight, const Scalar* bra,
                              const Scalar* ket, Index n) {
  const Scalar d = linalg::dot(bra, ket, n);
  if (weight.imag() == 0.0) {
    target += weight.real() * d;
  } else {
    target += weight * d;
  }
}

int PaddedMoments(int numMoments, int batchSize) {
  return ((numMoments + batchSize - 1) / batchSize) * batchSize;
}

Vector& SelectRegister(Register reg, Vector& ket, Vector& bra, Vector& aux) {
  switch (reg) {
    case Register::Ket:
      return ket;
    case Register::Bra:
      return bra;
    case Register::Aux:
      return aux;
  }
  return ket;
}

// ---------------------------------------------------------------------------
// Block (state_block > 1) helpers. Every block operation applies the identical
// per-column floating-point sequence of the corresponding serial routine, so a
// state blocked run reproduces the state_block = 1 result bit for bit.
// ---------------------------------------------------------------------------

VectorBlock MakeBlock(Index rows, int cols) {
  return VectorBlock(rows, cols);
}

void CopyBlock(const VectorBlock& src, VectorBlock& dst) {
  std::copy_n(src.Data(), static_cast<std::size_t>(src.Rows()) * src.Cols(), dst.Data());
}

void ScalBlock(Scalar a, VectorBlock& x) { x.Map() *= a; }

void AxpyBlock(Scalar a, const VectorBlock& x, VectorBlock& y) {
  y.Map() = a * x.Map() + y.Map();
}

// Three-term Chebyshev recursion on an `N x S` block: per column identical to
// `ChebyshevRecursion`, with one `MultiplyBlock` per polynomial order.
struct BlockRecursion {
  BlockRecursion(const LinearOperator& h, int cols)
      : H(&h), prev(MakeBlock(h.Dimension(), cols)), cur(MakeBlock(h.Dimension(), cols)),
        next(MakeBlock(h.Dimension(), cols)) {}

  void Start(const VectorBlock& phi) {
    CopyBlock(phi, cur);
    CopyBlock(phi, prev);
    H->MultiplyBlock(Scalar(1.0, 0.0), cur, Scalar(0.0, 0.0), next);
    order = 0;
  }

  void Step() {
    if (order == 0) {
      CopyBlock(cur, prev);
      cur.Swap(next);
      order = 1;
      return;
    }
    H->MultiplyBlock(Scalar(2.0, 0.0), cur, Scalar(-1.0, 0.0), prev);
    prev.Swap(cur);
    ++order;
  }

  const LinearOperator* H;
  VectorBlock prev;
  VectorBlock cur;
  VectorBlock next;
  int order = -1;
};

// Batched window of Chebyshev vectors for `states` states advancing together.
// Column layout is state-major: column `s * batch + m0` holds `T_{firstOrder+m0}` of state
// `s`, so the columns of one state form a contiguous sub-block.
struct BlockChebWindow {
  BlockChebWindow(Index dim, int batchSize, int numStates)
      : n(dim), batch(batchSize), states(numStates),
        vectors(MakeBlock(dim, batchSize * numStates)),
        scratch(MakeBlock(dim, batchSize * numStates)) {}

  void Fill(BlockRecursion& rec) {
    firstOrder = rec.order;
    for (int m = 0; m < batch; ++m) {
      for (int s = 0; s < states; ++s) {
        std::copy_n(rec.cur.Column(s), static_cast<std::size_t>(n),
                    vectors.Column(s * batch + m));
      }
      rec.Step();
    }
  }

  void ApplyToAll(const LinearOperator& A) {
    if (A.IsIdentity()) {
      return;
    }
    A.MultiplyBlock(Scalar(1.0, 0.0), vectors, Scalar(0.0, 0.0), scratch);
    vectors.Swap(scratch);
  }

  const Scalar* Column(int s, int m0) const { return vectors.Column(s * batch + m0); }
  const Scalar* StateBlockData(int s) const { return vectors.Data() + columnOffset(s); }
  std::size_t columnOffset(int s) const {
    return static_cast<std::size_t>(s) * batch * static_cast<std::size_t>(n);
  }

  Index n;
  int batch;
  int states;
  int firstOrder = 0;
  VectorBlock vectors;
  VectorBlock scratch;
};

void RunApplySeriesBlock(const LinearOperator& rescaledH, const std::vector<Scalar>& coeffs,
                         Scalar scale, VectorBlock& phi, VectorBlock& cheb0,
                         VectorBlock& cheb1) {
  if (coeffs.empty()) {
    ScalBlock(Scalar(0.0, 0.0), phi);
    return;
  }
  CopyBlock(phi, cheb0);
  rescaledH.MultiplyBlock(Scalar(1.0, 0.0), cheb0, Scalar(0.0, 0.0), cheb1);
  ScalBlock(Scalar(0.0, 0.0), phi);
  for (std::size_t n = 0; n < coeffs.size(); ++n) {
    AxpyBlock(coeffs[n], cheb0, phi);
    rescaledH.MultiplyBlock(Scalar(2.0, 0.0), cheb1, Scalar(-1.0, 0.0), cheb0);
    cheb0.Swap(cheb1);
  }
  ScalBlock(scale, phi);
}

void RunApplyDerivedSeriesBlock(const LinearOperator& rescaledH, const LinearOperator& commOp,
                                const std::vector<Scalar>& coeffs, Scalar scale,
                                VectorBlock& phi, VectorBlock& cheb0, VectorBlock& cheb1,
                                VectorBlock& chebx0, VectorBlock& chebx1) {
  if (coeffs.empty()) {
    ScalBlock(Scalar(0.0, 0.0), phi);
    return;
  }
  CopyBlock(phi, cheb0);
  CopyBlock(phi, chebx0);
  ScalBlock(Scalar(0.0, 0.0), chebx0);
  rescaledH.MultiplyBlock(Scalar(1.0, 0.0), cheb0, Scalar(0.0, 0.0), cheb1);
  commOp.MultiplyBlock(Scalar(1.0, 0.0), cheb0, Scalar(0.0, 0.0), chebx1);
  ScalBlock(Scalar(0.0, 0.0), phi);
  for (std::size_t n = 0; n < coeffs.size(); ++n) {
    const Scalar cn = (n == 0) ? Scalar(0.0, 0.0) : coeffs[n];
    AxpyBlock(cn, chebx0, phi);
    rescaledH.MultiplyBlock(Scalar(2.0, 0.0), cheb1, Scalar(-1.0, 0.0), cheb0);
    rescaledH.MultiplyBlock(Scalar(2.0, 0.0), chebx1, Scalar(-1.0, 0.0), chebx0);
    commOp.MultiplyBlock(Scalar(2.0, 0.0), cheb1, Scalar(1.0, 0.0), chebx0);
    cheb0.Swap(cheb1);
    chebx0.Swap(chebx1);
  }
  ScalBlock(scale, phi);
}

void ApplySeedOpsBlock(const LinearOperator& rescaledH,
                       const std::vector<const LinearOperator*>& operators,
                       const std::vector<SeedOp>& ops, const VectorBlock& init,
                       VectorBlock& out, VectorBlock& scratch, VectorBlock& cheb0,
                       VectorBlock& cheb1, VectorBlock& chebx0, VectorBlock& chebx1) {
  CopyBlock(init, out);
  for (const SeedOp& op : ops) {
    if (const auto* app = std::get_if<ApplyOperator>(&op)) {
      const LinearOperator& linOp = *operators[static_cast<std::size_t>(app->operatorIndex)];
      if (linOp.IsIdentity()) {
        continue;
      }
      CopyBlock(out, scratch);
      linOp.MultiplyBlock(Scalar(1.0, 0.0), scratch, Scalar(0.0, 0.0), out);
    } else if (const auto* ser = std::get_if<ApplySeries>(&op)) {
      RunApplySeriesBlock(rescaledH, ser->coeffs, ser->scale, out, cheb0, cheb1);
    } else if (const auto* dser = std::get_if<ApplyDerivedSeries>(&op)) {
      const LinearOperator& commOp = *operators[static_cast<std::size_t>(dser->operatorIndex)];
      RunApplyDerivedSeriesBlock(rescaledH, commOp, dser->coeffs, dser->scale, out, cheb0,
                                 cheb1, chebx0, chebx1);
    }
  }
}

void ApplySingleVectorOperatorsBlock(const std::vector<const LinearOperator*>& operators,
                                     const std::vector<int>& opIndices, const VectorBlock& init,
                                     VectorBlock& out, VectorBlock& scratch) {
  CopyBlock(init, out);
  for (int idx : opIndices) {
    const LinearOperator& op = *operators[static_cast<std::size_t>(idx)];
    if (op.IsIdentity()) {
      continue;
    }
    CopyBlock(out, scratch);
    op.MultiplyBlock(Scalar(1.0, 0.0), scratch, Scalar(0.0, 0.0), out);
  }
}

#ifdef LSQUANT_NATIVE_BUILD
// Apply the post-operator chain to `src` without an initial copy: the first non-identity
// operator reads `src` directly, later ones alternate between `bufA` and `bufB`. Returns a
// pointer to the block holding the result (`&src` when every operator is the identity).
const VectorBlock* ApplyPostOpsChain(const std::vector<const LinearOperator*>& operators,
                                     const std::vector<int>& opIndices, const VectorBlock& src,
                                     VectorBlock& bufA, VectorBlock& bufB) {
  const VectorBlock* in = &src;
  VectorBlock* out = &bufA;
  for (int idx : opIndices) {
    const LinearOperator& op = *operators[static_cast<std::size_t>(idx)];
    if (op.IsIdentity()) {
      continue;
    }
    op.MultiplyBlock(Scalar(1.0, 0.0), *in, Scalar(0.0, 0.0), *out);
    in = out;
    out = (out == &bufA) ? &bufB : &bufA;
  }
  return in;
}

// Copy column `s` of `src` into `dstColumnBase + s * dstColumnStride` for s < count, over a
// flat (state, row-chunk) task grid so all threads work even when `count` is small. Each
// element is written by exactly one task; copies carry no floating-point arithmetic, so the
// grid shape is free to target ~1024 tasks.
void CopyColumnsGrid(const VectorBlock& src, int count, Index n, Scalar* dstColumnBase,
                     Index dstColumnStride) {
  const Index rowChunks = std::max<Index>(1, 1024 / count);
#pragma omp parallel for schedule(static)
  for (Index t = 0; t < static_cast<Index>(count) * rowChunks; ++t) {
    const int s = static_cast<int>(t / rowChunks);
    const Index rc = t % rowChunks;
    const Index r0 = rc * n / rowChunks;
    const Index r1 = (rc + 1) * n / rowChunks;
    Scalar* dst = dstColumnBase + static_cast<std::size_t>(s) *
                                      static_cast<std::size_t>(dstColumnStride) *
                                      static_cast<std::size_t>(n) +
                  r0;
    std::copy_n(src.Column(s) + r0, static_cast<std::size_t>(r1 - r0), dst);
  }
}
#endif

VectorBlock& SelectRegisterBlock(Register reg, VectorBlock& ket, VectorBlock& bra,
                                 VectorBlock& aux) {
  switch (reg) {
    case Register::Ket:
      return ket;
    case Register::Bra:
      return bra;
    case Register::Aux:
      return aux;
  }
  return ket;
}

// Serial execution (state_block == 1): the legacy state-at-a-time path.
void RunSerialExecute(const LinearOperator& rescaledH,
                      const std::vector<const LinearOperator*>& operators,
                      StateGenerator& states, int batchSize, const ParsedPlan& pp,
                      ChebyshevMoments& mu, CheckpointHook& hook) {
  const bool hasTimeAxis = pp.hasTimeAxis;
  const TimeAxis& timeAxisSpec = pp.timeAxisSpec;
  const SideProgram& ketProg = pp.ketProg;
  const SideProgram& braProg = pp.braProg;
  const bool seenMirror = pp.seenMirror;
  const std::vector<Scalar>& traceWeights = pp.traceWeights;
  const bool scalarWeight = pp.scalarWeight;
  const bool validMomentWeights = pp.validMomentWeights;

  const Index dim = rescaledH.Dimension();

#ifdef LSQUANT_NATIVE_BUILD
  Vector r = MakeZeroVector(dim);
  Vector ketSeed = MakeZeroVector(dim);
  Vector braSeed = MakeZeroVector(dim);
  Vector auxSeed = MakeZeroVector(dim);
  Vector ketInit = MakeZeroVector(dim);
  Vector braInit = MakeZeroVector(dim);
  Vector scratch = MakeZeroVector(dim);
  Vector cheb0 = MakeZeroVector(dim);
  Vector cheb1 = MakeZeroVector(dim);
  Vector chebx0 = MakeZeroVector(dim);
  Vector chebx1 = MakeZeroVector(dim);
  Vector phiT = MakeZeroVector(dim);
  Vector postVec = MakeZeroVector(dim);
#else
  Vector r(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector ketSeed(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector braSeed(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector auxSeed(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector ketInit(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector braInit(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector scratch(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector cheb0(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector cheb1(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector chebx0(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector chebx1(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector phiT(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  Vector postVec(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
#endif

  if (hasTimeAxis) {
    ChebyshevRecursion rec(rescaledH);
    while (states.Next(r)) {
      linalg::copy(r, auxSeed);
      ApplySeedOps(rescaledH, operators, pp.seedKetOps, r, ketSeed, scratch, cheb0, cheb1,
                   chebx0, chebx1);
      const Vector& braSeedStart = pp.seedMirror ? ketSeed : r;
      ApplySeedOps(rescaledH, operators, pp.seedBraOps, braSeedStart, braSeed, scratch, cheb0,
                   cheb1, chebx0, chebx1);

      for (int stepIdx = 0; stepIdx < timeAxisSpec.numSteps; ++stepIdx) {
        ApplySeedOps(rescaledH, operators, ketProg.preOps, ketSeed, ketInit, scratch, cheb0,
                     cheb1, chebx0, chebx1);
        const Vector& braBodyStart = seenMirror ? ketInit : braSeed;
        ApplySeedOps(rescaledH, operators, braProg.preOps, braBodyStart, braInit, scratch,
                     cheb0, cheb1, chebx0, chebx1);

        if (ketProg.hasExpand) {
          const int numMoms = ketProg.numMoments;
          rec.Start(ketInit);
          for (int m = 0; m < numMoms; ++m) {
            const Vector* ketVec = &rec.Current();
            if (!ketProg.postOperators.empty()) {
              ApplySingleVectorOperators(operators, ketProg.postOperators, rec.Current(),
                                         postVec, scratch);
              ketVec = &postVec;
            }
            std::size_t wIdx = 0;
            if (validMomentWeights) {
              wIdx = static_cast<std::size_t>(m);
            } else if (!scalarWeight) {
              wIdx = static_cast<std::size_t>(m) *
                         static_cast<std::size_t>(timeAxisSpec.numSteps) +
                     static_cast<std::size_t>(stepIdx);
            }
            AccumulateWeightedDot(mu.At(m, stepIdx), traceWeights[wIdx], braInit, *ketVec);
            rec.Step();
          }
        } else {
          const std::size_t wIdx = scalarWeight ? 0 : static_cast<std::size_t>(stepIdx);
          AccumulateWeightedDot(mu.At(stepIdx), traceWeights[wIdx], braInit, ketInit);
        }

        for (const AxisStep& ev : timeAxisSpec.evolve) {
          if (const auto* er = std::get_if<EvolveRegister>(&ev)) {
            Vector& targetVec = SelectRegister(er->target, ketSeed, braSeed, auxSeed);
            RunApplySeries(rescaledH, er->coeffs, er->scale, targetVec, cheb0, cheb1);
          } else if (const auto* ed = std::get_if<EvolveDerived>(&ev)) {
            const Vector& srcVec = SelectRegister(ed->source, ketSeed, braSeed, auxSeed);
            Vector& targetVec = SelectRegister(ed->target, ketSeed, braSeed, auxSeed);
            linalg::copy(srcVec, phiT);
            const LinearOperator& commOp =
                *operators[static_cast<std::size_t>(ed->operatorIndex)];
            RunApplyDerivedSeries(rescaledH, commOp, ed->coeffs, ed->scale, phiT, cheb0,
                                  cheb1, chebx0, chebx1);
            linalg::axpy(Scalar(1.0, 0.0), phiT, targetVec);
          } else if (const auto* cp = std::get_if<CopyRegister>(&ev)) {
            const Vector& fromVec = SelectRegister(cp->from, ketSeed, braSeed, auxSeed);
            Vector& toVec = SelectRegister(cp->to, ketSeed, braSeed, auxSeed);
            linalg::copy(fromVec, toVec);
          }
        }
      }
      hook(mu, 1);
    }
    return;
  }

  if (!braProg.hasExpand && !ketProg.hasExpand) {
    // Rank 0: neither side opens a Chebyshev index.
    while (states.Next(r)) {
      ApplySeedOps(rescaledH, operators, ketProg.preOps, r, ketInit, scratch, cheb0, cheb1,
                   chebx0, chebx1);
      const Vector& braStart = r;
      ApplySeedOps(rescaledH, operators, braProg.preOps, braStart, braInit, scratch, cheb0,
                   cheb1, chebx0, chebx1);
      AccumulateWeightedDot(mu.Data()[0], traceWeights[0], braInit, ketInit);
      hook(mu, 1);
    }
    return;
  }

  if (!braProg.hasExpand && ketProg.hasExpand) {
    // Rank 1 (ket expanded): matches legacy SpectralMoments when postOperators is empty.
    const int numMoms = ketProg.numMoments;
    if (ketProg.postOperators.empty()) {
      ChebyshevRecursion rec(rescaledH);
      while (states.Next(r)) {
        ApplySeedOps(rescaledH, operators, ketProg.preOps, r, ketInit, scratch, cheb0, cheb1,
                     chebx0, chebx1);
        const Vector& braStart = r;
        ApplySeedOps(rescaledH, operators, braProg.preOps, braStart, braInit, scratch, cheb0,
                     cheb1, chebx0, chebx1);
        rec.Start(ketInit);
        for (int m = 0; m < numMoms; ++m) {
          const Scalar& w = scalarWeight ? traceWeights[0]
                                         : traceWeights[static_cast<std::size_t>(m)];
          AccumulateWeightedDot(mu.At(m), w, braInit, rec.Current());
          rec.Step();
        }
        hook(mu, 1);
      }
    } else {
      const int effBatch = std::min(batchSize, numMoms);
      const int paddedM = PaddedMoments(numMoms, effBatch);
      ChebyshevRecursion rec(rescaledH);
      ChebyshevVectors vecR(dim, effBatch);
      while (states.Next(r)) {
        ApplySeedOps(rescaledH, operators, ketProg.preOps, r, ketInit, scratch, cheb0, cheb1,
                     chebx0, chebx1);
        const Vector& braStart = r;
        ApplySeedOps(rescaledH, operators, braProg.preOps, braStart, braInit, scratch, cheb0,
                     cheb1, chebx0, chebx1);
        rec.Start(ketInit);
        for (int mR = 0; mR < paddedM; mR += effBatch) {
          vecR.Fill(rec);
          for (int opIdx : ketProg.postOperators) {
            vecR.ApplyToAll(*operators[static_cast<std::size_t>(opIdx)]);
          }
          for (int m0 = 0; m0 < effBatch; ++m0) {
            const int m = mR + m0;
            if (m < numMoms) {
              const Scalar& w = scalarWeight ? traceWeights[0]
                                             : traceWeights[static_cast<std::size_t>(m)];
              AccumulateWeightedDotPtr(mu.At(m), w, braInit.data(), vecR.Column(m0), dim);
            }
          }
        }
        hook(mu, 1);
      }
    }
    return;
  }

  if (braProg.hasExpand && !ketProg.hasExpand) {
    // Rank 1 (bra expanded).
    const int numMoms = braProg.numMoments;
    const int effBatch = std::min(batchSize, numMoms);
    const int paddedM = PaddedMoments(numMoms, effBatch);
    ChebyshevRecursion recL(rescaledH);
    ChebyshevVectors vecL(dim, effBatch);
    while (states.Next(r)) {
      ApplySeedOps(rescaledH, operators, ketProg.preOps, r, ketInit, scratch, cheb0, cheb1,
                   chebx0, chebx1);
      const Vector& braStart = r;
      ApplySeedOps(rescaledH, operators, braProg.preOps, braStart, braInit, scratch, cheb0,
                   cheb1, chebx0, chebx1);
      recL.Start(braInit);
      for (int mL = 0; mL < paddedM; mL += effBatch) {
        vecL.Fill(recL);
        for (int opIdx : braProg.postOperators) {
          vecL.ApplyToAll(*operators[static_cast<std::size_t>(opIdx)]);
        }
        for (int m0 = 0; m0 < effBatch; ++m0) {
          const int m = mL + m0;
          if (m < numMoms) {
            const Scalar& w = scalarWeight ? traceWeights[0]
                                           : traceWeights[static_cast<std::size_t>(m)];
            AccumulateWeightedDotPtr(mu.At(m), w, vecL.Column(m0), ketInit.data(), dim);
          }
        }
      }
      hook(mu, 1);
    }
    return;
  }

  // Rank 2 (both bra and ket expanded): matches legacy CorrelationExpansionMoments block order.
  const int numMomsL = braProg.numMoments;
  const int numMomsR = ketProg.numMoments;
  const int effBatchL = std::min(batchSize, numMomsL);
  const int effBatchR = std::min(batchSize, numMomsR);
  const int paddedML = PaddedMoments(numMomsL, effBatchL);
  const int paddedMR = PaddedMoments(numMomsR, effBatchR);

  ChebyshevRecursion recL(rescaledH);
  ChebyshevRecursion recR(rescaledH);
  ChebyshevVectors vecL(dim, effBatchL);
  ChebyshevVectors vecR(dim, effBatchR);

  while (states.Next(r)) {
    ApplySeedOps(rescaledH, operators, ketProg.preOps, r, ketSeed, scratch, cheb0, cheb1,
                 chebx0, chebx1);
    // Static words build both sides independently from the SAME trace state.
    // Any shared preparation must be explicit on each side of MirrorToBra.
    const Vector& braStart = r;
    ApplySeedOps(rescaledH, operators, braProg.preOps, braStart, braInit, scratch, cheb0,
                 cheb1, chebx0, chebx1);
    recL.Start(braInit);
    for (int mL = 0; mL < paddedML; mL += effBatchL) {
      vecL.Fill(recL);
      for (int opIdx : braProg.postOperators) {
        vecL.ApplyToAll(*operators[static_cast<std::size_t>(opIdx)]);
      }

      linalg::copy(ketSeed, ketInit);
      recR.Start(ketInit);
      for (int mR = 0; mR < paddedMR; mR += effBatchR) {
        vecR.Fill(recR);
        for (int opIdx : ketProg.postOperators) {
          vecR.ApplyToAll(*operators[static_cast<std::size_t>(opIdx)]);
        }

#ifdef LSQUANT_NATIVE_BUILD
        const int liveL = std::min(effBatchL, numMomsL - mL);
        const int liveR = std::min(effBatchR, numMomsR - mR);
        const auto table = linalg::MomentTableBlock(vecL.Block().Data(), liveL,
                                                    vecR.Block().Data(), liveR, dim);
        for (int m = 0; m < liveL; ++m) {
          for (int n = 0; n < liveR; ++n) {
            const std::size_t index = static_cast<std::size_t>(mL + m) * numMomsR + mR + n;
            const Scalar& weight = scalarWeight ? traceWeights[0] : traceWeights[index];
            const Scalar& value = table[static_cast<std::size_t>(m) * liveR + n];
            if (weight.imag() == 0.0) mu.At(mL + m, mR + n) += weight.real() * value;
            else mu.At(mL + m, mR + n) += weight * value;
          }
        }
#else
        // Each row has one owner; dot products and trace-state additions keep
        // their serial order, so thread count cannot change the moment values.
#pragma omp parallel for schedule(static)
        for (int m0 = 0; m0 < effBatchL; ++m0) {
          const int row = mL + m0;
          if (row >= numMomsL) {
            continue;
          }
          for (int m1 = 0; m1 < effBatchR; ++m1) {
            const int col = mR + m1;
            if (col >= numMomsR) {
              continue;
            }
            const std::size_t wIdx =
                static_cast<std::size_t>(row) * static_cast<std::size_t>(numMomsR) +
                static_cast<std::size_t>(col);
            const Scalar& w = scalarWeight ? traceWeights[0] : traceWeights[wIdx];
            AccumulateWeightedDotPtr(mu.At(row, col), w, vecL.Column(m0), vecR.Column(m1), dim);
          }
        }
#endif
      }
    }
    hook(mu, 1);
  }
}

// State-blocked execution (state_block > 1): `stateBlock` trace states advance together
// through one MultiplyBlock per operator application; per state the floating-point sequence
// is identical to RunSerialExecute and the tensor is accumulated in state order.
void RunBlockedExecute(const LinearOperator& rescaledH,
                       const std::vector<const LinearOperator*>& operators,
                       StateGenerator& states, int batchSize, int stateBlock,
                       const ParsedPlan& pp, ChebyshevMoments& mu, CheckpointHook& hook) {
  const bool hasTimeAxis = pp.hasTimeAxis;
  const TimeAxis& timeAxisSpec = pp.timeAxisSpec;
  const SideProgram& ketProg = pp.ketProg;
  const SideProgram& braProg = pp.braProg;
  const bool seenMirror = pp.seenMirror;
  const std::vector<Scalar>& traceWeights = pp.traceWeights;
  const bool scalarWeight = pp.scalarWeight;
  const bool validMomentWeights = pp.validMomentWeights;

  const Index dim = rescaledH.Dimension();
  const int S = stateBlock;

  Vector r(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  VectorBlock rBlock = MakeBlock(dim, S);
  VectorBlock ketSeed = MakeBlock(dim, S);
  VectorBlock braSeed = MakeBlock(dim, S);
  VectorBlock auxSeed = MakeBlock(dim, S);
  VectorBlock ketInit = MakeBlock(dim, S);
  VectorBlock braInit = MakeBlock(dim, S);
  VectorBlock scratch = MakeBlock(dim, S);
  VectorBlock cheb0 = MakeBlock(dim, S);
  VectorBlock cheb1 = MakeBlock(dim, S);
  VectorBlock chebx0 = MakeBlock(dim, S);
  VectorBlock chebx1 = MakeBlock(dim, S);
  VectorBlock phiT = MakeBlock(dim, S);
  VectorBlock postVec = MakeBlock(dim, S);

  if (hasTimeAxis) {
    BlockRecursion rec(rescaledH, S);
    while (true) {
      int count = 0;
      for (int s = 0; s < S; ++s) {
        if (!states.Next(r)) break;
        std::copy_n(r.data(), static_cast<std::size_t>(dim), rBlock.Column(s));
        ++count;
      }
      if (count == 0) break;

      CopyBlock(rBlock, auxSeed);
      ApplySeedOpsBlock(rescaledH, operators, pp.seedKetOps, rBlock, ketSeed, scratch, cheb0,
                        cheb1, chebx0, chebx1);
      const VectorBlock& braSeedStart = pp.seedMirror ? ketSeed : rBlock;
      ApplySeedOpsBlock(rescaledH, operators, pp.seedBraOps, braSeedStart, braSeed, scratch,
                        cheb0, cheb1, chebx0, chebx1);

      for (int stepIdx = 0; stepIdx < timeAxisSpec.numSteps; ++stepIdx) {
        ApplySeedOpsBlock(rescaledH, operators, ketProg.preOps, ketSeed, ketInit, scratch,
                          cheb0, cheb1, chebx0, chebx1);
        const VectorBlock& braBodyStart = seenMirror ? ketInit : braSeed;
        ApplySeedOpsBlock(rescaledH, operators, braProg.preOps, braBodyStart, braInit, scratch,
                          cheb0, cheb1, chebx0, chebx1);

        if (ketProg.hasExpand) {
          const int numMoms = ketProg.numMoments;
          rec.Start(ketInit);
          for (int m = 0; m < numMoms; ++m) {
            const VectorBlock* ketVec = &rec.cur;
            if (!ketProg.postOperators.empty()) {
              ApplySingleVectorOperatorsBlock(operators, ketProg.postOperators, rec.cur,
                                              postVec, scratch);
              ketVec = &postVec;
            }
            for (int s = 0; s < count; ++s) {
              std::size_t wIdx = 0;
              if (validMomentWeights) {
                wIdx = static_cast<std::size_t>(m);
              } else if (!scalarWeight) {
                wIdx = static_cast<std::size_t>(m) *
                           static_cast<std::size_t>(timeAxisSpec.numSteps) +
                       static_cast<std::size_t>(stepIdx);
              }
              AccumulateWeightedDotPtr(mu.At(m, stepIdx), traceWeights[wIdx],
                                       braInit.Column(s), ketVec->Column(s), dim);
            }
            rec.Step();
          }
        } else {
          const std::size_t wIdx = scalarWeight ? 0 : static_cast<std::size_t>(stepIdx);
          for (int s = 0; s < count; ++s) {
            AccumulateWeightedDotPtr(mu.At(stepIdx), traceWeights[wIdx], braInit.Column(s),
                                     ketInit.Column(s), dim);
          }
        }

        for (const AxisStep& ev : timeAxisSpec.evolve) {
          if (const auto* er = std::get_if<EvolveRegister>(&ev)) {
            VectorBlock& targetVec = SelectRegisterBlock(er->target, ketSeed, braSeed, auxSeed);
            RunApplySeriesBlock(rescaledH, er->coeffs, er->scale, targetVec, cheb0, cheb1);
          } else if (const auto* ed = std::get_if<EvolveDerived>(&ev)) {
            const VectorBlock& srcVec =
                SelectRegisterBlock(ed->source, ketSeed, braSeed, auxSeed);
            VectorBlock& targetVec = SelectRegisterBlock(ed->target, ketSeed, braSeed, auxSeed);
            CopyBlock(srcVec, phiT);
            const LinearOperator& commOp =
                *operators[static_cast<std::size_t>(ed->operatorIndex)];
            RunApplyDerivedSeriesBlock(rescaledH, commOp, ed->coeffs, ed->scale, phiT, cheb0,
                                       cheb1, chebx0, chebx1);
            AxpyBlock(Scalar(1.0, 0.0), phiT, targetVec);
          } else if (const auto* cp = std::get_if<CopyRegister>(&ev)) {
            const VectorBlock& fromVec =
                SelectRegisterBlock(cp->from, ketSeed, braSeed, auxSeed);
            VectorBlock& toVec = SelectRegisterBlock(cp->to, ketSeed, braSeed, auxSeed);
            CopyBlock(fromVec, toVec);
          }
        }
      }
      hook(mu, count);
    }
    return;
  }

  if (!braProg.hasExpand && !ketProg.hasExpand) {
    // Rank 0.
    while (true) {
      int count = 0;
      for (int s = 0; s < S; ++s) {
        if (!states.Next(r)) break;
        std::copy_n(r.data(), static_cast<std::size_t>(dim), rBlock.Column(s));
        ++count;
      }
      if (count == 0) break;
      ApplySeedOpsBlock(rescaledH, operators, ketProg.preOps, rBlock, ketInit, scratch, cheb0,
                        cheb1, chebx0, chebx1);
      ApplySeedOpsBlock(rescaledH, operators, braProg.preOps, rBlock, braInit, scratch, cheb0,
                        cheb1, chebx0, chebx1);
      for (int s = 0; s < count; ++s) {
        AccumulateWeightedDotPtr(mu.Data()[0], traceWeights[0], braInit.Column(s),
                                 ketInit.Column(s), dim);
      }
      hook(mu, count);
    }
    return;
  }

  if (!braProg.hasExpand && ketProg.hasExpand) {
    // Rank 1 (ket expanded).
    const int numMoms = ketProg.numMoments;
    if (ketProg.postOperators.empty()) {
      BlockRecursion rec(rescaledH, S);
      while (true) {
        int count = 0;
        for (int s = 0; s < S; ++s) {
          if (!states.Next(r)) break;
          std::copy_n(r.data(), static_cast<std::size_t>(dim), rBlock.Column(s));
          ++count;
        }
        if (count == 0) break;
        ApplySeedOpsBlock(rescaledH, operators, ketProg.preOps, rBlock, ketInit, scratch,
                          cheb0, cheb1, chebx0, chebx1);
        ApplySeedOpsBlock(rescaledH, operators, braProg.preOps, rBlock, braInit, scratch,
                          cheb0, cheb1, chebx0, chebx1);
        rec.Start(ketInit);
        for (int m = 0; m < numMoms; ++m) {
          const Scalar& w = scalarWeight ? traceWeights[0]
                                         : traceWeights[static_cast<std::size_t>(m)];
          for (int s = 0; s < count; ++s) {
            AccumulateWeightedDotPtr(mu.At(m), w, braInit.Column(s), rec.cur.Column(s), dim);
          }
          rec.Step();
        }
        hook(mu, count);
      }
    } else {
      const int effBatch = std::min(batchSize, numMoms);
      const int paddedM = PaddedMoments(numMoms, effBatch);
      BlockRecursion rec(rescaledH, S);
      BlockChebWindow vecR(dim, effBatch, S);
      while (true) {
        int count = 0;
        for (int s = 0; s < S; ++s) {
          if (!states.Next(r)) break;
          std::copy_n(r.data(), static_cast<std::size_t>(dim), rBlock.Column(s));
          ++count;
        }
        if (count == 0) break;
        ApplySeedOpsBlock(rescaledH, operators, ketProg.preOps, rBlock, ketInit, scratch,
                          cheb0, cheb1, chebx0, chebx1);
        ApplySeedOpsBlock(rescaledH, operators, braProg.preOps, rBlock, braInit, scratch,
                          cheb0, cheb1, chebx0, chebx1);
        rec.Start(ketInit);
        for (int mR = 0; mR < paddedM; mR += effBatch) {
          vecR.Fill(rec);
          for (int opIdx : ketProg.postOperators) {
            vecR.ApplyToAll(*operators[static_cast<std::size_t>(opIdx)]);
          }
          for (int s = 0; s < count; ++s) {
            for (int m0 = 0; m0 < effBatch; ++m0) {
              const int m = mR + m0;
              if (m < numMoms) {
                const Scalar& w = scalarWeight ? traceWeights[0]
                                               : traceWeights[static_cast<std::size_t>(m)];
                AccumulateWeightedDotPtr(mu.At(m), w, braInit.Column(s),
                                         vecR.Column(s, m0), dim);
              }
            }
          }
        }
        hook(mu, count);
      }
    }
    return;
  }

  if (braProg.hasExpand && !ketProg.hasExpand) {
    // Rank 1 (bra expanded).
    const int numMoms = braProg.numMoments;
    const int effBatch = std::min(batchSize, numMoms);
    const int paddedM = PaddedMoments(numMoms, effBatch);
    BlockRecursion recL(rescaledH, S);
    BlockChebWindow vecL(dim, effBatch, S);
    while (true) {
      int count = 0;
      for (int s = 0; s < S; ++s) {
        if (!states.Next(r)) break;
        std::copy_n(r.data(), static_cast<std::size_t>(dim), rBlock.Column(s));
        ++count;
      }
      if (count == 0) break;
      ApplySeedOpsBlock(rescaledH, operators, ketProg.preOps, rBlock, ketInit, scratch, cheb0,
                        cheb1, chebx0, chebx1);
      ApplySeedOpsBlock(rescaledH, operators, braProg.preOps, rBlock, braInit, scratch, cheb0,
                        cheb1, chebx0, chebx1);
      recL.Start(braInit);
      for (int mL = 0; mL < paddedM; mL += effBatch) {
        vecL.Fill(recL);
        for (int opIdx : braProg.postOperators) {
          vecL.ApplyToAll(*operators[static_cast<std::size_t>(opIdx)]);
        }
        for (int s = 0; s < count; ++s) {
          for (int m0 = 0; m0 < effBatch; ++m0) {
            const int m = mL + m0;
            if (m < numMoms) {
              const Scalar& w = scalarWeight ? traceWeights[0]
                                             : traceWeights[static_cast<std::size_t>(m)];
              AccumulateWeightedDotPtr(mu.At(m), w, vecL.Column(s, m0), ketInit.Column(s),
                                       dim);
            }
          }
        }
      }
      hook(mu, count);
    }
    return;
  }

  // Rank 2 (both bra and ket expanded).
  const int numMomsL = braProg.numMoments;
  const int numMomsR = ketProg.numMoments;
  const int effBatchL = std::min(batchSize, numMomsL);
  const int effBatchR = std::min(batchSize, numMomsR);
  const int paddedML = PaddedMoments(numMomsL, effBatchL);
  const int paddedMR = PaddedMoments(numMomsR, effBatchR);

  BlockRecursion recL(rescaledH, S);
  BlockRecursion recR(rescaledH, S);
  BlockChebWindow vecL(dim, effBatchL, S);
  BlockChebWindow vecR(dim, effBatchR, S);

  while (true) {
    int count = 0;
    for (int s = 0; s < S; ++s) {
      if (!states.Next(r)) break;
      std::copy_n(r.data(), static_cast<std::size_t>(dim), rBlock.Column(s));
      ++count;
    }
    if (count == 0) break;
    ApplySeedOpsBlock(rescaledH, operators, ketProg.preOps, rBlock, ketSeed, scratch, cheb0,
                      cheb1, chebx0, chebx1);
    ApplySeedOpsBlock(rescaledH, operators, braProg.preOps, rBlock, braInit, scratch, cheb0,
                      cheb1, chebx0, chebx1);
    recL.Start(braInit);
    for (int mL = 0; mL < paddedML; mL += effBatchL) {
      vecL.Fill(recL);
      for (int opIdx : braProg.postOperators) {
        vecL.ApplyToAll(*operators[static_cast<std::size_t>(opIdx)]);
      }

      recR.Start(ketSeed);
      for (int mR = 0; mR < paddedMR; mR += effBatchR) {
        vecR.Fill(recR);
        for (int opIdx : ketProg.postOperators) {
          vecR.ApplyToAll(*operators[static_cast<std::size_t>(opIdx)]);
        }

        for (int s = 0; s < count; ++s) {
#ifdef LSQUANT_NATIVE_BUILD
          const int liveL = std::min(effBatchL, numMomsL - mL);
          const int liveR = std::min(effBatchR, numMomsR - mR);
          const auto table = linalg::MomentTableBlock(
              vecL.vectors.Data() + vecL.columnOffset(s), liveL,
              vecR.vectors.Data() + vecR.columnOffset(s), liveR, dim);
          for (int m = 0; m < liveL; ++m) {
            for (int n = 0; n < liveR; ++n) {
              const std::size_t index = static_cast<std::size_t>(mL + m) * numMomsR + mR + n;
              const Scalar& weight = scalarWeight ? traceWeights[0] : traceWeights[index];
              const Scalar& value = table[static_cast<std::size_t>(m) * liveR + n];
              if (weight.imag() == 0.0) mu.At(mL + m, mR + n) += weight.real() * value;
              else mu.At(mL + m, mR + n) += weight * value;
            }
          }
#else
          // Each row has one owner; dot products and trace-state additions keep
          // their serial order, so thread count cannot change the moment values.
#pragma omp parallel for schedule(static)
          for (int m0 = 0; m0 < effBatchL; ++m0) {
            const int row = mL + m0;
            if (row >= numMomsL) {
              continue;
            }
            for (int m1 = 0; m1 < effBatchR; ++m1) {
              const int col = mR + m1;
              if (col >= numMomsR) {
                continue;
              }
              const std::size_t wIdx =
                  static_cast<std::size_t>(row) * static_cast<std::size_t>(numMomsR) +
                  static_cast<std::size_t>(col);
              const Scalar& w = scalarWeight ? traceWeights[0] : traceWeights[wIdx];
              AccumulateWeightedDotPtr(mu.At(row, col), w, vecL.Column(s, m0),
                                       vecR.Column(s, m1), dim);
            }
          }
#endif
        }
      }
    }
    hook(mu, count);
  }
}

#ifdef LSQUANT_NATIVE_BUILD
// Native rank-2 path: the ket side is stored per state in a state-major
// `N x (bR*S)` block (column `s*bR + k` holds `T_{rBase+k}` of state `s` after the ket
// post-operators), while the bra side is streamed in windows of `leftWindow` vectors
// (state-major `N x (bL*S)`, column `s*bL + m0`). Each window is contracted against the stored
// ket block with the deterministic split-K table (`linalg::MomentTableBlock`). With
// `rightWindow = 0` the ket side is stored once (`bR = M`, the 0012 path); with
// `rightWindow = B < M` the ket recursion is restarted once per block of B vectors and the bra
// side is re-streamed per block, bounding the ket store to `B*N*16 B` at
// `~M^2/B` extra matvecs per state (amendment 2). Per state the GEMM shapes and the chunk grid
// do not depend on `S`, so `state_block = 1` and `S > 1` apply the identical per-state
// arithmetic.
void RunRank2Native(const LinearOperator& rescaledH,
                    const std::vector<const LinearOperator*>& operators,
                    StateGenerator& states, int stateBlock, int leftWindow, int rightWindow,
                    const ParsedPlan& pp, ChebyshevMoments& mu, CheckpointHook& hook) {
  const SideProgram& ketProg = pp.ketProg;
  const SideProgram& braProg = pp.braProg;
  const std::vector<Scalar>& traceWeights = pp.traceWeights;
  const bool scalarWeight = pp.scalarWeight;
  const Index dim = rescaledH.Dimension();
  const int S = stateBlock;
  const int numMomsL = braProg.numMoments;
  const int numMomsR = ketProg.numMoments;
  const int bL = std::max(1, std::min(leftWindow, numMomsL));
  // rightWindow = 0 ("all") stores the whole ket side; B > 0 stores blocks of B vectors.
  const int bR = (rightWindow <= 0 || rightWindow >= numMomsR) ? numMomsR : rightWindow;

  Vector r(static_cast<std::size_t>(dim), Scalar(0.0, 0.0));
  VectorBlock rBlock = MakeBlock(dim, S);
  VectorBlock ketSeed = MakeBlock(dim, S);
  VectorBlock braInit = MakeBlock(dim, S);
  VectorBlock scratch = MakeBlock(dim, S);
  // The series scratch registers are only touched by ApplySeries / ApplyDerivedSeries seed
  // ops (never present in the Kubo-style static rank-2 plans); allocate them lazily so the
  // default right_window config's footprint stays at the window sizes plus live registers.
  bool needSeries = false, needDerived = false;
  for (const auto* side : {&ketProg.preOps, &braProg.preOps}) {
    for (const SeedOp& op : *side) {
      needSeries = needSeries || std::holds_alternative<ApplySeries>(op) ||
                   std::holds_alternative<ApplyDerivedSeries>(op);
      needDerived = needDerived || std::holds_alternative<ApplyDerivedSeries>(op);
    }
  }
  VectorBlock cheb0 = needSeries ? MakeBlock(dim, S) : VectorBlock(0, 0);
  VectorBlock cheb1 = needSeries ? MakeBlock(dim, S) : VectorBlock(0, 0);
  VectorBlock chebx0 = needDerived ? MakeBlock(dim, S) : VectorBlock(0, 0);
  VectorBlock chebx1 = needDerived ? MakeBlock(dim, S) : VectorBlock(0, 0);
  VectorBlock postVec = MakeBlock(dim, S);

  VectorBlock rightStore(dim, bR * S);
  VectorBlock leftWin(dim, bL * S);
  Vector tableScratch(linalg::TableScratchSize(dim, bL, numMomsR), Scalar(0.0, 0.0));

  BlockRecursion recR(rescaledH, S);
  BlockRecursion recL(rescaledH, S);

  // Exact allocation report for the memory accounting (EngineDiagnostics; native only).
  if (EngineDiagnosticsSetting().memoryReport) {
    const double mb = 1.0 / 1048576.0;
    const double vecMB = static_cast<double>(dim) * 16.0 * mb;
    const double regMB = (6.0 * S + (needSeries ? 2.0 * S : 0.0) +
                          (needDerived ? 2.0 * S : 0.0)) * vecMB + vecMB /* r */;
    std::fprintf(stderr,
                 "ENGINE_MEM dim=%lld S=%d bL=%d bR=%d right_store=%.1fMB left_window=%.1fMB "
                 "registers=%.1fMB recursions=%.1fMB table_scratch=%.1fMB tensor=%.1fMB\n",
                 static_cast<long long>(dim), S, bL, bR,
                 bR * S * vecMB, bL * S * vecMB, regMB, 6.0 * S * vecMB,
                 tableScratch.size() * 16.0 * mb,
                 static_cast<double>(mu.Data().size()) * 16.0 * mb);
  }

  // Phase timing for the throughput work (EngineDiagnostics; native only).
  const bool timing = EngineDiagnosticsSetting().phaseTiming;
  double tSeed = 0.0, tRightMatvec = 0.0, tRightPost = 0.0, tRightScatter = 0.0;
  double tLeftMatvec = 0.0, tLeftPost = 0.0, tLeftScatter = 0.0, tTable = 0.0,
         tAccum = 0.0;
  auto tick = []() { return std::chrono::steady_clock::now(); };
  auto since = [](std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  };

  while (true) {
    int count = 0;
    for (int s = 0; s < S; ++s) {
      if (!states.Next(r)) {
        break;
      }
      std::copy_n(r.data(), static_cast<std::size_t>(dim), rBlock.Column(s));
      ++count;
    }
    if (count == 0) {
      break;
    }

    if (timing) {
      auto t0 = tick();
      ApplySeedOpsBlock(rescaledH, operators, ketProg.preOps, rBlock, ketSeed, scratch, cheb0,
                        cheb1, chebx0, chebx1);
      ApplySeedOpsBlock(rescaledH, operators, braProg.preOps, rBlock, braInit, scratch, cheb0,
                        cheb1, chebx0, chebx1);
      tSeed += since(t0);
    } else {
      ApplySeedOpsBlock(rescaledH, operators, ketProg.preOps, rBlock, ketSeed, scratch, cheb0,
                        cheb1, chebx0, chebx1);
      ApplySeedOpsBlock(rescaledH, operators, braProg.preOps, rBlock, braInit, scratch, cheb0,
                        cheb1, chebx0, chebx1);
    }

    // The ket side is stored in blocks of bR vectors: column s*bR+k holds
    // postOps(T_{rBase+k}(H~) ketSeed_s). The right recursion is restarted from the seed once
    // per block and stepped to the block end (vectors before the block are not stored). With
    // bR = numMomsR this is exactly the store-once path; with bR < numMomsR the right store is
    // bounded to bR*N*16 B per state at the cost of M^2/(2 bR) + M/2 right matvecs per state
    // (amendment 2).
    for (int rBase = 0; rBase < numMomsR; rBase += bR) {
      const int liveR = std::min(bR, numMomsR - rBase);
      recR.Start(ketSeed);
      for (int n = 0; n < rBase + liveR; ++n) {
        auto t0 = timing ? tick() : std::chrono::steady_clock::time_point{};
        if (n >= rBase) {
          const VectorBlock* src = ketProg.postOperators.empty()
                                       ? &recR.cur
                                       : ApplyPostOpsChain(operators, ketProg.postOperators,
                                                           recR.cur, postVec, scratch);
          if (timing) {
            tRightPost += since(t0);
            t0 = tick();
          }
          CopyColumnsGrid(*src, count, dim, rightStore.Column(n - rBase), bR);
          if (timing) {
            tRightScatter += since(t0);
            t0 = tick();
          }
        }
        recR.Step();
        if (timing) {
          tRightMatvec += since(t0);
        }
      }

      // Stream the bra side in windows of bL vectors; contract each window against the
      // stored ket block.
      recL.Start(braInit);
      for (int mBase = 0; mBase < numMomsL; mBase += bL) {
        const int live = std::min(bL, numMomsL - mBase);
        for (int m0 = 0; m0 < live; ++m0) {
          auto t0 = timing ? tick() : std::chrono::steady_clock::time_point{};
          const VectorBlock* src = braProg.postOperators.empty()
                                       ? &recL.cur
                                       : ApplyPostOpsChain(operators, braProg.postOperators,
                                                           recL.cur, postVec, scratch);
          if (timing) {
            tLeftPost += since(t0);
            t0 = tick();
          }
          CopyColumnsGrid(*src, count, dim, leftWin.Column(m0), bL);
          if (timing) {
            tLeftScatter += since(t0);
            t0 = tick();
          }
          recL.Step();
          if (timing) {
            tLeftMatvec += since(t0);
          }
        }
        for (int s = 0; s < count; ++s) {
          auto t0 = timing ? tick() : std::chrono::steady_clock::time_point{};
          const auto table = linalg::MomentTableBlock(
              leftWin.Data() + static_cast<std::size_t>(s) * static_cast<std::size_t>(bL) *
                                   static_cast<std::size_t>(dim),
              live,
              rightStore.Data() + static_cast<std::size_t>(s) *
                                      static_cast<std::size_t>(bR) *
                                      static_cast<std::size_t>(dim),
              liveR, dim, tableScratch.data());
          if (timing) {
            tTable += since(t0);
            t0 = tick();
          }
          for (int m = 0; m < live; ++m) {
            for (int n = 0; n < liveR; ++n) {
              const std::size_t index =
                  static_cast<std::size_t>(mBase + m) * static_cast<std::size_t>(numMomsR) +
                  static_cast<std::size_t>(rBase + n);
              const Scalar& weight = scalarWeight ? traceWeights[0] : traceWeights[index];
              const Scalar& value = table[static_cast<std::size_t>(m) *
                                              static_cast<std::size_t>(liveR) +
                                          static_cast<std::size_t>(n)];
              if (weight.imag() == 0.0) {
                mu.At(mBase + m, rBase + n) += weight.real() * value;
              } else {
                mu.At(mBase + m, rBase + n) += weight * value;
              }
            }
          }
          if (timing) {
            tAccum += since(t0);
          }
        }
      }
    }
    hook(mu, count);
  }
  if (timing) {
    std::fprintf(stderr,
                 "ENGINE_TIMING seed=%.3f right_matvec=%.3f right_post=%.3f "
                 "right_scatter=%.3f left_matvec=%.3f left_post=%.3f left_scatter=%.3f "
                 "table=%.3f accum=%.3f\n",
                 tSeed, tRightMatvec, tRightPost, tRightScatter, tLeftMatvec, tLeftPost,
                 tLeftScatter, tTable, tAccum);
  }
}
#endif  // LSQUANT_NATIVE_BUILD

}  // namespace

MomentEngine::MomentEngine(const LinearOperator& rescaledH,
                           std::vector<const LinearOperator*> operators,
                           StateGenerator& states, int batchSize, int stateBlock,
                           int leftWindow, int rightWindow)
    : rescaledH_(&rescaledH),
      operators_(std::move(operators)),
      states_(&states),
      batchSize_(batchSize),
      stateBlock_(stateBlock),
      leftWindow_(leftWindow),
      rightWindow_(rightWindow) {
  if (rescaledH.Dimension() <= 0) {
    throw std::invalid_argument("MomentEngine: rescaledH dimension must be positive");
  }
  if (states.Dimension() != rescaledH.Dimension()) {
    throw std::invalid_argument("MomentEngine: states dimension mismatch");
  }
  if (batchSize <= 0) {
    throw std::invalid_argument("MomentEngine: batchSize must be positive");
  }
  if (stateBlock <= 0) {
    throw std::invalid_argument("MomentEngine: stateBlock must be positive");
  }
  if (leftWindow <= 0) {
    throw std::invalid_argument("MomentEngine: leftWindow must be positive");
  }
  if (rightWindow < 0) {
    throw std::invalid_argument("MomentEngine: rightWindow must be non-negative (0 = all)");
  }
  for (const LinearOperator* op : operators_) {
    if (op == nullptr) {
      throw std::invalid_argument("MomentEngine: null operator pointer");
    }
    if (op->Dimension() != rescaledH.Dimension()) {
      throw std::invalid_argument("MomentEngine: operator dimension mismatch");
    }
  }
}

ChebyshevMoments MomentEngine::Run(const Plan& plan) { return Run(plan, ""); }

ChebyshevMoments MomentEngine::Run(const Plan& plan, const std::string& checkpointPath) {
  const ParsedPlan pp = ParsePlan(plan, static_cast<int>(operators_.size()));

  ChebyshevMoments mu(pp.dims);
  mu.systemSize = rescaledH_->Dimension();
  mu.numberOfStates = static_cast<int>(states_->NumberOfStates());
  mu.stateLabel = states_->Label();
  mu.word = plan.word;
  if (pp.hasTimeAxis) {
    mu.numberOfTimeSteps = pp.timeAxisSpec.numSteps;
  }

  CheckpointHook hook;
  hook.path = checkpointPath;
  if (!checkpointPath.empty()) {
    if (EngineDiagnosticsSetting().checkpointAbortAfter >= 0) {
      hook.abortAfter = static_cast<Index>(EngineDiagnosticsSetting().checkpointAbortAfter);
    }
  }

  states_->Reset();
  if (!checkpointPath.empty()) {
    Index statesDone = 0;
    if (TryLoadCheckpoint(checkpointPath, mu, statesDone)) {
      hook.statesDone = statesDone;
      states_->Skip(statesDone);
    }
  }
#ifdef LSQUANT_NATIVE_BUILD
  const bool staticRank2 = !pp.hasTimeAxis && pp.braProg.hasExpand && pp.ketProg.hasExpand;
  if (staticRank2) {
    RunRank2Native(*rescaledH_, operators_, *states_, stateBlock_, leftWindow_, rightWindow_,
                   pp, mu, hook);
    return mu;
  }
#endif
  if (stateBlock_ <= 1) {
    RunSerialExecute(*rescaledH_, operators_, *states_, batchSize_, pp, mu, hook);
  } else {
    RunBlockedExecute(*rescaledH_, operators_, *states_, batchSize_, stateBlock_, pp, mu, hook);
  }
  return mu;
}

}  // namespace lsquant
