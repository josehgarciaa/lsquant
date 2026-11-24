/**
 * @file moment_engine.hpp
 * @brief Bytecode execution engine evaluating Chebyshev moment tensors from a Plan.
 *
 * Executes the seven-instruction `Plan` state machine over trace vectors yielded by a
 * `StateGenerator`, reproducing the exact recursion, batching, and dot-product accumulation
 * order of legacy `SpectralMoments` and `CorrelationExpansionMoments`.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_MOMENT_ENGINE_HPP_
#define LSQUANT_CORE_CHEBYSHEV_MOMENT_ENGINE_HPP_

#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/state_generator.hpp"
#include "core/chebyshev/chebyshev_moments.hpp"
#include "core/chebyshev/plan.hpp"

namespace lsquant {

/**
 * @brief Diagnostic switches of the moment engine, set once by the caller before a run.
 *
 * The engine never reads the environment; the Python package reads `LSQ_ENGINE_MEM`,
 * `LSQ_ENGINE_TIMING` and `LSQ_CHECKPOINT_ABORT_AFTER` at import and passes them here through
 * `ConfigureEngineDiagnostics` (bindings).
 * - `memoryReport`: print the exact allocation of the native rank-2 path to stderr.
 * - `phaseTiming`: print the phase timings of the native rank-2 path to stderr.
 * - `checkpointAbortAfter`: when >= 0, raise `std::runtime_error` right after that many
 *   checkpoint writes (fault injection for the resume test); -1 disables it.
 */
struct EngineDiagnostics {
  bool memoryReport = false;
  bool phaseTiming = false;
  long long checkpointAbortAfter = -1;
};

/** @brief Process-wide diagnostic switches of the moment engine (see EngineDiagnostics). */
inline EngineDiagnostics& EngineDiagnosticsSetting() {
  static EngineDiagnostics settings;
  return settings;
}


/**
 * @brief Execution engine evaluating `Plan` programs into `ChebyshevMoments` tensors.
 *
 * Solves the problem of computing arbitrary rank-0, rank-1, and rank-2 Chebyshev moment
 * traces with a single representation-independent state machine. For each trace vector `|r>`
 * produced sequentially by `states`, builds the ket and bra sides according to `plan.code`
 * and accumulates `<bra | ket>` inner products in the exact block order of the legacy
 * solvers. Complexity per state is `O(M * cost(H~))` for rank 1 and
 * `O((M / batchSize) * M * cost(H~) + M^2 * N)` for rank 2.
 *
 * @warning All operators and the state generator must share `rescaledH.Dimension()`, and
 *          `batchSize` must be positive.
 * @par Pitfall
 *      `Run` resets `states` at the start of execution so repeated calls on the same engine
 *      traverse the full state sequence deterministically.
 */
class MomentEngine {
public:
  /**
   * @brief Construct a `MomentEngine` bound to `rescaledH`, `operators`, `states`, and `batchSize`.
   *
   * @param rescaledH Rescaled Hamiltonian `H~` with spectrum in `[-1, 1]`.
   * @param operators Positional table of static `LinearOperator` pointers referenced by `ApplyOperator`.
   * @param states    Trace-vector generator yielding source states `|r>`.
   * @param batchSize Number of Chebyshev vectors buffered per batch during `Expand` (`B > 0`).
   * @param rightWindow Native rank-2 ket-side storage window: `0` stores all `M` right vectors
   *        per state (`"all"`, the fastest mode); `B > 0` stores `B` vectors and restarts the
   *        right recursion once per block (`M + M^2/B` matvecs per state instead of `2M`,
   *        memory `~B*N*16 B`). Native profile only.
   *
   * @warning Throws `std::invalid_argument` if `batchSize <= 0`, if any operator pointer is
   *          null, or if any dimension differs from `rescaledH.Dimension()`.
   * @par Pitfall
   *      When `Expand(M)` has `M` not divisible by `batchSize`, `M` is rounded up to the next
   *      multiple of `batchSize` during recursion and truncated on accumulation, matching
   *      legacy `CorrelationExpansionMoments`.
   */
  MomentEngine(const LinearOperator& rescaledH, std::vector<const LinearOperator*> operators,
               StateGenerator& states, int batchSize, int stateBlock = 1,
               int leftWindow = 32, int rightWindow = 0);

  /**
   * @brief Execute `plan` across all trace states and return the accumulated `ChebyshevMoments`.
   *
   * Validates `plan.code`, rewinds `states` via `Reset()`, and runs the ket/bra state machine
   * over the trace vectors. With `stateBlock == 1` the states are processed one at a time
   * (legacy order). With `stateBlock > 1`, `stateBlock` states advance together through one
   * `MultiplyBlock` per operator application; the moments of each state are computed with
   * exactly the per-state arithmetic of the `stateBlock == 1` path and accumulated into the
   * tensor in state order, so both settings produce bit-identical tensors. Populates
   * `systemSize`, `numberOfStates`, `stateLabel`, and `word` on the returned
   * `ChebyshevMoments`; the caller sets `window`.
   *
   * @param plan Instruction sequence ending with `Trace`.
   * @return     Accumulated `ChebyshevMoments` tensor of rank equal to the number of `Expand`
   *             instructions in `plan`.
   *
   * @warning Throws `std::runtime_error` if rank > 2 (a second `Expand` on the same side),
   *          and `std::invalid_argument` if `plan` is malformed.
   * @par Pitfall
   *      Rank-2 plans with `Trace(1)` return raw state sums; call `ApplyKpmWeights` afterwards
   *      when reproducing legacy `CorrelationExpansionMoments`.
   */
  ChebyshevMoments Run(const Plan& plan);

  /**
   * @brief Execute `plan` like `Run`, with trace-state-level checkpointing at `checkpointPath`.
   *
   * If a checkpoint file at `checkpointPath` matches this run (rank, dims, total states, state
   * label, word), the partial tensor is loaded, the generator is advanced by the recorded
   * number of finished states (`StateGenerator::Skip`, bit-exact), and execution continues;
   * the result is bit-identical to an uninterrupted run. After every finished trace state
   * (state block when `stateBlock > 1`) the partial tensor is rewritten atomically. An empty
   * `checkpointPath` disables checkpointing entirely. With
   * `EngineDiagnosticsSetting().checkpointAbortAfter = k >= 0` the run raises
   * `std::runtime_error` right after the k-th checkpoint write (fault-injection hook for the
   * resume test).
   *
   * @param plan           Instruction sequence ending with `Trace`.
   * @param checkpointPath Checkpoint file path, or empty to disable.
   * @return               Accumulated `ChebyshevMoments` tensor.
   *
   * @warning Throws `std::runtime_error` on a mismatched or corrupt checkpoint, and on any
   *          checkpoint I/O failure.
   */
  ChebyshevMoments Run(const Plan& plan, const std::string& checkpointPath);

private:
  const LinearOperator* rescaledH_;
  std::vector<const LinearOperator*> operators_;
  StateGenerator* states_;
  int batchSize_;
  int stateBlock_;
  int leftWindow_;  // streamed bra-side window of the native rank-2 path (native profile only)
  int rightWindow_;  // stored ket-side window of the native rank-2 path: 0 stores all M
                     // vectors ("all"); B > 0 bounds the stored right side to B vectors per
                     // state, restarting the right recursion once per block (amendment 2)
};

namespace chebyshev {
using lsquant::MomentEngine;
}  // namespace chebyshev

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_MOMENT_ENGINE_HPP_
