/**
 * @file test_core_stage1.cpp
 * @brief Unit tests of the stage-1 core: CSR algebra, CSR text I/O, spectral bounds,
 *        spectral window, Chebyshev recursion, state generators and the moment engine.
 *
 * Every test compares one component against an independent closed form or a dense
 * Eigen computation on a small matrix, so a failure points at one class:
 *
 *   - SparseMatrix: y = a A x + b y against the dense product; Adjoint, IsHermitian.
 *   - ReadCsrText / WriteCsrText: write, read back, same action on a vector.
 *   - GershgorinBounds / LanczosBounds: contain (Gershgorin) or approach (Lanczos) the
 *     dense extremal eigenvalues.
 *   - SpectralWindow: maps [lower, upper] onto [-1, 1].
 *   - ChebyshevRecursion: on a diagonal H~ the m-th vector is cos(m arccos d_i) phi_i.
 *   - RandomPhaseStates / LocalStates: unit norm, determinism for a fixed seed,
 *     different vectors for a different seed, basis vectors for LocalStates.
 *   - MomentEngine, rank 1: Tr[T_m(H~)] over the full local basis equals the sum over
 *     the dense eigenvalues of T_m.
 *   - MomentEngine, rank 2: Tr[A T_n(H~) A T_m(H~)] against a dense Chebyshev recursion,
 *     and the same tensor for every batch size (batching must not change the result).
 *
 * Exit code 0 when every check passes; the failing checks are printed otherwise.
 * Registered with CTest as `core_stage1_unit`.
 */
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/spectral_bounds.hpp"
#include "core/algebra/state_generator.hpp"
#include "core/chebyshev/chebyshev_moments.hpp"
#include "core/chebyshev/chebyshev_recursion.hpp"
#include "core/chebyshev/moment_engine.hpp"
#include "core/chebyshev/plan.hpp"
#include "core/chebyshev/rescaled_operator.hpp"
#include "core/chebyshev/spectral_window.hpp"
#include "io/csr_text.hpp"

namespace {

int g_failures = 0;
int g_checks = 0;

// Record one check; print the failing expression with its location.
void Check(bool ok, const std::string& what, int line) {
  ++g_checks;
  if (!ok) {
    ++g_failures;
    std::cerr << "FAIL line " << line << ": " << what << "\n";
  }
}
#define CHECK(cond) Check((cond), #cond, __LINE__)
#define CHECK_CLOSE(a, b, tol) \
  Check(std::abs((a) - (b)) <= (tol), std::string(#a " ~ " #b " (tol " #tol ")"), __LINE__)

using lsquant::Index;
using lsquant::Scalar;
using lsquant::Vector;

// Random complex Hermitian matrix of dimension `dim` with about `density` of the
// off-diagonal pairs filled; returns the CSR matrix and its dense copy.
lsquant::SparseMatrix RandomHermitian(Index dim, double density, unsigned seed,
                                      Eigen::MatrixXcd& dense) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  dense = Eigen::MatrixXcd::Zero(dim, dim);
  std::vector<Index> rows, cols;
  std::vector<Scalar> vals;
  for (Index i = 0; i < dim; ++i) {
    const Scalar d(u(rng), 0.0);
    dense(i, i) = d;
    rows.push_back(i); cols.push_back(i); vals.push_back(d);
    for (Index j = i + 1; j < dim; ++j) {
      if (std::abs(u(rng)) < density) {
        const Scalar v(u(rng), u(rng));
        dense(i, j) = v;
        dense(j, i) = std::conj(v);
        rows.push_back(i); cols.push_back(j); vals.push_back(v);
        rows.push_back(j); cols.push_back(i); vals.push_back(std::conj(v));
      }
    }
  }
  return lsquant::SparseMatrix::FromTriplets(dim, rows, cols, vals);
}

Eigen::VectorXcd ToEigen(const Vector& v) {
  Eigen::VectorXcd e(static_cast<Index>(v.size()));
  for (std::size_t i = 0; i < v.size(); ++i) e(static_cast<Index>(i)) = v[i];
  return e;
}

Vector RandomVector(Index dim, unsigned seed) {
  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  Vector v(static_cast<std::size_t>(dim));
  for (auto& x : v) x = Scalar(u(rng), u(rng));
  return v;
}

// ------------------------------------------------------------------ SparseMatrix
void TestSparseMultiply() {
  Eigen::MatrixXcd A;
  const auto csr = RandomHermitian(20, 0.3, 1, A);
  const Vector x = RandomVector(20, 2);
  Vector y = RandomVector(20, 3);
  const Eigen::VectorXcd y0 = ToEigen(y);
  const Scalar a(0.7, -0.2), b(-1.3, 0.4);
  csr.Multiply(a, x, b, y);  // y <- a A x + b y
  const Eigen::VectorXcd ref = a * (A * ToEigen(x)) + b * y0;
  CHECK((ToEigen(y) - ref).norm() < 1e-13 * ref.norm());
  CHECK(csr.Dimension() == 20);
  CHECK(csr.IsHermitian());
}

void TestSparseAdjointAndHermitian() {
  // non-Hermitian 2 x 2: [[0, 1], [0, 0]]
  const auto m = lsquant::SparseMatrix::FromTriplets(2, {0}, {1}, {Scalar(1.0, 0.0)});
  CHECK(!m.IsHermitian());
  const auto mt = m.Adjoint();
  Vector e0{Scalar(1.0, 0.0), Scalar(0.0, 0.0)}, y(2);
  mt.Multiply(Scalar(1.0, 0.0), e0, Scalar(0.0, 0.0), y);  // A^+ e0 = e1
  CHECK_CLOSE(y[1], Scalar(1.0, 0.0), 1e-15);
  CHECK_CLOSE(y[0], Scalar(0.0, 0.0), 1e-15);
}

// ------------------------------------------------------------------ CSR text I/O
void TestCsrTextRoundTrip() {
  Eigen::MatrixXcd A;
  const auto csr = RandomHermitian(15, 0.4, 4, A);
  const std::string path = "test_core_stage1_roundtrip.CSR";
  lsquant::WriteCsrText(path, csr);
  const auto back = lsquant::ReadCsrText(path);
  std::remove(path.c_str());
  CHECK(back.Dimension() == csr.Dimension());
  CHECK(back.NonZeros() == csr.NonZeros());
  const Vector x = RandomVector(15, 5);
  Vector y1(15), y2(15);
  csr.Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), y1);
  back.Multiply(Scalar(1.0, 0.0), x, Scalar(0.0, 0.0), y2);
  CHECK((ToEigen(y1) - ToEigen(y2)).norm() < 1e-12 * ToEigen(y1).norm());
}

// ------------------------------------------------------------------ spectral bounds
void TestSpectralBounds() {
  Eigen::MatrixXcd A;
  const auto csr = RandomHermitian(40, 0.2, 6, A);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);
  const double emin = es.eigenvalues().minCoeff(), emax = es.eigenvalues().maxCoeff();

  const auto g = lsquant::GershgorinBounds(csr, 0.10);
  CHECK(g.lower <= emin);
  CHECK(g.upper >= emax);

  const lsquant::SparseMatrixOperator op(csr);
  const auto l = lsquant::LanczosBounds(op, 40, 12345);
  // 40 Lanczos steps on a 40 x 40 matrix give the extremal eigenvalues essentially exactly
  CHECK_CLOSE(l.lower, emin, 1e-8);
  CHECK_CLOSE(l.upper, emax, 1e-8);
}

// ------------------------------------------------------------------ spectral window
void TestSpectralWindow() {
  const lsquant::SpectralWindow w{-3.0, 5.0};
  CHECK_CLOSE(w.HalfWidth(), 4.0, 1e-15);
  CHECK_CLOSE(w.Center(), 1.0, 1e-15);
  CHECK_CLOSE(w.Scale() * w.lower + w.Shift(), -1.0, 1e-15);
  CHECK_CLOSE(w.Scale() * w.upper + w.Shift(), 1.0, 1e-15);
}

// ------------------------------------------------------------------ Chebyshev recursion
void TestChebyshevRecursionDiagonal() {
  // H = diag(d_i) inside the window [-2, 2]: H~ = H / 2, T_m(H~) phi = cos(m acos(d_i/2)) phi_i
  const Index n = 9;
  std::vector<Index> idx;
  std::vector<Scalar> vals;
  std::vector<double> d;
  for (Index i = 0; i < n; ++i) {
    d.push_back(-1.9 + 3.8 * static_cast<double>(i) / static_cast<double>(n - 1));
    idx.push_back(i);
    vals.push_back(Scalar(d.back(), 0.0));
  }
  const auto H = lsquant::SparseMatrix::FromTriplets(n, idx, idx, vals);
  const lsquant::RescaledOperator Ht(H, lsquant::SpectralWindow{-2.0, 2.0});
  lsquant::ChebyshevRecursion rec(Ht);
  const Vector phi = RandomVector(n, 7);
  rec.Start(phi);
  CHECK(rec.Order() == 0);
  double worst = 0.0;
  for (int m = 0; m <= 60; ++m) {
    if (m > 0) rec.Step();
    CHECK(rec.Order() == m);
    for (Index i = 0; i < n; ++i) {
      const Scalar ref = std::cos(m * std::acos(d[i] / 2.0)) * phi[i];
      worst = std::max(worst, std::abs(rec.Current()[i] - ref));
    }
  }
  CHECK(worst < 1e-12);
}

// ------------------------------------------------------------------ state generators
void TestRandomPhaseStates() {
  lsquant::RandomPhaseStates a(50, 3, 12345), b(50, 3, 12345), c(50, 3, 999);
  CHECK(a.NumberOfStates() == 3);
  CHECK(a.Dimension() == 50);
  Vector va, vb, vc;
  int count = 0;
  while (a.Next(va)) {
    CHECK(b.Next(vb));
    CHECK(c.Next(vc));
    double norm2 = 0.0;
    for (const auto& x : va) norm2 += std::norm(x);
    CHECK_CLOSE(norm2, 1.0, 1e-13);          // unit norm
    CHECK(va == vb);                           // same seed -> identical
    CHECK((ToEigen(va) - ToEigen(vc)).norm() > 1e-3);  // other seed -> different
    ++count;
  }
  CHECK(count == 3);
  a.Reset();  // Reset restarts the same sequence
  Vector first;
  CHECK(a.Next(first));
  lsquant::RandomPhaseStates fresh(50, 3, 12345);
  Vector ref;
  fresh.Next(ref);
  CHECK(first == ref);
}

void TestLocalStates() {
  lsquant::LocalStates s(6, {0, 4});
  Vector v;
  CHECK(s.Next(v));
  CHECK(v.size() == 6 && v[0] == Scalar(1.0, 0.0) && v[4] == Scalar(0.0, 0.0));
  CHECK(s.Next(v));
  CHECK(v[4] == Scalar(1.0, 0.0) && v[0] == Scalar(0.0, 0.0));
  CHECK(!s.Next(v));
}

// ------------------------------------------------------------------ moment engine
std::vector<Index> AllSites(Index n) {
  std::vector<Index> s(static_cast<std::size_t>(n));
  for (Index i = 0; i < n; ++i) s[static_cast<std::size_t>(i)] = i;
  return s;
}

void TestEngineRank1ExactTrace() {
  Eigen::MatrixXcd A;
  const Index n = 30;
  const auto H = RandomHermitian(n, 0.25, 8, A);
  const lsquant::SpectralWindow w{-6.0, 6.0};
  const lsquant::RescaledOperator Ht(H, w);
  lsquant::LocalStates states(n, AllSites(n));
  lsquant::MomentEngine engine(Ht, {}, states, 7);
  const int M = 40;
  lsquant::Plan plan;
  plan.word = "1 delta";
  plan.code.push_back(lsquant::Expand{M});
  plan.code.push_back(lsquant::Trace{{Scalar(1.0, 0.0)}});
  const auto mu = engine.Run(plan);
  CHECK(mu.Rank() == 1);
  CHECK(mu.Dims().size() == 1 && mu.Dims()[0] == M);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);
  double worst = 0.0;
  for (int m = 0; m < M; ++m) {
    double ref = 0.0;  // Tr T_m(H~) = sum_k cos(m acos(e~_k))
    for (Index k = 0; k < n; ++k) {
      const double et = w.Scale() * es.eigenvalues()(k) + w.Shift();
      ref += std::cos(m * std::acos(et));
    }
    worst = std::max(worst, std::abs(mu.At(m) - Scalar(ref, 0.0)));
  }
  CHECK(worst < 1e-10);
}

// Dense reference of mu(n, m) = Tr[A T_n(H~) A T_m(H~)].
std::vector<Eigen::MatrixXcd> DenseChebyshev(const Eigen::MatrixXcd& Ht, int M) {
  std::vector<Eigen::MatrixXcd> T(static_cast<std::size_t>(M));
  const Index n = Ht.rows();
  T[0] = Eigen::MatrixXcd::Identity(n, n);
  if (M > 1) T[1] = Ht;
  for (int m = 2; m < M; ++m) T[m] = 2.0 * Ht * T[m - 1] - T[m - 2];
  return T;
}

lsquant::ChebyshevMoments RunRank2(const lsquant::RescaledOperator& Ht,
                                   const lsquant::SparseMatrixOperator& Aop, Index n, int M,
                                   int batch) {
  lsquant::LocalStates states(n, AllSites(n));
  lsquant::MomentEngine engine(Ht, {&Aop}, states, batch);
  lsquant::Plan plan;
  plan.word = "A delta A delta";
  plan.code.push_back(lsquant::Expand{M});
  plan.code.push_back(lsquant::ApplyOperator{0});
  plan.code.push_back(lsquant::MirrorToBra{});
  plan.code.push_back(lsquant::ApplyOperator{0});
  plan.code.push_back(lsquant::Expand{M});
  plan.code.push_back(lsquant::Trace{{Scalar(1.0, 0.0)}});
  return engine.Run(plan);
}

void TestEngineRank2AndBatching() {
  Eigen::MatrixXcd Hd, Ad;
  const Index n = 24;
  const auto H = RandomHermitian(n, 0.3, 9, Hd);
  const auto A = RandomHermitian(n, 0.3, 10, Ad);
  const lsquant::SpectralWindow w{-6.0, 6.0};
  const lsquant::RescaledOperator Ht(H, w);
  const lsquant::SparseMatrixOperator Aop(A);
  const int M = 12;

  const auto mu = RunRank2(Ht, Aop, n, M, M);
  CHECK(mu.Rank() == 2);

  const Eigen::MatrixXcd Htd =
      w.Scale() * Hd + w.Shift() * Eigen::MatrixXcd::Identity(n, n);
  const auto T = DenseChebyshev(Htd, M);
  double worst = 0.0, scale = 0.0;
  for (int a = 0; a < M; ++a) {
    for (int b = 0; b < M; ++b) {
      const Scalar ref = (Ad * T[a] * Ad * T[b]).trace();
      scale = std::max(scale, std::abs(ref));
      worst = std::max(worst, std::abs(mu.At(a, b) - ref));
    }
  }
  CHECK(worst < 1e-11 * scale);

  // batching (M = 12 with B = 5 is not a multiple; B = 1 is the extreme) must not change mu
  for (int batch : {1, 5, 12}) {
    const auto mb = RunRank2(Ht, Aop, n, M, batch);
    double diff = 0.0;
    for (std::size_t k = 0; k < mu.Data().size(); ++k) {
      diff = std::max(diff, std::abs(mb.Data()[k] - mu.Data()[k]));
    }
    Check(diff == 0.0, "rank-2 moments identical for batch " + std::to_string(batch), __LINE__);
  }
}

void TestEngineRejectsBadInput() {
  Eigen::MatrixXcd A;
  const auto H = RandomHermitian(5, 0.5, 11, A);
  const lsquant::RescaledOperator Ht(H, lsquant::SpectralWindow{-4.0, 4.0});
  lsquant::LocalStates states(5, AllSites(5));
  bool threw = false;
  try {
    lsquant::MomentEngine bad(Ht, {}, states, 0);  // batch size must be > 0
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CHECK(threw);

  lsquant::MomentEngine engine(Ht, {}, states, 4);
  lsquant::Plan noTrace;
  noTrace.code.push_back(lsquant::Expand{4});
  threw = false;
  try {
    engine.Run(noTrace);  // a plan must end with Trace
  } catch (const std::exception&) {
    threw = true;
  }
  CHECK(threw);
}

}  // namespace

int main() {
  TestSparseMultiply();
  TestSparseAdjointAndHermitian();
  TestCsrTextRoundTrip();
  TestSpectralBounds();
  TestSpectralWindow();
  TestChebyshevRecursionDiagonal();
  TestRandomPhaseStates();
  TestLocalStates();
  TestEngineRank1ExactTrace();
  TestEngineRank2AndBatching();
  TestEngineRejectsBadInput();
  std::cout << "core_stage1_unit: " << g_checks - g_failures << "/" << g_checks
            << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
