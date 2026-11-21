/**
 * @file linalg.hpp
 * @brief Dense complex vector linear algebra primitives backed by Eigen::Map.
 *
 * Provides thin, zero-copy wrappers around `Eigen::Map` that reproduce the exact arithmetic
 * and summation order of the legacy LSQUANT `linalg` namespace (`axpy`, `scal`, `copy`,
 * `dot`, `norm2`).
 */
#ifndef LSQUANT_CORE_ALGEBRA_LINALG_HPP_
#define LSQUANT_CORE_ALGEBRA_LINALG_HPP_

#include "core/algebra/types.hpp"

namespace lsquant {
namespace linalg {

/**
 * @brief Compute the in-place complex vector update `y = a * x + y`.
 *
 * Maps `x` and `y` as Eigen column vectors and evaluates `eig_y = a * eig_x + eig_y` in
 * `O(N)` time.
 *
 * @param a Multiplicative scalar factor for `x`.
 * @param x Input vector of length `N`.
 * @param y Input/output vector of length `N`.
 *
 * @warning Throws `std::invalid_argument` if `x.size() != y.size()`.
 * @par Pitfall
 *      Both vectors must be pre-allocated to the same length before calling `axpy`.
 */
void axpy(Scalar a, const Vector& x, Vector& y);

/**
 * @brief Scale a complex vector in place: `x = a * x`.
 *
 * Maps `x` as an Eigen vector and multiplies every entry by `a` in `O(N)` time.
 *
 * @param a Scalar multiplier.
 * @param x Input/output vector scaled in place.
 *
 * @warning Scaling by `0` overwrites existing values but can propagate NaNs if `x` is
 *          uninitialised.
 * @par Pitfall
 *      Prefer `std::fill` when zeroing a freshly allocated vector that may contain non-finite
 *      values.
 */
void scal(Scalar a, Vector& x);

/**
 * @brief Copy complex vector `x` into pre-allocated vector `y`.
 *
 * Maps `x` and `y` as Eigen vectors and assigns `eig_y = eig_x` in `O(N)` time, matching
 * legacy `linalg::copy`.
 *
 * @param x Source vector of length `N`.
 * @param y Destination vector of length `N`.
 *
 * @warning Throws `std::invalid_argument` if `x.size() != y.size()`.
 * @par Pitfall
 *      Unlike `std::vector::operator=`, `linalg::copy` requires `y` to be pre-sized to match
 *      `x` and throws on mismatch.
 */
void copy(const Vector& x, Vector& y);

/**
 * @brief Compute the sesquilinear inner product `conj(x) . y = sum_i conj(x[i]) * y[i]`.
 *
 * Evaluates `Eigen::Map(x).dot(Eigen::Map(y))` sequentially on a single thread in `O(N)`
 * time, matching legacy `linalg::vdot(x, y)` bit for bit.
 *
 * @param x Left-hand (bra) vector of length `N`, conjugated in the product.
 * @param y Right-hand (ket) vector of length `N`.
 * @return  Complex inner product `<x | y>`.
 *
 * @warning Throws `std::invalid_argument` if `x.size() != y.size()`.
 * @par Pitfall
 *      Do not parallelise this reduction with OpenMP; changing the summation tree alters
 *      floating-point rounding and breaks byte identity with the golden moment tables.
 */
Scalar dot(const Vector& x, const Vector& y);

/**
 * @brief Compute `conj(x) . y` on raw contiguous buffers of length `n`.
 *
 * Identical arithmetic to `dot(const Vector&, const Vector&)` (single-threaded
 * `Eigen::Map::dot` in index order), taking raw pointers so columns of a `VectorBlock` can
 * be dotted without copying.
 *
 * @param x Left-hand (bra) buffer of length `n`, conjugated in the product.
 * @param y Right-hand (ket) buffer of length `n`.
 * @param n Buffer length `n >= 0`.
 * @return  Complex inner product `<x | y>`.
 */
Scalar dot(const Scalar* x, const Scalar* y, Index n);

/**
 * @brief Compute the squared Euclidean norm `||x||_2^2 = sum_i |x[i]|^2`.
 *
 * Evaluates `Eigen::Map(x).squaredNorm()` in `O(N)` time, matching legacy `linalg::nrm2`.
 *
 * @param x Input vector of length `N`.
 * @return  Non-negative real squared norm `||x||_2^2`.
 *
 * @warning Returns the squared norm `||x||_2^2`, not `||x||_2`.
 * @par Pitfall
 *      Take `std::sqrt(linalg::norm2(x))` when the Euclidean length itself is needed (e.g.
 *      for state normalisation in Lanczos).
 */
double norm2(const Vector& x);

#ifdef LSQUANT_NATIVE_BUILD
/**
 * @brief Parallel, thread-count-invariant inner product `conj(x) . y` (native profile only).
 *
 * Splits the index range into fixed chunks of `2^17` elements (independent of the OpenMP
 * thread count), computes each chunk's partial product in parallel, and sums the partials in
 * chunk order. The result is identical for any thread count; it differs from Eigen's serial
 * `dot` at the last-ulp level because the summation tree changes.
 *
 * @param x Left-hand (bra) buffer of length `n`, conjugated in the product.
 * @param y Right-hand (ket) buffer of length `n`.
 * @param n Buffer length `n >= 0`.
 * @return  Complex inner product `<x | y>`.
 */
Scalar DotDeterministic(const Scalar* x, const Scalar* y, Index n);
#endif

/**
 * @brief Compute the inner-product table `<left[m] | right[n]>` in row-major order.
 */
std::vector<Scalar> MomentTable(const std::vector<const Vector*>& left,
                                const std::vector<const Vector*>& right);

/**
 * @brief Compute the inner-product table of two column blocks, `<left_m | right_n>`.
 *
 * `left` and `right` point to column-major blocks of `leftCount` / `rightCount` contiguous
 * columns of length `dimension`. Returns `leftCount * rightCount` values in row-major order
 * (`result[m * rightCount + n] = <left_m | right_n>`). Equivalent to `MomentTable` on the
 * column pointers, without assembling them as separate vectors.
 *
 * @param left       Pointer to the first column of the left block.
 * @param leftCount  Number of left columns `bL >= 0`.
 * @param right      Pointer to the first column of the right block.
 * @param rightCount Number of right columns `bR >= 0`.
 * @param dimension  Column length `N >= 0`.
 * @param scratch    Optional workspace for the native-profile split-K partials; when given
 *                   it must hold at least `TableScratchSize(dimension, bL, bR)` scalars.
 *                   Passing a persistent buffer avoids a fresh page-faulting allocation per
 *                   call (the engine issues one table per left window).
 * @return           Row-major `bL x bR` table of inner products.
 *
 * @warning `scratch` is written by OpenMP workers; it must not be shared between concurrent
 *          calls.
 */
std::vector<Scalar> MomentTableBlock(const Scalar* left, int leftCount, const Scalar* right,
                                     int rightCount, Index dimension, Scalar* scratch = nullptr);

/**
 * @brief Number of scalars the native split-K `MomentTableBlock` writes to `scratch`.
 *
 * `ceil(N / chunkRows)` partial tables of `leftCount * rightCount` entries, with the fixed
 * chunk grid of the native profile (`max(4096, ceil(N / 256))` rows per chunk).
 */
std::size_t TableScratchSize(Index dimension, int leftCount, int rightCount);

}  // namespace linalg
}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_LINALG_HPP_
