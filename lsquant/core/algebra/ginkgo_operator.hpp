/**
 * @file ginkgo_operator.hpp
 * @brief Optional Ginkgo-backed LinearOperator for the native build profile.
 *
 * Wraps a CSR `SparseMatrix` in a Ginkgo matrix (csr, sellp, ell, hybrid, or coo) on an
 * OpenMP or CUDA executor so that the Chebyshev engine can run its `Multiply` and
 * `MultiplyBlock` kernels through Ginkgo's SpMV/SpMM implementations.
 */
#ifdef LSQUANT_WITH_GINKGO
#ifndef LSQUANT_CORE_ALGEBRA_GINKGO_OPERATOR_HPP_
#define LSQUANT_CORE_ALGEBRA_GINKGO_OPERATOR_HPP_

#include <memory>
#include <string>

#include <ginkgo/ginkgo.hpp>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/sparse_matrix.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief `LinearOperator` backed by a Ginkgo sparse matrix on an `omp` or `cuda` executor.
 *
 * Solves the problem of running the Chebyshev recursion and static words through Ginkgo's
 * sparse kernels without changing the engine. Contiguous columns are wrapped as
 * zero-copy N x 1 Dense views because Ginkgo Dense storage is row-major. The matrix
 * index type is int32 when the dimension and non-zero count allow it, int64 otherwise.
 *
 * @warning Available only in the native build profile with `LSQUANT_WITH_GINKGO=ON`;
 *          arithmetic order differs from `SparseMatrix::Multiply`, so results agree with
 *          the golden profile only to floating-point reduction accuracy.
 * @par Pitfall
 *      Column-major host blocks use one sparse apply per column. The device engine
 *      uses executor-resident row-major blocks for a single SpMM across states.
 */
class GinkgoOperator : public LinearOperator {
public:
  /**
   * @brief Build a Ginkgo operator from a CSR `SparseMatrix`.
   *
   * @param matrix   Source sparse matrix (copied into Ginkgo storage).
   * @param executor Executor name: `"omp"` or `"cuda"` (cuda requires a CUDA-enabled
   *                 Ginkgo build).
   * @param format   Matrix format: `"csr"`, `"sellp"`, `"ell"`, `"hybrid"`, or `"coo"`.
   *
   * @warning Throws `std::invalid_argument` for unknown executor or format names, and
   *          `std::runtime_error` if the CUDA executor is requested in a build without
   *          CUDA support.
   * @par Pitfall
   *      The source matrix is converted once at construction; later changes to `matrix`
   *      do not propagate.
   */
  explicit GinkgoOperator(const SparseMatrix& matrix, std::string executor = "omp",
                          std::string format = "csr");

  /**
   * @brief Return the operator dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Compute `y = a * O * x + b * y` with zero-copy Ginkgo dense views.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  /**
   * @brief Compute `Y = a * O * X + b * Y` on an `N x b` block through zero-copy column views.
   */
  void MultiplyBlock(Scalar a, const VectorBlock& x, Scalar b, VectorBlock& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return cached Hermiticity of the source matrix.
   *
   * @return `true` if the source matrix is Hermitian within `1e-12`.
   */
  bool IsHermitian() const override;

  /**
   * @brief Return the executor name (`"omp"` or `"cuda"`).
   */
  const std::string& ExecutorName() const { return executorName_; }

  /**
   * @brief Return the matrix format name.
   */
  const std::string& FormatName() const { return formatName_; }

  /**
   * @brief Return the Ginkgo executor backing this operator.
   *
   * @return Shared pointer to the `omp` or `cuda` executor.
   */
  std::shared_ptr<gko::Executor> GinkgoExecutor() const;

  /**
   * @brief Apply `y = a * O * x + b * y` directly on executor-resident dense blocks.
   *
   * Zero-copy apply on `gko::matrix::Dense` blocks that already live on this operator's
   * executor (device-resident engine path). Dimensions must match `Dimension() x b`.
   */
  void ApplyDevice(Scalar a, const gko::matrix::Dense<Scalar>* x, Scalar b,
                   gko::matrix::Dense<Scalar>* y) const;

private:
  struct Impl;
  std::shared_ptr<Impl> impl_;
  Index dimension_ = 0;
  bool isHermitian_ = false;
  std::string executorName_;
  std::string formatName_;
};

/**
 * @brief Report whether this build includes the Ginkgo backend.
 *
 * @return `true` when compiled with `LSQUANT_WITH_GINKGO=ON`.
 */
constexpr bool GinkgoAvailable() { return true; }

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_GINKGO_OPERATOR_HPP_
#endif  // LSQUANT_WITH_GINKGO
