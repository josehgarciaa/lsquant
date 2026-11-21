/**
 * @file vector_block.hpp
 * @brief Contiguous column-major `N x b` block of complex state vectors.
 *
 * Stores `b` vectors of length `N` in one allocation with column-major layout (column `j`
 * occupies `data[j*N .. (j+1)*N)`), so that batched operator application (`MultiplyBlock`)
 * and block inner-product tables can stream one contiguous buffer and wrap it in Eigen or
 * BLAS-style views without per-column copies.
 */
#ifndef LSQUANT_CORE_ALGEBRA_VECTOR_BLOCK_HPP_
#define LSQUANT_CORE_ALGEBRA_VECTOR_BLOCK_HPP_

#include <algorithm>
#include <stdexcept>

#include <Eigen/Core>

#include "core/algebra/memory_placement.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Column-major dense block of `Cols()` state vectors of length `Rows()`.
 *
 * Solves the problem of holding a batch of Chebyshev or trace vectors in one contiguous
 * allocation so that sparse kernels stream the matrix once per row for all columns and the
 * moment table maps the columns as an Eigen matrix without copying. Column `j` starts at
 * `Data() + j * Rows()`, matching the layout of `Eigen::Matrix<Scalar, Dynamic, Dynamic>`
 * (column-major), so `Map()` views are zero-copy. Memory footprint is `16 * Rows() * Cols()`
 * bytes.
 *
 * @warning Views (`Column`, `Map`) are invalidated by any operation that reallocates the
 *          underlying storage (construction, `Swap` with a differently sized block).
 * @par Pitfall
 *      Columns are stored contiguously; algorithms that need row-major access must go
 *      through `Map()` (Eigen) rather than reinterpreting the raw buffer.
 */
class VectorBlock {
public:
  using DenseMap = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>;
  using ConstDenseMap =
      Eigen::Map<const Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>>;

  /**
   * @brief Construct an empty `0 x 0` block.
   */
  VectorBlock() = default;

  /**
   * @brief Allocate a zero-initialised `rows x cols` block.
   *
   * @param rows Vector length `N >= 0`.
   * @param cols Number of vectors `b >= 0`.
   *
   * @warning Throws `std::invalid_argument` if `rows < 0` or `cols < 0`.
   * @par Pitfall
   *      In native builds the pages are discarded and zeroed with a parallel first touch so
   *      each OpenMP thread faults its slice on its local NUMA node.
   */
  VectorBlock(Index rows, int cols) : rows_(rows), cols_(cols) {
    if (rows < 0 || cols < 0) {
      throw std::invalid_argument("VectorBlock: negative rows or cols");
    }
#ifdef LSQUANT_NATIVE_BUILD
    data_ = MakeZeroVector(rows * static_cast<Index>(cols));
#else
    data_.assign(static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols),
                 Scalar(0.0, 0.0));
#endif
  }

  /**
   * @brief Wrap an external column-major `rows x cols` buffer without taking ownership.
   *
   * Used by bindings and benchmark harnesses to apply operators directly on caller-owned
   * storage (e.g. a NumPy buffer) without a copy. The caller must guarantee the buffer
   * outlives this view and holds at least `rows * cols` scalars.
   *
   * @param rows Vector length `N >= 0`.
   * @param cols Number of vectors `b >= 0`.
   * @param data External column-major buffer (may be nullptr only when `rows * cols == 0`).
   */
  VectorBlock(Index rows, int cols, Scalar* data) : rows_(rows), cols_(cols), external_(data) {}

  /**
   * @brief Return the vector length `N`.
   *
   * @return Number of rows `N >= 0`.
   */
  Index Rows() const { return rows_; }

  /**
   * @brief Return the number of vectors `b` in the block.
   *
   * @return Number of columns `b >= 0`.
   */
  int Cols() const { return cols_; }

  /**
   * @brief Mutable access to the contiguous column-major storage.
   *
   * @return Pointer to the first of `Rows() * Cols()` scalars (nullptr when empty).
   */
  Scalar* Data() { return external_ != nullptr ? external_ : data_.data(); }

  /**
   * @brief Const access to the contiguous column-major storage.
   *
   * @return Pointer to the first of `Rows() * Cols()` scalars (nullptr when empty).
   */
  const Scalar* Data() const { return external_ != nullptr ? external_ : data_.data(); }

  /**
   * @brief Mutable access to column `j` (a contiguous vector of length `Rows()`).
   *
   * @param j Column index `0 <= j < Cols()`.
   * @return  Pointer to the first element of column `j`.
   *
   * @warning Throws `std::out_of_range` if `j < 0` or `j >= Cols()`.
   */
  Scalar* Column(int j) {
    if (j < 0 || j >= cols_) {
      throw std::out_of_range("VectorBlock::Column: index out of range");
    }
    return Data() + static_cast<std::size_t>(j) * static_cast<std::size_t>(rows_);
  }

  /**
   * @brief Const access to column `j` (a contiguous vector of length `Rows()`).
   *
   * @param j Column index `0 <= j < Cols()`.
   * @return  Const pointer to the first element of column `j`.
   *
   * @warning Throws `std::out_of_range` if `j < 0` or `j >= Cols()`.
   */
  const Scalar* Column(int j) const {
    if (j < 0 || j >= cols_) {
      throw std::out_of_range("VectorBlock::Column: index out of range");
    }
    return Data() + static_cast<std::size_t>(j) * static_cast<std::size_t>(rows_);
  }

  /**
   * @brief Mutable zero-copy Eigen view of the whole block as a column-major `N x b` matrix.
   *
   * @return `Eigen::Map` over the block storage.
   */
  DenseMap Map() {
    return DenseMap(Data(), static_cast<Eigen::Index>(rows_),
                    static_cast<Eigen::Index>(cols_));
  }

  /**
   * @brief Const zero-copy Eigen view of the whole block as a column-major `N x b` matrix.
   *
   * @return Const `Eigen::Map` over the block storage.
   */
  ConstDenseMap Map() const {
    return ConstDenseMap(Data(), static_cast<Eigen::Index>(rows_),
                         static_cast<Eigen::Index>(cols_));
  }

  /**
   * @brief Exchange the contents with `other` in `O(1)`.
   *
   * @param other Block to swap with (may have different dimensions).
   */
  void Swap(VectorBlock& other) {
    std::swap(rows_, other.rows_);
    std::swap(cols_, other.cols_);
    std::swap(external_, other.external_);
    data_.swap(other.data_);
  }

private:
  Index rows_ = 0;
  int cols_ = 0;
  Scalar* external_ = nullptr;
  Vector data_;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_VECTOR_BLOCK_HPP_
