/**
 * @file atom_block_cuda.hpp
 * @brief CUDA execution mode of the atom-block operator on the shared SoA data layout.
 *
 * Owns device copies of exactly the arrays that define `AtomBlockOperator` (atom offsets,
 * neighbour lists, block ids, the block dictionary, the on-site diagonal) uploaded verbatim
 * with `cudaMemcpy`, and evaluates `y = a * O * x + b * y` with a hand-written warp-per-atom
 * kernel: the dictionary is staged in dynamic shared memory when it fits, global memory
 * otherwise (where it stays L2-resident on crystals), and neighbour slices are staged through
 * shared memory with the block columns innermost. Per output element the accumulation order
 * (pair order, then block-column order) matches the OpenMP kernel, so CUDA and OpenMP results
 * differ only by device-level floating-point contraction (~1e-16).
 */
#ifndef LSQUANT_CORE_ALGEBRA_ATOM_BLOCK_CUDA_HPP_
#define LSQUANT_CORE_ALGEBRA_ATOM_BLOCK_CUDA_HPP_

#ifdef LSQUANT_WITH_CUDA

#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/** Wait for the current device's queued work; throws on asynchronous errors. */
void SynchronizeCuda();

/**
 * @brief Device-resident column-major `N x S` block of complex state vectors.
 *
 * Owns a `cudaMalloc`'d buffer with the same column-major layout as `VectorBlock` so the CUDA
 * matvec benchmark runs without per-call host transfers. Copy in/out is explicit
 * (`CopyFromHost` / `CopyToHost`).
 */
class CudaVectorBlock {
public:
  /**
   * @brief Allocate a zero-initialised device block of `rows x cols` scalars.
   *
   * @param rows Vector length `N >= 0`.
   * @param cols Number of vectors `S >= 0`.
   *
   * @warning Throws `std::runtime_error` on CUDA allocation or memset failure.
   */
  CudaVectorBlock(Index rows, int cols);

  /**
   * @brief Free the device buffer.
   */
  ~CudaVectorBlock();

  CudaVectorBlock(const CudaVectorBlock&) = delete;
  CudaVectorBlock& operator=(const CudaVectorBlock&) = delete;

  /**
   * @brief Return the vector length `N`.
   *
   * @return Row count `N >= 0`.
   */
  Index Rows() const { return rows_; }

  /**
   * @brief Return the number of vectors `S`.
   *
   * @return Column count `S >= 0`.
   */
  int Cols() const { return cols_; }

  /**
   * @brief Mutable device pointer to the contiguous column-major storage.
   *
   * @return Device pointer to `Rows() * Cols()` scalars.
   */
  Scalar* DeviceData() { return dev_; }

  /**
   * @brief Const device pointer to the contiguous column-major storage.
   *
   * @return Device pointer to `Rows() * Cols()` scalars.
   */
  const Scalar* DeviceData() const { return dev_; }

  /**
   * @brief Upload `Rows() * Cols()` scalars from a host buffer (column-major).
   *
   * @param host Host pointer to at least `Rows() * Cols()` scalars.
   */
  void CopyFromHost(const Scalar* host);

  /**
   * @brief Download `Rows() * Cols()` scalars into a host buffer (column-major).
   *
   * @param host Host pointer to at least `Rows() * Cols()` scalars.
   */
  void CopyToHost(Scalar* host) const;

private:
  Index rows_ = 0;
  int cols_ = 0;
  Scalar* dev_ = nullptr;
};

/**
 * @brief CUDA atom-block operator on the shared SoA layout (`y = a * O * x + b * y`).
 *
 * Construction validates and uploads the same arrays `AtomBlockOperator` consumes; no
 * reformatting or repacking happens on the device. `Multiply` / `MultiplyBlock` satisfy the
 * host `LinearOperator` interface by staging the vectors through host memory (correct for
 * engine correctness runs, deliberately not fast); `MultiplyBlockDevice` runs the kernel on
 * device-resident blocks without transfers (benchmarks, later the device engine).
 *
 * @warning Requires one visible CUDA device; construction throws `std::runtime_error` on any
 *          CUDA failure.
 * @par Pitfall
 *      State columns are chunked to keep `W_a * S_chunk <= 256` and shared scratch
 *      within 48 KiB. Individual atom widths above 256 are rejected by this kernel.
 */
class CudaAtomBlockOperator : public LinearOperator {
public:
  /**
   * @brief Validate and upload the atom-block arrays to the device.
   *
   * Parameters and validation are exactly those of `AtomBlockOperator`.
   *
   * @warning Throws `std::invalid_argument` on array violations, `std::runtime_error` on CUDA
   *          failures.
   */
  CudaAtomBlockOperator(std::vector<Index> atomOffset, std::vector<Index> pairOffset,
                        std::vector<std::int32_t> nbr, std::vector<std::int32_t> blockId,
                        std::vector<Index> dictOffset, std::vector<std::int32_t> dictRows,
                        std::vector<std::int32_t> dictCols, std::vector<Scalar> dictValues,
                        std::vector<Scalar> diagonal, bool hermitian);

  /**
   * @brief Free the device arrays.
   */
  ~CudaAtomBlockOperator() override;

  CudaAtomBlockOperator(const CudaAtomBlockOperator&) = delete;
  CudaAtomBlockOperator& operator=(const CudaAtomBlockOperator&) = delete;

  /**
   * @brief Return the Hilbert-space dimension `N`.
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Compute `y = a * O * x + b * y` via host staging (correct, not fast).
   *
   * @param a Multiplicative scalar for `O * x`.
   * @param x Host input vector of length `Dimension()`.
   * @param b Multiplicative scalar for `y`.
   * @param y Host input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  /**
   * @brief Compute `Y = a * O * X + b * Y` via host staging (correct, not fast).
   *
   * @param a Multiplicative scalar for `O * X`.
   * @param x Host input block of size `Dimension() x S`.
   * @param b Multiplicative scalar for `Y`.
   * @param y Host input/output block of size `Dimension() x S`.
   */
  void MultiplyBlock(Scalar a, const VectorBlock& x, Scalar b, VectorBlock& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Compute `Y = a * O * X + b * Y` on device-resident blocks without transfers.
   *
   * @param a Multiplicative scalar for `O * X`.
   * @param x Device input block of size `Dimension() x S`.
   * @param b Multiplicative scalar for `Y`.
   * @param y Device input/output block of size `Dimension() x S`.
   */
  void MultiplyBlockDevice(Scalar a, const CudaVectorBlock& x, Scalar b,
                           CudaVectorBlock& y) const;

  /**
   * @brief Return the Hermiticity flag supplied at construction.
   *
   * @return `true` if the caller declared the operator Hermitian.
   */
  bool IsHermitian() const override;

  /**
   * @brief Return the number of atoms `A`.
   *
   * @return Atom count `A >= 0`.
   */
  Index NumAtoms() const;

  /**
   * @brief Return the number of atom pairs `P`.
   *
   * @return Pair count `P >= 0`.
   */
  Index NumPairs() const;

  /**
   * @brief Return the number of distinct dictionary blocks `D`.
   *
   * @return Dictionary size `D >= 0`.
   */
  Index DictionarySize() const;

private:
  Index dimension_ = 0;
  Index numAtoms_ = 0;
  Index numPairs_ = 0;
  Index numBlocks_ = 0;
  Index maxBlockRows_ = 0;
  Index maxBlockCols_ = 0;
  bool hermitian_ = false;
  Index* atomOffset_ = nullptr;
  Index* pairOffset_ = nullptr;
  std::int32_t* nbr_ = nullptr;
  std::int32_t* blockId_ = nullptr;
  Index* dictOffset_ = nullptr;
  std::int32_t* dictRows_ = nullptr;
  std::int32_t* dictCols_ = nullptr;
  Scalar* dictValues_ = nullptr;
  Scalar* diagonal_ = nullptr;
  Index dictScalars_ = 0;
};

}  // namespace lsquant

#endif  // LSQUANT_WITH_CUDA
#endif  // LSQUANT_CORE_ALGEBRA_ATOM_BLOCK_CUDA_HPP_
