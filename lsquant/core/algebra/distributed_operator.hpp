/**
 * @file distributed_operator.hpp
 * @brief MPI wrapper over a local atom-block operator: halo exchange, never a third kernel.
 *
 * `DistributedOperator` presents a `LinearOperator` of the LOCAL row count to the engine: all
 * recursion vectors are distributed by rows across ranks. Each `MultiplyBlock` (a) starts a
 * non-blocking neighbour exchange (`MPI_Ineighbor_alltoallv` on a distributed-graph
 * communicator) of the boundary rows of X for all columns at once, (b) multiplies the
 * interior atoms — a contiguous prefix of the local atom range, guaranteed by the spatial
 * permutation — while the exchange is in flight, (c) waits, then (d) multiplies the boundary
 * atoms out of the received halo buffer. Every output row is accumulated by exactly one rank
 * in the same pair order as the OpenMP kernel, so the result is bit-identical to the OpenMP
 * atom-block operator row by row. No matrix data is duplicated across ranks; only the halo
 * (boundary atoms' vector slices) is exchanged per call.
 */
#ifndef LSQUANT_CORE_ALGEBRA_DISTRIBUTED_OPERATOR_HPP_
#define LSQUANT_CORE_ALGEBRA_DISTRIBUTED_OPERATOR_HPP_

#ifdef LSQUANT_WITH_MPI

#include <mpi.h>

#include <vector>

#include "core/algebra/atom_block_operator.hpp"
#include "core/algebra/linear_operator.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Halo description of one MPI rank: which local atoms are sent and which halo slots
 *        receive remote atoms.
 *
 * The local atom-block operator's neighbour indices are either local atoms `[0, localAtoms)`
 * or halo slots `[localAtoms, localAtoms + haloAtoms)`. For each neighbour rank (in the order
 * of the distributed-graph communicator), `sendOffset`/`sendAtoms` list the local atoms whose
 * vector slices this rank sends, and `recvOffset`/`recvAtoms` the halo atoms (in halo-slot
 * order) received from that rank.
 */
struct AtomBlockHalo {
  Index localAtoms = 0;               ///< atoms owned by this rank
  Index interiorAtoms = 0;            ///< prefix of local atoms with no remote neighbour
  std::vector<Index> atomOffset;      ///< local + halo atom orbital offsets (halo appended)
  std::vector<int> neighborRanks;     ///< peer ranks in communicator order
  std::vector<Index> sendOffset;      ///< per peer: offsets into sendAtoms (size peers + 1)
  std::vector<std::int32_t> sendAtoms;///< local atom indices sent to each peer
  std::vector<Index> recvOffset;      ///< per peer: offsets into recvAtoms (size peers + 1)
  std::vector<std::int32_t> recvAtoms;///< halo atom slots (localAtoms + k) filled per peer
};

/**
 * @brief `LinearOperator` over a local `AtomBlockOperator` with non-blocking halo exchange.
 *
 * Solves the problem of running the moment engine distributed by rows without a third matvec
 * kernel: the local operator (OpenMP, or the CUDA one behind the same interface in `mpi+cuda`
 * mode) does all arithmetic; this wrapper only moves boundary vector slices. `Dimension()` is
 * the local row count; `x`/`y` in every call are local vectors of that length. Interior rows
 * are computed during the exchange, boundary rows after it, and each output row is written by
 * exactly one rank, so global assembly needs no reduction inside the operator.
 *
 * @warning The halo exchange is collective across the ranks of the communicator: every rank
 *          must call `Multiply`/`MultiplyBlock` in the same sequence (the engine's plan loop
 *          guarantees this).
 * @par Pitfall
 *      The engine's dot tables are local partials; the global sum is the engine's reduction
 *      interface (`MomentReduction`), not this class.
 */
class DistributedOperator : public LinearOperator {
public:
  /**
   * @brief Construct the distributed operator over a local operator and its halo description.
   *
   * @param local  Local atom-block operator (OpenMP `AtomBlockOperator`); its atom neighbour
   *               indices may reference halo slots `>= halo.localAtoms`, resolved from the
   *               receive buffer after the exchange. Must outlive `*this` if non-owning is
   *               ever supported; this constructor takes ownership.
   * @param halo   Halo description (send/recv atom lists per peer rank).
   * @param comm   MPI communicator over the participating ranks (duplicated internally).
   *
   * @warning Throws `std::invalid_argument` on inconsistent halo arrays or MPI failures
   *          (as `std::runtime_error`).
   */
  DistributedOperator(std::unique_ptr<AtomBlockOperator> local, AtomBlockHalo halo,
                      MPI_Comm comm);

  /**
   * @brief Free the distributed-graph communicator and the duplicated communicator.
   */
  ~DistributedOperator() override;

  /**
   * @brief Return the LOCAL row count of this rank.
   *
   * @return Local dimension `N_local >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Compute the local rows of `y = a * O * x + b * y` with halo exchange.
   *
   * @param a Multiplicative scalar for `O * x`.
   * @param x Local input vector of length `Dimension()`.
   * @param b Multiplicative scalar for `y`.
   * @param y Local input/output vector of length `Dimension()`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  /**
   * @brief Compute the local rows of `Y = a * O * X + b * Y`, one exchange for all columns.
   *
   * @param a Multiplicative scalar for `O * X`.
   * @param x Local input block of size `Dimension() x S`.
   * @param b Multiplicative scalar for `Y`.
   * @param y Local input/output block of size `Dimension() x S`.
   */
  void MultiplyBlock(Scalar a, const VectorBlock& x, Scalar b, VectorBlock& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Return the local operator's Hermiticity flag.
   *
   * @return `true` if the local operator is Hermitian.
   */
  bool IsHermitian() const override;

private:
  std::unique_ptr<AtomBlockOperator> local_;
  AtomBlockHalo halo_;
  MPI_Comm comm_ = MPI_COMM_NULL;      ///< duplicated base communicator
  MPI_Comm graphComm_ = MPI_COMM_NULL; ///< distributed-graph communicator over the peers
  mutable std::vector<Scalar> sendBuf_;///< scratch: boundary slices of X per peer
  mutable std::vector<Scalar> recvBuf_;///< scratch: received halo slices of X
  mutable int maxCols_ = 0;            ///< scratch sizing watermark
};

}  // namespace lsquant

#endif  // LSQUANT_WITH_MPI
#endif  // LSQUANT_CORE_ALGEBRA_DISTRIBUTED_OPERATOR_HPP_
