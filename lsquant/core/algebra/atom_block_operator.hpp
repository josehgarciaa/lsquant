/**
 * @file atom_block_operator.hpp
 * @brief Atom-pair block sparse operator with a deduplicated dense-block dictionary.
 *
 * Stores a sparse operator whose rows and columns are grouped into atoms (contiguous orbital
 * ranges after a spatial permutation) as a neighbour list of atom pairs plus a dictionary of
 * distinct dense `W_a x W_b` blocks. On a crystal the dictionary degenerates to the handful of
 * distinct `H(R)` blocks (a stencil that stays in cache instead of streaming values); on an
 * amorphous system it holds one block per pair (a block-sparse matrix with spatial locality).
 * The class receives plain arrays and contains no physics. Python supplies the atom map and
 * permutation; AtomBlocksFromCsr builds neighbour lists, the dictionary and the diagonal.
 */
#ifndef LSQUANT_CORE_ALGEBRA_ATOM_BLOCK_OPERATOR_HPP_
#define LSQUANT_CORE_ALGEBRA_ATOM_BLOCK_OPERATOR_HPP_

#include <cstdint>
#include <vector>

#include "core/algebra/linear_operator.hpp"
#include "core/algebra/types.hpp"

namespace lsquant {

/** Plain atom-block arrays shared by construction and all execution modes. */
struct AtomBlockData {
  std::vector<Index> atomOffset, pairOffset, dictOffset;
  std::vector<std::int32_t> nbr, blockId, dictRows, dictCols;
  std::vector<Scalar> dictValues, diagonal;
};

/**
 * Build neighbour lists and deduplicate dense blocks by shape and exact bytes.
 * The permutation maps new orbitals to old ones and must keep each atom contiguous.
 * CSR duplicate entries are rejected; no Hermitian folding or physical convention is applied.
 * Work is parallel over atom rows, with private hash tables merged after construction.
 */
AtomBlockData AtomBlocksFromCsr(const Index* rowptr, const Index* columns,
                              const Scalar* values, Index dimension, Index nonzeros,
                              const Index* atomMap, const Index* permutation);

/**
 * @brief Matrix-free operator `y = a * O * x + b * y` over atom-pair dense blocks.
 *
 * Solves the problem of evaluating tight-binding matvecs without streaming one value per
 * non-zero from DRAM: each atom pair references a dictionary block by index, so on a crystal
 * only the pair indices (`8 B` per pair) and the vectors cross the memory system while the
 * block values stay resident in L1/L2. Atoms are contiguous orbital ranges of a permuted
 * Hilbert space (`atomOffset[a]` is the first orbital of atom `a`); for every atom the
 * neighbour list `nbr[pairOffset[a] .. pairOffset[a+1])` with block ids `blockId[...]` gives
 * the dense blocks to multiply against the neighbour's orbital slice. An optional on-site
 * diagonal (e.g. Anderson disorder) is applied separately so that disorder never multiplies
 * the dictionary. Thread parallelism partitions atoms into contiguous ranges balanced by pair
 * count (`RowBoundary` over `pairOffset`), each row-block reduction is owned by one thread and
 * summed in stored pair order, so results are bit-invariant across OpenMP thread counts.
 * `MultiplyBlock` performs, per column, exactly the arithmetic of `Multiply` (bit-identical
 * per state at any block width).
 *
 * @warning All arrays are validated at construction; the operator is immutable afterwards.
 * @par Pitfall
 *      The per-row accumulation order (pair order, then column order inside a block) differs
 *      from a column-sorted CSR row, so results match a CSR evaluation only up to floating-point
 *      reassociation (~1e-16 per row), not bit for bit.
 */
class AtomBlockOperator : public LinearOperator {
public:
  /**
   * @brief Construct the operator from plain atom-pair and dictionary arrays.
   *
   * @param atomOffset  Orbital offset per atom, length `A + 1`, non-decreasing, starting at 0;
   *                    `atomOffset[A]` is the Hilbert-space dimension `N`.
   * @param pairOffset  Pair offset per atom, length `A + 1`, non-decreasing, starting at 0;
   *                    `pairOffset[A]` is the pair count `P`.
   * @param nbr         Neighbour atom index per pair, length `P`, entries in `[0, A)`.
   * @param blockId     Dictionary block index per pair, length `P`, entries in `[0, D)`.
   * @param dictOffset  Value offset per dictionary block, length `D + 1`, non-decreasing,
   *                    starting at 0; block `d` occupies
   *                    `dictValues[dictOffset[d] .. dictOffset[d+1])` in row-major order with
   *                    `dictRows[d] * dictCols[d]` entries.
   * @param dictRows    Block row count per dictionary entry; must equal `W_a` of every atom
   *                    whose pair references it.
   * @param dictCols    Block column count per dictionary entry; must equal `W_b` of the
   *                    referenced neighbour atom.
   * @param dictValues  Contiguous complex values of all dictionary blocks.
   * @param diagonal    Optional on-site diagonal of length `N` (empty to disable); applied as
   *                    `acc[i] += diagonal[row] * x[row]` after the pair sum of each atom.
   * @param hermitian   Hermiticity flag reported by `IsHermitian()` (the caller knows the
   *                    physics; checking blockwise would cost `O(P * W^2)`).
   *
   * @warning Throws `std::invalid_argument` on any size, monotonicity, bounds, or block-shape
   *          violation.
   * @par Pitfall
   *      Neighbour atom slices are read through `atomOffset[nbr[p]]`; the permutation that makes
   *      atoms contiguous must be applied to vectors by the caller (`PermutedStates`).
   */
  AtomBlockOperator(std::vector<Index> atomOffset, std::vector<Index> pairOffset,
                    std::vector<std::int32_t> nbr, std::vector<std::int32_t> blockId,
                    std::vector<Index> dictOffset, std::vector<std::int32_t> dictRows,
                    std::vector<std::int32_t> dictCols, std::vector<Scalar> dictValues,
                    std::vector<Scalar> diagonal, bool hermitian);

  /**
   * @brief Return the Hilbert-space dimension `N` (total orbital count).
   *
   * @return Dimension `N >= 0`.
   */
  Index Dimension() const override;

  /**
   * @brief Compute `y = a * O * x + b * y`, row-blocks partitioned by pair count.
   *
   * For each atom `a` accumulates `acc[i] = sum_p sum_j D_p[i,j] * x[nbr_p slice j]` in stored
   * pair order, adds the optional diagonal term, then writes
   * `y[rowBegin + i] = a * acc[i] + b * y[rowBegin + i]`. Complexity is
   * `O(sum_p W_a * W_b)` work with `O(A)` parallel tasks.
   *
   * @param a Multiplicative scalar for `O * x`.
   * @param x Input vector of length `Dimension()`.
   * @param b Multiplicative scalar for `y`.
   * @param y Input/output vector of length `Dimension()`.
   *
   * @warning Throws `std::invalid_argument` on dimension mismatch or `&x == &y`.
   */
  void Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const override;

  /**
   * @brief Compute `Y = a * O * X + b * Y` on a contiguous column-major `N x S` block.
   *
   * Neighbour slices `x_b` (`W_b x S`, strided by `N` between columns) are staged into a
   * per-thread contiguous scratch and accumulated as small dense products with the `S`
   * direction innermost (SIMD-friendly); per `(row, column)` lane the accumulation order over
   * pairs and block columns is exactly the one of `Multiply`, so every output column is
   * bit-identical to a separate `Multiply` on that column. A single-column block dispatches to
   * the `Multiply` kernel directly.
   *
   * @param a Multiplicative scalar for `O * X`.
   * @param x Input block of size `Dimension() x S`.
   * @param b Multiplicative scalar for `Y`.
   * @param y Input/output block of size `Dimension() x S`.
   *
   * @warning Throws `std::invalid_argument` on dimension mismatch or aliased blocks.
   */
  void MultiplyBlock(Scalar a, const VectorBlock& x, Scalar b, VectorBlock& y) const override;

  using LinearOperator::Multiply;

  /**
   * @brief Compute `y = a * O * x + b * y` restricted to the atom rows `[atomBegin, atomEnd)`.
   *
   * Raw-pointer variant used by the MPI wrapper: interior and boundary atom ranges are
   * evaluated in separate phases around the halo exchange. `x` must cover every neighbour
   * slice the range references and `y` the full local row space; no size checks are made.
   * Rows are partitioned by pair count within the range, each row's reduction keeps the
   * stored pair order, so a `[0, interiorEnd) + [interiorEnd, A)` two-phase evaluation is
   * bit-identical to a single full `Multiply` call.
   *
   * @param a         Multiplicative scalar for `O * x`.
   * @param x         Input vector storage (full local row space).
   * @param b         Multiplicative scalar for `y`.
   * @param y         Input/output vector storage (full local row space).
   * @param atomBegin First atom row-block `0 <= atomBegin <= atomEnd`.
   * @param atomEnd   One past the last atom row-block `atomEnd <= NumAtoms()`.
   */
  void MultiplyRange(Scalar a, const Scalar* x, Scalar b, Scalar* y, Index atomBegin,
                     Index atomEnd) const;

  /**
   * @brief Compute `Y = a * O * X + b * Y` on atom rows `[atomBegin, atomEnd)` for `S` columns.
   *
   * Column-major blocks with leading dimension `Dimension()`; per column the arithmetic of
   * `MultiplyRange`. Same contract and use as `MultiplyRange`.
   *
   * @param a         Multiplicative scalar for `O * X`.
   * @param x         Column-major input block storage, leading dimension `Dimension()`.
   * @param b         Multiplicative scalar for `Y`.
   * @param y         Column-major input/output block storage, leading dimension `Dimension()`.
   * @param cols      Number of columns `S > 0`.
   * @param atomBegin First atom row-block.
   * @param atomEnd   One past the last atom row-block.
   */
  void MultiplyBlockRange(Scalar a, const Scalar* x, Scalar b, Scalar* y, int cols,
                          Index atomBegin, Index atomEnd) const;

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
  // Shared single-vector kernel on raw pointers; Multiply and the single-column
  // MultiplyBlock path both call it so their arithmetic is bit-identical.
  void MultiplyKernel(Scalar a, const Scalar* xp, Scalar b, Scalar* yp) const;
#ifdef LSQUANT_NATIVE_BUILD
  void PlaceOnThreads();
#endif
  std::vector<Index> atomOffset_;
  std::vector<Index> pairOffset_;
  std::vector<std::int32_t> nbr_;
  std::vector<std::int32_t> blockId_;
  std::vector<Index> dictOffset_;
  std::vector<std::int32_t> dictRows_;
  std::vector<std::int32_t> dictCols_;
  std::vector<Scalar> dictValues_;
  std::vector<Scalar> diagonal_;
  Index dimension_ = 0;
  Index maxBlockRows_ = 0;   // max W_a over atoms (scratch sizing)
  Index maxBlockCols_ = 0;   // max W_b over dictionary (scratch sizing)
  bool hermitian_ = false;
};

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_ATOM_BLOCK_OPERATOR_HPP_
