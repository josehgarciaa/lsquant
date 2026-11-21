/**
 * @file atom_block_validate.hpp
 * @brief Shared validation of the atom-block SoA arrays for the OpenMP and CUDA operators.
 *
 * Both `AtomBlockOperator` (OpenMP) and `CudaAtomBlockOperator` (CUDA) consume the same
 * plain arrays; this header keeps their construction-time validation and the derived
 * quantities (dimension, largest block row/column width) identical and in one place.
 */
#ifndef LSQUANT_CORE_ALGEBRA_ATOM_BLOCK_VALIDATE_HPP_
#define LSQUANT_CORE_ALGEBRA_ATOM_BLOCK_VALIDATE_HPP_

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/algebra/types.hpp"

namespace lsquant {

/**
 * @brief Derived, validated extents of an atom-block array set.
 */
struct AtomBlockExtents {
  Index dimension = 0;     ///< Hilbert-space dimension `N` (total orbitals)
  Index numAtoms = 0;      ///< Atom count `A`
  Index numPairs = 0;      ///< Pair count `P`
  Index numBlocks = 0;     ///< Dictionary size `D`
  Index maxBlockRows = 0;  ///< Largest atom width `W_a`
  Index maxBlockCols = 0;  ///< Largest dictionary block column count `W_b`
};

/**
 * @brief Validate the atom-block arrays and return their extents.
 *
 * Checks offset monotonicity and endpoints, index bounds, dictionary block shapes against
 * every referencing pair, and the diagonal length. Complexity is `O(A + P + D)`.
 *
 * @warning Throws `std::invalid_argument` on any violation.
 */
inline AtomBlockExtents ValidateAtomBlockArrays(
    const std::vector<Index>& atomOffset, const std::vector<Index>& pairOffset,
    const std::vector<std::int32_t>& nbr, const std::vector<std::int32_t>& blockId,
    const std::vector<Index>& dictOffset, const std::vector<std::int32_t>& dictRows,
    const std::vector<std::int32_t>& dictCols, const std::vector<Scalar>& dictValues,
    const std::vector<Scalar>& diagonal) {
  auto checkOffsets = [](const char* name, const std::vector<Index>& offsets,
                         Index expectedEnd) {
    if (offsets.empty() || offsets.front() != 0) {
      throw std::invalid_argument(std::string("AtomBlockOperator: ") + name +
                                  " must start at 0");
    }
    for (std::size_t i = 1; i < offsets.size(); ++i) {
      if (offsets[i] < offsets[i - 1]) {
        throw std::invalid_argument(std::string("AtomBlockOperator: ") + name +
                                    " must be non-decreasing");
      }
    }
    if (offsets.back() != expectedEnd) {
      throw std::invalid_argument(std::string("AtomBlockOperator: ") + name +
                                  " endpoint does not match the array lengths");
    }
  };

  AtomBlockExtents ext;
  if (atomOffset.size() < 2) {
    throw std::invalid_argument("AtomBlockOperator: at least one atom is required");
  }
  ext.numAtoms = static_cast<Index>(atomOffset.size()) - 1;
  ext.dimension = atomOffset.back();
  if (pairOffset.size() != atomOffset.size()) {
    throw std::invalid_argument("AtomBlockOperator: pairOffset must have one entry per atom");
  }
  if (nbr.size() != blockId.size()) {
    throw std::invalid_argument("AtomBlockOperator: nbr and blockId size mismatch");
  }
  checkOffsets("atomOffset", atomOffset, ext.dimension);
  checkOffsets("pairOffset", pairOffset, static_cast<Index>(nbr.size()));
  if (dictOffset.size() != dictRows.size() + 1 || dictRows.size() != dictCols.size()) {
    throw std::invalid_argument("AtomBlockOperator: dictionary array size mismatch");
  }
  checkOffsets("dictOffset", dictOffset, static_cast<Index>(dictValues.size()));
  if (!diagonal.empty() && static_cast<Index>(diagonal.size()) != ext.dimension) {
    throw std::invalid_argument("AtomBlockOperator: diagonal must be empty or have length N");
  }
  ext.numPairs = static_cast<Index>(nbr.size());
  ext.numBlocks = static_cast<Index>(dictRows.size());
  for (Index d = 0; d < ext.numBlocks; ++d) {
    if (dictRows[d] <= 0 || dictCols[d] <= 0) {
      throw std::invalid_argument("AtomBlockOperator: dictionary blocks must be non-empty");
    }
    const Index expect = static_cast<Index>(dictRows[d]) * dictCols[d];
    if (dictOffset[d + 1] - dictOffset[d] != expect) {
      throw std::invalid_argument("AtomBlockOperator: dictionary block size mismatch");
    }
    ext.maxBlockCols = std::max(ext.maxBlockCols, static_cast<Index>(dictCols[d]));
  }
  for (Index a = 0; a < ext.numAtoms; ++a) {
    const Index width = atomOffset[a + 1] - atomOffset[a];
    if (width <= 0) {
      throw std::invalid_argument("AtomBlockOperator: atoms must hold at least one orbital");
    }
    ext.maxBlockRows = std::max(ext.maxBlockRows, width);
  }
  for (Index a = 0; a < ext.numAtoms; ++a) {
    const Index width = atomOffset[a + 1] - atomOffset[a];
    for (Index p = pairOffset[a]; p < pairOffset[a + 1]; ++p) {
      const std::int32_t b = nbr[static_cast<std::size_t>(p)];
      const std::int32_t d = blockId[static_cast<std::size_t>(p)];
      if (b < 0 || b >= ext.numAtoms) {
        throw std::invalid_argument("AtomBlockOperator: neighbour atom index out of range");
      }
      if (d < 0 || d >= ext.numBlocks) {
        throw std::invalid_argument("AtomBlockOperator: dictionary index out of range");
      }
      if (dictRows[d] != width) {
        throw std::invalid_argument(
            "AtomBlockOperator: block row count does not match atom");
      }
      if (dictCols[d] != atomOffset[b + 1] - atomOffset[b]) {
        throw std::invalid_argument(
            "AtomBlockOperator: block column count does not match atom");
      }
    }
  }
  return ext;
}

}  // namespace lsquant

#endif  // LSQUANT_CORE_ALGEBRA_ATOM_BLOCK_VALIDATE_HPP_
