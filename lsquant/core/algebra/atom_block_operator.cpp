/**
 * @file atom_block_operator.cpp
 * @brief Implementation of AtomBlockOperator.
 */
#include "core/algebra/atom_block_operator.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <deque>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>

#include "core/algebra/atom_block_validate.hpp"
#include "core/algebra/memory_placement.hpp"
#include "core/algebra/work_partition.hpp"

namespace lsquant {

namespace {
constexpr std::size_t kDictionaryShards = 256;

// Finite operator/vector data use the elementary complex product. Keeping this explicit
// avoids the exceptional-value fallback in std::complex multiplication in the SIMD loop.
inline Scalar BlockProduct(Scalar left, Scalar right) {
  return Scalar(left.real() * right.real() - left.imag() * right.imag(),
                left.real() * right.imag() + left.imag() * right.real());
}

Index AtomRangeBoundary(const Index* offsets, Index begin, Index end, int part, int parts) {
  if (part == 0) return begin;
  if (part == parts) return end;
  const Index pairs = offsets[end] - offsets[begin];
  if (pairs == 0) return begin + (end - begin) * part / parts;
  const Index target = offsets[begin] + (pairs / parts) * part + (pairs % parts) * part / parts;
  return static_cast<Index>(std::lower_bound(offsets + begin, offsets + end, target) - offsets);
}
}  // namespace

AtomBlockData AtomBlocksFromCsr(const Index* rowptr, const Index* columns,
                              const Scalar* values, Index dimension, Index nonzeros,
                              const Index* atomMap, const Index* permutation) {
  if (dimension <= 0 || nonzeros < 0 || rowptr[0] != 0 || rowptr[dimension] != nonzeros) {
    throw std::invalid_argument("AtomBlocksFromCsr: invalid CSR extents");
  }
  std::vector<Index> inverse(dimension, -1), orbitalAtom(dimension), localIndex(dimension);
  for (Index i = 0; i < dimension; ++i) {
    if (rowptr[i] < 0 || rowptr[i] > rowptr[i + 1] || rowptr[i + 1] > nonzeros) {
      throw std::invalid_argument("AtomBlocksFromCsr: invalid row pointers");
    }
    Index previous = -1;
    for (Index p = rowptr[i]; p < rowptr[i + 1]; ++p) {
      if (columns[p] <= previous || columns[p] >= dimension) {
        throw std::invalid_argument("AtomBlocksFromCsr: columns must be sorted and unique");
      }
      previous = columns[p];
    }
    const Index old = permutation[i];
    if (old < 0 || old >= dimension || inverse[old] != -1) {
      throw std::invalid_argument("AtomBlocksFromCsr: invalid permutation");
    }
    inverse[old] = i;
  }
  AtomBlockData out;
  out.atomOffset.push_back(0);
  std::unordered_map<Index, Index> seenAtoms;
  Index previousAtom = -1;
  for (Index i = 0; i < dimension; ++i) {
    const Index atom = atomMap[permutation[i]];
    if (atom < 0) throw std::invalid_argument("AtomBlocksFromCsr: negative atom index");
    if (i == 0 || atom != previousAtom) {
      if (!seenAtoms.emplace(atom, static_cast<Index>(seenAtoms.size())).second) {
        throw std::invalid_argument("AtomBlocksFromCsr: permutation splits an atom");
      }
      if (i != 0) out.atomOffset.push_back(i);
      previousAtom = atom;
    }
    orbitalAtom[i] = static_cast<Index>(out.atomOffset.size()) - 1;
    localIndex[i] = i - out.atomOffset.back();
  }
  out.atomOffset.push_back(dimension);
  const Index atoms = static_cast<Index>(out.atomOffset.size()) - 1;
  if (atoms > std::numeric_limits<std::int32_t>::max()) {
    throw std::invalid_argument("AtomBlocksFromCsr: atom count exceeds int32");
  }
  out.pairOffset.resize(atoms + 1, 0);
  out.diagonal.resize(dimension, Scalar(0));
  // Only the small neighbour set of one atom is sorted, never an nnz-sized global key array.
#pragma omp parallel
  {
    std::vector<std::int32_t> neighbors;
#pragma omp for schedule(static)
    for (Index atom = 0; atom < atoms; ++atom) {
      neighbors.clear();
      for (Index i = out.atomOffset[atom]; i < out.atomOffset[atom + 1]; ++i) {
        const Index old = permutation[i];
        for (Index p = rowptr[old]; p < rowptr[old + 1]; ++p) {
          if (columns[p] != old) neighbors.push_back(static_cast<std::int32_t>(orbitalAtom[inverse[columns[p]]]));
        }
      }
      std::sort(neighbors.begin(), neighbors.end());
      out.pairOffset[atom + 1] = std::unique(neighbors.begin(), neighbors.end()) - neighbors.begin();
    }
  }
  for (Index a = 0; a < atoms; ++a) out.pairOffset[a + 1] += out.pairOffset[a];
  out.nbr.resize(out.pairOffset.back());
  out.blockId.resize(out.pairOffset.back());
  struct LocalDictionary {
    // Stable owned keys allow both hash tables to use views: no second/third key copy.
    std::unordered_map<std::string_view, std::int32_t> ids;
    std::deque<std::string> keys;
    std::vector<Index> firstPairs;
    std::array<std::vector<std::int32_t>, kDictionaryShards> shards;
    Index firstAtom = 0, lastAtom = 0;
  };
  std::vector<LocalDictionary> dictionaries(omp_get_max_threads());
#pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    const int threads = omp_get_num_threads();
    auto& dict = dictionaries[tid];
    dict.firstAtom = atoms * tid / threads;
    dict.lastAtom = atoms * (tid + 1) / threads;
    std::vector<std::int32_t> neighbors;
    std::vector<Index> offsets;
    std::vector<Scalar> blocks;
    std::string key;
    for (Index atom = dict.firstAtom; atom < dict.lastAtom; ++atom) {
      const Index begin = out.atomOffset[atom], end = out.atomOffset[atom + 1];
      const auto width = static_cast<std::int32_t>(end - begin);
      neighbors.clear();
      for (Index i = begin; i < end; ++i) {
        const Index old = permutation[i];
        for (Index p = rowptr[old]; p < rowptr[old + 1]; ++p) {
          if (columns[p] != old) neighbors.push_back(static_cast<std::int32_t>(orbitalAtom[inverse[columns[p]]]));
        }
      }
      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
      offsets.assign(neighbors.size() + 1, 0);
      for (std::size_t k = 0; k < neighbors.size(); ++k) {
        const Index nb = neighbors[k];
        offsets[k + 1] = offsets[k] + width * (out.atomOffset[nb + 1] - out.atomOffset[nb]);
      }
      blocks.assign(offsets.back(), Scalar(0));
      for (Index i = begin; i < end; ++i) {
        const Index old = permutation[i];
        for (Index p = rowptr[old]; p < rowptr[old + 1]; ++p) {
          if (columns[p] == old) { out.diagonal[i] = values[p]; continue; }
          const Index j = inverse[columns[p]], nb = orbitalAtom[j];
          const auto k = std::lower_bound(neighbors.begin(), neighbors.end(), nb) - neighbors.begin();
          const Index cols = out.atomOffset[nb + 1] - out.atomOffset[nb];
          blocks[offsets[k] + (i - begin) * cols + localIndex[j]] = values[p];
        }
      }
      for (std::size_t k = 0; k < neighbors.size(); ++k) {
        const Index nb = neighbors[k];
        const auto cols = static_cast<std::int32_t>(out.atomOffset[nb + 1] - out.atomOffset[nb]);
        key.resize(2 * sizeof(std::int32_t) + (offsets[k + 1] - offsets[k]) * sizeof(Scalar));
        std::memcpy(key.data(), &width, sizeof(width));
        std::memcpy(key.data() + sizeof(width), &cols, sizeof(cols));
        std::memcpy(key.data() + 2 * sizeof(width), blocks.data() + offsets[k], (offsets[k + 1] - offsets[k]) * sizeof(Scalar));
        auto entry = dict.ids.find(key);
        if (entry == dict.ids.end()) {
          const auto id = static_cast<std::int32_t>(dict.keys.size());
          dict.keys.push_back(key);
          dict.firstPairs.push_back(out.pairOffset[atom] + static_cast<Index>(k));
          entry = dict.ids.emplace(std::string_view(dict.keys.back()), id).first;
          const auto shard = std::hash<std::string_view>{}(dict.keys.back()) % kDictionaryShards;
          dict.shards[shard].push_back(id);
        }
        const Index p = out.pairOffset[atom] + static_cast<Index>(k);
        out.nbr[p] = static_cast<std::int32_t>(nb);
        out.blockId[p] = entry->second;
      }
    }
  }
  // An amorphous dictionary can contain millions of distinct blocks. Merge independent
  // hash shards in parallel instead of serialising every unique block through one table.
  std::vector<std::vector<Index>> remaps(dictionaries.size());
  for (std::size_t t = 0; t < dictionaries.size(); ++t) {
    remaps[t].resize(dictionaries[t].keys.size());
  }
  std::array<std::vector<std::string_view>, kDictionaryShards> shardKeys;
  std::array<std::vector<Index>, kDictionaryShards> shardFirstPairs;
  std::array<Index, kDictionaryShards + 1> dictionaryOffsets{}, valueOffsets{};
#pragma omp parallel for schedule(dynamic, 1)
  for (std::size_t shard = 0; shard < kDictionaryShards; ++shard) {
    std::size_t count = 0;
    for (const auto& dict : dictionaries) count += dict.shards[shard].size();
    std::unordered_map<std::string_view, Index> merged;
    merged.reserve(count);
    for (std::size_t t = 0; t < dictionaries.size(); ++t) {
      const auto& dict = dictionaries[t];
      for (const auto k : dict.shards[shard]) {
        const std::string_view key(dict.keys[k]);
        const auto inserted = merged.emplace(key, static_cast<Index>(merged.size()));
        remaps[t][k] = inserted.first->second;
        if (inserted.second) {
          shardKeys[shard].push_back(key);
          shardFirstPairs[shard].push_back(dict.firstPairs[k]);
          valueOffsets[shard + 1] += static_cast<Index>(
              (key.size() - 2 * sizeof(std::int32_t)) / sizeof(Scalar));
        }
      }
    }
    dictionaryOffsets[shard + 1] = static_cast<Index>(merged.size());
  }
  for (std::size_t shard = 0; shard < kDictionaryShards; ++shard) {
    dictionaryOffsets[shard + 1] += dictionaryOffsets[shard];
    valueOffsets[shard + 1] += valueOffsets[shard];
  }
  if (dictionaryOffsets.back() > std::numeric_limits<std::int32_t>::max()) {
    throw std::overflow_error("AtomBlocksFromCsr: dictionary exceeds int32");
  }
  out.dictRows.resize(dictionaryOffsets.back());
  out.dictCols.resize(dictionaryOffsets.back());
  out.dictOffset.resize(dictionaryOffsets.back() + 1);
  out.dictValues.resize(valueOffsets.back());
  std::vector<std::string_view> mergedKeys(dictionaryOffsets.back());
  std::vector<std::pair<Index, Index>> order(dictionaryOffsets.back());
#pragma omp parallel for schedule(static)
  for (std::size_t shard = 0; shard < kDictionaryShards; ++shard) {
    for (std::size_t k = 0; k < shardKeys[shard].size(); ++k) {
      const Index d = dictionaryOffsets[shard] + static_cast<Index>(k);
      mergedKeys[d] = shardKeys[shard][k];
      order[d] = {shardFirstPairs[shard][k], d};
    }
  }
  // Keep values in first-use order. Hash-bucket order would scatter amorphous dictionary
  // reads across DRAM and discard the locality supplied by the spatial atom permutation.
  if (order.size() > out.nbr.size() / 8) {
    // First-use positions are unique pair indices. Dense dictionaries can be ordered by
    // a linear index scan; comparison sorting is cheaper for a small crystal dictionary.
    std::vector<std::int32_t> firstUse(out.nbr.size(), -1);
#pragma omp parallel for schedule(static)
    for (std::size_t d = 0; d < order.size(); ++d) {
      firstUse[order[d].first] = static_cast<std::int32_t>(order[d].second);
    }
    std::size_t d = 0;
    for (std::size_t p = 0; p < firstUse.size(); ++p) {
      if (firstUse[p] >= 0) order[d++] = {static_cast<Index>(p), firstUse[p]};
    }
  } else {
    std::sort(order.begin(), order.end());
  }
  std::vector<std::int32_t> newId(order.size());
#pragma omp parallel for schedule(static)
  for (std::size_t d = 0; d < order.size(); ++d) {
    const auto key = mergedKeys[order[d].second];
    newId[order[d].second] = static_cast<std::int32_t>(d);
    std::memcpy(&out.dictRows[d], key.data(), sizeof(std::int32_t));
    std::memcpy(&out.dictCols[d], key.data() + sizeof(std::int32_t), sizeof(std::int32_t));
    out.dictOffset[d + 1] = static_cast<Index>(out.dictRows[d]) * out.dictCols[d];
  }
  out.dictOffset[0] = 0;
  for (std::size_t d = 0; d < order.size(); ++d) out.dictOffset[d + 1] += out.dictOffset[d];
#pragma omp parallel for schedule(static)
  for (std::size_t d = 0; d < order.size(); ++d) {
    const auto key = mergedKeys[order[d].second];
    std::memcpy(out.dictValues.data() + out.dictOffset[d], key.data() + 2 * sizeof(std::int32_t),
                (out.dictOffset[d + 1] - out.dictOffset[d]) * sizeof(Scalar));
  }
#pragma omp parallel for schedule(static)
  for (std::size_t t = 0; t < dictionaries.size(); ++t) {
    auto& dict = dictionaries[t];
    auto& remap = remaps[t];
    for (std::size_t k = 0; k < dict.keys.size(); ++k) {
      const auto shard = std::hash<std::string_view>{}(dict.keys[k]) % kDictionaryShards;
      remap[k] += dictionaryOffsets[shard];
    }
    for (Index p = out.pairOffset[dict.firstAtom]; p < out.pairOffset[dict.lastAtom]; ++p) {
      out.blockId[p] = newId[remap[out.blockId[p]]];
    }
    // Release the many amorphous keys on their worker threads too. Serial destruction of
    // thread-private hash nodes would otherwise become the construction bottleneck again.
    dict.ids.clear();
    dict.ids.rehash(0);
    dict.keys.clear();
  }
  return out;
}

AtomBlockOperator::AtomBlockOperator(std::vector<Index> atomOffset,
                                     std::vector<Index> pairOffset, std::vector<std::int32_t> nbr,
                                     std::vector<std::int32_t> blockId,
                                     std::vector<Index> dictOffset,
                                     std::vector<std::int32_t> dictRows,
                                     std::vector<std::int32_t> dictCols,
                                     std::vector<Scalar> dictValues,
                                     std::vector<Scalar> diagonal, bool hermitian)
    : atomOffset_(std::move(atomOffset)),
      pairOffset_(std::move(pairOffset)),
      nbr_(std::move(nbr)),
      blockId_(std::move(blockId)),
      dictOffset_(std::move(dictOffset)),
      dictRows_(std::move(dictRows)),
      dictCols_(std::move(dictCols)),
      dictValues_(std::move(dictValues)),
      diagonal_(std::move(diagonal)),
      hermitian_(hermitian) {
  const AtomBlockExtents ext =
      ValidateAtomBlockArrays(atomOffset_, pairOffset_, nbr_, blockId_, dictOffset_, dictRows_,
                              dictCols_, dictValues_, diagonal_);
  dimension_ = ext.dimension;
  maxBlockRows_ = ext.maxBlockRows;
  maxBlockCols_ = ext.maxBlockCols;
#ifdef LSQUANT_NATIVE_BUILD
  PlaceOnThreads();
#endif
}

Index AtomBlockOperator::Dimension() const { return dimension_; }

Index AtomBlockOperator::NumAtoms() const {
  return static_cast<Index>(atomOffset_.size()) - 1;
}

Index AtomBlockOperator::NumPairs() const { return static_cast<Index>(nbr_.size()); }

Index AtomBlockOperator::DictionarySize() const {
  return static_cast<Index>(dictRows_.size());
}

bool AtomBlockOperator::IsHermitian() const { return hermitian_; }

void AtomBlockOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index n = Dimension();
  if (static_cast<Index>(x.size()) != n || static_cast<Index>(y.size()) != n) {
    throw std::invalid_argument("AtomBlockOperator::Multiply: vector dimension mismatch");
  }
  if (&x == &y) {
    throw std::invalid_argument(
        "AtomBlockOperator::Multiply: input and output vectors must not alias");
  }
  MultiplyKernel(a, x.data(), b, y.data());
}

void AtomBlockOperator::MultiplyKernel(Scalar a, const Scalar* xp, Scalar b,
                                       Scalar* yp) const {
  MultiplyRange(a, xp, b, yp, 0, NumAtoms());
}

void AtomBlockOperator::MultiplyRange(Scalar a, const Scalar* xp, Scalar b, Scalar* yp,
                                      Index atomBegin, Index atomEnd) const {
  const Index numAtoms = NumAtoms();
  if (atomBegin < 0 || atomEnd > numAtoms || atomBegin > atomEnd) {
    throw std::invalid_argument("AtomBlockOperator::MultiplyRange: invalid atom range");
  }
  if (Dimension() == 0 || numAtoms == 0 || atomBegin == atomEnd) {
    return;
  }

  const Index* const atomOff = atomOffset_.data();
  const Index* const pairOff = pairOffset_.data();
  const std::int32_t* const nbr = nbr_.data();
  const std::int32_t* const blk = blockId_.data();
  const Index* const dictOff = dictOffset_.data();
  const std::int32_t* const dictCols = dictCols_.data();
  const Scalar* const dictVal = dictValues_.data();
  const Scalar* const diag = diagonal_.empty() ? nullptr : diagonal_.data();
  const Index maxRows = maxBlockRows_;

  if (maxRows == 1 && dictValues_.size() == dictRows_.size()) {
    // Every atom and dictionary block is scalar: offsets are exactly the identity.
    // Avoid three metadata gathers per pair (atom offset, dictionary offset and width).
#pragma omp parallel if(NumPairs() >= SpmvParallelThreshold())
    {
      const int tid = omp_get_thread_num(), threads = omp_get_num_threads();
      const Index first = AtomRangeBoundary(pairOff, atomBegin, atomEnd, tid, threads);
      const Index last = AtomRangeBoundary(pairOff, atomBegin, atomEnd, tid + 1, threads);
      for (Index row = first; row < last; ++row) {
        Scalar sum(0);
        for (Index p = pairOff[row]; p < pairOff[row + 1]; ++p) {
          sum += BlockProduct(dictVal[blk[p]], xp[nbr[p]]);
        }
        if (diag != nullptr) sum += BlockProduct(diag[row], xp[row]);
        yp[row] = BlockProduct(a, sum) + BlockProduct(b, yp[row]);
      }
    }
    return;
  }

  if (maxRows == 2 && Dimension() / numAtoms == 2 && Dimension() % numAtoms == 0) {
    // A uniform 2x2 block product keeps both accumulators in registers. The generic
    // variable-width path remains available for mixed atoms and all other block sizes.
#pragma omp parallel if(NumPairs() >= SpmvParallelThreshold())
    {
      const int tid = omp_get_thread_num(), threads = omp_get_num_threads();
      const Index first = AtomRangeBoundary(pairOff, atomBegin, atomEnd, tid, threads);
      const Index last = AtomRangeBoundary(pairOff, atomBegin, atomEnd, tid + 1, threads);
      for (Index atom = first; atom < last; ++atom) {
        Scalar sum0(0), sum1(0);
        for (Index p = pairOff[atom]; p < pairOff[atom + 1]; ++p) {
          const Scalar* block = dictVal + dictOff[blk[p]];
          const Index column = 2 * static_cast<Index>(nbr[p]);
          const Scalar x0 = xp[column], x1 = xp[column + 1];
          if (block[0] != Scalar(0)) sum0 += BlockProduct(block[0], x0);
          if (block[1] != Scalar(0)) sum0 += BlockProduct(block[1], x1);
          if (block[2] != Scalar(0)) sum1 += BlockProduct(block[2], x0);
          if (block[3] != Scalar(0)) sum1 += BlockProduct(block[3], x1);
        }
        const Index row = 2 * atom;
        if (diag != nullptr) {
          sum0 += BlockProduct(diag[row], xp[row]);
          sum1 += BlockProduct(diag[row + 1], xp[row + 1]);
        }
        yp[row] = BlockProduct(a, sum0) + BlockProduct(b, yp[row]);
        yp[row + 1] = BlockProduct(a, sum1) + BlockProduct(b, yp[row + 1]);
      }
    }
    return;
  }

#pragma omp parallel if(NumPairs() >= SpmvParallelThreshold())
  {
    const int tid = omp_get_thread_num();
    const int nthreads = omp_get_num_threads();
    const Index first = AtomRangeBoundary(pairOff, atomBegin, atomEnd, tid, nthreads);
    const Index last = AtomRangeBoundary(pairOff, atomBegin, atomEnd, tid + 1, nthreads);
    Vector accBuf(static_cast<std::size_t>(maxRows));
    for (Index atom = first; atom < last; ++atom) {
      const Index rowBegin = atomOff[atom];
      const Index rowEnd = atomOff[atom + 1];
      const int rows = static_cast<int>(rowEnd - rowBegin);
      std::fill_n(accBuf.begin(), rows, Scalar(0.0, 0.0));
      for (Index p = pairOff[atom]; p < pairOff[atom + 1]; ++p) {
        const std::int32_t d = blk[static_cast<std::size_t>(p)];
        const Scalar* const block = dictVal + dictOff[d];
        const int cols = dictCols[d];
        const Scalar* const xb = xp + atomOff[nbr[static_cast<std::size_t>(p)]];
        for (int i = 0; i < rows; ++i) {
          const Scalar* const blockRow = block + static_cast<std::size_t>(i) * cols;
          Scalar sum = accBuf[static_cast<std::size_t>(i)];
          for (int j = 0; j < cols; ++j) {
            sum += BlockProduct(blockRow[j], xb[j]);
          }
          accBuf[static_cast<std::size_t>(i)] = sum;
        }
      }
      if (diag != nullptr) {
        for (int i = 0; i < rows; ++i) {
          accBuf[static_cast<std::size_t>(i)] += BlockProduct(diag[rowBegin + i], xp[rowBegin + i]);
        }
      }
      for (int i = 0; i < rows; ++i) {
        yp[rowBegin + i] = BlockProduct(a, accBuf[static_cast<std::size_t>(i)]) +
                          BlockProduct(b, yp[rowBegin + i]);
      }
    }
  }
}

void AtomBlockOperator::MultiplyBlock(Scalar a, const VectorBlock& x, Scalar b,
                                      VectorBlock& y) const {
  const Index n = Dimension();
  if (&x == &y) {
    throw std::invalid_argument(
        "AtomBlockOperator::MultiplyBlock: input and output blocks must not alias");
  }
  if (x.Rows() != n || y.Rows() != n || x.Cols() != y.Cols()) {
    throw std::invalid_argument("AtomBlockOperator::MultiplyBlock: block dimension mismatch");
  }
  const int cols = x.Cols();
  const Index numAtoms = NumAtoms();
  if (n == 0 || cols == 0 || numAtoms == 0) {
    return;
  }
  if (cols == 1) {
    // Same per-element arithmetic as Multiply (bit-identical), directly on the column.
    MultiplyKernel(a, x.Data(), b, y.Data());
    return;
  }
  MultiplyBlockRange(a, x.Data(), b, y.Data(), cols, 0, NumAtoms());
}

void AtomBlockOperator::MultiplyBlockRange(Scalar a, const Scalar* xp, Scalar b, Scalar* yp,
                                           int cols, Index atomBegin, Index atomEnd) const {
  const Index n = Dimension();
  const Index numAtoms = NumAtoms();
  if (atomBegin < 0 || atomEnd > numAtoms || atomBegin > atomEnd) {
    throw std::invalid_argument("AtomBlockOperator::MultiplyBlockRange: invalid atom range");
  }
  if (n == 0 || cols == 0 || numAtoms == 0 || atomBegin == atomEnd) {
    return;
  }
  if (cols == 1) {
    MultiplyRange(a, xp, b, yp, atomBegin, atomEnd);
    return;
  }
  if (maxBlockRows_ <= 2 && maxBlockCols_ <= 2) {
    // Small blocks do too little arithmetic to amortise state-lane staging. Stream each
    // contiguous column through the same kernel; no transpose or alternate layout is used.
    for (int col = 0; col < cols; ++col) {
      MultiplyRange(a, xp + static_cast<std::size_t>(col) * n, b,
                    yp + static_cast<std::size_t>(col) * n, atomBegin, atomEnd);
    }
    return;
  }

  const Index* const atomOff = atomOffset_.data();
  const Index* const pairOff = pairOffset_.data();
  const std::int32_t* const nbr = nbr_.data();
  const std::int32_t* const blk = blockId_.data();
  const Index* const dictOff = dictOffset_.data();
  const std::int32_t* const dictCols = dictCols_.data();
  const Scalar* const dictVal = dictValues_.data();
  const Scalar* const diag = diagonal_.empty() ? nullptr : diagonal_.data();
  const Index maxRows = maxBlockRows_;
  const Index maxCols = maxBlockCols_;

#pragma omp parallel if(NumPairs() >= SpmvParallelThreshold())
  {
    const int tid = omp_get_thread_num();
    const int nthreads = omp_get_num_threads();
    const Index first = AtomRangeBoundary(pairOff, atomBegin, atomEnd, tid, nthreads);
    const Index last = AtomRangeBoundary(pairOff, atomBegin, atomEnd, tid + 1, nthreads);
    // Per-thread scratch: staged neighbour slice [j][s] and accumulator [i][s], state
    // contiguous so the inner loop over block columns vectorises per lane.
    Vector xBuf(static_cast<std::size_t>(maxCols) * static_cast<std::size_t>(cols));
    Vector accBuf(static_cast<std::size_t>(maxRows) * static_cast<std::size_t>(cols));
    for (Index atom = first; atom < last; ++atom) {
      const Index rowBegin = atomOff[atom];
      const int rows = static_cast<int>(atomOff[atom + 1] - rowBegin);
      std::fill_n(accBuf.begin(), static_cast<std::size_t>(rows) * cols, Scalar(0.0, 0.0));
      for (Index p = pairOff[atom]; p < pairOff[atom + 1]; ++p) {
        const std::int32_t d = blk[static_cast<std::size_t>(p)];
        const Scalar* const block = dictVal + dictOff[d];
        const int width = dictCols[d];
        const Scalar* const xb = xp + atomOff[nbr[static_cast<std::size_t>(p)]];
        for (int j = 0; j < width; ++j) {
          Scalar* const dst = xBuf.data() + static_cast<std::size_t>(j) * cols;
          const Scalar* const src = xb + j;
#pragma omp simd
          for (int s = 0; s < cols; ++s) {
            dst[s] = src[static_cast<std::size_t>(s) * n];
          }
        }
        for (int i = 0; i < rows; ++i) {
          const Scalar* const blockRow = block + static_cast<std::size_t>(i) * width;
          Scalar* const acc = accBuf.data() + static_cast<std::size_t>(i) * cols;
          for (int j = 0; j < width; ++j) {
            const Scalar value = blockRow[j];
            const Scalar* const xRow = xBuf.data() + static_cast<std::size_t>(j) * cols;
#pragma omp simd
            for (int s = 0; s < cols; ++s) acc[s] += BlockProduct(value, xRow[s]);
          }
        }
      }
      if (diag != nullptr) {
        for (int i = 0; i < rows; ++i) {
          const Scalar dval = diag[rowBegin + i];
          const Scalar* const xRow = xp + rowBegin + i;
          Scalar* const acc = accBuf.data() + static_cast<std::size_t>(i) * cols;
#pragma omp simd
          for (int s = 0; s < cols; ++s) {
            acc[s] += BlockProduct(dval, xRow[static_cast<std::size_t>(s) * n]);
          }
        }
      }
      for (int i = 0; i < rows; ++i) {
        const Scalar* const acc = accBuf.data() + static_cast<std::size_t>(i) * cols;
        Scalar* const yRow = yp + rowBegin + i;
#pragma omp simd
        for (int s = 0; s < cols; ++s) {
          Scalar& out = yRow[static_cast<std::size_t>(s) * n];
          out = BlockProduct(a, acc[s]) + BlockProduct(b, out);
        }
      }
    }
  }
}

#ifdef LSQUANT_NATIVE_BUILD
void AtomBlockOperator::PlaceOnThreads() {
  const Index numAtoms = NumAtoms();
  const Index numPairs = NumPairs();
  if (numAtoms == 0 || numPairs == 0) return;
  std::vector<std::int32_t> placedNbr(nbr_.size());
  std::vector<std::int32_t> placedBlk(blockId_.size());
  DiscardInteriorPages(placedNbr.data(), placedNbr.size() * sizeof(std::int32_t));
  DiscardInteriorPages(placedBlk.data(), placedBlk.size() * sizeof(std::int32_t));
#pragma omp parallel if(numPairs >= SpmvParallelThreshold())
  {
    const Index first =
        RowBoundary(pairOffset_.data(), numAtoms, omp_get_thread_num(), omp_get_num_threads());
    const Index last = RowBoundary(pairOffset_.data(), numAtoms, omp_get_thread_num() + 1,
                                   omp_get_num_threads());
    for (Index p = pairOffset_[first]; p < pairOffset_[last]; ++p) {
      placedNbr[static_cast<std::size_t>(p)] = nbr_[static_cast<std::size_t>(p)];
      placedBlk[static_cast<std::size_t>(p)] = blockId_[static_cast<std::size_t>(p)];
    }
  }
  nbr_.swap(placedNbr);
  blockId_.swap(placedBlk);
}
#endif

}  // namespace lsquant
