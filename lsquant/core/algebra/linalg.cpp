/**
 * @file linalg.cpp
 * @brief Implementation of dense vector linear algebra wrappers over Eigen::Map.
 */
#include "core/algebra/linalg.hpp"

#include <Eigen/Core>
#include <algorithm>
#include <stdexcept>
#include "core/algebra/memory_placement.hpp"

namespace lsquant {
namespace linalg {

#ifdef LSQUANT_NATIVE_BUILD
namespace {

// Fixed chunk length of the deterministic parallel dot reduction. The chunk grid depends
// only on the vector length, never on the OpenMP thread count, so partial sums are computed
// and combined in exactly the same order for 1 and 128 threads.
constexpr Index kReductionChunk = 1 << 17;

}  // namespace
#endif

// Row-chunk length of the native split-K moment table. Sized so the chunk count keeps all
// cores of a 128-core socket busy at N ~ 5e5 and the per-chunk partial table stays
// cache-resident; capped at 256 chunks so the partial storage stays small at large N.
// Depends only on N.
Index TableChunkRows(Index dimension) {
  const Index capped = (dimension + 255) / 256;
  return std::max<Index>(4096, capped);
}

std::size_t TableScratchSize(Index dimension, int leftCount, int rightCount) {
  if (dimension <= 0 || leftCount <= 0 || rightCount <= 0) {
    return 0;
  }
  const Index chunkRows = TableChunkRows(dimension);
  const Index numChunks = (dimension + chunkRows - 1) / chunkRows;
  return static_cast<std::size_t>(numChunks) * static_cast<std::size_t>(leftCount) *
         static_cast<std::size_t>(rightCount);
}

#ifdef LSQUANT_NATIVE_BUILD
Scalar DotDeterministic(const Scalar* x, const Scalar* y, Index n) {
  if (n <= 0) {
    return Scalar(0.0, 0.0);
  }
  const Index numChunks = (n + kReductionChunk - 1) / kReductionChunk;
  std::vector<Scalar> partials(static_cast<std::size_t>(numChunks));
#pragma omp parallel for schedule(static)
  for (Index c = 0; c < numChunks; ++c) {
    const Index begin = c * kReductionChunk;
    const Index len = std::min(kReductionChunk, n - begin);
    Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x + begin, static_cast<Eigen::Index>(len));
    Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_y(y + begin, static_cast<Eigen::Index>(len));
    partials[static_cast<std::size_t>(c)] = eig_x.dot(eig_y);
  }
  Scalar total(0.0, 0.0);
  for (const Scalar& p : partials) {
    total += p;
  }
  return total;
}
#endif

void axpy(Scalar a, const Vector& x, Vector& y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("linalg::axpy: vector size mismatch");
  }
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  Eigen::Map<Eigen::Vector<Scalar, -1>> eig_y(y.data(), static_cast<Eigen::Index>(y.size()));
  eig_y = a * eig_x + eig_y;
}

void scal(Scalar a, Vector& x) {
  Eigen::Map<Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  eig_x *= a;
}

void copy(const Vector& x, Vector& y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("linalg::copy: vector size mismatch");
  }
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  Eigen::Map<Eigen::Vector<Scalar, -1>> eig_y(y.data(), static_cast<Eigen::Index>(y.size()));
  eig_y = eig_x;
}

Scalar dot(const Vector& x, const Vector& y) {
  if (x.size() != y.size()) {
    throw std::invalid_argument("linalg::dot: vector size mismatch");
  }
  return dot(x.data(), y.data(), static_cast<Index>(x.size()));
}

Scalar dot(const Scalar* x, const Scalar* y, Index n) {
#ifdef LSQUANT_NATIVE_BUILD
  return DotDeterministic(x, y, n);
#else
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x, static_cast<Eigen::Index>(n));
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_y(y, static_cast<Eigen::Index>(n));
  return eig_x.dot(eig_y);
#endif
}

double norm2(const Vector& x) {
  Eigen::Map<const Eigen::Vector<Scalar, -1>> eig_x(x.data(), static_cast<Eigen::Index>(x.size()));
  return eig_x.squaredNorm();
}

std::vector<Scalar> MomentTable(const std::vector<const Vector*>& left,
                               const std::vector<const Vector*>& right) {
  if (left.empty() || right.empty()) return {};
  if (!left.front()) throw std::invalid_argument("MomentTable: null left vector");
  const Index dimension = static_cast<Index>(left.front()->size());
  for (const Vector* v : left) {
    if (!v || static_cast<Index>(v->size()) != dimension) throw std::invalid_argument("MomentTable: invalid left vector");
  }
  for (const Vector* v : right) {
    if (!v || static_cast<Index>(v->size()) != dimension) throw std::invalid_argument("MomentTable: invalid right vector");
  }
  using Dense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  Dense a(dimension, left.size()), b(dimension, right.size());
#ifdef LSQUANT_NATIVE_BUILD
  DiscardInteriorPages(a.data(), static_cast<std::size_t>(a.size()) * sizeof(Scalar));
  DiscardInteriorPages(b.data(), static_cast<std::size_t>(b.size()) * sizeof(Scalar));
#endif
#pragma omp parallel for schedule(static)
  for (Index row = 0; row < dimension; ++row) {
    for (std::size_t j = 0; j < left.size(); ++j) a(row, j) = (*left[j])[row];
    for (std::size_t j = 0; j < right.size(); ++j) b(row, j) = (*right[j])[row];
  }
  using RowDense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  std::vector<Scalar> values(left.size() * right.size());
#ifdef LSQUANT_NATIVE_BUILD
  DiscardInteriorPages(values.data(), values.size() * sizeof(Scalar));
#endif
  Eigen::Map<RowDense> result(values.data(), left.size(), right.size());
  constexpr Index tile = 16;
  const Index rows = static_cast<Index>(left.size()), cols = static_cast<Index>(right.size());
  const Index rowTiles = (rows + tile - 1) / tile, colTiles = (cols + tile - 1) / tile;
#pragma omp parallel for schedule(static)
  for (Index task = 0; task < rowTiles * colTiles; ++task) {
    const Index row = (task / colTiles) * tile, col = (task % colTiles) * tile;
    const Index nr = std::min(tile, rows - row), nc = std::min(tile, cols - col);
    result.block(row, col, nr, nc).noalias() = a.middleCols(row, nr).adjoint() * b.middleCols(col, nc);
  }
  return values;
}

std::vector<Scalar> MomentTableBlock(const Scalar* left, int leftCount, const Scalar* right,
                                     int rightCount, Index dimension, Scalar* scratch) {
  if (leftCount < 0 || rightCount < 0 || dimension < 0) {
    throw std::invalid_argument("MomentTableBlock: negative extent");
  }
  if (leftCount == 0 || rightCount == 0) return {};
  if (left == nullptr || right == nullptr) {
    throw std::invalid_argument("MomentTableBlock: null block pointer");
  }
  using Dense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using RowDense = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  std::vector<Scalar> values(static_cast<std::size_t>(leftCount) *
                             static_cast<std::size_t>(rightCount));
#ifdef LSQUANT_NATIVE_BUILD
  DiscardInteriorPages(values.data(), values.size() * sizeof(Scalar));
#else
  (void)scratch;  // the golden profile keeps the serial tiled GEMM and needs no partials
#endif
  Eigen::Map<RowDense> result(values.data(), leftCount, rightCount);
#ifdef LSQUANT_NATIVE_BUILD
  // Parallel, thread-count-invariant split-K reduction: fixed row chunks (independent of
  // the OpenMP thread count), one single-thread partial GEMM per chunk computed in
  // parallel, partial tables summed in chunk order (chunk loop outer, elements parallel).
  const Index chunkRows = TableChunkRows(dimension);
  const Index numChunks = (dimension + chunkRows - 1) / chunkRows;
  const std::size_t tileSize =
      static_cast<std::size_t>(leftCount) * static_cast<std::size_t>(rightCount);
  const std::size_t needed = static_cast<std::size_t>(numChunks) * tileSize;
  // Caller-provided workspace when given (the engine issues one table per left window; a
  // fresh 33 MB allocation per call costs ~10 ms of page faults at N ~ 5e5). The used prefix
  // is fully overwritten on every call.
  std::vector<Scalar> owned;
  if (scratch == nullptr) {
    owned.assign(needed, Scalar(0.0, 0.0));
    scratch = owned.data();
  }
  Scalar* const partials = scratch;
  using Strided = Eigen::Map<const Dense, 0, Eigen::OuterStride<>>;
  const Strided a(left, dimension, leftCount, Eigen::OuterStride<>(dimension));
  const Strided b(right, dimension, rightCount, Eigen::OuterStride<>(dimension));
#pragma omp parallel for schedule(static)
  for (Index c = 0; c < numChunks; ++c) {
    const Index begin = c * chunkRows;
    const Index len = std::min(chunkRows, dimension - begin);
    Eigen::Map<Dense> partial(partials + static_cast<std::size_t>(c) * tileSize,
                              leftCount, rightCount);
    partial.noalias() =
        a.middleRows(begin, len).adjoint() * b.middleRows(begin, len);
  }
  // Sum the partial tables in chunk order: chunk loop outer (sequential), element loop
  // inner (parallel; each partial tile is streamed contiguously). The element order of the
  // additions is fixed, so the result is thread-count invariant.
  for (Index c = 0; c < numChunks; ++c) {
    const Scalar* p = partials + static_cast<std::size_t>(c) * tileSize;
#pragma omp parallel for schedule(static)
    for (Index k = 0; k < static_cast<Index>(tileSize); ++k) {
      // partials are column-major (leftCount x rightCount); result is row-major.
      const std::size_t m = static_cast<std::size_t>(k) % static_cast<std::size_t>(leftCount);
      const std::size_t n = static_cast<std::size_t>(k) / static_cast<std::size_t>(leftCount);
      Scalar& out = result(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(n));
      if (c == 0) {
        out = p[k];
      } else {
        out += p[k];
      }
    }
  }
#else
  Eigen::Map<const Dense> a(left, dimension, leftCount);
  Eigen::Map<const Dense> b(right, dimension, rightCount);
  constexpr Index tile = 16;
  const Index rows = leftCount, cols = rightCount;
  const Index rowTiles = (rows + tile - 1) / tile, colTiles = (cols + tile - 1) / tile;
#pragma omp parallel for schedule(static)
  for (Index task = 0; task < rowTiles * colTiles; ++task) {
    const Index row = (task / colTiles) * tile, col = (task % colTiles) * tile;
    const Index nr = std::min(tile, rows - row), nc = std::min(tile, cols - col);
    result.block(row, col, nr, nc).noalias() =
        a.middleCols(row, nr).adjoint() * b.middleCols(col, nc);
  }
#endif
  return values;
}

}  // namespace linalg
}  // namespace lsquant
