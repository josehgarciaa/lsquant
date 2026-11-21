/**
 * @file chebyshev_vectors.cpp
 * @brief Implementation of ChebyshevVectors.
 */
#include "core/chebyshev/chebyshev_vectors.hpp"

#include <algorithm>
#include <stdexcept>

namespace lsquant {

ChebyshevVectors::ChebyshevVectors(Index dimension, int batchSize)
    : dimension_(dimension), batchSize_(batchSize), firstOrder_(0) {
  if (dimension <= 0) {
    throw std::invalid_argument("ChebyshevVectors: dimension must be positive");
  }
  if (batchSize <= 0) {
    throw std::invalid_argument("ChebyshevVectors: batchSize must be positive");
  }
  vectors_ = VectorBlock(dimension, batchSize);
  scratch_ = VectorBlock(dimension, batchSize);
}

void ChebyshevVectors::Fill(ChebyshevRecursion& recursion) {
  if (static_cast<Index>(recursion.Current().size()) != dimension_) {
    throw std::invalid_argument("ChebyshevVectors::Fill: recursion dimension mismatch");
  }
  firstOrder_ = recursion.Order();
  for (int m = 0; m < batchSize_; ++m) {
    std::copy_n(recursion.Current().data(), static_cast<std::size_t>(dimension_), Column(m));
    recursion.Step();
  }
}

void ChebyshevVectors::ApplyToAll(const LinearOperator& A) {
  if (A.Dimension() != dimension_) {
    throw std::invalid_argument("ChebyshevVectors::ApplyToAll: operator dimension mismatch");
  }
  if (A.IsIdentity()) {
    return;
  }
  A.MultiplyBlock(Scalar(1.0, 0.0), vectors_, Scalar(0.0, 0.0), scratch_);
  vectors_.Swap(scratch_);
}

Scalar* ChebyshevVectors::Column(int m) { return vectors_.Column(m); }

const Scalar* ChebyshevVectors::Column(int m) const { return vectors_.Column(m); }

VectorBlock& ChebyshevVectors::Block() { return vectors_; }

const VectorBlock& ChebyshevVectors::Block() const { return vectors_; }

int ChebyshevVectors::Size() const { return batchSize_; }

Index ChebyshevVectors::Dimension() const { return dimension_; }

int ChebyshevVectors::FirstOrder() const { return firstOrder_; }

}  // namespace lsquant
