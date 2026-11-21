/**
 * @file ginkgo_operator.cpp
 * @brief Implementation of GinkgoOperator.
 */
#ifdef LSQUANT_WITH_GINKGO
#include "core/algebra/ginkgo_operator.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

#include <ginkgo/ginkgo.hpp>

namespace lsquant {

namespace {

using GkoDense = gko::matrix::Dense<Scalar>;

// Wrap a contiguous host buffer as a Ginkgo dense view without copying. Ginkgo dense
// matrices are row-major: element (i, j) lives at ptr[i * stride + j].
std::unique_ptr<GkoDense> WrapHost(const std::shared_ptr<gko::Executor>& host, Scalar* ptr,
                                   Index rows, Index cols, Index stride) {
  const Index last = (rows - 1) * stride + cols;
  return GkoDense::create(
      host, gko::dim<2>(static_cast<gko::size_type>(rows), static_cast<gko::size_type>(cols)),
      gko::array<Scalar>::view(host, static_cast<gko::size_type>(last), ptr),
      static_cast<gko::size_type>(stride));
}

std::unique_ptr<GkoDense> WrapHostConst(const std::shared_ptr<gko::Executor>& host,
                                        const Scalar* ptr, Index rows, Index cols,
                                        Index stride) {
  return WrapHost(host, const_cast<Scalar*>(ptr), rows, cols, stride);
}

template <typename IT>
std::shared_ptr<gko::LinOp> BuildMatrix(const SparseMatrix& matrix,
                                        std::shared_ptr<gko::Executor> exec,
                                        const std::shared_ptr<gko::Executor>& host,
                                        const std::string& format) {
  using Csr = gko::matrix::Csr<Scalar, IT>;
  const auto& e = matrix.Eigen();
  const Index n = e.rows();
  const Index nnz = e.nonZeros();
  const Index* rp = e.outerIndexPtr();
  const Index* cl = e.innerIndexPtr();
  const Scalar* vl = e.valuePtr();

  gko::array<Scalar> values(host, static_cast<gko::size_type>(nnz));
  gko::array<IT> columns(host, static_cast<gko::size_type>(nnz));
  gko::array<IT> rowPointers(host, static_cast<gko::size_type>(n + 1));
  std::copy_n(vl, nnz, values.get_data());
  std::transform(cl, cl + nnz, columns.get_data(),
                 [](Index c) { return static_cast<IT>(c); });
  std::transform(rp, rp + n + 1, rowPointers.get_data(),
                 [](Index p) { return static_cast<IT>(p); });

  auto hcsr = Csr::create(host, gko::dim<2>(static_cast<gko::size_type>(n),
                                            static_cast<gko::size_type>(n)),
                          std::move(values), std::move(columns), std::move(rowPointers));

  if (format == "csr") {
    auto out = Csr::create(exec);
    out->copy_from(hcsr);
    return out;
  }
  if (format == "sellp") {
    auto h = gko::matrix::Sellp<Scalar, IT>::create(host);
    hcsr->convert_to(h);
    auto out = gko::matrix::Sellp<Scalar, IT>::create(exec);
    out->copy_from(h);
    return out;
  }
  if (format == "ell") {
    auto h = gko::matrix::Ell<Scalar, IT>::create(host);
    hcsr->convert_to(h);
    auto out = gko::matrix::Ell<Scalar, IT>::create(exec);
    out->copy_from(h);
    return out;
  }
  if (format == "hybrid") {
    auto h = gko::matrix::Hybrid<Scalar, IT>::create(host);
    hcsr->convert_to(h);
    auto out = gko::matrix::Hybrid<Scalar, IT>::create(exec);
    out->copy_from(h);
    return out;
  }
  if (format == "coo") {
    auto h = gko::matrix::Coo<Scalar, IT>::create(host);
    hcsr->convert_to(h);
    auto out = gko::matrix::Coo<Scalar, IT>::create(exec);
    out->copy_from(h);
    return out;
  }
  throw std::invalid_argument("GinkgoOperator: unknown format '" + format +
                              "'; expected csr, sellp, ell, hybrid, or coo");
}

}  // namespace

struct GinkgoOperator::Impl {
  std::shared_ptr<gko::Executor> host;
  std::shared_ptr<gko::Executor> exec;
  std::shared_ptr<gko::LinOp> matrix;
  bool onCuda = false;
};

GinkgoOperator::GinkgoOperator(const SparseMatrix& matrix, std::string executor,
                               std::string format)
    : impl_(std::make_shared<Impl>()),
      dimension_(matrix.Dimension()),
      isHermitian_(matrix.IsHermitian()),
      executorName_(std::move(executor)),
      formatName_(std::move(format)) {
  if (dimension_ <= 0) {
    throw std::invalid_argument("GinkgoOperator: matrix dimension must be positive");
  }
  impl_->host = gko::OmpExecutor::create();
  if (executorName_ == "omp") {
    impl_->exec = impl_->host;
    impl_->onCuda = false;
  } else if (executorName_ == "cuda") {
#ifdef LSQUANT_WITH_CUDA
    impl_->exec = gko::CudaExecutor::create(0, impl_->host);
    impl_->onCuda = true;
#else
    throw std::runtime_error(
        "GinkgoOperator: executor 'cuda' requested but Ginkgo was built without CUDA");
#endif
  } else {
    throw std::invalid_argument("GinkgoOperator: unknown executor '" + executorName_ +
                                "'; expected omp or cuda");
  }
  const Index nnz = matrix.NonZeros();
  const Index maxIndex = std::max(dimension_ + 1, nnz);
  if (maxIndex <= static_cast<Index>(std::numeric_limits<std::int32_t>::max())) {
    impl_->matrix =
        BuildMatrix<std::int32_t>(matrix, impl_->exec, impl_->host, formatName_);
  } else {
    impl_->matrix =
        BuildMatrix<std::int64_t>(matrix, impl_->exec, impl_->host, formatName_);
  }
}

Index GinkgoOperator::Dimension() const { return dimension_; }

void GinkgoOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index n = dimension_;
  if (static_cast<Index>(x.size()) != n || static_cast<Index>(y.size()) != n) {
    throw std::invalid_argument("GinkgoOperator::Multiply: vector dimension mismatch");
  }
  if (&x == &y) {
    throw std::invalid_argument(
        "GinkgoOperator::Multiply: input and output vectors must not alias");
  }
  if (n == 0) {
    return;
  }
  auto alpha = gko::initialize<GkoDense>({a}, impl_->exec);
  auto beta = gko::initialize<GkoDense>({b}, impl_->exec);
  if (!impl_->onCuda) {
    auto xv = WrapHostConst(impl_->exec, x.data(), n, 1, 1);
    auto yv = WrapHost(impl_->exec, y.data(), n, 1, 1);
    impl_->matrix->apply(alpha, xv, beta, yv);
    return;
  }
  // CUDA: stage through device copies (the device-resident engine avoids this path).
  auto hx = WrapHostConst(impl_->host, x.data(), n, 1, 1);
  auto dx = gko::clone(impl_->exec, hx);
  auto dy = GkoDense::create(impl_->exec, gko::dim<2>(static_cast<gko::size_type>(n), 1));
  if (b != Scalar(0.0, 0.0)) {
    auto hy = WrapHost(impl_->host, y.data(), n, 1, 1);
    dy->copy_from(hy);
  }
  impl_->matrix->apply(alpha, dx, beta, dy);
  impl_->exec->synchronize();
  auto hy = WrapHost(impl_->host, y.data(), n, 1, 1);
  hy->copy_from(dy);
}

void GinkgoOperator::MultiplyBlock(Scalar a, const VectorBlock& x, Scalar b,
                                   VectorBlock& y) const {
  const Index n = dimension_;
  if (&x == &y) {
    throw std::invalid_argument(
        "GinkgoOperator::MultiplyBlock: input and output blocks must not alias");
  }
  if (x.Rows() != n || y.Rows() != n || x.Cols() != y.Cols()) {
    throw std::invalid_argument("GinkgoOperator::MultiplyBlock: block dimension mismatch");
  }
  const int cols = x.Cols();
  if (n == 0 || cols == 0) {
    return;
  }
  // Ginkgo Dense is row-major; contiguous columns are zero-copy N x 1 views.
  auto alpha = gko::initialize<GkoDense>({a}, impl_->exec);
  auto beta = gko::initialize<GkoDense>({b}, impl_->exec);
  if (impl_->onCuda) {
    throw std::runtime_error("GinkgoOperator::MultiplyBlock: host blocks not supported on cuda");
  }
  for (int j = 0; j < cols; ++j) {
    auto xv = WrapHostConst(impl_->exec, x.Column(j), n, 1, 1);
    auto yv = WrapHost(impl_->exec, y.Column(j), n, 1, 1);
    impl_->matrix->apply(alpha, xv, beta, yv);
  }
}

bool GinkgoOperator::IsHermitian() const { return isHermitian_; }

std::shared_ptr<gko::Executor> GinkgoOperator::GinkgoExecutor() const { return impl_->exec; }

void GinkgoOperator::ApplyDevice(Scalar a, const gko::matrix::Dense<Scalar>* x, Scalar b,
                                 gko::matrix::Dense<Scalar>* y) const {
  auto alpha = gko::initialize<GkoDense>({a}, impl_->exec);
  auto beta = gko::initialize<GkoDense>({b}, impl_->exec);
  impl_->matrix->apply(alpha, x, beta, y);
}

}  // namespace lsquant
#endif  // LSQUANT_WITH_GINKGO
