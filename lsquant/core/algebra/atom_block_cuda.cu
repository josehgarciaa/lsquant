/**
 * @file atom_block_cuda.cu
 * @brief CUDA warp-per-atom kernel and host wrappers for the atom-block operator.
 *
 * One warp owns one atom row-block per iteration; neighbour slices `x_b` (`W_b x S`, columns
 * strided by `N`) are staged through shared memory coalesced (orbital index fastest across
 * lanes), the dictionary is staged in dynamic shared memory when it fits `kMaxSharedDict`
 * bytes (crystals) and read from global memory otherwise (amorphous; L2-resident). Each lane
 * owns the accumulator elements `e = lane, lane+32, ...` of the flattened `(row i, column s)`
 * space and sums serially over pairs and block columns, so the per-element operation order
 * matches the OpenMP kernel.
 */
#include "core/algebra/atom_block_cuda.hpp"

#ifdef LSQUANT_WITH_CUDA

#include "core/algebra/atom_block_validate.hpp"

#include <cuda_runtime.h>
#include <cuComplex.h>

#include <algorithm>
#include <stdexcept>
#include <string>

namespace lsquant {

namespace {

constexpr int kWarpSize = 32;
constexpr int kMaxAccumPerLane = 8;      // supports W_a * S <= 256 per atom
constexpr std::size_t kMaxSharedDict = 40 * 1024;  // bytes of dictionary staged in shared
constexpr int kBlockThreads = 256;
constexpr int kMaxStateChunk = 16;       // S chunked to keep the register bound

#define LSQ_CUDA_CHECK(call)                                                        \
  do {                                                                              \
    const cudaError_t err_ = (call);                                                \
    if (err_ != cudaSuccess) {                                                      \
      throw std::runtime_error(std::string("CUDA error at ") + __FILE__ + ":" +     \
                               std::to_string(__LINE__) + ": " +                    \
                               cudaGetErrorString(err_));                           \
    }                                                                               \
  } while (0)

__device__ __forceinline__ cuDoubleComplex cmul(cuDoubleComplex a, cuDoubleComplex b) {
  return cuCmul(a, b);
}

__device__ __forceinline__ cuDoubleComplex cadd(cuDoubleComplex a, cuDoubleComplex b) {
  return cuCadd(a, b);
}

__global__ void AtomBlockKernel(cuDoubleComplex a, const cuDoubleComplex* __restrict__ x,
                                cuDoubleComplex b, cuDoubleComplex* __restrict__ y,
                                int cols, Index n, Index numAtoms,
                                const Index* __restrict__ atomOffset,
                                const Index* __restrict__ pairOffset,
                                const std::int32_t* __restrict__ nbr,
                                const std::int32_t* __restrict__ blockId,
                                const Index* __restrict__ dictOffset,
                                const std::int32_t* __restrict__ dictCols,
                                const cuDoubleComplex* __restrict__ dictValues,
                                const cuDoubleComplex* __restrict__ diagonal,
                                int dictInShared, Index dictScalars, int maxCols) {
  extern __shared__ cuDoubleComplex smem[];
  cuDoubleComplex* sDict = smem;
  cuDoubleComplex* sX = smem + (dictInShared ? dictScalars : 0);
  const int warpsPerBlock = blockDim.x / kWarpSize;
  if (dictInShared) {
    for (Index i = threadIdx.x; i < dictScalars; i += blockDim.x) sDict[i] = dictValues[i];
  }
  __syncthreads();
  const int warp = threadIdx.x / kWarpSize;
  const int lane = threadIdx.x % kWarpSize;
  cuDoubleComplex* warpX = sX + static_cast<std::size_t>(warp) * maxCols * cols;

  for (Index atom = static_cast<Index>(blockIdx.x) * warpsPerBlock + warp; atom < numAtoms;
       atom += static_cast<Index>(gridDim.x) * warpsPerBlock) {
    const Index rowBegin = atomOffset[atom];
    const int rows = static_cast<int>(atomOffset[atom + 1] - rowBegin);
    const int nelem = rows * cols;
    cuDoubleComplex acc[kMaxAccumPerLane];
    for (int e = lane, k = 0; e < nelem; e += kWarpSize, ++k) {
      acc[k] = make_cuDoubleComplex(0.0, 0.0);
    }
    for (Index p = pairOffset[atom]; p < pairOffset[atom + 1]; ++p) {
      const std::int32_t d = blockId[p];
      const int width = dictCols[d];
      const cuDoubleComplex* block =
          (dictInShared ? sDict : dictValues) + dictOffset[d];
      const cuDoubleComplex* xb = x + atomOffset[nbr[p]];
      // stage x_b transposed ([s][j]) so consecutive lanes read consecutive orbitals
      for (int e = lane; e < width * cols; e += kWarpSize) {
        const int j = e % width;
        const int s = e / width;
        warpX[s * width + j] = xb[static_cast<std::size_t>(s) * n + j];
      }
      __syncwarp();
      for (int e = lane, k = 0; e < nelem; e += kWarpSize, ++k) {
        const int i = e / cols;
        const int s = e % cols;
        const cuDoubleComplex* blockRow = block + static_cast<std::size_t>(i) * width;
        const cuDoubleComplex* xRow = warpX + s * width;
        cuDoubleComplex sum = acc[k];
        for (int j = 0; j < width; ++j) sum = cadd(sum, cmul(blockRow[j], xRow[j]));
        acc[k] = sum;
      }
      __syncwarp();
    }
    if (diagonal != nullptr) {
      for (int e = lane, k = 0; e < nelem; e += kWarpSize, ++k) {
        const int i = e / cols;
        const int s = e % cols;
        acc[k] = cadd(acc[k], cmul(diagonal[rowBegin + i],
                                   x[(rowBegin + i) + static_cast<std::size_t>(s) * n]));
      }
    }
    for (int e = lane, k = 0; e < nelem; e += kWarpSize, ++k) {
      const int i = e / cols;
      const int s = e % cols;
      cuDoubleComplex& out = y[(rowBegin + i) + static_cast<std::size_t>(s) * n];
      out = cadd(cmul(a, acc[k]), cmul(b, out));
    }
  }
}

template <typename T>
T* DeviceAlloc(std::size_t count) {
  if (count == 0) return nullptr;
  void* ptr = nullptr;
  LSQ_CUDA_CHECK(cudaMalloc(&ptr, count * sizeof(T)));
  return static_cast<T*>(ptr);
}

template <typename T>
T* DeviceUpload(const std::vector<T>& host) {
  T* ptr = DeviceAlloc<T>(host.size());
  if (!host.empty()) {
    LSQ_CUDA_CHECK(
        cudaMemcpy(ptr, host.data(), host.size() * sizeof(T), cudaMemcpyHostToDevice));
  }
  return ptr;
}

void DeviceFree(void* ptr) {
  if (ptr != nullptr) {
    cudaFree(ptr);
  }
}

}  // namespace

void SynchronizeCuda() { LSQ_CUDA_CHECK(cudaDeviceSynchronize()); }

CudaVectorBlock::CudaVectorBlock(Index rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 0 || cols < 0) {
    throw std::invalid_argument("CudaVectorBlock: negative rows or cols");
  }
  dev_ = DeviceAlloc<Scalar>(static_cast<std::size_t>(rows) * cols);
  if (dev_ != nullptr) {
    LSQ_CUDA_CHECK(
        cudaMemset(dev_, 0, static_cast<std::size_t>(rows) * cols * sizeof(Scalar)));
  }
}

CudaVectorBlock::~CudaVectorBlock() { DeviceFree(dev_); }

void CudaVectorBlock::CopyFromHost(const Scalar* host) {
  LSQ_CUDA_CHECK(cudaMemcpy(dev_, host,
                            static_cast<std::size_t>(rows_) * cols_ * sizeof(Scalar),
                            cudaMemcpyHostToDevice));
}

void CudaVectorBlock::CopyToHost(Scalar* host) const {
  LSQ_CUDA_CHECK(cudaMemcpy(host, dev_,
                            static_cast<std::size_t>(rows_) * cols_ * sizeof(Scalar),
                            cudaMemcpyDeviceToHost));
}

CudaAtomBlockOperator::CudaAtomBlockOperator(
    std::vector<Index> atomOffset, std::vector<Index> pairOffset,
    std::vector<std::int32_t> nbr, std::vector<std::int32_t> blockId,
    std::vector<Index> dictOffset, std::vector<std::int32_t> dictRows,
    std::vector<std::int32_t> dictCols, std::vector<Scalar> dictValues,
    std::vector<Scalar> diagonal, bool hermitian)
    : hermitian_(hermitian) {
  const AtomBlockExtents ext =
      ValidateAtomBlockArrays(atomOffset, pairOffset, nbr, blockId, dictOffset, dictRows,
                              dictCols, dictValues, diagonal);
  dimension_ = ext.dimension;
  numAtoms_ = ext.numAtoms;
  numPairs_ = ext.numPairs;
  numBlocks_ = ext.numBlocks;
  maxBlockRows_ = ext.maxBlockRows;
  maxBlockCols_ = ext.maxBlockCols;
  dictScalars_ = static_cast<Index>(dictValues.size());

  atomOffset_ = DeviceUpload(atomOffset);
  pairOffset_ = DeviceUpload(pairOffset);
  nbr_ = DeviceUpload(nbr);
  blockId_ = DeviceUpload(blockId);
  dictOffset_ = DeviceUpload(dictOffset);
  dictCols_ = DeviceUpload(dictCols);
  dictValues_ = DeviceUpload(dictValues);
  dictRows_ = DeviceUpload(dictRows);
  diagonal_ = diagonal.empty() ? nullptr : DeviceUpload(diagonal);
}

CudaAtomBlockOperator::~CudaAtomBlockOperator() {
  DeviceFree(atomOffset_);
  DeviceFree(pairOffset_);
  DeviceFree(nbr_);
  DeviceFree(blockId_);
  DeviceFree(dictOffset_);
  DeviceFree(dictRows_);
  DeviceFree(dictCols_);
  DeviceFree(dictValues_);
  DeviceFree(diagonal_);
}

Index CudaAtomBlockOperator::Dimension() const { return dimension_; }

Index CudaAtomBlockOperator::NumAtoms() const { return numAtoms_; }

Index CudaAtomBlockOperator::NumPairs() const { return numPairs_; }

Index CudaAtomBlockOperator::DictionarySize() const { return numBlocks_; }

bool CudaAtomBlockOperator::IsHermitian() const { return hermitian_; }

void CudaAtomBlockOperator::MultiplyBlockDevice(Scalar a, const CudaVectorBlock& x, Scalar b,
                                                CudaVectorBlock& y) const {
  if (x.Rows() != dimension_ || y.Rows() != dimension_ || x.Cols() != y.Cols()) {
    throw std::invalid_argument(
        "CudaAtomBlockOperator::MultiplyBlockDevice: block dimension mismatch");
  }
  if (&x == &y) {
    throw std::invalid_argument(
        "CudaAtomBlockOperator::MultiplyBlockDevice: input and output blocks must not alias");
  }
  const int cols = x.Cols();
  if (dimension_ == 0 || cols == 0 || numAtoms_ == 0) return;

  const Index registerChunk = kWarpSize * kMaxAccumPerLane / maxBlockRows_;
  const Index sharedChunk = (48 * 1024) /
      ((kBlockThreads / kWarpSize) * std::max<Index>(1, maxBlockCols_) * sizeof(Scalar));
  const int chunk = static_cast<int>(std::min<Index>(
      kMaxStateChunk, std::min(registerChunk, sharedChunk)));
  if (chunk < 1) {
    throw std::invalid_argument(
        "CudaAtomBlockOperator::MultiplyBlockDevice: atom width exceeds kernel capacity");
  }
  for (int s0 = 0; s0 < cols; s0 += chunk) {
    const int sc = std::min(chunk, cols - s0);
    const std::size_t scratchBytes =
        static_cast<std::size_t>(kBlockThreads / kWarpSize) * maxBlockCols_ * sc * sizeof(Scalar);
    const std::size_t dictionaryBytes = static_cast<std::size_t>(dictScalars_) * sizeof(Scalar);
    const bool dictInShared = dictionaryBytes <= kMaxSharedDict &&
        dictionaryBytes + scratchBytes <= 48 * 1024;
    const std::size_t sharedBytes =
        (dictInShared ? static_cast<std::size_t>(dictScalars_) * sizeof(Scalar) : 0) +
        scratchBytes;
    const int warpsPerBlock = kBlockThreads / kWarpSize;
    const Index blocks = (numAtoms_ + warpsPerBlock - 1) / warpsPerBlock;
    const cuDoubleComplex av = make_cuDoubleComplex(a.real(), a.imag());
    const cuDoubleComplex bv = make_cuDoubleComplex(b.real(), b.imag());
    AtomBlockKernel<<<static_cast<unsigned>(blocks), kBlockThreads, sharedBytes>>>(
        av, reinterpret_cast<const cuDoubleComplex*>(x.DeviceData()) +
                static_cast<std::size_t>(s0) * dimension_,
        bv,
        reinterpret_cast<cuDoubleComplex*>(y.DeviceData()) +
            static_cast<std::size_t>(s0) * dimension_,
        sc, dimension_, numAtoms_, atomOffset_, pairOffset_, nbr_, blockId_,
        dictOffset_, dictCols_, reinterpret_cast<const cuDoubleComplex*>(dictValues_),
        reinterpret_cast<const cuDoubleComplex*>(diagonal_), dictInShared ? 1 : 0, dictScalars_,
        static_cast<int>(maxBlockCols_));
    LSQ_CUDA_CHECK(cudaGetLastError());
  }
}

void CudaAtomBlockOperator::Multiply(Scalar a, const Vector& x, Scalar b, Vector& y) const {
  const Index n = Dimension();
  if (static_cast<Index>(x.size()) != n || static_cast<Index>(y.size()) != n) {
    throw std::invalid_argument("CudaAtomBlockOperator::Multiply: vector dimension mismatch");
  }
  CudaVectorBlock dx(n, 1), dy(n, 1);
  dx.CopyFromHost(x.data());
  dy.CopyFromHost(y.data());
  MultiplyBlockDevice(a, dx, b, dy);
  dy.CopyToHost(y.data());
}

void CudaAtomBlockOperator::MultiplyBlock(Scalar a, const VectorBlock& x, Scalar b,
                                          VectorBlock& y) const {
  const Index n = Dimension();
  if (x.Rows() != n || y.Rows() != n || x.Cols() != y.Cols()) {
    throw std::invalid_argument("CudaAtomBlockOperator::MultiplyBlock: block dimension mismatch");
  }
  CudaVectorBlock dx(n, x.Cols()), dy(n, y.Cols());
  dx.CopyFromHost(x.Data());
  dy.CopyFromHost(y.Data());
  MultiplyBlockDevice(a, dx, b, dy);
  dy.CopyToHost(y.Data());
}

}  // namespace lsquant

#endif  // LSQUANT_WITH_CUDA
