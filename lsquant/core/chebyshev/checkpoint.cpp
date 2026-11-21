/**
 * @file checkpoint.cpp
 * @brief Implementation of the trace-state-level checkpoint format (version 1).
 */
#include "core/chebyshev/checkpoint.hpp"

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace lsquant {

namespace {

constexpr char kMagic[8] = {'L', 'S', 'Q', 'C', 'K', 'P', 'T', '\n'};
constexpr std::uint32_t kVersion = 1;

void WriteAll(FILE* fh, const void* data, std::size_t size, const char* what) {
  if (size > 0 && std::fwrite(data, 1, size, fh) != size) {
    throw std::runtime_error(std::string("WriteCheckpoint: failed to write ") + what);
  }
}

template <typename T>
void WriteScalar(FILE* fh, T value, const char* what) {
  WriteAll(fh, &value, sizeof(T), what);
}

void ReadAll(FILE* fh, void* data, std::size_t size, const char* what) {
  if (size > 0 && std::fread(data, 1, size, fh) != size) {
    throw std::runtime_error(std::string("TryLoadCheckpoint: truncated at ") + what);
  }
}

template <typename T>
T ReadScalar(FILE* fh, const char* what) {
  T value{};
  ReadAll(fh, &value, sizeof(T), what);
  return value;
}

std::string ReadString(FILE* fh, const char* what) {
  const std::uint32_t length = ReadScalar<std::uint32_t>(fh, what);
  if (length > (1u << 22)) {
    throw std::runtime_error(std::string("TryLoadCheckpoint: implausible string length at ") +
                             what);
  }
  std::string out(length, '\0');
  ReadAll(fh, out.data(), length, what);
  return out;
}

}  // namespace

void WriteCheckpoint(const std::string& path, const ChebyshevMoments& mu, Index statesDone) {
  const std::string tmp = path + ".tmp." + std::to_string(static_cast<long>(::getpid()));
  FILE* fh = std::fopen(tmp.c_str(), "wb");
  if (fh == nullptr) {
    throw std::runtime_error("WriteCheckpoint: cannot open " + tmp + ": " +
                             std::strerror(errno));
  }
  try {
    const auto& dims = mu.Dims();
    const auto& data = mu.Data();
    WriteAll(fh, kMagic, sizeof(kMagic), "magic");
    WriteScalar<std::uint32_t>(fh, kVersion, "version");
    WriteScalar<std::uint32_t>(fh, static_cast<std::uint32_t>(mu.Rank()), "rank");
    WriteScalar<std::uint32_t>(fh, static_cast<std::uint32_t>(dims.size()), "ndims");
    WriteScalar<std::int64_t>(fh, static_cast<std::int64_t>(data.size()), "element count");
    WriteScalar<std::int64_t>(fh, static_cast<std::int64_t>(mu.numberOfStates),
                              "number of states");
    WriteScalar<std::int64_t>(fh, static_cast<std::int64_t>(statesDone), "states done");
    WriteScalar<std::uint32_t>(fh, static_cast<std::uint32_t>(mu.stateLabel.size()),
                               "state label length");
    WriteAll(fh, mu.stateLabel.data(), mu.stateLabel.size(), "state label");
    WriteScalar<std::uint32_t>(fh, static_cast<std::uint32_t>(mu.word.size()), "word length");
    WriteAll(fh, mu.word.data(), mu.word.size(), "word");
    for (int d : dims) {
      WriteScalar<std::int64_t>(fh, static_cast<std::int64_t>(d), "dims");
    }
    WriteAll(fh, data.data(), data.size() * sizeof(Scalar), "tensor payload");
    if (std::fflush(fh) != 0) {
      throw std::runtime_error("WriteCheckpoint: flush failed on " + tmp);
    }
    if (::fsync(::fileno(fh)) != 0) {
      throw std::runtime_error("WriteCheckpoint: fsync failed on " + tmp);
    }
    if (std::fclose(fh) != 0) {
      throw std::runtime_error("WriteCheckpoint: close failed on " + tmp);
    }
  } catch (...) {
    std::fclose(fh);
    std::remove(tmp.c_str());
    throw;
  }
  if (std::rename(tmp.c_str(), path.c_str()) != 0) {
    std::remove(tmp.c_str());
    throw std::runtime_error("WriteCheckpoint: rename onto " + path + " failed: " +
                             std::strerror(errno));
  }
}

bool TryLoadCheckpoint(const std::string& path, ChebyshevMoments& mu, Index& statesDone) {
  FILE* fh = std::fopen(path.c_str(), "rb");
  if (fh == nullptr) {
    if (errno == ENOENT) {
      return false;
    }
    throw std::runtime_error("TryLoadCheckpoint: cannot open " + path + ": " +
                             std::strerror(errno));
  }
  auto fail = [&](const std::string& what) -> void {
    std::fclose(fh);
    throw std::runtime_error("TryLoadCheckpoint: " + what + " in " + path);
  };
  char magic[sizeof(kMagic)];
  ReadAll(fh, magic, sizeof(magic), "magic");
  if (std::memcmp(magic, kMagic, sizeof(kMagic)) != 0) {
    fail("bad magic");
  }
  if (ReadScalar<std::uint32_t>(fh, "version") != kVersion) {
    fail("unsupported format version");
  }
  const std::uint32_t rank = ReadScalar<std::uint32_t>(fh, "rank");
  const std::uint32_t ndims = ReadScalar<std::uint32_t>(fh, "ndims");
  const std::int64_t elements = ReadScalar<std::int64_t>(fh, "element count");
  const std::int64_t totalStates = ReadScalar<std::int64_t>(fh, "number of states");
  const std::int64_t done = ReadScalar<std::int64_t>(fh, "states done");
  const std::string label = ReadString(fh, "state label");
  const std::string word = ReadString(fh, "word");
  std::vector<int> dims(ndims);
  for (std::uint32_t k = 0; k < ndims; ++k) {
    const std::int64_t d = ReadScalar<std::int64_t>(fh, "dims");
    if (d <= 0 || d > (std::int64_t{1} << 40)) {
      fail("implausible tensor dimension");
    }
    dims[k] = static_cast<int>(d);
  }
  if (rank != static_cast<std::uint32_t>(mu.Rank()) || dims != mu.Dims() ||
      elements != static_cast<std::int64_t>(mu.Data().size()) ||
      totalStates != static_cast<std::int64_t>(mu.numberOfStates) || label != mu.stateLabel ||
      word != mu.word) {
    fail("checkpoint does not match this run (rank/dims/states/label/word)");
  }
  if (done < 0 || done > totalStates) {
    fail("states done out of range");
  }
  auto& data = mu.Data();
  ReadAll(fh, data.data(), data.size() * sizeof(Scalar), "tensor payload");
  std::fclose(fh);
  statesDone = static_cast<Index>(done);
  return true;
}

}  // namespace lsquant
