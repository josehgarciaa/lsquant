/**
 * @file checkpoint.hpp
 * @brief Trace-state-level checkpoint and resume for MomentEngine runs.
 *
 * A `MomentEngine` run accumulates the moment tensor state by state in a fixed order, so the
 * only state that must survive a crash is the partial tensor plus the number of finished
 * states: the state generator reconstructs the exact stream position from its seed
 * (`StateGenerator::Skip`). The file format is versioned (`LSQCKPT\\n`, version 1) and holds
 * the rank, dims, total state count, states done, state label, plan word and the raw
 * complex128 tensor payload. Writes are atomic (temporary file, `fsync`, `rename`). A resumed
 * run is bit-identical to an uninterrupted one.
 */
#ifndef LSQUANT_CORE_CHEBYSHEV_CHECKPOINT_HPP_
#define LSQUANT_CORE_CHEBYSHEV_CHECKPOINT_HPP_

#include <string>

#include "core/algebra/types.hpp"
#include "core/chebyshev/chebyshev_moments.hpp"

namespace lsquant {

/**
 * @brief Persist the partial moment tensor and the finished-state count, atomically.
 *
 * Serialises to `<path>.tmp.<pid>` and renames onto `path`, so a concurrent or crashed writer
 * never leaves a truncated file at `path`. Cost is one write of the tensor payload (a few MB
 * for rank 2) per call.
 *
 * @param path       Destination checkpoint path.
 * @param mu         Partial tensor (raw, pre-weighting accumulation).
 * @param statesDone Number of fully accumulated trace states `0 <= statesDone`.
 *
 * @warning Throws `std::runtime_error` on any I/O failure.
 */
void WriteCheckpoint(const std::string& path, const ChebyshevMoments& mu, Index statesDone);

/**
 * @brief Load a checkpoint into `mu` if a matching file exists at `path`.
 *
 * A missing file returns `false` (fresh run). An existing file must match `mu`'s rank, dims,
 * `numberOfStates`, `stateLabel` and `word` exactly and hold `0 <= statesDone <=
 * numberOfStates`; any mismatch, truncation or version difference throws rather than risking
 * a silent restart on the wrong problem. On success `mu`'s tensor data is replaced by the
 * stored partial accumulation and `statesDone` is set.
 *
 * @param path       Checkpoint path.
 * @param mu         Tensor whose metadata define the expected run; data overwritten on match.
 * @param statesDone Output: number of finished states recorded in the checkpoint.
 * @return           `true` if a matching checkpoint was loaded, `false` if none exists.
 *
 * @warning Throws `std::runtime_error` on any mismatch or corruption.
 */
bool TryLoadCheckpoint(const std::string& path, ChebyshevMoments& mu, Index& statesDone);

}  // namespace lsquant

#endif  // LSQUANT_CORE_CHEBYSHEV_CHECKPOINT_HPP_
