#ifndef LSQUANT_LEGACY_MOMENT_INPUT_HPP
#define LSQUANT_LEGACY_MOMENT_INPUT_HPP

#include <fstream>
#include <iostream>

namespace chebyshev {
/**
 * @brief Reject inaccessible moment inputs before legacy constructors run.
 *
 * Legacy readers assert that opening succeeded; release builds can instead use
 * uninitialised header values after an open failure. A separate, constant-space
 * probe lets command-line callers return status 1 before allocating moments or
 * creating output. Reading one byte also detects paths that open but cannot be
 * read (such as directories). The production reader and its arithmetic stay intact.
 *
 * @par Pitfall
 * This is an accessibility check, not format validation. The file must remain
 * unchanged between this probe and the constructor's subsequent open.
 */
inline bool LegacyMomentInputReadable(const char* filename) {
    std::ifstream input(filename);
    if (input.is_open()) input.peek();
    if (!input.is_open() || input.bad()) {
        std::cerr << "Cannot read moment file: " << filename << '\n';
        return false;
    }
    return true;
}
}  // namespace chebyshev
#endif
