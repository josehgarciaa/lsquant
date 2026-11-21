#!/usr/bin/env bash
# Chebyshev moments of a 1D chain with the stage-1 core, then the reconstructed DOS.
#
# usage: examples/core_engine/run.sh [build_dir]       (default: build)
#
# 1. stage1_moments evaluates the "dos" plan (Tr[T_m(H~)] over the full basis, M = 64)
#    on the window [-2.5, 2.5], and the "kubo" plan (Tr[VX T_m(H~) VX T_n(H~)], M = 64 x 64) on the 1D chain of
#    the golden set (N = 2048, hopping -1, periodic), with an exact trace;
# 2. dos_from_moments.py rebuilds the DOS from the .chebmom file with the Jackson
#    kernel and compares it with the analytic 1D DOS 1 / (pi sqrt(4 - E^2)).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
BUILD="$(cd "${1:-$REPO/build}" && pwd)"
EXE="$BUILD/lsquant_core/stage1_moments"
OUT="${OUT:-$HERE/out}"
mkdir -p "$OUT"
export KPM_SEED=12345
OPS="$REPO/goldens/stage0/chain1d/operators"
# kubo on the golden window [-2, 2]: the same tensor as G2
"$EXE" --operators "$OPS" --label chain1d --window -2 2 --trace exact --states 2048 \
       --moments 64 --word kubo --out "$OUT"
# dos on a window wider than the band. On exactly [-2, 2] the chain DOS
# 1/(pi sqrt(1 - x^2)) is the Chebyshev weight itself: every moment m >= 1 vanishes and
# the reconstruction is exact by construction, which tests nothing.
"$EXE" --operators "$OPS" --label chain1d --window -2.5 2.5 --trace exact --states 2048 \
       --moments 64 --word dos --out "$OUT"
python3 "$HERE/dos_from_moments.py" "$OUT/SpectralOp1chain1dKPM_M64_stateexact.chebmom" \
        --out "$OUT/dos_chain1d.dat"
