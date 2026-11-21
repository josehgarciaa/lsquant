#!/usr/bin/env bash
# The run-file workflow of the lsquant CLI on the Haldane model (examples/run_files/README.md).
# usage, from the repository root: examples/run_files/run.sh [out_dir]
set -euo pipefail
OUT="${1:-out/run_files}"
RF=examples/run_files/haldane_hall.json
echo "== 1. what does the word compute?"
python3 -m lsquant words describe "Tr[VX delta(E-H) VY dG-(E)/dE]" | head -20
echo "== 2. the run file, validated and resolved"
python3 -m lsquant inspect "$RF" | head -30
echo "== 3. run: moments and curves"
python3 -m lsquant run "$RF" --out "$OUT/haldane"
echo "== 4. reconstruct sigma_xy again at 150 meV from the stored VX-VY moments (no new moments)"
python3 -m lsquant reconstruct "$OUT/haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D" \
    --formula bastin --broadening-meV 150 --out "$OUT/haldane_150meV"
echo "== 5. Chern number from the plateau of sigma_xy at 10 and 150 meV"
python3 examples/run_files/chern.py "$OUT/haldane" "$OUT/haldane_150meV"
