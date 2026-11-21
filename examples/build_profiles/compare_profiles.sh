#!/usr/bin/env bash
# Build the golden and native profiles with the CMake presets, run the same run file with
# both, and compare the outputs and the run times (examples/build_profiles/README.md).
# usage, from the repository root:
#   examples/build_profiles/compare_profiles.sh [runfile] [out_dir]
set -euo pipefail
export LC_NUMERIC=C      # printf below uses a decimal point
RF="${1:-runs/goldens/g5.json}"
OUT="${2:-out/build_profiles}"
mkdir -p "$OUT"
for p in golden native; do
  if [ ! -x "build/$p/lsquant_core/stage1_moments" ]; then
    cmake --preset "$p" > "$OUT/configure_$p.log"
    cmake --build --preset "$p" > "$OUT/build_$p.log"
  fi
  s=$(date +%s.%N)
  PYTHONPATH="$PWD/build/$p/lsquant_core:$PWD/python" python3 -m lsquant run "$RF" --out "$OUT/$p" > "$OUT/run_$p.log"
  e=$(date +%s.%N)
  printf "%-7s run: %6.2f s\n" "$p" "$(echo "$e - $s" | bc)"
done
python3 examples/build_profiles/compare_outputs.py "$OUT/golden" "$OUT/native"
