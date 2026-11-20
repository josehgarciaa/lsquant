#!/usr/bin/env bash
# Golden-master gate. Runs the stage's programs on the committed golden inputs in a
# temporary directory and compares every produced file byte for byte with
# goldens/stage0/. It never writes into goldens/ and never modifies a golden.
# See docs/testing.md for the procedure, the meaning of each golden (G1..G5) and an
# example session.
#
# usage: scripts/check_goldens.sh <stage> [build_dir]
#
#   stage 0  legacy build (build_dir holds the inline_* drivers): reruns the G1, G2
#            (chain1d, exact trace, M = 64) and G5 (Haldane, exact trace, M = 128)
#            pipelines at OMP_NUM_THREADS=16, then G1 and G2 again at 1 and 64 threads.
#
#   Later stages (the new core and its Python layer) add their own case to the
#   dispatch at the end of this file when they are merged.
#
# Tolerance mode (GOLDEN_TOL=1; final validation on a machine that is not byte-identical to the
# machine where the goldens were written, e.g. one with gcc 9.4, Eigen 3.4, glibc 2.31).
# Every compare() that is not byte-identical is classified by a helper (tolcmp.py, written into
# $TMP) and reported as TOLPASS (counted as passing) or DIFFERS (counted as failing), with the
# scripts/golden_diff.py numbers either way:
#   moment files (*.chebmom*, 15 significant digits): the two header lines identical as text,
#     same layout (lines, tokens per line), and max_k |new_k - gold_k| <= 1e-12 * max_k |gold_k|
#     over the complex moments;
#   text curves (*.dat, 6 significant digits): same layout, the energy column (first column)
#     identical as text, every other value within one unit of the last printed significant digit
#     (unit = 10^(floor(log10|x|) - P + 1), P = the largest number of significant digits printed
#     in the two files, the larger unit of the two values; absolute, never relative to the value),
#     and never smaller than the noise floor 1e-12 * max |golden value| of that curve;
#   every other file (JSON, TSV, CSR, index files): non-numeric text identical, every number within
#     one unit of the last printed significant digit as above.
# Totals report byte-identical and tolerance-PASS counts separately.
#
# Output: one line per compared file, IDENTICAL or DIFFERS; for DIFFERS the first
# differing bytes (cmp -l | head -5) and the numeric localisation of
# scripts/golden_diff.py (first differing line, element index, both values, max
# relative difference over the file). Exit code 0 only if every file is IDENTICAL
# (or, with GOLDEN_TOL=1, IDENTICAL or TOLPASS).
#
# KPM_SEED=12345 always. Seeds, tolerances and goldens are never changed here.
set -uo pipefail

usage() { echo "usage: scripts/check_goldens.sh <stage> [build_dir]"; exit 2; }
[ $# -ge 1 ] || usage
STAGE="$1"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="$(cd "${2:-$REPO/build}" 2>/dev/null && pwd)" || { echo "ERROR: build dir '${2:-build}' not found"; exit 2; }
GOLD="$REPO/goldens/stage0"
export KPM_SEED=12345
[ -d "$GOLD" ] || { echo "ERROR: $GOLD missing"; exit 2; }

TMP="$(mktemp -d "${TMPDIR:-/tmp}/check_goldens.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
NFAIL=0
NCMP=0
GOLDEN_TOL="${GOLDEN_TOL:-0}"
NBYTE=0; NTOLP=0; CMP_LAST=""; FORCE_BYTE=0
TOLCMP="$TMP/tolcmp.py"
cat >"$TOLCMP" <<'TOLEOF'
# tolerance classifier of check_goldens.sh (GOLDEN_TOL=1); prints "TOLPASS ..." or "DIFFERS ...", exit 0/1
import math, os, re, sys
new_p, gold_p = sys.argv[1], sys.argv[2]
NUM = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|[-+]?(?:nan|inf)', re.I)
def toks(line):
    out, pos = [], 0
    for m in NUM.finditer(line):
        if m.start() > pos: out.append(("t", line[pos:m.start()]))
        out.append(("n", m.group(0))); pos = m.end()
    if pos < len(line): out.append(("t", line[pos:]))
    return out
def sig(tok):
    m = re.match(r'[-+]?(\d*)\.?(\d*)', tok)
    d = (m.group(1) + m.group(2)).lstrip("0") if m else ""
    return len(d)
def done(ok, msg):
    print(("TOLPASS " if ok else "DIFFERS ") + msg); sys.exit(0 if ok else 1)
new = open(new_p, errors="replace").read().split("\n")
gold = open(gold_p, errors="replace").read().split("\n")
if len(new) != len(gold):
    done(False, "layout: %d lines new vs %d golden" % (len(new), len(gold)))
ext = os.path.splitext(gold_p)[1]
if ext.startswith(".chebmom"):
    if new[:2] != gold[:2]:
        done(False, "moment header differs: new %r golden %r" % (new[:2], gold[:2]))
    a, b = [], []
    for i, (x, y) in enumerate(zip(new[2:], gold[2:])):
        tx, ty = x.split(), y.split()
        if len(tx) != len(ty):
            done(False, "layout: line %d has %d vs %d tokens" % (i + 3, len(tx), len(ty)))
        if not tx: continue
        a.append(complex(float(tx[0]), float(tx[1]) if len(tx) > 1 else 0.0))
        b.append(complex(float(ty[0]), float(ty[1]) if len(ty) > 1 else 0.0))
    g = max((abs(v) for v in b), default=0.0)
    d = [abs(x - y) for x, y in zip(a, b)]
    k = max(range(len(d)), key=lambda i: d[i]) if d else 0
    r = (d[k] / g) if g > 0 else (0.0 if not d or d[k] == 0 else math.inf)
    ok = all(math.isfinite(abs(v)) for v in a) and r <= 1e-12
    done(ok, "moments: max|d|/max|gold| = %.3e (<= 1e-12 required), max|d| = %.3e at moment %d, max|gold| = %.3e, header identical"
         % (r, d[k] if d else 0.0, k, g))
P = 1
for y in gold + new:
    for kind, t in toks(y):
        if kind == "n": P = max(P, sig(t))
def unit(x):
    if x == 0 or not math.isfinite(x): return 0.0
    return 10.0 ** (math.floor(math.log10(abs(x))) - P + 1)
curve = ext == ".dat"
# largest |value| of the golden (energy column of a curve excluded): sets the noise floor
gmax = 0.0
for y in gold:
    col = -1
    for kind, t in toks(y):
        if kind != "n": continue
        col += 1
        if curve and col == 0: continue
        try: v = abs(float(t))
        except ValueError: continue
        if math.isfinite(v): gmax = max(gmax, v)
# noise floor of a curve: a value at 1e-12 of the curve's largest value is numerically zero
# (the moments themselves are only required to 1e-12), so the unit is never smaller than that
floor = 1e-12 * gmax if curve else 0.0
worst, nd, nbad, where, maxabs = 0.0, 0, 0, "-", 0.0
for i, (x, y) in enumerate(zip(new, gold)):
    if x == y: continue
    tx, ty = toks(x), toks(y)
    if len(tx) != len(ty) or any(p[0] != q[0] for p, q in zip(tx, ty)):
        done(False, "layout: line %d tokenises differently: %r vs %r" % (i + 1, x[:120], y[:120]))
    col = -1
    for (kx, sx), (ky, sy) in zip(tx, ty):
        if kx == "n": col += 1
        if sx == sy: continue
        if kx == "t":
            done(False, "text differs at line %d: %r vs %r" % (i + 1, sx[:80], sy[:80]))
        if curve and col == 0:
            done(False, "energy column differs at line %d: %s vs %s" % (i + 1, sx, sy))
        a, b = float(sx), float(sy)
        if not (math.isfinite(a) and math.isfinite(b)):
            done(False, "non-finite value differs at line %d col %d: %s vs %s" % (i + 1, col + 1, sx, sy))
        u = max(unit(a), unit(b), floor)
        dd = abs(a - b); nd += 1; maxabs = max(maxabs, dd)
        q = dd / u if u > 0 else math.inf
        nbad += q > 1.0 + 1e-9
        if q > worst: worst, where = q, "line %d col %d new=%s golden=%s" % (i + 1, col + 1, sx, sy)
ok = worst <= 1.0 + 1e-9
done(ok, "%s: %d values differ, %d beyond one unit, max %.3g units of the last printed digit (P = %d significant digits; <= 1 required) at %s, max|d| = %.3e, max|d|/max|gold| = %.3e (informative)"
     % ("curve" if curve else "text", nd, nbad, worst, P, where, maxabs, maxabs / gmax if gmax > 0 else 0.0))
TOLEOF

# compare <produced file> <golden file> <tag>
compare() {
  local new="$1" gold="$2" tag="$3"
  NCMP=$((NCMP + 1))
  if [ ! -f "$new" ]; then
    echo "DIFFERS   $tag  (not produced: $(basename "$new"))"; NFAIL=$((NFAIL + 1)); CMP_LAST=DIFFERS; return
  fi
  if [ ! -f "$gold" ]; then
    echo "DIFFERS   $tag  (no golden: ${gold#$REPO/})"; NFAIL=$((NFAIL + 1)); CMP_LAST=DIFFERS; return
  fi
  if cmp -s "$new" "$gold"; then
    echo "IDENTICAL $tag"; NBYTE=$((NBYTE + 1)); CMP_LAST=IDENTICAL
  elif [ "$GOLDEN_TOL" = 1 ] && [ "$FORCE_BYTE" = 0 ]; then
    local t; t="$(python3 "$TOLCMP" "$new" "$gold" 2>&1)"
    if [ "${t%% *}" = TOLPASS ]; then
      echo "TOLPASS   $tag"; NTOLP=$((NTOLP + 1)); CMP_LAST=TOLPASS
    else
      echo "DIFFERS   $tag"; NFAIL=$((NFAIL + 1)); CMP_LAST=DIFFERS
    fi
    echo "          tolerance: ${t#* }"
    echo "          golden_diff.py: $(python3 "$REPO/scripts/golden_diff.py" "$new" "$gold" 2>&1)"
  else
    NFAIL=$((NFAIL + 1)); CMP_LAST=DIFFERS
    echo "DIFFERS   $tag"
    echo "          cmp: $(cmp "$new" "$gold" 2>&1 | head -1)"
    cmp -l "$new" "$gold" 2>/dev/null | head -5 | sed 's/^/          cmp -l: /'
    echo "          $(python3 "$REPO/scripts/golden_diff.py" "$new" "$gold" 2>&1)"
  fi
}

# exact-trace state file over the full basis of a CSR (legacy "local" state set)
write_exact() { local n; n=$(head -1 "$1" | awk '{print $1}'); { echo local; echo "$n"; seq 0 $((n - 1)); } > "$2"; }

# ---------------------------------------------------------------- stage 0 (legacy)
# G1 + G2: the chain1d part of scripts/regenerate_goldens.sh, on the golden operators.
legacy_chain1d() {   # <threads> <workdir>
  local T="$1" W="$2"
  mkdir -p "$W/operators"
  cp "$GOLD/chain1d/operators/"*.CSR "$W/operators/"
  cp "$GOLD/chain1d/BOUNDS" "$W/BOUNDS"
  write_exact "$W/operators/chain1d.HAM.CSR" "$W/exact"
  ( cd "$W"
    OMP_NUM_THREADS=$T "$BUILD/inline_compute-kpm-nonEqOp"   chain1d VX VX 64 exact >nonEq.log 2>&1
    OMP_NUM_THREADS=$T "$BUILD/inline_compute-kpm-spectralOp" chain1d 1 64 exact   >spectral.log 2>&1 )
  compare "$W/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" \
          "$GOLD/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D" \
          "[T=$T] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)"
  compare "$W/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
          "$GOLD/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D" \
          "[T=$T] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)"
}

# G5: Haldane, exact trace, M = 128, bounds [-4.5, 4.5] (test/hall_haldane.sh), 10 meV.
legacy_haldane() {   # <threads> <workdir>
  local T="$1" W="$2" f
  mkdir -p "$W/operators"
  cp "$GOLD/haldane/operators/"haldane.{HAM,VX,VY}.CSR "$W/operators/"
  cp "$GOLD/haldane/BOUNDS" "$W/BOUNDS"
  write_exact "$W/operators/haldane.HAM.CSR" "$W/exact"
  ( cd "$W"
    export OMP_NUM_THREADS=$T
    "$BUILD/inline_compute-kpm-nonEqOp" haldane VX VY 128 exact >xy.log 2>&1
    "$BUILD/inline_compute-kpm-nonEqOp" haldane VX VX 128 exact >xx.log 2>&1
    XY=NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D
    XX=NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D
    "$BUILD/inline_kuboBastinFromChebmom"     "$XY" 10 >bastin.log 2>&1
    "$BUILD/inline_kuboBastinSeaFromChebmom"  "$XY" 10 >sea.log 2>&1
    "$BUILD/inline_kuboBastinSurfFromChebmom" "$XY" 10 >surf.log 2>&1
    "$BUILD/inline_kuboGreenwoodFromChebmom"  "$XX" 10 >greenwood.log 2>&1 )
  for f in NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D \
           NonEqOpVX-VXhaldaneKPM_M128x128_stateexact.chebmom2D; do
    compare "$W/$f" "$GOLD/haldane/$f" "[T=$T] haldane/$f (G5)"
  done
  for f in "$GOLD"/haldane/*.dat; do
    compare "$W/$(basename "$f")" "$f" "[T=$T] haldane/$(basename "$f") (G5)"
  done
}

stage0() {
  for b in inline_compute-kpm-nonEqOp inline_compute-kpm-spectralOp inline_kuboBastinFromChebmom \
           inline_kuboBastinSeaFromChebmom inline_kuboBastinSurfFromChebmom inline_kuboGreenwoodFromChebmom; do
    [ -x "$BUILD/$b" ] || { echo "ERROR: legacy driver $BUILD/$b missing (build with the legacy tree)"; exit 2; }
  done
  echo "== stage 0: legacy pipelines vs goldens/stage0 (KPM_SEED=$KPM_SEED, build $BUILD)"
  legacy_chain1d 16 "$TMP/g12_T16"
  legacy_haldane 16 "$TMP/g5_T16"
  legacy_chain1d 1  "$TMP/g12_T1"
  legacy_chain1d 64 "$TMP/g12_T64"
}

case "$STAGE" in
  0) stage0 ;;
  *) echo "ERROR: stage $STAGE has no gate defined yet"; exit 2 ;;
esac

echo "== $NCMP files compared; failures: $NFAIL"
[ "$GOLDEN_TOL" = 1 ] && echo "== tolerance mode (GOLDEN_TOL=1): $NBYTE byte-identical, $NTOLP tolerance-PASS"
[ "$NFAIL" -eq 0 ] && { echo "GATE PASS (stage $STAGE)"; exit 0; }
echo "GATE FAIL (stage $STAGE)"; exit 1
