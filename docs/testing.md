\page testing Testing: golden masters and the regression gate

LSQUANT is tested against **golden masters**: output files written once by a known-good
build and committed under `goldens/` as permanent references. A change is accepted only
if the programs it touches still reproduce those files. This page explains what the
goldens contain, how to run the gate that compares against them, how to read its output,
and how to run the unit tests of the gate itself.

## What is in `goldens/stage0/`

Every file was produced by the legacy LSQUANT build (the `inline_*` drivers built from
`src/`) in one session, with `KPM_SEED=12345`. `goldens/stage0/PROVENANCE.md` records
the machine, compilers, library versions, repository commit and the exact command for
every file.

| Golden | Directory | Model | What it pins |
|---|---|---|---|
| G1 | `chain1d/` | 1D chain, N = 2048, window [-2, 2] | spectral (DOS) Chebyshev moments, exact trace, M = 64 |
| G2 | `chain1d/` | same | velocity-velocity moments ⟨VX δ VX δ⟩, M = 64 × 64 |
| G3 | `chain1d/` | same | reconstructed DOS curve (`DOS_reconstructed.dat`, Jackson kernel) |
| G4 | `graphene/` | graphene 8 × 8, 17 operators | Kubo-Bastin moments and curves, random-phase trace; FFT-grid, sea and surface variants |
| G5 | `haldane/` | Haldane model, dim 512, bounds [-4.5, 4.5] | VX-VY and VX-VX moments, M = 128; Kubo-Bastin (total, sea, surface) and Kubo-Greenwood curves at 10 meV |
| tutorials | `tutorials/01..06/` | the six `examples/` tutorials | every data file each tutorial script reads to draw its figure |

The `operators/` subdirectories hold the input Hamiltonians and operators in LSQUANT's CSR
text format, so the gate reads them directly and runs without a model generator.

## Running the gate

The gate is `scripts/check_goldens.sh <stage> [build_dir]`. Stage 0 checks the legacy
pipeline: it reruns G1 and G2 at 16 threads, G5 at 16 threads, and G1 and G2 again at 1
and 64 threads (the moments must remain identical across thread counts). Each run happens in a
temporary directory, leaving `goldens/` unchanged.

Example session, from a clean checkout:

```bash
# 1. build the legacy drivers (needs C++17, CMake 3.21+, Boost headers, Eigen, FFTW, GSL; see README.md)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLSQUANT_BUILD_INPUT_TOOLS=OFF
cmake --build build -j 16

# 2. run the stage-0 gate against that build
scripts/check_goldens.sh 0 build
```

Expected output (abridged):

```
== stage 0: legacy pipelines vs goldens/stage0 (KPM_SEED=12345, build /path/to/build)
IDENTICAL [T=16] chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D (G1)
IDENTICAL [T=16] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
IDENTICAL [T=16] haldane/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D (G5)
...
== 16 files compared; failures: 0
GATE PASS (stage 0)
```

The exit code is 0 only when every file is byte-identical, 1 when at least one differs and
2 for a usage error (unknown stage, missing build directory, missing driver).

## Running on another machine: tolerance mode

The goldens were written on one machine with one compiler (see `PROVENANCE.md`). On a
different compiler, C library or Eigen version the floating-point operations can be
reordered, and a few printed values change in their last digit while the physics stays the same.
For that case the gate has a tolerance mode:

```bash
GOLDEN_TOL=1 scripts/check_goldens.sh 0 build
```

The gate then classifies each file that differs byte for byte:

- **moment files** (`*.chebmom*`, 15 significant digits): the header lines must be
  identical and max_k |new_k − gold_k| ≤ 1e-12 · max_k |gold_k|;
- **curves** (`*.dat`, 6 significant digits): the energy column must be identical as text
  and every other value within one unit of its last printed digit, or within
  1e-12 × (largest value of that curve), whichever is larger. The second bound is a noise
  floor: a value a trillion times smaller than the curve maximum is numerically zero (the
  moments are only required to 1e-12), and its digits are rounding noise;
- **other text files**: non-numeric text identical, every number within one unit of its
  last printed digit.

A file that meets its rule is reported as `TOLPASS` and counts as passing; the summary line
gives the number of byte-identical and tolerance-passing files separately. Example from a
machine with g++ 9.4 (the goldens were written with g++ 11.4): every moment file is
identical, and two Haldane Kubo-Bastin curves differ by one unit in the sixth significant
digit, e.g. `5.91198e-06` against `5.91197e-06`:

```
TOLPASS   [T=16] haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5)
          tolerance: curve: 3 values differ, 0 beyond one unit, max 0.0929 units of the last
          printed digit ... new=5.91198e-06 golden=5.91197e-06, max|d| = 1.000e-11,
          max|d|/max|gold| = 9.285e-14 (informative)
TOLPASS   [T=16] haldane/KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat (G5)
...
== 12 files compared; failures: 0
== tolerance mode (GOLDEN_TOL=1): 10 byte-identical, 2 tolerance-PASS
GATE PASS (stage 0)
```

With the default strict byte comparison, the same run reports those two files as `DIFFERS`.
Use tolerance mode to check a build on a new machine; code changes require strict comparison.

## Reading a failure

When a file differs, the gate prints the first differing bytes and one line from
`scripts/golden_diff.py`, which locates the difference in the units of the file:

```
DIFFERS   [T=64] chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D (G2)
          cmp: ... differ: byte 812, line 40
          first-diff line 40 (m,n)=(0,37) re, new=0.0123 golden=0.0124;
          max abs diff 1.000e-04 at ...; max rel diff 8.065e-03 at ...; DEFECT level (>= 1e-12)
```

`golden_diff.py` understands the moment layouts (`.chebmom1D`: index m; `.chebmom2D`:
index (m, n); `.chebmomTD`: index (m, t)) and treats any other file as a table of
numbers (row, column). It classifies the largest relative difference
|a − b| / max(|a|, |b|) as:

- **identical values (formatting only)**: every number is the same, only its text differs;
- **rounding level (< 1e-12)**: a floating-point reordering (for example a different
  reduction order on another CPU). A maintainer must review it before it is accepted;
- **DEFECT level (≥ 1e-12)**: a real change of the result. The change is rejected.

It can also be run by hand on any two files:

```bash
python3 scripts/golden_diff.py my_output.chebmom2D goldens/stage0/chain1d/NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D
```

## Regenerating the tutorial goldens

`figures/stage_0.py` reruns every command block of the six tutorial READMEs with the
legacy build, captures each data file the tutorial scripts read, and compares them with
`goldens/stage0/tutorials/`. It also redraws the tutorial figures with a caption that
records model, size, number of moments and seed:

```bash
python3 figures/stage_0.py --build build --out figures/out/stage_0
```

It writes under `--out` and its work directory, leaving `examples/` and `goldens/` unchanged.

## Unit tests of the gate

The gate and the localiser have their own unit tests under `tests/unit/`. The gate tests
replace the six legacy drivers by small fake executables that copy (or corrupt, or omit)
the golden files, so they run without a build, in a few seconds:

```bash
python3 -m pytest tests/unit -v
```

They check that the gate passes on identical outputs, fails and localises a corrupted
moment, reports a missing output file, stops with exit code 2 when a driver is
missing, and leaves `goldens/` unchanged; and, in tolerance mode, that a curve off by one unit
of its last digit and moments off by 1e-13 pass while five units or 1e-9 fail.

## Rules

- Do not edit, regenerate or delete an existing golden. Add a new golden in a new
  directory, with its provenance, in the change that needs it.
- The gate fixes seeds (`KPM_SEED=12345`) and tolerances; the caller cannot change them.
- A pull request that changes numerical code states which gate stages it ran and pastes
  the final `GATE PASS` line.
