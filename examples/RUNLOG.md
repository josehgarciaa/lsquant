# Run log of the examples and the quickstart (2026-11-21)

All runs were on a shared 96-core Linux node that was oversubscribed
throughout (load 150–175), under `nice -n 19`, with the `golden` build preset of
the final code of this pull-request series (`_core.build_profile == "golden"`).
Each pull request carries the rows of the examples it adds; the numbers were
measured on the final code. Nothing was
compiled for these runs. Python 3.10.22 with numpy 2.2.6, scipy 1.15.2,
sisl 0.16.4 and mcp 2.3.0. The run directory was a copy of the
`examples/`, `python/` and the needed `goldens/stage0/{haldane,chain1d,nonorth}`.
The log of every command was kept with one line per command (command, threads,
wall time, exit code); the outputs were not copied into the example folders
(examples 01–08 do not commit outputs either).

Every command below was run from the repository root with
`PYTHONPATH=<build>/lsquant_core:<repo>/python`. `python3` is the
interpreter above. Wall times are whole-command times including Python start-up.

## Thread count: OMP_NUM_THREADS=1 on a busy node

The first round used `OMP_NUM_THREADS=4`. On this oversubscribed node (with
`OMP_PROC_BIND=spread`, `OMP_PLACES=cores` set in the environment) small
problems became 10–30 times slower than with one thread, because the parallel
regions wait for descheduled threads. Same command, same node, minutes apart:

| command | OMP 4 | OMP 1 |
|---|---|---|
| `lsquant run examples/11_nonorthogonal_toy/run_orthogonal.json` (N = 128, M = 256, exact) | 145.6 s | 5.6 s (1.8 s in the final run) |
| `lsquant run examples/01_dos_1d_chain/run_dos.json` | 357.7 s (first run), 2892 s (second run, killed by me after 48 min while 3 other OMP-4 jobs ran) | 148.7 s |

All final runs used `OMP_NUM_THREADS=1`. The OMP-4 runs that finished are kept
in the table where they are the only run of that example (examples 09, 10, 14)
and are marked "OMP 4".

## Results

| # | example | command | threads | wall | exit | key result / check |
|---|---|---|---|---|---|---|
| 01 | DOS 1D chain | `python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out out/ex01` | 4 | 357.7 s | 0 | 3 chebmom1D, eta{10,20,24,80}/ curves |
| 09 | inspect dir | `python3 -m lsquant inspect out/ex01` | 4 | 2 s | 0 | JSON, `kind: directory`, resolved run file + 6 moment files |
| 09 | inspect moments | `python3 -m lsquant inspect out/ex01/SpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256.chebmom1D` | 4 | 2 s | 0 | rank 1, dims [512], window [-2, 2]; legacy file: `word ""`, `number_of_states 0` |
| 09 | inspect run file | `python3 -m lsquant inspect examples/01_dos_1d_chain/run_dos.json` | 4 | 1 s | 0 | `kind: runfile`, defaults resolved |
| 09 | reconstruct (named) | `python3 -m lsquant reconstruct <N256 .chebmom1D> --formula dos --broadening-meV 60 --out out/ex09/eta60 --json` | 4 | 1 s | 0 | JSON summary, one curve |
| 09 | reconstruct (positional) | `python3 -m lsquant reconstruct <N256 .chebmom1D> dos 120 --out out/ex09/eta120` | 4 | 1 s | 0 | one curve |
| 09 | reconstruct (directory) | `python3 -m lsquant reconstruct out/ex01 -b 60 --out out/ex09/dir_eta60` | 4 | 3 s | 0 | 6 curves |
| 09 | check | `python3 examples/09_inspect_reconstruct/check_reconstruct.py out/ex01 out/ex09` | 4 | 0.2 s | 0 | **PASS**: 120 meV reconstruct byte-identical to run's curve; N=4096 weight on \|E\|<1.9 = 0.7978 (bulk 0.7978) |
| 09 | rank-2 formula | `python3 -m lsquant reconstruct out/omp1/ex13/bastin_split/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D sea 100 --out out/omp1/ex09_rank2` | 1 | 3.0 s | 0 | `KuboBastinSea_…JACKSON.dat` |
| 10 | real | `python3 -m lsquant run examples/10_bases_compare/run_real.json --out out/ex10/bases_real` | 4 | 41 s | 0 | DOS + σxx, 512 trace states |
| 10 | bloch | `… run_bloch.json --out out/ex10/bases_bloch` | 4 | ≈1 s | 0 | 2 trace states |
| 10 | eigen | `… run_eigen.json --out out/ex10/bases_eigen` | 4 | ≈1 s | 0 | 2 trace states |
| 10 | check | `python3 examples/10_bases_compare/compare_bases.py out/ex10` | 4 | ≈1 s | 0 | **PASS**: bases agree to 1.3e-14 (moments), curves identical; DOS vs analytic bands 1.1e-14 |
| 11 | generator | `python3 examples/11_nonorthogonal_toy/make_nonorth_chain.py 128` | 4 | 1.5 s | 0 | cond(S') = 5.958, nnz(S') = 640, max \|eig(H',S') − 2t cos\| = 3.2e-15 |
| 11 | orthogonal | `python3 -m lsquant run examples/11_nonorthogonal_toy/run_orthogonal.json --out out/omp1/ex11/nonorth_orthogonal` | 1 | 1.8 s | 0 | |
| 11 | Löwdin exact | `… run_lowdin.json --out out/omp1/ex11/nonorth_lowdin` | 1 | 4.4 s | 0 | |
| 11 | generalised exact | `… run_generalised.json --out out/omp1/ex11/nonorth_generalised` | 1 | 2.3 s | 0 | |
| 11 | Löwdin Chebyshev | `… run_lowdin_chebyshev.json --out out/omp1/ex11/nonorth_lowdin_chebyshev` | 1 | 9.4 s | 0 | |
| 11 | check | `python3 examples/11_nonorthogonal_toy/check_nonorth.py out/omp1/ex11` | 1 | 0.8 s | 0 | **PASS**: vs exact generalised diagonalisation 4.2e-14 / 7.2e-14 / 4.2e-14 / 1.3e-13 |
| 12 | real seed 7 | `python3 -m lsquant run examples/12_disorder_nonorthogonal/run_real_seed7.json --out out/omp1/ex12b/disorder_real_seed7` | 1 | 64.4 s | 0 | dim 1152, window [-68.5, 100.2] eV |
| 12 | real seed 8 | `… run_real_seed8.json --out out/omp1/ex12b/disorder_real_seed8` | 1 | 104.0 s | 0 | |
| 12 | bloch seed 7 | `… run_bloch_seed7.json --out out/omp1/ex12b/disorder_bloch_seed7` | 1 | 20.1 s | 0 | |
| 12 | check | `python3 examples/12_disorder_nonorthogonal/check_disorder.py out/omp1/ex12b` | 1 | 2.0 s | 0 | **PASS**: real vs bloch (same disorder) 1.3e-12; ∫DOS = 72.000 (72 spin-orbitals/cell) for all three |
| 13 | run | `python3 -m lsquant run examples/13_bastin_sea_surface/run_bastin_split.json --out out/omp1/ex13/bastin_split` | 1 | 23.9 s | 0 | 4 curves from one rank-2 tensor |
| 13 | check | `python3 examples/13_bastin_sea_surface/check_bastin_split.py out/omp1/ex13/bastin_split` | 1 | 0.2 s | 0 | **PASS**: \|sea+surf−bastin\|/max = 2.3e-3; C(bastin) = C(sea) = 1.017, C(surf) = 0; C(bastin_fft) = 1.020 |
| 14 | run | `python3 -m lsquant run examples/14_wannier/run_ssh_dos.json --out out/ex14/ssh_dos` | 4 | 25.5 s | 0 | no wannier2sparse needed |
| 14 | check | `python3 examples/14_wannier/check_ssh.py out/ex14/ssh_dos` | 4 | 1.5 s | 0 | **PASS**: moments vs analytic SSH bands 6.8e-15; gap DOS 1.4e-4 vs 0.30 in band |
| 15 | generators | `python3 examples/15_spatial_operators/make_bilayer_chain.py 256` and `… 4096` | 4 | 0.1 s, 0.2 s | 0 | |
| 15 | layer DOS | `python3 -m lsquant run examples/15_spatial_operators/run_layer_dos.json --out out/omp1/ex15/layer_dos` | 1 | 18.9 s | 0 | |
| 15 | MSD | `python3 -m lsquant run examples/15_spatial_operators/run_msd.json --out out/omp1/ex15/msd` | 1 | 8.7 s | 0 | |
| 15 | check | `python3 examples/15_spatial_operators/check_spatial.py out/omp1/ex15` | 1 | 0.2 s | 0 | **PASS**: ρA+ρB=ρ to 2.9e-6 (6-digit files); ⟨E⟩ = ±0.3000; MSD slope 2.00000; ΔX²/t² = 8.435 vs 8.448 (E=0), 5.572 vs 5.518 (E=1) |
| 16 | projection | `python3 -m lsquant run examples/16_projection_downfold/run_projection_shift.json --out out/omp1/ex16b/projection_shift` | 1 | 855.9 s | 0 | window [-2, 2] eV; spectrum [-24.69, 56.37] → [-2, 2] eV; 2–8 kept bands per k (mean 4.531) |
| 16 | check | `python3 examples/16_projection_downfold/check_projection.py out/omp1/ex16b/projection_shift` | 1 | 0.3 s | 0 | **PASS**: exact DOS projected vs unprojected 2.8e-14; ∫DOS = 71.998, ∫WINDOW-DOS = 4.5312 (mean kept 4.53125) |
| 17 | Python API | `python3 examples/17_python_api/python_api.py out/omp1/ex17/level1` | 1 | 2.8 s | 0 | **PASS**: level 1 (dict) vs level 2 (blocks): moments and curves identical (0.0) |
| 18 | MCP session | `python3 examples/18_mcp_session/mcp_session.py` | 1 | 4.9 s | 0 | all 5 tools answered; background run succeeded in 2.0 s; root guard refused `/etc` |

The first attempts that failed and were then fixed (all recorded in
`SUMMARY.txt`):

- ex11 first version on the committed 512-orbital Haldane non-orthogonal
  fixture: the exact Löwdin run had not finished after 38 min at OMP 4 (killed);
  replaced by the 128-orbital chain. The generator first failed with
  `make_local_transform: cond(S') = 4.776939 outside [5.0, 10.0]` for seed
  12345; seed 1 passes.
- ex12 real-space runs first failed with `metric probes exceed 1e-13`
  (probe residual 1.4e-13 for the exact applier, cond(S) ≈ 3550); fixed with
  `metric.tol = 1e-10`. The first check expected 52 states per cell (my
  mistake; Bi DZP has 18 orbitals per atom, 72 spin-orbitals per cell).
- ex13 check, first version: required `bastin_fft` to match `bastin` at all
  energies; it differs by 11 % of the maximum away from the gap (127 nodes,
  ends at 0.64 instead of 0). The check now compares the plateau and prints the
  full-range deviation.
- ex15: custom operator names `PTOP`/`PBOT` validated but failed at run time
  (`Unknown operator 'PTOP' …`), renamed `A`/`B`; layer DOS first set up on the
  8192-orbital system with an exact trace (about 4 TFlop, killed).
- ex18: `CallToolResult` has `is_error` (mcp 2.x), not `isError`.
- ex16 first failed with `Unknown operator 'WINDOW'`; fixed with
  `model.operators = ["WINDOW"]`.

## Quickstart

Fresh output directory `out/omp1/quickstart`, OMP 1:

| step | command | wall | exit | output |
|---|---|---|---|---|
| 1 | `python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out out/omp1/quickstart/chain_dos` | 148.7 s | 0 | prints written paths; 3 `SpectralOp1…chebmom1D`, 3 `moments_chain1d_N*_M512.chebmom1D`, `eta10/ eta20/ eta24/ eta80/` subfolders, `meanSpectralOp1chain1d_N256…JACKSON.dat`, `run.resolved.json` |
| 2 | `python3 -m lsquant inspect out/omp1/quickstart/chain_dos` | 1.3 s | 0 | 305 lines of JSON (`kind: directory`, `build_profile: golden`, runtime OMP variables, model `generator`/`make_chain`, 6 moment-file headers) |
| 3 | `python3 examples/01_dos_1d_chain/plot.py out/omp1/quickstart/chain_dos` | 3.3 s | 0 | `wrote …/fig_levels.png`, `wrote …/fig_moments.png`, `wrote …/fig_dos.png` |
| 4 | `python3 -m lsquant words` | 2.2 s | 0 | JSON list of 7 shapes (below) |
| 5 | `python3 -m lsquant words describe "Tr[VX delta(E-H) VY dG-(E)/dE]"` | 3.3 s | 0 | JSON (below) |

`plot.py` found everything it needs: it reads
`moments_chain1d_N{256,1024,4096}_M512.chebmom1D` from the run directory and
`eta20/meanmoments_chain1d_N{256,1024,4096}_M512JACKSON.dat`,
`eta{80,24,10}/meanmoments_chain1d_N256_M512JACKSON.dat`, all written by step 1
(it would rerun the run file itself if the moments were missing). It writes the
three figures into the run directory (`fig_levels.png` 43 kB,
`fig_moments.png` 27 kB, `fig_dos.png` 142 kB), not into the example folder.

`python3 -m lsquant words` (shape ids in order):
`dos_identity`, `dos_operator`, `kubo_greenwood`, `kubo_bastin`, `msd`,
`projected_evolution`, `time_correlation`. Each entry has `shape_id`,
`canonical_example`, `pattern_description`, `plan_builder`,
`default_reconstruction`, `legacy_moment_prefix`, `legacy_curve_prefix`,
`rank`, `time_dependent`. Full output: `out/logs/o1_qs4_words.log`.

`python3 -m lsquant words describe "Tr[VX delta(E-H) VY dG-(E)/dE]"` (abridged;
full output `out/logs/o1_qs5_describe.log`):

```json
{
  "word": {"canonical": "Tr[VX delta(E-H) VY dG-(E)/dE]", "sea": false, "rank": 2,
           "time_dependent": false, "slots": [ ...4 slots... ]},
  "catalogue": {"shape_id": "kubo_bastin", "plan_builder": "kubo_plan",
                "default_reconstruction": "kubo_bastin", "default_formulas": ["bastin"],
                "engine_word": "VX delta VY delta", "legacy_moment_prefix": "NonEqOpVX-VY",
                "legacy_moment_ext": ".chebmom2D", ...},
  "plan": {"word": "VX delta VY delta",
           "instructions": [{"op": "Expand", "num_moments": 64},
                            {"op": "ApplyOperator", "operator_index": 1, "operator_name": "VY"},
                            {"op": "MirrorToBra"},
                            {"op": "ApplyOperator", "operator_index": 0, "operator_name": "VX"},
                            {"op": "Expand", "num_moments": 64},
                            {"op": "Trace", "num_weights": 1}]},
  "operators": ["VX", "VY"],
  "reconstruction": {"shape_id": "kubo_bastin", "formulas": ["bastin"], "kernel": "jackson",
                     "broadening_meV": null, "use_fft_grid": false, ...},
  "dt": null
}
```

Without a run file, `describe` compiles with 64 moments.

## Example 07 (quick DOS)

Not rerun here. A separate run of `examples/07_spin_hall_bi_bilayer/run_dos_quick.json`
took 12 min with 8 threads on a less loaded node; at the throughput measured here
(one niced thread receiving 10–130 % of a core) it would take well over the
15 min limit per run used here.

A later quick run of example 07 (DOS and spin Hall) took 55 min with 8 threads
on a shared node that was also running two test-suite jobs. Its oracle
(`exact_shc.json`, unprojected, Gaussian of 20 meV, `24x24x1` mesh) gives
$\sigma^z_{xy}(E_F = 0) = 1.1442\ e^2/h$, and `spin_chern.json` gives
$C_\pm = \pm 1$, spin Chern number 2, $\mathbb{Z}_2 = 1$, a gap of 0.3044 eV and
30 occupied bands. The quick KPM curve gives 0.500 at $E_F = 0$: with
$M = 256$ on the SIESTA window it resolves about 0.6 eV, more than the gap.

## Known-limitation checks

### `shc` in a csr run file

Haldane operators, `Tr[VX delta(E-H) VY dG-(E)/dE]`, 32 moments, 16 exact
trace states, 200 meV, `num_cells = 256`, `cell_area = 5.196152422706632`.
The run files (not committed) were in `limits/`, the outputs in
`out/limits/`, OMP 4.

| case | command | exit | output |
|---|---|---|---|
| (a) direct: `formulas: ["shc"]` with `num_cells`/`cell_area` set | `python3 -m lsquant run limits/run_shc.json --out out/limits/shc_direct` | 1 | `lsquant error: shc reconstruction requires num_cells and cell_area in the resolved run file` (after the moments were computed and written) |
| (b1) workaround: same with `formulas: ["bastin"]` | `python3 -m lsquant run limits/run_bastin.json --out out/limits/bastin` | 0 | `.chebmom2D` + Bastin curve; `run.resolved.json` keeps `num_cells`, `cell_area` |
| (b2) reconstruct the file | `python3 -m lsquant reconstruct out/limits/bastin/NonEqOpVX-VYhaldaneKPM_M32x32_stateexact.chebmom2D --formula shc --broadening-meV 200 --out out/limits/bastin_file_rec` | 0 | `SpinHall_NonEqOpVX-VYhaldaneKPM_M32x32_stateexactJACKSON.dat`, header `# energy_eV sigma_z_xy [e^2/h, sigma_z Pauli]` |
| (b3) reconstruct the directory | `python3 -m lsquant reconstruct out/limits/bastin --formula shc` | 0 | same file, byte-identical to (b2) (broadening taken from `run.resolved.json`) |
| (c) no `num_cells`/`cell_area` in the run file | run `limits/run_bastin_nocells.json` (exit 0), then `reconstruct <file>.chebmom2D --formula shc --broadening-meV 200` | 1 | `lsquant error: shc reconstruction requires num_cells and cell_area in the resolved run file` |
| (d) directory with a rank-1 and a rank-2 file (`words`: Bastin + DOS, `formulas: ["bastin","dos"]`, cells set) | `python3 -m lsquant reconstruct out/limits/mixed --formula shc` | 1 | writes `SpinHall_NonEqOpVX-VY…JACKSON.dat` for the rank-2 file (processed first, alphabetical order), then fails on the rank-1 file: `lsquant error: rank-2 reconstruction requires moments.rank == 2`. Partial output is left in the directory. |

So (a) fails as expected: `execute_runfile` never passes `num_cells` and
`cell_area` to `_reconstruct_single_formula` (runfile.py lines 1703, 2228,
2307). The workaround (b) works in both file and directory form. Note that the
shc header always says `sigma_z_xy [e^2/h, sigma_z Pauli]`, also for this
charge-Hall word.

### Other limitations found while writing the examples

- **Custom operator names fail at run time.** `model.operators: ["PTOP",
  "PBOT"]` with the word `Tr[PTOP delta(E-H)]` passes validation (and
  `lsquant inspect`), but `lsquant run` fails with `Unknown operator 'PTOP' …
  (registered operators: ['1', 'A', 'B', 'JXSX', …, 'VZSZ'])`. In
  `execute_runfile` the per-word preview compile (runfile.py, around line 1536:
  `compile_word(parsed_w, basis={"name": rf.basis["name"]}, …)`) passes no
  operator list, so only the default catalogue names are accepted. Example 15
  uses `A` and `B`. The same preview exists in the Wannier path.
- **`WINDOW` (projected identity) needs `model.operators: ["WINDOW"]`** in a
  projection run file; without it validation refuses the word.
- **Default `metric.tol = 1e-13`** refuses the exact applier on the Bi DZP
  fixture (probe residual 1.4e-13; cond(S) ≈ 3550). The stage6c CLI itself uses
  `tol = 1e-10`.
- **Disorder window in a non-orthogonal basis** is padded by
  √cond(S)·max|w| (±30 eV for Bi DZP at W = 1 eV), so the window is
  [-68.5, 100.2] eV and 256 moments resolve about 1 eV.
- **`reconstruction.fft_grid`** is recorded in `run.resolved.json` and in the
  compiled reconstruction (`use_fft_grid`) but no executor reads it; the
  FFT-node curve is selected by the formula `bastin_fft` / `greenwood_fft`.
- **`bastin_fft`** has only M − 1 points and deviates from `bastin` by up to
  11 % of the maximum away from plateaus (Haldane, M = 128, 50 meV).
- **`surf` is identically zero** for the clean Haldane model (6e-14), so the
  sea/surface split is uninformative for clean isotropic systems.
- **`numerics.operator = "atom_block"`** (spatial.py) is refused in the golden
  profile; not exercised. `numerics.execution = "mpi"` is refused as "not wired
  yet" in runfile.py.
- **`inspect` has no `--json` flag** (it always prints JSON); the legacy text
  moment formats report `word: ""` and `number_of_states: 0`.

## Build presets (not built)

From `CMakePresets.json` and `CMakeLists.txt` / `lsquant/CMakeLists.txt`. All
presets use `Unix Makefiles`, `CMAKE_BUILD_TYPE=Release`,
`LSQUANT_BUILD_PYTHON=ON`, `LSQUANT_BUILD_DEV_CHECKS=ON`,
`LSQUANT_BUILD_BENCHMARKS=OFF`, binary dir `build/<preset>`, CMake ≥ 3.16.
Common requirements: a C++ compiler with OpenMP, Eigen3 ≥ 3.3, FFTW3 (double,
via pkg-config), GSL, Python with development headers, pybind11.

| preset | profile | extra options | extra requirements |
|---|---|---|---|
| `golden` | golden | — | — (portable IEEE build; the one used here) |
| `native` | native (`-O3 -march=native -mtune=native -fno-math-errno`, `LSQUANT_NATIVE_BUILD=1`) | — | none beyond the common ones; binaries are tied to the CPU they were built on |
| `native-mpi` | native | `LSQUANT_WITH_MPI=ON` | an MPI implementation with C++ bindings (`find_package(MPI REQUIRED COMPONENTS CXX)`, links `MPI::MPI_CXX`) |
| `native-cuda` | native | `LSQUANT_WITH_CUDA=ON` | CUDA toolkit and `nvcc` (`enable_language(CUDA)`, `find_package(CUDAToolkit REQUIRED)`, links `CUDA::cudart`); compiles `core/algebra/atom_block_cuda.cu` with **`CUDA_ARCHITECTURES "90"` hard-coded** (Hopper). Ampere GPUs such as the RTX A6000 are sm_86, so this preset as written does not produce code for them; the architecture needs to be changed or made a cache variable before such a build. |
| `full` | native | MPI + CUDA | both of the above |

`LSQUANT_WITH_GINKGO` (not in any preset) needs a Ginkgo installation, and a
CUDA-enabled one with `LSQUANT_WITH_CUDA`. CUDA, MPI and Ginkgo are refused
with the golden profile (CMake `FATAL_ERROR`).

## Example 04 check (added 2026-11-21)

The plot script of example 04 now checks the precession (`check_precession`
in `plot.py`). The operators `goldens/stage0/chain1d_mag/operators/*.CSR` were
copied into the run directory; their SHA-256 sums match this branch. Golden
build, OMP 1, `nice -n 19`:

| example | step | command | OMP | wall | exit |
|---|---|---|---|---|---|
| 04 | run | `python3 -m lsquant run examples/04_spin_precession_magnetic_chain/run_spin.json --out out/ex04_check` | 1 | 12.2 s | 0 |
| 04 | plot + check | `python3 examples/04_spin_precession_magnetic_chain/plot.py out/ex04_check` | 1 | 5.1 s | 0 |

The check prints a fitted frequency of 0.303853487060 rad/fs against
2J_ex/ħ = 0.303853487060 rad/fs (relative difference 3.7e-16), the largest
deviation from the analytic curves as 4.0e-15 (S_x), 3.2e-15 (S_y) and 0 (S_z),
and CHECK PASS.

## Not kept

None of the planned examples was dropped. Example 07 quick was not rerun
(see above). The build presets were documented, not built.
