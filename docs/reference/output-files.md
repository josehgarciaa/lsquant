# Output files

This page lists the files that [`python -m lsquant run`](cli.md#run) and
[`python -m lsquant reconstruct`](cli.md#reconstruct) write, with their names,
formats, columns and units. The sources are `python/lsquant/runfile.py`,
`python/lsquant/reconstruct.py`, the SIESTA executors in
`python/lsquant/apps/` and `python/lsquant/projection/runfile.py`, and the C++
output functions in `lsquant/io/` (bound in `lsquant._core` as `WriteChebmom`,
`WriteLegacyChebmom1D`, `WriteLegacyChebmom2D`, `WriteLegacyChebmomTD` and
`WriteCurveText`).

## Output directory

The output directory is `--out` if you give it, otherwise `output.dir` from
the run file, otherwise the current directory. LSQuant creates it if needed
and overwrites files of the same name. Two runs with the same output
directory, such as the DOS and spin Hall runs of
`examples/07_spin_hall_bi_bilayer/`, share one `run.resolved.json` and one
`basis.json`: the second run replaces them.

`run --json` prints the paths it wrote as `moment_files`, `curve_files` and
`written_files`.

## Files by model kind

| File | `csr` | `wannier` | `siesta` (clean) | `siesta` (disordered) | `siesta` (projection) | `generator` |
|---|---|---|---|---|---|---|
| `run.resolved.json` | yes | yes | yes | yes | yes | yes |
| Legacy moment file (`.chebmom1D`, `.chebmom2D`, `.chebmomTD`) | if `legacy_layout` or `legacy_text_round_trip` | if `legacy_layout` | if `legacy_layout` | if `legacy_layout` or `legacy_text_round_trip` | `word_<i>.*` if `legacy_layout` | driver decides |
| Unified moment file (`.chebmom`) | if `write_chebmom` or not `legacy_layout` | if `write_chebmom` or not `legacy_layout` | always | always | `word_<i>.chebmom`, always | driver decides |
| Legacy curve files (`*JACKSON.dat`) | if `legacy_names` | if `legacy_names` | — | — | — | driver decides |
| `curves.json` | if not `legacy_names` | if not `legacy_names` | — | — | — | — |
| `dos.dat`, `shc.dat` | — | — | yes | yes | — | — |
| `basis.json` | — | `<stem>.basis.json` if `write_basis_json` | yes | yes | — | — |
| Other files | see [`csr` extras](#csr-extras) | — | see [Clean SIESTA runs](#clean-siesta-runs) | see [Disordered SIESTA runs](#disordered-siesta-runs) | see [Projection runs](#projection-runs) | see [Generator runs](#generator-runs) |

The run-file keys in this table are described in the
[Run-file reference](run-file.md#output).

## `run.resolved.json`

The run file after validation, with every default filled in, written to the
output directory by every run. You can pass it back to `run` to repeat the
calculation, provided the loaded engine has the same build profile.

Compared with the input run file, it contains:

| Key | Content |
|---|---|
| `build_profile` | Build profile of the loaded engine, `golden` or `native`. |
| `runtime` | Values of `OMP_NUM_THREADS`, `OMP_PROC_BIND`, `OMP_PLACES`, `OMP_DYNAMIC`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `LSQUANT_SPMV_MIN_NNZ`, `LSQUANT_KERNEL_BATCH` and `LSQUANT_FFT_BACKEND` (`null` if unset). |
| `word`, `words` | Canonical forms. `words` is present only if it differs from `[word]`. |
| `model.dir` and other model paths | Absolute paths. |
| `output.dir` | Absolute path of the directory actually used. |
| `projection`, `numerics`, `reference` | Present only when set. |

SIESTA executors also write:

| Key | Content |
|---|---|
| `expansion.window` | Numeric bounds `[lower, upper]` of the window used, in eV. |
| `reconstruction.num_cells` | Number of k-points (or supercell cells). |
| `reconstruction.cell_area` | \|a1 × a2\| of the unit cell, in Å². |
| `receipt` | Phase list (`phase`, `status`, `wall_seconds`) and diagnostics. Clean runs add `receipt.reference` (`enabled`, `reason`). Disordered runs add `executor`, `basis` (a copy of `basis.json`), `total_wall_seconds` and `fermi_values`. Projection runs add `trace_basis`, `metric_residuals` and `projection`. |
| `model.params.nonorthogonal_disorder` | Set to `true` by the disordered executor, so that a replay uses the same executor. |

`csr` and `wannier` runs write `run.resolved.json` once, before computing.
SIESTA runs rewrite it at every phase boundary, so a failed run still records
which phase failed. `python -m lsquant reconstruct` reads
`reconstruction.broadening_meV`, `fermi_energies`, `num_cells` and `cell_area`
from the `run.resolved.json` next to the moments.

## Moment files

### Names

`csr` runs name moment files from the word's legacy prefix, the model label,
the moment count and the state label:

| Moments | Stem |
|---|---|
| Rank 1 (`dos_*` shapes, `projected_evolution`) | `<prefix><label>KPM_M<M>_state<state>` |
| Rank 2 (`kubo_*` shapes) | `<prefix><label>KPM_M<M>x<M>_state<state>` |
| Time-dependent `msd` and `time_correlation` | `<prefix><label>KPM_M<M>_state<tmax><state>` |

- `<label>` is `model.label`.
- `<state>` is `expansion.state_label` if set, otherwise `exact` for an exact
  trace and `default` for a random-phase trace.
- `<tmax>` is `time.tmax` (or `dt × (steps − 1)`) printed with `%g`, for
  example `150`.

Other executors use these stems:

| Executor | Stem |
|---|---|
| `wannier` | `<prefix><seed>_<basis>_KPM_M<M>_state<state>`, or `..._KPM_M<M>x<M>_state<state>` for rank-2 tensors (time-dependent tensors are rank 2) |
| SIESTA clean and disordered | `<prefix><label>KPM_M<M>_state<state>` or `<prefix><label>KPM_M<M>x<M>_state<state>` |
| SIESTA projection | `word_<i>`, where `<i>` is the index of the word in `words` |

Moment files keep these names whatever the value of `output.legacy_names`;
that key only chooses how curves are written (see [Curve files](#curve-files)).

The prefix and the legacy extension come from the word shape
(`words/catalogue.py`). `<A>`, `<B>`, `<V>` and `<P>` are operator names as
they appear in the word, for example `VX` or `PX` (from `P:PX`).

| Word shape | Prefix | Legacy extension |
|---|---|---|
| `dos_identity` | `SpectralOp1` | `.chebmom1D` |
| `dos_operator` | `SpectralOp<A>` | `.chebmom1D` |
| `kubo_greenwood` | `NonEqOp<A>-<B>` | `.chebmom2D` |
| `kubo_bastin` | `NonEqOp<A>-<B>` | `.chebmom2D` |
| `msd` | `Correlation<V>`, with `<V>` the velocity of the commutator axis | `.chebmomTD` |
| `projected_evolution` | `TimeProj<A>-<P>` | `.chebmomTD` |
| `time_correlation` | `TimeCorr<A>-<B>` | `.chebmomTD` |

The full name is the stem plus the legacy extension or `.chebmom`. Examples:

- `NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D`
  (`examples/03_haldane_hall/run_haldane.json`)
- `TimeProjSX-PXchain1d_magKPM_M32_stateexact.chebmomTD`
  (`examples/04_spin_precession_magnetic_chain/run_spin.json`)

### Legacy formats

The legacy files reproduce the moment files of the legacy code byte for byte.
Numbers are written with 15 significant digits. `bandwidth` is
`upper − lower` and `center` is `(upper + lower) / 2` of the spectral window,
in eV. A build other than `golden` appends ` # build_profile=<profile>` to the
first line.

| File | Line 1 | Line 2 | Data lines |
|---|---|---|---|
| `.chebmom1D` | `systemSize bandwidth center` | `M` | M lines `re im` |
| `.chebmom2D` | `systemSize bandwidth  center` (two spaces before `center`) | `M0 M1 ` (trailing space) | M0 × M1 lines `re im`, row-major |
| `.chebmomTD` | `systemSize bandwidth center N_t dt ` (trailing space) | `M N_t ` (trailing space) | M × N_t lines `re im`, moment index outer, time index inner |

`dt` is the time step in fs. When LSQuant reads a legacy file back, it takes
the state label from the file name after `_state`; the number of states and
the word are not stored.

### Unified `.chebmom` format

A text file with a header and full-precision data (17 significant digits):

```text
chebmom 1
rank <r>
dims <d0> [<d1>]
window <lower> <upper>
system_size <N>
states <label> <count>
word <word>
kpm_weighted <0|1>
symmetrised <0|1>
data
<re> <im>
...
```

Data are row-major. In a build other than `golden`, the first line ends with
the build-profile mark described under [Legacy formats](#legacy-formats). The
format has no time-step field, so `reconstruct` reads a time-dependent
`.chebmom` file as a rank-2 tensor without time steps. Reconstruct time curves from the `.chebmomTD` file.

`python -m lsquant inspect <file>` prints the metadata of any moment file.

## Curve files

### Names

With `output.legacy_names = true` (the default), `csr` and `wannier` runs
write one file per word and formula, named `<curve prefix><stem><suffix>`,
where `<stem>` is the moment-file stem. `reconstruct` uses the same rule with
the moment file name minus its extension as the stem.

| Formula | Curve prefix | Suffix |
|---|---|---|
| `dos` | `mean` | `JACKSON.dat` |
| `greenwood` | `KuboGreenWood_` | `JACKSON.dat` |
| `bastin` | `KuboBastin_` | `JACKSON.dat` |
| `sea` | `KuboBastinSea_` | `JACKSON.dat` |
| `surf` | `KuboBastinSurf_` | `JACKSON.dat` |
| `greenwood_fft` | `KuboGreenWood_` | `JACKSON_FFTgrid.dat` |
| `bastin_fft` | `KuboBastin_` | `JACKSON_FFTgrid.dat` |
| `time_correlation` | `mean` | `EF<E_F>JACKSON.dat`, one file per Fermi energy, `<E_F>` with 6 decimals |
| `shc` | `SpinHall_` | `JACKSON.dat` (`reconstruct` only; see [`shc` curves](#shc-curves)) |
| `row_m0` | — | no file |

Examples:

- `meanSpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256JACKSON.dat`
- `KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat`
- `meanCorrelationVX<label>KPM_M256_state150defaultEF2.000000JACKSON.dat`
  (an `msd` word with `tmax = 150`, a random-phase trace and `fermi_energies = [2.0]`)

The `legacy_curve_prefix` field that `python -m lsquant words list` shows is
catalogue metadata; the run uses the prefixes in this table.

### Format

`WriteCurveText` writes one point per line, in two columns separated by one
space. The file has no header line. Numbers use the default C++ stream format
(6 significant digits).

| Formula | Column 1 | Column 2 | Points |
|---|---|---|---|
| `dos` | Energy, eV | Extensive density of states N ρ(E), 1/eV (N is the dimension of H) | 30 × M_eff, uniform in the rescaled interval [−0.95, 0.95] |
| `greenwood` | Energy, eV | Kubo–Greenwood value; the code states no unit | M_eff − 1, uniform in [−0.95, 0.95] |
| `bastin`, `sea`, `surf` | Energy, eV | Kubo–Bastin value (total, Fermi sea, Fermi surface); the code states no unit | 30 × M_eff − 1, uniform in [−0.95, 0.95] |
| `greenwood_fft`, `bastin_fft` | Energy, eV | As above | M_eff − 1, on the Chebyshev nodes |
| `time_correlation` | Time t = n × dt, fs | Correlation at the Fermi energy; the code states no unit | N_t |

M_eff = min(⌈π/η⌉, M) with η = 2 × `broadening_meV` / 1000 / (E_max − E_min):
the number of moments the Jackson kernel keeps for the requested broadening.
For rank-2 moments M_eff is the larger of the two truncated dimensions.
Energies are in eV because the reconstruction converts `broadening_meV` with
a factor 1/1000; this assumes H is in eV.

### `curves.json`

With `output.legacy_names = false`, `csr` and `wannier` runs put every curve
into `curves.json` in place of the `.dat` files:

```json
{
  "<stem>:<formula>": {"grid": [...], "values": [...]}
}
```

The keys of each entry are the same columns as the `.dat` file. Curves of
[reconstruction targets](#reconstruction-targets) are still written as
`.dat` files.

### `shc` curves

`shc` gives the spin Hall conductivity σ^z_xy in units of e²/h, with the spin
current built from Pauli matrices (σ_z, not ħ/2 σ_z). It is the Kubo–Bastin
curve multiplied by 2π / (`num_cells` × `cell_area`).

- SIESTA runs write `shc.dat` (see below).
- `python -m lsquant reconstruct <file> shc` writes
  `SpinHall_<stem>JACKSON.dat` with the header line
  `# energy_eV sigma_z_xy [e^2/h, sigma_z Pauli]` and 17 significant digits.
  It reads `num_cells` and `cell_area` from `run.resolved.json`, as described
  above.
- A `csr` or `wannier` run with `formulas: ["shc"]` fails, because these
  executors do not pass `num_cells` and `cell_area` to the formula.

## Reconstruction targets

Each entry of `reconstruction.targets` (`csr` runs only) writes its curves to
`<output dir>/<subdir>/` (or the output directory if `subdir` is absent),
with the names and format of [Curve files](#curve-files). With
`also_top_dir = true`, LSQuant also copies each curve to the output directory,
unless a file of that name already exists there.

## Per-site DOS

With `reconstruction.per_site_dos = true`, `csr` runs also write
`DOS_reconstructed.dat` for every rank-1 word. The file has two columns,
energy in eV and the density of states divided by the dimension of H, in 1/eV,
printed with `%.10g`. It uses all M moments, whatever the broadening, and
drops every line whose first two values are `nan` or `inf`. With several
rank-1 words, each one overwrites the file.

## `csr` extras

| File | Written when | Content |
|---|---|---|
| `checkpoint_v2_<hash>.lsqckpt` | `expansion.checkpoint = true` | Trace-state checkpoint; `<hash>` is the full SHA-256 calculation fingerprint (settings and model input contents). Legacy and different fingerprints are reported and ignored. Deleted when the word finishes. |
| `spatial_receipt.json` | `numerics.operator = "atom_block"` | Source of the atom blocks, domain size, spatial ordering, block counts and construction times. |
| `permutation.npy` | `numerics.operator = "atom_block"` | Orbital permutation of the spatial ordering (NumPy array). |
| `captured/NNNN_<name>` | `output.capture_tutorial = true` | Numbered copies of the curve files and of `projected_evolution` moment files, for the tutorial golden check. |

## `wannier` extras

With `output.write_basis_json = true` (the default for `wannier`), each word
writes `<stem>.basis.json` with the keys `basis`, `seed`, `supercell`,
`mesh`, `dimension`, `trace`, `states`, `state_seed`, `window_source`
(`eigen`, `gershgorin` or `explicit`), `window`, `disorder_w`,
`disorder_seed` (`null` without disorder), `moments`, `broadening_meV`,
`word` (`dos`, `kuboxx` or `kuboxy`) and `operators` (the source of each
operator).

## Clean SIESTA runs

`apps/stage6b_siesta.py` writes, in addition to the moment files:

| File | Written when | Content |
|---|---|---|
| `dos.dat` | rank-1 word | Header `# energy_eV DOS [states/eV/cell]`. Energy in eV, DOS per unit cell (the extensive DOS divided by `num_cells`). |
| `shc.dat` | rank-2 word | Header `# energy_eV sigma_z_xy [e^2/h, sigma_z Pauli]`. Energy in eV, σ^z_xy in e²/h. |
| `basis.json` | always | Model metadata (sisl version, input files, energy and length units, printed Fermi level, sisl layout checks, orbital list, dipole quadrature certificate and cache), `metric_residuals`, `window`, `velocity_mode`, `route`, `metric_approximation`, `trace`, `trace_basis`, `units`, `num_cells`, `cell_area`. |
| `eigen_hermiticity.json` | always | Hermiticity residuals of the eigenbasis operators. |
| `eigen_comparison.json` | always | One entry per word: `status` (`completed` or `skipped`), the reason if skipped, and otherwise `max_moment_difference` and `max_curve_difference` between the production run and the exact eigenbasis reference. |
| `eigen_dos.chebmom`, `eigen_dos.dat`, `eigen_shc.chebmom`, `eigen_shc.dat` | reference enabled | Moments and curve of the exact eigenbasis reference, same format as `dos.dat` and `shc.dat`. |
| `exact_shc.json` | `model.params.oracle = true` | Exact per-k spin Hall oracle: `energies` (eV), `values` (e²/h), `units`, `trace`, `degeneracy`, `gaussian_standard_deviation_eV` (= `broadening_meV` / 1000). |
| `spin_chern.json` | `model.params.oracle = true` | Topology at the first Fermi energy: `chern_plus`, `chern_minus`, `spin_chern` (C₊ − C₋), `prodan_half_difference`, `z2`, `projected_spin_gap`, `occupied_bands`, `energy_gap`, `convention`; or `defined: false` with a `reason` when the Fermi energy is not in a gap. Always adds `fermi_energy_eV`, `degeneracy`, `sigma_xy` (e²/h) and `units`. |

`.dat` files from SIESTA runs use `numpy.savetxt` with 17 significant digits.
With `metric.route = "anticommutator"`, the headers end with
`; approximate anticommutator metric`. With several words, each rank-1 word
writes `dos.dat` and each rank-2 word writes `shc.dat`, so later words
overwrite earlier ones.

## Disordered SIESTA runs

`apps/stage6c_disorder.py` writes the moment files, `dos.dat` and `shc.dat`
(same headers as clean runs), and:

| File | Content |
|---|---|
| `basis.json` | Model metadata, `layout`, `route`, `metric_applier`, `metric_residuals`, `disorder` (`W_eV`, `seed`, `kind`, `sha256` of the disorder values, distribution, potential form), `window`, `window_source`, `overlap_bounds`, `num_cells`, `cell_area`, `trace`, `trace_basis`, `units`. |
| `disorder.npy` | The disorder values (NumPy array). |
| `fermi_values.json` | One entry per word and Fermi energy: `word`, `fermi_energy_eV`, `value` (the curve interpolated at that energy), `units`. |

## Projection runs

`projection/runfile.py` writes, for the word with index `<i>`:

| File | Content |
|---|---|
| `word_<i>.chebmom` | Moments. |
| `word_<i>.chebmom1D` or `word_<i>.chebmom2D` | Legacy moments, if `output.legacy_layout`. |
| `word_<i>.dat` | Header `# energy_eV value [<units>]`: e²/h with Pauli σ_z for the spin Hall word, `states/eV/cell` for rank-1 words. |
| `word_<i>_oracle.json` | Exact Gaussian-broadened curve of the word before and after projection on the oracle grid: `energies`, `unprojected`, `projected`, the oracle error and diagnostics. |
| `projection.json` | Projection receipt, off-grid interpolation errors, and per-word errors, including the KPM curve against the projected oracle. |

## Generator runs

For `model.kind = "generator"`, the tutorial driver decides the files and
ignores most `output` keys. Each example README lists the files; for example,
`examples/01_dos_1d_chain/run_dos.json` writes:

- `SpectralOp1chain1d_N<N>KPM_M512_stateexact_chain1d_N<N>.chebmom1D` and
  `moments_chain1d_N<N>_M512.chebmom1D` for each chain length;
- `meanSpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256JACKSON.dat`;
- `eta<meV>/meanmoments_chain1d_N<N>_M512JACKSON.dat` for each broadening;
- the operator files it builds, under `operators/`.

## Files outside the output directory

With `velocity_mode = "covariant"`, SIESTA runs cache the certified dipole
matrices as `dipoles-<key>.npz`. The cache directory is `LSQUANT_CACHE_DIR`
if set; otherwise `.lsquant-cache` next to the SIESTA directory if that
directory is named `siesta`, `.lsquant-cache` inside it if it is writable, or
`lsquant-cache` in the system temporary directory. A cache entry is reused
only if its key matches the geometry, overlap, radial files and quadrature
settings.

## See also

- [Run-file reference](run-file.md)
- [Command-line interface](cli.md)
- [Observable words](observable-words.md)
