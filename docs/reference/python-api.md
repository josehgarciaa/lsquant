# Python API

This page lists the public Python API of the package `lsquant` (`python/lsquant/`).
"Public" means the names exported by `lsquant/__init__.py` and
`lsquant/words/__init__.py`, the model loaders of `lsquant/models/`, and the
building blocks used by [example 17](../../examples/17_python_api/README.md).

To import the package, `PYTHONPATH` must contain the build directory of the
compiled module and the `python/` directory of the repository:

```bash
export PYTHONPATH="$PWD/build/<preset>/lsquant_core:$PWD/python"
```

The C++ types that appear in the signatures (`Plan`, `SparseMatrix`,
`LinearOperator`, `SpectralWindow`, `StateGenerator`, `ChebyshevMoments`) belong
to the compiled module `lsquant._core` and are listed in the
[core API reference](core-api.md). For the run-file keys, see the
[run-file reference](run-file.md). For a guided use, see
[How we run LSQuant](../how-to/run-lsquant.md) and
[example 17: LSQuant from Python](../../examples/17_python_api/README.md).

## Import-time side effects

`import lsquant` changes the process environment before it loads the compiled
module (`runtime.py`, called from `__init__.py:14` and `__init__.py:25`).

| Variable | Effect | Source |
|---|---|---|
| `OMP_PROC_BIND` | Set to `spread` if not already set. | `runtime.py:30` |
| `OMP_PLACES` | Set to `cores` if not already set. | `runtime.py:31` |
| `OMP_NUM_THREADS` | Set to the number of physical cores in the process CPU affinity if not already set. | `runtime.py:32–33` |
| `LSQ_ENGINE_MEM` | Any value turns on the allocation report of the native rank-2 path, on stderr. Passed to `_core.ConfigureEngineDiagnostics`. | `runtime.py:44–49` |
| `LSQ_ENGINE_TIMING` | Any value turns on the phase timings of the native rank-2 path, on stderr. Passed to `_core.ConfigureEngineDiagnostics`. | `runtime.py:44–49` |
| `LSQ_CHECKPOINT_ABORT_AFTER` | `k` stops a checkpointed run right after the k-th checkpoint write (default `-1`, never). Used by the resume test. Passed to `_core.ConfigureEngineDiagnostics`. | `runtime.py:44–49` |
| `LSQUANT_SPMV_MIN_NNZ` | Read once (default `16384`) and passed to `_core.ConfigurePerformance`, `native` builds only. | `runtime.py:50–54` |
| `LSQUANT_KERNEL_BATCH` | Read once (default `0`) and passed to `_core.ConfigurePerformance`, `native` builds only. | `runtime.py:50–54` |

Set these variables before the first `import lsquant`. A change after the
import has no effect on the OpenMP runtime.

> **Known limitation.** If `OMP_NUM_THREADS` is unset, the import starts one
> thread per physical core (96 on a 96-core node). On a shared or busy node this
> is much slower than one thread for small problems. Set `OMP_NUM_THREADS`
> yourself, before you import `lsquant`. To choose the value, see
> [Run on a cluster](../how-to/run-on-a-cluster.md).

Other environment variables read by the Python layer:

| Variable | Read by | Effect |
|---|---|---|
| `LSQUANT_FFT_BACKEND` | `kubo_fft_grid` (`reconstruct.py:588`) | Default `backend` (`tables` or `fft`) when the argument is `None`. |

## Top-level names

`from lsquant import *` exports the following names (`__init__.py:67–118`).

| Name | Kind | Defined in |
|---|---|---|
| `_core` | module | compiled extension, see [core-api.md](core-api.md) |
| `bases`, `cli`, `kernels`, `models`, `moments`, `reconstruct`, `runfile`, `spectral_functions`, `words` | modules | `python/lsquant/` |
| `RunFile`, `execute_runfile`, `reconstruct_target`, `inspect_target` | class, functions | `runfile.py` |
| `RealBasis`, `BlochBasis`, `EigenBasis`, `BlochMetricBasis`, `EigenMetricBasis`, `make_basis` | classes, function | `bases.py` |
| `HBAR`, `dos_plan`, `kubo_plan`, `msd_plan`, `projected_evolution_plan`, `time_correlation_plan`, `time_evolution_series`, `time_step_from_tmax`, `window_from_gershgorin`, `run_plan` | constant, functions | `moments.py` |
| `dos_curve`, `kubo_greenwood`, `kubo_bastin`, `kubo_bastin_sea`, `kubo_bastin_surf`, `kubo_fft_grid`, `time_correlation_curve`, `shc` | functions | `reconstruct.py` |
| `jackson_damping`, `jackson_damping_2d` | functions | `kernels.py` |
| `SpectralFunction`, `TimeEvolutionFunction`, `DeltaFunction`, `GaussianDelta`, `GreenRetarded`, `GreenRetardedDerivative`, `TimeEvolution`, `bessel_series_coefficients`, `legacy_bessel_coefficients`, `legacy_delta_table`, `legacy_green_tables` | classes, instances, functions | `spectral_functions.py` |

## `lsquant.runfile`

### `execute_runfile`

```python
execute_runfile(runfile, *, out_dir_override=None, repo_root=None) -> dict[str, Any]
```

Validates a run file, runs it, and writes moments, curves and
`run.resolved.json`. This is the function behind `python -m lsquant run`.

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `runfile` | `RunFile`, mapping or `str` | — | A `RunFile`, a run-file dictionary, or the path of a JSON run file. |
| `out_dir_override` | `str` or `None` | `None` | Output directory. Replaces `output.dir` of the run file. |
| `repo_root` | `str` or `None` | `None` | Root used to resolve relative model paths. `None` means the repository root (two levels above `python/lsquant/`). |

Returns a dictionary:

| Key | Type | Meaning |
|---|---|---|
| `output_dir` | `str` | Absolute output directory. |
| `resolved_runfile` | `str` | Path of `run.resolved.json`. |
| `moment_files` | `list[str]` | Moment files written (`.chebmom`, `.chebmom1D`, `.chebmom2D`, `.chebmomTD`). |
| `curve_files` | `list[str]` | Curve files written (`.dat`). |
| `written_files` | `list[str]` | Every file written. For SIESTA runs: every file in the output directory. |

Dispatch by `model.kind` (`runfile.py:1480–1491`): a SIESTA model with
`projection.kind != "none"` goes to the projection executor; a SIESTA model with
`basis.name = "real"`, `basis.disorder.W != 0` or
`model.params.nonorthogonal_disorder` goes to the disordered executor; any other
SIESTA model goes to the clean executor. See the
[run-file reference](run-file.md#siesta-executors).

Raises `ValueError` for an invalid run file. Raises `FileNotFoundError` when
`<label>.HAM.CSR` of a `csr` model is missing, or when `<label>.S.CSR` is
missing and the metric is non-orthogonal (`runfile.py:1502–1515`).

A relative `output.dir` (or `out_dir_override`) is taken relative to the current
working directory (`runfile.py:1465–1467`). Relative model paths are resolved
as described in [Path resolution](run-file.md#path-resolution).

> **Known limitation: operator names.** Before it runs each word,
> `execute_runfile` compiles it with `basis={"name": ...}` only
> (`runfile.py:1547`). The compiler then checks the operators against the
> built-in list `DEFAULT_KNOWN_OPERATORS` (`words/compiler.py:28`, fallback at
> `words/compiler.py:393–394`) and raises `ValueError("Unknown operator ...")`
> (`words/compiler.py:405`). A custom name listed in `model.operators` passes
> `RunFile.from_dict` and `inspect`, but fails here. Use the catalogue names
> (`A`, `B`, `O`, `P`, `PX`, `PY`, `PZ`, `SX`, …) for your own operators.

### `RunFile`

```python
@dataclass
class RunFile:
    model: dict; basis: dict; metric: dict; word: str; words: list[str]
    expansion: dict; time: dict | None; reconstruction: dict; output: dict
    source_path: str | None = None; reference: bool | None = None
    projection: dict = {"kind": "none", "window": None, "radius": None}
    numerics: dict = {}
```

A validated run file with every default filled in (`runfile.py:864`). Each
field holds one section of the run file, as listed in the
[run-file reference](run-file.md).

| Method | Signature | Meaning |
|---|---|---|
| `from_dict` | `RunFile.from_dict(raw, *, source_path=None, repo_root=None) -> RunFile` | Validate `raw` with `validate_runfile_dict`, fill defaults, resolve paths. `source_path` is the file the dictionary came from; relative paths are also tried against its directory. |
| `from_file` | `RunFile.from_file(path, *, repo_root=None) -> RunFile` | Read a JSON file and call `from_dict` with `source_path=path`. |
| `to_resolved_dict` | `to_resolved_dict() -> dict` | The dictionary written to `run.resolved.json`, including `build_profile` and the `runtime` environment record. |

Raises `ValueError` on any validation error.

### `reconstruct_target`

```python
reconstruct_target(target_path, *, formula=None, broadening_meV=None,
                   fermi_energies=None, out_dir=None) -> dict[str, Any]
```

Reconstructs curves from a moment file or from every moment file in a
directory. This is the function behind `python -m lsquant reconstruct`.

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `target_path` | `str` | — | A moment file, or a directory with moment files. |
| `formula` | `str` or `None` | `None` | Reconstruction formula. `None`: `time_correlation` for time-dependent moments, `bastin` for rank 2, `dos` for rank 1 (`runfile.py:2423–2429`). |
| `broadening_meV` | `float` or `None` | `None` | Broadening in meV. `None`: for a directory, `reconstruction.broadening_meV` of its `run.resolved.json`, else 10.0; for a file, 10.0. |
| `fermi_energies` | sequence of `float` or `None` | `None` | Fermi energies (eV) for `time_correlation`. `None`: for a directory, the value in `run.resolved.json`, else `[0.0]`; for a file, `[0.0]`. |
| `out_dir` | `str` or `None` | `None` | Where to write the curves. `None`: next to the moments. |

Returns, for a file: `moment_file`, `output_dir`, `formula`, `broadening_meV`,
`fermi_energies`, `curve_files`. For a directory: `output_dir`, `curve_files`.

For a directory, it reads `*.chebmom1D`, `*.chebmom2D` and `*.chebmomTD`, and
only if there are none, `*.chebmom` (`runfile.py:2364–2372`). For `shc`, it
takes `num_cells` and `cell_area` from the `run.resolved.json` next to the
moment file (`runfile.py:2431–2436`).

Raises `FileNotFoundError` if the path or the moment files are missing.

### `inspect_target`

```python
inspect_target(path) -> dict[str, Any]
```

Describes a moment file, a run file or an output directory. This is the
function behind `python -m lsquant inspect`.

| `kind` returned | Input | Other keys |
|---|---|---|
| `"directory"` | a directory | `path`, `resolved_runfile` (content of `run.resolved.json` or `None`), `defaults`, `moment_files` (one `inspect_target` result per `*.chebmom*` file) |
| `"runfile"` | a `.json` file | `path`, `resolved_runfile`, `defaults` |
| `"moments"` | any other file | `path`, `rank`, `dims`, `system_size`, `number_of_states`, `state_label`, `word`, `window`, `band_width`, `half_width`, `band_center`, `kpm_weighted`, `symmetrised`, `number_of_time_steps`, `time_step` |

Raises `FileNotFoundError` if `path` is missing, and `ValueError` if a
`.json` file is an invalid run file.

### Other functions

| Signature | Meaning |
|---|---|
| `validate_runfile_dict(raw) -> None` | Validate a raw run-file dictionary before any file is read. Raises `ValueError`. |
| `defaults_for_basis(basis_name='real') -> dict` | Copy of `DEFAULT_RUNFILE_DEFAULTS` with `expansion.window` set for the basis. This is the `defaults` block printed by `inspect`. |
| `default_window_for_basis(basis_name) -> str` | Default spectral-window mode for a basis. |
| `resolve_generator_path(name, repo_root=None) -> str` | Path of a generator script from its catalogue name. |

## `lsquant.bases`

A basis object holds the Hamiltonian, the static operators and the trace
states of one model on one supercell or k mesh. It is the input of `run_plan`.

### `make_basis`

```python
make_basis(basis_name, model, supercell, *, disorder_w=0.0, seed=12345,
           disorder_seed=None, disorder=None, csr_dir=None)
    -> RealBasis | BlochBasis | EigenBasis
```

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `basis_name` | `str` | — | `"real"`, `"bloch"` or `"eigen"` (case-insensitive). |
| `model` | `TightBindingModel` | — | A model with `num_wann`, `uc`, `tau`, `Hk`, `dHk`, `k_mesh`, `bloch_phases`, `hamiltonian_csr`, `velocity_csr`, for example a `WannierModel`. |
| `supercell` | sequence of 3 `int` | — | Supercell (`real`) or k mesh (`bloch`, `eigen`) N₁ × N₂ × N₃. |
| `disorder_w` | `float` | `0.0` | Anderson box disorder width W in eV, on-site values in [−W/2, W/2]. |
| `seed` | `int` | `12345` | Seed of random-phase states, and of disorder if `disorder_seed` is `None`. |
| `disorder_seed` | `int` or `None` | `None` | Seed of the disorder values. |
| `disorder` | array or `None` | `None` | Explicit on-site values, length `dimension`. Overrides `disorder_w`. |
| `csr_dir` | `str` or `None` | `None` | Directory of `wannier2sparse` CSR files (`real` only). |

Raises `ValueError` for an unknown name. For `"eigen"`, `make_basis` passes
only `disorder_w` and `seed` on and drops `disorder_seed`, `disorder` and
`csr_dir`; `disorder_w != 0` raises `ValueError` (`bases.py:925–928`).

### `RealBasis`, `BlochBasis`, `EigenBasis`

```python
RealBasis(model, supercell, *, disorder_w=0.0, seed=12345, disorder_seed=None,
          disorder=None, csr_dir=None, use_w2s_csr=True)
BlochBasis(model, supercell, *, disorder_w=0.0, seed=12345, disorder_seed=None, disorder=None)
EigenBasis(model, supercell, *, disorder_w=0.0, seed=12345)
```

| Class | Representation | Disorder |
|---|---|---|
| `RealBasis` | Real-space supercell CSR. Reads `<seed>.HAM.CSR` and `<seed>.V*.CSR` from `csr_dir` or `model.model_dir` if present and `use_w2s_csr` is true; otherwise builds them. | yes (diagonal operator added to H) |
| `BlochBasis` | One block H(k) per k-point (block-diagonal operator). | yes (applied between Bloch transforms) |
| `EigenBasis` | Exact eigenstates of each block H(k). | no (`ValueError` if `disorder_w != 0`) |

Common attributes and methods:

| Member | Meaning |
|---|---|
| `dimension` | Hilbert-space size N₁N₂N₃ × `num_wann`. |
| `supercell`, `num_cells`, `num_wann`, `seed` | As constructed. |
| `hamiltonian()` | H as a `_core.LinearOperator`. |
| `operator(name)` | Static operator: `"1"` (also `I`, `ID`, `IDENTITY`), `"H"` (`HAM`), `"VX"`, `"VY"`, `"VZ"`. Other names raise `ValueError`. |
| `states(trace_mode='exact', num_states=None, seed=None)` | A `_core.StateGenerator`. `"exact"`: local states (all `dimension` states by default). `"random"`: random-phase states (1 by default) with seed `seed` or the basis seed. |
| `window(pad=0.1)` | `_core.SpectralWindow`. `RealBasis`: Gershgorin bound widened by `pad` and by W/2·(1 + `pad`) with disorder. `BlochBasis`, `EigenBasis`: block-eigenvalue range widened by `pad` × width. |
| `operator_sources`, `operator_source(name)` | Where each operator came from (CSR file or built). |
| `EigenBasis.exact_dos(energies, eta)`, `EigenBasis.exact_sigma_xx(energies, eta)` | Reference DOS and σ_xx from the exact eigenvalues with broadening `eta`. |

For `BlochBasis` at W = 0 and `EigenBasis`, `states("exact")` returns
`_core.BlockLocalStates(num_cells, num_wann)`: `num_wann` state vectors that
together give the exact trace over the whole k mesh. `execute_runfile` then sets
`mu.number_of_states = basis.dimension` on the moments
(`runfile.py:1680–1681`). Do the same when you call `run_plan` yourself, as
example 17 does.

### `BlochMetricBasis`, `EigenMetricBasis`

```python
BlochMetricBasis(model, supercell, *, velocity_mode='berry_connection', route='lowdin', seed=12345)
EigenMetricBasis(bloch)
```

Clean non-orthogonal bases used by the SIESTA executors. `BlochMetricBasis`
builds S⁻¹ and S^(−1/2) per k block once.

| Parameter | Values | Meaning |
|---|---|---|
| `model` | — | A model with `Hk`, `Sk`, their derivatives, `k_mesh`, `tau`, `uc`, and `dipoles` for `covariant`. If the model has a method `spin_blocks(k_points, component)` that returns a matrix, `SX`, `SY`, `SZ` and the spin currents `J?S?` use it instead of the Pauli matrices on the spinor labels (see [Orthogonal models](#orthogonal-models-for-the-metric-executors)). If `model.spinful` is `False`, spin operators and spin currents raise `ValueError`. |
| `velocity_mode` | `bare`, `berry_connection`, `covariant` | How velocities are built. `covariant` raises `ValueError` if `model.dipoles` is `None`. |
| `route` | `lowdin`, `generalised`, `anticommutator` | Metric route. |
| `bloch` | a `BlochMetricBasis` | `EigenMetricBasis` reuses its cached eigenvectors. |

See the [`basis` and `metric` keys of the run file](run-file.md#basis).

### Other names

| Signature | Meaning |
|---|---|
| `make_anderson_disorder(dimension, disorder_w, seed=12345)` | Deterministic box disorder in [−W/2, W/2], complex array of length `dimension`. |
| `validate_w2s_csr_dir(csr_dir, seed_name, supercell, *, model_dir=None, operators=('VX',))` | Check a `wannier2sparse` directory. Raises `W2sCsrError` (a `ValueError` and `FileNotFoundError`). |
| `format_w2s_command(model_dir, seed_name, supercell, csr_dir, operators=('VX',))` | The `wannier2sparse` command that produces `csr_dir`. |
| `TightBindingModel` | `typing.Protocol` of the model interface. |

## `lsquant.moments`

### `run_plan`

```python
run_plan(plan, hamiltonian, window, states, operators=(), batch_size=None,
         symmetrise=True, dt=None, state_block=1, left_window=16, right_window=16)
    -> _core.ChebyshevMoments
```

Runs a compiled plan through `_core.MomentEngine` and fills the metadata of the
moments.

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `plan` | `_core.Plan` | — | Compiled plan, for example `words.compile(...).plan`. |
| `hamiltonian` | `_core.SparseMatrix` or `_core.LinearOperator` | — | Physical H; rescaled to `window` inside. |
| `window` | `_core.SpectralWindow` | — | Spectral window; must contain the spectrum. |
| `states` | `_core.StateGenerator` | — | Trace states. |
| `operators` | sequence of `_core.LinearOperator` | `()` | Operators in the order the plan refers to them (`CompiledWord.operators`). |
| `batch_size` | `int` or `None` | `None` | Chebyshev batch size. `None`: the moment count of the first `Expand` instruction. |
| `symmetrise` | `bool` | `True` | Hermitian-symmetrise rank-2 moments. |
| `dt` | `float` or `None` | `None` | Time step recorded on time-dependent moments. |
| `state_block` | `int` | `1` | Trace states advanced together. |
| `left_window` | `int` | `16` | Bra-side streaming window of the native rank-2 path. |
| `right_window` | `int` or `None` | `16` | Ket-side storage window of the native rank-2 path. `None` stores all M vectors (fastest, most memory). |

Returns the moments with `window`, `system_size`, `number_of_states`,
`state_label`, `word`, `time_step`, `kpm_weighted` and `symmetrised` set.
Rank-2 moments get the KPM weights (`_core.ApplyKpmWeights`).

Raises `ValueError` if `state_block`, `left_window` or `right_window` is below 1.
With `state_block > 1`, issues a `RuntimeWarning` when the estimated working set
exceeds 80 % of `MemAvailable` (`moments.py:568–582`).

### Plan builders

The word compiler calls these functions. Use them directly only if you build a
plan by hand.

| Signature | Plan |
|---|---|
| `dos_plan(num_moments, num_states, word, *, operator_index=None)` | Rank 1, `Tr[A delta(E-H)]`; `A` is omitted if `operator_index` is `None`. |
| `kubo_plan(num_moments, *, left_index, right_index, word)` | Rank 2, `Tr[O_left T_m O_right T_n]`. |
| `msd_plan(num_moments, num_steps, num_states, coeffs, *, velocity_index=0, phase=1+0j, word='[X, U(t)]^+ delta [X, U(t)]')` | Mean-square displacement. |
| `projected_evolution_plan(num_moments, num_steps, num_states, coeffs, *, operator_index, projector_index, phase=1+0j, word='P U(t)^+ O U(t) delta P')` | Time-evolved projected operator. |
| `time_correlation_plan(num_moments, num_steps, num_states, coeffs, *, left_index, right_index, phase=1+0j, word='A delta U(t)^+ B U(t)')` | Two-operator time correlation. |

### Helpers

| Signature | Returns |
|---|---|
| `HBAR` | ħ = 0.6582119624 eV·fs (`moments.py:22`). Time is in fs. |
| `time_step_from_tmax(num_steps, tmax) -> float` | `tmax / (num_steps - 1)`. Raises `ValueError` if `num_steps <= 1`. |
| `time_evolution_series(window, dt, tol=1e-15, use_legacy_bessel=False) -> (coeffs, phase)` | Bessel coefficients of U(dt) and the band-centre phase. |
| `window_from_gershgorin(hamiltonian, pad=0.1) -> _core.SpectralWindow` | Symmetric window [−ρ(1+pad), ρ(1+pad)], ρ = max row sum of \|H\|. |
| `estimate_engine_bytes(dim, plan, batch_size, state_block=1, left_window=16, right_window=None) -> int` | Estimated peak working set of the engine in bytes. |
| `choose_state_block(dim, num_moments, num_states, left_window=16, fraction=0.6, right_window=None) -> (int, str)` | `state_block` from available memory, at most 16, and a reason string. |

## `lsquant.reconstruct`

All functions take moments that carry `window` and `system_size` (as set by
`run_plan` or read with `_core.ReadChebmom`) and return a pair
`(x, y)` of 1D `float64` arrays. `broadening_meV` must be positive. Each
function first truncates the moments to M_eff = min(⌈π/η⌉, M), with
η = 2 · `broadening_meV` / 1000 / (window width). It then applies the Jackson
kernel, the only kernel the reconstruction uses. `alpha` (0 < `alpha` < 1,
default 0.95) keeps the uniform grid inside [−`alpha`, `alpha`] in rescaled
units. Units and normalisation are those of the curve files; see
[Output files](output-files.md#curve-files).

| Signature | Moments | Grid (points) | Returns |
|---|---|---|---|
| `dos_curve(moments, broadening_meV, alpha=0.95)` | rank 1 | uniform, 30 M_eff | E (eV), extensive DOS N·ρ(E) (1/eV). Divide by `moments.system_size` for the DOS per orbital. |
| `kubo_greenwood(moments, broadening_meV, alpha=0.95)` | rank 2 | uniform, M_eff − 1 | E (eV), Kubo–Greenwood σ_xx. |
| `kubo_bastin(moments, broadening_meV, alpha=0.95, *, backend='tables')` | rank 2 | uniform, 30 M_eff − 1 | E (eV), Kubo–Bastin σ_xy (cumulative integral). `backend='folded'` gives the same series with less memory. |
| `kubo_bastin_sea(moments, broadening_meV, alpha=0.95)` | rank 2 | uniform, 30 M_eff − 1 | E (eV), Fermi-sea term. |
| `kubo_bastin_surf(moments, broadening_meV, alpha=0.95)` | rank 2 | uniform, 30 M_eff − 1 | E (eV), Fermi-surface term. |
| `kubo_fft_grid(moments, broadening_meV, formula='bastin', *, backend=None)` | rank 2 | Chebyshev nodes, M_eff − 1 | E (eV), Bastin (`'bastin'`, `'bastin_fft'`, `True`) or Greenwood (`'greenwood'`, `'greenwood_fft'`, `False`). `backend`: `'tables'` or `'fft'`; `None` reads `LSQUANT_FFT_BACKEND`, default `'tables'`. |
| `time_correlation_curve(moments, broadening_meV, fermi_energy)` | rank 2 (M × N_t) | time steps | t (fs, `n · moments.time_step`), correlation at `fermi_energy` (eV). |
| `shc(moments, broadening_meV, num_cells, cell_area, *, backend='tables')` | rank 2 | as `kubo_bastin` | E (eV), `kubo_bastin` × 2π / (`num_cells` · `cell_area`), in e²/h for Pauli spin. |

Raises `ValueError` for the wrong rank, a non-positive broadening, `alpha`
outside (0, 1), an unknown `formula` or `backend`, or (for `shc`) a
non-positive `num_cells` or `cell_area`.

> **Known limitation.** M_eff − 1 points on the Chebyshev nodes (`kubo_fft_grid`)
> is a coarse grid: for `bastin_fft` with M = 128 the curve is off by up to 11 % away from
> plateaus, and ends at 0.64 instead of 0 above the bands. The end value is the
> quadrature error of the cumulative integral and falls as 1/M_eff (0.158 with
> 511 nodes); the largest deviation falls more slowly (5.9 % with 511 nodes).
> Prefer `kubo_bastin` on the uniform grid when you compare curves.

## `lsquant.kernels`

| Signature | Returns |
|---|---|
| `jackson_damping(num_moments) -> ndarray` | Jackson weights g_m, m = 0 … M − 1, with g_0 = 1. |
| `jackson_damping_2d(num_moments_0, num_moments_1) -> ndarray` | Outer product g(m₀, M₀) · g(m₁, M₁). |

The reconstruction functions of `lsquant.reconstruct` (above) apply these
weights.

## `lsquant.spectral_functions`

Symbolic Chebyshev coefficients of functions of H, built with SymPy. The
rescaled energy is x ∈ (−1, 1).

| Name | Signature or value | Meaning |
|---|---|---|
| `SpectralFunction` | `SpectralFunction(expr, x=None, m=None)` | Wraps a SymPy expression of `x` and the moment index `m`. Methods: `table(x_grid, num_moments)` returns the (M, N) table; `evaluate(x_values, param)`. |
| `DeltaFunction` | instance | cos(m arccos x) / (π √(1 − x²)). |
| `GaussianDelta` | instance | Gaussian of width η (here `m` is η). |
| `GreenRetarded` | instance | −i e^(−i m arccos x) / √(1 − x²). |
| `GreenRetardedDerivative` | instance | d/dx of `GreenRetarded`. |
| `TimeEvolutionFunction` | `TimeEvolutionFunction()` | exp(−i x t); `series_coefficients(x, tol=1e-15)` gives the Bessel coefficients. |
| `TimeEvolution` | instance | `TimeEvolutionFunction()`. |
| `bessel_series_coefficients(x, tol=1e-15)` | function | Jacobi–Anger coefficients with SciPy. |
| `legacy_bessel_coefficients(x, tol=1e-15)` | function | Same with GSL, through `_core.compat`. |
| `legacy_delta_table(x, num_moments)` | function | Delta table in the legacy order (golden byte compatibility). |
| `legacy_green_tables(x, num_moments_0, num_moments_1=None)` | function | `(greenR, DgreenR)` tables through `_core.compat`. |

## `lsquant.words`

The observable-word language. For the shapes, operators and formulas, see the
[observable-word reference](observable-words.md).

### `compile`

```python
compile(word, basis=None, params=None) -> CompiledWord
```

Parses a word, matches it against the catalogue and builds its plan.

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `word` | `str` or `Word` | — | Observable word, for example `"Tr[VX delta(E-H) VX delta(E-H)]"`. |
| `basis` | `str`, mapping or `None` | `None` | Basis name (`"real"`, `"bloch"`, `"eigen"`; `None` means `"real"`), or a mapping with `name`, `operators`, `operator_bases`, `metric`. |
| `params` | mapping or `None` | `None` | `expansion` (`moments`, default 64; `states`, default 1; `window`, a `SpectralWindow` or `[lower, upper]`, default [−1, 1]), `time` (`steps`, `tmax` or `dt`, `bessel`), `metric`, `reconstruction`. The `expansion` keys may also be given at the top level. |

`basis["operators"]` sets which operator names the word may use:

| `basis["operators"]` | Allowed names | Operator objects in the result |
|---|---|---|
| absent, or `basis` is a string or `None` | `DEFAULT_KNOWN_OPERATORS` (`words/compiler.py:28`) | none |
| list of names | `"1"` and those names | none |
| mapping name → `_core.LinearOperator` | `"1"` and those names | the mapped operators, in the order of the word |

Raises `ValueError` for a syntax error, a word that matches no shape, an unknown
basis or operator, or a non-Hermitian operator. It also raises `ValueError` for
a time-dependent word that lacks `time.steps`, or lacks both `time.tmax` and
`time.dt`, and for a formula outside those allowed for the shape.

### `CompiledWord`

| Field | Type | Meaning |
|---|---|---|
| `word` | `Word` | Parsed word. |
| `match` | `CatalogueMatch` | Matched shape. |
| `plan` | `_core.Plan` | Bytecode plan for `run_plan`. |
| `operator_names` | `tuple[str, ...]` | Operator names in plan order. |
| `operators` | `list[_core.LinearOperator]` | Operator objects in plan order (empty if `basis` gave no objects). |
| `reconstruction` | `ReconstructionSpec` | Default formula and allowed formulas. |
| `dt` | `float` or `None` | Time step of time-dependent words. |
| `to_dict()` | method | `word`, `catalogue`, `plan`, `operators`, `reconstruction`, `dt` as JSON-ready data. |

### `ReconstructionSpec`

| Field | Type | Default |
|---|---|---|
| `shape_id` | `str` | — |
| `default_reconstruction` | `str` | — |
| `formulas` | `tuple[str, ...]` | — |
| `kernel` | `str` | `"jackson"` |
| `broadening_meV` | `float` or `None` | `None` |
| `fermi_energies` | `tuple[float, ...]` | `()` |
| `num_points` | `int` or `None` | `None` |
| `use_fft_grid` | `bool` | `False` |
| `legacy_moment_prefix` | `str` | `""` |
| `legacy_moment_ext` | `str` | `".chebmom1D"` |
| `legacy_text_round_trip` | `bool` | `False` |

### Other names

| Signature | Meaning |
|---|---|
| `describe(word, basis=None, params=None) -> dict` | `compile(...).to_dict()`. With `params=None`, a time-dependent word compiles with a placeholder time axis. Behind `python -m lsquant words describe`. |
| `list_catalogue() -> list[dict]` | Every catalogue shape as a dictionary. Behind `python -m lsquant words`. |
| `parse(text) -> Word` | Parse a word into its syntax tree. Raises `ValueError`. |
| `tokenize(text) -> list[Token]` | Lexical tokens (`kind`, `value`, `pos`). |
| `match_catalogue(word) -> CatalogueMatch` | Match a parsed word to a shape. Raises `ValueError` with the nearest shape. |
| `nearest_catalogue_entry(word) -> CatalogueEntry` | Nearest shape, for error messages. |
| `validate_word_bases(word, calculation_basis=None, operator_bases=None) -> None` | Check that all operators are in one basis or explicitly transformed. |
| `format_plan_instructions(plan, operator_names=()) -> list[dict]` | Plan instructions as JSON-ready data. |
| `CATALOGUE`, `CATALOGUE_BY_ID` | Tuple of `CatalogueEntry`, and a dictionary by `shape_id`. |
| `Word` | Parsed word: `slots`, `sea`; properties `rank`, `is_time_dependent`; `canonical()`, `to_dict()`. |
| `OperatorSlot`, `SpectralSlot`, `CommutatorSlot` | Slot types of a `Word`. |
| `CatalogueEntry` | `shape_id`, `canonical_example`, `pattern_description`, `plan_builder`, `default_reconstruction`, `legacy_moment_prefix`, `legacy_curve_prefix`, `rank`, `time_dependent`. |
| `CatalogueMatch` | `entry`, `operators`, `engine_word`, `legacy_moment_prefix`, `legacy_moment_ext`, `default_formulas`. |

## `lsquant.models`

Model readers. `lsquant.models` exports the Wannier90 reader and the CSR and
basis-rotation helpers (`models/__init__.py`).

### `load_wannier_model`

```python
load_wannier_model(model_dir, seed=None) -> WannierModel
```

| Parameter | Type | Default | Meaning |
|---|---|---|---|
| `model_dir` | `str` | — | Directory with `<seed>_hr.dat`, `<seed>.uc`, `<seed>.xyz`, or the path of the `_hr.dat` file. |
| `seed` | `str` or `None` | `None` | File prefix. `None`: taken from the only `*_hr.dat` in the directory. |

Returns a `WannierModel` with `num_wann`, `uc` (3 × 3, Å, rows a₁ a₂ a₃),
`tau` (W × 3, Å), `r_int`, `r_cart`, `hoppings` (n_R × W × W, `1/ndegen`
folded in), `seed_name`, `model_dir`, and the methods `Hk(k_points)`,
`dHk(k_points, direction)`, `k_mesh(supercell)`, `bloch_phases(k_points)`,
`hamiltonian_csr(supercell)`, `velocity_csr(supercell, direction)`. A
`WannierModel` can be passed to `make_basis`.

```python
WannierModel.velocity_blocks(k_points, direction, mode="berry_connection") -> ndarray
```

Velocity blocks of shape (N_k, W, W) at Cartesian `k_points` (1/Å), in the
atom gauge of `Hk`, with the sign of `dHk`. A Hall word has two velocities,
so an overall sign of the velocity cancels.

| `mode` | Velocity | Extra file |
|---|---|---|
| `"bare"` | Derivative of the cell-gauge $H(\mathbf k)$, phases on $\mathbf R$ only. | — |
| `"berry_connection"` | Derivative of the atom-gauge $H(\mathbf k)$, phases on $\mathbf R + \tau_j - \tau_i$; the same matrix as `dHk`. Uses the Wannier centres only. | — |
| `"covariant"` | The `bare` derivative plus the commutator term $\propto [H, A]$, with $A$ the full position matrix $\langle 0m\vert\mathbf r\vert\mathbf Rn\rangle$, folded Hermitian with the degeneracies of `_hr.dat`. It adds the off-diagonal position elements between Wannier functions, which the other two modes take as zero. | `<seed>_r.dat` in `model_dir`, written by Wannier90 with `write_rmn = true`; read once per model. |

In a run file, `model.params.velocity_mode` builds the same three modes
through `models/orthogonal.py`, with the sign +∂H (see
[Wannier90 and PAOFLOW models](run-file.md#wannier90-and-paoflow-models-in-the-metric-executors)).
[Example 20](../../examples/20_pdse2_spin_hall_routes/README.md) uses the
covariant mode.

### Orthogonal models for the metric executors

`lsquant.models.orthogonal` (not re-exported by `lsquant.models`) gives a
Wannier90 or a PAOFLOW Hamiltonian the model protocol of `BlochMetricBasis`
with S = 1. The run-file executors use it for `model.kind = "paoflow"` and
for `model.kind = "wannier"` with `params.velocity_mode`, `spin` or
`fermi_energy`.

| Name | Meaning |
|---|---|
| `OrthogonalModel(uc, tau, r_int, hoppings, *, seed_name="model", positions=None, spin=None, metadata=None)` | Atom-gauge H(R) (the convention of `hk_blocks`). `positions = (r_int, A)` from `read_position_matrix` enables `covariant`, with the connection A − diag(τ). `spin` is `None` (Pauli matrices on interleaved spinor labels), a `(3, W, W)` on-site array, or an object with `blocks(k_points, component)`. Methods `Hk`, `Sk` (identity), `derivatives(k, axis, cell=False)` (+∂H, atom or cell derivative), `Ak(k, axis)`, `k_mesh`, `spin_blocks(k, component)`. |
| `read_wannier_spin(model_dir, seed, *, unitarity_tolerance=1e-6)` | Returns `(k_frac, S_W)`, S_W of shape (3, N_k, W, W) in Pauli units: V†S_B V from `<seed>.spn` (formatted or unformatted), with V = U_dis U from `<seed>_u.mat` and `<seed>_u_dis.mat`. For a disentangled model it reads `<seed>.eig` and the outer window of `<seed>.win`, because the rows of `_u_dis.mat` list only the bands inside the window, packed at the top. Raises `ValueError` if V†V differs from the identity. |
| `WannierSpin(model_dir, seed, uc, tau)` | Fourier transform of S_W(k) to S(R) on the centred R range of the Wannier90 mesh; `blocks(k_points, component)` returns the Hermitian atom-gauge matrix at Cartesian k. |
| `load_wannier_metric_model(model_dir, seed, *, spin="labels", positions=None, fermi_energy=0.0)` | `OrthogonalModel` of a Wannier90 run. `spin = "exact"` uses `WannierSpin`. `positions=None` reads `<seed>_r.dat` if present. `fermi_energy` (eV) is subtracted from H(R = 0). |
| `write_paoflow_model(paoflow, path)` | Call in a PAOFLOW script after `pao_hamiltonian()`. Writes an `.npz` with H(R) (stored as H(−R), because PAOFLOW sums with e^{−ik·R}; on an even FFT grid the Nyquist block is split between +n/2 and −n/2 to keep H(k) Hermitian), the orbital positions from PAOFLOW's `Dnm`, and twice PAOFLOW's on-site spin operator `Sj` (Pauli units). Energies stay measured from the DFT Fermi energy. |
| `load_paoflow_model(path)` | `OrthogonalModel` of that file, with its on-site spin and no position matrix. |
| `load_runfile_model(model, velocity_mode, *, phase=None)` | Model of a run-file `model` mapping for the metric executors (`siesta`, `wannier`, `paoflow` and the machine-learning kinds of [Machine-learning Hamiltonians](#machine-learning-hamiltonians)). |

### Machine-learning Hamiltonians

`lsquant.models.lcao` turns the atom-pair blocks of a machine-learning
Hamiltonian into a model of the metric executors; one module per code reads
its files. These modules are not re-exported by `lsquant.models`.

| Name | Module | Meaning |
|---|---|---|
| `load_deeph_model(directory, *, hamiltonian_file=None, fermi_energy=None)` | `lsquant.models.deeph` | Model of a DeepH-format directory: the processed-structure layout of DeepH-pack (before its JAX version), DeepH-E3, xDeepH and MACE-H (`hamiltonians_pred.h5`, `overlaps.h5`, `lat.dat`, `site_positions.dat`, `orbital_types.dat`, `info.json`), or the DeepH-pack layout (`POSCAR`, `info.json`, chunked `hamiltonian_pred.h5`, `overlap.h5`, optional `hkb.h5`). Needs `h5py`. `fermi_energy` (eV) defaults to the value of `info.json`. |
| `read_poscar(path)` | `lsquant.models.deeph` | `(lattice, elements, cartesian)` of a VASP 5 POSCAR, in Å. |
| `assemble_blocks(pairs, norb, *, spinful, name="H")` | `lsquant.models.lcao` | `(r_int, M)` from `(R, i, j, block)` items, where `block` couples atom `i` in the home cell with atom `j` in cell `R`; a spinful block is `[[uu, ud], [du, dd]]` and is stored with interleaved spinors. |
| `hermitian_pairs(r_int, mats)` | `lsquant.models.lcao` | Adds missing −R blocks as M(R)† and replaces each pair by (M(R) + M(−R)†)/2; returns the relative residual before the mean. |
| `build_lcao_model(uc, positions, norb, r_int, hamiltonian, overlap=None, *, spinful, seed_name, metadata, fermi_energy=0.0, overlap_r_int=None)` | `lsquant.models.lcao` | `OrthogonalModel` without `overlap`, `SiestaModel` with it; removes H − E_F S and sets the attribute `spinful`, which `BlochMetricBasis` reads to refuse spin operators on a spinless model. |
| `load_runfile_model(model, velocity_mode, *, phase=None)` | Model of a run-file `model` mapping for the metric executors (`siesta`, `openmx`, `wannier`, `paoflow`). |

`lsquant.models.openmx` (not re-exported by `lsquant.models`):

| Name | Meaning |
|---|---|
| `read_scfout(path)` | Read an OpenMX 3.9 scfout file (version 3) into a dictionary: atoms, orbitals per atom, neighbour pairs with lattice vectors, cell, positions, Hks, iHks, overlap, position matrix, chemical potential. Units of the file (Bohr, Hartree). |
| `load_openmx_model(path, *, dipoles=False)` | `SiestaModel` of a non-collinear spin–orbit scfout file: H − μS in eV, S lifted to interleaved spinors, and with `dipoles=True` the connection A = D − (τ_j + R − τ_i) S/2 from OpenMX's row-centred position matrix D. H, S and A are folded Hermitian (`metadata["hermitian_fold_relative"]`). |

### Other functions

| Name | Module | Meaning |
|---|---|---|
| `read_hr(path, wsvec_path=None)` | `models.wannier` | Read a `_hr.dat` file (dense or sparse); returns `(num_wann, r_int, hoppings)`. Reads `<seed>_wsvec.dat` next to it if present. |
| `read_uc(path)`, `read_xyz(path)` | `models.wannier` | Lattice vectors and orbital positions. |
| `reciprocal_lattice(uc)` | `models.wannier` | 2π (A⁻¹)ᵀ. |
| `k_mesh(uc, supercell)` | `models.wannier` | N₁N₂N₃ k-points, m₁ fastest. |
| `bloch_phases(k_points, tau)` | `models.wannier` | exp(−i k·τ), shape (N_k, W). |
| `hk_blocks(r_cart, tau, hoppings, k_points)`, `dhk_blocks(r_cart, tau, hoppings, k_points, direction)` | `models.wannier` | H(k) and ∂H/∂k blocks. |
| `supercell_csr(num_wann, r_int, hoppings, supercell)`, `supercell_velocity_csr(num_wann, r_int, r_cart, tau, hoppings, supercell, direction)` | `models.wannier` | Real-space supercell H and velocity as `_core.SparseMatrix`. |
| `read_csr_text(path, fold_hermitian=True, *, legacy_fold=False)`, `write_csr_text(path, matrix)` | `models.nonorthogonal` | Read and write LSQuant `.CSR` text files as `scipy.sparse.csr_matrix`. |
| `hermitian_fold(matrix)` | `models.nonorthogonal` | Hermitian part of a CSR matrix. |
| `make_local_transform(H, alpha=0.55, complex_entries=True, seed=12345)`, `rotate_basis(ops, R)` | `models.nonorthogonal` | Build a local transformation R and rotate operators to a non-orthogonal basis (returns `S = R†R` and `R† O R`). Used to make test fixtures. |
| `load_siesta_model(directory, label=None, *, dipoles=False, quadrature_order=24, quadrature_tolerance=1e-5, quadrature_max_order=384, phase=None, cache_dir=None)` | `models.siesta` (not re-exported by `lsquant.models`) | Read a SIESTA calculation (`.fdf` with `.TSHS` or `.HSX`, optional `.ion.xml`) with `sisl`. Needs the `sisl` package. `label=None` requires exactly one `.fdf` in the directory. See [Use a SIESTA calculation](../how-to/use-siesta.md). |

## See also

- [Core API](core-api.md): the compiled module `lsquant._core`.
- [Run-file reference](run-file.md) and [CLI reference](cli.md).
- [Observable-word reference](observable-words.md).
- [Example 17: LSQuant from Python](../../examples/17_python_api/README.md).
