# `lsquant` — Python API, Run-File Schema, and MCP Reference

For an overview of LSQUANT, installation instructions, the observable word
catalogue, and the worked examples `01`–`08`, see the main repository
[**README.md**](../../README.md). This page documents the Python package API,
the JSON run-file schema accepted by `python -m lsquant run`, and the MCP server
surface.

---

## Python API Summary

- **High-level execution** (`lsquant.runfile`):
  - `execute_runfile(runfile_path_or_dict, *, out_dir_override=None)` $\to$ `dict` (`output_dir`, `resolved_runfile`, `moment_files`, `curve_files`, `written_files`)
  - `reconstruct_target(target_path, *, formula=None, broadening_meV=None, fermi_energies=None, out_dir=None)`
  - `inspect_target(path)`
- **Bases** (`lsquant.bases`, `lsquant.disorder_bases`):
  - `RealBasis(model, supercell=(1, 1, 1), *, disorder_w=0.0, seed=12345, disorder_seed=7, csr_dir=None)`
  - `BlochBasis(model, supercell=(1, 1, 1), *, disorder_w=0.0, seed=12345, disorder_seed=7)`
  - `EigenBasis(model, supercell=(1, 1, 1))`
  - `BlochMetricBasis(model, supercell, *, route="lowdin", applier="exact", velocity_mode="covariant", ...)`
  - `EigenMetricBasis(model, supercell, *, velocity_mode="covariant")`
  - `DisorderedBlochMetricBasis` and `RealSpaceMetricBasis`, the disordered SIESTA layouts, take `disorder_kind="scalar"` (`"scalar"`, `"magnetic"` or `"independent"`). A model that defines `physical_spin` accepts `"independent"` only.
  - `spinor_disorder_values(dimension, width, seed, kind="scalar")`: the disorder values of an interleaved spinor basis.
  - `make_basis(basis_name, model, supercell=(1, 1, 1), ...)`
- **Models** (`lsquant.models`):
  - `load_wannier_model(directory, seed)` — reads Wannier90 `<seed>_hr.dat`, `.uc`, and `.xyz`.
  - `load_siesta_model(directory, label=None, *, velocity_mode="covariant", ...)` — reads SIESTA `.fdf` + `.TSHS`/`.HSX` and `.ion.xml` via `sisl`, with certified Gauss–Legendre midpoint quadrature for orbital dipole matrices $A_\alpha(\mathbf{R})$.
- **Words & Compilation** (`lsquant.words`):
  - `lsquant.words.parse(word_str)` $\to$ `Word` AST
  - `lsquant.words.compile(word, basis=None, params=None)` $\to$ `CompiledWord` (`plan`, `operator_names`, `operators`, `reconstruction`, `dt`)
  - `lsquant.words.describe(word, basis=None, params=None)` $\to$ `dict`
  - `lsquant.words.list_catalogue()` $\to$ `list[dict]`
- **Engine Execution & Reconstruction** (`lsquant.moments`, `lsquant.reconstruct`):
  - `dos_plan(num_moments, num_states, word, *, operator_index=None)` $\to$ `_core.Plan`
  - `kubo_plan(num_moments, *, left_index, right_index, word)` $\to$ `_core.Plan`
  - `run_plan(plan, hamiltonian, window, states, operators=(), batch_size=None, symmetrise=True, dt=None)` $\to$ `_core.ChebyshevMoments`
  - `dos_curve(moments, broadening_meV, alpha=0.95)`
  - `kubo_greenwood(moments, broadening_meV, alpha=0.95)`
  - `kubo_bastin(moments, broadening_meV, alpha=0.95, *, backend="tables")` (plus `kubo_bastin_sea`, `kubo_bastin_surf`, `kubo_fft_grid`)
  - `shc(moments, broadening_meV, num_cells, cell_area, *, backend="tables")`
  - `time_correlation_curve(moments, broadening_meV, fermi_energy)`
- **Spectral Projection & Downfolding** (`lsquant.projection`):
  - `spectral_projection(hamiltonian, overlap, window, *, kind="shift")` (`kind` in `"hard"` | `"shift"`)
  - `downfold_orthogonal(hamiltonian, overlap, window)` and `LocalMatrix`
  - `compare_windows(source, words, windows, energies, *, cell_area, ...)`

---

## Command-Line Interface (`python -m lsquant`)

```bash
python3 -m lsquant words
python3 -m lsquant words describe "Tr[VX delta(E-H) VY dG-(E)/dE]"
python3 -m lsquant run <run.json> [--out <output_dir>] [--json]
python3 -m lsquant reconstruct <moment_file_or_dir> [--formula <name>] [--broadening-meV <val>] [--fermi-energies <E...>] [--out <dir>] [--json]
python3 -m lsquant inspect <run_json_or_moment_file_or_dir>
```

---

## Run-File Schema (`python/lsquant/runfile.py`)

Every calculation is configured by a strictly validated JSON run file. Unknown keys
at the top level or inside any sub-object are rejected before any operator file is
loaded.

See the maintained [run-file reference](../../docs/reference/run-file.md) for
all accepted keys, defaults, model sources and executor restrictions. The main
sections select the model, basis, metric, observable word(s), expansion, time
grid, reconstruction and output; `numerics` selects engine options.

The metric routes are `lowdin`, `generalised`, and the approximate
`anticommutator` route, available only in the clean SIESTA executor. A `time`
section requires `steps > 1` and at least one of `tmax` or `dt` (in fs); `dt`
takes precedence when both are given. There is no default `tmax`.

`reconstruction.num_points` sets only the default SIESTA oracle energy grid
(1001 points when null), including models using that executor. Explicit
`energies` or `energy_grid` override it. It does not set KPM curve grids.
`reconstruction.fft_grid` does not select a grid: the formula names
`greenwood_fft` and `bastin_fft` select the Chebyshev-node grid. The flag is
recorded and checked for incompatibility with `bastin_backend = "folded"`.
Validation emits one `RuntimeWarning` when `fft_grid` is true without a
requested `*_fft` formula, or `num_points` is set for a model outside the
metric executors, where it has no effect. Both keys remain accepted.

---

## MCP Server (`python -m lsquant.mcp_server`)

`lsquant.mcp_server` exposes five tools over `stdio` using `FastMCP` (`pip install mcp`):

- `words_describe(word)` — parses and compiles `word`, returning the same JSON structure as `python -m lsquant words describe "<word>"`.
- `run(runfile_json)` — validates `runfile_json` and dispatches `python -m lsquant run` in a background subprocess, returning `{"run_id", "log", "output_dir"}`.
- `status(run_id)` — polls a background calculation started by `run`.
- `reconstruct(output_dir, formula, broadening_meV, fermi_energies)` — reconstructs curves from moment files in `output_dir`.
- `inspect(path)` — inspects a moment file, run file, or output directory.

All file operations are restricted to `LSQUANT_MCP_ROOT` (defaults to the repository root or `--root <dir>`).

```json
{
  "mcpServers": {
    "lsquant": {
      "command": "python3",
      "args": ["-m", "lsquant.mcp_server"],
      "env": {
        "PYTHONPATH": "/path/to/lsquant/build/golden/lsquant_core:/path/to/lsquant/python",
        "LSQUANT_MCP_ROOT": "/path/to/lsquant",
        "OMP_NUM_THREADS": "16"
      }
    }
  }
}
```
