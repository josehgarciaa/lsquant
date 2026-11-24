# Run-file reference

A run file is the JSON file that defines one LSQuant calculation: the model,
the basis, the metric, the observable word, the Chebyshev expansion, the
reconstruction, and the output. You pass it to
[`python -m lsquant run`](cli.md#run).

This page lists every key the validator accepts. The source of truth is
`python/lsquant/runfile.py` (`validate_runfile_dict`, `RunFile.from_dict`,
`execute_runfile`), together with the executors it dispatches to.

## How LSQuant reads a run file

1. `validate_runfile_dict` checks the whole file before any model file is
   read. An unknown key at the top level or inside any section stops the run
   with `Unknown key '<key>' in run file (<section>); allowed keys are [...]`.
2. `RunFile.from_dict` fills in every default. The **Default** columns on this
   page show the value after this step.
3. `execute_runfile` writes `run.resolved.json` and runs the executor for
   `model.kind`. The SIESTA executors apply further checks at this stage; see
   [SIESTA executors](#siesta-executors).

To print the resolved run file and stop before any computation, run
[`python -m lsquant inspect run.json`](cli.md#inspect).

A key set to `null` is treated as absent: the default applies.

### Path resolution

`model.dir`, `model.csr_dir`, `model.file`, `model.hr_file`, and `model.positions` and
`model.atom_map` when they are strings, are resolved in this order:

1. An absolute path (after `~` expansion) is used as is.
2. A relative path is tried against the repository root (the directory two
   levels above `python/lsquant/runfile.py`).
3. If that does not exist, it is tried against the directory of the run file.
4. If neither exists, the repository-root path is used, and the executor
   reports the missing file.

`output.dir` follows its own rule: it is relative to the current working
directory.

For MCP runs, the configured server root replaces the repository root and is
also the subprocess working directory. After resolution, the server checks all
explicit model paths, generator scripts and `output.dir` with symbolic links
resolved. ABACUS's `model.params.overlap_file`, `stru_file`, `orbital_dir` and
`position_file` follow the input resolution above. GPAW's `dipole_cache` is
relative to the server root and is checked for reads and writes; DeepH's
`hamiltonian_file` is checked relative to `model.dir`. Any escape is refused
before dispatch; the error names the field and its supplied value without
disclosing the resolved outside path. The default `output.dir: "."` remains
allowed and writes into the root itself. These restrictions apply only to MCP;
the command line can still read paths outside the repository.

## Top-level keys

| Key | Type | Required | Meaning |
|---|---|---|---|
| `model` | object | yes | Where the Hamiltonian and operators come from. See [`model`](#model). |
| `word` | string | yes | The observable word. See [`word` and `words`](#word-and-words). |
| `expansion` | object | yes | Chebyshev expansion settings. See [`expansion`](#expansion). |
| `words` | list of strings | no | Several words evaluated in one run. |
| `basis` | object | no | Representation of H. See [`basis`](#basis). |
| `metric` | object | no | Treatment of the overlap matrix S. See [`metric`](#metric). |
| `projection` | object | no | Spectral projection of a SIESTA model. See [`projection`](#projection). |
| `time` | object | conditional | Time grid. Required for time-dependent words. See [`time`](#time). |
| `reconstruction` | object | no | How moments become curves. See [`reconstruction`](#reconstruction). |
| `output` | object | no | Output directory and file layout. See [`output`](#output). |
| `numerics` | object | no | Engine back end and operator layout (native build only). See [`numerics`](#numerics). |
| `reference` | boolean | no | SIESTA only. Forces or suppresses the exact eigenbasis reference. See [Provenance keys](#provenance-keys). |
| `build_profile` | string | no | Must equal the build profile of the loaded engine. See [Provenance keys](#provenance-keys). |
| `receipt` | any | no | Accepted and ignored on input. SIESTA executors write it into `run.resolved.json`. |
| `runtime` | any | no | Accepted and ignored on input. Every run writes it into `run.resolved.json`. |

## `model`

| Key | Type | Allowed values | Default | Used by | Meaning |
|---|---|---|---|---|---|
| `kind` | string | `"csr"`, `"wannier"`, `"siesta"`, `"openmx"`, `"pyscf"`, `"paoflow"`, `"generator"`, `"deeph"`, `"deeptb"`, `"hamgnn"`, `"threebodytb"`, `"abacus"`, `"gpaw"` | `"csr"` | all | Input source. |
| `dir` | path | — | `null` | `csr`, `wannier`, `siesta`, `deeph`, `deeptb`, `generator` (tutorials 3 and 4) | Directory that holds the model files. |
| `label` | string | — | `null` | `csr`, `siesta` | File prefix: `<label>.HAM.CSR` (`csr`) or `<label>.fdf` (`siesta`). |
| `seed` | string | — | `null` | `wannier` | Wannier90 seed: LSQuant reads `<seed>_hr.dat`, `<seed>.uc` and `<seed>.xyz`. |
| `file` | path | — | `null` | `openmx`, `pyscf`, `paoflow`, `hamgnn`, `threebodytb`, `abacus`, `gpaw` | `openmx`: the `<System.Name>.scfout` file. `pyscf`: the model file written by `write_pyscf_model` (`.npz`). `paoflow`: the model file written by `write_paoflow_model` (`.npz`). `hamgnn`: the `.npz` of `write_hamgnn_model`. `threebodytb`: the XML of ThreeBodyTB's `write_tb_crys` (`.xml`, `.xml.gz`). `abacus`: H(R) in the ABACUS sparse format. `gpaw`: the `.npz` of `export_gpaw_lcao` (script `models/gpaw.py`). |
| `operators` | list of strings | operator names | `null` | `csr` | Operators to load as `<label>.<OP>.CSR`. If absent, every `<label>.*.CSR` except `HAM` and `S` is loaded. A listed operator without a file raises `FileNotFoundError` naming the missing file, for every CSR backend. |
| `csr_dir` | path | — | `null` | `wannier`, `basis.name = "real"` | Directory with a prebuilt supercell `<seed>.HAM.CSR` and its `wannier2sparse` receipt `<seed>.out`. LSQuant checks it against `seed` and `basis.supercell` before it loads anything. |
| `name` | string | see [Generator models](#generator-models) | `null` | `generator` | Generator from the catalogue. |
| `tutorial` | integer | `1`–`6` | from `name` | `generator` | Tutorial driver to run. Not validated; an unknown value fails at execution with `Unsupported generator`. |
| `params` | object | depends on `kind` | `{}` | `siesta`, `openmx`, `pyscf`, `wannier`, `paoflow`, `deeph`, `deeptb`, `hamgnn`, `threebodytb`, `abacus`, `gpaw`, `generator` | Model parameters. See [SIESTA models](#siesta-models), [OpenMX models](#openmx-models), [PySCF models](#pyscf-models), [Wannier90 and PAOFLOW models in the metric executors](#wannier90-and-paoflow-models-in-the-metric-executors), [Machine-learning Hamiltonians in the metric executors](#machine-learning-hamiltonians-in-the-metric-executors) and [Generator models](#generator-models). For `csr` the object is not checked and not used; for `wannier` it is checked only when it holds `velocity_mode`, `spin` or `fermi_energy`. |
| `positions` | path or array | file of `x y z` lines, `.xyz`-like file, or `(N, 2)`/`(N, 3)` array | `null` | `csr` with `atom_block` | Orbital positions. |
| `orbitals_per_atom` | integer | positive | `null` | `csr` with `atom_block` | Uniform number of orbitals per atom. Give this or `atom_map`, not both. Required with `hr_file`. |
| `atom_map` | path or array | one integer per orbital | `null` | `csr` with `atom_block` | Atom index of every orbital. |
| `hr_file` | path | — | `null` | `csr` with `atom_block` | Hopping list `H(R)`: header `W nhops`, then rows `r1 r2 i j re im`. Replaces `dir` and `label`. |
| `supercell` | list of 2 integers | `[L1, L2]`, both ≥ 2 | `null` | `csr` with `hr_file` | Supercell built from `hr_file`. |

### Requirements by kind

| `model.kind` | Required keys | Validator error if missing |
|---|---|---|
| `csr` | `dir` and `label`, unless `hr_file` is set | `Run file with model.kind='csr' requires 'dir' and 'label'` |
| `csr` with `hr_file` | `numerics.operator = "atom_block"` | `model.hr_file requires numerics.operator='atom_block'` |
| `wannier` | `dir` and `seed` | `Run file with model.kind='wannier' requires 'dir' and 'seed'` |
| `siesta` | `dir` and `label` | `SIESTA model requires dir and label` |
| `openmx` | `file` | `Run file with model.kind='openmx' requires 'file' (the .scfout file)` |
| `pyscf` | `file` | `Run file with model.kind='pyscf' requires 'file' (see write_pyscf_model)` |
| `paoflow` | `file` | `Run file with model.kind='paoflow' requires 'file' (see write_paoflow_model)` |
| `deeph` | `dir` | `Run file with model.kind='deeph' requires 'dir'` |
| `deeptb` | `dir` | `Run file with model.kind='deeptb' requires 'dir'` |
| `hamgnn` | `file` | `Run file with model.kind='hamgnn' requires 'file'` |
| `threebodytb` | `file` | `Run file with model.kind='threebodytb' requires 'file'` |
| `abacus` | `file`, `params.overlap_file`, `params.stru_file` | `Run file with model.kind='abacus' requires 'file'`; `model.params.overlap_file is required for model.kind='abacus' ...` |
| `gpaw` | `file` | `Run file with model.kind='gpaw' requires 'file'` |
| `generator` | `name` | `Run file with model.kind='generator' requires 'name' from [...]` |

For `csr`, the executor also needs `<dir>/<label>.HAM.CSR` and, when
`metric.policy` is `"nonorthogonal"`, `<dir>/<label>.S.CSR`.

### SIESTA models

`model.kind = "siesta"` reads `<dir>/<label>.fdf` and `<dir>/<label>.TSHS`
(or `<label>.HSX`) with `sisl`. `models/siesta.py` accepts only a
non-orthogonal Hamiltonian with spin–orbit coupling. For any other Hamiltonian
it stops with
`load_siesta_model requires a nonorthogonal SIESTA SOC Hamiltonian`. With
`velocity_mode = "covariant"` it also reads `<dir>/<species>.ion.xml` for every
species, to integrate the orbital dipole matrices.

LSQuant uses the energies as stored in the TSHS file, where E_F = 0. It records
the printed Fermi level in `basis.json` and never subtracts it.

Allowed keys in `model.params`:

| Key | Type | Allowed values | Default | Meaning |
|---|---|---|---|---|
| `velocity_mode` | string | `"bare"`, `"berry_connection"`, `"covariant"` | `"berry_connection"` (clean and projected runs); `"covariant"` (disordered and real-space runs) | How velocity operators are built from H, S and the dipoles. See [Non-orthogonal bases](../explanation/nonorthogonal-bases.md). |
| `oracle` | boolean | — | `false` | Clean runs: also compute the exact per-k spin Hall oracle (`exact_shc.json`, `spin_chern.json`). Disordered runs refuse `true`. |
| `oracle_degeneracy_tolerance` | number | positive, finite | `1e-6` | Degeneracy tolerance of the oracle, in eV. |
| `bloch_memory_mb` | number | positive, finite | `256.0` | Workspace limit of the Bloch sums, in MiB (clean runs). It sets how many k points are summed at once; a model too large for one k point within the limit is summed one k point at a time. |
| `quadrature_order` | integer | ≥ 2 | `24` | Initial Gauss–Legendre order of the dipole quadrature. |
| `quadrature_max_order` | integer | ≥ 2 × `quadrature_order` | `384` | Maximum quadrature order. |
| `quadrature_tolerance` | number | positive | `1e-5` | Refinement tolerance of the dipole quadrature. |
| `nonorthogonal_disorder` | boolean | — | `false` | Sends a clean Bloch run to the disordered executor. The disordered executor sets it to `true` in `run.resolved.json`. |

### OpenMX models

`model.kind = "openmx"` reads the `<System.Name>.scfout` file that OpenMX 3.9
writes with `HS.fileout on` (`models/openmx.py`). It needs a non-collinear
calculation with spin–orbit coupling and stops otherwise. The model is
non-orthogonal and runs through the clean and projection executors of SIESTA
with their rules (`metric.policy = "nonorthogonal"`, the same words). Only
`basis.name` `bloch` or `eigen` with `basis.disorder.W = 0` is accepted: the
disordered executor and the `real` basis are SIESTA only.
Energies are measured from the chemical potential stored in the file, so
E_F = 0. `velocity_mode = "covariant"` uses the position matrix of the basis
orbitals written in the same file. H, S and that position matrix are folded
Hermitian; the removed relative parts are recorded in the model metadata
(`hermitian_fold_relative`, copied to `basis.json` by the clean executor).

Allowed keys in `model.params`: `velocity_mode` (`"bare"`,
`"berry_connection"`, `"covariant"`; default `"berry_connection"`), `oracle`,
`oracle_degeneracy_tolerance`, `bloch_memory_mb`. `spin` and `fermi_energy`
are refused.

### PySCF models

`model.kind = "pyscf"` reads the `.npz` file that `write_pyscf_model`
(`models/pyscf_lcao.py`) writes from a converged periodic two-component PySCF
calculation (`pyscf.pbc.dft.KGKS`, for example with `.x2c1e()`). The model is
non-orthogonal and runs like an OpenMX model: `metric.policy =
"nonorthogonal"`, `basis.name` `bloch` or `eigen`, `basis.disorder.W = 0`.
Energies are measured from the `fermi_energy` stored in the file (by default
the middle of the gap). `velocity_mode = "covariant"` uses the analytic
position matrix of the Gaussian basis.

Allowed keys in `model.params`: `velocity_mode` (`"bare"`,
`"berry_connection"`, `"covariant"`; default `"berry_connection"`), `oracle`,
`oracle_degeneracy_tolerance`, `bloch_memory_mb`. `spin` and `fermi_energy`
are refused.

### Wannier90 and PAOFLOW models in the metric executors

A `wannier` model whose `params` hold `velocity_mode`, `spin` or
`fermi_energy`, and every `paoflow` model, run through the clean or the
projection executor of SIESTA (see [SIESTA executors](#siesta-executors)).
`models/orthogonal.py` gives them the model protocol of `BlochMetricBasis`
with S = 1. They therefore get the three velocity modes, a physical spin
operator in `SX`, `SY`, `SZ` and `JXSZ`, spectral projection, the exact
oracle and the KPM words of a SIESTA model. A `wannier` model without these
keys keeps the original Wannier executor and its `real` basis.

These runs need `basis.name` `bloch` or `eigen` and `basis.disorder.W = 0`.
`metric.policy = "orthogonal"` is accepted and run as the exact Löwdin route,
which is the same with S = 1.

| Key | Kinds | Allowed values | Default | Meaning |
|---|---|---|---|---|
| `velocity_mode` | `wannier`, `paoflow` | `"bare"`, `"berry_connection"`, `"covariant"` (`wannier` only) | `"berry_connection"` | `bare`: cell derivative of H. `berry_connection`: derivative in the atom gauge, H(k) = Σ_R e^{ik·(R+τ_j−τ_i)} H(R). `covariant`: ∂H_cell + i[H, A], with the full position matrix A of `<seed>_r.dat`. Velocities have the sign +∂H; `WannierModel.velocity_blocks` returns −∂H, and a Hall word gives the same result with either sign. |
| `spin` | `wannier` | `"labels"`, `"exact"` | `"labels"` | `labels`: Pauli matrices on the spinor labels (Wannier90 orders each spinor projection up, then down). `exact`: S_W(k) = V†S_B V from `<seed>.spn`, with V = U_dis U from `<seed>_u.mat` and `<seed>_u_dis.mat`. With spin–orbit coupling the Wannier functions mix spin, so only `exact` is the spin operator. |
| `fermi_energy` | `wannier` | number (eV) | `0.0` | Subtracted from H(R = 0). Energies, `projection.window` and `reconstruction.fermi_energies` are then measured from it. |
| `oracle`, `oracle_degeneracy_tolerance`, `bloch_memory_mb` | both | as for SIESTA | as for SIESTA | As in [SIESTA models](#siesta-models). |

For `spin = "exact"` on a disentangled model LSQuant also reads `<seed>.eig`
and `dis_win_min`/`dis_win_max` from `<seed>.win`: the rows of `_u_dis.mat`
list only the bands inside the outer window, packed at the top. It stops with
`V^+V differs from the identity` when the window does not match the files.

A `paoflow` file holds H(R), the orbital positions and PAOFLOW's on-site spin
operator (in Pauli units), with energies measured from the DFT Fermi energy.
Write it in the PAOFLOW script, after `pao_hamiltonian()`:

```python
from lsquant.models.orthogonal import write_paoflow_model
write_paoflow_model(paoflow, "pdse2_paoflow.npz")
```

PAOFLOW has no position matrix, so `velocity_mode = "covariant"` is refused;
`berry_connection` is the velocity of PAOFLOW's `gradient_and_momenta`. `spin`
and `fermi_energy` are refused for `paoflow`.

### Machine-learning Hamiltonians in the metric executors

Machine-learning models of the DFT Hamiltonian run through the clean or the
projection executor of SIESTA, like `paoflow`. Each kind has a reader in
`python/lsquant/models/` that converts the Fourier sign, the units (eV, Å),
the orbital positions and the spinor order of its code. A non-orthogonal
model (with S(R)) takes the exact Löwdin route; `metric.policy` can be left
out. These runs need `basis.name` `bloch` or `eigen` and
`basis.disorder.W = 0`.

| `kind` | Codes | Path key | Reader |
|---|---|---|---|
| `deeph` | DeepH-pack, DeepH-E3, xDeepH, MACE-H | `dir`: the directory of one structure with the prediction, the overlap and the structure files | `models/deeph.py`, see [Use a DeepH-format model](../how-to/use-deeph.md) |
| `deeptb` | DeePTB-SK, DeePTB-E3 (spinless) | `dir`: the directory written by `TBSystem.to_pardiso_json` or `dptb pdso` | `models/deeptb.py`, see [Use a DeePTB model](../how-to/use-deeptb.md) |
| `hamgnn` | HamGNN (OpenMX, ABACUS and SIESTA/HONPAS bases) | `file`: the `.npz` of `write_hamgnn_model` (script `models/hamgnn.py`) | `models/hamgnn.py`, see [Use a HamGNN model](../how-to/use-hamgnn.md) |
| `threebodytb` | ThreeBodyTB (TB3) | `file`: the XML of `write_tb_crys` (`.xml` or `.xml.gz`) | `models/threebodytb.py`, see [Use a ThreeBodyTB model](../how-to/use-threebodytb.md) |
| `abacus` | ABACUS (LCAO), NextHAM | `file`: H(R) in the ABACUS sparse format (ABACUS `data-HR-sparse_SPIN0.csr` or `hrs1_nao.csr`, NextHAM `hr_tot`), text or binary | `models/abacus.py`, see [Use an ABACUS calculation](../how-to/use-abacus.md) and [Use a NextHAM model](../how-to/use-nextham.md) |
| `gpaw` | GPAW, LCAO mode (not machine learning: the same input, atomic-orbital H(R) and S(R); the position matrix is built from the radial functions of the export) | `file`: the `.npz` of `export_gpaw_lcao` (script `models/gpaw.py`) | `models/gpaw.py`, see [Use a GPAW calculation](../how-to/use-gpaw.md) |

| Key | Kinds | Allowed values | Default | Meaning |
|---|---|---|---|---|
| `velocity_mode` | all | `"bare"`, `"berry_connection"`; for `abacus` with `position_file` and for `gpaw` also `"covariant"` | `"berry_connection"` | As for SIESTA. Without a position matrix `"covariant"` is refused. |
| `fermi_energy` | all | number (eV) | `deeph`: `fermi_level` or `fermi_energy_eV` of `info.json`, else `0.0`; `deeptb`: `0.0`; `hamgnn`: `0.0`; `threebodytb`: `efermi` of the file; `abacus`: `0.0`; `gpaw`: the Fermi level of GPAW in the file | Removed as H(R) − E_F S(R). Energies, `projection.window` and `reconstruction.fermi_energies` are then measured from it. |
| `hamiltonian_file` | `deeph` | file name in `dir` | `hamiltonians_pred.h5` or `hamiltonian_pred.h5`, then the DFT file `hamiltonians.h5` or `hamiltonian.h5` | The H(R) file to read. |
| `orbitals_per_element` | `threebodytb` | object, element → positive integer | ThreeBodyTB's table (`Atomdata.jl`) | Orbitals of one atom of an element; for elements missing from the table or a modified database. |
| `overlap_file`, `stru_file` | `abacus` | paths | required | S(R) in the ABACUS sparse format and the ABACUS `STRU` of the same structure. Resolved like `model.file`. |
| `orbital_dir` | `abacus` | path | directory of `stru_file` | Where the numerical orbital files named in `STRU` are. They give the orbitals per atom and their angular momenta. |
| `orbitals_per_element` | `abacus` | object, element → positive integer | from the orbital files | Orbitals of one atom of an element, when the orbital files are missing. |
| `position_file` | `abacus` | path | none | ABACUS's position matrix `data-rR-sparse.csr` (`out_mat_r 1`), text or binary, of the same run. Gives the connection A = Q − (τ_i + τ_j + R) S/2 of `"covariant"`; read only for that velocity. The removed anti-Hermitian part goes to `basis.json` as `connection_hermiticity_residual`. Resolved like `model.file`. |
| `soc` | `gpaw` | `true`, `false` | `true` | Add GPAW's spin–orbit term (spinful model) or not (spinless). |
| `velocity_mode` | `gpaw` | also `"covariant"` | `"berry_connection"` | The position matrix is integrated from the radial functions of the export. |
| `quadrature_order`, `dipole_cache` | `gpaw` | integer; path | `48`; none | Order of that integration; file that stores it for reuse. |
| `oracle`, `oracle_degeneracy_tolerance`, `bloch_memory_mb` | all | as for SIESTA | as for SIESTA | As in [SIESTA models](#siesta-models). |

LSQuant makes H(R) Hermitian as the mean of H(R) and H(−R)† and writes the
largest difference, relative to the largest element, to `basis.json` as
`hamiltonian_hermiticity_residual`; above 10⁻³ it warns. A model without
spin–orbit coupling has no spin operator: the words with `SX`, `SY`, `SZ` or
`JXSZ` stop with `... needs spin, but the model is spinless`.

### Generator models

Only the executable tutorial generators listed below are accepted. The former
`gen_models` catalogue entry has been removed because it had no executor.

`model.kind = "generator"` runs a tutorial driver that builds its own models.
The validator accepts only `basis.name = "real"` and
`metric.policy = "orthogonal"`. It also requires `basis.disorder` to be absent
or `null`; otherwise it stops with `pass disorder parameters in model.params`.
Generator receipts omit the inert default `basis.disorder`, so
`run.resolved.json` can be passed to `run` again.
The drivers decide which files are written, and they ignore most `output` and
`reconstruction` keys.

| `name` | Tutorial | Allowed `params` (default) | Other keys the driver reads |
|---|---|---|---|
| `make_chain` | 1 | `sizes` (`[256, 1024, 4096]`), `hopping_t` (`-1.0`) | `expansion.moments`, `expansion.window`, `word`, `reconstruction.broadening_meV` (default 20) |
| `make_disordered_chain` | 2 | `sizes` (`[256, 1024, 4096]`), `W` (`1.0`), `realisations` (`16`), `readme_size` (`512`), `readme_disorder_seed` (`1`) | `expansion.seed`, `moments`, `states` (default 1), `window`; `time.steps`, `tmax`, `bessel`; `reconstruction.broadening_meV` (default 10), `fermi_energies` (default `[0, 1, -1, 1.9]`), `legacy_text_round_trip` |
| `make_haldane` | 3 | none | `model.dir` (operator directory), `expansion.moments`, `reconstruction.broadening_meV` (default 10) |
| `make_chain_mag` | 4 | none | `model.dir` (operator directory) |
| `make_square` | 5 | `L` (`200`), `W` (`[0.5, 1.0, 2.0, 4.0]`), `disorder_seed` (`1`), `dos_moments` (`1024`) | `expansion.seed`, `moments`, `states` (default 1), `window`; `time`; `words`; `reconstruction.broadening_meV` (default 80), `fermi_energies` (default `[2.0]`), `legacy_text_round_trip` |
| `make_ey`, `make_rashba` | 6 | `ey_l` (`20`), `ey_ds` (`[0.05, 0.07, 0.10]`), `dp_l` (`20`), `dp_lambda` (`0.2`), `dp_ws` (`[0.0, 6.0]`), `dp_moments` (`160`), `dp_steps` (`120`), `dp_tmax` (`45.0`), `ratelaw_l` (`50`), `ratelaw_ds` (`[0.04, 0.06, 0.08, 0.10, 0.14]`), `ratelaw_nr` (`20`), `crosscheck_l` (`8`), `crosscheck_d` (`0.10`), `crosscheck_nr` (`20`), `crosscheck_moments` (`96`), `crosscheck_steps` (`60`), `crosscheck_tmax` (`220.0`), `disorder_seed` (`1`) | `expansion`, `time`, `words`, `reconstruction.kernel` (only `"jackson"`), `broadening_meV` (default 140), `fermi_energies` (default `[2.0]`), `formulas`, `legacy_text_round_trip`, `output.capture_tutorial` |

`RunFile` resolves `time.bessel`, `legacy_text_round_trip` and
`capture_tutorial` to `"scipy"`, `false` and `false` before the driver starts.
The driver defaults (`"legacy"`, `true` and `true`) therefore never apply. Set
these keys explicitly, as the example run files do.

The defaults printed by `inspect` match a minimal resolved run: in particular,
`model.operators` is `null`, `expansion.left_window` is `32`, and for the
64-moment example `expansion.batch` is `64`. When omitted, `batch` follows the
requested moment count.

## `numerics`

`numerics` selects alternative engine paths. Every value other than the default
needs the `native` [build preset](build-presets.md).

| Key | Type | Allowed values | Default | Meaning |
|---|---|---|---|---|
| `backend.name` | string | `"native"`, `"ginkgo"` | `"native"` | Sparse matrix–vector back end. |
| `backend.executor` | string | `"omp"`, `"cuda"` | `"omp"` | Ginkgo executor. Checked only when `name` is `"ginkgo"`. |
| `backend.format` | string | `"csr"`, `"sellp"`, `"ell"`, `"hybrid"`, `"coo"` | `"csr"` | Ginkgo matrix format. Checked only when `name` is `"ginkgo"`. |
| `operator` | string | `"csr"`, `"atom_block"` | `"csr"` | Operator storage: plain CSR, or atom blocks ordered in spatial domains. |
| `domain_atoms` | integer | positive | computed (64–16384, from the largest atom width and `state_block`) | Atoms per spatial domain for `atom_block`. |
| `execution` | string | `"openmp"`, `"cuda"`, `"mpi"`, `"mpi+cuda"` | `"openmp"` | Execution mode for `atom_block`. `"mpi"` and `"mpi+cuda"` pass validation but fail at execution with `is not wired yet`. |

Constraints the validator enforces:

| Condition | Error |
|---|---|
| `backend.name = "ginkgo"` in a non-native build | `numerics.backend ginkgo is only available in the native build profile` |
| `backend.name = "ginkgo"` without Ginkgo support | `numerics.backend ginkgo requested but this build has LSQUANT_WITH_GINKGO=OFF` |
| `backend.executor = "cuda"` without CUDA | `executor cuda requested but LSQUANT_WITH_CUDA=OFF` |
| `ginkgo` with `model.kind` other than `"csr"` or `basis.name` other than `"real"` | `numerics.backend ginkgo applies to real-space csr models` |
| `ginkgo` with `cuda` and a non-orthogonal metric | `metric composites not supported on cuda` |
| `operator = "atom_block"` in a non-native build | `numerics.operator='atom_block' is only available in the native build profile` |
| `atom_block` with `model.kind` other than `"csr"` | `numerics.operator='atom_block' requires model.kind='csr'` |
| `atom_block` with `backend.name = "ginkgo"` | `numerics.operator='atom_block' cannot be combined with numerics.backend ginkgo` |
| `atom_block` with `basis.name` other than `"real"` | `numerics.operator='atom_block' requires basis.name='real'` |
| `atom_block` with a non-orthogonal metric | `numerics.operator='atom_block' requires metric.policy='orthogonal'` |
| `atom_block` without `model.positions` | `numerics.operator='atom_block' requires model.positions` |
| `atom_block` without `orbitals_per_atom` or `atom_map` | `numerics.operator='atom_block' requires model.orbitals_per_atom or model.atom_map` |
| `hr_file` without a valid `model.supercell` | `model.hr_file requires model.supercell = [L1, L2] of positive ints` |
| `hr_file` without an explicit `expansion.window` | `model.hr_file requires an explicit expansion.window` |
| `execution` other than `"openmp"` in a non-native build, or without `atom_block` | `numerics.execution='<mode>' ... requires numerics.operator='atom_block'` (or the build-profile message) |
| `execution` with `cuda` or `mpi` not compiled in | `... this build has LSQUANT_WITH_CUDA=OFF` / `LSQUANT_WITH_MPI=OFF` |

With `hr_file`, only the operators `VX`, `VY` and `VZ` can be built (as
i[H, X_a]).

## `basis`

| Key | Type | Allowed values | Default | Meaning |
|---|---|---|---|---|
| `name` | string | `"real"`, `"bloch"`, `"eigen"` | `"real"` | Representation of H. See [The representation](../how-to/run-lsquant.md#the-representation-supercell-or-k-mesh). |
| `supercell` | list of 3 integers | positive | `[1, 1, 1]` | Supercell size (`real`) or k-point mesh (`bloch`, `eigen`) for `wannier` and `siesta`. Not used by `csr` or `generator`. |
| `disorder.W` | number | ≥ 0 | `0.0` | Width of the uniform Anderson disorder, in eV. Used by `wannier` and `siesta`; nonzero W is refused for `csr`; refused for `generator`. |
| `disorder.seed` | integer | — | `7` | Seed of the disorder draw. |
| `disorder.kind` | string | `"scalar"`, `"magnetic"`, `"independent"` | `"scalar"` | How the draw acts on spin. `csr`, `wannier` and `generator` models take `"scalar"` only, one value per basis function. See below. |

For `siesta`, the disorder potential is $V_{ij} = (w_i + w_j) S_{ij}/2$, with
draws uniform in `[-W/2, W/2]`. The spinor basis is interleaved: index $2i$ is
orbital $i$ with spin up and $2i+1$ the same orbital with spin down. The three
kinds set $w$ as follows:

- `scalar`: one value $u_i$ per orbital, on both spin components. The potential
  is spin-independent and keeps time-reversal symmetry.
- `magnetic`: one value $h_i$ per orbital, $+h_i$ on spin up and $-h_i$ on spin
  down. This is a random Zeeman field $h_i \sigma_z$, which breaks time-reversal
  symmetry.
- `independent`: one value per spinor coefficient. This is a scalar part plus a
  magnetic part, each of variance $W^2/24$. It is the layout of earlier
  versions, kept to reproduce their results.

For the same seed, `scalar` and `magnetic` draw the same values and differ only
in the sign on spin down. `basis.json` and `run.resolved.json` record the
kind and the distribution.

For `csr`, `wannier` and `generator` models, `"magnetic"` and `"independent"`
stop the run with `basis.disorder.kind='<k>' needs the spinor layout of a SIESTA
model; csr, wannier and generator models take kind 'scalar' (one value per
basis function)`.
<!-- runfile.py:183 ALLOWED_DISORDER_KEYS; :289 defaults; :501-510 validation and refusals; :965 resolution; disorder_bases.py:29 spinor_disorder_values; apps/stage6c_disorder.py:137-156 provenance -->

For `csr`, the CSR matrices are used as given. A `bloch` or `eigen` basis
there only changes the default window to `"eigen"`, which `csr` refuses.

## `metric`

| Key | Type | Allowed values | Default | Meaning |
|---|---|---|---|---|
| `policy` | string | `"orthogonal"`, `"nonorthogonal"` | `"orthogonal"` | Whether the basis has an overlap matrix S ≠ 1. |
| `route` | string | `"lowdin"`, `"generalised"`, `"anticommutator"` | `"lowdin"` | Treatment of S: S^(-1/2) H S^(-1/2), S^(-1) H, or an approximate anticommutator form (SIESTA only). |
| `applier` | string | `"exact"`, `"chebyshev"`, `"ischo"` | `"exact"` | How S^(-1) or S^(-1/2) is applied: dense factorisation, Chebyshev expansion, or preconditioned conjugate gradients (S^(-1) only). |
| `order` | integer | — | `60` | Chebyshev order of the `chebyshev` applier. |
| `coeff_tol` | number | — | `1e-12` | Coefficient truncation threshold of the `chebyshev` applier. |
| `tol` | number | — | `1e-13` | Conjugate-gradient tolerance of the `ischo` applier. Also the tolerance of the metric probe of real-space non-orthogonal runs with disorder, for every applier (`disorder_bases.py:315–319`); raise it for an ill-conditioned S (example 12 uses `1e-10`). |
| `max_iter` | integer | — | `200` | Maximum conjugate-gradient iterations of the `ischo` applier. |
| `max_dim` | integer | — | `4096` | Largest matrix dimension the `exact` applier accepts. |

| Condition | Error |
|---|---|
| `route = "anticommutator"` and `model.kind` is not `"siesta"` | `anticommutator metric is an explicitly approximate SIESTA route` |
| `policy = "nonorthogonal"`, `route = "lowdin"`, `applier = "ischo"` | `metric.applier='ischo' provides S^-1 only and cannot be used with route='lowdin'` |

## `projection`

Spectral projection works with the models of the metric executors: `siesta`,
`openmx`, `pyscf`, `paoflow`, the machine-learning kinds (`deeph`, `deeptb`, `hamgnn`, `threebodytb`, `abacus`), `gpaw`, and `wannier` with
`params.velocity_mode`, `spin` or `fermi_energy`. For any other model the validator stops with
`Expected a metric-executor model (...) when projection is enabled`.

| Key | Type | Allowed values | Default | Meaning |
|---|---|---|---|---|
| `kind` | string | `"none"`, `"hard"`, `"shift"`, `"orthogonal"` | `"none"` | Projection type. `orthogonal` reduces the dimension. |
| `window` | list of 2 numbers | finite, increasing | `null` | Energy window `[E_min, E_max]` in eV, inclusive, in stored energies (E_F = 0 for SIESTA). Required unless `kind` is `"none"`. |
| `radius` | number | finite, ≥ 0 | `null` | Real-space cutoff of the projected hoppings, in Å. `null` keeps every hopping. |

With `kind = "none"`, leave out `window` and `radius`; otherwise the validator
stops with `projection none does not accept window or radius`.
`run.resolved.json` contains `projection` only when `kind` is `"hard"`,
`"shift"` or `"orthogonal"`.

## `word` and `words`

| Key | Type | Required | Meaning |
|---|---|---|---|
| `word` | string | yes | An observable word, for example `Tr[delta(E-H)]`. Must be a non-empty string. |
| `words` | list of strings | no | Words evaluated on the same model in one run. Must be non-empty. When present, LSQuant runs these words and not `word`; `word` is still required and parsed. |

Each word must match one of the shapes in the [observable-word
catalogue](observable-words.md); run `python -m lsquant words list` to see
them. Operators in a word must exist in `model.operators` (or the discovered
CSR files) for `csr`. `run.resolved.json` stores the canonical form of every
word.

Custom operator names listed in `model.operators` are accepted by validation and
execution. For CSR models, provide each corresponding `<label>.<OP>.CSR` file;
a missing listed file raises `FileNotFoundError` naming that file.

## `expansion`

| Key | Type | Allowed values | Default | Meaning |
|---|---|---|---|---|
| `moments` | integer | positive | **required** | Number of Chebyshev moments M per expanded slot. |
| `window` | string or list | `"gershgorin"`, `"eigen"`, `[lo, hi]` with `hi > lo` | `"gershgorin"` for `basis.name = "real"`; `"eigen"` for `"bloch"` and `"eigen"` | Spectral window onto which H is rescaled: a Gershgorin bound, the exact eigenvalue range, or explicit bounds in eV. |
| `trace` | string | `"exact"`, `"random"` | `"exact"` | Trace over basis states, or a random-phase estimate. |
| `states` | integer | positive | `csr`: the dimension of H (`exact`) or 1 (`random`); `siesta`: 16 (`random`) | Number of trace states. For `csr` with `exact`, the first `states` basis vectors; must not exceed the dimension. |
| `seed` | integer | — | `12345` | Seed of the random-phase states. |
| `batch` | integer | positive | `moments` | Number of Chebyshev vectors held per batch. |
| `state_label` | string | — | `"exact"` (`exact`) or the generator label (`random`, `"default"`) | `csr` only: the label written into moment file names after `_state`. |
| `state_block` | integer or string | positive integer, `"auto"` | `1` | Number of trace states propagated together. `"auto"` picks a value from available memory in the native build and resolves to 1 in the golden build. The clean and disordered SIESTA executors treat `"auto"` as 1; the projection executor requires an integer and rejects `"auto"` during validation. |
| `left_window` | integer | positive | `32` | Bra-side streaming window of the native rank-2 path. |
| `right_window` | integer or string | positive integer, `"all"` | `16` | Ket-side storage window of the native rank-2 path; `"all"` stores all M vectors (resolves to `null`). |
| `checkpoint` | boolean | — | `false` | `csr` only: checkpoint after every trace state to `checkpoint_v2_<hash>.lsqckpt` in the output directory and resume from it on restart only when the calculation fingerprint matches. The file is deleted when the word finishes. |
| `state_block_rule` | string | — | — | Output only. Ignored on input. The CSR and ordinary Wannier executors record the effective numeric block size and its selection reason in `run.resolved.json`. |
| `trace_coverage` | object | — | — | Output only. Ignored on input. Written to `run.resolved.json` when an exact trace uses fewer states than the dimension: `{"states": n, "dimension": N}`, an average over the first n of N basis functions. |

| Condition | Error |
|---|---|
| `window = "eigen"` with `model.kind` `"csr"` or `"generator"` | `expansion.window='eigen' is not supported for model.kind='<kind>' (no k blocks; use gershgorin or an explicit window)` |
| `basis.name` `"bloch"` or `"eigen"` with `csr` and no `window` | Same message, raised by the executor (the default resolves to `"eigen"`). |
| `window = "gershgorin"` with `model.kind = "siesta"` | Raised by the SIESTA executors. |
| `window = "gershgorin"` with `model.hr_file` | `model.hr_file requires an explicit expansion.window` |

With an exact trace, a truncated local seed set still uses the legacy file
label `exact` and normalises by `states`. It emits a `RuntimeWarning`:
the result is an average over the first n of N basis functions, not the
trace. `expansion.trace_coverage = {"states": n, "dimension": N}` in
`run.resolved.json` records this case. Full block-local Bloch/eigen traces
use all orbitals at every k point and are not truncated local traces.
The coverage field is output only and ignored on input. Generator suites
with several dimensions also retain each distinct coverage in
`receipt.trace_coverages`.

Checkpoint filenames use a versioned, full SHA-256 calculation fingerprint.
It covers the engine word and expansion sizes, resolved numerical window,
expansion settings (including trace kind, seed, resolved state count and state
block), basis and disorder settings, metric, numerics, resolved atom-domain
size, build profile and time settings. Model settings and the size and SHA-256
of every model file read are included: Hamiltonian or H(R), loaded operators,
overlap matrix when used, and file-based atom-block positions and atom maps.
Inline geometry is covered by the model settings. Output and reconstruction
settings do not affect the fingerprint.

Files are streamed in 1 MiB chunks once per checkpoint-enabled execution, even
when multiple words run. This adds one sequential read of the input bytes,
with bounded memory; no mtime cache is used, so same-size edits with preserved
timestamps are detected. Keep model inputs unchanged during loading and execution.
Other fingerprints and old word-only checkpoints are reported with an
`Ignoring checkpoint` warning, left in place, and never passed to the loader.
An identical calculation can still resume byte-identically. This changes neither
the binary checkpoint format nor the C++ loader. Wannier and SIESTA executors
do not currently checkpoint.

## `time`

Required when a word contains `U(t)`. If the section is missing, compilation
fails with
`Word '<word>' contains time evolution U(t) without 'time' parameters (steps and tmax required)`.
With a static word, the validator accepts a `time` section and the run ignores
it.

| Key | Type | Allowed values | Default | Meaning |
|---|---|---|---|---|
| `steps` | integer | > 1 | **required** | Number of time points N_t, including t = 0. |
| `tmax` | number | — | `null` | Last time, in fs. The time step is `tmax / (steps - 1)`. |
| `dt` | number | — | `null` | Time step, in fs. Takes precedence over `tmax` when both are set. |
| `bessel` | string | `"scipy"`, `"legacy"` | `"scipy"` | Source of the Bessel coefficients of U(dt): SciPy, or the GSL values of the legacy code (golden byte compatibility). |

One of `tmax` or `dt` is required: `time requires 'tmax' or 'dt'`.
Time is in fs because LSQuant uses ħ = 0.6582119624 eV fs with energies in eV.
The SIESTA projection executor refuses a `time` section.

## `reconstruction`

| Key | Type | Allowed values | Default | Meaning |
|---|---|---|---|---|
| `broadening_meV` | number | positive | `null` | Target energy resolution of the Jackson kernel, in meV. With `null`, `csr` and `wannier` runs write moments but no curves; SIESTA runs fail. |
| `formulas` | list of strings | see [Formulas by word shape](#formulas-by-word-shape) | default formula of each word | Reconstruction formulas. Matching is case-insensitive. |
| `fermi_energies` | list of numbers | — | `[]` | Fermi energies in eV. `time_correlation` and `time_correlation_curve` write one curve per entry and require a nonempty list when broadening requests a curve (including the default MSD and time-correlation formulas). `row_m0` and `row0` need no list; null broadening keeps moments only. Explicit time-correlation targets also require their own nonempty list. SIESTA runs use the first entry (or 0.0) for the oracle and `spin_chern.json`, and every entry (or 0.0) for `fermi_values.json`. |
| `kernel` | string | any; only `"jackson"` is implemented | `"jackson"` | Kernel name. The value is lower-cased and recorded; every reconstruction uses the Jackson kernel. |
| `targets` | list of objects | see [Reconstruction targets](#reconstruction-targets) | `[]` | Extra reconstructions of the same moments, `csr` only. |
| `per_site_dos` | boolean | — | `false` | `csr` only: also write `DOS_reconstructed.dat`, the DOS divided by the dimension of H, from all M moments. |
| `legacy_text_round_trip` | boolean | — | `false` | Golden/legacy byte compatibility: write the legacy moment file and reconstruct from the copy read back (15 significant digits). `csr` and SIESTA runs only; the clean SIESTA and projection executors write the legacy file only when `output.legacy_layout` is `true`. |
| `fft_grid` | boolean | — | `false` | Recorded only. It does not change the curves. To reconstruct on the Chebyshev-node grid, use the formulas `greenwood_fft` or `bastin_fft`. With `true`, `bastin_backend = "folded"` is refused. |
| `num_points` | integer | — | `null` | Number of points of the default SIESTA oracle energy grid (1001 when `null`). Not used by KPM curves. |
| `energies` | list of numbers | non-empty, finite, strictly increasing | absent | SIESTA oracle energy grid, in eV. |
| `energy_grid` | list of 3 numbers | `[E_min, E_max, n]`, `E_min < E_max`, integer `n ≥ 2` | absent | SIESTA oracle energy grid as a uniform grid. Give `energies` or `energy_grid`, not both. |
| `num_cells` | integer | — | absent | Number of unit cells for the `shc` normalisation. SIESTA runs overwrite it with the k-mesh size. |
| `cell_area` | number | — | absent | Unit-cell area for the `shc` normalisation. SIESTA runs overwrite it with \|a1 × a2\| in Å². |
| `bastin_backend` | string | `"tables"`, `"folded"` | `"tables"` | Kubo–Bastin contraction used by `shc` in SIESTA runs: legacy tables, or a folded series with O(M² + MN) work. `"folded"` requires `model.kind = "siesta"` and `fft_grid = false`. |

When `energies` and `energy_grid` are both absent, the SIESTA oracle grid
spans the gap around the first Fermi energy, extended by
`max(0.25 × gap, 5 × broadening)`.
If the Fermi energy lies outside the spectrum, the run fails with
`Fermi level is outside the spectrum; supply reconstruction.energies or energy_grid`.

For ordinary `csr` and `wannier` runs, `formulas: ["shc"]` requires both
`reconstruction.num_cells` and `reconstruction.cell_area`. Missing either is
refused during validation, before computing moments. Reconstruction targets
also use these two keys. The metric executors obtain geometry from the model.
With a non-null broadening, the run writes the `SpinHall_…` curve directly.

To reconstruct it again, keep the moments next to `run.resolved.json`, which
supplies the geometry, and pass the same broadening explicitly:

```bash
python -m lsquant reconstruct <dir>/<file>.chebmom2D --formula shc --broadening-meV <same as the run>
```

A single file otherwise defaults to 10 meV and Fermi energies `[0.0]`.
For a directory and an explicit formula, incompatible ranks or word shapes
are skipped and reported in `skipped_files` (file and reason) and by the CLI.
If no file admits the formula, reconstruction fails before writing curves.
Legacy moment files take precedence over unified `.chebmom` files. Legacy
files have no word, and static engine words do not distinguish Greenwood
from Bastin, so filtering uses the union of catalogue formulas for their
rank and time axis; full observable words retain their shape restrictions.

The `shc` header uses the stored word: `JXSZ` on the left and `VY` on the
right give `sigma_z_xy [e^2/h, sigma_z Pauli]`; other spin components and
directions change the indices. A word without a left spin current is
labelled `normalised_bastin[A,B] [2*pi/(num_cells*cell_area)]`, naming its
operators without claiming spin Hall conductivity. Missing words are labelled
with unknown operators. The numbers and normalisation do not change. The
SIESTA executor continues to write its existing `shc.dat` header.

### Formulas by word shape

From `SHAPE_ALLOWED_FORMULAS` in `words/catalogue.py`. Aliases are in
parentheses.

| Word shape | Example word | Allowed formulas | Default |
|---|---|---|---|
| `dos_identity` | `Tr[delta(E-H)]` | `dos` (`dos_curve`) | `dos` |
| `dos_operator` | `Tr[A delta(E-H)]` | `dos` (`dos_curve`) | `dos` |
| `kubo_greenwood` | `Tr[VX delta(E-H) VX delta(E-H)]` | `greenwood` (`kubo_greenwood`), `greenwood_fft` (`kubo_greenwood_fft`) | `greenwood` |
| `kubo_bastin` | `Tr[VX delta(E-H) VY dG-(E)/dE]` | `bastin` (`kubo_bastin`), `sea` (`kubo_bastin_sea`), `surf` (`kubo_bastin_surf`), `bastin_fft` (`kubo_bastin_fft`), `shc` | `bastin`; `sea` for a `sea[...]` word |
| `msd` | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | `time_correlation` (`time_correlation_curve`) | `time_correlation` |
| `projected_evolution` | `Tr[P U(t)^+ SX U(t) delta(E-H) P]` | `row_m0` (`row0`), `time_correlation` (`time_correlation_curve`) | `row_m0` |
| `time_correlation` | `Tr[A delta(E-H) U(t)^+ B U(t)]` | `time_correlation` (`time_correlation_curve`) | `time_correlation` |

MSD words accept centred and off-centre enclosing spectral windows. Supply the
physical velocity `V = i[H, X]` in the operator file; the executor scales it by
`1/a`, with `a = (upper - lower)/2`, without adding an identity shift. The
Hamiltonian retains the requested window. Direct Python callers should use
`RescaledOperator(V, velocity_window(window))`; see [time evolution](../time-evolution.md).

With several words, LSQuant applies each formula to every word whose shape
allows it. Each formula must match at least one word; otherwise the validator
stops with
`Reconstruction formula '<f>' in reconstruction.formulas does not match any word in [...]`.
`row_m0` writes no curve file.

### Reconstruction targets

Each entry of `reconstruction.targets` reconstructs the moments of every word
again. Only the `csr` executor processes targets. The `wannier` executor
ignores them and the SIESTA executors refuse a non-empty list.

| Key | Type | Required | Default | Meaning |
|---|---|---|---|---|
| `formulas` | list of strings | yes | — | Formulas to apply. They are not checked against the word shapes by the validator; a formula that does not fit a word's moments fails at execution. |
| `broadening_meV` | number | yes | — | Broadening of this target, in meV. |
| `fermi_energies` | list of numbers | no | `[]` | Fermi energies for `time_correlation`. |
| `subdir` | string | no | none | Subdirectory of the output directory for the curves. Without it, curves go to the output directory. |
| `also_top_dir` | boolean | no | `false` | Also copy each curve to the output directory, if no file of that name exists there. |
| `legacy_text_round_trip` | boolean | no | `reconstruction.legacy_text_round_trip` | Reconstruct from the legacy moment file read back. |

Target curves are always written as `.dat` files, whatever
`output.legacy_names` says.

Finite-M Jackson reconstruction estimates the effective broadening as
`max(broadening_meV, 1000*pi*(upper-lower)/(2*M))`, in meV.
A request below this floor emits a `RuntimeWarning` with requested and
effective values; moment truncation and curve numbers are unchanged.
`receipt.curve_resolution` records the requested value and a list of
`effective_broadening_meV` values, one per spectral axis (time is excluded).
When `curves.json` is written, each curve also carries these fields.
For rectangular tensors the two spectral axes may have different floors.
The estimate describes KPM curves, not the separate Gaussian eigen oracle.

For the real-space nonorthogonal metric probe, the acceptance rule and default
`metric.tol` are unchanged. Exact-applier failures report measured residuals,
the rounding floor of an ill-conditioned overlap, an estimate of cond(S) from
existing overlap bounds, and advice to raise `metric.tol`. Chebyshev failures
continue to suggest increasing order or reducing `coeff_tol`.

## `output`

| Key | Type | Default | Meaning |
|---|---|---|---|
| `dir` | path | `"."` | Output directory, relative to the current working directory. `--out` on the command line overrides it. |
| `legacy_names` | boolean | `true` | `csr` and `wannier`: write each curve as a `.dat` file with the legacy name. With `false`, all curves go into one `curves.json`. Moment file names do not change. |
| `legacy_layout` | boolean | `true` | Write the legacy moment files `.chebmom1D`, `.chebmom2D` or `.chebmomTD`. With `false`, `csr` and `wannier` write the unified `.chebmom` instead. |
| `write_chebmom` | boolean | `true` if `model.kind = "wannier"`, `metric.policy = "nonorthogonal"`, or `legacy_layout = false`; otherwise `false` | `csr` and `wannier`: also write the unified `.chebmom` file. SIESTA runs always write it. |
| `write_basis_json` | boolean | `true` if `model.kind = "wannier"`, otherwise `false` | `wannier` only: write `<stem>.basis.json` with the basis provenance. SIESTA runs always write `basis.json`. |
| `capture_tutorial` | boolean | `false` | Golden/legacy byte compatibility: copy curve files (and `projected_evolution` moment files) to `captured/NNNN_<name>` for the tutorial golden check. |

See [Output files](output-files.md) for the names and formats.

## Provenance keys

These keys let a `run.resolved.json` be read back as a run file.

| Key | Type | Default | Meaning |
|---|---|---|---|
| `build_profile` | string | the loaded profile | If present, it must equal the profile of the loaded `_core` module (`golden` or `native`), otherwise: `Run file build_profile=... differs from loaded engine profile=...`. Every run writes it. |
| `runtime` | object | — | Ignored on input. Every run writes the values of `OMP_NUM_THREADS`, `OMP_PROC_BIND`, `OMP_PLACES`, `OMP_DYNAMIC`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `LSQUANT_SPMV_MIN_NNZ`, `LSQUANT_KERNEL_BATCH` and `LSQUANT_FFT_BACKEND`. |
| `receipt` | object | — | Ignored on input. SIESTA executors write phase timings and diagnostics. |
| `reference` | boolean | automatic | Models of the metric executors only (`Expected a metric-executor model (...) when 'reference' is set`). `true` forces and `false` suppresses the exact eigenbasis reference of a clean SIESTA run. By default it runs unless the production run already computes the exact eigenbasis trace. The disordered and projection executors ignore it. |

## SIESTA executors

`execute_runfile` sends a `siesta` run to one of three executors. An `openmx`
or `paoflow` run, or a `wannier` run with `params.velocity_mode`, `spin` or
`fermi_energy`, goes to the projection or the clean executor by the same rule. Each checks
the resolved run file before it reads the SIESTA files and fails with the
quoted message.

| Executor | Selected when |
|---|---|
| Projection (`projection/runfile.py`) | `projection.kind` is not `"none"`. |
| Disordered (`apps/stage6c_disorder.py`) | otherwise, if `basis.name = "real"`, `basis.disorder.W ≠ 0`, or `model.params.nonorthogonal_disorder = true`. |
| Clean (`apps/stage6b_siesta.py`) | otherwise. |

| Check | Clean | Disordered | Projection |
|---|---|---|---|
| `basis.name` | `bloch`, `eigen` | `real`, `bloch` | `real`, `bloch`, `eigen` |
| `metric.policy` | `nonorthogonal` | `nonorthogonal` | `nonorthogonal` |
| `metric.route` | any (`anticommutator` allowed) | `lowdin`, `generalised` | `lowdin`, `generalised` |
| `metric.applier` | `exact` | `exact` (`bloch`); `exact`, `chebyshev` (`real`) | `exact`; also `chebyshev` for `real` |
| `basis.disorder.W` | must be 0 | any | must be 0 |
| Words | `Tr[delta(E-H)]` or `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` | same | `Tr[delta(E-H)]`, `Tr[A delta(E-H)]`, or `Tr[JXSZ delta(E-H) VY dG-(E)/dE]`; no `time` section |
| `reconstruction.formulas` | subset of `shc`, `dos`, `dos_curve` | same | same |
| `reconstruction.targets` | must be empty | must be empty | must be empty |
| `reconstruction.broadening_meV` | required | required, positive, finite | required, positive, finite |
| `expansion.window` | `eigen` or `[lo, hi]` enclosing the generalised spectrum | `eigen` or `[lo, hi]` enclosing the clean spectrum plus √cond(S) × max\|w\| | `eigen` or `[lo, hi]` enclosing the projected spectrum |
| `model.params.oracle` | allowed | refused | not used |
| Other | — | — | `projection.kind = "orthogonal"` refused for `basis.name = "real"` |

The clean executor reconstructs every rank-2 word with `shc` and every rank-1
word with `dos`, whatever `formulas` lists.

## Keys for golden and legacy byte compatibility

These keys exist so that LSQuant reproduces the reference outputs of the
legacy code byte for byte. New calculations do not need them.

| Key | Effect |
|---|---|
| `output.legacy_names` | Legacy curve file names (default `true`). |
| `output.legacy_layout` | Legacy moment file layout (default `true`). |
| `output.capture_tutorial` | Copies outputs into `captured/` for the tutorial golden check. |
| `reconstruction.legacy_text_round_trip` | Reconstructs from the 15-digit legacy moment file. |
| `time.bessel = "legacy"` | GSL Bessel coefficients. |

## Examples

### Minimal run file

This is the smallest run file that writes a curve. It reads
`<dir>/haldane.HAM.CSR` and computes the density of states with the default
Gershgorin window and an exact trace. The block also sets
`reconstruction.broadening_meV`, which every curve needs (see
[`reconstruction`](#reconstruction)):

```json
{
  "model": {"kind": "csr", "dir": "goldens/stage0/haldane/operators", "label": "haldane"},
  "word": "Tr[delta(E-H)]",
  "expansion": {"moments": 256},
  "reconstruction": {"broadening_meV": 10.0}
}
```

The first tutorial file, `examples/01_dos_1d_chain/run_dos.json`, is a
generator run file: it builds the chain itself instead of reading model files.

### Annotated SIESTA run file

`examples/07_spin_hall_bi_bilayer/run_shc_quick.json`:

```json
{
  "model": {
    "kind": "siesta",
    "dir": "examples/07_spin_hall_bi_bilayer/siesta",
    "label": "bi",
    "params": {"velocity_mode": "covariant", "oracle": true}
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "metric": {"policy": "nonorthogonal", "route": "lowdin", "applier": "exact"},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 256, "window": "eigen", "trace": "exact", "seed": 12345},
  "reconstruction": {
    "kernel": "jackson",
    "broadening_meV": 20.0,
    "formulas": ["shc"],
    "fermi_energies": [0.0]
  },
  "output": {"dir": "out/bi_quick", "legacy_layout": true, "write_chebmom": true}
}
```

| Key | Value | Effect |
|---|---|---|
| `model.dir` | `examples/07_spin_hall_bi_bilayer/siesta` | Resolved against the repository root. The directory holds `bi.fdf`, `bi.TSHS` and `Bi.ion.xml`. |
| `model.params.velocity_mode` | `"covariant"` | Velocities include the orbital dipoles, so `Bi.ion.xml` is read and the dipole quadrature runs (or its cache is reused). |
| `model.params.oracle` | `true` | Also writes `exact_shc.json` and `spin_chern.json`. |
| `basis` | `eigen`, `[24, 24, 1]` | Exact eigenstates on a 24 × 24 × 1 k-mesh; `num_cells` = 576. No disorder, so the clean executor runs. |
| `metric` | `nonorthogonal`, `lowdin`, `exact` | The only applier the clean executor accepts. |
| `word` | spin Hall word | Shape `kubo_bastin`; rank-2 moments. |
| `expansion.moments` | `256` | M = 256 in each direction. |
| `expansion.window` | `"eigen"` | Window from the generalised eigenvalues; `run.resolved.json` records the numeric bounds. |
| `expansion.trace` | `"exact"` | Exact block trace in the eigenbasis. Because production already computes this trace, the eigen reference is skipped (`reference` unset). |
| `expansion.seed` | `12345` | Recorded; random phases are not used with an exact trace. |
| `reconstruction.broadening_meV` | `20.0` | Jackson resolution of `shc.dat`, and the Gaussian width (0.020 eV) of the oracle. |
| `reconstruction.formulas` | `["shc"]` | Spin Hall conductivity in e²/h (Pauli σ_z). |
| `reconstruction.fermi_energies` | `[0.0]` | Fermi energy of the default oracle grid and of `spin_chern.json`. |
| `output.dir` | `out/bi_quick` | Relative to the current working directory. |
| `output.legacy_layout` | `true` | Also writes the `.chebmom2D` file next to the `.chebmom` file. |
| `output.write_chebmom` | `true` | Has no effect in SIESTA runs, which always write `.chebmom`. |

Defaults that this file leaves to the resolver: `expansion.batch` = 256,
`expansion.state_block` = 1, `reconstruction.legacy_text_round_trip` = `false`,
`output.legacy_names` = `true`.

## See also

- [Command-line interface](cli.md)
- [Output files](output-files.md)
- [Observable words](observable-words.md)
- [Non-orthogonal bases](../explanation/nonorthogonal-bases.md)
- [Build presets](build-presets.md)
