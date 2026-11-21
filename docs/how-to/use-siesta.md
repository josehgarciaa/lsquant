# Use a SIESTA calculation

Start from your own SIESTA calculation and end with a converged density of
states and spin Hall conductivity. Save the Hamiltonian in SIESTA, write the
run file, then converge the number of moments, the k-mesh or supercell, and
the broadening. You should know SIESTA and have read
[Running LSQuant](run-lsquant.md).

For a SIESTA crystal with spin–orbit coupling, clean or with Anderson
disorder, LSQuant computes two quantities: the density of states and the
intrinsic spin Hall conductivity. The longitudinal and Hall conductivities
σ_xx and σ_xy, time-dependent observables, checkpoints, the atom-block layout
and the GPU code are not available for SIESTA models. Two use cases show the
full workflow:
<!-- source: quantities and time-dependent: reference/run-file.md (SIESTA executors, Words row; "no time section"); checkpoint: reference/run-file.md (expansion.checkpoint "csr only"); atom-block: reference/run-file.md (numerics: "atom_block requires model.kind='csr'"); GPU: how-to/run-on-a-cluster.md (CUDA: "SIESTA and Wannier models run on the CPU"), reference/run-file.md (numerics: ginkgo "applies to real-space csr models") -->

- [Example 07: spin Hall conductivity of a Bi(111) bilayer](../../examples/07_spin_hall_bi_bilayer/README.md);
- [Example 08: spin Hall conductivity of a PdSe₂ monolayer](../../examples/08_spin_hall_pdse2_monolayer/README.md).

The [table in the README](../../README.md#quantities-by-model-source) compares
what each model source computes.

Every key used here is listed in the [run-file reference](../reference/run-file.md#siesta-models).

## Requirements

You need:

- LSQuant built and importable; see [Install LSQuant](install.md).
- The Python package `sisl`, which reads the SIESTA files. The examples were
  run with `sisl` 0.16.4.

  Install it with pip:

  ```bash
  pip install sisl
  ```

## The SIESTA calculation

LSQuant reads only non-orthogonal SIESTA Hamiltonians with spin–orbit
coupling. The atomic-orbital basis of SIESTA is always non-orthogonal, so the
condition on your calculation is `Spin spin-orbit`. For any other SIESTA
output, collinear and spin-unpolarised runs included, a run with
`model.kind = "siesta"` stops with this message:

```text
load_siesta_model requires a nonorthogonal SIESTA SOC Hamiltonian
```
<!-- source: python/lsquant/models/siesta.py:311-312 (refuses orthogonal or not is_spinorbit) -->

1. Turn on spin–orbit coupling in the `.fdf` file:

   ```text
   Spin  spin-orbit
   ```

2. Ask SIESTA to save the Hamiltonian and the overlap:

   ```text
   TS.HS.Save  true     # writes <SystemLabel>.TSHS
   SaveHS      true     # writes <SystemLabel>.HSX
   ```

   LSQuant needs one of the two files. If both are present, it reads the
   `.TSHS` file. The input of example 07,
   `examples/07_spin_hall_bi_bilayer/siesta/bi.fdf`, sets both.

3. Converge the self-consistent cycle with a k-grid that is dense enough for
   your system. LSQuant uses H(R) and S(R) exactly as SIESTA wrote them and
   runs no SIESTA step itself.

4. Collect these files in one directory:

   | File | Needed for |
   |---|---|
   | `<label>.fdf` | Always. `sisl` reads the geometry, the species and the orbitals from it. |
   | `<label>.TSHS` or `<label>.HSX` | Always. H(R) and S(R). |
   | `<species>.ion.xml`, one per species | `velocity_mode = "covariant"` only. The radial functions of the basis orbitals. `<species>` is the label in `ChemicalSpeciesLabel`, for example `Bi.ion.xml`. |

   `<label>` is the `SystemLabel` of the SIESTA run. LSQuant reads only these
   files, not the pseudopotential files.

Give every energy in the run file (windows, Fermi energies, projection
windows) relative to the Fermi level, E_F = 0. LSQuant takes the energies
stored in the `.TSHS` or `.HSX` file as already measured from E_F and uses
them unshifted. The Fermi level printed by SIESTA is copied to `basis.json`
for information only.
<!-- source: python/lsquant/models/siesta.py:4-5 (docstring), 310 (one read path for .TSHS and .HSX), 343 (energy_reference "stored H, E_F already zero"), 345-347 (printed_fermi_level_eV recorded only) -->

## The run file

Start from a run file of example 07. This one computes the spin Hall
conductivity of the clean crystal on a 24 × 24 k-mesh:

```json
{
  "model": {
    "kind": "siesta",
    "dir": "/path/to/my_siesta_run",
    "label": "bi",
    "params": {"velocity_mode": "covariant", "oracle": true}
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "metric": {"policy": "nonorthogonal", "route": "lowdin", "applier": "exact"},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 256, "window": "eigen", "trace": "exact"},
  "reconstruction": {"broadening_meV": 20.0, "formulas": ["shc"], "fermi_energies": [0.0]},
  "output": {"dir": "out/shc_quick"}
}
```

This template is the quick run of example 07, a test of the workflow rather
than a converged result. Its broadening is set to 20 meV, but with 256
moments the curve is resolved only to about 0.6 eV; see
[Number of moments and broadening](#number-of-moments-and-broadening).

Adapt the template:

1. Set `model.dir` to the directory with the SIESTA files and `model.label`
   to the `SystemLabel`. Use an absolute path. LSQuant tries a relative path
   against the LSQuant source tree first and only then against the directory
   of the run file, so a relative path can pick up the wrong directory.
2. Set `basis.supercell` to the k-mesh. For a 2D system the third entry is 1.
   The mesh is Γ-centred: k = (i/N1, j/N2, l/N3) · B, with B the reciprocal
   lattice vectors.
3. Keep `metric.policy = "nonorthogonal"`. Every SIESTA run requires it.
4. Set `reconstruction.broadening_meV`. A SIESTA run fails when it is
   missing.
5. `model.params.oracle = true` adds an exact reference: LSQuant also
   diagonalises H(k) at every k-point and computes the spin Hall
   conductivity exactly, for comparison with the KPM curve.
<!-- source: path resolution python/lsquant/runfile.py:900-915; k mesh python/lsquant/models/bloch.py:12-18, models/siesta.py:121-123 -->

Check the file, then run it:

```bash
python -m lsquant inspect my_run.json
python -m lsquant run my_run.json --out out/shc_quick
```

The `--out` option overrides `output.dir` of the run file.
<!-- source: python/lsquant/cli.py:133 (out_dir_override) -->

`inspect` checks the keys and the observable. Further checks run when the
run starts, most of them before the SIESTA files are read, so a file that
`inspect` accepts can still stop at once in `run`. The next section lists
those checks.

## Clean, disordered and projected calculations

LSQuant distinguishes three kinds of SIESTA calculation, a clean crystal, a
disordered supercell and a projected Hamiltonian. It takes the kind from the
contents of the run file:

| Your run file has | Kind of run | Typical use |
|---|---|---|
| `projection.kind` other than `"none"` | projection | Keep only the states in an energy window. |
| `basis.name = "real"`, or `basis.disorder.W` ≠ 0 | disordered | Anderson disorder, real-space supercells. |
| anything else | clean | Clean crystal on a k-mesh. |

Each kind of run accepts a different set of settings. The full table is in
the [run-file reference](../reference/run-file.md#siesta-executors); the
rules you meet most often follow.

### Settings for every kind

- The `word` is the DOS, `Tr[delta(E-H)]`, or the spin Hall conductivity,
  `Tr[JXSZ delta(E-H) VY dG-(E)/dE]`. Projection runs also accept
  `Tr[A delta(E-H)]` with one operator, for example `WINDOW`. A
  time-dependent observable stops the run.
- `reconstruction.formulas` is a subset of `shc`, `dos` and `dos_curve`;
  `dos_curve` is another name for `dos`. For the units of `shc`, see
  [Results and units](#results-and-units).
- `expansion.window` is `"eigen"` or an explicit `[lo, hi]` that encloses the
  spectrum. `"gershgorin"` stops the run. Because the default window of
  `basis.name = "real"` is `"gershgorin"`, set `"window": "eigen"` explicitly
  for a real-space run.
- `reconstruction.targets` holds extra reconstructions of the same moments
  with other broadenings. Leave it empty, because a non-empty list stops the
  run, and reconstruct again afterwards with
  [`python -m lsquant reconstruct`](../reference/cli.md#reconstruct).
<!-- source: targets defined in reference/run-file.md (reconstruction table: "Extra reconstructions of the same moments"; Reconstruction targets: SIESTA executors refuse a non-empty list) -->

### Clean crystals

- `basis.name` is `"bloch"` or `"eigen"`.
- `metric.applier` is `"exact"`. `metric.route` can be `"lowdin"`,
  `"generalised"`, or the approximate `"anticommutator"`.
- `basis.disorder.W` is 0. A non-zero value makes the run a disordered one.
- With `model.params.oracle = true`, the run also writes the exact per-k spin
  Hall conductivity (`exact_shc.json`) and the spin Chern number
  (`spin_chern.json`).
- A clean run reconstructs every spin Hall observable with `shc` and every
  DOS observable with `dos`, whatever `formulas` lists.

### Disordered supercells

- `basis.name` is `"real"` (a supercell of `basis.supercell` unit cells) or
  `"bloch"`.
- `metric.route` is `"lowdin"` or `"generalised"`. `metric.applier` is
  `"exact"`, or also `"chebyshev"` for `"real"`.
- `"applier": "exact"` stops the run when the supercell dimension exceeds
  `metric.max_dim` (default 4096). For larger supercells, use
  `"applier": "chebyshev"`.
- `model.params.oracle` must be `false`, because an exact per-k reference
  has no meaning with disorder.
- The disorder potential is (D S + S D)/2, with draws uniform in
  [−W/2, W/2]. To average over disorder, run one run file per
  `basis.disorder.seed`.
- Choose how the disorder acts on spin with `basis.disorder.kind`:
  - `"scalar"` (default): one energy per orbital on both spin components, a
    spin-independent potential that keeps time-reversal symmetry;
  - `"magnetic"`: a random Zeeman field along z, +h on spin up and −h on spin
    down of each orbital. It breaks time-reversal symmetry and acts as
    magnetic disorder on the spin Hall conductivity;
  - `"independent"`: the layout of earlier versions, kept only to reproduce
    their results.

  The [run-file reference](../reference/run-file.md#basis) gives the exact
  draws of each kind.
<!-- source: max_dim python/lsquant/disorder_bases.py:264-265, runfile.py:299; kind: disorder_bases.py:29, runfile.py:501-510 -->

[Example 12](../../examples/12_disorder_nonorthogonal/README.md) shows this
case for Anderson disorder in the non-orthogonal SIESTA basis of the Bi
bilayer.

### Projected Hamiltonians

- `projection.kind` is `"hard"`, `"shift"` or `"orthogonal"`, and
  `projection.window` is the energy window in eV relative to E_F.
  `"orthogonal"` stops the run with `basis.name = "real"`.
- `basis.disorder.W` is 0, and the run file has no `time` section.
- An observable with the operator `WINDOW` also needs
  `"operators": ["WINDOW"]` in `model`. When `model.operators` is set,
  validation accepts only the listed names, so list every operator of every
  observable.

[Example 16](../../examples/16_projection_downfold/README.md) projects a Bi
Hamiltonian onto the states in [−2, 2] eV.

## Basis

| Basis | What it does | Use it for |
|---|---|---|
| `eigen` | Diagonalises S⁻¹H at every k and works in the eigenbasis. | Clean crystals. With `trace: "exact"` it is the cheapest route and the reference for the others. |
| `bloch` | Keeps one H(k), S(k) block per k-point. | Clean crystals; disorder without a supercell matrix, the potential acting on the H(k), S(k) blocks through a Fourier transform. |
| `real` | Builds a real-space supercell of `basis.supercell` unit cells. | Disorder. |
<!-- source: bloch with disorder: examples/12_disorder_nonorthogonal/README.md (V = ½{B D B†, S(k)}, no supercell matrix stored) -->

For a clean crystal, `eigen` and `bloch` give the same moments.

## Velocity operators

`model.params.velocity_mode` sets how LSQuant builds the velocity operators,
and with them the spin current `JXSZ`, from H(k), S(k) and the orbital
positions:

| Mode | Includes | Needs |
|---|---|---|
| `bare` | ∂H(k)/∂k only. | — |
| `berry_connection` | ∂H/∂k, the intracell position phases, and the overlap connection −½{H, S⁻¹ ∂S/∂k}. | — |
| `covariant` | Everything in `berry_connection`, plus the orbital dipole matrix elements, integrated from the radial functions. | `<species>.ion.xml` for every species. |

The default depends on the kind of SIESTA run: `berry_connection` for clean
and projection runs, `covariant` for disordered runs. Set the mode
explicitly, so that clean and disordered runs of the same system use the same
velocity. Examples 07 and 08 use `covariant`.
[Non-orthogonal bases](../explanation/nonorthogonal-bases.md) explains the
modes.
<!-- source: python/lsquant/apps/stage6b_siesta.py:203, python/lsquant/apps/stage6c_disorder.py:111 -->

With `covariant`, the first run integrates the dipoles and caches them. Later
runs on the same SIESTA files reuse the cache, whatever the k-mesh, metric or
observable; a change in the input files invalidates it. To keep the cache in
a fixed place, set the environment variable `LSQUANT_CACHE_DIR`. Without it,
the cache goes to `.lsquant-cache` in the parent of `model.dir` when that
directory is named `siesta`. Otherwise it goes to `lsquant-cache` in the
system temporary directory when `model.dir` is not writable, and to
`.lsquant-cache` in `model.dir` itself in every other case.
<!-- source: python/lsquant/models/siesta.py:287-289 (docstring: cache outside siesta/ input directories, content-based invalidation), 382-391 (order: LSQUANT_CACHE_DIR, parent of siesta/, temp dir if not writable, model.dir) -->

The accuracy of the dipole integrals is set by `quadrature_order`,
`quadrature_max_order` and `quadrature_tolerance` in `model.params`.
Examples 07 and 08 use the defaults, which the
[run-file reference](../reference/run-file.md#siesta-models) lists.
<!-- source: no quadrature_* key in examples/07_*/run_*.json or examples/08_*/run_*.json; defaults 24, 384, 1e-5 in reference/run-file.md (SIESTA models) -->

## Number of moments and broadening

Choose the number of moments M from the broadening you need. M sets the
energy resolution of every curve, as the smearing does for a DFT density of
states. The finest broadening of a run is πa/M, with a the half-width of the
energy window. A finer value in `broadening_meV` uses all M moments, and
LSQuant applies πa/M instead, without a warning.

> **Rule of thumb.** For a broadening η, use M ≈ πa/η moments, where a is
> the half-width of the energy window in eV.

SIESTA spectra are wide, so a SIESTA run needs more moments than a
tight-binding model for the same broadening. With `"window": "eigen"`, LSQuant
takes the window from the eigenvalues of the generalised problem and pads it
by 10 %. The spectrum of the Bi bilayer spans about [−24.7, 56.4] eV, so the
window becomes [−32.8, 64.5] eV, with a half-width a ≈ 48.6 eV. With 256
moments, πa/M ≈ 0.6 eV, and LSQuant uses this value in place of the requested
20 meV. The DOS file of example 07's quick run, which has the same settings,
confirms it: its 7681 energy points show that all 256 moments were used. To
resolve 20 meV, you need about 7600 moments.
<!-- source: 7681 points = grid of 30 points per moment for all 256 moments (previous page; generic grid num_div = 30*max_mom in python/lsquant/reconstruct.py:101 gives 7680, the extra point of the SIESTA writer not checked); window margin runfile.py:1603-1604 (0.10 * width) -->

To need fewer moments, narrow the window with a projection run instead:
[example 16](../../examples/16_projection_downfold/README.md) keeps the states
in [−2, 2] eV, so the expansion covers 4 eV instead of 81 eV.

In a disordered run, the window is also padded by √cond(S) × max|w|, where w is
the disorder potential. For the Bi bilayer at W = 1 eV, the window is
[−68.5, 100.2] eV, and 256 moments resolve only about 1 eV.

The spin Hall curve is an integral over energy that starts at the bottom of
the reconstruction grid, at E = b − 0.95 a, with b the centre of the window.
All occupied bands must lie inside this inner 95 % of the window. With
`"window": "eigen"`, the 10 % padding keeps them there; check this when you
give an explicit window.

## k-mesh and supercell

For a clean crystal, `basis.supercell` is the k-mesh and sets the
Brillouin-zone sampling error, as in SIESTA. Start with `eigen` and a coarse
k-mesh, as example 07 does with 24 × 24, then increase the mesh. The
production setting of example 07 uses a 48 × 48 mesh with M = 1024 and a
requested broadening of 10 meV. With a ≈ 48.6 eV, 1024 moments resolve about
πa/M ≈ 0.15 eV, so this setting is also limited by M. Its converged value at
E_F = 0 is still under study.
<!-- source: examples/07_spin_hall_bi_bilayer/README.md (Step 2; Step 3: "converged KPM value of the production setting at E_F = 0 is still pending"), run_*_production.json; π × 48.6 / 1024 = 0.149 eV -->

For a disordered run, `basis.supercell` is the number of unit cells of the
real-space supercell and sets the finite-size error, as the k-mesh does for a
crystal. Average the result over several `basis.disorder.seed` values.

## Overlap tolerance `metric.tol`

A real-space disordered run checks the inverse overlap S⁻¹ and S^(−1/2) with a
probe, and stops if the residual exceeds `metric.tol`:

```text
metric probes exceed 1e-13; increase Chebyshev order or reduce coeff_tol: {...}
```

The default `metric.tol` is 10⁻¹³. In double precision, with
ε ≈ 2.2 × 10⁻¹⁶, the rounding errors of an exact factorisation of S can reach
cond(S) × ε. For the DZP basis of the Bi bilayer, cond(S) ≈ 3550, so rounding
errors can reach about 8 × 10⁻¹³. The measured residual, 1.4 × 10⁻¹³, is of
rounding size but exceeds the default, a tolerance tighter than double
precision can guarantee for this overlap.

Raise the tolerance in the `metric` section:

```json
"metric": {"policy": "nonorthogonal", "route": "lowdin", "applier": "exact", "tol": 1e-10}
```

To see the residuals of a run, look at `metric_residuals` in its
`basis.json`. Example 12 shows the case and checks the result against the
`bloch` basis.

> **Known limitation: `metric.tol` and its error message.** The default
> tolerance stays at 10⁻¹³ whatever cond(S) is. The message suggests changing
> the Chebyshev order or `coeff_tol` even with `"applier": "exact"`, where
> those keys have no effect. Set `tol` as shown above.
<!-- source: python/lsquant/disorder_bases.py:315–319, python/lsquant/runfile.py:297 -->

## Results and units

After the run, look at these files in the output directory:

| File | What to check |
|---|---|
| `basis.json` | Number of orbitals, k-mesh, window, `metric_residuals`, velocity mode, the Fermi level printed by SIESTA. |
| `dos.dat` | Two columns: energy (eV) and DOS per unit cell (states/eV/cell). The generic curve files of LSQuant hold the extensive DOS Nρ(E) of the whole system; `dos.dat` is divided by the number of unit cells. Its integral up to E_F should equal the number of valence electrons per cell. Example 07 gives 30.02 for the 30 occupied bands of the Bi bilayer. |
| `shc.dat` | Two columns: energy (eV) and σ^z_xy(E) in units of e²/h. |
| `exact_shc.json`, `spin_chern.json` | Clean runs with `oracle: true`: the exact per-k result to compare with `shc.dat`. |
| `run.resolved.json` | The run file with every default filled in, and the phase timings under `receipt`. |
<!-- source: python/lsquant/apps/stage6b_siesta.py:286-296 (dos divided by num_cells; headers) -->

The spin current in `shc.dat` is defined with the Pauli matrix σ_z rather
than the spin ħσ_z/2: J = ½{v_x, σ_z}. To compare with values quoted for
spin ħ/2 in units of (ħ/e) e²/h, multiply `shc.dat` by 1/2. The header of the
file calls this convention `sigma_z Pauli`.
<!-- source: python/lsquant/reconstruct.py:731-743; factor 1/2 derived -->

The [output-file reference](../reference/output-files.md) describes every file.

## Cost

These times were measured on a shared 96-core node with the `golden` build:

| Run | Threads | Time | Memory |
|---|---|---|---|
| Example 07, DOS, 24 × 24 k-mesh, `eigen`, M = 256 | 8 | 12 min | 1.2 GB |
| Example 12, DOS, 4 × 4 × 1 real-space supercell (dimension 1152), W = 1 eV, 4 random vectors, M = 256 | 1 | 64 s | not recorded |
| Example 16, projection onto [−2, 2] eV, 8 × 8 Bloch mesh, M = 256 | 1 | 14 min | not recorded |

On one thread of a heavily loaded shared node, the quick DOS and spin Hall
runs of example 07 took 18 min each.
<!-- source: examples/07_spin_hall_bi_bilayer/README.md (Runtime) -->

The spin Hall conductivity needs a double Chebyshev expansion, with an M × M
table of moments. For the same k-mesh and number of random vectors, its cost
grows as M², against M for the DOS, so at large M expect a spin Hall run to
take much longer than a DOS run with the same M.
<!-- source: docs/explanation/kpm.md (Cost: R M² inner products; "for large M the inner products dominate") -->

See [Run on a cluster](run-on-a-cluster.md) for threads and memory.
