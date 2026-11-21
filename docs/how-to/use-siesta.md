# Use a SIESTA calculation

This guide shows you how to prepare a SIESTA calculation for LSQuant, write the
run file, and stay within the settings that LSQuant accepts for SIESTA models.
It assumes you know SIESTA and have read [How we run LSQuant](run-lsquant.md).

LSQuant computes the density of states and the intrinsic spin Hall
conductivity of a SIESTA crystal with spin–orbit coupling, clean or with
Anderson disorder. Two use cases show it end to end:

- [Example 07: spin Hall conductivity of a Bi(111) bilayer](../../examples/07_spin_hall_bi_bilayer/README.md);
- [Example 08: spin Hall conductivity of a PdSe₂ monolayer](../../examples/08_spin_hall_pdse2_monolayer/README.md).

Every key used here is listed in the [run-file reference](../reference/run-file.md#siesta-models).

## Before you begin

You need:

- LSQuant built and importable; see [Install LSQuant](install.md).
- The Python package `sisl`, which reads the SIESTA files. The examples were
  run with `sisl` 0.16.4.

  Install it with pip:

  ```bash
  pip install sisl
  ```

## Check what LSQuant accepts

LSQuant reads only non-orthogonal SIESTA Hamiltonians with spin–orbit
coupling. For any other SIESTA output, the run stops with this message:

```text
load_siesta_model requires a nonorthogonal SIESTA SOC Hamiltonian
```
<!-- source: python/lsquant/models/siesta.py:311-312 -->

SIESTA's atomic-orbital basis is always non-orthogonal, so in practice the
condition is `Spin spin-orbit`. LSQuant cannot load a collinear or
spin-unpolarised calculation with `model.kind = "siesta"`; the run stops with
this message.

## Prepare the SIESTA run

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

LSQuant treats the energies stored in the `.TSHS` file as already measured
from the Fermi level, E_F = 0, and uses them unshifted. The Fermi level that
SIESTA prints goes into `basis.json` as information only. Give every energy
in the run file (windows, Fermi energies, projection windows) relative to
E_F.
<!-- source: python/lsquant/models/siesta.py:4-5, 340-345 -->

## Write the run file

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

The template asks for 20 meV, but the run delivers a broadening of about
0.6 eV. The finest broadening LSQuant can deliver is πa/M, with a the
half-width of the window and M the number of moments. The spectrum of the Bi
bilayer spans about [−24.7, 56.4] eV. With the default 10 % padding, the
window becomes [−32.8, 64.5] eV, a half-width of a ≈ 48.6 eV. With 256
moments, πa/M ≈ 0.6 eV, and LSQuant uses this value instead of the requested
20 meV, without a warning. The DOS of example 07's quick run, which has the
same settings, shows the cap: its 7681 points are the grid of 30 points per
moment for all 256 moments. To resolve 20 meV, you need about 7600 moments.
Alternatively, narrow the window with a projection run:
[example 16](../../examples/16_projection_downfold/README.md) keeps the states
in [−2, 2] eV, so the expansion covers 4 eV instead of 81 eV.

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

`inspect` checks the keys and the word. LSQuant applies further checks when
the run starts, most of them before it reads the SIESTA files. A file that
`inspect` accepts can therefore still stop at once in `run`. The next
sections describe those checks.

## Know which kind of SIESTA run you have

LSQuant handles a SIESTA run file in one of three ways. The kind of run
follows from the contents of the run file:

| Your run file has | Kind of run | Typical use |
|---|---|---|
| `projection.kind` other than `"none"` | projection | Keep only the states in an energy window. |
| `basis.name = "real"`, or `basis.disorder.W` ≠ 0 | disordered | Anderson disorder, real-space supercells. |
| anything else | clean | Clean crystal on a k-mesh. |

Each kind of run accepts a different set of settings. The full table is in
[SIESTA executors](../reference/run-file.md#siesta-executors). The rules you
meet most often are below.

### Rules for every SIESTA run

- **Words.** Only the DOS, `Tr[delta(E-H)]`, and the spin Hall word
  `Tr[JXSZ delta(E-H) VY dG-(E)/dE]`. Projection runs also accept
  `Tr[A delta(E-H)]` with one operator, for example `WINDOW`.
  Time-dependent words are refused.
- **Formulas.** A subset of `shc`, `dos` and `dos_curve`; `dos_curve` is
  another name for `dos`. For the units of `shc`, see
  [Check the result](#check-the-result).
- **Window.** `expansion.window` is `"eigen"` or an explicit `[lo, hi]` that
  encloses the spectrum. `"gershgorin"` is refused. Because the default
  window of `basis.name = "real"` is `"gershgorin"`, set `"window": "eigen"`
  explicitly for a real-space run.
- **Targets.** Leave `reconstruction.targets` empty; a non-empty list stops
  the run. Reconstruct again afterwards with
  [`python -m lsquant reconstruct`](../reference/cli.md#reconstruct).

### Clean runs

- `basis.name` is `"bloch"` or `"eigen"`.
- `metric.applier` is `"exact"`. `metric.route` can be `"lowdin"`,
  `"generalised"`, or the approximate `"anticommutator"`.
- `basis.disorder.W` is 0. A non-zero value makes the run a disordered one.
- With `model.params.oracle = true`, the run also writes the exact per-k spin
  Hall conductivity (`exact_shc.json`) and the spin Chern number
  (`spin_chern.json`).
- A clean run reconstructs every spin Hall word with `shc` and every DOS word
  with `dos`, whatever `formulas` lists.

### Disordered runs

- `basis.name` is `"real"` (a supercell of `basis.supercell` unit cells) or
  `"bloch"`.
- `metric.route` is `"lowdin"` or `"generalised"`. `metric.applier` is
  `"exact"`, or also `"chebyshev"` for `"real"`.
- The exact applier is refused when the supercell dimension exceeds
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
  - `"independent"`: one energy per spin-orbital, the layout of earlier
    versions. Use it only to reproduce earlier results.

  The [run-file reference](../reference/run-file.md#basis) gives the exact
  draws.
<!-- source: max_dim python/lsquant/disorder_bases.py:264-265, runfile.py:299; kind: disorder_bases.py:29, runfile.py:501-510 -->

Use case: [Example 12: Anderson disorder in a non-orthogonal SIESTA basis](../../examples/12_disorder_nonorthogonal/README.md).

### Projection runs

- `projection.kind` is `"hard"`, `"shift"` or `"orthogonal"`, and
  `projection.window` is the energy window in eV relative to E_F.
  `"orthogonal"` is refused with `basis.name = "real"`.
- `basis.disorder.W` is 0, and the run file has no `time` section.
- A word with the operator `WINDOW` also needs `"operators": ["WINDOW"]` in
  `model`. When `model.operators` is set, validation accepts only the listed
  names, so list every operator of every word.

Use case: [Example 16: spectral projection of a Bi Hamiltonian onto [−2, 2] eV](../../examples/16_projection_downfold/README.md).

## Choose the basis

| Basis | What it does | Use it for |
|---|---|---|
| `eigen` | Diagonalises S⁻¹H at every k and works in the eigenbasis. | Clean crystals. With `trace: "exact"` it is the cheapest route and the reference for the others. |
| `bloch` | Keeps one H(k), S(k) block per k-point. | Clean crystals; the disordered Bloch layout. |
| `real` | Builds a real-space supercell of `basis.supercell` unit cells. | Disorder. |

For a clean crystal, `eigen` and `bloch` give the same moments. Start with
`eigen` and a coarse k-mesh, as example 07 does with 24 × 24, then increase
the mesh.

## Choose the velocity mode

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
<!-- source: python/lsquant/apps/stage6b_siesta.py:203, python/lsquant/apps/stage6c_disorder.py:111 -->

[Bases and metrics](../reference/bases-and-metrics.md) and
[Non-orthogonal bases](../explanation/nonorthogonal-bases.md) explain the
modes.

With `covariant`, the first run integrates the dipoles and caches them. Later
runs on the same SIESTA files reuse the cache, whatever the k-mesh, metric or
word; a change in the input files invalidates it. LSQuant writes the cache to
the first of these locations that applies:

1. the directory in the environment variable `LSQUANT_CACHE_DIR`;
2. `.lsquant-cache` in the parent of `model.dir`, if `model.dir` is named
   `siesta`;
3. `lsquant-cache` in the system temporary directory, if `model.dir` is not
   writable;
4. `.lsquant-cache` in `model.dir` itself.

The order keeps the cache out of the SIESTA input directory where it can. To
keep the cache in a fixed place, set `LSQUANT_CACHE_DIR`.
<!-- source: python/lsquant/models/siesta.py:287-289 (docstring: cache outside siesta/ input directories, content-based invalidation), 382-397 -->

The quadrature is controlled by `quadrature_order`, `quadrature_max_order` and
`quadrature_tolerance` in `model.params`.

## Set the energy window and the number of moments

With `"window": "eigen"`, LSQuant takes the window from the eigenvalues of the
generalised problem. SIESTA spectra are wide, so a SIESTA run needs more
moments than a tight-binding model for the same broadening. A broadening finer
than πa/M silently uses all M moments, with a the half-width of the window and
M the number of moments. For the numbers of the Bi bilayer, see
[Write the run file](#write-the-run-file).

The spin Hall curve is an integral over energy that starts at the bottom of
the reconstruction grid, at E = b − 0.95 a, with b the centre of the window.
All occupied bands must lie inside this inner 95 % of the window. The 10 %
padding of `"window": "eigen"` keeps them there; check it when you give an
explicit window.

In a disordered run, the window is also padded by √cond(S) × max|w|, where w is
the disorder potential. For the Bi bilayer at W = 1 eV, the window is
[−68.5, 100.2] eV, and 256 moments resolve only about 1 eV. A projection run
narrows the spectrum, as shown in [Write the run file](#write-the-run-file).

## Relax `metric.tol` for an ill-conditioned overlap

A real-space disordered run checks the inverse overlap S⁻¹ and S^(−1/2) with a
probe, and stops if the residual exceeds `metric.tol`:

```text
metric probes exceed 1e-13; increase Chebyshev order or reduce coeff_tol: {...}
```

The default `metric.tol` is 10⁻¹³. In double precision, with
ε ≈ 2.2 × 10⁻¹⁶, the rounding errors of an exact factorisation of S can reach
cond(S) × ε. For the DZP basis of the Bi bilayer, cond(S) ≈ 3550, so rounding
errors can reach about 8 × 10⁻¹³. The measured residual, 1.4 × 10⁻¹³, is of
rounding size, yet the default refuses it: 10⁻¹³ asks for more than double
precision guarantees for this overlap.

Raise the tolerance in the `metric` section:

```json
"metric": {"policy": "nonorthogonal", "route": "lowdin", "applier": "exact", "tol": 1e-10}
```

To see the residuals of a run, look at `metric_residuals` in its
`basis.json`. Example 12 shows the case and checks the result against the
Bloch layout.

> **Known limitation: `metric.tol` and its error message.** The default
> tolerance does not scale with cond(S)
> (`python/lsquant/disorder_bases.py:315–319`, `python/lsquant/runfile.py:297`).
> The message suggests changing the Chebyshev order or `coeff_tol` even with
> the exact applier, where those keys have no effect. The
> [run-file reference](../reference/run-file.md#metric) describes `tol` only
> as the tolerance of the `ischo` applier. Set `tol` as shown above.

## Check the result

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

The spin Hall conductivity is a rank-2 word. For the same mesh and number of
random vectors, its cost grows as M², against M for the DOS, so expect it to
take much longer than the DOS run with the same M.

See [Run on a cluster](run-on-a-cluster.md) for threads and memory.
