# Use a DeepH-format model

Start from a Hamiltonian predicted by DeepH-pack, DeepH-E3, xDeepH or MACE-H
and compute its density of states, or, for a model with spin–orbit coupling,
its spin Hall conductivity. These codes predict the DFT Hamiltonian H(R) in
the atomic-orbital basis of the DFT code and use the DFT overlap S(R). LSQuant
reads the prediction as the code wrote it. Check the training set, the
transferability and the agreement with DFT in your code. You should know your
code and have read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

The DeepH reader has been tested only on synthetic files in both layouts, not
yet on a real prediction; compare your first results with DFT.
<!-- source: formats read from the sources of DeepH-pack 141faaa, DeepH-E3 036bf95, MACE-H 3d6248b, xDeepH 99e6440 and DeepH-dock 60ad934; unit tests use synthetic files in both layouts, no real prediction yet -->

LSQuant computes these quantities for the clean crystal on a k-mesh; a
real-space supercell or disorder stops the run. Every key is listed in the
[run-file reference](../reference/run-file.md#machine-learning-hamiltonians-in-the-metric-executors).
<!-- source: python/lsquant/runfile.py:1576–1579 -->

## Requirements

You need:

- LSQuant built and importable; see [Install LSQuant](install.md).
- `h5py`, which reads the HDF5 files of the prediction.

## 1. Run the inference

LSQuant reads the directory of one structure after inference. It accepts two
layouts and detects which one the directory uses.

### DeepH-E3, MACE-H, xDeepH and DeepH-pack before JAX

Run the inference so that the prediction is written to the processed-structure
directory:

- DeepH-E3 and MACE-H: in `eval.ini`, set `inference = True` and
  `output_dir` equal to the processed data directory, then run
  `deephe3-eval.py eval.ini`. Each structure directory then contains
  `hamiltonians_pred.h5`.
- xDeepH: prepare the overlap with `deeph-inference` (task `[1]`), then run
  `xdeeph-eval.py` with your inference config, as in the xDeepH README.
- DeepH-pack: run `deeph-inference` with the default tasks; the result is
  written to the work directory.

The directory must contain `hamiltonians_pred.h5`, `overlaps.h5`,
`orbital_types.dat`, `site_positions.dat`, `lat.dat` (or `rlat.dat`) and
`info.json`. These files come from the preprocessing of the DFT code
(OpenMX, ABACUS, SIESTA or FHI-aims); keep them next to the prediction.

### DeepH-pack (JAX) with DeepH-dock

The directory of a structure contains `POSCAR`, `info.json`, `overlap.h5` and
the prediction `hamiltonian_pred.h5`. For a new structure, DeepH-dock computes
`overlap.h5` from `POSCAR` and the basis files (`dock compute overlap -h`).

### Spin–orbit coupling

For spin–orbit coupling, train and predict a spinful model, one whose
Hamiltonian acts on both spin components of every orbital (`isspinful` or
`spinful` is `true` in `info.json`).

## 2. Transfer the model

There is nothing to run. LSQuant reads the directory:

| File | Needed for |
|---|---|
| `hamiltonians_pred.h5` or `hamiltonian_pred.h5` | H(R), one block per atom pair |
| `overlaps.h5` or `overlap.h5` | S(R) |
| `hkb.h5` | the spin–orbit part, when `info.json` sets `split_hkb`: the Hamiltonian file then contains the spinless part H₀, and `hkb.h5` the spinful Kleinman–Bylander part H_KB, with H = H₀ on both spins + H_KB |
| `lat.dat`, `site_positions.dat`, `orbital_types.dat` or `POSCAR` | cell, atomic positions, orbitals per atom |
| `info.json` | spin, Fermi energy |
<!-- source: python/lsquant/models/deeph.py:31–33 (split_hkb: H = I_2 ⊗ H_0 + H_KB), 175 (split_hkb needs spinful) -->

Without a prediction file, LSQuant reads the DFT file (`hamiltonians.h5` or
`hamiltonian.h5`), so you can run the same calculation on the DFT Hamiltonian
and compare. `model.params.hamiltonian_file` selects another file in the
directory.

## 3. Run LSQuant

This run file computes the spin Hall conductivity of a spinful model on a
24 × 24 k-mesh, with the spectrum shifted inside [−3, 3] eV:

```json
{
  "model": {
    "kind": "deeph",
    "dir": "/path/to/structure",
    "params": {"velocity_mode": "berry_connection", "fermi_energy": -3.85}
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "projection": {"kind": "shift", "window": [-3, 3]},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 6.0,
                     "energy_grid": [-0.5, 1.0, 751]},
  "output": {"dir": "out/deeph_shc"}
}
```

For the density of states, use `"word": "Tr[delta(E-H)]"` and
`"formulas": ["dos"]`; this also works for a spinless model.

The template has the settings of the projected SIESTA run of
[example 20](../../examples/20_pdse2_spin_hall_routes/README.md), with a
narrower energy grid. It is the exact reference path of LSQuant, and each
setting does the following:

- `word` is the observable, here the spin Hall conductivity σ^z_xy: the spin
  current with spin along z and flow along x, in response to an electric
  field along y. See [Observable words](../reference/observable-words.md).
- `basis.name = "eigen"` diagonalises H(k) at every point of the k-mesh, and
  `basis.supercell` is that mesh, Γ-centred, with 1 as the third entry for a
  2D system. With `expansion.trace = "exact"`, LSQuant sums over all states,
  so the curve has no random-vector noise.
- Leave out the `metric` section: LSQuant treats the model as non-orthogonal
  and orthogonalises the basis with S^(−1/2) (Löwdin), computed exactly.
- `projection` with `"kind": "shift"` moves every state outside [−3, 3] eV,
  relative to E_F, to the nearest edge of the window. The expansion then
  covers at most 7.2 eV instead of the full spectrum, and needs fewer moments.
- `expansion.moments` sets the energy resolution, as the smearing does for a
  DFT density of states. With `"window": "eigen"`, LSQuant pads the projected
  spectrum by 10 %, to at most [−3.6, 3.6] eV. With 1024 moments, the finest
  broadening is πa/M ≈ 11 meV, with a = 3.6 eV the half-width of the window.
  LSQuant applies this value in place of the requested 6 meV; to resolve
  6 meV, use about 1900 moments. See
  [Number of moments and broadening](use-siesta.md#number-of-moments-and-broadening).
- `reconstruction.energy_grid`, `[E_min, E_max, n]` in eV relative to E_F,
  is the grid of the exact reference: here −0.5 to 1.0 eV in 751 points,
  2 meV apart. The KPM curve has its own energy grid.

Whether the 24 × 24 mesh is converged depends on your material: increase it
until the curve stops changing.
<!-- source: template close to examples/20_pdse2_spin_hall_routes/run_siesta_projected.json (energy_grid differs); metric: reference/run-file.md (Machine-learning Hamiltonians: exact Löwdin route, metric can be left out), runfile.py:1580–1583; shift: projection/spectral.py:119–125; window pad: bases.py:1410–1416; resolution: reconstruct.py:35–36; energy_grid: apps/stage6b_siesta.py:132–143, projection/runfile.py:214 -->

Check the file, then run it:

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Values and meaning |
|---|---|
| `model.params.velocity_mode` | `"berry_connection"` (the default: derivative of H(k) and S(k) with the orbital positions) or `"bare"`. There is no position matrix, so `"covariant"` stops the run. |
| `model.params.fermi_energy` | In eV. LSQuant shifts the Hamiltonian to H − E_F S. Default: the value in `info.json` (`fermi_level` or `fermi_energy_eV`), or 0 if there is none; with 0, energies stay on the absolute scale of the DFT code. |
| `model.params.hamiltonian_file` | File name of H(R) in `model.dir`. Default: the prediction, then the DFT file. |
| `model.params.oracle` | `true` adds an exact reference in a run without a projection. A projection run writes its exact reference in any case. |
| `projection.window` | In eV, measured from `model.params.fermi_energy`. |
<!-- source: reference/run-file.md (Machine-learning Hamiltonians table); absolute scale: models/deeph.py docstring ("H is stored on the absolute scale of the DFT code"); oracle: apps/stage6b_siesta.py:322, projection/runfile.py:214–228 -->

With a spinless model, `JXSZ` and the other spin observables stop the run
with `... needs spin, but the model is spinless`.
<!-- source: reference/run-file.md (Machine-learning Hamiltonians, last paragraph) -->

With a projection, LSQuant writes the KPM curve to `word_0.dat` and the exact
result to `word_0_oracle.json`, which contains the exact curve both before and
after the projection; `projection.json` gives the difference between the KPM
curve and the exact one. Without a projection, it writes `dos.dat` or
`shc.dat`. The spin Hall curve is σ^z_xy in units of e²/h, with the Pauli
matrix σ_z in the spin current; divide by 2 to compare with values for spin
ħ/2. The DOS is per unit cell, in states/eV/cell. See
[Output files](../reference/output-files.md#projection-runs) and
[Results and units](use-siesta.md#results-and-units).
<!-- source: reference/output-files.md (Projection runs; Clean SIESTA runs); docs/explanation/kubo-formulas.md -->

## 4. Recommendations

- Check the prediction before transport: compare its bands with DFT on a
  path, with the band tools of your code.
- Set `model.params.fermi_energy` from a calculation with the predicted
  Hamiltonian (for example the Fermi level of a DOS on a dense mesh). The
  value in `info.json` comes from the DFT run of the training structure or is
  absent.
- LSQuant replaces H(R) by its Hermitian part, the mean of H(R) and H(−R)†,
  and records the relative size of the removed part in `basis.json`
  (`hamiltonian_hermiticity_residual`). It warns above 10⁻³ of the largest
  element; a large value points to a prediction problem.
- The overlap must belong to the same structure and basis as the prediction.
  LSQuant stops the run when a block does not match the orbitals of its
  atoms.
