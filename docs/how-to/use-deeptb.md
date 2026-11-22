# Use a DeePTB model

Start from a DeePTB tight-binding model (DeePTB-SK or DeePTB-E3) and compute
its density of states. DeePTB predicts H(R), and S(R) for a non-orthogonal
model, from a trained network and a structure. LSQuant reads the matrices
that DeePTB exports. Check the training data, the transferability and the
agreement with DFT in DeePTB. You should know DeePTB and have read
[Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

> **Known limitation: no spin Hall conductivity with DeePTB.** DeePTB-SK adds
> its on-site spin–orbit coupling when it builds H(k), and the export contains
> only the spinless blocks. LSQuant stops the run on an export marked
> `spinful` in `structure.json`, so DeePTB models are currently restricted to
> spinless runs and the density of states.
<!-- source: python/lsquant/models/deeptb.py:32–34, 96–98 -->

The DeePTB reader has been tested only on synthetic export files, not yet on
a real model; compare your first results with DFT.
<!-- format read from the sources of DeePTB 597a036 (postprocess/unified/system.py, pardiso backend); unit tests use synthetic export files, no real model yet -->

LSQuant computes the density of states of the clean crystal on a k-mesh; a
real-space supercell or disorder stops the run. Every key is listed in the
[run-file reference](../reference/run-file.md#machine-learning-hamiltonians-in-the-metric-executors).
<!-- source: python/lsquant/runfile.py:1576–1579 -->

## Requirements

You need:

- LSQuant built and importable; see [Install LSQuant](install.md).
- `h5py`, which reads the exported HDF5 files.

## 1. Prepare the model in DeePTB

Use your trained model and structure as they are; the export needs no change
in the training. DeePTB exports the real-space matrices of a structure for its
sparse band-structure code, its Julia/Pardiso band code,
and LSQuant reads that export.
<!-- source: python/lsquant/models/deeptb.py:3–5 ("for its Julia/Pardiso band code") -->

## 2. Transfer the model

In the Python environment of DeePTB, compute the Fermi energy and export the
matrices:

```python
from dptb.postprocess.unified.system import TBSystem

tbsys = TBSystem(data="structure.vasp", calculator="model.pth")
print(tbsys.get_efermi(kmesh=[24, 24, 1]))   # Fermi energy for the run file
tbsys.to_pardiso_json(output_dir="deeptb_export")
```

The command `dptb pdso band.json -i model.pth -stu structure.vasp -o
deeptb_export` writes the same files before it runs DeePTB's Julia band code.
The directory contains:

| File | Needed for |
|---|---|
| `predicted_hamiltonians.h5` | H(R), one block per atom pair |
| `predicted_overlaps.h5` | S(R); absent for an orthogonal model |
| `structure.json` | cell, atomic positions, orbitals per element |

The older `to_pardiso` export (`cell.dat`, `positions.dat`,
`atomic_numbers.dat`, `basis.dat`) is also read.

## 3. Run LSQuant

This run file computes the density of states on a 48 × 48 k-mesh:

```json
{
  "model": {
    "kind": "deeptb",
    "dir": "/path/to/deeptb_export",
    "params": {"velocity_mode": "berry_connection", "fermi_energy": -4.21}
  },
  "basis": {"name": "eigen", "supercell": [48, 48, 1]},
  "word": "Tr[delta(E-H)]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["dos"], "broadening_meV": 10.0,
                     "energy_grid": [-3.0, 3.0, 601]},
  "output": {"dir": "out/deeptb_dos"}
}
```

The template is the exact reference path of LSQuant, and each setting does
the following:

- `word` is the observable, here the density of states `Tr[delta(E-H)]`. See
  [Observable words](../reference/observable-words.md).
- `basis.name = "eigen"` diagonalises H(k) at every point of the k-mesh, and
  `basis.supercell` is that mesh, Γ-centred, with 1 as the third entry for a
  2D system. With `expansion.trace = "exact"`, LSQuant sums over all states,
  so the curve has no random-vector noise.
- Leave out the `metric` section: LSQuant orthogonalises a non-orthogonal
  model with S^(−1/2) (Löwdin), computed exactly.
- `expansion.moments` sets the energy resolution, as the smearing does for a
  DFT density of states. With `"window": "eigen"`, LSQuant takes the window
  from the eigenvalues and pads it by 10 % of the bandwidth on each side. The
  finest broadening is πa/M, with a the half-width of the window: with
  M = 1024, the requested 10 meV is reached only when a is below about
  3.3 eV. For a wider spectrum, LSQuant applies πa/M instead; see
  [Number of moments and broadening](use-siesta.md#number-of-moments-and-broadening).
- `reconstruction.energy_grid` sets only the grid of the exact reference,
  which this run does not compute, so it has no effect here. The KPM curve
  has its own energy grid.
- `model.params.velocity_mode` has no effect on the density of states.

Whether the 48 × 48 mesh is converged depends on your material: increase it
until the curve stops changing.
<!-- source: eigen + exact: use-siesta.md (Basis table); metric: reference/run-file.md (Machine-learning Hamiltonians), runfile.py:1580–1583; window pad: bases.py:1410–1416; resolution: reconstruct.py:35–36, a = M η / π = 1024 × 0.010 / π = 3.26 eV; energy_grid only in oracle_energies (apps/stage6b_siesta.py:132–143, 322–329), no oracle and no projection in this template; the DOS word uses no velocity operator -->

Check the file, then run it:

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Values and meaning |
|---|---|
| `model.params.velocity_mode` | `"berry_connection"` (the default) or `"bare"`. There is no position matrix, so `"covariant"` stops the run. |
| `model.params.fermi_energy` | In eV. LSQuant shifts the Hamiltonian to H − E_F S, or H − E_F for an orthogonal model. The export contains no Fermi energy, so take it from `get_efermi`; the default is 0, which leaves energies on DeePTB's scale. |
<!-- source: reference/run-file.md (Machine-learning Hamiltonians table); models/deeptb.py:29–30 ("Energies are on DeePTB's scale; the files hold no Fermi energy") -->

LSQuant writes the density of states to `dos.dat`: energy in eV and the DOS
per unit cell, in states/eV/cell. See
[Output files](../reference/output-files.md#clean-siesta-runs).
<!-- source: reference/output-files.md (Clean SIESTA runs, dos.dat); apps/stage6b_siesta.py writes dos.dat for metric-executor models -->

## 4. Recommendations

- Check the bands of the export with DeePTB's own band code before
  transport. DeePTB's Python band code uses the phase e^(−ik·R) and its Julia
  code e^(+ik·R); for the real spinless blocks of the export the two give the
  same bands, so either serves as a check. LSQuant follows the Julia code.
