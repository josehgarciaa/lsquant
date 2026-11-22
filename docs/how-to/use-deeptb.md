# Use a DeePTB model

Connect a DeePTB tight-binding model (DeePTB-SK or DeePTB-E3) to LSQuant and
compute its density of states. DeePTB predicts H(R), and S(R) for a
non-orthogonal model, from a trained network and a structure. LSQuant reads
the matrices that DeePTB exports; the quality of the model (training data,
transferability, agreement with DFT) is yours. You should know DeePTB and have
read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

Every key is listed in the
[run-file reference](../reference/run-file.md#machine-learning-hamiltonians-in-the-metric-executors).

## Requirements

- LSQuant built and importable; see [Install LSQuant](install.md).
- `h5py`, which reads the exported HDF5 files.
<!-- format read from the sources of DeePTB 597a036 (postprocess/unified/system.py, pardiso backend); unit tests use synthetic export files, no real model yet -->

## 1. Turn on in DeePTB

Nothing to change in the training. DeePTB exports the real-space matrices of
a structure for its sparse band-structure code; LSQuant reads that export.

## 2. Transfer the model

In the Python environment of DeePTB:

```python
from dptb.postprocess.unified.system import TBSystem

tbsys = TBSystem(data="structure.vasp", calculator="model.pth")
print(tbsys.get_efermi(kmesh=[24, 24, 1]))   # Fermi energy for the run file
tbsys.to_pardiso_json(output_dir="deeptb_export")
```

The command `dptb pdso band.json -i model.pth -stu structure.vasp -o
deeptb_export` writes the same files before it runs DeePTB's Julia band code.
The directory holds:

| File | Used for |
|---|---|
| `predicted_hamiltonians.h5` | H(R), one block per atom pair |
| `predicted_overlaps.h5` | S(R); absent for an orthogonal model |
| `structure.json` | cell, atomic positions, orbitals per element |

The older `to_pardiso` export (`cell.dat`, `positions.dat`,
`atomic_numbers.dat`, `basis.dat`) is also read.

## 3. Run LSQuant

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

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Value |
|---|---|
| `velocity_mode` | `"berry_connection"` or `"bare"`. There is no position matrix, so `"covariant"` is refused. |
| `fermi_energy` | In eV, removed as H − E_F S (H − E_F for an orthogonal model). The export holds no Fermi energy: take it from `get_efermi`. Default 0. |

## 4. Recommendations

- DeePTB-SK adds its on-site spin–orbit coupling when it builds H(k); the
  export holds only the spinless blocks. LSQuant refuses an export marked
  `spinful` in `structure.json`, so spin Hall runs are not possible with
  DeePTB models for now.
- Check the bands of the export with DeePTB's own band code before transport.
- DeePTB's Python band code uses the phase e^{−ik·R}, its Julia code
  e^{+ik·R}. For the real spinless blocks of the export the two give the same
  bands; LSQuant follows the Julia code.
