# Use a DeepH-format model

Connect a Hamiltonian predicted by DeepH-pack, DeepH-E3, xDeepH or MACE-H to
LSQuant and compute its density of states or, with spin–orbit coupling, its
spin Hall conductivity. These codes predict the DFT Hamiltonian H(R) in the
atomic-orbital basis of the DFT code and use the DFT overlap S(R). LSQuant
reads the prediction as the code wrote it; the quality of the prediction
(training set, transferability, agreement with DFT) is yours. You should know
your code and have read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

Every key is listed in the
[run-file reference](../reference/run-file.md#machine-learning-hamiltonians-in-the-metric-executors).

## Requirements

- LSQuant built and importable; see [Install LSQuant](install.md).
- `h5py`, which reads the HDF5 files of the prediction.
<!-- formats read from the sources of DeepH-pack 141faaa, DeepH-E3 036bf95, MACE-H 3d6248b, xDeepH 99e6440 and DeepH-dock 60ad934; unit tests use synthetic files in both layouts, no real prediction yet -->

## 1. Turn on in your code

LSQuant reads the directory of one structure after inference. Two layouts
exist, and LSQuant detects which one it has.

**DeepH-E3, MACE-H, xDeepH and DeepH-pack before its JAX version.** Run the
inference so that the prediction lands in the processed-structure directory:

- DeepH-E3 and MACE-H: in `eval.ini`, set `inference = True` and
  `output_dir` equal to the processed data directory, then run
  `deephe3-eval.py eval.ini`. Each structure folder then holds
  `hamiltonians_pred.h5`.
- xDeepH: prepare the overlap with `deeph-inference` (task `[1]`), then run
  `xdeeph-eval.py` with your inference config, as in the xDeepH README.
- DeepH-pack: `deeph-inference` with the default tasks; the work directory
  holds the result.

The directory must hold `hamiltonians_pred.h5`, `overlaps.h5`,
`orbital_types.dat`, `site_positions.dat`, `lat.dat` (or `rlat.dat`) and
`info.json`. These files come from the preprocessing of the DFT code
(OpenMX, ABACUS, SIESTA or FHI-aims); keep them next to the prediction.

**DeepH-pack (JAX version) with DeepH-dock.** The directory of a structure
holds `POSCAR`, `info.json`, `overlap.h5` and the prediction
`hamiltonian_pred.h5`. For a new structure, DeepH-dock computes `overlap.h5`
from `POSCAR` and the basis files (`dock compute overlap -h`).

For spin–orbit coupling, train and predict a spinful model
(`isspinful` or `spinful` is `true` in `info.json`).

## 2. Transfer the model

Nothing to run. LSQuant reads the directory:

| File | Used for |
|---|---|
| `hamiltonians_pred.h5` or `hamiltonian_pred.h5` | H(R), one block per atom pair |
| `overlaps.h5` or `overlap.h5` | S(R) |
| `hkb.h5` | the spin–orbit part, when `info.json` sets `split_hkb` |
| `lat.dat`, `site_positions.dat`, `orbital_types.dat` or `POSCAR` | cell, atomic positions, orbitals per atom |
| `info.json` | spin, Fermi energy |

Without a prediction file, LSQuant reads the DFT file (`hamiltonians.h5` or
`hamiltonian.h5`), which helps to compare the prediction with DFT.
`params.hamiltonian_file` selects another file in the directory.

## 3. Run LSQuant

Spin Hall conductivity of a spinful model:

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

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Value |
|---|---|
| `velocity_mode` | `"berry_connection"` (derivative of H(k) and S(k) with the orbital positions) or `"bare"`. There is no position matrix, so `"covariant"` is refused. |
| `fermi_energy` | In eV, removed as H − E_F S. Default: the value in `info.json` (`fermi_level` or `fermi_energy_eV`), or 0 if there is none. |
| `hamiltonian_file` | File name of H(R) in `dir`. Default: the prediction, then the DFT file. |
| `projection.window` | Measured from `fermi_energy`. |

The metric section can be left out: the model is non-orthogonal and LSQuant
uses the exact Löwdin route. A spinless model refuses `JXSZ` and the other
spin words.

## 4. Recommendations

- Check the prediction before transport: compare its bands with DFT on a
  path, with the band tools of your code. LSQuant only checks that the files
  are complete and consistent.
- Set `fermi_energy` from a calculation with the predicted Hamiltonian
  (for example the Fermi level of a DOS on a dense mesh). The value in
  `info.json` comes from the DFT run of the training structure or is absent.
- LSQuant makes H(R) Hermitian as the mean of H(R) and H(−R)†, records the
  difference in `basis.json` (`hamiltonian_hermiticity_residual`) and warns
  above 10⁻³ of the largest element. A large value points to a prediction
  problem.
- The overlap must belong to the same structure and basis as the prediction.
  LSQuant stops when a block does not match the orbitals of its atoms.
