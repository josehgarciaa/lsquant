# Use a HamGNN model

Connect a Hamiltonian predicted by HamGNN to LSQuant and compute its density
of states or its spin Hall conductivity. HamGNN predicts H(R) in the
atomic-orbital basis of OpenMX, ABACUS or SIESTA/HONPAS, with spin–orbit
coupling when the model was trained with it; the overlap S(R) comes from the
DFT code. LSQuant reads the prediction as HamGNN wrote it; the quality of the
model (training data, transferability, agreement with DFT) is yours. You
should know HamGNN and have read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

Every key is listed in the
[run-file reference](../reference/run-file.md#machine-learning-hamiltonians-in-the-metric-executors).

## Requirements

- In the Python environment of HamGNN: NumPy (already there). The export step
  opens `graph_data.npz`, which needs HamGNN's `torch` and `torch_geometric`.
  LSQuant does not need to be installed there.
- LSQuant built and importable where you run it; see [Install LSQuant](install.md).
<!-- format read from the sources of HamGNN 64dba1d (DFT_interfaces/openmx/band_cal.py, graph_data_gen.py, openmx_postprocess/read_openmx.c, abacus/band_cal.py, siesta/read_siesta.py, models/Model.py); unit tests use synthetic files, no real model yet -->

## 1. Turn on in HamGNN

Nothing to change in the training. Build the graph of the structure you want
with HamGNN's usual tools (for OpenMX, `openmx_postprocess` and
`graph_data_gen.py`; the structure needs the overlap of a non-self-consistent
DFT step, as for HamGNN's own band structure). Then run a prediction with
`stage: test` on that `graph_data.npz`. HamGNN writes
`prediction_hamiltonian.npy` in the log directory of the run.

Note the three settings of your model: `nao_max` (14, 19 or 26 for OpenMX;
27 or 40 for ABACUS), and whether it uses `soc_switch` (spin–orbit coupling)
or `spin_colinear` (collinear spin).

## 2. Transfer the model

In the HamGNN environment, run the export script that comes with LSQuant:

```bash
python /path/to/lsquant/python/lsquant/models/hamgnn.py \
    graph_data.npz version_1/prediction_hamiltonian.npy mymodel.npz \
    --nao-max 26 --spin soc --structure 0 --dft-code openmx
```

| Option | Value |
|---|---|
| `--nao-max` | `nao_max` of the HamGNN config. |
| `--spin` | `none`, `collinear` (`spin_colinear: True`) or `soc` (`soc_switch: True`). |
| `--structure` | Index of the structure in `graph_data.npz`, from 0. |
| `--dft-code` | `openmx`, `abacus` or `siesta`: the DFT code of the training data. It records the angular momentum of every orbital, from HamGNN's padding tables, in the file and in `basis.json`. |

The script writes one structure to `mymodel.npz`: the lattice, the atomic
positions, H(R) and S(R) of the real orbitals of each atom (HamGNN pads every
atom to `nao_max` orbitals; the script keeps those with S_ii = 1). Replace the
prediction file by `dft` to export the DFT Hamiltonian stored in
`graph_data.npz` instead; this is how you compare the prediction with DFT in
LSQuant. If LSQuant is importable in the HamGNN environment, the same step is
`from lsquant.models.hamgnn import write_hamgnn_model`.

## 3. Run LSQuant

```json
{
  "model": {
    "kind": "hamgnn",
    "file": "/path/to/mymodel.npz",
    "params": {"velocity_mode": "berry_connection", "fermi_energy": -3.82}
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "projection": {"kind": "shift", "window": [-3, 3]},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 10.0,
                     "energy_grid": [-1.0, 1.0, 401]},
  "output": {"dir": "out/hamgnn_shc"}
}
```

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Value |
|---|---|
| `velocity_mode` | `"berry_connection"` or `"bare"`. There is no position matrix, so `"covariant"` is refused. |
| `fermi_energy` | In eV, removed as H − E_F S. HamGNN's Hamiltonian is on the absolute energy scale of the DFT code, and the export holds no Fermi energy: take the Fermi energy of the DFT run of the same structure (OpenMX prints it in Hartree as `Chemical potential`), or, for an insulator, the valence-band maximum that HamGNN's `band_cal.py` prints as `max_val`. Default 0. |

LSQuant converts Hartree to eV and Bohr to Å, puts every orbital on its atom,
and stores the spinors interleaved, so `SZ` is the Pauli matrix on the spinor
labels. HamGNN's convention for the blocks and the Bloch phase is the one of
LSQuant; no sign changes.

## Basis, positions and training data

The velocity and the spin that LSQuant builds are only as good as the basis
and the orbital positions. With HamGNN both come from the training data:

- **Basis.** The atomic orbitals of the DFT code that produced the training
  set (OpenMX pseudo-atomic orbitals, ABACUS numerical orbitals or SIESTA/HONPAS
  orbitals), with the radial cut-offs and the number of shells per element of
  that set. A prediction is meaningful only for elements and bases that were
  in the training set. Give `--dft-code`; the export then records l for every
  orbital.
- **Positions.** These orbitals are centred on the atoms, so LSQuant puts every
  orbital on its atom. The position operator is then diagonal: the
  intra-atomic dipoles (for example between s and p orbitals of one atom) are
  left out. This is the usual approximation of `berry_connection`; there is no
  `covariant` route without the position matrix of the DFT code.
- **Training data.** Note the dataset of your model (DFT code, functional,
  pseudopotentials, with or without spin–orbit coupling) next to your results.
  The energy scale, the Fermi energy and the strength of the spin–orbit
  coupling are those of that dataset.

## What the files give and what LSQuant builds

| Quantity | From HamGNN | Built by LSQuant |
|---|---|---|
| H(R) | yes, predicted | — |
| S(R) | yes, from the DFT step stored in `graph_data.npz`; a graph without overlap gives an orthogonal model (then pass HamGNN's `basis_def` to `write_hamgnn_model`) | — |
| Velocity | no | from ∂H/∂k and ∂S/∂k with every orbital on its atom (`bare` or `berry_connection`); `covariant` needs a position matrix and is refused |
| Spin | no operator; only the spinor layout of a `soc` or `collinear` model | σ on the spinor labels (interleaved after the export); a `none` model has no spin and the spin words stop with an error |
| Spin current | no | from the velocity and the spin, so only for a spinful model |

## 4. Recommendations

- Export the DFT Hamiltonian of a test structure with `dft` and run the same
  run file on both files. The difference is the error of the prediction in the
  quantity you compute.
- Use `"oracle": true` on a small supercell to compare the kernel polynomial
  result with exact diagonalization.
- A model without `soc_switch` has no spin–orbit coupling. With `--spin none`
  it has no spin operator, and the spin Hall word stops with an error; compute
  the DOS instead.
- LSQuant makes H(R) Hermitian as the mean of H(R) and H(−R)†. If it warns
  that they differ by more than 10⁻³, the prediction is not Hermitian: check
  the model before transport.
- The structure must be the one of the prediction: the script reads the
  prediction rows of each structure in the order of `graph_data.npz`, as
  HamGNN's `band_cal.py` does, and stops if the sizes do not match `nao_max`
  and the spin mode.

## References

- Y. Zhong, H. Yu, M. Su, X. Gong and H. Xiang, "Transferable equivariant graph neural networks
  for the Hamiltonians of molecules and solids", npj Comput. Mater. 9, 182
  (2023).
- HamGNN source and documentation: <https://github.com/QuantumLab-ZY/HamGNN>.
