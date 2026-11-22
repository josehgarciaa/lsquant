# Use a Wannier90 model

Connect a Quantum ESPRESSO + Wannier90 calculation with spin–orbit coupling to
LSQuant and compute its spin Hall conductivity. LSQuant reads the Wannier90
files as they are; the Wannier model itself (windows, convergence, agreement
with the DFT bands) is yours. You should know Wannier90 and have read
[Running LSQuant](run-lsquant.md) and [Bring a model to LSQuant](bring-a-model.md).

Every key is listed in the
[run-file reference](../reference/run-file.md#wannier90-and-paoflow-models-in-the-metric-executors).

## Requirements

- LSQuant built and importable; see [Install LSQuant](install.md).
- `scipy`, which reads the unformatted `.spn` file.
<!-- tested with Quantum ESPRESSO 7.4, Wannier90 3.1.0 and 4.0.3, scipy 1.15.2 -->

## 1. Turn on in Wannier90

In `<seed>.win`, next to your usual keys (`spinors = true`, projections,
windows, `unit_cell_cart`):

```text
write_hr         = true    ! <seed>_hr.dat: H(R)
write_xyz        = true    ! <seed>_centres.xyz: Wannier centres
write_rmn        = true    ! <seed>_r.dat: position matrix (covariant velocity)
write_u_matrices = true    ! <seed>_u.mat, <seed>_u_dis.mat (exact spin)
```

In the input of `pw2wannier90.x`:

```text
write_spn = .true.         ! <seed>.spn: DFT spin matrices (exact spin)
```

Then run `wannier90.x -pp`, `pw2wannier90.x` and `wannier90.x` as usual.

## 2. Transfer the model

Nothing to run. LSQuant reads the Wannier90 directory:

| File | Used for |
|---|---|
| `<seed>_hr.dat` | Hamiltonian |
| `<seed>.win` | cell (`unit_cell_cart`), outer window (`dis_win_min`, `dis_win_max`) |
| `<seed>_centres.xyz` | Wannier centres |
| `<seed>_r.dat` | `velocity_mode = "covariant"` |
| `<seed>.spn`, `<seed>.eig`, `<seed>_u.mat`, `<seed>_u_dis.mat` | `spin = "exact"` |

On the first load it writes `<seed>.uc` and `<seed>.xyz` (cell and centres)
next to them.

## 3. Run LSQuant

```json
{
  "model": {
    "kind": "wannier",
    "dir": "/path/to/w90",
    "seed": "myseed",
    "params": {
      "velocity_mode": "covariant",
      "spin": "exact",
      "fermi_energy": -1.31,
      "oracle": true
    }
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "metric": {"policy": "orthogonal"},
  "projection": {"kind": "shift", "window": [-3, 3]},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 6.0,
                     "energy_grid": [-0.5, 2.0, 1251]},
  "output": {"dir": "out/w90_shc"}
}
```

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Value |
|---|---|
| `velocity_mode` | `"covariant"`, `"berry_connection"` (Wannier centres only) or `"bare"`. |
| `spin` | `"exact"` (spin projected from the DFT bands) or `"labels"` (Pauli matrix on the spinor labels). |
| `fermi_energy` | The Fermi energy printed by `pw.x`, in eV. Energies, `projection.window` and `fermi_energies` are measured from it. |
| `oracle` | `true` also writes the exact diagonalisation on the same k mesh. |

With a projection the KPM curve is `word_0.dat` and the exact one
`word_0_oracle.json`; without it, `shc.dat` and `exact_shc.json`.

LSQuant stops when a file it needs is missing or inconsistent:

| Message | Cause |
|---|---|
| `V^+V differs from the identity` | The outer window in `<seed>.win` or the bands of `<seed>.eig` do not match the run that wrote `_u_dis.mat`. |
| `covariant velocity requires the position matrix (_r.dat)` | No `write_rmn = true`. |
| `<seed>_centres.xyz not found` | No `write_xyz = true` and no `<seed>.xyz`. |

## 4. Recommendations

- For the spin Hall conductivity use `velocity_mode = "covariant"` and
  `spin = "exact"`. With spin–orbit coupling the Wannier functions mix spin,
  so the Pauli matrix on the spinor labels is not the spin operator, and the
  velocity from the centres alone misses the position matrix between Wannier
  functions.
- Keep the frozen window over every band that carries the effect you study.
- Compare the KPM curve with the exact one (`oracle`). In a gap they must
  agree; between bands both depend on the k mesh and on the broadening.
- For an insulator, the Fermi energy of a DFT code lies somewhere in the gap,
  at a position set by the smearing. To compare with another code, measure
  energies from the valence-band maximum.
