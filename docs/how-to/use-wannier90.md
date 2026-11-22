# Use a Wannier90 model

Feed a Quantum ESPRESSO + Wannier90 calculation with spin–orbit coupling to
LSQuant and compute its spin Hall conductivity with the kernel polynomial
method (KPM). Turn on four outputs in Wannier90 and one in `pw2wannier90.x`,
then point a run file at the directory. LSQuant reads the files as they are;
the quality of the Wannier model (windows, convergence, band agreement) is
yours to check. You should know Wannier90 and have read
[Running LSQuant](run-lsquant.md).

The run file sends the Wannier model through the same executor as a SIESTA
model: same KPM engine, projection and exact check. Every key is listed in the
[run-file reference](../reference/run-file.md#wannier90-and-paoflow-models-in-the-metric-executors).
[Example 20](../../examples/20_pdse2_spin_hall_routes/README.md) uses it for a
PdSe₂ monolayer.

## Requirements

- LSQuant built and importable; see [Install LSQuant](install.md).
- `scipy`, which reads the unformatted `.spn` file.
<!-- tested with Quantum ESPRESSO 7.4, Wannier90 3.1.0 and 4.0.3, scipy 1.15.2 -->

## What to turn on

In `<seed>.win`, next to your usual keys (`spinors = true`, the projections,
the windows and `unit_cell_cart`):

```text
write_hr         = true    ! <seed>_hr.dat: H(R)
write_xyz        = true    ! <seed>_centres.xyz: Wannier centres
write_rmn        = true    ! <seed>_r.dat: position matrix, for the covariant velocity
write_u_matrices = true    ! <seed>_u.mat, <seed>_u_dis.mat: for the exact spin
```

In the input of `pw2wannier90.x`:

```text
write_spn = .true.         ! <seed>.spn: DFT spin matrices, for the exact spin
```

Run `wannier90.x -pp`, `pw2wannier90.x` and `wannier90.x` as usual. LSQuant
reads, from the directory of the model:

| File | Used for |
|---|---|
| `<seed>_hr.dat` | Hamiltonian |
| `<seed>.win` | cell (`unit_cell_cart`) and outer window (`dis_win_min`, `dis_win_max`) |
| `<seed>_centres.xyz` | Wannier centres |
| `<seed>_r.dat` | `velocity_mode = "covariant"` |
| `<seed>.spn`, `<seed>.eig`, `<seed>_u.mat`, `<seed>_u_dis.mat` | `spin = "exact"` |

The first time it loads the model, LSQuant writes `<seed>.uc` and
`<seed>.xyz` from `<seed>.win` and `<seed>_centres.xyz`.

## The run file

```json
{
  "model": {
    "kind": "wannier",
    "dir": "path/to/w90",
    "seed": "pdse2",
    "params": {
      "velocity_mode": "covariant",
      "spin": "exact",
      "fermi_energy": -1.3106,
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

| Key | Value |
|---|---|
| `velocity_mode` | `"covariant"`: uses the position matrix between Wannier functions. `"berry_connection"` uses the centres only, `"bare"` neither. |
| `spin` | `"exact"`: the spin operator projected from the DFT bands. `"labels"` uses the Pauli matrix on the spinor labels, which is not the spin operator once spin–orbit coupling mixes spin (1.58 instead of 1.95 e²/h in Example 20). |
| `fermi_energy` | The Fermi energy printed by `pw.x`, in eV. Energies, `projection.window` and `fermi_energies` are measured from it. |
| `oracle` | `true` also writes the exact diagonalisation on the same k mesh. |

Run it with `python -m lsquant run run.json`. With a projection, the KPM
curve is `word_0.dat` and the exact one `word_0_oracle.json`; without it,
`shc.dat` and `exact_shc.json`.

## What LSQuant checks

LSQuant stops when a file it needs is missing or inconsistent:

| Message | Cause |
|---|---|
| `V^+V differs from the identity` | The outer window in `<seed>.win` or the bands of `<seed>.eig` do not match the run that wrote `_u_dis.mat`. |
| `covariant velocity requires the position matrix (_r.dat)` | No `write_rmn = true`. |
| `<seed>_centres.xyz not found` | No `write_xyz = true` and no `<seed>.xyz`. |
| `model.params.spin must be 'labels' or 'exact'` | Any other value. |

It does not check the Wannier model itself: the convergence of the
disentanglement and of the spreads, or the agreement of the Wannier bands with
the DFT bands.
