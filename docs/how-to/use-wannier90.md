# Use a Wannier90 model

Start from a Quantum ESPRESSO + Wannier90 calculation with spin–orbit coupling
and end with the spin Hall conductivity of the Wannier model, computed with the
kernel polynomial method (KPM) and checked against exact diagonalisation. You
should know Wannier90 and have read [Running LSQuant](run-lsquant.md).

The run file of this guide sends the Wannier model through the same executor
as a SIESTA model. That executor builds the velocity from the full position
matrix of Wannier90 (`covariant`), the spin from the DFT spin matrices
(`exact`), and computes the density of states and the spin Hall conductivity
in the `bloch` or `eigen` basis. Every key is listed in the
[run-file reference](../reference/run-file.md#wannier90-and-paoflow-models-in-the-metric-executors).

[Example 20](../../examples/20_pdse2_spin_hall_routes/README.md) follows this
guide for a PdSe₂ monolayer and compares the result with SIESTA and PAOFLOW.

## Requirements

- LSQuant built and importable; see [Install LSQuant](install.md).
- `scipy`, which reads the unformatted `.spn` file.
- Quantum ESPRESSO with `pw2wannier90.x`, and Wannier90.
<!-- tested with Quantum ESPRESSO 7.4, Wannier90 3.1.0 and 4.0.3, scipy 1.15.2 -->

## The DFT and Wannier90 calculation

1. Run `pw.x` with spin–orbit coupling (`noncolin = .true.`,
   `lspinorb = .true.`): an scf calculation, then an nscf calculation on the
   full uniform k mesh of Wannier90, with `nosym = .true.` and
   `noinv = .true.`.

2. Ask Wannier90 for the files LSQuant reads. In `<seed>.win`:

   ```text
   spinors          = true
   write_hr         = true    ! H(R): <seed>_hr.dat
   write_xyz        = true    ! Wannier centres: <seed>_centres.xyz
   write_rmn        = true    ! position matrix: <seed>_r.dat (covariant velocity)
   write_u_matrices = true    ! U and U_dis: <seed>_u.mat, <seed>_u_dis.mat (exact spin)
   ```

   Keep `unit_cell_cart`, `dis_win_min` and `dis_win_max` in the same
   `<seed>.win`: LSQuant reads the cell and the outer window from it.

3. Ask `pw2wannier90.x` for the spin matrices as well as the usual files:

   ```text
   &INPUTPP
     write_mmn = .true., write_amn = .true., write_eig = .true.
     write_spn = .true.
   /
   ```

4. Run `wannier90.x -pp`, `pw2wannier90.x` and `wannier90.x`. Check in the
   `.wout` file that the disentanglement and the spread minimisation
   converged, and compare the Wannier bands with the DFT bands in the energy
   range you need.

The model directory then holds:

| File | Written by | Used for |
|---|---|---|
| `<seed>_hr.dat` | Wannier90 | Hamiltonian H(R) |
| `<seed>.win`, `<seed>_centres.xyz` | you, Wannier90 | cell, outer window, Wannier centres |
| `<seed>_r.dat` | Wannier90 | `velocity_mode = "covariant"` |
| `<seed>.spn`, `<seed>.eig` | `pw2wannier90.x` | `spin = "exact"` |
| `<seed>_u.mat`, `<seed>_u_dis.mat` | Wannier90 | `spin = "exact"` |

LSQuant writes `<seed>.uc` and `<seed>.xyz` (cell and centres) from
`<seed>.win` and `<seed>_centres.xyz` the first time it loads the model. To
write them yourself, call
`lsquant.models.orthogonal.write_wannier_geometry(model_dir, seed)`.

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

| Key | What to set |
|---|---|
| `velocity_mode` | `"covariant"` for the spin Hall conductivity. It adds the position matrix between Wannier functions, which `berry_connection` (centres only) and `bare` leave out. |
| `spin` | `"exact"`. With spin–orbit coupling the Wannier functions mix spin, so the Pauli matrix on the spinor labels (`"labels"`) is not the spin operator. In Example 20 it gives 1.58 instead of 1.95 e²/h in the topological gap. |
| `fermi_energy` | The Fermi energy printed by `pw.x`, in eV. Energies, `projection.window` and `fermi_energies` are then measured from it. |
| `projection` | Optional. A `shift` window around E_F keeps the Chebyshev expansion on the bands that matter. |
| `oracle` | `true` writes the exact diagonalisation on the same mesh next to the KPM curve. |

Run it with `python -m lsquant run run.json`. With a projection the KPM curve
is in `word_0.dat` and the exact one in `word_0_oracle.json`; without it they
are `shc.dat` and `exact_shc.json`.

## Check the result

- **KPM against exact.** In a gap the two must agree. With 1024 moments the
  KPM plateau of Example 20 is 1 % above the exact one (2.020 against 1.998
  e²/h); with 2048 moments it is 2.014.
- **k mesh.** In a gap the value converges fast. Between bands it depends on
  the mesh and on the broadening; compare codes only with the same
  broadening.
- **Energy origin.** For an insulator, the Fermi energy of a DFT code lies
  somewhere in the gap, at a position set by the smearing. To compare with
  another code, measure energies from the valence-band maximum.

## Common errors

| Message | Cause |
|---|---|
| `V^+V differs from the identity` | `dis_win_min`/`dis_win_max` in `<seed>.win`, or the band count of `<seed>.eig`, do not match the run that wrote `_u_dis.mat`. |
| `covariant velocity requires the position matrix (_r.dat)` | Wannier90 ran without `write_rmn = true`. |
| `<seed>_centres.xyz not found` | Wannier90 ran without `write_xyz = true`, and there is no `<seed>.xyz`. |
| `model.params.spin must be 'labels' or 'exact'` | Any other value. |
