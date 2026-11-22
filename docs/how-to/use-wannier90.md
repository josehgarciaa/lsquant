# Use a Wannier90 model

Start from your own Quantum ESPRESSO + Wannier90 calculation with spin–orbit
coupling and end with its spin Hall conductivity. LSQuant reads the Wannier90
files as they are. Choose the windows, converge the Wannier functions and
check the agreement with the DFT bands in Wannier90, as usual. You should know
Wannier90 and have read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

With `model.params.velocity_mode`, `spin` or `fermi_energy` in the run file,
LSQuant computes the density of states and the spin Hall conductivity of the
clean crystal on a k-mesh; a real-space supercell or disorder stops the run.
Every key is listed in the
[run-file reference](../reference/run-file.md#wannier90-and-paoflow-models-in-the-metric-executors).
<!-- source: python/lsquant/runfile.py:1486–1499 (uses_metric_executor), 1576–1579 (clean bloch or eigen only); words: reference/run-file.md (SIESTA executors, Words row) -->

## Requirements

You need:

- LSQuant built and importable; see [Install LSQuant](install.md).
- `scipy`, which reads the unformatted `.spn` file.
<!-- tested with Quantum ESPRESSO 7.4, Wannier90 3.1.0 and 4.0.3, scipy 1.15.2 -->

## 1. Turn on the outputs in Wannier90

Add these keys to `<seed>.win`, next to your usual ones (`spinors = true`,
projections, windows, `unit_cell_cart`):

```text
write_hr         = true    ! <seed>_hr.dat: H(R)
write_xyz        = true    ! <seed>_centres.xyz: Wannier centres
write_rmn        = true    ! <seed>_r.dat: position matrix (covariant velocity)
write_u_matrices = true    ! <seed>_u.mat, <seed>_u_dis.mat (exact spin)
```

Add this key to the input of `pw2wannier90.x`:

```text
write_spn = .true.         ! <seed>.spn: DFT spin matrices (exact spin)
```

Then run `wannier90.x -pp`, `pw2wannier90.x` and `wannier90.x` as usual.
The inputs of [example 20](../../examples/20_pdse2_spin_hall_routes/README.md)
(`qe/scf.in`, `qe/nscf.in`, `qe/pdse2_proj.win`, `qe/pw2wan_spn.in`) show a
complete calculation; its `pw.x` runs set `noncolin = .true.` and
`lspinorb = .true.`, and the nscf run also `nosym = .true.` and
`noinv = .true.`.
<!-- source: examples/20_pdse2_spin_hall_routes/qe/scf.in:20–21, qe/nscf.in:20–22 -->

## 2. Transfer the model

There is nothing to run. LSQuant reads these files from the Wannier90
directory:

| File | Needed for |
|---|---|
| `<seed>_hr.dat` | Hamiltonian |
| `<seed>.win` | cell (`unit_cell_cart`), outer window (`dis_win_min`, `dis_win_max`) |
| `<seed>_centres.xyz` | Wannier centres |
| `<seed>_r.dat` | `model.params.velocity_mode = "covariant"` |
| `<seed>.spn`, `<seed>.eig`, `<seed>_u.mat`, `<seed>_u_dis.mat` | `model.params.spin = "exact"` |

The first time LSQuant reads the directory, it writes the cell and the centres
to `<seed>.uc` and `<seed>.xyz` next to these files.

## 3. Run LSQuant

This run file computes the spin Hall conductivity on a 24 × 24 k-mesh, with
exact spin and the covariant velocity:

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

The template has the settings of the projected SIESTA run of
[example 20](../../examples/20_pdse2_spin_hall_routes/README.md). It is the
exact reference path of LSQuant, and each setting does the following:

- `word` is the observable, here the spin Hall conductivity σ^z_xy: the spin
  current with spin along z and flow along x, in response to an electric
  field along y. See [Observable words](../reference/observable-words.md).
- `basis.name = "eigen"` diagonalises H(k) at every point of the k-mesh, and
  `basis.supercell` is that mesh, Γ-centred, with 1 as the third entry for a
  2D system. With `expansion.trace = "exact"`, LSQuant sums over all states,
  so the curve has no random-vector noise.
- `metric.policy = "orthogonal"`: Wannier functions are orthonormal.
- `projection` with `"kind": "shift"` moves every state outside [−3, 3] eV,
  relative to E_F, to the nearest edge of the window. The expansion then
  covers at most 7.2 eV, and needs fewer moments.
- `expansion.moments` sets the energy resolution, as the smearing does for a
  DFT density of states. With `"window": "eigen"`, LSQuant pads the projected
  spectrum by 10 %, to at most [−3.6, 3.6] eV. With 1024 moments, the finest
  broadening is πa/M ≈ 11 meV or finer, with a ≤ 3.6 eV the half-width of the
  window. When πa/M exceeds the requested 6 meV, LSQuant applies πa/M
  instead; see
  [Number of moments and broadening](use-siesta.md#number-of-moments-and-broadening).
- `reconstruction.energy_grid`, `[E_min, E_max, n]` in eV relative to E_F,
  is the grid of the exact reference: here −0.5 to 2.0 eV in 1251 points,
  2 meV apart. The KPM curve has its own energy grid.

Whether the 24 × 24 mesh is converged depends on your material: increase it
until the curve stops changing. Example 20 computes the spin Hall
conductivity of its Wannier90 model on a 300 × 300 mesh.
<!-- source: template = examples/20_pdse2_spin_hall_routes/run_siesta_projected.json; eigen + exact: use-siesta.md (Basis table); shift: python/lsquant/projection/spectral.py:119–125; window pad: python/lsquant/bases.py:1410–1416; resolution: python/lsquant/reconstruct.py:35–36. A 12-band Wannier model may span less than [−3, 3] eV, so a can be below 3.6 eV; hence "or finer". 300 × 300: examples/20 README (The routes) -->

Check the file, then run it:

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Values and meaning |
|---|---|
| `model.params.velocity_mode` | `"covariant"`: the derivative of H(k) plus the full position matrix of `<seed>_r.dat`. `"berry_connection"` (the default): the derivative of H(k) with the Wannier centres as orbital positions, without the position matrix between different Wannier functions. `"bare"`: the derivative of H(k) alone. |
| `model.params.spin` | `"exact"` (spin projected from the DFT bands) or `"labels"` (the default: Pauli matrix on the spinor labels). |
| `model.params.fermi_energy` | The Fermi energy printed by `pw.x`, in eV; the default is 0. LSQuant subtracts it from H, and energies, `projection.window` and `reconstruction.fermi_energies` are then measured from it. |
| `model.params.oracle` | `true` adds an exact reference in a run without a projection (`exact_shc.json`). A projection run writes its exact reference in any case. |
<!-- source: reference/run-file.md (Wannier90 and PAOFLOW models: velocity_mode, spin, fermi_energy defaults and meaning); oracle: apps/stage6b_siesta.py:322, projection/runfile.py:214–228 -->

With a projection, LSQuant writes the KPM curve to `word_0.dat` and the exact
result to `word_0_oracle.json`, which contains the exact curve both before and
after the projection; `projection.json` gives the difference between the KPM
curve and the exact one. Without a projection, it writes `shc.dat` and
`exact_shc.json`. The curves are σ^z_xy in units of e²/h, with the Pauli
matrix σ_z in the spin current; divide by 2 to compare with values for spin
ħ/2. See [Output files](../reference/output-files.md#projection-runs) and
[Results and units](use-siesta.md#results-and-units).
<!-- source: reference/output-files.md (Projection runs; Clean SIESTA runs); docs/explanation/kubo-formulas.md -->

The template run on the 12-band Wannier90 model of PdSe₂ in example 20 takes
about 5 min on 8 threads. The cost grows with the number of bands at each
k-point, so a small Wannier90 model is much cheaper than an atomic-orbital
model of the same material.
<!-- source: measured, examples/20_pdse2_spin_hall_routes README Runtime / data/wannier_covariant_spn.json, AMD EPYC 9374F -->

LSQuant stops the run when a file it needs is missing or inconsistent:

| Message | Cause |
|---|---|
| `V^+V differs from the identity` | The outer window in `<seed>.win` or the bands of `<seed>.eig` do not match the run that wrote `_u_dis.mat`. V = U_dis U is the matrix from the Bloch bands to the Wannier functions, read from `<seed>_u_dis.mat` and `<seed>_u.mat`. |
| `covariant velocity requires the position matrix (_r.dat)` | No `write_rmn = true`. |
| `<seed>_centres.xyz not found` | No `write_xyz = true` and no `<seed>.xyz`. |
<!-- source: V = U_dis U: reference/run-file.md (spin row) and python/lsquant/models/orthogonal.py:137 -->

## 4. Recommendations

- For the spin Hall conductivity, use
  `model.params.velocity_mode = "covariant"` and
  `model.params.spin = "exact"`. With spin–orbit coupling the Wannier
  functions mix spin, so the Pauli matrix on the spinor labels is not the
  spin operator, and the velocity from the centres alone misses the position
  matrix between Wannier functions. On the PdSe₂ model of example 20, at the
  centre of the topological gap, the spin operator matters more than the
  velocity:

  | Velocity | Spin projected from DFT (e²/h) | Pauli σ_z on the labels (e²/h) |
  |---|---|---|
  | `bare` | 1.824 | 1.583 |
  | `berry_connection` | 1.947 | 1.507 |
  | `covariant` | 1.946 | 1.577 |

  With the projected spin, `covariant` and `berry_connection` agree to 0.1 %
  and `bare` is 6 % lower.
  <!-- source: measured, examples/20_pdse2_spin_hall_routes README (The spin operator and the Wannier90 plateau) / data/wannier_covariant_spn.json key operators_at_upper_gap; 300 × 300 mesh -->
- Keep the frozen window over every band that carries the effect you study.
- Compare the KPM curve with the exact reference. In a gap they must agree;
  between bands both depend on the k-mesh and on the broadening.
- For an insulator, the Fermi energy of a DFT code lies somewhere in the gap,
  at a position set by the smearing. To compare with another code, measure
  energies from the valence-band maximum.
