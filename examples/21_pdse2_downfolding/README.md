# Example 21: set up and run the PdSe₂ downfolding comparison

Compute the spin Hall conductivity of pentagonal PdSe₂ from five models of
different size: the full SIESTA Hamiltonian, the same Hamiltonian restricted
to an energy window, SIESTA reduced to 12 bands, and the Wannier90 and
PAOFLOW models of Quantum ESPRESSO. The [example page](PAGE.md) shows the
result.

The figure is made from the curves in [`data/`](data/) and in the `data/`
folder of [example 20](../20_pdse2_spin_hall_routes/README.md), so you can
plot it without running anything (step 3).

## Requirements

- LSQuant built and importable; see [Install LSQuant](../../docs/how-to/install.md).
- The SIESTA files of [example 8](../08_spin_hall_pdse2_monolayer/README.md),
  committed in the repository.
- For the Wannier90 and PAOFLOW models: steps 1 and 2 of
  [example 20](../20_pdse2_spin_hall_routes/README.md).

## Step 1: the three SIESTA models

From the repository root:

```bash
python -m lsquant run examples/20_pdse2_spin_hall_routes/run_siesta_direct.json
python -m lsquant run examples/20_pdse2_spin_hall_routes/run_siesta_projected.json
python -m lsquant run examples/21_pdse2_downfolding/run_siesta_reduced12.json
```

The three runs differ only in the `projection` section:

| Run file | `projection` | Model |
|---|---|---|
| `run_siesta_direct.json` | none | all 220 spinor orbitals |
| `run_siesta_projected.json` | `{"kind": "shift", "window": [-3, 3]}` | 220 orbitals; states outside the window moved to its edges |
| `run_siesta_reduced12.json` | `{"kind": "orthogonal", "window": [-1.77, 3.17]}` | 12 orthonormal orbitals spanning the 12 bands inside the window |

Window energies are in eV from E_F, the energy origin of the SIESTA files.

The `orthogonal` projection needs a window that holds the same number of
bands at every k point. For PdSe₂ the 12 bands from −1.75 to 2.05 eV (from
E_F) are separated from the bands below by 56 meV and from the bands above by
2.2 eV; the window edges lie in those two gaps. See
[`projection`](../../docs/reference/run-file.md#projection) in the run-file
reference.

## Step 2: the Wannier90 and PAOFLOW models

Run steps 1 and 2 of [example 20](../20_pdse2_spin_hall_routes/README.md),
then the `run_wannier_kpm.json` and `run_paoflow_kpm.json` lines of its
step 4.

## Step 3: the figure

```bash
python examples/21_pdse2_downfolding/plot.py
```

It writes `fig_downfolding.png`.

## Runtime

Total wall time of each run on an AMD EPYC 9374F:

| Run | Threads | Wall time |
|---|---|---|
| `run_siesta_reduced12.json` | 16 | 1 h 17 min, of which 58 min for the orbital dipoles of the covariant velocity |
| `run_siesta_projected.json` | 16 | 3 h 10 min |
| `run_wannier_kpm.json` | 8 | 5 min |
| `run_paoflow_kpm.json` | 8 | 2 h 30 min |

The orbital dipoles are cached in a `.lsquant-cache` directory after the first SIESTA run
with the covariant velocity, and later runs on the same SIESTA files reuse
them. The cost of the Chebyshev expansion grows with the number of bands per
k point: 220 for the window run, 12 for the reduced model.

## Files

| File | Content |
|---|---|
| `run_siesta_reduced12.json` | LSQuant run file of the reduced SIESTA model |
| `data/siesta_reduced12_kpm.dat`, `data/siesta_reduced12_exact.dat` | σᶻ_xy of the reduced model against E − E_F, KPM and exact |
| `plot.py` | the figure |
