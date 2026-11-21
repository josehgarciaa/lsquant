# Kubo conductivities: the Haldane model across its topological transition

`haldane_phase.py` builds the Haldane model in NumPy and computes its Hall and longitudinal
conductivities with the rank-2 Kubo plan and the Python reconstructions.

The model:
- an L × L honeycomb supercell with periodic boundaries; L = 16 gives 512 sites;
- nearest-neighbour hopping t₁ = −1 eV;
- next-nearest-neighbour hopping t₂ = 0.15 eV with phase +φ on sublattice A and −φ on B,
  φ = π/2;
- a Semenoff mass +m on A and −m on B;
- velocities V_a = i[H, X_a], with ħ = 1 and lengths in units of the C–C distance.

For each mass the script:
1. computes μ_mn = Tr[VX T_m VY T_n] and Tr[VX T_m VX T_n] with an exact trace over all
   sites (M = 128, window [−4.5, 4.5]);
2. reconstructs σ_xy(E) with `lsquant.kubo_bastin` and σ_xx(E) with `lsquant.kubo_greenwood`
   at 10 meV;
3. takes the mean of σ_xy over |E| < 0.4 × half-gap and returns C = 2π⟨σ_xy⟩/A, where A is
   the supercell area.

## Run

```bash
PYTHONPATH=build/lsquant_core:python python3 examples/kubo_python/haldane_phase.py
```

Options: `--L` (16), `--moments` (128), `--broadening-meV` (10), `--masses` (0 1.2),
`--out` (`out/kubo_python`).

Output on a 64-core workstation (Ubuntu 22.04, g++ 11.4), 8 threads, 21 s wall time:

```
Haldane model, 16 x 16 cells (512 sites), M = 128, eta = 10.0 meV; transition at |m| = 0.779 eV
m = 0.00 eV: half-gap 0.779 eV, C = -1.017 (|C| expected 1); sigma_xx/A in the gap 3.1e-07, maximum in the bands 1.99
m = 1.20 eV: half-gap 0.421 eV, C = +0.008 (|C| expected 0); sigma_xx/A in the gap 2.5e-06, maximum in the bands 2.08
checks: PASS
```

The script writes `sigma_xy_m<m>.dat` and `sigma_xx_m<m>.dat` (E in eV, σ/A). It exits with
status 1 when |C| differs from its expected value by more than 0.05.

## Reading the result

- **Phase diagram.** The gap at K and K′ is 2·|m ∓ 3√3 t₂ sin φ|. It closes at
  |m| = 0.779 eV, where C jumps from ±1 to 0. With m = 0 the model is a Chern insulator, and
  with m = 1.2 eV a trivial insulator with a half-gap of 0.421 eV.
- **Quantisation.** |C| = 1.017 at m = 0: the plateau is quantised to 1.7 %, the
  finite-size and finite-M accuracy of this 512-site supercell. The golden Haldane model,
  with the same size and M, gives C = +1.017 (`python3 -m lsquant.apps.tutorial03`).
- **Sign.** The sign of C follows the orientation of the t₂ phases. Here the +φ hops run
  along a₁, a₂ − a₁ and −a₂, which gives C = −1; the opposite orientation gives +1.
- **σ_xx.** It vanishes in the gap of both phases, about 1e-6 of its maximum in the bands,
  as it must for an insulator.

The window is centred on 0, as the uniform-grid reconstructions require at this stage (see
*Energy axis* in `docs/kubo.md`).
