# Bases: graphene in real space and in k space

`graphene_bases.py` computes the density of states of graphene twice, in two of the stage-5
bases, from the same Wannier model: `goldens/stage0/bases/graphene/model`, with
nearest-neighbour hopping t = −1 eV and a C–C distance of 1 Å.

1. **Bloch basis, exact trace.** The Hamiltonian is block diagonal in k, with one 2 × 2
   block per k-point. `BlockLocalStates` traces over all L² k-points with W = 2 vectors.
2. **Real-space basis, stochastic trace.** The same supercell as one sparse matrix of
   2L² orbitals, traced with R random-phase vectors.

Both use M = 2048 moments on the window [−3.3, 3.3] eV and the Jackson kernel at 15 meV.

## Run

```bash
PYTHONPATH=build/lsquant_core:python python3 examples/bases/graphene_bases.py
```

Options: `--L` (512), `--moments` (2048), `--states` (4), `--seed` (12345),
`--broadening-meV` (15), `--out` (`out/bases`).

Output on a 64-core workstation (Ubuntu 22.04, g++ 11.4), 8 threads, 24 s wall time:

```
graphene 512 x 512 cells (524288 orbitals), M = 2048, eta = 15.0 meV
  Bloch basis, exact trace (2 vectors):    8.5 s
  real basis, 4 random vectors:                       14.7 s
integral of rho per cell (Bloch) = 2.0000 (2 orbitals per cell)
Dirac slope s0 from rho/|E| = s0 + s2 E^2 on 0.1 < |E| < 0.4 eV: 0.3673 /eV^2 (exact 2/(sqrt(3) pi) = 0.3676)
van Hove peaks at E = -1.000 and +1.000 eV (exact -+1)
real (random trace) vs Bloch (exact trace), mean relative deviation in the band: 6.20e-03
checks: PASS
```

The script writes `graphene_dos.dat` (E in eV, ρ per cell from both bases) and exits with
status 1 if one of the four checks fails.

## Reading the result

- **Sum rule.** ρ per cell integrates to 2, one state per orbital.
- **Dirac point.** Near E = 0, two Dirac cones give ρ(E) = 2|E|/(√3π t²) per cell to
  leading order. The fit ρ/|E| = s₀ + s₂E² on 0.1 < |E| < 0.4 eV returns s₀ = 0.3673,
  0.1 % from 0.3676. The quadratic term is the trigonal warping of the cones. The fit stays
  above 0.1 eV, because the broadening rounds the cusp at E = 0.
- **Van Hove singularities.** The saddle points at the M points of the Brillouin zone give
  logarithmic peaks at |E| = |t| = 1 eV, found at ±1.000 eV.
- **Two bases, one result.** The real-space curve differs from the Bloch curve by 0.6 % on
  average: this is the stochastic error of 4 random vectors on 524 288 orbitals. With the
  same random vectors in both bases, the moments agree to rounding (see
  `tests/unit/test_python_stage5.py`). Here the Bloch basis gives the exact trace in less
  time than four random vectors take in real space.
