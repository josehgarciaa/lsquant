# Disorder kinds: scalar and magnetic Anderson disorder in the Bi bilayer

`kramers.py` shows what `basis.disorder.kind` changes. It uses the SIESTA Bi(111) bilayer of
`examples/07_spin_hall_bi_bilayer/siesta`: spin-orbit coupling, a non-orthogonal basis, 72
spinor orbitals per cell. It builds a 3 × 3 × 1 supercell (648 spinor orbitals) and adds
Anderson disorder of width W = 0.5 eV, seed 7, as V_ij = (w_i + w_j) S_ij / 2, once with each
kind:

| Kind | Values on the spinor pair of each orbital |
|---|---|
| `scalar` | (u, u): a spin-independent potential |
| `magnetic` | (+h, −h): a Zeeman field along z |
| `independent` | two independent values: scalar and magnetic parts together |

For each kind it then solves H c = E S c densely and prints the largest splitting
E₂ₖ₊₁ − E₂ₖ of the eigenvalue pairs.

With spin-orbit coupling, time-reversal symmetry makes every level doubly degenerate
(Kramers pairs). A real scalar potential keeps time-reversal symmetry; a Zeeman field breaks
it.

## Run

```bash
PYTHONPATH=build/golden/lsquant_core:python python3 examples/disorder_kinds/kramers.py
```

Output on a 64-core workstation, 8 threads, 2.1 s:

```
Bi bilayer, supercell (3, 3, 1): 648 spinor orbitals, W = 0.5 eV, seed 7
               largest Kramers-pair splitting: all levels | levels with |E| < 5 eV
  clean        3.26e-04 eV | 3.26e-04 eV
  scalar       3.25e-04 eV | 3.25e-04 eV
  magnetic     1.54e-01 eV | 5.00e-02 eV
  independent  9.94e-02 eV | 4.11e-02 eV
checks: PASS
```

The script exits with status 1 if the scalar kind changes the pairing of the clean spectrum
by more than a factor of 2, or if a kind with a magnetic part splits the pairs by less than
100 times the clean splitting.

## Reading the result

- **Clean.** The pairs of the stored SIESTA Hamiltonian agree to 3.3e-4 eV, the residual of
  the stored matrices.
- **Scalar.** The disorder moves the levels, but the pairs stay paired at the clean residual
  of 3.3e-4 eV: the potential is time-reversal symmetric.
- **Magnetic.** The Zeeman field splits pairs by up to 50 meV near E_F and 0.15 eV over the
  whole spectrum, a fraction of W = 0.5 eV.
- **Independent.** This was the only layout before `kind` existed. It contains a magnetic
  part of variance W²/24, so it also splits the pairs; it is kept to reproduce earlier
  results.

In a run file:

```json
"basis": {"name": "bloch", "supercell": [24, 24, 1], "disorder": {"W": 0.5, "seed": 7, "kind": "magnetic"}}
```

`scalar` is the default. csr, wannier and generator models accept `scalar` only. The manual
page *Materials* (`docs/materials.md`) describes the three kinds.
