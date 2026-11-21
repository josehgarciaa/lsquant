\page bases Bases: real space, Bloch and eigenbasis

The engine of the \ref core_engine sees only linear operators and trace vectors. Stage 5
adds three **basis adapters** that build them from one periodic tight-binding model:
- a Wannier90 `_hr.dat` file;
- the unit cell `.uc`;
- the orbital positions `.xyz`.

The same plan then runs in real space, in k space or in the per-k eigenbasis, and gives
the same moments. The k-space bases evaluate an exact trace over every k-point of the
supercell with only W vectors, W being the number of orbitals per cell.

## Model

```python
from lsquant.models.wannier import load_wannier_model
model = load_wannier_model("goldens/stage0/bases/graphene/model", seed="graphene")
```

`load_wannier_model` reads `<seed>_hr.dat` (sparse or dense layout), `<seed>.uc` (lattice
vectors, Å) and `<seed>.xyz` (orbital positions, Å). It divides each hopping by its
`ndegen` weight, and splits it over the `<seed>_wsvec.dat` shifts when that file exists. The
`WannierModel` holds `num_wann`, `uc`, `tau`, the lattice vectors `r_int`/`r_cart` and the
hopping matrices H(R). `lsquant.models.wannier` also provides H(k) and ∂H/∂k blocks
(`hk_blocks`, `dhk_blocks`) and supercell CSR builders.

## The three bases

```python
from lsquant.bases import make_basis
basis = make_basis("bloch", model, (512, 512, 1))        # "real", "bloch" or "eigen"
states = basis.states("exact")                           # or ("random", num_states=4, seed=12345)
mu = lsq.run_plan(lsq.dos_plan(M, states.number_of_states, word="1 delta"),
                  basis.hamiltonian(), window, states)
```

Every basis provides `hamiltonian()`, `operator(name)` (`"VX"`, `"VY"`, …), `states(...)`
and `window(pad=0.10)`.

| Basis | Hamiltonian | Exact trace | Disorder |
|---|---|---|---|
| `real` | the supercell CSR in the wannier2sparse cell numbering i_R = n₁ + N₁n₂ + N₁N₂n₃. It is read from `<seed>.HAM.CSR` when that file exists for this supercell (`csr_dir=`), and built from the model otherwise | all N_k·W basis vectors | diagonal V added as `OperatorSum(H, DiagonalOperator(V))` |
| `bloch` | `BlockDiagonalOperator` of the N_k blocks H(k) | W vectors (`BlockLocalStates`) for W = 0 | V in real space, sandwiched by the FFTW `BlochTransform` B: `OperatorSum(H_k, SandwichOperator(B, V))` |
| `eigen` | the eigenvalues ε_n(k), with operators rotated to U(k)† O(k) U(k) | W vectors | W = 0 only: `make_basis("eigen", ..., disorder_w=1.0)` raises `ValueError` |

**Random-phase traces.** `bloch.states("random", ...)` applies B to the same random-phase
vectors that `real.states("random", ...)` uses, so both bases trace over the same physical
vectors and give the same moments to rounding. The same holds with disorder: for
`disorder_w=1.0, disorder_seed=7`, both bases build the same box disorder in [−W/2, W/2].

**Windows.** `real.window()` is the Gershgorin bound, the legacy default. `bloch.window()`
and `eigen.window()` use the range of the block eigenvalues, which is much tighter for
models with long-range hoppings. Each window is widened by W/2 when there is disorder.

**Exact references.** `EigenBasis.exact_dos(E, eta)` and `EigenBasis.exact_sigma_xx(E, eta)`
evaluate the Gaussian-broadened DOS per cell and the Kubo-Greenwood σ_xx directly from the
eigenvalues and velocity matrix elements, with no Chebyshev expansion. σ_xx has the
normalisation of `kubo_greenwood`.

## Example: graphene in two bases

`examples/bases/graphene_bases.py` computes the DOS of graphene on a 512 × 512 supercell
(524 288 orbitals, M = 2048, 15 meV) in the Bloch basis with an exact trace (2 vectors),
and in real space with 4 random vectors:

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

The Bloch result is the exact trace over all 262 144 k-points. The real-space result
carries the stochastic error of 4 vectors, 0.6 % on average. Near the Dirac point
ρ(E)/|E| = s₀ + s₂E², and s₀ matches 2/(√3π t²) to 0.1 %. `examples/bases/README.md`
explains the output.

## Tests

| Test | What it checks | Run |
|---|---|---|
| `tests/unit/test_python_stage5.py` | on the G9 graphene model, 6 × 4 × 1 supercell: `BlochTransform` unitary; exact-trace DOS and VX–VX moments identical in the three bases (< 1e-12); random-phase moments identical in real and Bloch; disorder moments identical in real and Bloch; `exact_dos` sums to 2 per cell; `EigenBasis` rejects disorder; unknown basis names rejected. It needs no wannier2sparse binary | `PYTHONPATH=build/lsquant_core:python python3 -m pytest tests/unit` |
| `tests/dev/stage_5_check.py` | unitarity, cross-basis identities and the wannier2sparse convention, as a standalone script | `PYTHONPATH=... W2S_BIN=... python3 tests/dev/stage_5_check.py` |
| stage-5 gate | stages 1–4, then `python3 -m lsquant.apps.stage5_bases` on graphene and Haldane, supercell 48 × 48 × 1, M = 256, 40 meV, at 16 and 1 threads (details below) | `W2S_BIN=... W2S_TOOLS=... scripts/check_goldens.sh 5 build` |

The stage-5 gate runs:
- the real and Bloch bases with random phases, and the Bloch and eigen bases with an exact trace, for DOS and σ_xx;
- with disorder W = 1 eV, the real and Bloch bases;
- the FFTW convention of `BlochTransform` against a dense transform;
- the real basis against the wannier2sparse supercell CSR;
- a tolerance table (`figures/stage_5_table.py`) with three tolerances: real vs Bloch moments 1e-9, Bloch vs eigen curves 1e-8, `EigenBasis` exact values vs golden G9 1e-8.

## Golden G9

G9 is the one golden not written by the legacy LSQUANT build. It is the exact
diagonalisation of the graphene and Haldane models on the Γ-centred 48 × 48 k mesh,
computed by the independent tool `hr_exactdiag.py` of wannier2sparse. The broadening is
Gaussian, 40 meV. `goldens/stage0/PROVENANCE.md` gives the commands, versions and checksums;
`goldens/stage0/bases/gen_g9.sh` regenerates it.

The gate needs a wannier2sparse build, set in two environment variables:
- `W2S_BIN`: the executable (default `$HOME/wannier2sparse/build/wannier2sparse`);
- `W2S_TOOLS`: its `tools/` directory, which holds `hr_exactdiag.py` (default
  `$HOME/wannier2sparse/tools`).

When the tolerance table stops before its last row, the gate counts it as a failure.
`gen_g9.sh` reads the checkout root from `W2S_ROOT`.
