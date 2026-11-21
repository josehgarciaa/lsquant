\page nonorthogonal Non-orthogonal bases: the overlap metric

Atomic-orbital bases from DFT codes such as SIESTA are not orthogonal. Each one comes
with an overlap matrix S_ij = ⟨φ_i|φ_j⟩, and the Schrödinger equation reads H c = E S c.
Stage 6 adds the overlap to the engine of the \ref core_engine. It enters only through
linear operators, so the Chebyshev recursion and the trace machinery stay as they are.

## Coefficient representation

A state |ψ⟩ = Σ_i c_i |φ_i⟩ is a coefficient vector c. The matrix that acts on c for an
operator O is O_c = S⁻¹O, so:
- the Hamiltonian acting on coefficients is S⁻¹H;
- functions of it are f(S⁻¹H);
- the rescaled operator is S⁻¹(aH + bS) = aS⁻¹H + b.

A trace of a product written entirely in coefficient representation does not depend on
the basis. A plain random vector r with E[rr†] = I estimates it, and no S-weighting is
needed.

For the Kubo moments,

μ_mn = Tr[S⁻¹V_x T_m(S⁻¹H̃) S⁻¹V_y T_n(S⁻¹H̃)],

with V the Hermitian matrix of the velocity operator in the basis.

## Two exact routes

| Route | Recursion on | Operators | Trace |
|---|---|---|---|
| **Löwdin** | H_L = S^{−1/2} H S^{−1/2} | O_L = S^{−1/2} O S^{−1/2} (`LowdinOperator(isqrt, O)`) | the orthogonal plans, unchanged |
| **generalised** | S⁻¹(aH + bS) (`GeneralisedRescaledOperator(H, S, inv, window)`) | S⁻¹O (`MetricInverseProduct(inv, O)`) | ket: Expand, S⁻¹V_y; bra: S⁻¹, S⁻¹V_x, Expand, S |

S⁻¹V and S⁻¹H̃ are not Hermitian, so in the generalised route the bra side is the adjoint
of the ket side:

⟨bra_m| = [S T_m(S⁻¹H̃) S⁻¹V_x S⁻¹ r]†,   |ket_n⟩ = S⁻¹V_y T_n(S⁻¹H̃) r,   μ_mn = ⟨bra_m|ket_n⟩.

The density of states needs no S on the bra: Tr[T_m(S⁻¹H̃)] = E[r† T_m r]. The two routes
evaluate the same μ_mn and differ in cost:
- the Löwdin route applies S^{−1/2} twice per product;
- the generalised route applies S⁻¹ once per recursion step.

### The naive route

The legacy non-orthogonal drivers recursed on S⁻¹H̃ with S⁻¹V on both sides and a plain
dot product: the stage-1 plan with no S on the bra. That estimator computes
Tr[V T_m S⁻² V T_n] instead of Tr[V T_m S⁻¹ V T_n S⁻¹]. For a basis change T it gives
Tr[V T_m (TT†)⁻¹ V T_n (TT†)], which depends on the basis. It is correct only for S = I.
It gives the right density of states and wrong Kubo moments. `stage6_nonorth --route naive`
reproduces it as evidence, and the example below shows the effect on σ_xy.

## Overlap appliers

| Applier | Computes | Cost and accuracy |
|---|---|---|
| `ExactInverseOverlap(S)`, `ExactInverseSqrtOverlap(S)` | S⁻¹, S^{−1/2} by dense factorisation | exact; dimension ≤ `max_dimension` (4096 by default) |
| `IschoInverseOverlap(S, tolerance=1e-13, max_iterations=200)` | S⁻¹v iteratively | any size; `last_iterations` and `last_error` report each solve |
| `ChebyshevInverseOverlap(S, bounds, order)`, `ChebyshevInverseSqrtOverlap` | Chebyshev series of 1/s and 1/√s on the spectrum of S | any size; accuracy set by `order` and cond(S); `effective_order` after coefficient truncation |

`bounds` is a `SpectralBounds` with 0 < lower ≤ λ_min(S) and upper ≥ λ_max(S). The
Gershgorin bound of S can reach below 0, so use known eigenvalue bounds instead: the
fixtures record λ_min and λ_max in `manifest.json`. A `MetricPolicy` groups S, S⁻¹ and S^{−1/2} for later stages:
- `OrthogonalMetric(dim)` for S = I;
- `NonOrthogonalMetric(S, inv, inv_sqrt)` otherwise.

## Fixtures

`goldens/stage0/nonorth/` holds the golden chain (2048 sites) and Haldane (512 sites)
models rewritten in a non-orthogonal basis: H' = R†HR, V' = R†VR, S' = R†R. Here R = I + αK
is a local transformation on the sparsity graph of H, with α = 0.55, which gives cond(S')
between 5.2 and 6.7 over the four fixtures. Each
model comes in a real and a complex variant of K. Every trace in coefficient representation
equals the orthogonal golden (G1, G5). Each `manifest.json` records cond(S'), the spectrum
check max|eig(H', S') − eig(H)| (below 2e-14), and the checksums.
`python3 -m lsquant.apps.make_nonorth_fixtures` regenerates the fixtures with the functions
of `lsquant.models.nonorthogonal` (`make_local_transform`, `rotate_basis`, CSR I/O).

## Example: the Haldane model in a non-orthogonal basis

`examples/nonorthogonal/haldane_rotated.py` computes σ_xy of the rotated Haldane model by
each route and extracts the Chern number from the gap plateau:

```bash
PYTHONPATH=build/lsquant_core:python python3 examples/nonorthogonal/haldane_rotated.py
```

The output and the timings are in `examples/nonorthogonal/README.md`:
- the Löwdin and generalised routes reproduce the orthogonal moments and C = +1.017;
- the naive estimator does not.

## The C++ driver

`build/lsquant_core/stage6_nonorth` runs the DOS (`--word dos`) and the Kubo moments
(`--word kubo|kuboxy`) on a fixture, for any route (`--route lowdin|generalised|naive`) and
applier (`--applier exact|ischo|chebyshev`, `--order`, `--tol`). It writes the legacy
`.chebmom1D`/`.chebmom2D` file and the unified `.chebmom` file;
`--check-appliers` prints the applier residuals instead. `stage6_nonorth --help` lists every
option.

## Tests

| Test | What it checks | Run |
|---|---|---|
| `tests/unit/test_python_stage6.py` | exact appliers (S·S⁻¹ = 1, (S^{−1/2})²S = 1); Ischo and Chebyshev S⁻¹ against the exact applier (< 1e-10); Löwdin and generalised DOS on the rotated chain = G1 (< 1e-12); Löwdin and covariant generalised Kubo on the rotated Haldane model = G5; naive Kubo differs from G5 by more than 1 %; metric policies | `PYTHONPATH=build/lsquant_core:python python3 -m pytest tests/unit` |
| `tests/dev/stage_6_check.py` | spectrum invariance of the fixtures; `stage6_nonorth` DOS by both routes = G1; approximate appliers against order and tolerance; Chebyshev against exact S⁻¹ on random vectors (< 1e-10) | `PYTHONPATH=... python3 tests/dev/stage_6_check.py` |

Stage 6 adds no gate stage. Its identities run in the unit tests and in the stage-6 check,
and the gate keeps covering stages 0–5.
