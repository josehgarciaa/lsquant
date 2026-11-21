\page kubo Kubo conductivities from rank-2 moments

This page covers the rank-2 part of the \ref python_layer. The Kubo plan of the
\ref core_engine computes the moments

μ_mn = Tr[A T_m(H̃) B T_n(H̃)]

for two Hermitian operators A and B, usually velocities. Python then reconstructs
from them the longitudinal (Kubo-Greenwood) and the full (Kubo-Bastin) conductivity, and the
Fermi-sea and Fermi-surface parts of Kubo-Bastin. Every reconstruction reproduces the legacy
`.dat` files byte for byte: goldens G4 (graphene) and G5 (Haldane).

## Moments

```python
plan = lsq.kubo_plan(M, left_index=0, right_index=1, word="VX delta VY delta")
mu = lsq.run_plan(plan, ham, window, states, operators=[vx, vy])
```

`left_index` and `right_index` select A and B in the `operators` list. `run_plan` applies
the legacy weights (2 − δ_{m0})(2 − δ_{n0})/R once. With `symmetrise=True`, the default, it
also sets μ_mn = (μ_mn + μ̄_nm)/2, so the tensor satisfies μ_mn = μ̄_nm.

`window_from_gershgorin(H, pad=0.10)` returns the window [−ρ(1 + pad), ρ(1 + pad)], with
ρ = max_i Σ_j |H_ij|. It is the legacy fallback when bounds are omitted.

## Reconstructions

All functions take the rank-2 moments and a broadening in meV. They return (E in eV,
σ(E)) in units e = ħ = 1, extensive: divide by the area (2D) or the volume (3D), in the
length unit of the position operator, to get a conductivity. The moments are damped with the
2D Jackson kernel, `jackson_damping_2d`.

| Function | Quantity | Grid |
|---|---|---|
| `kubo_greenwood` | σ_aa(E) ∝ Tr[A δ(E − H) A δ(E − H)], Fermi surface only | M − 1 uniform points |
| `kubo_bastin` | σ_ab(E_F), the Bastin formula with ∂G/∂E, integrated over the Fermi sea up to E_F | 30·M − 1 uniform points |
| `kubo_bastin_sea` | the Fermi-sea part, built from G^R − G^A and ∂G^R/∂E + ∂G^A/∂E, integrated | 30·M − 1 points |
| `kubo_bastin_surf` | the Fermi-surface part, Tr[A δ(E − H) B δ(E − H)] | 30·M − 1 points |
| `kubo_fft_grid(mu, b, formula="greenwood"\|"bastin")` | Greenwood or Bastin on the M Chebyshev nodes cos((2k + 1)π/(2M)) | M points |

`kubo_bastin` and `kubo_bastin_sea` + `kubo_bastin_surf` agree up to a quadrature error.
`kubo_bastin` integrates the Fermi sea with the forward rectangle rule and `kubo_bastin_sea`
with the trapezoid rule, both legacy choices kept for byte identity. On the Haldane golden
the difference is at most 0.23 % of the maximum of σ_xy, near |E| = 1 eV where σ_xy changes
fastest, and 8e-8 of it at E = 0. The unit test checks both bounds.

The Green-function coefficient tables are SymPy expressions, `GreenRetarded` and
`GreenRetardedDerivative`. `legacy_green_tables` gives the same tables in the legacy
operation order (< 1e-12 apart), and the reconstructions use it for byte identity.

**Energy axis.** All functions map the energy exactly, E = a·x + c for a window of
half-width a and centre c. The off-centre golden G4c (graphene, window [−8.5, 7.5]) checks
this in the gate. The legacy uniform-grid drivers wrote a·x + c/a instead (legacy defect D2,
evidence in `goldens/stage0/legacy_defects/offcentre_axis/`), so the G4c curve golden was
written by the new code. Its value column is byte-identical to the legacy output.

## Example: the Haldane model across its transition

`examples/kubo_python/haldane_phase.py` builds the Haldane model in NumPy on a 16 × 16
honeycomb supercell (512 sites, periodic). The parameters are t₁ = −1 eV, t₂ = 0.15 eV,
φ = π/2 and a Semenoff mass ±m, with the velocities V = i[H, X]. For two masses the script
computes the VX–VY and VX–VX moments (M = 128, exact trace) and extracts the Chern number
from the plateau of σ_xy in the gap, C = 2π⟨σ_xy⟩/A:

```bash
PYTHONPATH=build/lsquant_core:python python3 examples/kubo_python/haldane_phase.py
```
```
Haldane model, 16 x 16 cells (512 sites), M = 128, eta = 10.0 meV; transition at |m| = 0.779 eV
m = 0.00 eV: half-gap 0.779 eV, C = -1.017 (|C| expected 1); sigma_xx/A in the gap 3.1e-07, maximum in the bands 1.99
m = 1.20 eV: half-gap 0.421 eV, C = +0.008 (|C| expected 0); sigma_xx/A in the gap 2.5e-06, maximum in the bands 2.08
checks: PASS
```

Below the transition |m| = 3√3 t₂ sin φ = 0.779 eV the model is a Chern insulator, and
above it a trivial insulator. The 1.7 % excess of |C| is the finite-size and finite-M
error, the same as for the golden Haldane model (C = +1.017). The sign of C follows the
orientation of the t₂ phases. In both phases σ_xx vanishes in the gap. More details are in
`examples/kubo_python/README.md`.

`python3 -m lsquant.apps.tutorial03` rebuilds tutorial 03 (`examples/03_haldane_hall`) on the
golden Haldane operators and prints C = +1.017.

## Tests

| Test | What it checks | Run |
|---|---|---|
| `tests/unit/test_python_stage3.py` | Kubo plan with A ≠ B on a random 20 × 20 model against Tr[A T_m B T_n] from dense diagonalisation, with the legacy weights (< 1e-12); the symmetrised tensor is exactly Hermitian; the Gershgorin window; the 2D Jackson kernel; SymPy Green tables against the legacy tables; Haldane Greenwood, Bastin, sea and surface curves and graphene FFT-grid curves byte-identical to the goldens; Bastin = sea + surface within 0.3 %; σ_xx ≥ 0; input validation | `PYTHONPATH=build/lsquant_core:python python3 -m pytest tests/unit` |
| `tests/dev/stage_3_check.py` | the G4/G5 identities as a standalone script, and `from_csr` with SciPy arrays | `PYTHONPATH=... python3 tests/dev/stage_3_check.py` |
| stage-3 gate | stages 1 and 2; then `python3 -m lsquant.apps.stage3_kubo` at 16 and 1 threads. G5: Haldane moments and the four curves. G4: graphene moments with the Gershgorin window and a random-phase trace, the Bastin curve at 64 meV, sea and surface at 10 meV, both FFT-grid curves. All byte-identical. A grep keeps `lsquant/compat` out of `lsquant/core` and `lsquant/io` | `scripts/check_goldens.sh 3 build` |

`lsquant/compat/legacy_kernels.*` holds legacy formulas that only the gate programs and the
bindings use; the core stays independent of them.
