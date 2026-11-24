\page time_evolution Time evolution with Chebyshev series

Stage 4 adds a time axis to the plans of the \ref core_engine. The evolution operator
U(dt) = exp(−iH dt/ħ) is applied as a Chebyshev series in H̃ with Bessel-function
coefficients. A `TimeAxis` instruction repeats the body of a plan at t_n = n·dt and
evolves the trace vectors between steps. The result is a moment tensor μ[m, n] with one
Chebyshev index and one time index. Three plans use it:
- the mean-square displacement (MSD);
- a time-evolved operator between projectors;
- a two-operator time correlation.

All three reproduce the legacy `.chebmomTD` files byte for byte: goldens G6 (ballistic
chain, MSD) and G7 (magnetic chain, spin precession). Units: energies in eV, times in fs,
ħ = 0.6582119624 eV fs (`lsquant.moments.HBAR`).

## The series

```python
dt = lsq.time_step_from_tmax(num_steps, tmax)          # dt = tmax / (num_steps - 1)
coeffs, phase = lsq.time_evolution_series(window, dt)
```

`time_evolution_series` returns two things:
- the coefficients c_m = (2 − δ_{m0})(−i)^m J_m(a·dt/ħ), computed until |c_m| ≤ 1e-15. By
  default the Bessel functions come from SciPy (`scipy.special.jv`). With
  `use_legacy_bessel=True` they come from GSL, as in the legacy code; the gate programs
  and tutorial apps use this option for byte identity. The two agree to 5.6e-16;
- the band-centre phase exp(−i c·dt/ħ).

Together they give phase · Σ_m c_m T_m(H̃) = U(dt) for every window, centred or not. The unit
test checks this at every eigenvalue of a random model to 1e-13. `TimeEvolutionFunction`
holds the same expansion as a SymPy expression; `legacy_bessel_coefficients` gives the
legacy values.

## Plans

With w_m = (2 − δ_{m0})/R and t_n = n·dt:

| Builder | μ[m, n] | Operators |
|---|---|---|
| `msd_plan(M, S, R, coeffs, velocity_index=0, phase=...)` | w_m Tr[[X, U(t_n)]† T_m(H̃) [X, U(t_n)]] | the rescaled velocity Ṽ = V/a, V = i[H, X] (below) |
| `projected_evolution_plan(M, S, R, coeffs, operator_index, projector_index, phase=...)` | w_m Tr[P U(t_n)† O U(t_n) T_m(H̃) P] | O and the projector P |
| `time_correlation_plan(M, S, R, coeffs, left_index, right_index, phase=...)` | w_m Tr[A T_m(H̃) U(t_n)† B U(t_n)] | A and B |

Each formula is checked against a dense computation with matrix exponentials (MSD to 1e-10), on a centred
and an off-centre window (`tests/unit/test_python_stage4.py`). Pass `dt=dt` to `run_plan`, so
that the tensor records the time step. Then
`time_correlation_curve(mu, broadening_meV, fermi_energy)` returns (t_n in fs, the
Jackson-damped series at E_F), that is Tr[… δ(E_F − H) …]/N for each t_n.

**The MSD velocity.** `msd_plan` builds [X, U(t)] from a derived Chebyshev series, which
needs the commutator in rescaled units, [X, H̃] = [X, H]/a. Pass the velocity as
`RescaledOperator(V, velocity_window(window))`, with V = i[H, X] in eV·(length unit).
Import `velocity_window` from `lsquant`; it returns `SpectralWindow(-a, a)`,
where a = (upper − lower)/2, so the velocity is V/a without an identity shift.
The Hamiltonian and time-evolution series keep the original physical window.
Both centred and off-centre enclosing windows give the MSD formula above.

Pitfall: `RescaledOperator(V, window)` adds −(c/a)·I for an off-centre
window and corrupts the MSD. The unit tests check agreement with dense matrix
exponentials and retain a check that this shifted velocity is wrong.

## Energy axis

From this stage on, every reconstruction maps E = a·x + c, the exact inverse of the
rescaling, for any window. The off-centre golden G4c (graphene, window [−8.5, 7.5]) checks
it in the gate. The legacy uniform-grid drivers wrote a·x + c/a (legacy defect D2), so the
G4c curve golden was written by the new code. Its value column is byte-identical to the
legacy output, and the legacy curve is kept as evidence in
`goldens/stage0/legacy_defects/offcentre_axis/`.

## Example: ballistic spreading in a clean chain

`examples/time_evolution/ballistic_msd.py` builds the periodic chain (t = −1 eV, 4096 sites)
and its velocity in NumPy. It computes the MSD and the density of states with R = 4 random
vectors (M = 256, 41 steps up to 100 fs), and divides the two, which gives
MSD(E, t) = Tr[[X, U]† δ(E − H) [X, U]]/Tr[δ(E − H)]. In a clean chain the motion is
ballistic, MSD = v(E)² t²/ħ² with v(E)² = 4 − E²:

```
   E (eV)   MSD/t^2 (a0^2/fs^2)   (4 - E^2)/hbar^2   ratio
    0.00               9.2312             9.2327   0.9998
    0.50               8.6523             8.6557   0.9996
    1.00               6.9154             6.9245   0.9987
    1.50               4.0499             4.0393   1.0026
checks: PASS
```

The run takes 1.8 s at 8 threads; `examples/time_evolution/README.md` explains the output.

The apps `python3 -m lsquant.apps.tutorial02`, `tutorial04` and `tutorial05` rebuild
tutorials 02 (localisation), 04 (spin precession) and 05 (charge transport) of `examples/`
from Python.

## Tests

| Test | What it checks | Run |
|---|---|---|
| `tests/unit/test_python_stage4.py` | the Bessel series against exp(−iE dt/ħ); the projected, correlation and MSD plans against dense matrix exponentials (centred and off-centre windows; shift-free MSD and an exact-trace run-file regression); the exact Kubo energy axis on G4c; `time_correlation_curve` against its definition | `PYTHONPATH=build/lsquant_core:python python3 -m pytest tests/unit` |
| `tests/dev/stage_4_check.py` | G6 moments and curves; GSL Bessel coefficients against SciPy (5.6e-16); `time_correlation_plan` against the legacy `inline_compute-kpm-CorrOp` when its binary is found (`LSQUANT_LEGACY_BUILD`) | `PYTHONPATH=... python3 tests/dev/stage_4_check.py` |
| stage-4 gate | stages 1–3 (now including G4c); then `python3 -m lsquant.apps.stage4_time` at 16 and 1 threads. G6: MSD moments, and its four curves from the golden and from the new moments. G7: SX, SY, SZ projected on PX. All byte-identical. A grep keeps the words time, velocity, spin, projector, bessel and gsl out of the code lines of `lsquant/core` and `lsquant/io` | `scripts/check_goldens.sh 4 build` |

The `_core` module now links GSL, for the legacy Bessel functions, and the Python layer
needs SciPy.
