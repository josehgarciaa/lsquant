\page python_layer Python layer: driving the core from Python

The `lsquant` Python package (`python/lsquant/`) exposes the C++ core of \ref core_engine
through a pybind11 module, `lsquant._core`, and adds three things on top of it:
- plan builders for the density of states and the velocity correlation;
- the Jackson kernel and the delta-function coefficients, written in SymPy;
- the reconstruction of the density of states from rank-1 moments.

The rank-1 reconstruction reproduces the legacy `.dat` output byte for byte (golden G3).

## Building

The module is built with the core when `LSQUANT_BUILD_PYTHON=ON`, which is the default. It needs
pybind11, which CMake finds through its config file:

```bash
pip install pybind11 numpy sympy
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -Dpybind11_DIR=$(python3 -m pybind11 --cmakedir)
cmake --build build -j 16
export PYTHONPATH=$PWD/build/lsquant_core:$PWD/python
python3 -c "import lsquant; print(lsquant._core.SpectralWindow(-2, 2).scale)"   # 0.5
```

Build with the same Python interpreter you run with. CMake takes the first `python3` on
`PATH`; pass `-DPython_EXECUTABLE=...` to choose another one.

## What the package contains

| Module | Content |
|---|---|
| `lsquant._core` | the C++ classes of the core: `SparseMatrix` (`from_csr`, `is_hermitian`), `SparseMatrixOperator`, `SpectralWindow`, `LocalStates`, `RandomPhaseStates`, `UserStates`, the plan instructions, `Plan`, `MomentEngine`, `ChebyshevMoments`, `ApplyKpmWeights`, `ReadCsrText`, `ReadChebmom`, `WriteChebmom`, `ContractRank1`, `WriteCurveText`. Type stubs are in `_core.pyi` |
| `lsquant.moments` | `dos_plan(M, R, word)`, `kubo_plan(M, *, left_index, right_index, word)`, `window_from_gershgorin(H, pad=0.10)` and `run_plan(plan, H, window, states, operators=(), batch_size=None, symmetrise=True)` |
| `lsquant.kernels` | `jackson_damping(M)`: g_m = [(M − m + 1) cos(mφ) + sin(mφ)/tan φ] φ/π, φ = π/(M + 1) |
| `lsquant.spectral_functions` | `SpectralFunction`: a SymPy expression c_m(x), turned into an (M, N) table. `DeltaFunction` is c_m(x) = cos(m arccos x)/(π√(1 − x²)). `legacy_delta_table` gives the same table in the legacy operation order, for byte identity |
| `lsquant.reconstruct` | `dos_curve(moments, broadening_meV, alpha=0.95)` → (E in eV, N·ρ(E) in 1/eV) |
| `lsquant.apps.stage2_dos` | the gate program (below) |
| `lsquant.apps.tutorial01` | tutorial 01 of `examples/` rebuilt in Python |

**Normalisation.** `dos_plan` folds the KPM weights (2 − δ_{m0})/R into its `Trace`, so
μ₀ = 1 for normalised trace vectors. `run_plan` applies `ApplyKpmWeights` to rank-2 tensors
when the plan has yet to apply them.

**Reconstruction.** `dos_curve` converts the broadening into an effective number of moments:
η = 2·broadening/(E_max − E_min) and M_eff = min(⌈π/η⌉, M). It then:
- damps the first M_eff moments with the Jackson kernel;
- evaluates the series on 30·M_eff uniform points in x ∈ [−α, α];
- returns E = x/scale + shift together with N·ρ(E), the extensive density of states of
  the legacy files.

Divide by `moments.system_size` to get ρ per site.

**Energy axis.** The window [E_min, E_max] has half-width a and centre c, so
scale = 1/a and shift = −c/a. `dos_curve` computes E = x/scale + shift = a·x − c/a:
- for a window centred on 0 (c = 0, as in every golden and in the example) this is the
  exact energy a·x, and the output is byte-identical to the legacy drivers (G3);
- for an off-centre window the axis is displaced by −c(1 + 1/a). The legacy drivers write
  a·x + c/a, displaced by −c(1 − 1/a) (legacy defect D2), so off-centre output matches
  neither the exact axis nor the legacy one.

The time-evolution PR corrects the axis to E = a·x + c and adds the off-centre golden G4c.
Until then, use a window centred on 0. `test_dos_curve_energy_axis_for_centred_and_offset_windows`
pins both cases.

## Example

`examples/python_layer/dos_square_lattice.py` builds the L × L square lattice in NumPy
(t = −1 eV, periodic) and passes it directly to the engine with `SparseMatrix.from_csr`. It then:
- computes M = 512 moments with R = 8 random-phase vectors;
- reconstructs ρ(E) at 50 meV;
- compares the result with the exact density of states of the infinite lattice,
  ρ(E) = K(√(1 − E²/16))/(2π²).

```bash
PYTHONPATH=build/lsquant_core:python python3 examples/python_layer/dos_square_lattice.py
```

The core of the script:

```python
import lsquant as lsq
from lsquant import _core

ham = _core.SparseMatrix.from_csr(n, row_pointers, columns, values)
window = _core.SpectralWindow(-4.2, 4.2)                 # must contain the spectrum
states = _core.RandomPhaseStates(n, 8, 12345)            # R = 8, fixed seed
mu = lsq.run_plan(lsq.dos_plan(512, 8, word="1 delta"), ham, window, states)
energies, values = lsq.dos_curve(mu, 50.0)               # values = N * rho(E)
rho = values / mu.system_size
```

The output and timing are in `examples/python_layer/README.md`.

## Tests

| Test | What it checks | Run |
|---|---|---|
| `tests/unit/test_python_stage2.py` | Jackson weights against the closed form; SymPy delta table against the legacy table (< 1e-14); plan structure and input validation; dos moments of a random 24 × 24 Hermitian matrix against dense diagonalisation (< 1e-12); both plans on the golden chain identical to G1 and G2; `dos_curve` on G1 byte-identical to G3; the weight inside the grid, (2/π) asin 0.95; skipped when the build omits `_core` | `PYTHONPATH=build/lsquant_core:python python3 -m pytest tests/unit` |
| `tests/dev/stage_2_check.py` | the same identities as a standalone script | `PYTHONPATH=... python3 tests/dev/stage_2_check.py` |
| stage-2 gate | stage 1 in full; then `stage2_dos` at 1 and 16 threads, (a) from G1 to G3 and (b) from the golden operators to G1 and G3, byte-identical; a grep that keeps physics names (`delta`, `DOS`, `Kubo`, `basis`) out of the code lines of `lsquant/core` and `lsquant/io` | `scripts/check_goldens.sh 2 build` |

The gate program can also be used by hand:

```bash
# reconstruct from a moment file
python3 -m lsquant.apps.stage2_dos --moments goldens/stage0/chain1d/SpectralOp1chain1dKPM_M64_stateexact.chebmom1D \
    --broadening-meV 120 --out out/
# compute the moments from CSR operators, then reconstruct
python3 -m lsquant.apps.stage2_dos --operators goldens/stage0/chain1d/operators --label chain1d \
    --window -2 2 --trace exact --states 2048 --moments 64 --broadening-meV 120 --out out/
```

See \ref testing for the gate and its tolerance mode.
