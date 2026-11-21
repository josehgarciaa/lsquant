\page core_engine Core engine: Chebyshev moments as small programs

The new C++ core (`lsquant/`) computes every kernel-polynomial (KPM) quantity of LSQUANT
in one engine. A quantity is written as a short **plan**: a list of instructions that
expand a vector in Chebyshev polynomials, apply operators, switch from the ket to the bra
side, and finally take the trace. The engine runs the plan over a set of trace vectors and
returns a moment tensor. This page describes the pieces, how the `dos` and `kubo` plans
are built, an example run, and the tests.

## Building

The core is built next to the legacy drivers by default (`LSQUANT_BUILD_NEW_CORE=ON`):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 16
```

It produces the static library `lsquant_core`, the executable
`build/lsquant_core/stage1_moments`, and two test executables (see *Tests* below). It
needs a C++17 compiler, Eigen 3.4 and OpenMP.

## The pieces

| Class / file | Role |
|---|---|
| `SparseMatrix` (`core/algebra/sparse_matrix.hpp`) | complex CSR matrix: `y ← a A x + b y`, adjoint, Hermiticity check |
| `LinearOperator` (`core/algebra/linear_operator.hpp`) | the interface the engine uses: `Dimension()`, `Multiply(a, x, b, y)`, `IsHermitian()` |
| `GershgorinBounds`, `LanczosBounds` (`core/algebra/spectral_bounds.hpp`) | spectral bounds of H: guaranteed (Gershgorin, padded) or estimated (Lanczos) |
| `SpectralWindow` (`core/chebyshev/spectral_window.hpp`) | the interval [E_min, E_max] mapped onto [−1, 1]: H̃ = (H − c)/a |
| `RescaledOperator` | applies H̃ without forming it |
| `ChebyshevRecursion` | T₀v = v, T₁v = H̃v, T_{m+1}v = 2H̃T_mv − T_{m−1}v |
| `StateGenerator`: `LocalStates`, `RandomPhaseStates`, `UserStates` | the trace vectors \|r⟩: basis vectors (exact trace), random phases e^{iφ}/√N with a fixed seed (stochastic trace), or given vectors |
| `Plan` (`core/chebyshev/plan.hpp`) | the instruction list and a descriptive word |
| `MomentEngine` | runs a plan over all trace vectors and accumulates the moment tensor |
| `ChebyshevMoments` | the result: rank, dimensions, data, window, number of states, word |
| `io/csr_text`, `io/chebmom_file` | LSQUANT's CSR text format; the unified `.chebmom` format and the legacy `.chebmom1D/2D` layouts |

The core never prints, never reads an environment variable and knows nothing about a
basis or a model: the stage-1 gate checks this with a grep (`basis`, `std::cout`, `getenv`
must not appear under `lsquant/core`).

## The plan language

A plan has five instructions:

| Instruction | Effect |
|---|---|
| `Expand{M}` | replaces the live vector v by the M vectors T₀(H̃)v … T_{M−1}(H̃)v; opens one tensor index |
| `ApplyOperator{k}` | applies operator k (from the engine's operator table) to every live vector |
| `MirrorToBra{}` | the ket side is done; continue on the bra side, starting again from \|r⟩ |
| `ApplySeries{c}` | applies Σ_m c_m T_m(H̃) on the fly (used by time evolution; not executed in stage 1) |
| `Trace{w}` | μ[i_bra, i_ket] += w · ⟨bra_{i_bra} \| ket_{i_ket}⟩, summed over trace vectors |

The rank of the result equals the number of `Expand` instructions. Two plans are used in
this stage:

```cpp
// density of states: mu_m = (2 - delta_m0) / N * sum_r <r| T_m(H~) |r>
Plan dos;
dos.word = "1 delta";
dos.code = {Expand{M}, Trace{weights}};          // weights[m] = (m == 0 ? 1 : 2) / N

// velocity correlation: mu_mn = sum_r <r| VX T_n(H~) VX T_m(H~) |r>
Plan kubo;
kubo.word = "VX delta VX delta";
kubo.code = {Expand{M}, ApplyOperator{0}, MirrorToBra{}, ApplyOperator{0}, Expand{M},
             Trace{{1.0}}};
MomentEngine engine(rescaledH, {&vx}, states, batchSize);
ChebyshevMoments mu = engine.Run(kubo);
ApplyKpmWeights(mu, states.NumberOfStates(), /*symmetrise=*/true);
```

`batchSize` sets how many Chebyshev vectors are kept in memory at once on the expanded
side. It changes memory use only, not the result: the batched rank-2 tensor is
byte-identical to the unbatched one. That is checked by the gate and by the unit tests.
The legacy driver had a defect here: with `BATCH_SIZE` smaller than M it did not carry the
recursion across batches. Its wrong output is kept as evidence in
`goldens/stage0/chain1d/legacy_defects/batch16/` and is not used as a reference.

## Example

`examples/core_engine/` runs both plans on the 1D chain of the golden set (N = 2048,
M = 64, exact trace) and rebuilds the DOS:

```bash
examples/core_engine/run.sh build
```

The commands it runs are:

```bash
build/lsquant_core/stage1_moments --operators goldens/stage0/chain1d/operators \
    --label chain1d --window -2 2 --trace exact --states 2048 --moments 64 \
    --word kubo --out out/
build/lsquant_core/stage1_moments --operators goldens/stage0/chain1d/operators \
    --label chain1d --window -2.5 2.5 --trace exact --states 2048 --moments 64 \
    --word dos --out out/
```

| Option | Meaning |
|---|---|
| `--operators DIR --label L` | reads `DIR/L.HAM.CSR` (and `DIR/L.VX.CSR` for `kubo`) |
| `--window lo hi` | the spectral window; must contain the spectrum |
| `--trace exact\|random` | basis vectors 0 … states−1, or random-phase vectors |
| `--states n`, `--seed s` | number of trace vectors; seed of the random phases (default 12345) |
| `--moments M` | number of Chebyshev moments per index |
| `--word dos\|kubo` | which plan |
| `--batch B` | Chebyshev vectors per batch (default: M) |
| `--out DIR` | output directory |

Each run writes a `.chebmom` file and its legacy twin. `dos_from_moments.py` reads the
`.chebmom`, applies the Jackson kernel and compares with the analytic DOS
1/(π√(4 − E²)) of the chain: the integral is 1.0000 and the largest deviation for
|E| < 1.5 is 0.75 %. The DOS run uses the window [−2.5, 2.5]. On exactly [−2, 2] the chain
DOS is the Chebyshev weight function and all moments with m ≥ 1 vanish, so the comparison
would be exact by construction.

## Tests

| Test | What it checks | Run |
|---|---|---|
| `core_stage1_unit` (`tests/unit/test_core_stage1.cpp`) | each class against a closed form or a dense Eigen computation: CSR product and adjoint, CSR text round trip, Gershgorin/Lanczos bounds against dense eigenvalues, window mapping, T_m on a diagonal H̃ = cos(m arccos d), random-phase norm and seed determinism, local states, rank-1 trace = Σ_k T_m(ε̃_k), rank-2 Tr[A T_n A T_m] against a dense recursion, identical results for batch sizes 1, 5 and 12, rejection of malformed input | `ctest -R core` |
| `core_stage1_dense_check` (`tests/dev/stage_1_check.cpp`) | `dos` and `kubo` plans on a random 64 × 64 Hermitian model against dense diagonalisation; bit-for-bit identical at 1 and 8 threads | `ctest -R core` |
| stage-1 gate | G1 and G2 byte-identical at 1, 8 and 64 threads, G2 with `--batch 16`, two random-phase runs identical, the grep | `scripts/check_goldens.sh 1 build` |

See \ref testing for the gate and the tolerance mode.
