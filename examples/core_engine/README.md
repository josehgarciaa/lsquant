# Core engine: Chebyshev moments of a 1D chain

This example runs the new C++ core directly, through its stage-1 executable
`stage1_moments`, on the 1D chain of the golden set (N = 2048 sites, hopping −1,
periodic). It computes two moment tensors with an exact trace over all sites:

| word | plan | output |
|---|---|---|
| `dos` (window [−2.5, 2.5]) | `Expand(64) → Trace` | μ_m = Tr[T_m(H̃)]/N (×2 for m ≥ 1), `SpectralOp1chain1dKPM_M64_stateexact.chebmom` |
| `kubo` | `Expand(64) → Apply(VX) → MirrorToBra → Apply(VX) → Expand(64) → Trace` | μ_mn = Tr[VX T_m(H̃) VX T_n(H̃)], `NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom` |

Each word is written twice: in the unified `.chebmom` format (header with rank, dimensions,
window, number of states, word) and in the legacy `.chebmom1D`/`.chebmom2D` layout, which is
byte-identical to the goldens G1 and G2 in `goldens/stage0/chain1d/`.

`dos_from_moments.py` then rebuilds the density of states from the `.chebmom` file with the
Jackson kernel and compares it with the analytic DOS of the chain, ρ(E) = 1/(π√(4 − E²)).

## Run

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 16 --target stage1_moments
examples/core_engine/run.sh build
```

Output on a 64-core node (Ubuntu 22.04, g++ 11.4; each `stage1_moments` call takes about
1 s at 8 threads):

```
wrote examples/core_engine/out/dos_chain1d.dat: 401 points, M = 64
integral of rho over the band = 1.0000 (1 expected, minus the weight beyond |x| = 0.999)
max relative deviation from 1/(pi sqrt(4 - E^2)) for |E| < 1.5: 7.52e-03
```

The `dos` word runs on the window [−2.5, 2.5], wider than the band [−2, 2]. On exactly
[−2, 2] the chain DOS 1/(π√(1 − x²)) is the Chebyshev weight function itself, so every
moment with m ≥ 1 vanishes and the reconstruction is exact by construction. The wider
window gives non-trivial moments (μ₂ = −0.72, μ₄ = −0.66, …), and the 0.75 % deviation is
the Jackson broadening at M = 64. The `kubo` word stays on [−2, 2], so its legacy file is
byte-identical to G2.

The outputs go to `examples/core_engine/out/` (set `OUT=...` to change it). The same moments
are what `scripts/check_goldens.sh 1` compares byte for byte with G1 and G2.

See the manual page *Core engine* (`docs/core-engine.md`) for the plan language and the
classes involved.
