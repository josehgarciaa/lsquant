# Example 13: the Hall conductivity split into Fermi-sea and Fermi-surface terms

## The physical question

The Kubo–Bastin formula gives the Hall conductivity of a non-interacting
system at zero temperature as

$$ \sigma_{xy}(E_F) \;\propto\; \int_{-\infty}^{E_F}\! dE\;
\mathrm{Tr}\Big[v_x\,\delta(E-H)\,v_y\,\frac{dG^+}{dE}
- v_x\,\frac{dG^-}{dE}\,v_y\,\delta(E-H)\Big] , $$

an integral over all occupied states. Středa split this expression into a
Fermi-surface term $\sigma^\mathrm{I}$, which involves only states at $E_F$,
and a Fermi-sea term $\sigma^\mathrm{II}$, an integral over the occupied
states [1]. In a gap only the sea term survives, and it carries the quantised
Chern part. In a metal the surface term is in general finite.

LSQuant reconstructs a split of its own, which must not be read as Středa's
$\sigma^\mathrm{I}$ and $\sigma^\mathrm{II}$:

- `surf` is the $\delta$–$\delta$ contraction of the moments, the part of
  $\mathrm{Tr}[v_x\,\delta(E-H)\,v_y\,\delta(E-H)]$ that is symmetric under
  $x \leftrightarrow y$. It needs states at the energy $E$ and vanishes in a
  gap;
- `sea` integrates over the occupied energies the part of the Bastin integrand
  that pairs $G^+ - G^-$ with $2\,\mathrm{Re}\,dG^+/dE$.
<!-- source: python/lsquant/reconstruct.py:424-429 (sea kernel), 460-520 (surf) -->

LSQuant computes one rank-2 tensor of Chebyshev moments $\mu_{mn}$ for the word
`Tr[VX delta(E-H) VY dG-(E)/dE]`. From these moments alone it rebuilds the
full Bastin curve (`bastin`), the two parts (`sea` and `surf`), and the full
curve on the Chebyshev nodes (`bastin_fft`). We do this for the Haldane model
in its Chern-insulator phase. We check that the parts add up, and we find which
part carries the quantised plateau.

## The model

We use the Haldane model on a honeycomb lattice, with nearest-neighbour
hopping $t_1 = -1$ eV, complex next-nearest-neighbour hopping
$t_2 e^{\pm i\phi}$ with $t_2 = 0.15$ eV and $\phi = \pi/2$, and zero sublattice
mass. With these parameters it is a Chern insulator with $C = +1$ and a gap
around $E = 0$. The operators are the committed ones of example 3, for 256
cells and 512 orbitals: `haldane.HAM.CSR`, `haldane.VX.CSR` and
`haldane.VY.CSR` in `goldens/stage0/haldane/operators`.

## The run

[`run_bastin_split.json`](run_bastin_split.json) evaluates the Bastin word with
128 × 128 moments, an exact trace over all 512 orbitals and the window
$[-4.5, 4.5]$ eV. It asks for four reconstructions at 50 meV broadening:

```json
"reconstruction": {
  "kernel": "jackson",
  "broadening_meV": 50.0,
  "formulas": ["bastin", "sea", "surf", "bastin_fft"],
  "fft_grid": true
}
```

The run is started with:

```bash
python3 -m lsquant run examples/13_bastin_sea_surface/run_bastin_split.json --out /tmp/ex13/bastin_split
```

The output directory holds the moments,
`NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D` and `.chebmom`, and four
curves. All four give the extensive $\sigma_{xy}(E)$ in units $e = \hbar = 1$:

| File | Formula | Grid |
|---|---|---|
| `KuboBastin_NonEqOpVX-VY…JACKSON.dat` | `bastin` | 3839 uniform points on 95 % of the window |
| `KuboBastinSea_NonEqOpVX-VY…JACKSON.dat` | `sea` | same grid |
| `KuboBastinSurf_NonEqOpVX-VY…JACKSON.dat` | `surf` | same grid |
| `KuboBastin_NonEqOpVX-VY…JACKSON_FFTgrid.dat` | `bastin_fft` | 127 Chebyshev nodes |

We obtain the Chern number by dividing $\sigma_{xy}$ by the area
$A = 256 \times 3\sqrt3/2$, in units of the nearest-neighbour distance, and
multiplying by $2\pi$: $C = 2\pi\,\sigma_{xy}/A$. Any other part can be
rebuilt later from the same moments. The sea term at 100 meV, for example,
comes from:

```bash
python3 -m lsquant reconstruct /tmp/ex13/bastin_split/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D sea 100
```

> **Known limitation: `bastin_fft` is coarse, and `fft_grid` has no effect.**
> The `*_fft` formulas evaluate the curve on the $M$ Chebyshev nodes, and the
> integrated curve keeps $M - 1$ of them, 127 here. Away from the gap it
> deviates from `bastin` by up to 11 % of the maximum, and above the bands it
> ends at 0.64 instead of 0. We reran the same model with more moments. The
> value above the bands falls exactly as $1/M_\mathrm{eff}$, so it is the
> quadrature error of the cumulative energy integral on the nodes. The largest
> deviation away from the gap falls more slowly:
>
> | $M$ | Broadening (meV) | Nodes | Value above the bands | Largest deviation | $C$ from `bastin_fft` |
> |---|---|---|---|---|---|
> | 128 | 50 | 127 | 0.637 | 11.3 % | 1.020 |
> | 256 | 50 | 255 | 0.317 | 7.4 % | 1.006 |
> | 512 | 50 | 282 | 0.287 | 7.0 % | 1.005 |
> | 512 | 25 | 511 | 0.158 | 5.9 % | 1.002 |
>
> At 50 meV the 512-moment run keeps only 283 moments, the most this broadening
> allows. The curve is selected by the formula name
> `bastin_fft`. The `fft_grid` flag is recorded in `run.resolved.json` and in
> the compiled reconstruction, but the run does not read it.
> <!-- source: python/lsquant/reconstruct.py:544-550, 605-610 (nodes cos[pi(2i+1/2)/M], sorted); words/compiler.py:721 (use_fft_grid compiled; not read in runfile.py) -->

## The result

In the gap the full Bastin curve forms a plateau with $C = 1.017$. The 1.7 %
excess over $C = 1$ comes from the finite resolution of 50 meV and 128
moments. The plateau comes entirely from the Fermi sea, with
$C_\mathrm{sea} = 1.017$.

The surface term is zero not only in the gap but at every energy. Its largest
value inside the bands is $6\times10^{-14}$, against 107.7 for the full curve.
LSQuant's `surf` keeps only the $x \leftrightarrow y$ symmetric $\delta$–$\delta$
part, and in a clean, isotropic crystal like this one that part vanishes. The
whole Hall response is therefore carried by LSQuant's `sea` term. This
statement is about `surf` only: Středa's $\sigma^\mathrm{I}$ need not vanish in
the bands. The split becomes informative when `surf` is finite, for example
with disorder or anisotropy. The run file and the check stay the same in that
case.

## Checks

[`check_bastin_split.py`](check_bastin_split.py) compares the four curves. We
run it with:

```bash
python3 examples/13_bastin_sea_surface/check_bastin_split.py /tmp/ex13/bastin_split
```

Output on the reference run:

```text
grids: bastin/sea/surf 3839 points (same grid: True), bastin_fft 127 points
[1] max |sea + surf - bastin| / max |bastin| = 2.345e-03
    (information) max |bastin_fft - bastin| / max |bastin| = 1.125e-01 over 102 FFT nodes; bastin_fft at the top of the band = 0.6369
[2] gap |E| < 0.5 eV: C(bastin_fft) = +1.0198 vs C(bastin) = +1.0170
[3] gap |E| < 0.5 eV: C(bastin) = +1.0170, C(sea) = +1.0170, C(surf) = +0.0000
    check 1: PASS
    check 2: PASS
    check 3: PASS
CHECK PASS
```

1. Sea plus surface reproduces the full Bastin curve to $2.3\times10^{-3}$ of
   its maximum, everywhere on the grid. Since `surf` is zero here, this
   compares `sea` with `bastin`. The two are different evaluations of the same
   energy integral, with the rectangle rule for `bastin` and the trapezoid rule
   for `sea`. Their difference is the quadrature error.
2. The Chebyshev-node curve gives the same plateau, $C = 1.020$ against 1.017.
   Its deviation outside the gap is the limitation described above.
3. Line [3] prints the Chern numbers in the gap $|E| < 0.5$ eV: $C = 1.0170$
   for `bastin` and for `sea`, and $C = 0$ for `surf`. The Fermi sea carries
   the whole plateau and the surface term adds nothing, as
   [the result](#the-result) explains.

## Runtime

On a shared 96-core node with one OpenMP thread, the run took 24 s, and the
check and each extra `reconstruct` took under 1 s.

## References

[1] P. Středa, J. Phys. C 15, L717 (1982).
