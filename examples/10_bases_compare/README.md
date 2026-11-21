# Example 10: one Hamiltonian, three bases (`real`, `bloch`, `eigen`)

## The physical question

A clean crystal is invariant under lattice translations, and Bloch's theorem
reduces the cost of a trace accordingly. We consider an $N_1 \times N_2$
supercell with periodic boundaries and $N_\mathrm{orb}$ orbitals per cell. Its Hilbert
space has dimension $N_k N_\mathrm{orb}$, with $N_k = N_1 N_2$. The
translations commute with $H$, so $H$ is block-diagonal in the Bloch states of
the Monkhorst–Pack mesh

$$ \mathbf k = \frac{m_1}{N_1}\,\mathbf b_1 + \frac{m_2}{N_2}\,\mathbf b_2 ,
\qquad m_i = 0, \dots, N_i - 1 , $$

with one $N_\mathrm{orb} \times N_\mathrm{orb}$ block $H(\mathbf k)$ per
$\mathbf k$. A trace takes the same value in every basis, so for any function $F$

$$ \mathrm{Tr}\,F(H) = \sum_{\mathbf k} \mathrm{tr}\,F\big(H(\mathbf k)\big)
= \sum_{\mathbf k}\sum_{n=1}^{N_\mathrm{orb}} F\big(\varepsilon_n(\mathbf k)\big) . $$

The velocity operators are also translation invariant, so the same reduction
holds for the Kubo–Greenwood trace
$\mathrm{Tr}[v_x\,\delta(E-H)\,v_x\,\delta(E-H)]$. LSQuant offers the three
forms of this identity as three bases, set by `basis.name` in the run file:

| Basis | What the engine multiplies | Unit vectors for an exact trace |
|---|---|---|
| `real` | the $N_k N_\mathrm{orb} \times N_k N_\mathrm{orb}$ real-space supercell matrix (sparse) | $N_k N_\mathrm{orb}$ local orbitals |
| `bloch` | the $N_k$ blocks $H(\mathbf k)$, $N_\mathrm{orb} \times N_\mathrm{orb}$ each | $N_\mathrm{orb}$ block-local vectors |
| `eigen` | the eigenvalues $\varepsilon_n(\mathbf k)$ of each block, with operators rotated into the eigenbasis | $N_\mathrm{orb}$ block-local vectors |

The three bases must give the same moments, up to round-off. What differs is
the cost. An exact trace in the `real` basis needs $N_k N_\mathrm{orb}$ unit
vectors. The `bloch` and `eigen` bases need only $N_\mathrm{orb}$. We check
both statements on the DOS and the Kubo–Greenwood $\sigma_{xx}$ of graphene.

The `real` basis is the one that accepts disorder and defects, which break the
translation symmetry. The `bloch` and `eigen` bases are the fast routes for
clean crystals, and the reference against which we test the `real` one. Only
the `wannier` and `siesta` model kinds build all three bases. A `csr` model is
`real` only.

## The model

We take graphene with one $p_z$ orbital per carbon atom and nearest-neighbour
hopping $t = -2.7$ eV, on a triangular lattice with $a = 2.46$ Å and
$N_\mathrm{orb} = 2$ orbitals per cell. Its Bloch Hamiltonian is

$$ H(\mathbf k) = \begin{pmatrix} 0 & t\,f(\mathbf k) \\ t\,f^*(\mathbf k) & 0 \end{pmatrix},
\qquad f(\mathbf k) = 1 + e^{-i\mathbf k\cdot\mathbf a_1} + e^{-i\mathbf k\cdot\mathbf a_2} , $$

with the two bands $E_\pm(\mathbf k) = \pm|t|\,|f(\mathbf k)|$. They span
$[-3|t|, 3|t|] = [-8.1, 8.1]$ eV, touch at the Dirac points, and have van Hove
singularities at $E = \pm|t|$.

[`make_graphene.py`](make_graphene.py) writes this model in Wannier90 format
under [`graphene/`](graphene/). The hoppings go to `graphene_hr.dat`, the
lattice vectors to `graphene.uc` and the two orbital centres to
`graphene.xyz`. The files are committed, and `python3 make_graphene.py`, run
inside this folder, regenerates them.
[Example 14](../14_wannier/README.md) explains the file format.

## The runs

All three run files use a $16 \times 16$ supercell or k mesh, so $N_k = 256$.
They share $M = 128$ moments, an exact trace and the explicit spectral window
$[-9, 9]$ eV. A moment-by-moment comparison needs the same window in all three.
By default the `real` basis takes a Gershgorin window, and the `bloch` and
`eigen` bases take the band extrema. Only `basis.name` changes between the
files.

In the `real` basis the engine works with the $512 \times 512$ supercell matrix
and traces over 512 local orbitals:

```bash
python3 -m lsquant run examples/10_bases_compare/run_real.json  --out /tmp/ex10/bases_real
```

In the `bloch` basis it works with the 256 blocks $H(\mathbf k)$ and 2 unit
vectors:

```bash
python3 -m lsquant run examples/10_bases_compare/run_bloch.json --out /tmp/ex10/bases_bloch
```

In the `eigen` basis it works with the eigenvalues of each block, again with 2
unit vectors:

```bash
python3 -m lsquant run examples/10_bases_compare/run_eigen.json --out /tmp/ex10/bases_eigen
```

Each run computes the DOS through the word `Tr[delta(E-H)]` and $\sigma_{xx}$
through `Tr[VX delta(E-H) VX delta(E-H)]`. In `bases_bloch/`, for example, it
writes:

| File | Content |
|---|---|
| `SpectralOp1graphene_bloch_KPM_M128_stateexact.chebmom1D` and `.chebmom` | DOS moments, text and binary |
| `meanSpectralOp1graphene_bloch_KPM_M128_stateexactJACKSON.dat` | $N\rho(E)$, 150 meV Jackson broadening |
| `NonEqOpVX-VXgraphene_bloch_KPM_M128x128_stateexact.chebmom2D` and `.chebmom` | rank-2 Kubo moments $\mu_{mn}$ |
| `KuboGreenWood_NonEqOpVX-VXgraphene_bloch_KPM_M128x128_stateexactJACKSON.dat` | $\sigma_{xx}(E)$, extensive, LSQuant units |
| `*.basis.json` | basis, mesh, dimension, window and operator provenance |
| `run.resolved.json` | the run file with every default filled in |

## The result

The three bases give the same DOS and the same $\sigma_{xx}$. The moments agree
to $1.3\times10^{-14}$, which is round-off, and the curve files are identical as
written. This is Bloch's theorem at work: the three runs evaluate one trace in
three bases.

The `real` basis costs more because it needs more unit vectors. With four
OpenMP threads on a shared 96-core node, the `real` run took 41 s with
512 unit vectors. The `bloch` and `eigen` runs took about 1 s each with 2, a
factor $N_k = 256$ fewer vectors. The ratio of the run times is about 40, not
256, because for runs this short the run time does not scale with the vector
count.

## Checks

[`compare_bases.py`](compare_bases.py) compares the moments and curves of the
three runs. As an independent check, it also computes the DOS moments directly
from the analytic bands on the same mesh,

$$ \mu_m = (2-\delta_{m0})\,\frac{1}{2N_k}\sum_{\mathbf k,\pm} T_m\!\Big(\frac{E_\pm(\mathbf k)-b}{a}\Big) , $$

where $a$ and $b$ are the half-width and centre of the window. We run it with:

```bash
python3 examples/10_bases_compare/compare_bases.py /tmp/ex10
```

Output on the reference run:

```text
DOS      bloch vs real: max rel |d mu| = 6.33e-15, max rel |d curve| = 0.00e+00  PASS
DOS      eigen vs real: max rel |d mu| = 7.44e-15, max rel |d curve| = 0.00e+00  PASS
sigma_xx bloch vs real: max rel |d mu| = 1.26e-14, max rel |d curve| = 0.00e+00  PASS
sigma_xx eigen vs real: max rel |d mu| = 1.23e-14, max rel |d curve| = 0.00e+00  PASS
DOS real vs analytic bands: max |d mu| = 1.13e-14  PASS
CHECK PASS
```

1. The four comparison lines compare the `bloch` and `eigen` runs with the
   `real` one. The moments differ by at most $1.26\times10^{-14}$ in relative
   terms, which is round-off, and the curves do not differ at all. The three
   bases therefore give the same DOS and the same $\sigma_{xx}$, as
   [the result](#the-result) states.
2. The last line shows that the `real` DOS moments match the analytic graphene
   bands to $1.1\times10^{-14}$. The supercell, the blocks and the closed-form
   bands describe the same spectrum.

## Runtime

The run times are in [the result](#the-result). The check took about 1 s.
