# The kernel polynomial method

The density of states, the conductivity and most other linear-response
quantities of a non-interacting system are functions of its Hamiltonian $H$.
For a tight-binding model with $N$ orbitals, exact diagonalisation costs
$O(N^3)$ operations and $O(N^2)$ memory. This rules out the disordered
supercells of $10^6$ orbitals and more that transport calculations need. The
kernel polynomial method (KPM) avoids diagonalisation altogether [1]: each
spectral function is expanded in Chebyshev polynomials of $H$, and the
expansion is evaluated with sparse matrix–vector products only.

In what follows, we go through the steps of a calculation. The spectrum of
$H$ is mapped onto the interval on which Chebyshev polynomials are defined,
the Chebyshev moments are computed and damped with the Jackson kernel, and
the trace is evaluated exactly or with random vectors. Every quantity that
LSQuant computes is written as a trace of this kind. The two Kubo formulas for
the conductivity are explained in [Kubo formulas](kubo-formulas.md), and the
time-dependent quantities in [Time evolution](time-evolution.md).

Three parameters control the accuracy of a KPM result, and they play the
roles of familiar DFT convergence parameters. The number of moments $M$ sets
the energy resolution, like the smearing of a DFT density of states, and the
number of random vectors $R$ sets the statistical noise. The supercell size
sets the finite-size error, like the k-point mesh; in the `bloch` basis the
same key, `basis.supercell`, is the k-point mesh.
<!-- source: run-file.md basis.supercell: "Supercell size (real) or k-point mesh (bloch, eigen)" -->

## Rescaling the spectrum

Chebyshev polynomials $T_m(x) = \cos(m \arccos x)$ are defined on
$-1 \le x \le 1$, so we first map the spectrum of $H$ into that interval. With
a spectral window $[E_\mathrm{min}, E_\mathrm{max}]$ that contains every
eigenvalue, we define

$$
\tilde H = \frac{H - b}{a}, \qquad
a = \frac{E_\mathrm{max} - E_\mathrm{min}}{2}, \qquad
b = \frac{E_\mathrm{max} + E_\mathrm{min}}{2},
$$

where $a$ is the half-width of the window and $b$ its centre. An energy $E$
maps to $x = (E - b)/a$.

Every eigenvalue has to lie inside the window. At an eigenvalue outside it,
$|x| > 1$, the polynomials $T_m(x)$ grow exponentially with $m$, and the
recursion below diverges.

> **Known limitation: the window is not checked against the spectrum.**
> A window that cuts the spectrum gives diverging moments and no error
> message, because an enclosing window is assumed and never tested. A window
> given by hand in `expansion.window` is safest when it is taken from the
> band structure or, for an orthogonal basis, from the `"gershgorin"` bound,
> with some margin.
<!-- source: lsquant/core/chebyshev/rescaled_operator.hpp:44–46 -->

The window also sets the resolution: the same number of moments is spread
over the whole window, so a wider window gives a coarser curve. The key
`expansion.window` takes one of three values:

- an explicit interval `[lo, hi]`, in eV;
- `"gershgorin"`, the row-sum bound $\pm\rho$ of $H$, padded by 10 %;
- `"eigen"`, the extreme eigenvalues of the Bloch blocks $H(k)$, padded by
  10 % of the band width.

For an orthogonal basis the Gershgorin bound always encloses the spectrum.
For long-range Wannier hoppings, however, it can be several times wider than
the true spectrum, and the resolution drops by the same factor. The `bloch`
and `eigen` bases therefore use `"eigen"` by default, and the `real` basis
uses `"gershgorin"`.

In a `csr` run with `metric.policy = "nonorthogonal"`, the `"gershgorin"`
bound is computed from $H$ alone, whereas the expansion acts on
$S^{-1/2} H S^{-1/2}$ or $S^{-1} H$ (see
[Non-orthogonal bases](nonorthogonal-bases.md)). The eigenvalues of these
operators can lie outside the bound of $H$, so such a run needs an explicit
window, as the run files of
[example 11](../../examples/11_nonorthogonal_toy/README.md) give.
<!-- source: python/lsquant/runfile.py:1804–1810 (window_from_gershgorin(ham) on H only), 1960–1996 (LowdinOperator / GeneralisedRescaledOperator built on that window); examples/11_nonorthogonal_toy/run_*.json all set an explicit "window" -->

Anderson disorder with on-site energies drawn from $[-W/2, W/2]$ moves the
eigenvalues by at most $W/2$. In an orthogonal basis, both estimates are
widened by this amount before the padding is applied. In a non-orthogonal
basis the disorder acts together with the overlap matrix $S$, and $W/2$ is no
longer a bound. The window is then widened by
$\sqrt{\mathrm{cond}(S)}\,\max|w|$, where $w$ are the on-site energies drawn
and $\mathrm{cond}(S)$ is the ratio of the largest to the smallest eigenvalue
of $S$.
<!-- source: moments.py:376-407 window_from_gershgorin pad 0.10; bases.py:708-720 (RealBasis.window, extra = W/2*(1+pad)), 877-902 (BlochBasis: emin-W/2, emax+W/2, then pad*width), 1047-1067 (EigenBasis); disorder_bases.py:116-129 (_window: sqrt(smax/smin)*max|w|); runfile.py:266-274 default_window_for_basis -->

## Chebyshev moments

The density of states per orbital, $\rho(E) = \frac{1}{N}\mathrm{Tr}\,\delta(E - H)$,
has the Chebyshev expansion

$$
\rho(x) = \frac{1}{\pi\sqrt{1 - x^2}} \sum_{m=0}^{\infty} \mu_m\, T_m(x),
\qquad
\mu_m = \frac{2 - \delta_{m0}}{N}\,\mathrm{Tr}\, T_m(\tilde H),
$$

where the moments $\mu_m$ carry all the information about the system and the
factor $1/\sqrt{1-x^2}$ is the Chebyshev weight. The moments depend on $H$ and
on the window through $a$ and $b$, so a different window gives different
moments. The factor $(2 - \delta_{m0})/N$ is included when the moments are
computed, so the moment file on disk holds $\mu_m$ exactly as written above.
<!-- source: moments.py:25-32 _kpm_weights_1d: w[0]=1/Ns, w[m]=2/Ns; spectral_functions.py:233 legacy_delta_table has no (2-delta) factor -->

The moments are obtained from the three-term recursion

$$
T_0(\tilde H)|r\rangle = |r\rangle, \qquad
T_1(\tilde H)|r\rangle = \tilde H|r\rangle, \qquad
T_{m+1}(\tilde H)|r\rangle = 2\tilde H\, T_m(\tilde H)|r\rangle - T_{m-1}(\tilde H)|r\rangle .
$$

Each new moment costs one product of the sparse $\tilde H$ with a vector, and
only two vectors are stored at a time. For a sparse $H$ this product costs
$O(N)$ operations, which is the origin of the linear scaling of KPM with
system size.

A conductivity needs two Chebyshev indices. The Kubo formulas contain two
spectral functions separated by operators, so the two-index moments

$$
\mu_{mn} = \frac{1}{N}\mathrm{Tr}\big[A\, T_m(\tilde H)\, B\, T_n(\tilde H)\big]
$$

are computed, where $A$ and $B$ are, in the Kubo formulas, velocity or current
operators. One Chebyshev series is expanded on the ket and the other on the
bra, and the two are contracted. The same procedure works for any operator
that can be applied to a vector, so a spin operator or a projector can take
the place of $A$ or $B$.
<!-- source: moments.py kubo_plan -->

## Truncation and the Jackson kernel

A calculation keeps only $M$ moments. Cutting a Chebyshev series sharply
produces Gibbs oscillations: ripples around every sharp feature, and negative
values of a density that must be positive. The KPM remedy is to multiply each
moment by a damping factor $g_m$ that goes smoothly from 1 to 0 [1]. The
Jackson kernel is used,

$$
g_m = \frac{(M - m + 1)\cos\phi_m + \sin\phi_m \cot\frac{\pi}{M+1}}{M + 1},
\qquad \phi_m = \frac{\pi m}{M + 1},
$$

which keeps the reconstructed density positive and turns each delta peak into
a near-Gaussian of width about $\pi/M$ in the rescaled variable $x$.
<!-- source: kernels.py:15-65 jackson_damping -->

> **Known limitation: `reconstruction.kernel` has no effect.** The run file
> accepts the key, but every reconstruction uses the Jackson kernel. The
> default value, `"jackson"`, describes what is computed.
<!-- source: python/lsquant/reconstruct.py:18 imports only the Jackson damping functions; no reconstruction routine reads the key -->

In energy units, the width $\pi/M$ becomes a broadening $\eta \approx \pi a / M$.
This broadening is the energy resolution of every curve, and it plays the
role of the smearing in a DFT density of states. In the run file the relation
is used the other way round: the broadening is given in meV in
`reconstruction.broadening_meV`, and the reconstruction keeps

$$
M_\mathrm{eff} = \min\!\left(\left\lceil \frac{\pi a}{\eta} \right\rceil,\, M\right)
$$

moments and applies the Jackson kernel of order $M_\mathrm{eff}$. The finest
broadening available is therefore $\pi a/M$. Here $\eta$ is a nominal width
that fixes $M_\mathrm{eff}$: the peaks have the shape of the Jackson kernel,
not that of a Lorentzian of half-width $\eta$. A larger broadening uses fewer
moments and gives a smoother curve. The moments are independent of $\eta$, so
one moment file can be reconstructed at several broadenings from a single run
of the expansion.
<!-- source: reconstruct.py:89-96 (dos_curve), 122-150 (_truncate_and_damp_2d), 699-703 (time_correlation_curve): eta = 2*Gamma/(Emax-Emin) = Gamma/a, M_eff = min(ceil(pi/eta), M) -->

> **Rule of thumb.** Resolving a broadening $\eta$ takes at least
> $M = \pi a / \eta$ moments. For example, a 20 eV window ($a = 10$ eV) and
> $\eta = 10$ meV need $M \ge 3142$.

> **Known limitation: a broadening finer than $\pi a / M$ is silently
> replaced.** When $\eta < \pi a / M$, the reconstruction uses all $M$
> moments without a warning, and the actual broadening is $\pi a / M$, not
> the value in `broadening_meV`. A finer curve then needs more moments,
> following the rule of thumb above.
<!-- source: python/lsquant/reconstruct.py:89–96 and :136–145 -->

> **Rule of thumb.** In a finite system of $N$ orbitals, keep the number of
> moments used below $N$. Such a system has $N$ discrete levels, and as soon
> as more than $N$ moments contribute, the curve resolves them as separate
> peaks instead of a smooth density of states. In
> [Getting started](../tutorials/getting-started.md), the 256-site ring
> breaks into peaks at 20 meV (about 315 moments with $a = 2$ eV) and gives
> the smooth bulk curve at 60 meV (about 105 moments).
<!-- source: docs/tutorials/getting-started.md, "Density of states" and "Resolution from the same moments" sections (N levels, 315 moments at 20 meV, 105 at 60 meV, a = 2 eV) -->

The reconstructed curve covers only the inner 95 % of the window,
$-0.95 \le x \le 0.95$. The Chebyshev weight $1/\sqrt{1-x^2}$ diverges at
$x = \pm 1$, and stopping short of the edges keeps every grid point finite.
With the default 10 % padding of the window, the band edges lie inside this
range. The curve file has two columns: the energy in eV and the extensive
density of states $N\rho(E)$ in 1/eV. Dividing the second column by $N$ gives
the density of states per orbital. The number of points of each curve is
listed in the [output-file reference](../reference/output-files.md#format).
<!-- source: reconstruct.py:25 alpha=0.95; :101 num_div = 30*max_mom; :111 prefactor = system_size*scale -->

## Exact and random-phase traces

The trace in $\mu_m$ is a sum over $N$ basis vectors. An exact trace
(`expansion.trace = "exact"`) runs the recursion once from each unit vector
$|i\rangle$, which costs $N \times M$ matrix–vector products. In the `bloch`
and `eigen` bases without disorder, every operator is block-diagonal in $k$,
and the exact trace over the whole $k$ mesh is then obtained from only
$N_\mathrm{orb}$ unit vectors, one per orbital of the unit cell.
<!-- source: state_generator.hpp:336-351 BlockLocalStates; bases.py:856-860, 1018-1020 -->

For large systems we estimate the trace stochastically [1]. A random-phase
vector has components $r_i = e^{i\varphi_i}/\sqrt{N}$, with independent phases
$\varphi_i$ uniform in $[0, 2\pi)$. Its average satisfies
$\langle\!\langle \langle r|O|r\rangle \rangle\!\rangle = \mathrm{Tr}\,O / N$,
so $R$ such vectors (`expansion.trace = "random"`, `expansion.states = R`)
replace the $N$ unit vectors.
<!-- source: state_generator.hpp:98-112 RandomPhaseStates: x_i = (cos phi, sin phi)/sqrt(N), mt19937(seed) -->

The random vectors set the noise of the curve. The random-phase estimate is
self-averaging, and its relative error is proportional to $1/\sqrt{R N}$. A
large disordered supercell therefore often needs only one or a few random
vectors, while a small cell needs many. The error is largest where the
density of states is small, near band edges and in gaps. The phases are drawn
from a Mersenne Twister generator seeded with `expansion.seed` (default
12345), so a run with the same seed reproduces the same moments exactly.
<!-- source: runfile.py:306, 987 seed default 12345; state_generator.hpp:106-110 -->

An exact trace with `expansion.states` smaller than $N$ sums over the first
`states` unit vectors and divides by their number. The result is a
deterministic local average over these orbitals, not a stochastic estimate of
the full trace.
<!-- source: bases.py:692-698 (RealBasis.states), runfile.py:1819-1831 -->

## Cost

The work is counted in matrix–vector products with $\tilde H$. With $R$
random vectors, a quantity with one Chebyshev index, such as the density of
states, needs $R\,M$ products. A conductivity, with two indices, needs $2RM$
products and $R\,M^2$ inner products of length $N$, and for large $M$ the
inner products dominate. Memory grows with the number of Chebyshev vectors
held at once, which the keys `expansion.left_window`,
`expansion.right_window` and `expansion.state_block` control (see the
[run-file reference](../reference/run-file.md#expansion)).
<!-- source: moments.py:287-318 estimate_engine_bytes; kubo_plan two Expand(M) + Trace -->

As a scale, the recursion of $M = 512$ density-of-states moments for a
graphene model of $2 \times 10^6$ sites takes 17.9 s on one core and 1.77 s on
256 cores of a dual-socket AMD EPYC 9754 node, with a build tuned to that CPU
([performance measurements](../performance/README.md)).
<!-- source: docs/performance/README.md, "Result 1 — thread scaling" table (full recursion 17930 ms at 1 thread, 1767 ms at 256 threads; N = 2,000,000, M = 512, DOS; gcc -O3 -march=native) -->

The number of moments and the window enter the cost together, through the
resolution $\pi a / M$. Halving the window halves the number of moments
needed for the same broadening, and halving the broadening doubles it. A
conductivity needs $R\,M^2$ inner products, so a tight window pays off most
for conductivities.

[Example 01](../../examples/01_dos_1d_chain/README.md) compares the KPM
density of states of a chain with the analytic result, and
[example 09](../../examples/09_inspect_reconstruct/README.md) reconstructs one
set of moments at several broadenings. [Running LSQuant](../how-to/run-lsquant.md)
shows how to set `expansion.moments`, the window and the trace in a run file.

## References

1. A. Weiße, G. Wellein, A. Alvermann and H. Fehske, The kernel polynomial
   method, Rev. Mod. Phys. **78**, 275 (2006).
