# The kernel polynomial method

The density of states, the conductivity and most other linear-response
quantities of a non-interacting system are functions of its Hamiltonian $H$.
For a tight-binding model with $N$ orbitals, exact diagonalisation costs
$O(N^3)$ operations and $O(N^2)$ memory, which rules out the disordered
supercells of $10^6$ orbitals and more that transport calculations need. The
kernel polynomial method (KPM) avoids diagonalisation [1]. It expands each
spectral function in Chebyshev polynomials of $H$ and evaluates the expansion
with sparse matrix–vector products only.

This page follows the steps LSQuant takes. It maps $H$ onto the interval
where Chebyshev polynomials live, computes the Chebyshev moments, damps them
with the Jackson kernel, and estimates the trace exactly or with random
vectors. [How LSQuant works](unifying-idea.md) explains why every observable
reduces to such a trace. The two Kubo formulas are explained in
[Kubo formulas](kubo-formulas.md), and time-dependent words in
[Time evolution](time-evolution.md).

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
$|x| > 1$, the polynomials $T_m(x)$ grow exponentially with $m$ and the
recursion below diverges.

> **Known limitation: the window is not checked against the spectrum.**
> A window that cuts the spectrum gives diverging moments and no error
> message. The C++ engine assumes an enclosing window and does not test it
> (`lsquant/core/chebyshev/rescaled_operator.hpp:44–46`). A window given by
> hand in `expansion.window` is safest when taken from the `"gershgorin"`
> bound or from the band structure, with some margin.

A wide window costs resolution, because the same number of moments is spread
over a wider energy range. LSQuant therefore offers three sources for the
window (`expansion.window`):

- an explicit interval `[lo, hi]`;
- `"gershgorin"`, the row-sum bound $\pm\rho$ of $H$, padded by 10 %;
- `"eigen"`, the extreme eigenvalues of the Bloch blocks $H(k)$, padded by
  10 % of the band width.

The Gershgorin bound is safe, but for long-range Wannier hoppings it can be
several times wider than the true spectrum. The `bloch` and `eigen` bases
therefore use `"eigen"` by default, and the `real` basis uses `"gershgorin"`.

Anderson disorder with on-site energies drawn from $[-W/2, W/2]$ moves the
eigenvalues by at most $W/2$. In an orthogonal basis both estimates are
widened by this amount before the padding is applied. In a non-orthogonal
basis the disorder acts together with the overlap matrix $S$, and $W/2$ is no
longer a bound. The window is then widened by
$\sqrt{\mathrm{cond}(S)}\,\max|w|$, where $w$ are the on-site energies drawn
and $\mathrm{cond}(S)$ is the ratio of the largest to the smallest eigenvalue
of $S$.
<!-- moments.py:376-407 window_from_gershgorin pad 0.10; bases.py:708-720 (RealBasis.window, extra = W/2*(1+pad)), 877-902 (BlochBasis: emin-W/2, emax+W/2, then pad*width), 1047-1067 (EigenBasis); disorder_bases.py:116-129 (_window: sqrt(smax/smin)*max|w|); runfile.py:266-274 default_window_for_basis -->

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
moments. LSQuant includes the factor $(2 - \delta_{m0})/N$ when it computes
the moments, and the moment file on disk holds $\mu_m$ exactly as written
above.
<!-- moments.py:25-32 _kpm_weights_1d: w[0]=1/Ns, w[m]=2/Ns; spectral_functions.py:233 legacy_delta_table has no (2-delta) factor -->

The moments come from the three-term recursion

$$
T_0(\tilde H)|r\rangle = |r\rangle, \qquad
T_1(\tilde H)|r\rangle = \tilde H|r\rangle, \qquad
T_{m+1}(\tilde H)|r\rangle = 2\tilde H\, T_m(\tilde H)|r\rangle - T_{m-1}(\tilde H)|r\rangle .
$$

Each new moment costs one product of the sparse $\tilde H$ with a vector, and
only two vectors need to be stored at a time. This is the origin of the linear
scaling: for a sparse $H$ a matrix–vector product costs $O(N)$.

Transport needs two Chebyshev indices. The Kubo formulas contain two spectral
functions separated by operators, so LSQuant computes the rank-2 tensor

$$
\mu_{mn} = \frac{1}{N}\mathrm{Tr}\big[A\, T_m(\tilde H)\, B\, T_n(\tilde H)\big],
$$

where $A$ and $B$ are, in the Kubo formulas, velocity or current operators.
One Chebyshev series is expanded on the ket and the other on the bra, and the
two are contracted. The same procedure works for any operator that can be
applied to a vector, so a spin operator or a projector can take the place of
$A$ or $B$.
<!-- moments.py kubo_plan -->

## Truncation and the Jackson kernel

A computer keeps only $M$ moments. Cutting a Chebyshev series sharply produces
Gibbs oscillations: ripples around every sharp feature, and negative values of
a density that must be positive. The KPM remedy is to multiply each moment by
a damping factor $g_m$ that goes smoothly from 1 to 0 [1]. LSQuant uses the
Jackson kernel,

$$
g_m = \frac{(M - m + 1)\cos\phi_m + \sin\phi_m \cot\frac{\pi}{M+1}}{M + 1},
\qquad \phi_m = \frac{\pi m}{M + 1},
$$

which keeps the reconstructed density positive and turns each delta peak into
a near-Gaussian of width about $\pi/M$ in the rescaled variable $x$.
<!-- kernels.py:15-65 jackson_damping -->

> **Known limitation: `reconstruction.kernel` has no effect.** The run file
> accepts the key, but every reconstruction uses the Jackson kernel. The
> reconstruction module imports only the Jackson damping functions
> (`python/lsquant/reconstruct.py:18`), and no reconstruction routine reads
> the key. The default value, `"jackson"`, describes what is computed.

In physical units, the width $\pi/M$ becomes an energy broadening
$\Gamma \approx \pi a / M$. LSQuant uses this relation the other way round.
The run file gives a broadening, `broadening_meV`, and the reconstruction
keeps

$$
M_\mathrm{eff} = \min\!\left(\left\lceil \frac{\pi a}{\Gamma} \right\rceil,\, M\right)
$$

moments and applies the Jackson kernel of order $M_\mathrm{eff}$. The finest
broadening available is therefore $\pi a/M$. Here $\Gamma$ is a nominal width
that fixes $M_\mathrm{eff}$: the peaks have the shape of the Jackson kernel,
not of a Lorentzian of half-width $\Gamma$. A larger broadening uses fewer
moments and gives a smoother curve. The moments are independent of $\Gamma$,
so we can reconstruct one moment file at several broadenings from a single
run of the expansion.
<!-- reconstruct.py:89-96 (dos_curve), 122-150 (_truncate_and_damp_2d), 699-703 (time_correlation_curve): eta = 2*Gamma/(Emax-Emin) = Gamma/a, M_eff = min(ceil(pi/eta), M) -->

> **Known limitation: a broadening finer than $\pi a / M$ is silently
> replaced.** When $\Gamma < \pi a / M$, the reconstruction uses all $M$
> moments without a warning, and the actual broadening is $\pi a / M$, not
> the value in `broadening_meV` (`python/lsquant/reconstruct.py:89–96` and
> `:136–145`). Resolving $\Gamma$ takes at least $M = \pi a / \Gamma$
> moments. For example, a 20 eV window ($a = 10$ eV) and $\Gamma = 10$ meV
> need $M \ge 3142$.

The reconstructed curve is evaluated only on the inner 95 % of the window,
$-0.95 \le x \le 0.95$. The Chebyshev weight $1/\sqrt{1-x^2}$ diverges at
$x = \pm 1$, and stopping short of the edges keeps every grid point finite.
With the default 10 % padding of the window, the band edges lie inside this
range. The density of states is written on $30\,M_\mathrm{eff}$ uniform
points. The curve file has two columns: the energy in eV and the extensive
density of states $N\rho(E)$ in 1/eV. Dividing the second column by $N$ gives
the density of states per orbital.
<!-- reconstruct.py:25 alpha=0.95; :101 num_div = 30*max_mom; :111 prefactor = system_size*scale -->

## Exact and random-phase traces

The trace in $\mu_m$ is a sum over $N$ basis vectors. An exact trace
(`expansion.trace = "exact"`) runs the recursion once from each unit vector
$|i\rangle$, which costs $N \times M$ matrix–vector products. In the `bloch`
and `eigen` bases without disorder, every operator is block-diagonal in $k$.
LSQuant then obtains the exact trace over the whole $k$ mesh from only
$N_\mathrm{orb}$ unit vectors, one per orbital of the unit cell.
<!-- state_generator.hpp:336-351 BlockLocalStates; bases.py:856-860, 1018-1020 -->

For large systems we estimate the trace stochastically [1]. A random-phase
vector has components $r_i = e^{i\varphi_i}/\sqrt{N}$, with independent phases
$\varphi_i$ uniform in $[0, 2\pi)$. Its average satisfies
$\langle\!\langle \langle r|O|r\rangle \rangle\!\rangle = \mathrm{Tr}\,O / N$,
so $R$ such vectors (`expansion.trace = "random"`, `expansion.states = R`)
replace the $N$ unit vectors.
<!-- state_generator.hpp:98-112 RandomPhaseStates: x_i = (cos phi, sin phi)/sqrt(N), mt19937(seed) -->

The random-phase estimate is self-averaging. Its relative error is
proportional to $1/\sqrt{R N}$, so a large disordered supercell often needs
only one or a few random vectors, while a small cell needs many. The error is
largest where the density of states is small, near band edges and in gaps.
The phases come from a Mersenne Twister seeded with `expansion.seed` (default
12345), so the same seed reproduces the same moments bit for bit.
<!-- runfile.py:306, 987 seed default 12345; state_generator.hpp:106-110 -->

An exact trace with `expansion.states` smaller than $N$ gives a deterministic
local average over the first `states` orbitals, not a stochastic estimate of
the full trace. LSQuant sums over the first `states` unit vectors and divides
by their number.
<!-- bases.py:692-698 (RealBasis.states), runfile.py:1819-1831 -->

## Cost

The work is counted in matrix–vector products with $\tilde H$. For $R$ random
vectors, a rank-1 word such as the density of states needs $R\,M$ products. A
rank-2 word needs $2RM$ products and $R\,M^2$ inner products of length $N$.
For large $M$ the inner products dominate. Memory grows with the number of
Chebyshev vectors held at once, which the `expansion.left_window`,
`expansion.right_window` and `expansion.state_block` keys control
(see the [run-file reference](../reference/run-file.md)).
<!-- moments.py:287-318 estimate_engine_bytes; kubo_plan two Expand(M) + Trace -->

Both $M$ and the window enter through the resolution $\pi a / M$: halving the
window halves the number of moments needed for the same broadening, and
halving the broadening doubles it. A rank-2 word needs $R\,M^2$ inner
products, so the window matters most for conductivities.

[Example 01](../../examples/01_dos_1d_chain/README.md) compares the KPM
density of states of a chain with the analytic result, and
[example 09](../../examples/09_inspect_reconstruct/README.md) reconstructs one
set of moments at several broadenings. [How we run LSQuant](../how-to/run-lsquant.md)
shows how to set `expansion.moments`, the window and the trace in a run file.

## References

1. A. Weiße, G. Wellein, A. Alvermann and H. Fehske, The kernel polynomial
   method, Rev. Mod. Phys. **78**, 275 (2006).
