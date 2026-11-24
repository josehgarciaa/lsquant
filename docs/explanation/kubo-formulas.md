# Kubo formulas

The conductivity tensor of a non-interacting electron system follows from
the Kubo formula of linear response. Written with spectral functions of $H$,
it splits into two cases. The longitudinal conductivity $\sigma_{xx}$ is a
property of the states at the Fermi energy. The Hall conductivity
$\sigma_{xy}$ also has a contribution from all occupied states, which carries
the Berry-curvature and topological physics.

In what follows, we evaluate $\sigma_{xx}$ with the Kubo–Greenwood formula
and $\sigma_{xy}$ with the Kubo–Bastin formula. Both start from the same
two-index Chebyshev moments, explained in
[The kernel polynomial method](kpm.md), and differ only in how the moments
are combined into a curve. All curves are at zero temperature, with a
step-function Fermi distribution.

## Two-index moments

Both formulas contain two spectral functions of $H$ separated by velocity or
current operators. With the rescaled Hamiltonian $\tilde H$, one table of
moments is computed per observable,

$$
\mu_{mn} = \frac{1}{N}\,\mathrm{Tr}\big[A\, T_m(\tilde H)\, B\, T_n(\tilde H)\big],
$$

where $A$ and $B$ are the operators of the observable, for example $v_x$ and
$v_y$. The observable is the quantity to compute, written as a trace in the
`word` key of the run file. The Greenwood observable
`Tr[VX delta(E-H) VX delta(E-H)]` and the Bastin observable
`Tr[VX delta(E-H) VY dG-(E)/dE]` have their moments computed in the same way,
with different operators in place of $A$ and $B$. Each observable then allows
its own reconstruction formulas: Greenwood formulas for the first, Bastin
formulas for the second.
<!-- source: words/compiler.py:516-611 (same plan for both shapes); words/catalogue.py:129-155 SHAPE_ALLOWED_FORMULAS -->

The moments are the expensive part, and every formula below reads the same
`.chebmom2D` file. We therefore compute the moments once and reconstruct the
total, Fermi-sea and Fermi-surface curves at any broadening with
`python -m lsquant reconstruct` (see the [CLI reference](../reference/cli.md)).

## Kubo–Greenwood: the longitudinal conductivity

For the diagonal response, the Kubo formula reduces to the Kubo–Greenwood
expression

$$
\sigma_{xx}(E) = \frac{\pi \hbar e^2}{\Omega}\,
\mathrm{Tr}\big[v_x\, \delta(E - H)\, v_x\, \delta(E - H)\big],
$$

where $\Omega$ is the volume (or area, or length) of the sample. Only states
at the energy $E$ contribute. Inserting the Chebyshev series of both delta
functions gives a double sum $\sum_{mn} g_m g_n\, \mu_{mn}\, \delta_m(x)\, \delta_n(x)$,
with Jackson factors $g_m$ and the delta-function series coefficients
$\delta_m(x) = T_m(x)/(\pi\sqrt{1-x^2})$.

In the evaluation, one delta function is taken from its Chebyshev series and
the other as $-\mathrm{Im}\,G^+(E)/\pi$, which is the same function. The
`greenwood` curve is written on a coarser energy grid than the density of
states, with $M_\mathrm{eff} - 1$ uniform points on the inner 95 % of the
window; the [output-file reference](../reference/output-files.md#format)
gives the grid of every formula.
<!-- source: reconstruct.py:159-239 kubo_greenwood -->

## Kubo–Bastin: the Hall conductivity

For the off-diagonal response we use the Kubo–Bastin formula,

$$
\sigma_{xy}(E_F) = \frac{i \hbar e^2}{\Omega} \int_{-\infty}^{E_F} \! dE\;
\mathrm{Tr}\!\left[v_x\, \delta(E - H)\, v_y\, \frac{dG^+(E)}{dE}
- v_x\, \frac{dG^-(E)}{dE}\, v_y\, \delta(E - H)\right],
$$

where $G^\pm(E) = (E - H \pm i0)^{-1}$ are the retarded and advanced Green's
functions. The integral runs over all occupied states below $E_F$, not only
those at $E_F$, which is why the formula captures the quantised Hall plateau
in a gap, where no state exists at $E_F$.

The run-file observable `Tr[VX delta(E-H) VY dG-(E)/dE]` is the name of this
Bastin form, not a literal instruction. The `bastin` formula evaluates the
full Fermi-sea expression above, with both terms of the trace.
<!-- source: python/lsquant/reconstruct.py:289–291 (docstring "evaluating the full Fermi-sea Kubo-Bastin conductivity"), contraction delta(x, m0) * DgreenR(x, m1) in the docstring at :298–299 -->

The Chebyshev expansion of the Green's function is

$$
G^+(x) = -\frac{i}{\sqrt{1 - x^2}} \sum_{n} (2 - \delta_{n0})\, e^{-in\arccos x}\, T_n(\tilde H),
$$

and its energy derivative follows term by term. For Hermitian $v_x$ and $v_y$
the two terms in the trace are complex conjugates of each other, so only the
first is evaluated and twice its imaginary part is taken. The result is
integrated cumulatively in energy from the bottom of the grid at $x = -0.95$
up to each $E_F$. The lower limit of the integral is thus the edge of the
grid, and every occupied state must lie above it; the default 10 % padding of
the window leaves the band bottom inside the grid.
<!-- source: reconstruct.py:283-366 kubo_bastin -->

The `bastin` curve is extensive and uses units with $e = \hbar = 1$. For a
two-dimensional system, we divide by the sample area $A$ and multiply by
$2\pi$ to obtain $\sigma_{xy}$ in units of $e^2/h$. On a plateau of a Chern
insulator this gives the Chern number,

$$
C = \frac{2\pi\, \sigma_{xy}}{A},
$$

as in [example 03](../../examples/03_haldane_hall/README.md) for the Haldane
model.

## Fermi sea and Fermi surface

Středa split the Hall conductivity into a Fermi-surface term $\sigma^I$, built
from states at $E_F$, and a Fermi-sea term $\sigma^{II}$, built from all
occupied states [1]. The formulas `sea` and `surf` split the Bastin integrand
from the same moments, but this split is not Středa's: its two parts must not
be read as $\sigma^I$ and $\sigma^{II}$. They are defined as follows.

- `surf` is the delta–delta contraction $\pi \sum_{mn} \mathrm{Re}\,\mu_{mn}\, \delta_m(x)\, \delta_n(x)$
  (with Jackson factors). It keeps only the part of
  $\mathrm{Tr}[v_x\, \delta(E-H)\, v_y\, \delta(E-H)]$ that is symmetric under
  $x \leftrightarrow y$, and it vanishes in a gap.
- `sea` contracts $G^+ - G^-$ on one index with $2\,\mathrm{Re}\,dG^+/dE$ on
  the other and integrates over energy with the trapezoid rule.
<!-- source: reconstruct.py:460-532 kubo_bastin_surf; 369-457 kubo_bastin_sea, a_tbl = Gr - Ga and Re dGr at 424-432 -->

Středa's $\sigma^I$ contains more than this symmetric product, and in a clean
metal it is in general finite.
[Example 13](../../examples/13_bastin_sea_surface/README.md) applies the
`sea`/`surf` split to the Haldane model. The sea term carries the whole
plateau, $C = 1.017$. The surface term is zero to $6 \times 10^{-14}$ at every
energy, because in a clean isotropic crystal the $x \leftrightarrow y$-symmetric
product that `surf` keeps vanishes. The sum `sea + surf` matches `bastin` to
$2.3 \times 10^{-3}$ of its maximum; the residual is the difference between
the trapezoid rule of `sea` and the rectangle rule of `bastin`. With disorder
or anisotropy, `surf` becomes finite.

## Spin currents and the spin Hall conductivity

The spin Hall conductivity $\sigma^z_{xy}$ is the Kubo–Bastin response of a
spin current to an electric field. In the observable, $v_x$ is replaced by the
spin current operator

$$
J^z_x = \tfrac{1}{2}\left(v_x \sigma_z + \sigma_z v_x\right),
$$

where $\sigma_z$ is the Pauli matrix, without the factor $\hbar/2$. The
observable is `Tr[JXSZ delta(E-H) VY dG-(E)/dE]`, and its moments are the
same two-index table as above.
<!-- source: models/nonorthogonal.py:364-374 spin_current; metric_forms.py:125-127 spin_current_forms -->

The `shc` formula multiplies the `bastin` curve by
$2\pi/(N_\mathrm{cells} A_\mathrm{cell})$. The result is $\sigma^z_{xy}$ in
units of $e^2/h$ for the spin current $J^z_x$ defined above. In the
convention where the spin is $\hbar/2$ and the spin Hall conductivity is
quoted in units of $(\hbar/e)\,e^2/h$, the values are half as large. This
factor is not applied, so we multiply the `shc` curve by $1/2$ before
comparing with such values.
[Example 07](../../examples/07_spin_hall_bi_bilayer/README.md) computes the
spin Hall conductivity of a Bi(111) bilayer from SIESTA. In a non-orthogonal
basis, the spin and velocity operators need the overlap matrix; see
[Non-orthogonal bases](nonorthogonal-bases.md).
<!-- source: reconstruct.py:731-743 shc = 2*pi*bastin/(num_cells*cell_area); factor 1/2 derived from J = (hbar/2) * (1/2){v, sigma_z} -->

Ordinary `csr` and `wannier` runs accept `formulas: ["shc"]` when both
`reconstruction.num_cells` and `reconstruction.cell_area` are supplied.
Validation refuses missing geometry before computing moments. The metric
executors derive it from the model. The curve header follows the stored word:
a left spin current selects its Pauli component and tensor indices; other
operators are named as a normalised Bastin response without a spin Hall
label. This changes no numerical values.

## The FFT-node variants

The formulas `greenwood_fft` and `bastin_fft` evaluate the same sums on the
$M$ Chebyshev nodes instead of a uniform grid. These nodes are the
Chebyshev–Gauss points $\cos[\pi(n + 1/2)/M]$; they are generated in another
order, as $x_i = \cos[\pi(2i + 1/2)/M]$ with $i = 0, \dots, M-1$, and sorted
before use. On these nodes the double sum can be computed with a fast cosine
transform. By default the sums are evaluated with coefficient tables, and the
environment variable `LSQUANT_FFT_BACKEND=fft` selects an FFT contraction of
the same sums.
<!-- source: reconstruct.py:534-662; grid at 544-550 and 605-610; angles above pi reflect onto the odd-offset nodes, no repeats (checked numerically) -->

The price is a coarse energy axis. The curve has only $M_\mathrm{eff} - 1$
points, dense near the window edges and sparse in the centre, and for
Kubo–Bastin the cumulative integral over energy is taken on these few nodes.

> **Known limitation: `bastin_fft` is coarse.** With $M = 128$ the curve has
> only 127 points. In
> [example 13](../../examples/13_bastin_sea_surface/README.md) it reproduces
> the Hall plateau ($C = 1.020$ against 1.017), but away from the gap it
> deviates from `bastin` by up to 11 % of the maximum, and above the bands it
> ends at 0.64 instead of 0. This offset is the quadrature error of the
> cumulative energy integral on the nodes, and it decreases as
> $1/M_\mathrm{eff}$: 0.317, 0.287 and 0.158 with 255, 282 and 511 nodes. The
> largest deviation away from the gap decreases more slowly, to 5.9 % with
> 511 nodes, while $C$ converges to 1.002. `bastin_fft` is therefore suited to
> plateau values, and `bastin` to full curves. The run-file key
> `reconstruction.fft_grid` does not select this grid; only the formula names
> `greenwood_fft` and `bastin_fft` do. Validation emits a `RuntimeWarning`
> when the flag is true without a requested `*_fft` formula. The flag
> remains accepted and does not change the curves.
<!-- source: python/lsquant/reconstruct.py:534–662 -->

[Running LSQuant](../how-to/run-lsquant.md) lists the observable and formula
for each quantity, and the
[observable reference](../reference/observable-words.md) lists the observables
that LSQuant accepts.

## References

1. P. Středa, Theory of quantised Hall conductivity in two dimensions,
   J. Phys. C **15**, L717 (1982).
