# Non-orthogonal bases

Local-orbital DFT codes such as SIESTA expand the Kohn–Sham states in atomic
orbitals that overlap. The Hamiltonian then comes with an overlap matrix
$S_{ij} = \langle i | j \rangle \neq \delta_{ij}$, and the band energies solve
the generalised eigenvalue problem

$$
H c = E\, S c .
$$

The overlap changes the spectrum and every operator that appears in a Kubo
formula, so a calculation has to keep it throughout. In what follows, we see
how $S$ enters the Chebyshev expansion, why the velocity and spin-current
operators need it, and how the conditioning of $S$ sets the cost and the
energy resolution of a run.

## Matrices and operators

In a non-orthogonal basis, the matrix $O_{ij} = \langle i | \hat O | j \rangle$
of an operator differs from the map that the operator applies to coefficient
vectors. For a state $|\psi\rangle = \sum_j c_j |j\rangle$ with coefficients
$c$, the state $\hat O|\psi\rangle$ has coefficients $S^{-1} O c$. The
identity operator has the matrix $S$, not the unit matrix, and the product of
two operators has the matrix $A S^{-1} B$.

A trace of a product of operators and spectral functions therefore picks up
a factor $S^{-1}$ between neighbouring matrices. The physical matrices ($H$,
the velocity $V$ and the spin matrices) are kept Hermitian, and the overlap is
inserted only when the operators that act on vectors during the expansion are
built.
<!-- source: metric_forms.py:94-97 -->

## The Löwdin and generalised routes

KPM needs one operator whose spectrum is the generalised spectrum. Two are
available, selected by `metric.route`:

- the Löwdin route (`"lowdin"`, the default) uses
  $H_L = S^{-1/2} H S^{-1/2}$. It is Hermitian, and every static operator
  becomes $S^{-1/2} O S^{-1/2}$;
- the generalised route (`"generalised"`) uses $H_c = S^{-1} H$, and every
  static operator becomes $S^{-1} O$.
<!-- source: run-file.md metric.route default "lowdin" -->

The two Hamiltonians are similar, $H_L = S^{1/2} H_c S^{-1/2}$, so any
function of them obeys $f(H_L) = S^{1/2} f(H_c) S^{-1/2}$. Inside a trace the
factors $S^{\pm 1/2}$ cancel cyclically, and every observable LSQuant accepts
gives the same value on both routes.
[Example 11](../../examples/11_nonorthogonal_toy/README.md) checks this on a
chain with a known generalised spectrum.
<!-- source: metric_forms.py:148-175; rescaled_operator.hpp: GeneralisedRescaledOperator owns scale*H + shift*S and applies S^-1 on every multiplication; runfile.py:1960-1996 -->

The routes differ in what they need: the Löwdin route needs $S^{-1/2}$,
whereas the generalised route needs only $S^{-1}$, which an iterative method
can apply. On the generalised route, however, $S^{-1}H$ is not Hermitian in
the ordinary inner product, so the bra side of a two-index expansion cannot
reuse the ket recursion. $S^{-1}$ is therefore applied before the bra
expansion and $S$ after it.
<!-- source: words/compiler.py:246-266 -->

SIESTA models also accept a third route, `"anticommutator"`, which replaces
each operator by $\tfrac12(S^{-1} O + O S^{-1})$ and is labelled as an
approximation. It is exact for a single operator, but for products only if
$S$ commutes with $H$ and with every operator of the observable.
<!-- source: metric_forms.py:176-186; runfile.py:523-528 -->

## Applying $S^{-1}$ and $S^{-1/2}$

$S^{-1}$ and $S^{-1/2}$ are applied to vectors and never formed as a matrix
product with $H$. The method is chosen with `metric.applier`:

| `metric.applier` | How it applies the metric | Scaling |
|---|---|---|
| `exact` | Dense factorisation of $S$, for $S^{-1}$ and $S^{-1/2}$ | $O(N^3)$ setup, accepted only up to dimension `metric.max_dim` (default 4096) |
| `chebyshev` | Chebyshev expansion of $1/s$ or $s^{-1/2}$ on the spectrum of $S$, truncated when the coefficients fall below `metric.coeff_tol` | Linear in $N$; `metric.order` products with the sparse $S$ |
| `ischo` | Conjugate gradients preconditioned with an incomplete Cholesky factor; provides $S^{-1}$ only | Linear in $N$ per iteration |
<!-- source: lsquant/core/algebra/overlap_appliers.hpp -->

The `ischo` method provides only $S^{-1}$, while the Löwdin route needs
$S^{-1/2}$. A run file that combines them stops with the error
`metric.applier='ischo' provides S^-1 only and cannot be used with route='lowdin'`.
<!-- source: runfile.py:534-537 -->

The cost of the linear-scaling methods is set by the condition number
$\kappa = \mathrm{cond}(S) = s_\mathrm{max}/s_\mathrm{min}$. Large, diffuse
basis sets have a large $\kappa$, and every application of $\tilde H$ then
costs many products with $S$.

> **Rule of thumb.** A polynomial approximation of $1/s$ on
> $[s_\mathrm{min}, s_\mathrm{max}]$ needs a degree of order $\sqrt{\kappa}$
> times the number of digits required, and conjugate gradients need a similar
> number of iterations.

## Velocities need the overlap

The velocity operator is $\hat v = i[\hat H, \hat R]/\hbar$, with $\hat R$ the
position operator. In an orthogonal tight-binding model this reduces to the
$k$-derivative of $H(k)$. In a non-orthogonal basis the orbitals move with
$k$ as well, and the derivative of the overlap enters. With $\hbar = 1$, the
physical velocity matrix along direction $a$ is

$$
V_a = \partial_a H - \tfrac12\left(\partial_a S\, S^{-1} H + H S^{-1}\, \partial_a S\right)
+ i\left(H S^{-1} A_a - A_a S^{-1} H\right),
$$

where $\partial_a$ is the derivative with respect to $k_a$. The matrix $A_a$
is the intracell dipole matrix of the orbitals, centred at the midpoint of
each bond.
<!-- source: models/nonorthogonal.py:306-329; metric_forms.py:105-122 -->

The three SIESTA velocity modes (`model.params.velocity_mode`) keep different
parts of this expression:

- `bare` takes the derivatives of the lattice phase $e^{ik\cdot R}$ only and
  drops $A$;
- `berry_connection` takes the derivatives of the full phase
  $e^{ik\cdot(R + \tau_j - \tau_i)}$, which includes the orbital centres
  $\tau$, and drops $A$;
- `covariant` adds the dipole term, computed from the radial orbitals in the
  SIESTA `.ion.xml` files.
<!-- source: models/siesta.py:116-118 -->

The default depends on the SIESTA calculation: clean-crystal and projected
runs use `berry_connection`, and disordered-supercell runs use `covariant`.
The modes agree for the density of states. They can differ for the Hall and
spin Hall conductivities, which depend on interband matrix elements of $V$, so
a comparison of these conductivities between runs needs the same velocity
mode.
<!-- source: apps/stage6b_siesta.py:203 (default berry_connection); apps/stage6c_disorder.py:111 (default covariant) -->

The physical matrix $V_a$ is built first and only then transformed for the
chosen route, like any other operator. The derivative of the Löwdin
Hamiltonian $S^{-1/2} H S^{-1/2}$ is never used as a velocity.
<!-- source: models/nonorthogonal.py:315-318 -->

## Spin and spin currents

The spin matrix of a non-orthogonal spinor basis is $\Sigma_z = S\,\sigma_z$,
where $\sigma_z$ is the Pauli matrix acting on the spin index. The overlap of
SIESTA spinor orbitals is the same for both spins, so $S$ commutes with
$\sigma_z$, and the operator that acts on coefficients is $\sigma_z$ itself.
Only an overlap that commutes with all three Pauli matrices is accepted. The
spin current is $J^z_a = \tfrac12(V_a \sigma_z + \sigma_z V_a)$, the operator
used in the spin Hall observable of [Kubo formulas](kubo-formulas.md).
<!-- source: models/nonorthogonal.py:350-361 -->

## Disorder and conditioning

An on-site Anderson potential $w_i$ in a non-orthogonal basis has the matrix
$V_{ij} = \tfrac12 (w_i + w_j) S_{ij}$, the Hermitian average of assigning
the potential to the bra orbital and to the ket orbital. Its effect on the
spectrum is bounded by $\sqrt{\kappa}\, \max_i |w_i|$ rather than by $W/2$,
and the spectral window is padded by that bound, at the cost of energy
resolution. For the Bi DZP basis of
[example 12](../../examples/12_disorder_nonorthogonal/README.md),
$\kappa \approx 3550$ and $W = 1$ eV; the window grows to $[-68.5, 100.2]$ eV,
and 256 moments resolve only about 1 eV.
<!-- source: metric_forms.py:197-218 anderson_forms; disorder_bases.py:116-129 _window -->

In a model with spin–orbit coupling the orbital index $i$ carries spin, and
we choose how $w_i$ treats it with `basis.disorder.kind`. The default,
`scalar`, puts the same energy on both spin components of an orbital, so the
potential is spin-independent and time-reversal symmetric. `magnetic` puts
$+h_i$ on spin up and $-h_i$ on spin down, a random Zeeman field $h_i\sigma_z$.
Because $S$ is the same for both spins, this field flips no spin and conserves
$s_z$. In-plane spins dephase while $s_z$ is conserved, and for the spin Hall
conductivity the field acts as magnetic disorder. A third kind,
`independent`, reproduces earlier versions of LSQuant; the
[run-file reference](../reference/run-file.md#basis) describes it.
<!-- source: disorder_bases.py:29 spinor_disorder_values; runfile.py:965 default scalar; run-file.md basis section (independent: scalar and magnetic part, each of variance W^2/24) -->

A disordered SIESTA run checks the exact $S^{-1}$ and $S^{-1/2}$ with a probe
vector and stops if the residual exceeds `metric.tol`, whose default is
$10^{-13}$. For $\kappa \approx 3550$, the rounding errors of an exact
factorisation in double precision can reach
$\kappa\,\varepsilon \approx 8 \times 10^{-13}$. The residual measured in
example 12, $1.4 \times 10^{-13}$, is of rounding size.
<!-- source: disorder_bases.py:315-319 probe vs tol; runfile.py:297 default tol 1e-13 -->

> **Known limitation: the default `metric.tol` ignores cond(S).** For an
> overlap with $\kappa$ in the thousands, the default $10^{-13}$ can be smaller than the
> rounding error of double precision, and the probe can fail on rounding
> noise. Example 12 sets `metric.tol` to $10^{-10}$, which avoids the
> failure. Exact-applier failures now report the residual, estimated cond(S),
> and rounding floor, with advice to raise `metric.tol`; Chebyshev failures
> retain their order/coefficient advice.
<!-- source: python/lsquant/disorder_bases.py:315–319, python/lsquant/runfile.py:297 -->

The [run-file reference](../reference/run-file.md#metric) lists the `metric`
keys and the [`basis` keys](../reference/run-file.md#basis), and
[Running LSQuant](../how-to/run-lsquant.md#atomic-orbital-bases-and-the-overlap-matrix)
points to the use cases.
