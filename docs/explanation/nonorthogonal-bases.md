# Non-orthogonal bases

Local-orbital DFT codes such as SIESTA expand the Kohn–Sham states in atomic
orbitals that overlap. The Hamiltonian then comes with an overlap matrix
$S_{ij} = \langle i | j \rangle \neq \delta_{ij}$, and the band energies solve
the generalised eigenvalue problem

$$
H c = E\, S c .
$$

The overlap changes the spectrum and every operator that appears in a Kubo
formula, so a calculation has to keep it throughout. This page explains how
LSQuant brings $S$ into the Chebyshev expansion and why the velocity and
spin-current operators need it. It also explains what the conditioning of $S$
costs.

## Matrices and operators

In a non-orthogonal basis, the matrix $O_{ij} = \langle i | \hat O | j \rangle$
of an operator differs from the map that the operator applies to coefficient
vectors. For a state $|\psi\rangle = \sum_j c_j |j\rangle$ with coefficients
$c$, the state $\hat O|\psi\rangle$ has coefficients $S^{-1} O c$. The
identity operator has the matrix $S$, not the unit matrix, and the product of
two operators has the matrix $A S^{-1} B$.

A trace of a word of operators and spectral functions therefore picks up a
factor $S^{-1}$ between neighbouring matrices. LSQuant keeps the physical
matrices Hermitian: $H$, the velocity $V$ and the spin matrices. It inserts
the metric only when it builds the operators that the engine applies.
<!-- metric_forms.py:94-97 -->

## The Löwdin and generalised routes

KPM needs one operator whose spectrum is the generalised spectrum. LSQuant
offers two, selected by `metric.route`:

- the **Löwdin route** (`"lowdin"`) uses $H_L = S^{-1/2} H S^{-1/2}$. It is
  Hermitian, and every static operator becomes $S^{-1/2} O S^{-1/2}$;
- the **generalised route** (`"generalised"`) uses $H_c = S^{-1} H$, and every
  static operator becomes $S^{-1} O$.

The two Hamiltonians are similar, $H_L = S^{1/2} H_c S^{-1/2}$, so any
function of them obeys $f(H_L) = S^{1/2} f(H_c) S^{-1/2}$. Inside a trace the
factors $S^{\pm 1/2}$ cancel cyclically. Every word in the catalogue therefore
gives the same value on both routes.
[Example 11](../../examples/11_nonorthogonal_toy/README.md) checks this on a
chain with a known generalised spectrum.
<!-- metric_forms.py:148-175; rescaled_operator.hpp: GeneralisedRescaledOperator owns scale*H + shift*S and applies S^-1 on every multiplication; runfile.py:1960-1996 -->

The routes differ in what they need. The Löwdin route needs $S^{-1/2}$. The
generalised route needs only $S^{-1}$, which an iterative method can apply.
On the generalised route, however, $S^{-1}H$ is not Hermitian in the ordinary
inner product. The bra side of a rank-2 expansion then cannot reuse the ket
recursion. LSQuant handles this by applying $S^{-1}$ before the bra expansion
and $S$ after it.
<!-- words/compiler.py:246-266 -->

SIESTA models also accept a third route, `"anticommutator"`, which replaces
each operator by $\tfrac12(S^{-1} O + O S^{-1})$. It is labelled as an
approximation. For a single operator it is exact. For products it is exact
only if $S$ commutes with $H$ and with every operator of the word.
<!-- metric_forms.py:176-186; runfile.py:523-528 -->

## Applying $S^{-1}$ and $S^{-1/2}$

The engine applies $S^{-1}$ and $S^{-1/2}$ to vectors, and never forms them
as a matrix product with $H$. The applier is chosen with `metric.applier`:

| Applier | How it applies the metric | Scaling |
|---|---|---|
| `exact` | Dense factorisation of $S$, for $S^{-1}$ and $S^{-1/2}$ | $O(N^3)$ setup, refused above `metric.max_dim` (default 4096) |
| `chebyshev` | Chebyshev expansion of $1/s$ or $s^{-1/2}$ on the spectrum of $S$, truncated when the coefficients fall below `metric.coeff_tol` | Linear in $N$; `metric.order` products with the sparse $S$ |
| `ischo` | Conjugate gradients preconditioned with an incomplete Cholesky factor; provides $S^{-1}$ only | Linear in $N$ per iteration |
<!-- lsquant/core/algebra/overlap_appliers.hpp -->

The `ischo` applier provides only $S^{-1}$, while the Löwdin route needs
$S^{-1/2}$, so the validator refuses that combination.
<!-- runfile.py:534-537 -->

The cost of the linear-scaling appliers is set by the condition number
$\kappa = \mathrm{cond}(S) = s_\mathrm{max}/s_\mathrm{min}$. As a rule of
thumb, a polynomial approximation of $1/s$ on $[s_\mathrm{min}, s_\mathrm{max}]$
needs a degree of order $\sqrt{\kappa}$ times the number of digits required.
Conjugate gradients need a similar number of iterations. Large, diffuse basis
sets have a large $\kappa$. Every application of $\tilde H$ then costs many
products with $S$.

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
<!-- models/nonorthogonal.py:306-329; metric_forms.py:105-122 -->

The three SIESTA velocity modes (`model.params.velocity_mode`) keep different
parts of this expression:

- `bare` takes the derivatives of the lattice phase $e^{ik\cdot R}$ only and
  drops $A$;
- `berry_connection` takes the derivatives of the full phase
  $e^{ik\cdot(R + \tau_j - \tau_i)}$, which includes the orbital centres
  $\tau$, and drops $A$;
- `covariant` adds the dipole term, computed from the radial orbitals in the
  SIESTA `.ion.xml` files.
<!-- models/siesta.py:116-118 -->

The default depends on the kind of SIESTA run. Clean and projected runs use
`berry_connection`, and disordered runs use `covariant`. The modes agree for
the density of states. They can differ for the Hall and spin Hall
conductivities, which depend on interband matrix elements of $V$.
<!-- apps/stage6b_siesta.py:203 (default berry_connection); apps/stage6c_disorder.py:111 (default covariant) -->

LSQuant builds the physical matrix $V_a$ first, and only then transforms it
to the route, like any other operator. It never uses the derivative of the
Löwdin Hamiltonian $S^{-1/2} H S^{-1/2}$ as a velocity.
<!-- models/nonorthogonal.py:315-318 -->

## Spin and spin currents

The spin matrix of a non-orthogonal spinor basis is $\Sigma_z = S\,\sigma_z$,
where $\sigma_z$ is the Pauli matrix acting on the spin index. The overlap of
SIESTA spinor orbitals is the same for both spins, so $S$ commutes with
$\sigma_z$. The operator that acts on coefficients is then $\sigma_z$ itself.
LSQuant accepts only an overlap that commutes with all three Pauli matrices.
The spin current is
$J^z_a = \tfrac12(V_a \sigma_z + \sigma_z V_a)$, the operator used in the
spin Hall word of [Kubo formulas](kubo-formulas.md).
<!-- models/nonorthogonal.py:350-361 -->

## Disorder and conditioning

An on-site Anderson potential $w_i$ in a non-orthogonal basis has the matrix
$V_{ij} = \tfrac12 (w_i + w_j) S_{ij}$. This is the Hermitian average of
assigning the potential to the bra orbital and to the ket orbital. Its effect
on the spectrum is bounded by $\sqrt{\kappa}\, \max_i |w_i|$, not by $W/2$,
and LSQuant pads the spectral window by that bound. For the Bi DZP basis of
[example 12](../../examples/12_disorder_nonorthogonal/README.md),
$\kappa \approx 3550$ and $W = 1$ eV. The window grows to $[-68.5, 100.2]$ eV,
and 256 moments resolve only about 1 eV.
<!-- metric_forms.py:197-218 anderson_forms; disorder_bases.py:116-129 _window -->

In a model with spin–orbit coupling the orbital index $i$ carries spin, and
we choose how $w_i$ treats it with `basis.disorder.kind`. The default,
`scalar`, puts the same energy on both spin components of an orbital, so the
potential is spin-independent and time-reversal symmetric. `magnetic` puts
$+h_i$ on spin up and $-h_i$ on spin down, a random Zeeman field $h_i\sigma_z$.
Because $S$ is the same for both spins, this field flips no spin and conserves
$s_z$. In-plane spins dephase while $s_z$ does not, and for the spin Hall
conductivity it acts as magnetic disorder. `independent` draws one value per
spin-orbital, a scalar and a magnetic part of variance $W^2/24$ each; it
reproduces earlier versions of LSQuant.
<!-- disorder_bases.py:29 spinor_disorder_values; runfile.py:965 default scalar -->

Disordered SIESTA runs check the exact $S^{-1}$ and $S^{-1/2}$ with a probe
vector. They refuse the run if the residual exceeds `metric.tol`, whose
default is $10^{-13}$. For $\kappa \approx 3550$, the rounding errors of an
exact factorisation in double precision can reach
$\kappa\,\varepsilon \approx 8 \times 10^{-13}$. The residual measured in
example 12, $1.4 \times 10^{-13}$, is of rounding size.
<!-- disorder_bases.py:315-319 probe vs tol; runfile.py:297 default tol 1e-13 -->

> **Known limitation: the default `metric.tol` ignores cond(S).** The
> default $10^{-13}$ asks for more than double precision guarantees for an
> overlap with $\kappa$ in the thousands, and the probe can fail on
> rounding noise (`python/lsquant/disorder_bases.py:315–319`,
> `python/lsquant/runfile.py:297`). Example 12 sets `metric.tol` to
> $10^{-10}$, which avoids the failure.

The [bases-and-metrics reference](../reference/bases-and-metrics.md) lists the
keys, and [How we run LSQuant](../how-to/run-lsquant.md#a-basis-with-an-overlap-matrix)
points to the use cases.
