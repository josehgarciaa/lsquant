# Time evolution

Several properties of a disordered or spin–orbit-coupled system are most
naturally expressed in time rather than in energy. Among them are the
distance a wave packet spreads in a given time, the saturation of that
spreading when the states are localised, and the lifetime of a spin
polarisation. They are obtained from the time-evolution operator of a
time-independent Hamiltonian,

$$
U(t) = e^{-iHt/\hbar},
$$

applied to vectors with the same Chebyshev expansion that produces the
density of states ([The kernel polynomial method](kpm.md)). In what follows,
we see how $U(t)$ is expanded, how the mean-square displacement and the spin
dynamics are built from it, and how the results are read.

## The Chebyshev–Bessel expansion of $U(t)$

With the rescaled Hamiltonian $\tilde H = (H - b)/a$, the propagator over one
time step $\Delta t$ factorises into a phase and a function of $\tilde H$. The
Jacobi–Anger identity expands that function in Chebyshev polynomials:

$$
U(\Delta t) = e^{-ib\Delta t/\hbar} \sum_{n=0}^{\infty} c_n\!\left(\frac{a\Delta t}{\hbar}\right) T_n(\tilde H),
\qquad
c_0(z) = J_0(z), \quad c_n(z) = 2(-i)^n J_n(z),
$$

where $J_n$ are Bessel functions of the first kind, $a$ is the half-width and
$b$ the centre of the spectral window.
<!-- source: moments.py:42-69 -->

The Bessel functions decay faster than exponentially once $n$ exceeds their
argument $z = a\Delta t/\hbar$. The series is truncated at the first $n$ where
$|J_n(z)| \le 10^{-15}$, so one step costs about $z$ plus a few tens of
matrix–vector products. A wider spectral window therefore makes each time
step more expensive, as it coarsens the energy resolution of a static
calculation.
<!-- source: spectral_functions.py:186-213 -->

Times are in femtoseconds, with $\hbar = 0.6582119624$ eV fs, so with $H$ in
eV the keys `time.dt` and `time.tmax` are in fs. The time step is `time.dt` if
it is given, and otherwise $t_\mathrm{max}/(N_t - 1)$ for $N_t$ = `time.steps`.
The expansion converges for any $\Delta t$, so the step size sets only the
time resolution of the output; a large step does not make the evolution
unstable.
<!-- source: moments.py:22 HBAR; words/compiler.py:446-456 dt -->

## Moments along a time axis

A time-dependent observable has one Chebyshev index in energy and one time
index. At each time step $t_n = n\,\Delta t$, the evolved vectors are expanded
in energy with $M$ moments and then advanced by one step of $U(\Delta t)$. The
result is a table of moments $\mu_m(t_n)$ of size $M \times N_t$.
<!-- source: lsquant/core/chebyshev/plan.hpp, TimeAxis instruction -->

The reconstruction combines the energy moments at a fixed Fermi energy
$E_F$, with the Jackson kernel and the broadening of the static case:

$$
C(E_F, t_n) = \frac{1}{a} \sum_{m=0}^{M_\mathrm{eff}-1} g_m\, \mu_m(t_n)\, \delta_m(x_F),
\qquad x_F = \frac{E_F - b}{a} .
$$

The formula `time_correlation` writes one curve per entry of
`reconstruction.fermi_energies`. The broadening sets the energy window of the
states that contribute to each curve, so a finer broadening follows a
narrower group of states around $E_F$.
<!-- source: reconstruct.py:665-728 time_correlation_curve -->

## Mean-square displacement

The spreading of a wave packet at energy $E$ is measured by the mean-square
displacement,

$$
\Delta X^2(E, t) = \frac{\mathrm{Tr}\big[[X, U(t)]^\dagger\, \delta(E - H)\, [X, U(t)]\big]}{\mathrm{Tr}\,\delta(E - H)},
$$

where $X$ is the position operator along $x$. Its time dependence classifies
transport. Ballistic motion gives $\Delta X^2 \propto t^2$, diffusion gives
$\Delta X^2 = 2Dt$ in one direction with diffusion coefficient $D$, and
Anderson localisation makes $\Delta X^2$ saturate at a value set by the
localisation length. In the diffusive regime, the plateau of
$D(E, t) = \Delta X^2/2t$ gives the conductivity through the Einstein
relation $\sigma(E) = e^2 \rho(E) D$, with $\rho(E)$ the density of states;
[example 05](../../examples/05_charge_transport_square/README.md) obtains
$\sigma_{xx}$ of a disordered square lattice in this way.
<!-- source: examples/05_charge_transport_square/README.md, "The physics" (D = ΔX²/2t, plateau D_max, σ(E) = e² ρ(E) D_max) and Step 3 -->

The position operator is ill-defined in a periodic supercell, but its
commutator with $U(t)$ is well defined. $X$ itself is never built; instead,
$[X, U(t)]$ is updated step by step with

$$
[X, U(t + \Delta t)] = [X, U(\Delta t)]\, U(t) + U(\Delta t)\, [X, U(t)],
$$

and $[X, U(\Delta t)]$ is obtained from a second Chebyshev recursion driven
by the velocity operator $v_x \propto i[H, X]$. The commutator `[X,U(t)]` in
the observable therefore needs the velocity matrix `VX` of the model, not a
position matrix. Hopping across the supercell boundary enters through this
velocity matrix, so the supercell keeps every bond.
<!-- source: words/catalogue.py:381-391 ([X,U(t)] -> VX); lsquant/core/chebyshev/plan.hpp:82-100 ApplyDerivedSeries -->

The curve written to file is the numerator only,
$\mathrm{Tr}\big[[X, U]^\dagger \delta(E-H) [X, U]\big]/N$. This is the
product $\rho(E_F)\,\Delta X^2(E_F, t)$ per orbital, in Å²/eV. To obtain
$\Delta X^2$ in Å², we divide by the density of states per orbital computed
with the same random vectors, as
[example 15](../../examples/15_spatial_operators/README.md) does.
[Example 02](../../examples/02_localization_1d_chain/README.md) uses the
saturation of the mean-square displacement to extract the localisation length
of a disordered chain.

## Spin dynamics

Spin precession and relaxation follow from the time evolution of a spin
operator $s$ in the Heisenberg picture, $s(t) = U^\dagger(t)\, s\, U(t)$. We
prepare a spin-polarised ensemble at the Fermi energy with a projector $P$,
for example onto spin along $+x$, and follow its average spin [1]:

$$
S(E_F, t) = \frac{\mathrm{Tr}\big[P\, U^\dagger(t)\, s\, U(t)\, \delta(E_F - H)\, P\big]}
{\mathrm{Tr}\big[P\, \delta(E_F - H)\big]} .
$$

The observable `Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]` gives the numerator
for spin along $x$. The prefix `P:` marks `PX` as the projector, and `SX` is
the spin operator that is followed. Because the projector acts on both sides,
the trace counts only states that start with the chosen spin.

The numerator of $S(E_F, t)$ is written to file, and the formula
`time_correlation` gives its energy-resolved curve at each Fermi energy. The
default formula of this observable, `row_m0`, gives instead the spin
integrated over all energies, which is the $m = 0$ moment of each time step;
it writes no curve file, so its values are read from the `.chebmomTD` moment
file (see the [run-file reference](../reference/run-file.md#formulas-by-word-shape)).
[Example 04](../../examples/04_spin_precession_magnetic_chain/README.md)
normalises by the value at $t = 0$ and recovers Larmor precession in a
Zeeman field. [Example 06](../../examples/06_spin_relaxation_square/README.md)
fits the decay of the spin and compares the relaxation times with the
Elliott–Yafet and Dyakonov–Perel limits.
<!-- source: run-file.md "Formulas by word shape" (projected_evolution: allowed row_m0, time_correlation; default row_m0; row_m0 writes no curve file); runfile.py:1277-1278 row_m0 returns no curve; m = 0 moment and .chebmomTD from the previous version of this page -->

The same expansion evaluates a general time correlation
$\mathrm{Tr}[A\, \delta(E - H)\, U^\dagger(t)\, B\, U(t)]$ for any two static
operators $A$ and $B$.

## Cost

Each random vector carries a few vectors through time: two for the spin and
time-correlation observables, and three coupled vectors for the mean-square
displacement. Per time step, each evolved vector costs about $a\Delta t/\hbar$
matrix–vector products for the Bessel series, and the energy expansion adds
$M$ products. A run with $N_t$ steps and $R$ random vectors therefore costs
of order $R\, N_t\, (M + k\, a\Delta t/\hbar)$ products, with $k$ the number
of evolved vectors. Memory stays at a few vectors of length $N$ per random
vector, so real-space supercells of millions of orbitals can be followed in
time.

[Running LSQuant](../how-to/run-lsquant.md) shows the `time` block of a run
file, and the [run-file reference](../reference/run-file.md#time) lists its
keys.

## References

1. A. W. Cummings, J. H. Garcia, J. Fabian and S. Roche, Giant spin lifetime
   anisotropy in graphene induced by proximity effects, Phys. Rev. Lett.
   **119**, 206601 (2017).
