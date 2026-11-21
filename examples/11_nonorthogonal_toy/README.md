# Example 11: a non-orthogonal basis, Löwdin and generalised routes

## The physical question

Local-orbital DFT codes such as SIESTA expand the Kohn–Sham states in atomic
orbitals that overlap. The Hamiltonian $H$ then comes with an overlap matrix
$S \neq 1$, and the energies solve the generalised eigenvalue problem

$$ H\,c = E\,S\,c . $$

The density of states counts these generalised eigenvalues,
$\rho(E) = \frac1N\sum_i \delta(E - E_i)$. To expand it in Chebyshev
polynomials, KPM needs a single matrix with the same spectrum that is Hermitian,
or at least diagonalisable. LSQuant offers two such matrices:

- the **Löwdin route** (`metric.route = "lowdin"`) works with
  $S^{-1/2} H S^{-1/2}$, a Hermitian matrix in an orthonormalised basis;
- the **generalised route** (`metric.route = "generalised"`) works with
  $S^{-1} H$, which is not Hermitian but is similar to the Löwdin matrix.

Both matrices have exactly the eigenvalues $E_i$. Both must therefore give the
same moments $\mu_m = (2-\delta_{m0})\,\frac1N\,\mathrm{Tr}\,T_m(\tilde H)$.
We test this on a toy model whose generalised spectrum we know in closed form.
We also compare the two ways of applying $S^{-1/2}$ or $S^{-1}$. The `exact`
applier uses a dense factorisation and keeps an $N\times N$ factor of $S$. It
is therefore meant for small systems, and `metric.max_dim`, 4096 by default,
sets the limit. The `chebyshev` applier expands the inverse as a polynomial in
$S$. It only multiplies by the sparse $S$, so it scales linearly with the
system size.

## The model

[`make_nonorth_chain.py`](make_nonorth_chain.py) starts from the periodic chain
of $N = 128$ sites with hopping $t = -1$ eV, whose eigenvalues are
$E_k = 2t\cos(2\pi k/N)$. It then changes basis with a local, non-unitary
transformation $R = 1 + \alpha K$, where $K$ is random on the
nearest-neighbour graph and $\alpha = 0.55$:

$$ S = R^\dagger R, \qquad H' = R^\dagger H R . $$
<!-- source: lsquant.models.nonorthogonal.make_local_transform -->

The pair $(H', S)$ is a genuine non-orthogonal tight-binding model. Its overlap
is banded with a condition number of about 6, and its generalised eigenvalues
are still exactly $2t\cos(2\pi k/N)$. We generate it with:

```bash
python3 examples/11_nonorthogonal_toy/make_nonorth_chain.py 128
```

The script prints the spectrum of the overlap and the deviation of the
generalised eigenvalues from the closed form:

```text
N = 128: spec(S') = [0.3453, 2.0572], cond(S') = 5.958, nnz(S') = 640
max |eig(H', S') - 2t cos(2 pi k/N)| = 3.21e-15
```

It writes the orthogonal chain to `operators/chain.HAM.CSR` and the
non-orthogonal model to `operators/nonorth_chain.{HAM,S,VX}.CSR`. A `csr` model
with `metric.policy = "nonorthogonal"` reads `<label>.S.CSR` next to
`<label>.HAM.CSR`.

## The runs

We evaluate the same word, `Tr[delta(E-H)]`, four times. Every run uses 256
moments, an exact trace and the window $[-2.2, 2.2]$ eV. The orthogonal chain
gives the reference:

```bash
python3 -m lsquant run examples/11_nonorthogonal_toy/run_orthogonal.json --out /tmp/ex11/nonorth_orthogonal
```

Next, we run the non-orthogonal model on the Löwdin route with the exact
$S^{-1/2}$. The metric section of [`run_lowdin.json`](run_lowdin.json) reads:

```json
"metric": {"policy": "nonorthogonal", "route": "lowdin", "applier": "exact"}
```

The run is started with:

```bash
python3 -m lsquant run examples/11_nonorthogonal_toy/run_lowdin.json --out /tmp/ex11/nonorth_lowdin
```

We then switch to the generalised route with the exact $S^{-1}$, as set in
[`run_generalised.json`](run_generalised.json):

```json
"metric": {"policy": "nonorthogonal", "route": "generalised", "applier": "exact"}
```

The run is started with:

```bash
python3 -m lsquant run examples/11_nonorthogonal_toy/run_generalised.json --out /tmp/ex11/nonorth_generalised
```

Finally, we return to the Löwdin route and expand $S^{-1/2}$ as a Chebyshev
polynomial of $S$. The expansion is truncated when its coefficients drop below
`coeff_tol`. The metric section of
[`run_lowdin_chebyshev.json`](run_lowdin_chebyshev.json) reads:

```json
"metric": {"policy": "nonorthogonal", "route": "lowdin", "applier": "chebyshev",
           "order": 200, "coeff_tol": 1e-14}
```

The run is started with:

```bash
python3 -m lsquant run examples/11_nonorthogonal_toy/run_lowdin_chebyshev.json --out /tmp/ex11/nonorth_lowdin_chebyshev
```

Each run writes the moments as `SpectralOp1*.chebmom1D` and as the binary
`.chebmom`, which a non-orthogonal metric writes by default. It also writes the
DOS at 30 meV as `meanSpectralOp1*JACKSON.dat`, and `run.resolved.json`.

## The result

The four DOS curves are the same chain DOS $1/\pi\sqrt{4-E^2}$, resolved into
the discrete levels of the 128-site ring. The overlap changes the basis, not
the physics. The metric costs time, and
[the physical question](#the-physical-question) describes the cost of each
applier. For the Chebyshev applier, `order`, `coeff_tol` and the condition
number of $S$ set the number of products.

## Checks

[`check_nonorth.py`](check_nonorth.py) solves $H'c = ESc$ with
`scipy.linalg.eigh`, builds the exact moments from the eigenvalues and compares
them with each run. We run it with:

```bash
python3 examples/11_nonorthogonal_toy/check_nonorth.py /tmp/ex11
```

Output on the reference run:

```text
overlap S: N = 128, spectrum [0.3453, 2.0572], cond(S) = 5.958
generalised eigenvalues vs 2t cos(2 pi k/N): max deviation 2.44e-15
orthogonal         max |mu_m - mu_m(exact)| = 4.216e-14  (tol 1e-10)  PASS
lowdin             max |mu_m - mu_m(exact)| = 7.178e-14  (tol 1e-10)  PASS
generalised        max |mu_m - mu_m(exact)| = 4.199e-14  (tol 1e-10)  PASS
lowdin_chebyshev   max |mu_m - mu_m(exact)| = 1.303e-13  (tol 1e-08)  PASS
CHECK PASS
```

Both routes reproduce the exact moments of the generalised problem to
$10^{-13}$, with either applier.

## Runtime

On a shared 96-core node with one OpenMP thread, the orthogonal run took 2 s,
the exact Löwdin run 4 s, the exact generalised run 2 s, the Chebyshev Löwdin
run 9 s and the check 1 s. On such a busy node a small system runs much faster
with `OMP_NUM_THREADS=1`: the orthogonal run took 146 s with four threads and
6 s with one.
