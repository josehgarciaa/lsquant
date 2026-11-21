# Example 12: Anderson disorder in a non-orthogonal SIESTA basis

## The physical question

In an orthogonal tight-binding model, on-site disorder adds a diagonal matrix to
the Hamiltonian: $H \to H + D$, with $D = \mathrm{diag}(w_i)$ and $w_i$ drawn
uniformly from $[-W/2, W/2]$. A non-orthogonal basis needs more than this
diagonal. An orbital $|i\rangle$ overlaps its neighbours, so a potential that
shifts orbital $i$ by $w_i$ has off-diagonal matrix elements as well. LSQuant
uses the Hermitian form

$$ V = \tfrac12\,(D\,S + S\,D), \qquad H \to H + V, \qquad S \text{ unchanged}, $$

which reduces to $D$ when $S = 1$. It then solves the generalised problem
$(H + V)\,c = E\,S\,c$.

We compute the DOS of the Bi(111) bilayer of example 7 with $W = 1$ eV of
disorder, for two disorder realisations. We also check that two different
layouts of the same calculation agree:

- in the `real` basis, the supercell $H$, $S$ and $V$ are sparse CSR matrices,
  and $S^{-1/2}$ or $S^{-1}$ is applied exactly or by a Chebyshev expansion.
  This is the general route, and the one to use for defects or any
  non-periodic structure;
- in the `bloch` basis, the clean blocks $H(\mathbf k)$ and $S(\mathbf k)$
  are kept, and the potential acts through a Fourier transform $B$ as
  $V = \frac12\{B D B^\dagger, S(\mathbf k)\}$. The metric is applied exactly
  at each $\mathbf k$, and no supercell matrix is stored.
<!-- source: real layout = RealSpaceMetricBasis, bloch layout = DisorderedBlochMetricBasis, both in python/lsquant/disorder_bases.py -->

## The model

The model is the committed SIESTA calculation of example 7, with label `bi` in
`examples/07_spin_hall_bi_bilayer/siesta`. It has two Bi atoms with 18 DZP
orbitals each and spin–orbit coupling, so a cell holds 72 spin-orbitals. On a
$4\times4$ supercell the Hamiltonian has dimension 1152. The run file sets the
disorder in the `basis` section:

```json
"basis": {"name": "real", "supercell": [4, 4, 1], "disorder": {"W": 1.0, "seed": 7, "kind": "independent"}}
```

With `"kind": "independent"` each spin-orbital coefficient gets its own value,
the layout of earlier versions of LSQuant; we keep it so that the numbers below
are reproduced exactly. Without the key, the run uses the default `"scalar"`:
one value per orbital on both spin components, a spin-independent potential.
`"magnetic"` gives a random Zeeman field along $z$ instead. The
[run-file reference](../../docs/reference/run-file.md#basis) gives the exact
draws.

The independent draw is a scalar part plus a random Zeeman field along $z$ on
each orbital: $(w_\uparrow - w_\downarrow)/2$ has variance $W^2/24$. The
overlap is spin-scalar, so this field flips no spin; it breaks time-reversal
symmetry and conserves $s_z$. None of this changes the checks below, which
compare layouts and seeds and test the sum rule.
<!-- source: python/lsquant/disorder_bases.py:29 (spinor_disorder_values); python/lsquant/runfile.py:501-510 -->

The metric uses the Löwdin route with the exact applier. The applier is checked
against a tolerance `metric.tol`, whose default is $10^{-13}$, and we raise it
to $10^{-10}$. For this DZP basis $\mathrm{cond}(S) \approx 3550$, so in double
precision the rounding errors of an exact factorisation can reach
$\mathrm{cond}(S)\,\varepsilon \approx 8\times10^{-13}$. The residual we
measure, $1.4\times10^{-13}$, is of that rounding size. The default
$10^{-13}$ therefore asks for more than double precision guarantees for this
overlap.
<!-- source: python/lsquant/disorder_bases.py:315-319 (metric probe vs tol) -->

```json
"metric": {"policy": "nonorthogonal", "route": "lowdin", "applier": "exact", "tol": 1e-10}
```

## The runs

The trace is stochastic, with 4 random-phase vectors and 256 moments. The
window comes from the spectrum (`"eigen"`). We make three runs: seed 7 in both
layouts, and seed 8 in the real-space layout:

```bash
python3 -m lsquant run examples/12_disorder_nonorthogonal/run_real_seed7.json  --out /tmp/ex12/disorder_real_seed7
python3 -m lsquant run examples/12_disorder_nonorthogonal/run_real_seed8.json  --out /tmp/ex12/disorder_real_seed8
python3 -m lsquant run examples/12_disorder_nonorthogonal/run_bloch_seed7.json --out /tmp/ex12/disorder_bloch_seed7
```

Each run directory holds `dos.dat`, the DOS in states/eV/cell, and
`fermi_values.json`, the DOS at $E_F = 0$. It also holds the moments,
`disorder.npy` with the disorder vector actually used, and `basis.json`. The
last file records the disorder seed and its SHA-256 hash, the potential form
`(D S+S D)/2`, the window, the overlap bounds and the metric residuals.

## The result

The two seeds give DOS curves that differ by at most 0.14 states/eV/cell, for
a DOS whose maximum is 5.05. They are two realisations of the same disordered
system. At $E_F = 0$ the seed-7 DOS is 1.12 states/eV/cell.

The window needs a comment. For a non-orthogonal basis the disorder can shift
the spectrum of $S^{-1/2}(H+V)S^{-1/2}$ by up to
$\sqrt{\mathrm{cond}(S)}\,\max|w_i|$, not only by $W/2$. LSQuant pads the
window by that bound so that the expansion is guaranteed to converge. With
$\mathrm{cond}(S) \approx 3550$ the padding is $\pm 30$ eV, and the window
becomes $[-68.5, 100.2]$ eV. With 256 moments the resolution is then about
$\pi \times 84\ \mathrm{eV} / 256 \approx 1$ eV, not the requested 50 meV. A
production calculation needs proportionally more moments in
`expansion.moments`, or a basis with a better-conditioned overlap.

## Checks

[`check_disorder.py`](check_disorder.py) compares the three runs. We run it
with:

```bash
python3 examples/12_disorder_nonorthogonal/check_disorder.py /tmp/ex12
```

Output on the reference run:

```text
dimension 1152, window [-68.48572453353346, 100.16258869984046], overlap bounds [0.0009599617497710542, 3.409242598313619]
disorder sha256: seed7 real 2d42a49aab5f, seed7 bloch 2d42a49aab5f, seed8 real 4eb408e220b9
[1] real vs bloch, seed 7: same disorder vector = True; max |d mu| = 1.27e-12  PASS
[2] seed 7 vs seed 8: max |d DOS| = 0.136 states/eV/cell (DOS max 5.05); differ = True
[3] seed7 real: integral DOS dE = 72.000 (expected 72)  PASS
[3] seed7 bloch: integral DOS dE = 72.000 (expected 72)  PASS
[3] seed8 real: integral DOS dE = 72.000 (expected 72)  PASS
CHECK PASS
```

1. The real-space layout and the Bloch layout draw the same disorder vector
   from the same seed. They also use the same random-phase vectors, mapped by
   $B^\dagger$, so their moments agree to $10^{-12}$.
2. Seeds 7 and 8 draw different potentials, with different hashes.
3. The DOS per cell integrates to the 72 spin-orbitals per cell, as the sum
   rule requires.

## Runtime

On a shared 96-core node with one OpenMP thread, the two real-space runs took
64 s and 104 s, the Bloch run 20 s and the check 7 s.
