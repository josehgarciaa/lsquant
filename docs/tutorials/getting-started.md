# Getting started: the density of states of a chain

As a first example, we compute the density of states (DOS) of a
one-dimensional tight-binding chain and compare it with the exact result. A
chain with one orbital per site and nearest-neighbour hopping $t$ has the
band $E(k) = 2t\cos k$. Its DOS per site is

$$
\rho(E) = \frac{1}{\pi\sqrt{4t^2 - E^2}}, \qquad |E| < 2|t|,
$$

which diverges at the band edges $E = \pm 2|t|$. These divergences are the
van Hove singularities of a one-dimensional band.

LSQuant replaces the diagonalisation of the chain by an expansion. It expands
$\rho(E)$ in Chebyshev polynomials of the Hamiltonian and computes the
expansion coefficients, the **Chebyshev moments**, from repeated products of
$H$ with vectors. A second, cheap step, the **reconstruction**, turns the
moments into a curve. A JSON file, the **run file**, describes the whole
calculation, as an `.fdf` file describes a SIESTA run. Every LSQuant
calculation follows these three steps: we write a run file, compute the
moments, and reconstruct the curve.

## Build

We assume that the compiler, the libraries and the Python packages are
installed as described in [Install LSQuant](../how-to/install.md). From the
top of the repository, we build LSQuant with the `golden` preset, which gives
the same numbers on every machine:

```bash
cmake --preset golden
cmake --build --preset golden
```

The presets need CMake 3.21 or later; [Install LSQuant](../how-to/install.md)
gives the commands for an older CMake. We make the compiled module and the
LSQuant package visible to Python through `PYTHONPATH` and check that the
import works:

```bash
export PYTHONPATH="$PWD/build/golden/lsquant_core:$PWD/python"
python3 -c "import lsquant._core; print('LSQuant core loaded')"
```

## Run file

For the chain we use `examples/01_dos_1d_chain/run_dos.json`, whose main
sections are:

```json
{
  "model": { "kind": "generator", "name": "make_chain", "tutorial": 1,
             "params": { "sizes": [256, 1024, 4096], "hopping_t": -1.0 } },
  "basis": { "name": "real" },
  "word": "Tr[delta(E-H)]",
  "expansion": { "moments": 512, "window": [-2.0, 2.0], "trace": "exact" },
  "reconstruction": { "kernel": "jackson", "broadening_meV": 20.0, "formulas": ["dos"] }
}
```

The `model` section builds three periodic chains, or rings, of 256, 1024 and
4096 sites with $t = -1$ eV, using the generator
`examples/01_dos_1d_chain/make_chain.py`; the key `"tutorial": 1` selects
the calculation this generator runs for the present tutorial, which is also
its default. The ring length plays the role of the supercell size: it sets the
finite-size error, as the k-mesh does in a DFT calculation. The `basis`
section stores the Hamiltonian as a sparse matrix in real space.
<!-- source: tutorial key: python/lsquant/runfile.py:91 (make_chain -> 1), 1292-1300 (tutorial 1 runs apps/tutorial01); docs/reference/run-file.md#model -->

The quantity we compute, the **observable**, is the trace
$\mathrm{Tr}\,\delta(E - H)$, which is $N$ times the DOS per site for a ring
of $N$ sites. The run file gives it in the key `word`, as `Tr[delta(E-H)]`.
The conductivities and the spin dynamics are traces of the same kind, with
velocity, spin-current or spin operators between the functions of $H$.
<!-- source: quantity table of docs/how-to/run-lsquant.md; python/lsquant/words/catalogue.py:45-123 -->

The `expansion` section sets 512 moments, the number that fixes the finest
energy resolution the curve can reach. The Chebyshev polynomials are defined
on $[-1, 1]$, so LSQuant first maps the energy interval $[-2, 2]$ eV, the
**spectral window**, onto $[-1, 1]$. The window must contain the whole
spectrum. Here it ends exactly at the band edges, which is safe, because
$|T_m(\pm 1)| = 1$ for every polynomial. With $a = 2$ eV, the half-width of
the window, $M = 512$ moments resolve about $\pi a/M \approx 12$ meV. The
trace is exact: the recursion starts from every site, so the result carries
no statistical noise from random vectors.
<!-- source: pi*2/512 = 0.0123 eV; reconstruct.py:89-92 (max_mom = ceil(pi/eta), eta = 2*broadening/width) -->

The `reconstruction` section turns the moments into a curve with the Jackson
kernel and an energy resolution, the **broadening** $\eta$, of 20 meV. The
broadening is comparable with the smearing of a DFT DOS, and 20 meV lies
above the 12 meV limit set by the 512 moments. The
[run-file reference](../reference/run-file.md) lists every key.

## Moments

We run the calculation on one thread. On a shared machine, more threads can
make this small problem much slower (see [Run on a cluster](../how-to/run-on-a-cluster.md)):

```bash
export OMP_NUM_THREADS=1
python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/chain_dos
```

The run took about 2.5 minutes on one core of a shared node. It takes
minutes for so small a system because the exact trace costs $N \times M$
matrix–vector products, one recursion of $M$ steps from each of the $N$ sites
(see [Exact and random-phase traces](../explanation/kpm.md#exact-and-random-phase-traces)).
Large systems use $R$ random vectors instead, which cost $R \times M$
products. The command prints the path of each file it writes. The output
directory now holds the moments of each ring, in files ending in
`.chebmom1D`, and `run.resolved.json`, the run file with every default filled
in.

The DOS of the three rings at 20 meV is in `eta20/`, in files ending in
`JACKSON.dat`. Each curve file has two columns: the energy in eV and
$N\rho(E)$ in 1/eV, so we divide by the number of sites $N$ to obtain the DOS
per site.

The chain generator also writes extra curves for the figures of this
tutorial: the 256-site ring at 80, 24 and 10 meV, in `eta80/`, `eta24/` and
`eta10/`, and at 120 meV at the top level. These extra curves belong to this
example only. With any other kind of model, a run writes one curve for each
observable and formula, at the broadening of the run file.

## Output directory

We then ask LSQuant what the directory contains:

```bash
python3 -m lsquant inspect /tmp/chain_dos
```

The command prints JSON: the resolved run file and, for each moment file, the
observable, the number of moments, the spectral window and the system size.
We can use it whenever we need to know what produced a set of files.

## Density of states

Finally, we plot the result:

```bash
python3 examples/01_dos_1d_chain/plot.py /tmp/chain_dos
```

The script writes three figures to the output directory; Fig. 1 is
`fig_dos.png`. The reconstructed curve covers the inner 95 % of the window,
$|E| \le 1.9$ eV, and the van Hove singularities show up as the steep rise of
the DOS towards both ends of it.

![Density of states per site of rings of 256, 1024 and 4096 sites](../../examples/01_dos_1d_chain/fig_dos.png)

*Figure 1. DOS per site $\rho(E)$ of periodic chains with $t = -1$ eV,
reconstructed with the Jackson kernel. (a) Broadening $\eta = 20$ meV for
rings of $N = 256$ (solid black), 1024 (dashed blue) and 4096 sites
(dash-dotted orange). (b) The 256-site ring at $\eta = 80$ meV (solid
black), 24 meV (dashed blue) and 10 meV (dash-dotted orange); the legend
gives $\pi a/\eta$ for each broadening.*
<!-- source: figure as committed and as drawn by examples/01_dos_1d_chain/plot.py:138-169 -->

Fig. 1(a) shows that the DOS per site of the 1024- and 4096-site rings follows
the infinite-chain result, while the curve of the 256-site ring breaks into
peaks, the discrete levels of the finite ring. A ring of $N$ sites has $N$
discrete levels, and the curve resolves them as soon as more than $N$ moments
contribute. At 20 meV about $\pi a/\eta \approx 315$ moments contribute:
more than the 256 levels of the smallest ring, and fewer than the levels of
the two larger rings.
Fig. 1(b) shows the 256-site ring at 80, 24 and 10 meV, on both sides of this
threshold.

> **Rule of thumb.** A broadening $\eta$ keeps about $\pi a/\eta$ moments,
> at most the $M$ computed, where $a$ is the half-width of the window. Keep
> this number below the number of levels of the system, $N$ for a ring of $N$
> sites, or the curve shows the discrete levels instead of the bulk DOS.
<!-- source: reconstruct.py:89-92 (max_mom = ceil(pi/eta), capped at the stored moments) -->

## Resolution from the same moments

The moments contain all the information needed for the DOS of the rings. We
can therefore reuse them to look at the DOS with a coarser resolution of
60 meV:

```bash
python3 -m lsquant reconstruct /tmp/chain_dos --formula dos --broadening-meV 60 --out /tmp/chain_dos_60
```

This takes a few seconds, against 2.5 minutes for the moments, because the
reconstruction works from the stored moments and skips the products with the
Hamiltonian. It writes one curve per moment file to `/tmp/chain_dos_60`. At
60 meV the reconstruction uses about 105 moments, fewer than the 256 levels of
the smallest ring, and all three rings give the smooth bulk DOS.

As an exercise, reconstruct the 256-site ring at broadenings between 60 and
10 meV and find where its discrete levels appear.

## Next steps

- [Example 01](../../examples/01_dos_1d_chain/README.md) studies the same chain
  in more depth, including what the moments of a finite ring look like.
- [Example 09](../../examples/09_inspect_reconstruct/README.md) reconstructs one
  set of moments at several broadenings and checks the result against the
  exact DOS.
- [Running LSQuant](../how-to/run-lsquant.md) goes through the choices of a
  calculation: model, basis, observable and number of moments.
- [The kernel polynomial method](../explanation/kpm.md) explains the expansion
  and the kernel.
