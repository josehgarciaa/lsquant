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

LSQuant replaces the diagonalisation of the chain by an expansion. Roughly
speaking, it expands $\rho(E)$ in Chebyshev polynomials of the Hamiltonian and
computes the expansion coefficients, the **Chebyshev moments**, from repeated
products of $H$ with vectors. A second, cheap step, the **reconstruction**,
turns the moments into a curve. A JSON file, the **run file**, describes the
whole calculation. Every LSQuant calculation follows these three steps: we
write a run file, compute the moments, and reconstruct the curve.

## Building LSQuant

We assume that the compiler, the libraries and the Python packages are
installed as described in [Install LSQuant](../how-to/install.md). From the
top of the repository, we build LSQuant with the `golden` preset, which gives
the same numbers on every machine:

```bash
cmake --preset golden
cmake --build --preset golden
```

The presets need CMake 3.21 or later; [Install LSQuant](../how-to/install.md)
gives the commands for an older CMake. Python finds the compiled module and the
LSQuant package through `PYTHONPATH`, and we check that the import works:

```bash
export PYTHONPATH="$PWD/build/golden/lsquant_core:$PWD/python"
python3 -c "import lsquant._core; print('LSQuant core loaded')"
```

## The run file

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
`examples/01_dos_1d_chain/make_chain.py`. The `basis` section stores the
Hamiltonian as a sparse matrix in real space.

The quantity we compute is the trace $\mathrm{Tr}\,\delta(E - H)$, which is $N$
times the DOS per site for a ring of $N$ sites. In the run file it is written
as `Tr[delta(E-H)]`. LSQuant calls such a trace an **observable word**;
conductivities and spin dynamics are other words of the same kind.

The `expansion` section asks for 512 moments. The Chebyshev polynomials are
defined on $[-1, 1]$, so LSQuant first maps the energy interval $[-2, 2]$ eV,
the **spectral window**, onto $[-1, 1]$. The window must contain the whole
spectrum. Here it ends exactly at the band edges, which is safe, because
$|T_m(\pm 1)| = 1$ for every polynomial. The trace is exact: it runs over all
sites.

The `reconstruction` section turns the moments into a curve with the Jackson
kernel and an energy resolution, the **broadening**, of 20 meV. The
[run-file reference](../reference/run-file.md) lists every key.

## Computing the moments

We run the calculation on one thread. On a shared machine, more threads can
make this small problem much slower (see [Run on a cluster](../how-to/run-on-a-cluster.md)):

```bash
export OMP_NUM_THREADS=1
python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/chain_dos
```

The run took about 2.5 minutes on one core of a shared node. The command prints
the path of each file it writes. The output directory now holds the moments of
each ring, in files ending in `.chebmom1D`, and `run.resolved.json`, the run
file with every default filled in. The DOS of the three rings at 20 meV is in
`eta20/`, in files ending in `JACKSON.dat`. Each curve file has two columns:
the energy in eV and $N\rho(E)$ in 1/eV, so we divide by the number of sites
$N$ to obtain the DOS per site.

For the figures of this tutorial, the chain generator also writes curves of the
256-site ring at 80, 24 and 10 meV (in `eta80/`, `eta24/` and `eta10/`) and
at 120 meV (at the top level). A run file with any other kind of model writes
one curve per word and formula, at the broadening we set.

## Inspecting the result

We then ask LSQuant what the directory contains:

```bash
python3 -m lsquant inspect /tmp/chain_dos
```

The command prints JSON: the resolved run file and, for each moment file, the
word, the number of moments, the spectral window and the system size. We can
use it whenever we need to know what produced a set of files.

## Plotting the density of states

Finally, we plot the result:

```bash
python3 examples/01_dos_1d_chain/plot.py /tmp/chain_dos
```

The script writes three figures to the output directory. The reconstructed
curve covers the inner 95 % of the window, $|E| \le 1.9$ eV, and the van Hove
singularities show up as the steep rise of the DOS towards both ends of it.

In the upper panel of `fig_dos.png` we see that the DOS per site of the 1024-
and 4096-site rings follows the infinite-chain result. For the 256-site ring,
we see the curve break into peaks. A broadening $\eta$ uses about $\pi a/\eta$
moments, where $a = 2$ eV is the half-width of the window, so 20 meV uses
about 315 moments. A ring of $N$ sites has $N$ discrete levels, and the
moments resolve them as soon as more than $N$ of them contribute. The lower
panel shows the 256-site ring at 80, 24 and 10 meV, on both sides of this
threshold.

## Changing the resolution from the same moments

The moments contain all the information about the rings. We can therefore
reuse them to look at the DOS with a coarser resolution of 60 meV:

```bash
python3 -m lsquant reconstruct /tmp/chain_dos --formula dos --broadening-meV 60 --out /tmp/chain_dos_60
```

This takes a few seconds, because the reconstruction works from the stored
moments and skips the products with the Hamiltonian. It writes one curve per
moment file to `/tmp/chain_dos_60`. At 60 meV the reconstruction uses about
105 moments, fewer than the 256 levels of the smallest ring, and all three
rings give the smooth bulk DOS.

As an exercise, reconstruct the 256-site ring at broadenings between 60 and
10 meV and find where its discrete levels appear.

## Next steps

- [Example 01](../../examples/01_dos_1d_chain/README.md) studies the same chain
  in more depth, including what the moments of a finite ring look like.
- [Example 09](../../examples/09_inspect_reconstruct/README.md) reconstructs one
  set of moments at several broadenings and checks the result against the
  exact DOS.
- [Running LSQuant](../how-to/run-lsquant.md) goes through the choices of a
  calculation: model, basis, word and number of moments.
- [The kernel polynomial method](../explanation/kpm.md) explains the expansion
  and the kernel.
