# The Hall conductivity of the Haldane model

In this tutorial we compute the Hall conductivity of a Chern insulator and see
it locked to $e^2/h$ inside the gap, at zero magnetic field. The Haldane model
puts a real nearest-neighbour hopping $t_1$ and a complex next-nearest-neighbour
hopping $t_2 e^{\pm i\phi}$ on a honeycomb lattice. At $t_1 = -1$ eV,
$t_2 = 0.15$ eV and $\phi = \pi/2$ the lower band has Chern number $C = 1$, and
with the Fermi level in the gap

$$
\sigma_{xy} = C\,\frac{e^2}{h}.
$$

LSQuant evaluates the Hall conductivity with the Kubo–Bastin formula, written
as the observable word `Tr[VX delta(E-H) VY dG-(E)/dE]`, where `VX` and `VY`
are the velocity operators and `dG-(E)/dE` is the energy derivative of the
advanced Green's function. In the same run it computes the longitudinal
conductivity σ_xx with the Kubo–Greenwood word, which vanishes in the gap.

We assume LSQuant is built as in [Getting started](getting-started.md). We run
the example, inspect the output and plot it:

```bash
export OMP_NUM_THREADS=1
python3 -m lsquant run examples/03_haldane_hall/run_haldane.json --out /tmp/haldane
python3 -m lsquant inspect /tmp/haldane
python3 examples/03_haldane_hall/plot.py /tmp/haldane
```

The run file asks for 128 moments on the window [−4.5, 4.5] eV and for four
reconstructions: the full Kubo–Bastin curve (`bastin`), its Fermi-sea and
Fermi-surface parts (`sea`, `surf`) and σ_xx (`greenwood`).

The full lesson, with the plot read against the Chern number and the reason
the plateau survives where no states sit at the Fermi level, is in
[example 03](../../examples/03_haldane_hall/README.md).
[Example 13](../../examples/13_bastin_sea_surface/README.md) checks the
sea and surface split on the same model, and
[Kubo formulas](../explanation/kubo-formulas.md) explains the two formulas.
