# The Hall conductivity of the Haldane model

The Haldane model [1] is the simplest Chern insulator: a honeycomb lattice
that carries a quantised Hall current at zero magnetic field. In this tutorial
we compute its Hall conductivity with LSQuant and see it locked to $e^2/h$
inside the gap.

## Model

The model puts a real nearest-neighbour hopping $t_1$ and a complex
next-nearest-neighbour hopping $t_2 e^{\pm i\phi}$ on a honeycomb lattice. At
$t_1 = -1$ eV, $t_2 = 0.15$ eV and $\phi = \pi/2$ the lower band has Chern
number $C = 1$, and with the Fermi level in the gap

$$
\sigma_{xy} = C\,\frac{e^2}{h}.
$$

LSQuant computes the Hall conductivity with the Kubo–Bastin formula. In the
run file this observable is the trace `Tr[VX delta(E-H) VY dG-(E)/dE]`, where
`VX` and `VY` are the velocity operators along x and y, and `dG-(E)/dE` is the
energy derivative of the advanced Green's function. In the same run, the
Kubo–Greenwood trace gives the longitudinal conductivity σ_xx, which vanishes
in the gap.

## Run

We assume LSQuant is built as in [Getting started](getting-started.md). We run
the example, inspect the output and plot it:

```bash
export OMP_NUM_THREADS=1
python3 -m lsquant run examples/03_haldane_hall/run_haldane.json --out /tmp/haldane
python3 -m lsquant inspect /tmp/haldane
python3 examples/03_haldane_hall/plot.py /tmp/haldane
```

The run file expands both conductivities in 128 moments on the energy window
[−4.5, 4.5] eV. The number of moments sets the energy resolution, like the
smearing of a DFT density of states. The run file requests a broadening of
10 meV, but 128 moments resolve only about πa/M ≈ 0.11 eV, with a = 4.5 eV the
half-width of the window and M = 128. LSQuant therefore uses all 128 moments
and the coarser resolution.
<!-- source: examples/03_haldane_hall/run_haldane.json (moments 128, broadening_meV 10); resolution πa/M from docs/explanation/kpm.md and python/lsquant/reconstruct.py:89-96, 122-147 (min(ceil(π/η), M) moments kept); legacy_text_round_trip keeps a = 4.5 eV: runfile.py:2192-2233, lsquant/io/chebmom_file.cpp:236-237, 261-264 -->

The trace runs exactly over all 512 states of the honeycomb supercell, so the
curves carry no random-sampling noise. This exact trace costs N × M
matrix–vector products, here 512 × 128, which is why even this small model
takes minutes. Large systems replace the N basis states with R random vectors;
see [The kernel polynomial method](../explanation/kpm.md#exact-and-random-phase-traces).
<!-- source: examples/03_haldane_hall/run_haldane.json ("states": 512, "trace": "exact"); cost N × M from docs/explanation/kpm.md (Exact and random-phase traces) -->

From these moments the run builds four curves: the full Kubo–Bastin Hall
conductivity (`bastin`), its Fermi-sea and Fermi-surface parts (`sea`,
`surf`), and σ_xx (`greenwood`). On one thread of a heavily loaded shared
node, the run took 334 s.
<!-- source: examples/03_haldane_hall/README.md (Runtime) -->

## Result

![Hall conductivity σ_xy (total, Fermi sea, Fermi surface) and longitudinal conductivity σ_xx of the Haldane model, in units of e²/h](../../examples/03_haldane_hall/fig_haldane_hall.png)

*Fig. 1. Conductivities of the Haldane model ($t_1 = -1$ eV, $t_2 = 0.15$ eV,
$\phi = \pi/2$, 512 states, $M = 128$), in units of $e^2/h$. (a) The Hall
conductivity $\sigma_{xy}$ forms a flat plateau at one $e^2/h$ through the
gap; the Fermi-sea part carries all of it, and the Fermi-surface part stays
at zero there. (b) The longitudinal conductivity $\sigma_{xx}$ is zero in the
gap and finite in the bands.*
<!-- source: caption from examples/03_haldane_hall/README.md (Step 2: figure alt text and panel description) -->

Fig. 1 shows the Hall conductivity on a plateau at $+1\,e^2/h$ through the
gap. From this plateau the plot script prints $C = +1.017$, 1.7 % above the
Chern number $C = 1$. The plateau persists where there are no states at the
Fermi level, because the Fermi sea carries all of it.
<!-- source: examples/03_haldane_hall/README.md (Step 2, Runtime: C = +1.017); the README gives no cause for the deviation -->

[Example 03](../../examples/03_haldane_hall/README.md) also compares the
plateau with the Chern number computed in k-space across the sublattice mass.
[Example 13](../../examples/13_bastin_sea_surface/README.md) checks the split
into sea and surface parts on the same model, and
[Kubo formulas](../explanation/kubo-formulas.md) explains the two formulas.
<!-- source: examples/03_haldane_hall/README.md (fig_phase_diagram.png: real-space plateau vs Fukui–Hatsugai Chern number across M) -->

## References

[1] F. D. M. Haldane, *Model for a quantum Hall effect without Landau levels*,
Phys. Rev. Lett. **61**, 2015 (1988).
<!-- source: examples/03_haldane_hall/README.md (References) -->
