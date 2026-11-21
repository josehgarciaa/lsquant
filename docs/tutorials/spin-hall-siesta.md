# The spin Hall conductivity of a Bi(111) bilayer from SIESTA

The free-standing Bi(111) bilayer is a two-dimensional quantum spin Hall
insulator with strong spin–orbit coupling. In this tutorial we start from a
SIESTA calculation of the bilayer and compute its intrinsic spin Hall
conductivity $\sigma^z_{xy}(E)$ and its density of states. We follow the
workflow of [Getting started](getting-started.md), now with a
first-principles Hamiltonian in the non-orthogonal atomic-orbital basis of
SIESTA.
<!-- source: examples/07_spin_hall_bi_bilayer/README.md (spin Chern number 2, Z2 = 1) -->

## SIESTA calculation

The SIESTA input and output files of the bilayer are in
`examples/07_spin_hall_bi_bilayer/siesta/`. The input `bi.fdf` describes two Bi
atoms per cell with a DZP basis, the PBE functional, `Spin spin-orbit` and a
24 × 24 Monkhorst–Pack k-grid. SIESTA saved the Hamiltonian H and the overlap
matrix S in `bi.TSHS`, and the basis orbitals in `Bi.ion.xml`. LSQuant reads
these files directly; it runs no SIESTA step itself.
<!-- source: examples/07_spin_hall_bi_bilayer/siesta/bi.fdf -->

## Spin Hall conductivity in LSQuant

LSQuant computes the spin Hall conductivity with the Kubo–Bastin formula, with
a spin current in place of one velocity. In the run file this observable is
the trace `Tr[JXSZ delta(E-H) VY dG-(E)/dE]`, where `JXSZ` is the spin current
$J^z_x = \frac12\{v_x, \sigma_z\}$. The result is in units of $e^2/h$ for a spin
current defined with the Pauli matrix $\sigma_z$, without the factor
$\hbar/2$.

## Run

We run the quick versions of the DOS and spin Hall calculations into one
directory and plot them:

```bash
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_dos_quick.json --out /tmp/bi_quick
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_shc_quick.json --out /tmp/bi_quick
python3 examples/07_spin_hall_bi_bilayer/plot.py /tmp/bi_quick
```

The run files read H and S from the SIESTA files and work in the `eigen`
basis on a 24 × 24 k-mesh, the same mesh as the SIESTA run. They remove the
overlap with the Löwdin route, which turns the generalised eigenvalue problem
with S into an ordinary one. Both expand the trace in M = 256 moments with a
requested broadening of 20 meV. The quick DOS run took 12 minutes with 8
threads on a lightly loaded node. On one thread of a heavily loaded node, the
DOS and spin Hall runs took 18 minutes each.
<!-- source: examples/07_spin_hall_bi_bilayer/run_dos_quick.json, run_shc_quick.json (basis eigen, supercell [24,24,1], metric.route lowdin, moments 256, broadening_meV 20); times: examples/RUNLOG.md (12 min, 8 threads), examples/07_spin_hall_bi_bilayer/README.md (Runtime: 18 min each, one thread) -->

## Result

Both run files also diagonalise H(k) at every k-point, with `oracle: true`,
which gives an exact reference for the KPM curve. With a Gaussian of 20 meV on
the quick 24 × 24 mesh, this exact calculation gives
$\sigma^z_{xy}(E_F = 0) = 1.144\ e^2/h$. It also finds 30 occupied bands below
a gap of 0.304 eV, and Chern numbers $C_\pm = \pm 1$ of the two spin sectors,
so the spin Chern number is 2.
<!-- source: examples/07_spin_hall_bi_bilayer/README.md (Step 3: exact_shc.json, Gaussian of 20 meV, 24x24x1; spin_chern.json); "oracle": true in run_*_quick.json -->

The exact value is not 2 because spin–orbit coupling mixes the two spin
directions. The spin component $s_z$ is then not conserved, and
$\sigma^z_{xy}$ is not quantised even though the spin Chern number is.
<!-- source: REVISION_BRIEF_3 (fact correction for this page); consistent with examples/07 README, which gives 1.144 and 1.0798 e²/h for the exact value, not 2 -->

The KPM curve of the quick run gives $\sigma^z_{xy}(E_F = 0) \approx 0.50\
e^2/h$, below the exact 1.144 $e^2/h$. With M = 256 on the wide SIESTA window,
the curve resolves only about 0.6 eV, twice the 0.304 eV gap, so it cannot
show the gap. Resolving 20 meV on this window needs about 7600 moments.
<!-- source: examples/RUNLOG.md (Example 07: "The quick KPM curve gives 0.500 at E_F = 0 ... resolves about 0.6 eV, more than the gap"); docs/how-to/use-siesta.md (Number of moments and broadening: a ≈ 48.6 eV, 0.6 eV, about 7600 moments) -->

[Example 07](../../examples/07_spin_hall_bi_bilayer/README.md) plots the KPM
curve against this reference, and
[example 08](../../examples/08_spin_hall_pdse2_monolayer/README.md) repeats the
calculation for a PdSe₂ monolayer.
[Use a SIESTA calculation](../how-to/use-siesta.md) explains how to prepare a
SIESTA calculation of another material and how to choose the number of
moments for the wide spectrum of a DFT Hamiltonian.
