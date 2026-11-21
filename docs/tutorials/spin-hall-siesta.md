# The spin Hall conductivity of a Bi(111) bilayer from SIESTA

In this tutorial we start from a SIESTA calculation of a two-dimensional
Bi(111) bilayer with spin–orbit coupling and compute its intrinsic spin Hall
conductivity $\sigma^z_{xy}(E)$. This is the second task of most readers: the
workflow of [Getting started](getting-started.md), applied to a
first-principles Hamiltonian in a non-orthogonal basis.

LSQuant evaluates the Kubo–Bastin formula with a spin current in place of one
velocity, the observable word `Tr[JXSZ delta(E-H) VY dG-(E)/dE]`. Here
`JXSZ` is the spin current $J^z_x = \frac12\{v_x, \sigma_z\}$. The result is in
units of $e^2/h$ for a spin current defined with $\sigma_z$, without the factor
$\hbar/2$.

The SIESTA input and output files of the bilayer are in
`examples/07_spin_hall_bi_bilayer/siesta/`. We run the quick versions of the DOS
and spin Hall calculations into one directory and plot them:

```bash
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_dos_quick.json --out /tmp/bi_quick
python3 -m lsquant run examples/07_spin_hall_bi_bilayer/run_shc_quick.json --out /tmp/bi_quick
python3 examples/07_spin_hall_bi_bilayer/plot.py /tmp/bi_quick
```

The quick DOS run took 12 minutes with 8 threads on a lightly loaded node. The
run files read the Hamiltonian and the overlap matrix $S$ from the SIESTA
files, work in the `eigen` basis on a k mesh, and treat $S$ with the Löwdin
route.

The full lesson is in [example 07](../../examples/07_spin_hall_bi_bilayer/README.md),
and [example 08](../../examples/08_spin_hall_pdse2_monolayer/README.md) repeats
it for a PdSe₂ monolayer. [Use a SIESTA calculation](../how-to/use-siesta.md)
shows how to prepare a new SIESTA calculation. It also explains how the
number of moments limits the energy resolution for the wide spectrum of a DFT
Hamiltonian.
