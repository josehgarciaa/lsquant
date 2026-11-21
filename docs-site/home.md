---
hide:
  - navigation
  - toc
---

# LSQuant

LSQuant computes the electronic and transport properties of materials models
too large to diagonalise. It reads a SIESTA calculation, a Wannier90
Hamiltonian or a tight-binding model, and returns the density of states, the
longitudinal and Hall conductivity, the spin Hall conductivity, and the time
evolution of spin and charge. The cost grows linearly with the number of
orbitals, so supercells with disorder or defects of millions of orbitals are
within reach of one workstation.

[Get started](quickstart.md){ .md-button .md-button--primary }
[Install](docs/how-to/install.md){ .md-button }
[Source on GitHub](https://github.com/josehgarciaa/lsquant){ .md-button }

## What you can compute

<div class="grid cards" markdown>

-   ![Hall conductivity plateau of the Haldane model](examples/03_haldane_hall/fig_haldane_hall.png)

    **Quantized Hall conductivity**

    A Chern insulator carries a Hall plateau of one $e^2/h$ through its gap,
    computed with the Kubo–Bastin formula in real space.

    [:octicons-arrow-right-24: Example 3](examples/03_haldane_hall/README.md)

-   ![Conductivity and mean free path against disorder](examples/05_charge_transport_square/fig_transport.png)

    **Charge transport with disorder**

    Diffusion coefficient, conductivity and mean free path of a disordered
    metal, from the spreading of wave packets.

    [:octicons-arrow-right-24: Example 5](examples/05_charge_transport_square/README.md)

-   ![Spin relaxation under disorder](examples/06_spin_relaxation_square/fig_dp.png)

    **Spin relaxation**

    Dyakonov–Perel and Elliott–Yafet relaxation of a spin in a
    spin–orbit-coupled lattice, and the anisotropy of the spin lifetimes.

    [:octicons-arrow-right-24: Example 6](examples/06_spin_relaxation_square/README.md)

-   ![Spin precession in a magnetic chain](examples/04_spin_precession_magnetic_chain/fig_precession.png)

    **Spin dynamics**

    A spin precesses at the Larmor frequency set by the exchange field. The
    computed curves match the analytic result to machine precision.

    [:octicons-arrow-right-24: Example 4](examples/04_spin_precession_magnetic_chain/README.md)

-   ![Localization length against energy](examples/02_localization_1d_chain/fig_localization.png)

    **Anderson localisation**

    The localization length of a disordered wire, compared with the Born
    approximation.

    [:octicons-arrow-right-24: Example 2](examples/02_localization_1d_chain/README.md)

-   ![Density of states reconstructed from Chebyshev moments](examples/01_dos_1d_chain/fig_dos.png)

    **Density of states**

    The density of states from a Chebyshev expansion, and how the number of
    moments sets the energy resolution.

    [:octicons-arrow-right-24: Example 1](examples/01_dos_1d_chain/README.md)

</div>

## From your DFT calculation

- **SIESTA:** LSQuant reads the Hamiltonian and overlap that SIESTA saves,
  with spin–orbit coupling and a non-orthogonal basis, and computes the spin
  Hall conductivity of the material. Examples
  [7, Bi(111) bilayer](examples/07_spin_hall_bi_bilayer/README.md) and
  [8, PdSe₂ monolayer](examples/08_spin_hall_pdse2_monolayer/README.md).
- **Wannier90:** LSQuant reads a `_hr.dat` file and builds a supercell, a
  k-point mesh or the exact eigenstates from it
  ([example 14](examples/14_wannier/README.md)).
- **Disorder and defects:** Anderson disorder, scalar or magnetic, in
  real-space supercells, averaged over realisations
  ([example 12](examples/12_disorder_nonorthogonal/README.md),
  [disorder kinds](examples/disorder_kinds/README.md)).

## How a calculation runs

You describe the calculation in one JSON run file: the model, the quantity,
the number of Chebyshev moments, which sets the energy resolution, and the
number of random vectors, which sets the statistical noise. You then run
`python -m lsquant run`. LSQuant stores the moments, so you can rebuild the
curve at another resolution in seconds, without recomputing.
[Run a calculation](docs/how-to/run-lsquant.md) explains each choice.
