# LSQuant: linear-scaling quantum transport

**Documentation, examples and gallery: <https://josehgarciaa.github.io/lsquant/>**

LSQuant computes the spectral and transport properties of systems too large to
diagonalise. It gives the total and projected density of states, the
longitudinal conductivity σ_xx from the Kubo–Greenwood formula, and the Hall
and spin Hall conductivities σ_xy and σ^z_xy from the Kubo–Bastin formula. It
also follows the mean-square displacement and the real-time spin dynamics.
The method is the kernel polynomial method (KPM): each quantity is expanded in
Chebyshev polynomials of the Hamiltonian, so no diagonalisation is needed and
the cost grows linearly with the number of orbitals.

The Hamiltonian comes from a SIESTA calculation with spin–orbit coupling, a
Wannier90 `_hr.dat` file, sparse-matrix files you write yourself, or the model
generators in `examples/`. LSQuant is aimed at 2D materials, spin–orbit and
topological materials, and disordered crystals. Small systems can also be
solved exactly, so you can check a KPM result against an exact reference.

## Quantities by model source

The quantities you can compute depend on where the Hamiltonian comes from:

| Hamiltonian from | DOS | σ_xx | σ_xy | Spin Hall σ^z_xy | Time-dependent (spreading, spin dynamics) |
|---|---|---|---|---|---|
| SIESTA with spin–orbit coupling (`siesta`) | yes | no | no | yes | no |
| OpenMX with spin–orbit coupling (`openmx`) | yes | no | no | yes | no |
| Wannier90 `_hr.dat` (`wannier`) | yes | yes | not documented | no | not documented |
| Wannier90 with `params.velocity_mode`, `spin` or `fermi_energy` (`wannier`) | yes | no | no | yes | no |
| PAOFLOW model file (`paoflow`) | yes | no | no | yes | no |
| Sparse-matrix files (`csr`) | yes, also projected | yes | yes | yes¹ | yes |
| Model generators of examples 01–06 (`generator`) | yes | yes | yes | not documented | yes |

<!-- source, siesta row: the clean and disordered executors accept only Tr[delta(E-H)] and Tr[JXSZ delta(E-H) VY dG-(E)/dE] (python/lsquant/apps/stage6b_siesta.py:24-25, 179-182; apps/stage6c_disorder.py:86-89); the projection executor refuses a time section (python/lsquant/projection/runfile.py:136); docs/reference/run-file.md#siesta-executors; examples 07, 08 -->
<!-- source, wannier row: DOS and sigma_xx in examples/10_bases_compare/run_*.json and examples/14_wannier/run_ssh_dos.json (dos, greenwood); the model provides only VX, VY, VZ (python/lsquant/runfile.py:1637; docs/reference/observable-words.md, "Names each model provides"), so no spin-current operator; no documented sigma_xy or time-dependent run -->
<!-- source, openmx, wannier with params and paoflow rows: they run through the SIESTA clean and projection executors (python/lsquant/runfile.py, uses_metric_executor; models/orthogonal.py), which accept only Tr[delta(E-H)] and Tr[JXSZ delta(E-H) VY dG-(E)/dE]; docs/how-to/use-wannier90.md, docs/how-to/use-paoflow.md -->
<!-- source, csr row: projected DOS examples/15_spatial_operators/run_layer_dos.json; sigma_xx and sigma_xy examples/03_haldane_hall/run_haldane.json (greenwood, bastin) and examples/13_bastin_sea_surface; time-dependent examples/04 run_spin.json and examples/15 run_msd.json; spin Hall through reconstruct: docs/reference/run-file.md, known limitation "shc in csr and wannier runs" (runfile.py:1714, 2228, 2307, 2438) -->
<!-- source, generator row: make_chain (tutorial 1, DOS), make_haldane (tutorial 3, kubo_bastin and kubo_greenwood in python/lsquant/apps/tutorial03.py:28-31, 255), make_square (tutorial 5, sigma_xx, examples/05 README), make_disordered_chain, make_chain_mag, make_ey/make_rashba (tutorials 2, 4, 6, time-dependent); docs/reference/run-file.md#generator-models -->

¹ Give the spin current as an operator file `<label>.JXSZ.CSR`, run with the
`bastin` formula, and apply `shc` afterwards with `python -m lsquant
reconstruct`, as the [run-file reference](docs/reference/run-file.md#reconstruction)
describes.

## Requirements

You need a C++17 compiler with OpenMP, CMake 3.21 or later for the presets
(3.16 or later for a manual build), Eigen, FFTW3, GSL and Python 3.10 or later
with NumPy, SciPy, SymPy and pybind11. [Install LSQuant](docs/how-to/install.md)
gives the full list and the commands with and without root access.

## Quickstart

Build LSQuant with the `golden` preset, the portable build that gives the same
numbers on every machine. Then compute the density of states of a
one-dimensional chain and plot it:

```bash
git clone https://github.com/josehgarciaa/lsquant.git
cd lsquant
pip install numpy scipy sympy pybind11 matplotlib
cmake --preset golden
cmake --build --preset golden
export PYTHONPATH="$PWD/build/golden/lsquant_core:$PWD/python"

export OMP_NUM_THREADS=1
python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/chain_dos
python3 -m lsquant inspect /tmp/chain_dos
python3 examples/01_dos_1d_chain/plot.py /tmp/chain_dos
```

The run computes the density of states of three rings of 256, 1024 and 4096
sites with 512 Chebyshev moments, and takes about 2.5 minutes on one core.
It uses the exact trace, which costs N × M matrix–vector products for N
orbitals and M moments; large systems use R random vectors instead, and their
cost grows linearly with N. When the run finishes, `/tmp/chain_dos` holds the
Chebyshev moments, the density-of-states curves and `run.resolved.json` (the
run file with every default filled in). The last command draws the figures
into the same directory. [Getting started](docs/tutorials/getting-started.md)
explains each step, and [Install LSQuant](docs/how-to/install.md) covers a
cluster without root access.
<!-- source: ring sizes, moments, exact trace: examples/01_dos_1d_chain/run_dos.json; N x M and R M products: docs/explanation/kpm.md, "Exact and random-phase traces" and "Cost" -->

## Documentation

The full documentation is at <https://josehgarciaa.github.io/lsquant/>, and the
same pages are in [`docs/`](docs/index.md). Start with
[Getting started](docs/tutorials/getting-started.md) and the tutorials on the
[Haldane model](docs/tutorials/haldane-hall.md) and on the
[Bi(111) bilayer from SIESTA](docs/tutorials/spin-hall-siesta.md). The guides
[Running LSQuant](docs/how-to/run-lsquant.md),
[Use a SIESTA calculation](docs/how-to/use-siesta.md) and
[Run on a cluster](docs/how-to/run-on-a-cluster.md) cover your own
calculations. The [run-file reference](docs/reference/run-file.md) lists every
key, and [The kernel polynomial method](docs/explanation/kpm.md) gives the
background on the method.

## Examples

Each folder in [`examples/`](examples/) holds a run file, a model generator or
input files, a plotting or checking script, and a README that walks through
the physics.

| Example | Model | Quantity |
|---|---|---|
| `01_dos_1d_chain` | 1D chain | DOS |
| `02_localization_1d_chain` | disordered 1D chain | mean-square displacement, Anderson localisation length |
| `03_haldane_hall` | Haldane model | quantised σ_xy |
| `04_spin_precession_magnetic_chain` | magnetic chain | spin precession |
| `05_charge_transport_square` | disordered square lattice | σ_xx, diffusion |
| `06_spin_relaxation_square` | Rashba and Elliott–Yafet models | spin relaxation |
| `07_spin_hall_bi_bilayer` | Bi(111) bilayer, SIESTA | DOS, spin Hall conductivity |
| `08_spin_hall_pdse2_monolayer` | PdSe₂ monolayer, SIESTA | DOS, spin Hall conductivity |
| `09`–`18` | 1D, SSH, non-orthogonal and bilayer chains, graphene, Haldane model, Bi(111) bilayer | use cases of [Running LSQuant](docs/how-to/run-lsquant.md): bases, overlap matrices, disorder, Wannier90 input, Fermi-sea and Fermi-surface terms of σ_xy, layer-resolved DOS, projections, Python and MCP |
| `disorder_kinds` | Bi(111) bilayer, SIESTA, 3×3 supercell | Kramers splitting under scalar, magnetic and independent disorder |
<!-- source: example 02 reads the localisation length from the plateau of the mean-square displacement (examples/02_localization_1d_chain/README.md:40-42) -->

## Check your build

To check a build, run `scripts/check_goldens.sh`. The script runs the build on
inputs stored in the repository and compares its results with the reference
outputs in `goldens/`. With the reference toolchain, a `golden` build passes
every stage of the check. On other machines, run it with `GOLDEN_TOL=1`: a few
`DIFFERS` lines can remain there without the build being broken.
[Check a build against the reference outputs](docs/how-to/testing.md) explains
how to read each line and when such differences are acceptable.

## License

See [LICENSE](LICENSE).
