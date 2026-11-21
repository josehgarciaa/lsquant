# LSQuant: linear-scaling quantum transport

LSQuant computes spectral and transport properties of large tight-binding and
first-principles Hamiltonians with the kernel polynomial method (KPM). It
expands each quantity in Chebyshev polynomials of the Hamiltonian, so it
avoids diagonalisation and its cost grows linearly with the number of orbitals.
You can compute:

- the density of states, total or projected;
- the longitudinal conductivity σ_xx (Kubo–Greenwood);
- the Hall and spin Hall conductivities σ_xy and σ^z_xy (Kubo–Bastin);
- the mean-square displacement and real-time spin dynamics.

Models come from SIESTA (with spin–orbit coupling), from Wannier90, from sparse
matrices you build yourself, or from the generators in `examples/`. Small
systems can also be solved exactly, so you can check a KPM result against an
exact reference.

## Quickstart

You need a C++17 compiler, CMake 3.21 or later for the presets (3.16 or later for a manual build), Eigen, FFTW3, GSL, OpenMP and
Python 3.10 or later. [Install LSQuant](docs/how-to/install.md) gives the full
list and the commands for each system.

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

The `golden` preset is the portable build, which gives the same numbers on
every machine. The run computes the density of states of a one-dimensional
chain in about 2.5 minutes on one core, and writes the Chebyshev moments, the
reconstructed curves and `run.resolved.json` to `/tmp/chain_dos`. The last
command draws the figures. [Getting started](docs/tutorials/getting-started.md)
explains each step.

## Documentation

The documentation follows the [Diátaxis](https://diataxis.fr) layout. Start at
the [documentation index](docs/index.md).

- **Tutorials:** [Getting started](docs/tutorials/getting-started.md),
  [Haldane Hall conductivity](docs/tutorials/haldane-hall.md),
  [Spin Hall conductivity from SIESTA](docs/tutorials/spin-hall-siesta.md).
- **How-to guides:** [Install](docs/how-to/install.md),
  [Running LSQuant](docs/how-to/run-lsquant.md),
  [Use a SIESTA calculation](docs/how-to/use-siesta.md),
  [Run on a cluster](docs/how-to/run-on-a-cluster.md),
  [Check a build](docs/how-to/testing.md),
  [Connect an AI assistant (MCP)](docs/how-to/mcp-server.md).
- **Reference:** [Run file](docs/reference/run-file.md),
  [Command line](docs/reference/cli.md),
  [Observable words](docs/reference/observable-words.md),
  [Output files](docs/reference/output-files.md),
  [Build presets](docs/reference/build-presets.md),
  [Python API](docs/reference/python-api.md).
- **Explanation:** [How LSQuant works](docs/explanation/unifying-idea.md),
  [The kernel polynomial method](docs/explanation/kpm.md),
  [Kubo formulas](docs/explanation/kubo-formulas.md),
  [Non-orthogonal bases](docs/explanation/nonorthogonal-bases.md).

## Examples

Each folder in [`examples/`](examples/) holds a run file, a model generator or
input files, a plotting or checking script, and a README that walks through the
physics. Examples 09–18 are the use cases of [Running LSQuant](docs/how-to/run-lsquant.md):

| Example | Model | Quantity |
|---|---|---|
| `01_dos_1d_chain` | 1D chain | DOS |
| `02_localization_1d_chain` | disordered 1D chain | Anderson localisation |
| `03_haldane_hall` | Haldane model | quantised σ_xy |
| `04_spin_precession_magnetic_chain` | magnetic chain | spin precession |
| `05_charge_transport_square` | disordered square lattice | σ_xx, diffusion |
| `06_spin_relaxation_square` | Rashba and Elliott–Yafet models | spin relaxation |
| `07_spin_hall_bi_bilayer` | Bi(111) bilayer, SIESTA | DOS, spin Hall conductivity |
| `08_spin_hall_pdse2_monolayer` | PdSe₂ monolayer, SIESTA | DOS, spin Hall conductivity |
| `09_inspect_reconstruct` | 1D chain | one set of moments at several broadenings |
| `10_bases_compare` | graphene, Wannier90 format | the same DOS and σ_xx in three bases |
| `11_nonorthogonal_toy` | non-orthogonal chain | Löwdin and generalised routes |
| `12_disorder_nonorthogonal` | Bi(111) bilayer, SIESTA | Anderson disorder with an overlap matrix |
| `13_bastin_sea_surface` | Haldane model | σ_xy split into Fermi-sea and Fermi-surface terms |
| `14_wannier` | SSH chain, Wannier90 | DOS from a `_hr.dat` file |
| `15_spatial_operators` | biased bilayer chain | layer-resolved DOS, mean-square displacement |
| `16_projection_downfold` | Bi(111) bilayer, SIESTA | spectral projection onto an energy window |
| `17_python_api` | graphene | the same calculation from Python |
| `18_mcp_session` | graphene | a session driven through the MCP server |
| `disorder_kinds` | Bi(111) bilayer, SIESTA, 3×3 supercell | Kramers splitting under scalar, magnetic and independent disorder |

## Check your build

`scripts/check_goldens.sh` runs a build on committed inputs and compares the
outputs with the reference outputs in `goldens/`. On the reference machine, a
preset build reports two failures. Both come from checks of the C++ source,
not from the comparison of the numbers. On other machines, a few `DIFFERS`
lines can also appear, so a higher failure count alone does not mean the build
is broken. [Check a build](docs/how-to/testing.md) explains how to read each
line and when such differences are acceptable.

## License

See [LICENSE](LICENSE).
