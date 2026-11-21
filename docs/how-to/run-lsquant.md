# Running LSQuant: from the run file to the curve

In this guide we follow one calculation from the run file to the curve. We
name each choice on the way and leave the worked details to the examples in
`examples/`. Each of them makes the choice on a physical system and checks the
result against an exact one. We assume that LSQuant is built and on
`PYTHONPATH`, as in [Getting started](../tutorials/getting-started.md).

## The three steps of a calculation

Every LSQuant calculation has three steps:

1. We write a **run file**, a JSON file that names the model, the basis, the
   quantity to compute and how to expand it.
2. `python -m lsquant run` computes the **Chebyshev moments** of that quantity
   and writes them to disk. This is the expensive step.
3. A **reconstruction** turns the moments into a curve, such as a DOS or a
   conductivity σ(E), with the Jackson kernel and an energy broadening. `run`
   does this at the end, and `python -m lsquant reconstruct` repeats it later
   from the saved moments, without applying the Hamiltonian again.

For the DOS of a tight-binding model stored as sparse matrices, a minimal run
file is:

```json
{
  "model": {"kind": "csr", "dir": "goldens/stage0/haldane/operators", "label": "haldane"},
  "word": "Tr[delta(E-H)]",
  "expansion": {"moments": 256},
  "reconstruction": {"broadening_meV": 10.0}
}
```

It reads `haldane.HAM.CSR`, the Hamiltonian of the Haldane model of
[example 03](../../examples/03_haldane_hall/README.md), which the repository
keeps next to its reference outputs. It expands the DOS on the default spectral
window, a Gershgorin bound of the spectrum, with an exact trace. Every key that
we leave out takes its default value; the
[run-file reference](../reference/run-file.md) lists them all. Relative paths in
a run file are tried against the LSQuant source tree first, so for our own data
we give absolute paths.

Before we spend compute time on a run file, we let LSQuant validate it and fill
in the defaults:

```bash
python -m lsquant inspect my_run.json
```

A run file that `inspect` refuses also fails in `run`. Operator names outside
the standard set, however, pass `inspect` and fail in `run` (see the
[run-file reference](../reference/run-file.md#word-and-words)).

## Choosing the model

The key `model.kind` selects the source of the Hamiltonian:

| Starting point | `model.kind` | Example |
|---|---|---|
| A model we write in Python | `generator` (the drivers of examples 01–06) or `csr` (sparse-matrix files we write ourselves) | [01](../../examples/01_dos_1d_chain/README.md), [03](../../examples/03_haldane_hall/README.md) |
| A Wannier90 `_hr.dat` file | `wannier` | [14](../../examples/14_wannier/README.md) |
| A SIESTA calculation with spin–orbit coupling | `siesta` | [07](../../examples/07_spin_hall_bi_bilayer/README.md), [08](../../examples/08_spin_hall_pdse2_monolayer/README.md) |

[Example 14](../../examples/14_wannier/README.md) reads the
Su–Schrieffer–Heeger chain from a hand-written `_hr.dat` file and reproduces its
two analytic bands. For a SIESTA model, [Use a SIESTA calculation](use-siesta.md)
lists what the SIESTA run must save and which settings LSQuant accepts.

## Choosing the basis

The key `basis.name` sets how LSQuant represents $H$: as a supercell in real
space (`real`), as one block per k-point of the unit cell (`bloch`), or by the
exact eigenstates of each block (`eigen`). For a clean periodic system the three
give the same moments, and `bloch` or `eigen` are much cheaper. We use `real`
when disorder or a defect breaks the translation symmetry.

[Example 10](../../examples/10_bases_compare/README.md) computes one graphene
Hamiltonian in the three bases; their moments agree to 10⁻¹⁴.

## Working with an overlap matrix

Local-orbital DFT bases are non-orthogonal, with an overlap matrix $S$. For
them we set `metric.policy` to `"nonorthogonal"`, and LSQuant works either with
$S^{-1/2} H S^{-1/2}$ (the Löwdin route) or with $S^{-1} H$ (the generalised
route).

[Example 11](../../examples/11_nonorthogonal_toy/README.md) solves a
non-orthogonal chain on both routes and reproduces the exact generalised
eigenvalue problem to 10⁻¹³.
[Non-orthogonal bases](../explanation/nonorthogonal-bases.md) explains why $S$
matters for velocities and currents.

## Adding disorder

For Wannier and SIESTA models we set the Anderson disorder strength
`basis.disorder.W`, in eV, and its seed `basis.disorder.seed`. For SIESTA
models with spin–orbit coupling, `basis.disorder.kind` chooses a
spin-independent potential (`"scalar"`, the default) or a random Zeeman field
(`"magnetic"`). To average over disorder we run one run file per seed. Sparse-matrix models take their disorder
from the files themselves, and the generators build their own.

[Example 02](../../examples/02_localization_1d_chain/README.md) follows the
Anderson localisation of a one-dimensional chain, and
[example 05](../../examples/05_charge_transport_square/README.md) the resistance
of a disordered square lattice.
[Example 12](../../examples/12_disorder_nonorthogonal/README.md) adds disorder
to a SIESTA basis with an overlap matrix and explains why the default
`metric.tol` must then be relaxed.

## Choosing the quantity and the reconstruction

Each quantity is a trace of operators and functions of $H$, which LSQuant calls
the **observable word**. LSQuant expands the word into Chebyshev moments, and a
reconstruction formula turns the moments into the physical curve. In the words
of the table, `VX` and `VY` are the velocity operators $v_x$ and $v_y$, `JXSZ`
is the spin current $J^z_x$, `dG-(E)/dE` is the energy derivative of the
advanced Green's function, `U(t)` is the time-evolution operator and `^+` the
adjoint. In `P:PX`, the operator `PX` of example 04 is the projector
$P_x = (1 + S_x)/2$ onto spin along $+x$; the
[observable-word reference](../reference/observable-words.md) explains this
syntax.

| Quantity | Word | Formula | Example |
|---|---|---|---|
| DOS | `Tr[delta(E-H)]` | `dos` | [01](../../examples/01_dos_1d_chain/README.md) |
| DOS projected on a region or an orbital set | `Tr[A delta(E-H)]` | `dos` | [15](../../examples/15_spatial_operators/README.md) |
| σ_xx (Kubo–Greenwood) | `Tr[VX delta(E-H) VX delta(E-H)]` | `greenwood` | [05](../../examples/05_charge_transport_square/README.md), [10](../../examples/10_bases_compare/README.md) |
| σ_xy (Kubo–Bastin) | `Tr[VX delta(E-H) VY dG-(E)/dE]` | `bastin` | [03](../../examples/03_haldane_hall/README.md) |
| σ_xy split into Fermi-sea and Fermi-surface terms | as above | `sea`, `surf` | [13](../../examples/13_bastin_sea_surface/README.md) |
| Spin Hall conductivity | `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` | `shc` | [07](../../examples/07_spin_hall_bi_bilayer/README.md) |
| Spreading of a wave packet (mean-square displacement) | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | `time_correlation` | [02](../../examples/02_localization_1d_chain/README.md), [15](../../examples/15_spatial_operators/README.md) |
| Spin dynamics | `Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]` and related | `row_m0`, `time_correlation` | [04](../../examples/04_spin_precession_magnetic_chain/README.md), [06](../../examples/06_spin_relaxation_square/README.md) |

The formula `row_m0` keeps the $m = 0$ moment of each time step, which is the
spin integrated over energy, and writes no curve file. The command `words`
lists the word shapes that LSQuant accepts, and `words describe` shows how one
word is compiled:

```bash
python -m lsquant words
python -m lsquant words describe "Tr[VX delta(E-H) VY dG-(E)/dE]"
```

> **Known limitation.** With `model.kind` set to `csr` or `wannier`, the
> formula `shc` fails at the end of `run`. We run with `bastin` and then apply
> `shc` with `reconstruct`, as described in the
> [run-file reference](../reference/run-file.md#reconstruction).

## Restricting the energy window

For a SIESTA model, the `projection` section keeps only the states in an energy
window, so the Chebyshev expansion covers a much narrower spectrum.
[Example 16](../../examples/16_projection_downfold/README.md) projects a Bi
bilayer Hamiltonian onto [−2, 2] eV around the Fermi level and recovers the
exact DOS of the kept bands.

## Choosing the number of moments and the trace

The finest broadening we can resolve is $\pi a/M$, where $a$ is the half-width
of the spectral window and $M$ the number of moments. We therefore choose $M$
from the resolution we need. A broadening finer than $\pi a/M$ is replaced by
$\pi a/M$ without a warning, so a wide window, as in most DFT Hamiltonians,
needs many moments. [The kernel polynomial method](../explanation/kpm.md)
explains the kernel and the broadening.

With `expansion.trace = "exact"`, LSQuant sums over every basis state. For large
real-space systems we use `"random"` with a few random-phase vectors, set by
`expansion.states`. With $R$ vectors in a system of dimension $N$, the relative
error of this estimate falls as $1/\sqrt{RN}$, and it is larger where the DOS
is small.

## Running the calculation

We run the file with:

```bash
python -m lsquant run my_run.json --out out/my_run
```

The option `--out` overrides `output.dir` of the run file. The output directory
holds the moment files, the curves and `run.resolved.json`, the run file with
every default filled in. The [output-file reference](../reference/output-files.md)
describes each file.

On a shared node we set `OMP_NUM_THREADS=1`: threads that wait for each other
can make a small run 10 to 30 times slower than one thread.
[Run on a cluster](run-on-a-cluster.md) gives the details.

## Inspecting and reconstructing a finished run

To see what produced a directory and to draw a new curve from its moments, we
use:

```bash
python -m lsquant inspect out/my_run
python -m lsquant reconstruct out/my_run --formula dos --broadening-meV 120
```

`reconstruct` gives a new curve in seconds, with another broadening or with
another formula that the same word allows.
[Example 09](../../examples/09_inspect_reconstruct/README.md) reconstructs one
set of moments at several broadenings and checks each curve against the exact
DOS.

## Working without a run file

A run file covers the standard calculations. To change one step of a
calculation, we drive LSQuant from Python, as in
[example 17](../../examples/17_python_api/README.md). An AI assistant can also
run LSQuant through its MCP server:
[example 18](../../examples/18_mcp_session/README.md) shows a session, and
[Connect an AI assistant](mcp-server.md) shows how to set it up.
