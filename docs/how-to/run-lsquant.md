# Running LSQuant

Here we go through the choices of an LSQuant calculation in the order in which
they arise in a materials calculation: the model, the quantity, the
representation, then the convergence parameters. For each choice we point to
an example in `examples/` in which it is made for a physical system and the
result is checked against an exact one. We assume that LSQuant is built and on
`PYTHONPATH`, as in [Getting started](../tutorials/getting-started.md).

## The three steps of a calculation

Every LSQuant calculation has three steps:

1. We write a **run file**, the JSON input of one calculation, much as an
   `.fdf` file is the input of a SIESTA run. It names the model, the
   representation, the quantity and the convergence parameters.
2. `python -m lsquant run` computes the **Chebyshev moments** of the quantity
   and writes them to disk. These numbers come from the Hamiltonian, and every
   curve is built from them. This is the expensive step.
3. A **reconstruction** turns the moments into a curve, such as a DOS or a
   conductivity σ(E), with the Jackson kernel and a chosen energy broadening.
   `run` does this at the end, and `python -m lsquant reconstruct` repeats it later from the saved
   moments, without applying the Hamiltonian again.

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

This file reads `haldane.HAM.CSR`, the Hamiltonian of the Haldane model of
[example 03](../../examples/03_haldane_hall/README.md), stored in the
repository with the reference outputs. The DOS is computed with an exact trace
on the default spectral window. The **spectral window** is an energy range that
encloses every eigenvalue; here it is a Gershgorin bound, computed from the
rows of the Hamiltonian. Every key that we leave out takes its default value,
and the [run-file reference](../reference/run-file.md) lists them all.
<!-- source: python/lsquant/runfile.py:266-274, :303 (default window "gershgorin" for real), :1804-1814; python/lsquant/moments.py:376-385 (rho = max row sum, padded 10 %) -->

The requested 10 meV is finer than these 256 moments can resolve. The
Gershgorin bound of this model is ±4.3 eV, so the curve comes out with a
broadening of about π × 4.3 eV / 256 ≈ 50 meV, as explained in
[Number of moments and energy resolution](#number-of-moments-and-energy-resolution).
A 10 meV curve would need about 1350 moments, more than the 512 states of
this model, and would then show its discrete levels (see the upper bound on
$M$ in the same section).
<!-- source: derived from goldens/stage0/haldane/operators/haldane.HAM.CSR (512 rows, 4608 entries: per row three hoppings |t1| = 1 eV and six |t2| = 0.15 eV, as in examples/03_haldane_hall/README.md), so rho = 3.9 eV and the padded window is ±4.29 eV; pi*4.29/256 = 0.053 eV; pi*4.29/0.010 = 1348. Clamp: python/lsquant/reconstruct.py:89–96, :136–145 -->

Relative paths in a run file are tried against the LSQuant source tree first
([path resolution](../reference/run-file.md#path-resolution)). For our own
data we therefore give absolute paths.

Before we spend compute time on a run file, we let LSQuant validate it and fill
in the defaults:

```bash
python -m lsquant inspect my_run.json
```

A run file that fails `inspect` also fails `run`. Operator names outside the
standard set, however, pass `inspect` and fail in `run` (see the
[run-file reference](../reference/run-file.md#word-and-words)).

## The model: SIESTA, Wannier90 or tight-binding

The key `model.kind` selects where the Hamiltonian comes from:

| Starting point | `model.kind` | Example |
|---|---|---|
| A SIESTA calculation with spin–orbit coupling | `siesta` | [07](../../examples/07_spin_hall_bi_bilayer/README.md), [08](../../examples/08_spin_hall_pdse2_monolayer/README.md) |
| A Wannier90 `_hr.dat` file | `wannier` | [14](../../examples/14_wannier/README.md) |
| A tight-binding model we write in Python | `generator` (the model generators of examples 01–06) or `csr` (sparse-matrix files we write ourselves) | [01](../../examples/01_dos_1d_chain/README.md), [03](../../examples/03_haldane_hall/README.md) |

For a SIESTA model, [Use a SIESTA calculation](use-siesta.md) lists what the
SIESTA run must save and which settings LSQuant accepts.
[Example 14](../../examples/14_wannier/README.md) reads the
Su–Schrieffer–Heeger chain from a hand-written `_hr.dat` file and reproduces its
two analytic bands.

### Atomic-orbital bases and the overlap matrix

The atomic-orbital basis of SIESTA is non-orthogonal, with an overlap matrix
$S$. For such a basis we set `metric.policy` to `"nonorthogonal"`. The
calculation then works either with $S^{-1/2} H S^{-1/2}$ (the Löwdin route)
or with $S^{-1} H$ (the generalised route).

[Example 11](../../examples/11_nonorthogonal_toy/README.md) treats a
non-orthogonal chain on both routes and reproduces the exact generalised
eigenvalue problem to 10⁻¹³.
[Non-orthogonal bases](../explanation/nonorthogonal-bases.md) explains why $S$
changes the velocities and currents.

## The quantity

The quantity to compute, the **observable**, is written as a trace of
operators and functions of $H$ in the key `word`. Its Chebyshev moments are
turned into the physical curve by a reconstruction formula, the key
`reconstruction.formulas`. The table gives, for each quantity, the trace, the
formula and the unit of the curve. All curves are at zero temperature, with a
step-function Fermi distribution.
<!-- source: zero temperature from docs/explanation/kubo-formulas.md:13-14; units from docs/reference/output-files.md:208-212, :237-240 and docs/explanation/kubo-formulas.md:93-95, :147-153 -->

In the traces, `VX` and `VY` are the velocity operators $v_x$ and $v_y$,
`JXSZ` is the spin current $J^z_x$, `dG-(E)/dE` is the energy derivative of the
advanced Green's function, `U(t)` is the time-evolution operator and `^+` the
adjoint. In `P:PX`, the operator `PX` of example 04 is the projector
$P_x = (1 + S_x)/2$ onto spin along $+x$; the
[reference on observables](../reference/observable-words.md) explains this
syntax and lists every formula that each form of `word` allows.

| Quantity | Word | Formula | Unit of the curve | Example |
|---|---|---|---|---|
| DOS | `Tr[delta(E-H)]` | `dos` | 1/eV, summed over all N states: N ρ(E) | [01](../../examples/01_dos_1d_chain/README.md) |
| DOS projected on a region or an orbital set | `Tr[A delta(E-H)]` | `dos` | as the DOS | [15](../../examples/15_spatial_operators/README.md) |
| σ_xx (Kubo–Greenwood) | `Tr[VX delta(E-H) VX delta(E-H)]` | `greenwood` | not stated by the code | [05](../../examples/05_charge_transport_square/README.md), [10](../../examples/10_bases_compare/README.md) |
| σ_xy (Kubo–Bastin) | `Tr[VX delta(E-H) VY dG-(E)/dE]` | `bastin` | summed over the sample, e = ħ = 1; times 2π/A gives e²/h in 2D | [03](../../examples/03_haldane_hall/README.md) |
| σ_xy split into Fermi-sea and Fermi-surface terms | as above | `sea`, `surf` | as `bastin` | [13](../../examples/13_bastin_sea_surface/README.md) |
| Spin Hall conductivity | `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` | `shc` | e²/h, spin as Pauli σ_z (halve for the ħ/2 convention) | [07](../../examples/07_spin_hall_bi_bilayer/README.md) |
| Spreading of a wave packet (mean-square displacement) | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | `time_correlation` | not stated by the code; time in fs | [02](../../examples/02_localization_1d_chain/README.md), [15](../../examples/15_spatial_operators/README.md) |
| Spin dynamics | `Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]` and related | `time_correlation` | not stated by the code; time in fs | [04](../../examples/04_spin_precession_magnetic_chain/README.md), [06](../../examples/06_spin_relaxation_square/README.md) |

The command `words` lists the forms of `word` that LSQuant accepts, and
`words describe` shows how one of them is read:

```bash
python -m lsquant words
python -m lsquant words describe "Tr[VX delta(E-H) VY dG-(E)/dE]"
```

With an ordinary `csr` or `wannier` model, set `reconstruction.num_cells`
and `reconstruction.cell_area` when requesting `shc`. Validation requires
both keys before moments are computed; the run writes the curve directly.
See the [run-file reference](../reference/run-file.md#reconstruction).

## The representation: supercell or k-mesh

The key `basis.name` sets how $H$ is represented:

- `real`: a supercell in real space, for disorder and defects.
- `bloch`: one block per k-point of the unit cell, for clean crystals.
- `eigen`: the exact eigenstates of each block, for small reference checks.

For a clean periodic system the three give the same moments, and `bloch` or
`eigen` are much cheaper. We use `real` when disorder or a defect breaks the
translation symmetry. [Example 10](../../examples/10_bases_compare/README.md)
computes one graphene Hamiltonian in the three bases; their moments agree to
10⁻¹⁴.

### Disorder

For Wannier and SIESTA models we set the Anderson disorder strength
`basis.disorder.W`, in eV, and its seed `basis.disorder.seed`. The on-site
energies are drawn uniformly from $[-W/2, W/2]$. For SIESTA models with
spin–orbit coupling, `basis.disorder.kind` chooses a spin-independent
potential (`"scalar"`, the default) or a random Zeeman field (`"magnetic"`).
Sparse-matrix models take their disorder from the files themselves, and the
generators build their own.
<!-- source: python/lsquant/bases.py:472-494 (uniform in [-W/2, W/2]); python/lsquant/disorder_bases.py:32-34; docs/reference/run-file.md:189-194 -->

[Example 02](../../examples/02_localization_1d_chain/README.md) follows the
Anderson localisation of a one-dimensional chain, and
[example 05](../../examples/05_charge_transport_square/README.md) the resistance
of a disordered square lattice.
[Example 12](../../examples/12_disorder_nonorthogonal/README.md) adds disorder
to a SIESTA basis with an overlap matrix and explains why the default
`metric.tol` must then be relaxed.

## Convergence parameters

We converge the parameters in the order of the sections below. The window
fixes how many moments a given resolution costs, so we set it first; the noise
of the random vectors and the disorder average come last.

### Spectral window

For a given number of moments, the energy resolution scales with the width of
the spectral window, so a narrower window makes every curve cheaper. For a
SIESTA model, the `projection` section keeps only the states in an energy
window, so the expansion covers a much narrower spectrum.
[Example 16](../../examples/16_projection_downfold/README.md) projects a Bi
bilayer Hamiltonian onto [−2, 2] eV around the Fermi level and recovers the
exact DOS of the kept bands.

### Number of moments and energy resolution

The number of moments $M$, the key `expansion.moments`, sets the energy
resolution of the curve, much as the smearing does for a DFT DOS. The finest
broadening we can resolve is $\pi a/M$, where $a$ is the half-width of the
spectral window. We therefore choose $M$ from the resolution we need.

A broadening finer than $\pi a/M$ is replaced by $\pi a/M$ with a
`RuntimeWarning`; the effective estimate is recorded in the curve metadata or run receipt.
A wide window, as in most DFT Hamiltonians, therefore needs many moments: a
20 eV window and a 10 meV broadening need $M \ge 3142$.
[The kernel polynomial method](../explanation/kpm.md) explains the kernel and
the broadening.
<!-- source: python/lsquant/reconstruct.py:89–96, :136–145 (silent clamp; numbers from docs/explanation/kpm.md) -->

$M$ also has an upper bound in a finite system. A supercell with $N$ states
has $N$ discrete levels, and the curve resolves them as separate peaks as soon
as more than $N$ moments contribute. [Getting started](../tutorials/getting-started.md)
shows this for a ring of 256 sites.
<!-- source: docs/tutorials/getting-started.md:134-141 -->

[Example 09](../../examples/09_inspect_reconstruct/README.md) reconstructs one
set of moments at several broadenings and checks each curve against the exact
DOS.

### Supercell or k-mesh and finite-size error

For Wannier and SIESTA models, the key `basis.supercell` sets the supercell
size in the `real` representation and the k-point mesh in `bloch` and `eigen`.
It sets the finite-size error, as the k-mesh does in a DFT calculation. We
enlarge it, at fixed broadening, until the curve changes by less than the
accuracy we need. [Example 10](../../examples/10_bases_compare/README.md) uses
a 16 × 16 supercell or k-mesh, and [example 14](../../examples/14_wannier/README.md)
a 256-point k-mesh, `basis.supercell = [256, 1, 1]`.

### Random vectors and statistical noise

With `expansion.trace = "exact"`, the trace is summed over every basis state.
For large real-space systems we use `"random"`, which estimates the trace with
a few random-phase vectors; their number $R$, the key `expansion.states`, sets
the noise of the curve. In a system of dimension $N$, the relative error falls
as $1/\sqrt{RN}$, and it is larger where the DOS is small. Large systems
therefore need few vectors.
[Example 06](../../examples/06_spin_relaxation_square/README.md) uses
$R = 20$ random vectors on an $L = 50$ lattice.

### Disorder realisations

To average over disorder, we run one run file per seed and average the
curves; a single run computes one realisation.
[Example 02](../../examples/02_localization_1d_chain/README.md) averages the
mean-square displacement over disorder samples, and
[example 12](../../examples/12_disorder_nonorthogonal/README.md) compares two
realisations.

### Cost

The run time is counted in sparse matrix–vector products with $H$. A DOS
needs $R\,M$ products with $R$ random vectors, and $N\,M$ with an exact
trace in a real-space supercell.
A conductivity needs $2RM$ products and $R\,M^2$ inner products of length
$N$, and for large $M$ the inner products dominate. The memory grows with the
number of vectors of length $N$ held at once;
[Run on a cluster](run-on-a-cluster.md#estimate-and-limit-the-memory) gives
the memory of each setting.
<!-- source: docs/explanation/kpm.md, sections "Exact and random-phase traces" and "Cost" (moments.py:287-318) -->

## Calculation

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

## Inspection and reconstruction

To see what produced a directory and to draw a new curve from its moments, we
use:

```bash
python -m lsquant inspect out/my_run
python -m lsquant reconstruct out/my_run --formula dos --broadening-meV 120
```

`reconstruct` gives a new curve in seconds, with another broadening or with
another formula that the same `word` allows. A convergence test in the
broadening therefore costs one `run` and several reconstructions.

## Calculations without a run file

A run file covers the standard calculations. To change one step of a
calculation, we drive LSQuant from Python, as in
[example 17](../../examples/17_python_api/README.md). An AI assistant can also
run LSQuant through its MCP server:
[example 18](../../examples/18_mcp_session/README.md) shows a session, and
[Connect an AI assistant](mcp-server.md) shows how to set it up.
