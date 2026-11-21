# LSQuant documentation

LSQuant computes the density of states, the conductivity, the Hall and spin
Hall conductivities and spin dynamics of systems too large to diagonalise. It
uses the kernel polynomial method, at a cost that grows linearly with the
number of orbitals. The pages below follow the order in which you meet them:
install, a first run, your own material, disorder, then production runs on a
cluster.

## Installation and build check

- [Install LSQuant](how-to/install.md): the dependencies with or without root
  access, the build presets and `PYTHONPATH`.
- [Getting started](tutorials/getting-started.md): build LSQuant, compute the
  density of states of a chain in about 2.5 minutes on one core, and plot it
  against the exact result.
- [Check a build against the reference outputs](how-to/testing.md).

## Your own material

- [Running LSQuant](how-to/run-lsquant.md): the choices of a calculation
  (model, basis, quantity, number of moments), each with a use case from
  examples 09–18.
- [Use a SIESTA calculation](how-to/use-siesta.md): prepare a SIESTA
  calculation with spin–orbit coupling and write its run file.
- [The spin Hall conductivity of a Bi(111) bilayer from SIESTA](tutorials/spin-hall-siesta.md):
  the workflow of Getting started applied to a first-principles Hamiltonian.
- [The Hall conductivity of the Haldane model](tutorials/haldane-hall.md):
  σ_xy locked to e²/h inside the gap.
- [Example 14](../examples/14_wannier/README.md): the DOS from a Wannier90
  `_hr.dat` file.

## Physics examples

Examples 01–08 in [`examples/`](../examples/) each address one physical
question: the DOS of a chain, Anderson localisation, the quantised Hall
conductivity of the Haldane model, spin precession, the resistance of a
disordered metal, spin relaxation, and the spin Hall conductivity of two
SIESTA materials.

## Disorder and defects

- [Adding disorder](how-to/run-lsquant.md#disorder) in Running LSQuant,
  with [example 12](../examples/12_disorder_nonorthogonal/README.md): Anderson
  disorder in a SIESTA model with an overlap matrix.
- [Scalar and magnetic disorder](../examples/disorder_kinds/README.md): the
  Kramers splitting of a Bi(111) bilayer in a 3×3 supercell under scalar,
  magnetic and independent disorder.

## Production runs on a cluster

- [Run on a cluster](how-to/run-on-a-cluster.md): threads, memory, resuming a
  killed run and a Slurm batch script.

## Reference

- [Run file](reference/run-file.md), including the [basis](reference/run-file.md#basis)
  and [metric](reference/run-file.md#metric) keys
- [Command line](reference/cli.md)
- [Observables (the `word` key)](reference/observable-words.md): how to
  write the quantity to compute as a trace
- [Output files](reference/output-files.md)
- [Build presets](reference/build-presets.md)
- [Connect an AI assistant through MCP](how-to/mcp-server.md) and the
  [MCP tools](reference/mcp-tools.md)
- [Python API](reference/python-api.md)

## Background on the method

- [The kernel polynomial method](explanation/kpm.md): the Chebyshev moments,
  the Jackson kernel, exact and random-phase traces, and the cost.
- [Kubo formulas](explanation/kubo-formulas.md)
- [Time evolution](explanation/time-evolution.md)
- [Non-orthogonal bases](explanation/nonorthogonal-bases.md)

## For developers

- [Architecture](explanation/architecture.md)
- [Core module `lsquant._core`](reference/core-api.md)
