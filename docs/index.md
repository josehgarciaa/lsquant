# LSQuant documentation

LSQuant computes spectral and transport properties of large tight-binding and
first-principles Hamiltonians with the kernel polynomial method, at a cost that
grows linearly with the number of orbitals. The pages are grouped by what you
want to do.

## Tutorials: learn by doing

- [Getting started](tutorials/getting-started.md): build LSQuant and compute the
  density of states of a chain.
- [The Hall conductivity of the Haldane model](tutorials/haldane-hall.md).
- [The spin Hall conductivity of a Bi(111) bilayer from SIESTA](tutorials/spin-hall-siesta.md).
- The examples 01–08 in [`examples/`](../examples/) are full physics lessons:
  the DOS of a chain, Anderson localisation, the Haldane model, spin
  precession, the resistance of a disordered metal, spin relaxation, and the
  spin Hall conductivity of two SIESTA materials.

## How-to guides: get a task done

- [Running LSQuant](how-to/run-lsquant.md): the choices of a calculation, each
  with a use case from examples 09–18.
- [Install LSQuant](how-to/install.md)
- [Use a SIESTA calculation](how-to/use-siesta.md)
- [Run on a cluster](how-to/run-on-a-cluster.md)
- [Check a build against the reference outputs](how-to/testing.md)
- [Connect an AI assistant through MCP](how-to/mcp-server.md)

## Reference: look something up

- [Run file](reference/run-file.md)
- [Command line](reference/cli.md)
- [Observable words](reference/observable-words.md)
- [Output files](reference/output-files.md)
- [Bases and metrics](reference/bases-and-metrics.md)
- [Build presets](reference/build-presets.md)
- [MCP tools](reference/mcp-tools.md)
- [Python API](reference/python-api.md)
- [Core module `lsquant._core`](reference/core-api.md)

## Explanation: understand how it works

- [How LSQuant works](explanation/unifying-idea.md): observables as traces of
  operator words.
- [The kernel polynomial method](explanation/kpm.md)
- [Kubo formulas](explanation/kubo-formulas.md)
- [Time evolution](explanation/time-evolution.md)
- [Non-orthogonal bases](explanation/nonorthogonal-bases.md)
- [Architecture](explanation/architecture.md)
