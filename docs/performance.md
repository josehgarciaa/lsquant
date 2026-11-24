\page performance Performance: thread placement, state blocks, atom blocks, checkpoints

This page covers the performance features of the native build profile (see \ref build) and
the checkpointing available in both profiles. None of them changes what a run file computes:
- the golden profile keeps the serial reference path;
- the native profile reorders the arithmetic and agrees with it to rounding.

## Thread placement

When the process is the `lsquant` command (`lsquant run` or `python -m lsquant run`) and it
has its own cores (a CPU set from the batch system, smaller than the machine), the package sets
`OMP_PROC_BIND=spread` and `OMP_PLACES=cores` at import. Without a CPU set it does not bind:
OpenMP numbers the places from the first core for every process, so two commands on one node
would bind to the same cores. `LSQUANT_BIND_THREADS=1` binds anyway (a node you have to
yourself) and `LSQUANT_BIND_THREADS=0` never binds. A plain `import lsquant`
in a script or notebook does not set them: with binding, OpenMP pins the importing process to
one core. If `OMP_NUM_THREADS` is unset, the package sets it to the number of physical cores in
the process's CPU affinity (`lsquant.runtime`). Each variable is set only when it is unset, so a
value from the caller wins.

In the native profile two kernel settings reach the engine once, at import:
- `LSQUANT_SPMV_MIN_NNZ`: below this number of non-zeros a matrix-vector product runs on one
  thread; default 16384;
- `LSQUANT_KERNEL_BATCH`: the batch size of the batched kernels; default 0, automatic.

The engine itself never reads the environment.

## State blocks

`expansion.state_block` processes several trace vectors at once in the native profile, which
turns matrix-vector products into matrix-matrix products:

```json
"expansion": {"moments": 1024, "trace": "random", "states": 32, "state_block": "auto"}
```

`"auto"` picks S = min(states, 16, ⌊0.6 · MemAvailable / ((right_window + left_window) · N · 16 B)⌋),
never below 1 (`lsquant.moments.choose_state_block`). `left_window` and `right_window` set how
many Chebyshev vectors of each side stay in memory. The golden profile always uses S = 1. The
resolved value and the reason are written back into the run's `run.resolved.json`
(`state_block_rule`).

## Atom-block operators

```json
"numerics": {"operator": "atom_block"},
"model": {"kind": "csr", "dir": "...", "label": "...", "positions": "positions.xyz", "orbitals_per_atom": 2}
```

In the native profile, `numerics.operator = "atom_block"` stores the Hamiltonian as a
dictionary of dense atom blocks. Before that, the atoms are reordered:
- atoms grouped in spatial domains (`numerics.domain_atoms`);
- Morton order inside and across domains;
- interior atoms before boundary atoms (`lsquant.spatial`).

The reordering is a permutation, so traces are unchanged. The example below checks this
against the CSR operator. The layout needs a real-space basis, an orthogonal metric and the
atom positions. A crystal can also be given as hoppings H(R) (`model.hr_file`) with a
`model.supercell`.

## Checkpoints

`"expansion": {"checkpoint": true}` makes the engine write the moment tensor to
`<out>/checkpoint_<hash>.lsqckpt` after every fully accumulated trace state:
- the file is versioned and written atomically (temporary file, `fsync`, `rename`);
- a run restarted with the same run file resumes after the last state in the file;
- the result is bit-identical to an uninterrupted run;
- the file is removed when the run completes.

This works in both profiles.

## Optional backends

These need the native profile and are off by default:

| Option | What it adds | Run file |
|---|---|---|
| `-DLSQUANT_WITH_CUDA=ON` | the CUDA atom-block kernel | `numerics.execution = "cuda"` with `numerics.operator = "atom_block"` |
| `-DLSQUANT_WITH_GINKGO=ON` | Ginkgo sparse kernels (formats `csr`, `sellp`, `ell`, `hybrid`, `coo`; executors `omp`, `cuda`) | `numerics.backend = {"name": "ginkgo", "executor": "omp", "format": "csr"}` |
| `-DLSQUANT_WITH_MPI=ON` | the distributed operator (`DistributedOperator`) | `numerics.execution = "mpi"` is refused by the run file: the MPI path is not wired into run files yet |

Configuring with any of these and `LSQUANT_BUILD_PROFILE=golden` stops with an error. The
golden profile is OpenMP-only by construction.

## Example: atom blocks and checkpoints from run files

`examples/performance/atom_block.py` builds a square lattice in Python and runs the same
density-of-states word with the CSR operator, with the atom-block operator, and with the
atom-block operator plus checkpointing. On 64 × 64 sites with an exact trace the three agree;
on 1024 × 1024 sites it reports the run times. `examples/performance/README.md` gives the
results.
