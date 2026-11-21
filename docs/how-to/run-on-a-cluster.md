# Run on a cluster

This guide shows you how to set the thread count, choose a build, size the
memory and write a batch script for LSQuant runs on a cluster or a shared
compute node. It assumes you can already run a calculation; if not, read
[How we run LSQuant](run-lsquant.md) first.

## Set the number of threads

LSQuant parallelises over OpenMP threads inside one process. Set the thread
count with `OMP_NUM_THREADS` before Python starts.

When Python imports `lsquant`, it sets three variables, each one only if you
have left it unset:
<!-- source: python/lsquant/runtime.py -->

| Variable | Value set by `import lsquant` |
|---|---|
| `OMP_NUM_THREADS` | The number of physical cores in the CPU affinity of the process (cpusets included). |
| `OMP_PROC_BIND` | `spread` |
| `OMP_PLACES` | `cores` |

The default suits a node that you have to yourself. On a node with a CPU
allocation (cgroups or cpusets, as most batch systems use), the affinity
already contains only your cores. On a shared node without an allocation, the
default is every physical core of the machine, which is almost never what you
want.

Every run records the values it used in `run.resolved.json`, under
`runtime`, so you can check afterwards what a run did.

### On a node you have to yourself

Use one thread per physical core of your allocation:

```bash
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK   # or the number of cores you requested
```

On a node with two sockets, a large system can gain from memory placement.
Run on one socket, or interleave the memory over both:

```bash
numactl --cpunodebind=0 --membind=0 python -m lsquant run run.json --out out/run
numactl --interleave=all python -m lsquant run run.json --out out/run
```

[Performance and scaling](../performance/README.md) shows the thread scaling of
the sparse product on a 256-core node.

### On a shared or busy node

On a node where other jobs also run, set `OMP_NUM_THREADS=1`, or let idle
threads sleep with `OMP_WAIT_POLICY=passive`:

```bash
export OMP_NUM_THREADS=1
# or, to keep several threads:
export OMP_WAIT_POLICY=passive
```

An OpenMP parallel region ends only when its last thread finishes. On an
oversubscribed node the operating system deschedules some threads, and the
others wait for them. With the default wait policy, waiting threads also spin
and take CPU time from the threads they wait for.

The following times were measured on a shared 96-core node with the `golden`
build, with the same command run minutes apart:

| Run | `OMP_NUM_THREADS=4` | `OMP_NUM_THREADS=1` |
|---|---|---|
| Example 11, orthogonal chain (N = 128, M = 256, exact trace) | 145.6 s | 5.6 s |
| Example 01, DOS of three chains (M = 512) | 357.7 s; 2892 s when three other 4-thread jobs ran (stopped after 48 min) | 148.7 s |

For these small problems, 4 threads were 10 to 30 times slower than one.
`OMP_WAIT_POLICY=passive` was not measured on that node. The binding that LSQuant
sets by default (`OMP_PROC_BIND=spread`) also pins threads to cores that other
jobs may be using; set `OMP_PROC_BIND=false` if you share cores with other
jobs.

### Choose a build for small matrices

| | `golden` | `native` |
|---|---|---|
| Sparse matrix–vector product | Every product opens an OpenMP parallel region, whatever the size of the matrix. | Products with fewer non-zeros than `LSQUANT_SPMV_MIN_NNZ` (default 16384) run on one thread. |
| Results and thread count | The moments do not depend on the thread count. The golden-file check tests this at 1, 8 and 64 threads. | Results can differ from `golden` in the last digits. |
<!-- source: lsquant/core/algebra/sparse_matrix.cpp:148 (golden), :137 (native threshold) -->

For small matrices, extra threads in a `golden` build add only
synchronisation cost. Use one thread, or a `native` build.

In a `native` build, `lsquant` reads two tuning variables when it is imported:
<!-- source: python/lsquant/runtime.py:50-54 -->

| Variable | Default | Meaning |
|---|---|---|
| `LSQUANT_SPMV_MIN_NNZ` | `16384` | Smallest number of non-zeros for which the sparse product runs in parallel. |
| `LSQUANT_KERNEL_BATCH` | `0` | Batch width of the sparse matrix–matrix kernel; 0 lets the engine choose. |

The `golden` build ignores both.

## Use the MPI and CUDA builds

- **MPI.** The `native-mpi` and `full` presets build the MPI distributed
  operator, but no run-file path uses it yet:
  `numerics.execution = "mpi"` stops with `is not wired yet`.
  <!-- source: python/lsquant/runfile.py:1896-1899 -->
  Run one process per node with
  OpenMP threads, and spread independent runs (seeds, energies, k-meshes)
  over nodes with your batch system. To average over disorder seeds, submit
  one job per seed, or a job array with one run file per seed, rather than
  one job with many threads.
- **CUDA.** The `native-cuda` preset builds the GPU kernel for compute
  capability 9.0 only; see the [known limitation](../reference/build-presets.md).
  A run uses the GPU only with a real-space `csr` model and
  `numerics.operator = "atom_block"`, `numerics.execution = "cuda"` (see
  [the atom-block layout](#use-the-atom-block-layout-for-large-real-space-models)
  and [`numerics`](../reference/run-file.md#numerics)). SIESTA and Wannier
  models run on the CPU.

## Use the atom-block layout for large real-space models

For a large real-space model, a `native` build can store the Hamiltonian in an
atom-block layout instead of a plain CSR matrix. LSQuant groups the orbitals
of each atom into one block, gathers neighbouring atoms into spatial domains,
and orders atoms and domains along a space-filling (Morton) curve. Vectors
then stay in cache during the sparse product. The CUDA execution uses the
same layout. The moments are traces, so they are invariant under the
reordering.

To turn the layout on, add these keys to the run file:

```json
"numerics": {"operator": "atom_block"},
"model": {
  "kind": "csr",
  "positions": "/path/to/positions.xyz",
  "orbitals_per_atom": 2
}
```

The keys mean:

- `model.positions`: the position of every orbital, as a path or as an array
  of shape (N, 3). A file holds one `x y z` line per orbital, or uses the
  `.xyz` format: a count line, a comment line, then `tag x y z` lines.
- `model.orbitals_per_atom`: the number of consecutive orbitals that form one
  atom. For atoms with different numbers of orbitals, give `model.atom_map`
  instead.
- `numerics.domain_atoms` (optional): the number of atoms per spatial domain.
  LSQuant chooses a value if you leave it out.

The layout accepts only a `csr` model in the `real` basis with the orthogonal
metric. The `golden` build refuses it with
`numerics.operator='atom_block' is only available in the native build profile`.
<!-- source: python/lsquant/spatial.py:1-9 (permutation, Morton order), :33-60 (read_positions); python/lsquant/runfile.py:690-734 (validation), :1768-1900 (execution) -->

## Estimate and limit the memory

The engine stores Chebyshev vectors of dimension N, the dimension of the
Hamiltonian. One vector of complex doubles takes 16 N bytes: 16 MB for
N = 10⁶. The expansion keys of the run file set how many vectors are held at
once.

| Key | Default | Effect on memory | Cost |
|---|---|---|---|
| `expansion.batch` | `moments` | In the `golden` build, a rank-2 word (σ_xx, σ_xy) holds two blocks of `batch` vectors, about 2 × batch × 16 N bytes per random vector. | With `batch` < M, the right-hand recursion is repeated for every block of the left side, so the time grows about M / batch times. The moments do not change. |
| `expansion.state_block` | `1` | Random or unit vectors propagated together. Memory grows in proportion. `"auto"` picks a value from the available memory in the `native` build, and 1 in `golden`. | Larger blocks use the cache better. |
| `expansion.left_window`, `expansion.right_window` | `32`, `16` | `native` build, rank-2 words: about (left_window + right_window) × 16 N bytes per random vector; `right_window: "all"` stores all M vectors of the right side. | Not used by the `golden` build. |
| `expansion.checkpoint` | `false` | No effect on memory. `csr` models only: saves the accumulated moments after every random or unit vector, so that a run killed by the wall-time limit resumes. | One file write per vector. |

Example: σ_xx with M = 1024 on N = 10⁶ in the `golden` build holds about
2 × 1024 × 16 MB = 32 GB with the default `batch`. With `batch: 128` it holds
4 GB and takes up to 8 times longer.

SIESTA runs also hold dense k-blocks and, for non-orthogonal real-space
runs with the exact applier, dense factorisations of S, up to
`metric.max_dim` (default 4096). `model.params.bloch_memory_mb` (default 256)
limits the work space of the Bloch sums in clean SIESTA runs. See
[Use a SIESTA calculation](use-siesta.md).

### Resume a killed run

With `expansion.checkpoint: true`, a `csr` run writes
`checkpoint_<hash>.lsqckpt` into the output directory after every random or
unit vector. Run the same command again with the same `--out` directory, and
the engine continues after the last finished vector. The file is deleted when
the word finishes.

> **Known limitation: the checkpoint does not record the model or the window.**
> The file name depends only on the word (`python/lsquant/runfile.py:2148`),
> and on loading the engine checks only the rank, the dimensions, the number of
> trace vectors, the state label and the word
> (`lsquant/core/chebyshev/checkpoint.cpp:144–148`). If you change the
> Hamiltonian, the window or the seed and run again into a directory that still
> holds a checkpoint from a killed run, the new run continues from the old
> moments. Delete stale `checkpoint_*.lsqckpt` files, or use a new output
> directory, whenever you change the run file.

## Write a batch script

This Slurm script runs one calculation on one node. Adapt the directives to
your batch system.

```bash
#!/bin/bash
#SBATCH --job-name=lsq-shc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00

set -euo pipefail

# Environment: the same one that built LSQuant.
conda activate lsq                       # or: source .venv/bin/activate
LSQ=$HOME/lsquant                         # the source tree
export PYTHONPATH="$LSQ/build/native/lsquant_core:$LSQ/python"

# Threads: one per allocated core; idle threads sleep.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_WAIT_POLICY=passive
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

cd "$SLURM_SUBMIT_DIR"
python -m lsquant inspect run_shc.json > /dev/null   # stop early if the run file is invalid
python -m lsquant run run_shc.json --out "out/shc_$SLURM_JOB_ID"
```

Notes on the script:

- Run `inspect` before `run`. It validates the run file in seconds, so you
  find a typo before the queue wait. Some errors still pass `inspect` and stop
  `run`; see the known limitations in the
  [run-file reference](../reference/run-file.md).
- Relative `model.dir` paths are tried against the LSQuant source tree first,
  then against the directory of the run file. Use absolute paths in batch
  jobs. `--out` is relative to the working directory.
- Set `OPENBLAS_NUM_THREADS=1` and `MKL_NUM_THREADS=1` so that NumPy and SciPy
  keep to one thread each and leave the cores to the OpenMP threads.
- For a `native` build, the job must run on the CPU type that built it. See
  [Install LSQuant](install.md#build-native).
- For disorder averages, use one job per seed, as described in
  [Use the MPI and CUDA builds](#use-the-mpi-and-cuda-builds).

## See also

- [Install LSQuant](install.md)
- [Check a build against the golden files](testing.md)
- [Run-file reference: `expansion`](../reference/run-file.md#expansion)
- [Build presets](../reference/build-presets.md)
