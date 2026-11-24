# Run on a cluster

To run LSQuant on a cluster or a shared compute node, set the thread count,
choose a build, estimate the memory and write a batch script, as described
below. Start from a run file that already works on a small case; if you have
none, read [Running LSQuant](run-lsquant.md) first.

## Set the number of threads

LSQuant runs in parallel on OpenMP threads inside one process. Set the thread
count with `OMP_NUM_THREADS` before Python starts.

When Python imports `lsquant`, it sets three variables, each one only if you
have left it unset:
<!-- source: python/lsquant/runtime.py -->

| Variable | Value set by `import lsquant` |
|---|---|
| `OMP_NUM_THREADS` | The number of physical cores in the CPU affinity of the process (cpusets included). |
| `OMP_PROC_BIND` | `spread`, only for the `lsquant` command with its own cores (see below) |
| `OMP_PLACES` | `cores`, under the same condition |

The default suits a node that you have to yourself. On a node with a CPU
allocation (cgroups or cpusets, as most batch systems use), the affinity
already contains only your cores. On a shared node without an allocation, the
default is every physical core of the machine, so set the thread count
yourself there.

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
the sparse matrix–vector product on a 256-core node. For the DOS of a graphene
model with N = 2 × 10⁶ and M = 512, the whole recursion takes 17.9 s on one
thread, 3.9 s on 32 and 1.8 s on 128. The gain stops near 10 times because
the inner products of the moments run on one thread.
<!-- source: docs/performance/README.md, "Result 1 — thread scaling" table (full recursion 17930, 3924, 1786 ms) and the note on serial vdot/AXPY -->

### On a shared or busy node

On a node where other jobs also run, set `OMP_NUM_THREADS=1`, or let idle
threads sleep with `OMP_WAIT_POLICY=passive`:

```bash
export OMP_NUM_THREADS=1
# or, to keep several threads:
export OMP_WAIT_POLICY=passive
```

Each parallel step ends only when its last thread finishes. On an
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
`OMP_WAIT_POLICY=passive` was not measured on that node.

The `lsquant` command binds its threads to cores (`OMP_PROC_BIND=spread`,
`OMP_PLACES=cores`) only when the batch system gave the job its own cores: a
CPU set smaller than the machine. On a node without one, every process would
bind to the same cores, counted from the first, and two runs side by side
would share them while the other cores idle. There the command does not bind.
Set `LSQUANT_BIND_THREADS=1` to bind anyway, on a node you have to yourself,
or `LSQUANT_BIND_THREADS=0` never to bind. Python scripts that import
`lsquant` are not bound; set `OMP_PROC_BIND` yourself if you want binding
there.

### Choose a build for small matrices

| | `golden` | `native` |
|---|---|---|
| Sparse matrix–vector product | Every product runs on all threads, whatever the size of the matrix. | Products with fewer non-zeros than `LSQUANT_SPMV_MIN_NNZ` (default 16384) run on one thread. |
| Results and thread count | The moments do not depend on the thread count. The reference check of the build tests this at 1, 8 and 64 threads. | Results can differ from `golden` in the last digits. |
<!-- source: lsquant/core/algebra/sparse_matrix.cpp:148 (golden), :137 (native threshold) -->

For small matrices, extra threads in a `golden` build only add waiting time.
Use one thread, or a `native` build.

In a `native` build, `lsquant` reads two tuning variables when it is imported:
<!-- source: python/lsquant/runtime.py:50-54; lsquant/core/algebra/work_partition.hpp:27-29, 48-53 (min(setting, vectors available); 0 = half a 512 KiB L2) -->

| Variable | Default | Meaning |
|---|---|---|
| `LSQUANT_SPMV_MIN_NNZ` | `16384` | Smallest number of non-zeros for which the sparse product runs in parallel. |
| `LSQUANT_KERNEL_BATCH` | `0` | Largest number of vectors multiplied together in one sparse matrix–matrix product; 0 lets LSQuant choose from the cache size. |

The `golden` build ignores both.

## Use the MPI and CUDA builds

The `native-mpi` and `full` presets build the distributed (MPI) operator,
but run files cannot use it yet: `numerics.execution = "mpi"` stops with
`is not wired yet`.
<!-- source: python/lsquant/runfile.py:1896-1899 -->
Run one process per node with threads, and spread independent runs (seeds,
energies, k-meshes) over nodes with your batch system. To average over
disorder seeds, submit one job per seed, or a job array with one run file per
seed, rather than one job with many threads.

The `native-cuda` preset builds the GPU code for compute capability 9.0 by
default; set `LSQUANT_CUDA_ARCHITECTURES` for other GPUs (see
[CUDA architecture](../reference/build-presets.md#presets)). A run uses the GPU
only with a real-space `csr` model and `numerics.operator = "atom_block"`,
`numerics.execution = "cuda"` (see
[the atom-block layout](#use-the-atom-block-layout-for-large-real-space-models)
and [`numerics`](../reference/run-file.md#numerics)). SIESTA and Wannier
models run on the CPU.

## Use the atom-block layout for large real-space models

For a large real-space model, a `native` build can speed up the sparse
products with the atom-block layout, which keeps the vectors in cache during
each product. The Hamiltonian is then stored not as a plain sparse (CSR)
matrix but in blocks: the orbitals of each atom form one block, neighbouring
atoms are gathered into spatial domains, and atoms and domains are ordered
along a space-filling (Morton) curve. The GPU runs use the same layout. The
moments are traces, so the reordering leaves them unchanged.

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
metric. The `golden` build stops with
`numerics.operator='atom_block' is only available in the native build profile`.
<!-- source: python/lsquant/spatial.py:1-9 (permutation, Morton order), :33-60 (read_positions); python/lsquant/runfile.py:690-734 (validation), :1768-1900 (execution) -->

## Estimate and limit the memory

LSQuant stores Chebyshev vectors of dimension N, the dimension of the
Hamiltonian. One vector of complex doubles takes 16 N bytes: 16 MB for
N = 10⁶. The `expansion` keys of the run file set how many vectors are held
at once. In the table, a random vector is one sample of a random trace, and a
unit vector one basis state of an exact trace.

| Key | Default | Effect on memory | Cost |
|---|---|---|---|
| `expansion.batch` | `moments` | In the `golden` build, a conductivity (σ_xx, σ_xy), which has two Chebyshev expansions, holds two blocks of `batch` vectors, about 2 × batch × 16 N bytes per random vector. | With `batch` < M, the second expansion is repeated for every block of the first, so the time grows about M / batch times. The moments do not change. |
| `expansion.state_block` | `1` | Random or unit vectors propagated together. Memory grows in proportion. `"auto"` picks a value from the available memory in the `native` build, and 1 in `golden`. | Larger blocks use the cache better. |
| `expansion.left_window`, `expansion.right_window` | `32`, `16` | `native` build, conductivities: about (left_window + right_window) × 16 N bytes per random vector; `right_window: "all"` stores all M vectors of the second expansion. | Not used by the `golden` build. |
| `expansion.checkpoint` | `false` | No effect on memory. `csr` models only: saves the accumulated moments after every random or unit vector, so that a run killed by the wall-time limit resumes. | One file write per vector. |
<!-- "conductivity" here means a word with two expansion slots (rank-2 word); batch, left_window and right_window act on that path -->

Example: σ_xx with M = 1024 on N = 10⁶ in the `golden` build holds about
2 × 1024 × 16 MB = 32 GB with the default `batch`. With `batch: 128` it holds
4 GB and takes up to 8 times longer.

SIESTA runs also hold dense k-blocks. Non-orthogonal real-space runs with
`metric.applier = "exact"`, the default, also hold dense factorisations of S,
up to `metric.max_dim` (default 4096). `model.params.bloch_memory_mb`
(default 256) limits the work space of the Bloch sums in clean SIESTA runs.
See [Use a SIESTA calculation](use-siesta.md).

### Resume a killed run

With `expansion.checkpoint: true`, a `csr` run writes
`checkpoint_v2_<hash>.lsqckpt` into the output directory after every random or
unit vector. Run the same command again with the same `--out` directory, and
LSQuant continues after the last finished vector. The file is deleted when the
calculation of the quantity finishes.

> **Checkpoint identity:** the filename fingerprints the resolved window, state
> generator and seed, expansion and time settings, basis, metric, numerics and
> model input contents. Changing the window, seed or matrix bytes starts a fresh
> calculation. Other fingerprints and legacy word-only checkpoints are explicitly
> reported and ignored, and left in place for their original calculations. No
> manual deletion is required. See [the run-file reference](../reference/run-file.md)
> for the fingerprint and hashing cost.

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
set +u                                   # conda's scripts read unset variables
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate lsq                       # or: source .venv/bin/activate
set -u
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

- A batch job runs in a non-interactive shell, which does not read the conda
  setup of your login shell. The `source` line loads it, so that
  `conda activate` works.
- Run `inspect` before `run`. It validates the run file in seconds, so you
  find a typo before the queue wait. Some errors still pass `inspect` and stop
  `run`; see the known limitations in the
  [run-file reference](../reference/run-file.md).
- Relative `model.dir` paths are tried against the LSQuant source tree first,
  then against the directory of the run file. Use absolute paths in batch
  jobs. `--out` is relative to the working directory.
- Set `OPENBLAS_NUM_THREADS=1` and `MKL_NUM_THREADS=1` so that NumPy and SciPy
  keep to one thread each and leave the cores to LSQuant.
- For a `native` build, the job must run on the CPU type that built it. See
  [Install LSQuant](install.md#build-native).
- For disorder averages, use one job per seed, as described in
  [Use the MPI and CUDA builds](#use-the-mpi-and-cuda-builds).

## See also

- [Install LSQuant](install.md)
- [Check a build against the reference outputs](testing.md)
- [Run-file reference: `expansion`](../reference/run-file.md#expansion)
- [Build presets](../reference/build-presets.md)
