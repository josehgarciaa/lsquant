# Check a build against the golden files

This guide shows you how to check that your LSQuant build reproduces the
reference outputs. These reference outputs are the golden files under
`goldens/stage0/`: outputs of the earlier codes, generated once and never
regenerated. The script that compares a build with them,
`scripts/check_goldens.sh`, is called the gate. Run it after you build
LSQuant on a new machine, and after you change numerical code.

The gate is organised in stages, which follow the order in which LSQuant was
developed. Each stage checks one part of the code, from the C++ moment
engine (stage 1) to the run files, the command line and the MCP server
(stage 7).

> **Note.** On every preset build, stage 4 ends with `GATE FAIL`, even when
> every number is right. Two checks of the C++ source report `FOUND` whatever
> the numbers are. Judge your build by the per-file lines instead: if the only
> failures are the two `GREP … FOUND` lines, the numerical results are
> correct. For the cause, see the known limitations under
> [Read the result](#read-the-result).

[Testing: golden masters and the regression gate](../testing.md) explains what
the golden files contain and the rules of the tolerance mode. It also shows how
to read a `DIFFERS` line and how to compare two files by hand with
`scripts/golden_diff.py`. This page covers the task: which stage to run on a
preset build, what it needs, and how long it takes.

## Before you begin

You need:

- A build from a [preset](install.md), normally `golden`. The gate uses the
  build directory `build/golden/` by default.
- `python3` on your `PATH` must be the interpreter you built for, with NumPy,
  SciPy and SymPy. The gate runs `python3` and sets `PYTHONPATH` itself,
  ignoring the value in your environment.
- For stages 5, 6 and 7: the `wannier2sparse` binary. See
  [Get wannier2sparse](#get-wannier2sparse).
- For stage 7: the `mcp` Python package, because the stage also tests the MCP
  server.
- A compute node. The gate fixes its own thread counts, up to
  `OMP_NUM_THREADS=64` in stage 1, whatever you set. Do not run it on a login
  node.

## Choose a stage

The script takes a stage and, optionally, the build directory:

```text
scripts/check_goldens.sh <stage> [build_dir]
```

Stages `1` to `6` are cumulative: each one first runs every stage below it,
starting from stage 1. Stage `6`, for example, runs stages 1 to 5 first. The
stages `5only`, `6only` and `7runs` run one part alone.

The times below are for the `golden` build on an otherwise idle 96-core node.

| Stage | What it checks | Needs | Time |
|---|---|---|---|
| `4` | The C++ moments of a chain at 1, 8 and 64 threads; the Python binding and the DOS reconstruction; Kubo–Bastin and Kubo–Greenwood moments and curves of graphene and the Haldane model; time evolution of a chain. Also greps of the source in `lsquant/core` and `lsquant/io`. | build | 14 min |
| `7runs` | Every run file under `runs/goldens/` through `python -m lsquant run`, compared with its golden file. With `S7_SET=quick`, a subset of 15 run files at 16 threads. | build, `wannier2sparse` | 8 min with `S7_SET=quick` |
| `5only` | Real, Bloch and eigen bases of graphene and the Haldane model. | build, `wannier2sparse` | not measured |
| `6only` | Non-orthogonal metric: Löwdin and generalised routes, `exact`, `ischo` and `chebyshev` appliers, and the Hall plateau of the Haldane model. | build | not measured |
| `7` | All run files at 16 and 1 threads, the MCP server, the command-line refusals and the source greps. | build, `wannier2sparse`, `mcp` | about 10 h for `6` and `7` together |

The header of `scripts/check_goldens.sh` describes every stage and every
environment variable in detail.

For a new machine, run `4` and then `7runs` with `S7_SET=quick`. Run the full
`6only` and `7` before you accept a change to numerical code.

## Run the gate

1. Go to the top of the source tree.

2. For a `golden` build on the machine where the golden files were written,
   run the stage:

   ```bash
   scripts/check_goldens.sh 4
   ```

3. For a `native` build, or for a `golden` build on any other machine, turn on
   the tolerance mode and pass the build directory:

   ```bash
   GOLDEN_TOL=1 scripts/check_goldens.sh 4 build/native
   ```

   Another compiler, C library or Eigen version reorders floating-point
   operations, so some files differ in their last printed digit. The tolerance
   mode accepts those differences and reports the file as `TOLPASS`.

4. Run the run-file check:

   ```bash
   GOLDEN_TOL=1 S7_SET=quick scripts/check_goldens.sh 7runs build/golden
   ```

The gate writes only to a temporary directory under `$TMPDIR` (or `/tmp`) and
deletes it when it ends. It never writes into `goldens/`.

## Read the result

The gate prints one line per compared file, then a summary. This is the end
of `GOLDEN_TOL=1 scripts/check_goldens.sh 4` on a `golden` build on a node
with gcc 9.4. The known limitations below explain why it fails:

```text
== 71 files compared; failures: 8
== tolerance mode (GOLDEN_TOL=1): 49 byte-identical, 16 tolerance-PASS (compare() and index.tsv rows)
GATE FAIL (stage 4)
```

Each line starts with one of these words:

| Line starts with | Meaning |
|---|---|
| `IDENTICAL` | The file is byte-identical to its golden file. |
| `TOLPASS` | Only with `GOLDEN_TOL=1`: the file differs, but within the tolerance. It counts as a pass. |
| `DIFFERS` | The file differs beyond the tolerance, or was not produced. The next lines locate the difference. |
| `GREP` | A check on the C++ source. `FOUND` counts as a failure. |

The exit code of the script tells you the verdict:

| Exit code | Meaning |
|---|---|
| 0 | `GATE PASS`: every comparison and every grep passed. |
| 1 | `GATE FAIL`: at least one `DIFFERS` or one `GREP … FOUND`. |
| 2 | The gate could not start: unknown stage, missing build directory, or a missing program such as `stage1_moments` or `wannier2sparse`. |

To read a `DIFFERS` line, see
[Testing: golden masters and the regression gate](../testing.md).

> **Known limitation: stages 1 to 4 and 7 end in `GATE FAIL` on source
> greps.** The `GREP` lines for `getenv` report `FOUND`, because the C++
> engine reads three environment variables
> (`lsquant/core/chebyshev/moment_engine.cpp:1258`, `1273`, `1473`;
> `scripts/check_goldens.sh:389`, `:420–430`). The stages `7runs`, `5only`
> and `6only` skip the greps.

> **Known limitation: stage 0 does not run on a preset build.** It stops with
> `legacy driver … missing`, because it needs programs of the earlier code
> that the presets do not build (`scripts/check_goldens.sh:327–331`). Run
> stage 7 instead, which builds them itself in a temporary worktree,
> `/tmp/legacy-opus7` by default; set `S7_LEGACY_WT` to move it. Stage 7
> deletes that directory first (`scripts/check_goldens.sh:1506`).

On a machine other than the one where the golden files were written, a few
`DIFFERS` can remain even with `GOLDEN_TOL=1`. The tolerance, one unit in the
last printed digit, fails by construction on curves that are numerically
zero. On a node with gcc 9.4, a `golden` build gave 6 `DIFFERS` in stage 4
and 16 in `7runs` with `S7_SET=quick`. All were at a relative level below
10⁻¹², and identical at 1 and 16 threads. A byte-for-byte verdict needs the
reference machine and toolchain, recorded in `goldens/stage0/PROVENANCE.md`.

## Get wannier2sparse

Stages 5 and 7 build real-space supercells with `wannier2sparse`, an external
program, at a fixed commit. The script `scripts/fetch_wannier2sparse.sh`
clones and builds it into `.cache/wannier2sparse/`. It needs network access to
GitHub the first time.

Run the script and point the gate to the binary:

```bash
scripts/fetch_wannier2sparse.sh
export W2S_BIN="$PWD/.cache/wannier2sparse/build/wannier2sparse"
```

If `W2S_BIN` is unset, the gate looks for
`$HOME/wannier2sparse/build/wannier2sparse`. If Eigen is installed outside the
standard locations, set `EIGEN_DIR` before you run the fetch script.

## Environment variables

| Variable | Stage | Effect |
|---|---|---|
| `GOLDEN_TOL=1` | all | Tolerance mode. |
| `W2S_BIN` | 5, 7 | Path of the `wannier2sparse` binary. |
| `S7_SET` | 7, `7runs` | `all` (default), `quick`, or `list` with the names in `S7_LIST`. |
| `S7_THREADS` | 7, `7runs` | Thread counts of the run files (default `16 1`). |
| `S7_KEEP=<dir>` | 7, `7runs` | Keep the outputs of the run files in `<dir>`. |
| `S7_LEGACY_WT` | 7 | Directory of the temporary worktree for the earlier code (default `/tmp/legacy-opus7`). The gate deletes this directory before it starts, so point it to a path that does not exist. |
| `S6_THREADS` | 6 | Thread count of the stage-6 runs (default 16). |
| `S6_KEEP=<dir>`, `S5_KEEP=<dir>` | 6, 5 | Keep the outputs of these stages. |

The gate fixes the seed (`KPM_SEED=12345`), the tolerances and the golden
files. You cannot change them from the command line.

## Run the gate as a batch job

The gate takes minutes to hours, so submit it as a job. A minimal Slurm script:

```bash
#!/bin/bash
#SBATCH --job-name=lsq-gate
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00

conda activate lsq            # or: source .venv/bin/activate
cd "$HOME/lsquant"
GOLDEN_TOL=1 scripts/check_goldens.sh 4 build/golden > gate4.log 2>&1
echo "gate exit code: $?" >> gate4.log
```

See [Run on a cluster](run-on-a-cluster.md) for thread settings and batch
scripts.
