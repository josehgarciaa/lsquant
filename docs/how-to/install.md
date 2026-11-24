# Install LSQuant

The steps below install the dependencies, build LSQuant with a CMake preset
and make the Python package importable. The main route is for a cluster
without root access; the
[route with root access](#install-the-dependencies-with-root-access) follows
it. To build LSQuant together with a first calculation, follow
[Getting started](../tutorials/getting-started.md) instead.

LSQuant has two parts:

- a compiled C++ engine with a Python binding, the module `lsquant._core`,
  which CMake builds;
- the Python package `lsquant` in `python/lsquant/`, which is used directly
  from the source tree.

Both are made importable with `PYTHONPATH`; LSQuant has no `pip install`
step.

The [build-preset reference](../reference/build-presets.md) lists every preset,
cache option and dependency version in tables.

## Before you begin

You need:

- A C++17 compiler with OpenMP, to build the engine. The reference build
  used the system compiler, GCC 9.4.
- CMake 3.21 or later for both preset and manual builds.
- Eigen 3.4 or later (the legacy drivers do not compile with 3.3), Boost headers (for legacy Boost.Math Bessel functions),
  FFTW3 (double precision) and `pkg-config`, for the engine.
- GSL, Python 3.10 or later with its development headers, and pybind11, for
  the Python binding.
- NumPy, SciPy and SymPy, for the Python package.
- Optional: `sisl` for SIESTA models, `mcp` for the MCP server, an MPI
  library or the CUDA Toolkit for the MPI and GPU presets.

The [build-preset reference](../reference/build-presets.md#requirements)
gives the details of how CMake finds each dependency.

Choose the preset you need:

| You want to | Preset |
|---|---|
| Build for the first time, check the build, or get results that are identical on every machine | `golden` |
| Run production calculations on one node | `native` |
| Build the distributed (MPI) operator, which run files do not use yet (see [below](#build-native-mpi)) | `native-mpi` |
| Run on an NVIDIA GPU | `native-cuda` |
| Build both MPI and CUDA | `full` |

Start with `golden` even if you plan to use another preset. It is the build
that the [check against the reference outputs](testing.md) compares byte for
byte.

## Install the dependencies without root access

On a cluster you usually lack root access, so `apt-get` is not available.
Install the libraries, CMake and Python into a conda environment instead. The
commands are the same with `mamba` or `micromamba`.

1. Create and activate the environment:

   ```bash
   conda create -n lsq -c conda-forge python=3.10 cmake pkg-config eigen boost-cpp fftw gsl pybind11 \
       numpy scipy sympy matplotlib
   conda activate lsq
   ```

   Boost headers are explicitly included for the legacy drivers.

2. Install the optional packages you need:

   ```bash
   pip install sisl mcp
   ```

3. Check that a C++17 compiler with OpenMP is on your `PATH`, for example
   with `g++ --version`. If your cluster uses environment modules, load a
   compiler module first, for example `module load gcc`. The module names
   depend on your site.

When you configure in the next section, add these options so that CMake uses
the libraries and the interpreter of the environment:

```bash
-DPython_EXECUTABLE="$(which python)" -DGSL_ROOT_DIR="$CONDA_PREFIX"
```

CMake finds Eigen in `$CONDA_PREFIX/include/eigen3` on its own, and
`pkg-config` from the environment finds FFTW3.

If you prefer the libraries of your site's environment modules, load them and
point CMake at them with `-DEIGEN3_INCLUDE_DIR=…` and `-DGSL_ROOT_DIR=…`.
CMake always finds FFTW3 through `pkg-config`, so add its `lib/pkgconfig`
directory to `PKG_CONFIG_PATH`.

Continue with [Build the golden preset](#build-the-golden-preset).

## Install the dependencies with root access

On a workstation with Ubuntu or Debian, use the system packages:

1. Install the compiler, CMake and the libraries:

   ```bash
   sudo apt-get install g++ cmake pkg-config libeigen3-dev libfftw3-dev libgsl-dev python3-dev
   ```

2. Create a virtual environment and install the Python packages:

   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   pip install numpy scipy sympy pybind11 matplotlib
   ```

3. Install the optional packages you need: `sisl` for SIESTA models, `mcp`
   for the MCP server.

   ```bash
   pip install sisl mcp
   ```

## Build the golden preset

Run these commands from the top of the source tree.

1. Configure:

   ```bash
   cmake --preset golden
   ```

   Add the options of the conda section if you use a conda environment.

2. Build:

   ```bash
   cmake --build --preset golden
   ```

   The build preset runs 8 parallel jobs. On a shared login node, use fewer:
   `cmake --build --preset golden -j 2`.

3. Check that the module exists:

   ```bash
   ls build/golden/lsquant_core/
   ```

   If the build has succeeded, the directory contains the Python module
   `_core.cpython-<version>-<platform>.so`. The other executables in it are
   test programs, used by the check against the reference outputs.
   <!-- source: build/golden/lsquant_core also holds stage1_moments, stage6_nonorth and stage_1_check (docs/build.md:21, 49; docs/core-engine.md:20) -->

To configure without a preset (also requires CMake 3.21 or later):

```bash
cmake -S . -B build/golden -DCMAKE_BUILD_TYPE=Release -DLSQUANT_BUILD_PROFILE=golden
cmake --build build/golden -j 8
```

## Set PYTHONPATH

Python needs two directories: the build directory that holds `_core`, and
`python/`, which holds the package.

```bash
export PYTHONPATH="$PWD/build/golden/lsquant_core:$PWD/python"
```

Run this from the top of the source tree, or write absolute paths. Add the
line to your shell start-up file (such as `.bashrc` in your home directory)
or to your batch scripts. For another preset, replace `golden` with the
preset name.

Keep one build on `PYTHONPATH` at a time, because Python loads the first
`_core` it finds.

> **Known limitation: `pip install -e python/` fails.** The package metadata
> points to the top-level `README.md`, a file outside the package directory,
> and recent versions of setuptools refuse it. With setuptools 84.0 the install
> stops with `Cannot access '…/python/../README.md' (or anything outside
> '…/python')`. Use `PYTHONPATH` as shown above, and run the command-line
> interface as `python -m lsquant`.
<!-- source: python/pyproject.toml:9 (readme = "../README.md") -->

## Check the installation

1. Import the package and print the build profile:

   ```bash
   python -c "import lsquant; from lsquant import _core; print(_core.build_profile)"
   ```

   The command prints `golden` (or `native` for the other presets).

2. List the quantities the engine accepts:

   ```bash
   python -m lsquant words list
   ```

   The command prints, as JSON, the seven forms of trace that LSQuant can
   evaluate. Each form is a template for the *observable*, the quantity to
   compute, which a run file gives in the key `word` (see the
   [observable reference](../reference/observable-words.md)).
   <!-- source: python/lsquant/words/catalogue.py:45-123 (7 entries); compiler.py:764-766; cli.py:171-174 -->

On import, `lsquant` sets `OMP_NUM_THREADS` (the number of physical cores the
process may use) and, for the `lsquant` command with its own cores only
(a CPU set from the batch system, or `LSQUANT_BIND_THREADS=1`), `OMP_PROC_BIND=spread`
and `OMP_PLACES=cores`; each only if it is unset. Set `OMP_NUM_THREADS` yourself on a shared node; see
[Run on a cluster](run-on-a-cluster.md).
<!-- source: python/lsquant/runtime.py configure_threads, binding_requested -->

If the import fails:

| Message | Cause | Fix |
|---|---|---|
| `ModuleNotFoundError: No module named 'lsquant'` | `python/` is missing from `PYTHONPATH`. | Set `PYTHONPATH` as above. |
| `ModuleNotFoundError: No module named '_core'` | The build directory is missing from `PYTHONPATH`, or the build failed. | Check that `build/<preset>/lsquant_core/_core*.so` exists. |
| The same message, although `_core*.so` exists | The module was built for another Python version: the version is part of the file name, for example `cpython-310`. | Configure again with `-DPython_EXECUTABLE=$(which python)` and rebuild. |
| `ImportError: libgsl.so…` or `libfftw3.so…: cannot open shared object file` | The libraries of a conda environment are missing from the library search path at run time. | Activate the environment, or add `$CONDA_PREFIX/lib` to `LD_LIBRARY_PATH`. |

## Build the other presets

Each preset builds into its own directory, `build/<preset>/`, so you can keep
several builds side by side.

### Build `native`

```bash
cmake --preset native
cmake --build --preset native
export PYTHONPATH="$PWD/build/native/lsquant_core:$PWD/python"
```

The `native` preset compiles with `-march=native`, so the module runs only on
CPUs with the instruction set of the machine that built it. On a cluster,
build on a compute node of the partition where you will run, rather than on
the login node. For a heterogeneous cluster, keep one `native` build per CPU
type.

A `native` build gives numbers that differ from `golden` in the last digits.
Check it with `GOLDEN_TOL=1`; see [Check a build against the reference outputs](testing.md).

### Build `native-mpi`

You also need an MPI library with C++ support, for example Open MPI or MPICH.
Load it, or install `openmpi` into the conda environment, then run:

```bash
cmake --preset native-mpi
cmake --build --preset native-mpi
```

After the import, `_core.mpi_available` is `True`.

> **Known limitation: run files cannot currently use MPI.** An MPI build accepts
> `numerics.execution = "mpi"` and `"mpi+cuda"` in a run file, but the run
> stops with `numerics.execution='mpi' is not wired yet`. A `native-mpi` build
> runs the same multithreaded calculations as `native`.
<!-- source: python/lsquant/runfile.py:1896–1899 -->

### Build `native-cuda`

You also need the CUDA Toolkit (`nvcc`) and an NVIDIA GPU.

```bash
cmake --preset native-cuda
cmake --build --preset native-cuda
```

After the import, `_core.cuda_available` is `True`.

The GPU code is compiled for compute capability 9.0 (H100 class) by default.
Only 90 is tested. For another GPU, set `LSQUANT_CUDA_ARCHITECTURES`; see
[CUDA architecture](../reference/build-presets.md#presets).

The GPU is used by two run-file settings, both for real-space `csr` models.
The first is the atom-block layout, `numerics.operator = "atom_block"` with
`numerics.execution = "cuda"`. The second is the Ginkgo back end,
`numerics.backend.executor = "cuda"`, which also needs a build with
`-DLSQUANT_WITH_GINKGO=ON`. See [`numerics`](../reference/run-file.md#numerics)
in the run-file reference.

### Build `full`

```bash
cmake --preset full
cmake --build --preset full
```

`full` needs both MPI and the CUDA Toolkit, and has the limitations of both
`native-mpi` and `native-cuda`.

## If the configuration fails

CMake stops with one of these messages when a dependency is missing:

| Message starts with | Fix |
|---|---|
| `Eigen3 (>= 3.3) not found` | Install Eigen, or pass `-DEIGEN3_INCLUDE_DIR=/path/to/eigen3`. |
| `OpenMP for C++ not found` | Use a compiler with OpenMP support (GCC, or Clang with `libomp`). |
| `pkg-config not found` | Install `pkg-config`. |
| `FFTW3 (double precision, 'fftw3') not found via pkg-config` | Install FFTW3, or add its `lib/pkgconfig` directory to `PKG_CONFIG_PATH`. |
| `GSL (GNU Scientific Library) not found` | Install GSL, or pass `-DGSL_ROOT_DIR=/path/to/gsl`. |
| `pybind11 not found` | `pip install pybind11` in the Python you build for, or pass `-Dpybind11_DIR=$(python -m pybind11 --cmakedir)`. |
| `LSQUANT_WITH_CUDA/LSQUANT_WITH_MPI require LSQUANT_BUILD_PROFILE=native` | You enabled MPI or CUDA on top of the golden profile. Use `native-mpi`, `native-cuda` or `full`. |

## Next steps

- Check your build: [Check a build against the reference outputs](testing.md).
- Run your first calculation: [Getting started](../tutorials/getting-started.md).
- Run on a cluster: [Run on a cluster](run-on-a-cluster.md).
