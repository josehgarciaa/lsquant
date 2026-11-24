# Build presets and options

LSQuant builds with CMake presets defined in `CMakePresets.json`. Each preset
has a configure preset and a build preset of the same name, and writes to
`build/<preset>/`.

```bash
cmake --preset <name>
cmake --build --preset <name>          # or: cmake --build build/<name> -j
export PYTHONPATH="$PWD/build/<name>/lsquant_core:$PWD/python"
```

## Presets

| Preset | Profile | Extra components | Compiler flags | Use it for |
|---|---|---|---|---|
| `golden` | `golden` | — | Portable IEEE arithmetic, no machine-specific or fast-math flags. | First build, checks against the golden files, results that must be reproducible across machines. |
| `native` | `native` | — | `-O3 -march=native -mtune=native -fno-math-errno`, `LSQUANT_NATIVE_BUILD=1` (SIMD, NUMA-aware, cache-blocked kernels). | Production runs on one node. |
| `native-mpi` | `native` | MPI distributed operator (`LSQUANT_WITH_MPI=ON`) | as `native` | Systems that do not fit in the memory of one node. |
| `native-cuda` | `native` | CUDA atom-block kernel (`LSQUANT_WITH_CUDA=ON`) | as `native` | Runs on an NVIDIA GPU. |
| `full` | `native` | MPI and CUDA | as `native` | Several GPU nodes. |

All presets set `CMAKE_BUILD_TYPE=Release`, build the Python module
(`LSQUANT_BUILD_PYTHON=ON`) and the developer self-check executable
(`LSQUANT_BUILD_DEV_CHECKS=ON`), and skip the benchmarks.

**CUDA architecture.** `native-cuda` and `full` compile the CUDA kernel for the
compute capabilities in the cache variable `LSQUANT_CUDA_ARCHITECTURES`. The default
is `90` (Hopper, for example H100 or H200). Only 90 is tested; other values are
passed to the CUDA compiler unchanged. For another GPU, give its compute capability
(for example `86`), a list such as `"80;86;90"`, or `native` for the GPU of the build
machine.
`CMAKE_CUDA_ARCHITECTURES` does not change the kernel; use the variable above.

A `native` build produces numbers that differ from the `golden` build in the
last digits, so check a `native` build against the golden files with
`GOLDEN_TOL=1` (see [Check a build](../how-to/testing.md)).

## Requirements

| Dependency | Version | Needed by | How CMake finds it |
|---|---|---|---|
| C++ compiler | C++17 | all | `CMAKE_CXX_COMPILER` |
| CMake | ≥ 3.21 | all | — |
| Eigen | ≥ 3.4 (3.4.0 wrote the goldens) | all | CMake package `Eigen3`, otherwise the header directory in `$CONDA_PREFIX/include/eigen3`, `/usr/include/eigen3` or `/usr/local/include/eigen3`. Override with `-DEIGEN3_INCLUDE_DIR=…`. |
| Boost headers | — | legacy drivers (`kpm_lib`) | `find_package(Boost)`; override with `-DBoost_ROOT=…`. No compiled Boost library is needed. |
| OpenMP | — | all | compiler support |
| FFTW3, double precision | — | all | `pkg-config fftw3` (so `pkg-config` itself is required) |
| GSL | — | Python module | `find_package(GSL)`; override with `-DGSL_ROOT_DIR=…` |
| Python, with development headers | ≥ 3.10 | Python module | `-DPython_EXECUTABLE=…` |
| pybind11 | — | Python module | `python -m pybind11 --cmakedir`; override with `-Dpybind11_DIR=…` |
| MPI (C++) | — | `native-mpi`, `full` | `find_package(MPI)` |
| CUDA Toolkit | — | `native-cuda`, `full` | `find_package(CUDAToolkit)` |

On Ubuntu or Debian:

```bash
sudo apt-get install g++ cmake pkg-config libeigen3-dev libboost-dev libfftw3-dev libgsl-dev python3-dev
```

## Cache options

| Option | Default | Meaning |
|---|---|---|
| `LSQUANT_BUILD_PROFILE` | set by the preset | `golden` or `native`. |
| `LSQUANT_BUILD_PYTHON` | `ON` | Build the pybind11 module `lsquant._core`. |
| `LSQUANT_BUILD_DEV_CHECKS` | `ON` | Build the developer self-check executable `stage_1_check`. |
| `LSQUANT_BUILD_BENCHMARKS` | `OFF` | Build the kernel benchmark `stage8_kernels`. |
| `LSQUANT_WITH_MPI` | `OFF` | Build the MPI distributed operator (native profile only). |
| `LSQUANT_WITH_CUDA` | `OFF` | Build the CUDA atom-block kernel (native profile only). |
| `LSQUANT_WITH_GINKGO` | `OFF` | Build the optional Ginkgo linear-operator backend (native profile only). |

## Python packages

The Python front end needs `numpy`, `scipy` and `sympy`. The examples and some
model loaders also use:

| Package | Needed for |
|---|---|
| `matplotlib` | the `plot.py` scripts of the examples |
| `sisl` | `model.kind = "siesta"` |
| `mcp` | the MCP server |

```bash
pip install numpy scipy sympy pybind11 matplotlib sisl mcp
```
