\page build Building: presets and build profiles

The new core builds with two numerical **profiles**:
- `golden`: the reference;
- `native`: tuned for the host CPU.

`CMakePresets.json` names both. Preset and manual builds require CMake ≥ 3.21,
a C++17 compiler, and Boost headers for the legacy Boost.Math Bessel functions.

For CUDA builds, an explicit or cached `LSQUANT_CUDA_ARCHITECTURES` wins.
Otherwise, a supplied `CMAKE_CUDA_ARCHITECTURES` initializes the project value;
with neither set, the fallback is `90`. See the
[CUDA architecture reference](reference/build-presets.md#presets).

## Presets

```bash
cmake --preset golden && cmake --build --preset golden       # build/golden
cmake --preset native && cmake --build --preset native       # build/native
export PYTHONPATH=$PWD/build/golden/lsquant_core:$PWD/python
```

| Cache variable | `golden` preset | `native` preset | Meaning |
|---|---|---|---|
| `LSQUANT_BUILD_PROFILE` | `golden` | `native` | numerical profile (below) |
| `LSQUANT_BUILD_PYTHON` | `ON` | `ON` | the pybind11 module `lsquant._core` (needs pybind11 and GSL) |
| `LSQUANT_BUILD_DEV_CHECKS` | `ON` | `ON` | `stage_1_check`, used by CTest and the gate |
| `LSQUANT_BUILD_BENCHMARKS` | `OFF` | `OFF` | the kernel benchmark `stage8_kernels` |
| `LSQUANT_BUILD_INPUT_TOOLS` | `OFF` | `OFF` | fetch and build wannier2sparse at configure time |
| `CMAKE_BUILD_TYPE` | `Release` | `Release` | |

Each preset builds into its own directory (`build/golden`, `build/native`), so both can coexist.
When `build/` holds no build, `scripts/check_goldens.sh` uses `build/golden`.
pybind11 is found from the running Python (`python3 -m pybind11 --cmakedir`). Pass
`-Dpybind11_DIR=...` to choose another one.

## Profiles

| Profile | Flags | Use |
|---|---|---|
| `golden` | portable IEEE arithmetic: no `-march=native`, no fast-math | the reference; the goldens are byte-identical in this profile on the golden toolchain |
| `native` | `-O3 -march=native -mtune=native -fno-math-errno`, `LSQUANT_NATIVE_BUILD=1`; batched kernels, NUMA page placement, cache-blocked partitions | production runs on the build machine |

The native profile changes the order of floating-point operations, so its results differ from
the golden ones in the last digits. The example below measures the difference.
`lsquant._core.build_profile` reports the profile of the loaded module.

## Install

```bash
cmake --install build/golden --prefix /opt/lsquant
```

The install places `liblsquant_core.a` and the headers under `lib/` and `include/lsquant/`,
and `stage1_moments` and `stage6_nonorth` under `bin/`. With `LSQUANT_BUILD_PYTHON=ON` it also
installs the `lsquant` package with its `_core` module under `python/lsquant/`.

## Benchmark

With `-DLSQUANT_BUILD_BENCHMARKS=ON` the build adds `stage8_kernels`:

```bash
build/native/lsquant_core/stage8_kernels spmv|spmm|block|dot|contract SIZE BATCH REPEATS
```

For each kernel it reports the wall time, the modelled FLOP and byte rates, and an output
hash. Run it at several `OMP_NUM_THREADS` values to see the thread scaling.

## Example: golden and native profiles on the same run file

`examples/build_profiles/compare_profiles.sh` builds both presets, runs the same run file with
each (by default `runs/goldens/g5.json`, the Haldane Kubo moments and curves) and compares
the outputs file by file. `examples/build_profiles/README.md` gives the result.
