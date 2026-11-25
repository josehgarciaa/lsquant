# Reduced goldens of the SIESTA run files (G11)

One JSON file per run file `runs/goldens/g11_<material>_<velocity>_<route>_<dos|shc>.json`, written
by `scripts/siesta_golden.py make <run output directory> <file>` from a run of
`python -m lsquant run runs/goldens/<name>.json`: every two-column curve sampled on at most 256
rows with its number of rows, largest absolute value and integral, and the numbers of the exact
oracle (`spin_chern.json`, `exact_shc.json`).

Produced at commit `f41c384` with the `golden` preset: GCC 11.4.0 (Ubuntu 22.04), Eigen 3.4.0,
FFTW 3.3.8, GSL 2.7.1, Python 3.10.12, NumPy 2.2.6, SciPy 1.15.2, sisl 0.16.4, 20 threads.

These are compared with a relative tolerance, not byte for byte (stage 7 of
`scripts/check_goldens.sh`; `S7_SIESTA_RTOL`, `S7_SIESTA_ORACLE_RTOL`). Measured between two machines
with different BLAS libraries: reconstructed density of states 8e-14, reconstructed spin Hall curve
7e-10, exact oracle 7e-7. The tolerances are 1e-7 for the curves and 1e-5 for the oracle.

To regenerate one after an intended change: run the run file, then
`python scripts/siesta_golden.py make <output directory> goldens/stage7_siesta/<name>.json`, and
review the difference in the pull request.
