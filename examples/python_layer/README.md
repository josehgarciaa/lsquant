# Python layer: density of states of a square lattice

`dos_square_lattice.py` drives the C++ core from Python with the matrix built in memory:

1. builds the nearest-neighbour Hamiltonian of an L × L square lattice in NumPy
   (t = −1 eV, periodic boundaries, CSR arrays);
2. hands it to the engine with `lsquant._core.SparseMatrix.from_csr`;
3. computes M Chebyshev moments of the density of states with `lsquant.dos_plan` and
   `lsquant.run_plan`, using a stochastic trace over R random-phase vectors;
4. reconstructs ρ(E) with `lsquant.dos_curve` (Jackson kernel);
5. compares it with the exact density of states of the infinite lattice,
   ρ(E) = K(√(1 − E²/16))/(2π²) for |E| < 4, with K the complete elliptic integral of
   the first kind.

## Run

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -Dpybind11_DIR=$(python3 -m pybind11 --cmakedir)
cmake --build build -j 16
PYTHONPATH=build/lsquant_core:python python3 examples/python_layer/dos_square_lattice.py
```

Options: `--L` (default 256), `--moments` (512), `--states` (8), `--seed` (12345),
`--broadening-meV` (50), `--out` (`out/python_layer`).

Output with the defaults on a 64-core workstation (Ubuntu 22.04, g++ 11.4), 8 threads,
2.2 s wall time:

```
N = 65536 sites, M = 512, R = 8 random vectors, eta = 50.0 meV (M_eff = 264)
integral of rho over the grid = 0.9980
rho(0) = 0.3258 1/eV, mean rho(+-1) = 0.1411 1/eV (van Hove peak resolved: True)
mean |rho - rho_exact| for 1 < |E| < 3: 1.03e-03 1/eV (rho_exact there: 0.112 1/eV)
checks: PASS
```

## Reading the result

- **Sum rule.** The grid covers x ∈ [−0.95, 0.95], that is |E| < 3.99 eV for the window
  [−4.2, 4.2]. The missing 0.2 % is the weight of the band edges beyond that.
- **Van Hove singularity.** The exact ρ(E) diverges logarithmically at E = 0. The
  reconstruction shows it as a peak of finite height, 0.326 1/eV, set by the 50 meV
  broadening (M_eff = 264 of the 512 moments).
- **Accuracy.** Away from the peak, the mean deviation is 1.0e-3 1/eV, 0.9 % of ρ there.
  Two errors contribute: the 50 meV Jackson broadening, and the stochastic trace, whose
  relative error is of order 1/√(R·N) ≈ 0.14 %.

The script exits with status 1 if one of the three checks fails. It writes `dos_square.dat`
(E in eV, ρ per site in 1/eV) and, when matplotlib is installed, `dos_square.png`.

The window must contain the spectrum, and its centre must be 0 at this stage: see the
*Energy axis* paragraph of the manual page *Python layer* (`docs/python-layer.md`).
