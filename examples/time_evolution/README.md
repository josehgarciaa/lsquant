# Time evolution: ballistic spreading in a clean chain

`ballistic_msd.py` computes the energy-resolved mean-square displacement (MSD) of a clean
periodic chain from Python and compares it with the ballistic law.

The script:
1. builds the chain Hamiltonian (hopping t = −1 eV, lattice constant a₀ = 1, N = 4096
   sites, periodic) and its velocity V = i[H, X] in NumPy;
2. computes the MSD moments Tr[[X, U(t)]† T_m(H̃) [X, U(t)]] with `lsquant.msd_plan`, over
   41 time steps up to 100 fs, M = 256 moments and R = 4 random-phase vectors;
3. computes the density-of-states moments with the same vectors;
4. reconstructs both at 10 meV and divides them:
   MSD(E, t) = Tr[[X, U]† δ(E − H) [X, U]]/Tr[δ(E − H)].

In a clean chain every state moves at its group velocity v(E) = dE/dk = 2 sin k, so
MSD(E, t) = v(E)² t²/ħ² with v(E)² = 4 − E² (eV a₀)² and ħ = 0.6582 eV fs.

## Run

```bash
PYTHONPATH=build/lsquant_core:python python3 examples/time_evolution/ballistic_msd.py
```

Options: `--N` (4096), `--moments` (256), `--states` (4), `--seed` (12345), `--steps` (41),
`--tmax` (100 fs), `--broadening-meV` (10), `--out` (`out/time_evolution`).

Output on a 64-core workstation (Ubuntu 22.04, g++ 11.4), 8 threads, 1.8 s wall time:

```
chain N = 4096, M = 256, R = 4, 41 steps up to 100.0 fs (dt = 2.500 fs); ballistic window t < 337 fs
   E (eV)   MSD/t^2 (a0^2/fs^2)   (4 - E^2)/hbar^2   ratio
    0.00               9.2312             9.2327   0.9998
    0.50               8.6523             8.6557   0.9996
    1.00               6.9154             6.9245   0.9987
    1.50               4.0499             4.0393   1.0026
checks: PASS
```

The script writes `msd_EF<E>.dat` (t in fs, MSD in a₀²) and exits with status 1 if a ratio
differs from 1 by more than 2 %.

## Reading the result

- **Ballistic law.** MSD/t² is constant in time and equals v(E)²/ħ² within 0.3 % at all
  four energies. The deviation is the stochastic trace (R·N = 16 384) plus the 10 meV
  broadening, which averages v² over neighbouring energies.
- **Time window.** The fastest states (v = 2 eV a₀ at E = 0) cross a quarter of the ring,
  N/4 sites, after 337 fs. The script only uses t ≤ 100 fs, where the periodic boundary
  plays no role.
- **Velocity input.** `msd_plan` takes the velocity of the rescaled Hamiltonian, passed as
  `RescaledOperator(vel, window)`, on a window centred on 0: [−2.2, 2.2] here. The manual
  page *Time evolution* (`docs/time-evolution.md`) explains why the window must be centred.
- **Time step.** `run_plan(..., dt=dt)` records the time step in the tensor, and
  `time_correlation_curve` uses it to return the time axis.
