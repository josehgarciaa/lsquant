#!/usr/bin/env python3
"""Ballistic spreading in a clean chain: the mean-square displacement from Python.

Builds the periodic 1D chain (hopping t = -1 eV, lattice constant a0 = 1) with N sites
in NumPy, and its velocity V = i[H, X] (bond vectors +-1). It then computes, with a
stochastic trace over R random-phase vectors:

  * the time-resolved moments Tr[[X, U(t)]^+ T_m(H~) [X, U(t)]]  (lsquant.msd_plan),
  * the density-of-states moments                                 (lsquant.dos_plan),

and the energy-resolved mean-square displacement

    MSD(E, t) = Tr[[X, U(t)]^+ delta(E - H) [X, U(t)]] / Tr[delta(E - H)].

In a clean chain the motion is ballistic, MSD(E, t) = v(E)^2 t^2 / hbar^2, with the group
velocity v(E) = dE/dk = 2 sin k, so v(E)^2 = 4 - E^2 (eV a0)^2 and hbar = 0.6582 eV fs.
The script prints MSD / t^2 at four energies against (4 - E^2) / hbar^2, for times short
enough that no wave packet has crossed half of the ring (v t / hbar < N / 4).

Usage (repository root, after building the core with Python bindings):
    PYTHONPATH=build/lsquant_core:python python3 examples/time_evolution/ballistic_msd.py \
        [--N 4096] [--moments 256] [--states 4] [--steps 41] [--tmax 100] [--out out/time_evolution]
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

import lsquant as lsq
from lsquant import _core
from lsquant.moments import HBAR


def chain(n, t=-1.0):
    """CSR arrays of H and of V = i[H, X] for the periodic chain."""
    rows = np.arange(0, 2 * n + 1, 2, dtype=np.int64)
    cols = np.empty(2 * n, dtype=np.int64)
    h = np.empty(2 * n, dtype=np.complex128)
    v = np.empty(2 * n, dtype=np.complex128)
    for i in range(n):
        left, right = (i - 1) % n, (i + 1) % n
        pair = sorted([(left, -1.0), (right, +1.0)])          # (column, x_j - x_i)
        for k, (j, dx) in enumerate(pair):
            cols[2 * i + k] = j
            h[2 * i + k] = t
            v[2 * i + k] = 1j * t * dx                         # V_ij = i (x_j - x_i) H_ij
    return rows, cols, h, v


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--N", type=int, default=4096)
    ap.add_argument("--moments", type=int, default=256)
    ap.add_argument("--states", type=int, default=4)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--steps", type=int, default=41)
    ap.add_argument("--tmax", type=float, default=100.0, help="fs")
    ap.add_argument("--broadening-meV", type=float, default=10.0, dest="eta")
    ap.add_argument("--out", default="out/time_evolution")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    rows, cols, h, v = chain(a.N)
    ham = _core.SparseMatrix.from_csr(a.N, rows, cols, h)
    vel = _core.SparseMatrix.from_csr(a.N, rows, cols, v)
    window = _core.SpectralWindow(-2.2, 2.2)
    dt = lsq.time_step_from_tmax(a.steps, a.tmax)
    coeffs, phase = lsq.time_evolution_series(window, dt)

    plan = lsq.msd_plan(a.moments, a.steps, a.states, coeffs, velocity_index=0, phase=phase)
    mu_msd = lsq.run_plan(plan, ham, window, _core.RandomPhaseStates(a.N, a.states, a.seed),
                          operators=[_core.RescaledOperator(vel, lsq.velocity_window(window))], dt=dt)
    mu_dos = lsq.run_plan(lsq.dos_plan(a.moments, a.states, word="1 delta"), ham, window,
                          _core.RandomPhaseStates(a.N, a.states, a.seed))
    e_grid, n_rho = lsq.dos_curve(mu_dos, a.eta)                    # N * rho(E)

    v_max = 2.0
    t_ok = (a.N / 4) * HBAR / v_max                                 # fs
    print(f"chain N = {a.N}, M = {a.moments}, R = {a.states}, {a.steps} steps up to {a.tmax} fs "
          f"(dt = {dt:.3f} fs); ballistic window t < {t_ok:.0f} fs")
    print("   E (eV)   MSD/t^2 (a0^2/fs^2)   (4 - E^2)/hbar^2   ratio")
    ok = True
    rows_out = []
    for ef in (0.0, 0.5, 1.0, 1.5):
        t, msd_rho = lsq.time_correlation_curve(mu_msd, a.eta, ef)      # Tr[... delta(E_F - H) ...] / N
        rho = float(np.interp(ef, e_grid, n_rho)) / a.N
        msd = np.asarray(msd_rho) / rho
        sel = (np.asarray(t) > 0) & (np.asarray(t) < min(a.tmax, t_ok))
        slope = float(np.mean(msd[sel] / np.asarray(t)[sel] ** 2))
        exact = (4.0 - ef * ef) / HBAR ** 2
        ok &= abs(slope / exact - 1.0) < 0.02
        print(f"  {ef:6.2f}   {slope:18.4f}   {exact:16.4f}   {slope / exact:6.4f}")
        rows_out.append(np.column_stack([t, msd]))
        np.savetxt(os.path.join(a.out, f"msd_EF{ef:g}.dat"), np.column_stack([t, msd]),
                   header=f"t (fs), MSD (a0^2) at E_F = {ef} eV")
    print("checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
