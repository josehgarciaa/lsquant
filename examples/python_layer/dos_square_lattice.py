#!/usr/bin/env python3
"""Density of states of a square lattice from Python, without input files.

Builds the nearest-neighbour tight-binding Hamiltonian of an L x L square lattice
(hopping t = -1 eV, periodic boundaries) in NumPy, hands it to the C++ engine through
``lsquant._core.SparseMatrix.from_csr``, computes M Chebyshev moments of the density of
states with a stochastic trace (R random-phase vectors), and reconstructs rho(E) with
the Jackson kernel. The result is compared with the exact density of states of the
infinite lattice,

    rho(E) = K(sqrt(1 - (E/4)^2)) / (2 pi^2)        for |E| < 4,

with K the complete elliptic integral of the first kind (modulus k, parameter m = k^2).

Usage (repository root, after building the core with Python bindings):
    PYTHONPATH=build/lsquant_core:python python3 examples/python_layer/dos_square_lattice.py \
        [--L 256] [--moments 512] [--states 8] [--broadening-meV 50] [--out out/python_layer]

It writes dos_square.dat (E in eV, rho per site in 1/eV) and, if matplotlib is
installed, dos_square.png, and prints three checks:
    integral of rho over the grid (1 minus the weight beyond |x| = 0.95),
    the van Hove peak at E = 0 is resolved (rho(0) > rho(+-1)),
    mean |rho - rho_exact| for 1 < |E| < 3, away from the logarithmic peak.
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

import lsquant as lsq
from lsquant import _core


def square_lattice_csr(L: int, t: float = -1.0):
    """CSR arrays of the L x L periodic square lattice, site index i = x + L*y."""
    n = L * L
    rows = np.arange(0, 4 * n + 1, 4, dtype=np.int64)
    cols = np.empty(4 * n, dtype=np.int64)
    for y in range(L):
        for x in range(L):
            i = x + L * y
            nb = sorted({(x + 1) % L + L * y, (x - 1) % L + L * y,
                         x + L * ((y + 1) % L), x + L * ((y - 1) % L)})
            cols[4 * i:4 * i + 4] = nb
    vals = np.full(4 * n, complex(t, 0.0))
    return n, rows, cols, vals


def exact_square_dos(e: np.ndarray) -> np.ndarray:
    """rho(E) of the infinite square lattice with |t| = 1, by the AGM form of K."""
    out = np.zeros_like(e)
    for i, en in enumerate(e):
        if abs(en) >= 4.0 or en == 0.0:
            out[i] = np.nan if en == 0.0 else 0.0
            continue
        k = math.sqrt(1.0 - (en / 4.0) ** 2)
        a, b = 1.0, math.sqrt(1.0 - k * k)          # K(k) = pi / (2 AGM(1, k'))
        for _ in range(60):
            a, b = 0.5 * (a + b), math.sqrt(a * b)
        out[i] = (math.pi / (2.0 * a)) / (2.0 * math.pi ** 2)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--L", type=int, default=256)
    ap.add_argument("--moments", type=int, default=512)
    ap.add_argument("--states", type=int, default=8)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--broadening-meV", type=float, default=50.0, dest="eta")
    ap.add_argument("--out", default="out/python_layer")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    n, rows, cols, vals = square_lattice_csr(a.L)
    ham = _core.SparseMatrix.from_csr(n, rows, cols, vals)
    assert ham.is_hermitian()

    window = _core.SpectralWindow(-4.2, 4.2)            # contains the band [-4, 4]
    states = _core.RandomPhaseStates(n, a.states, a.seed)
    mu = lsq.run_plan(lsq.dos_plan(a.moments, a.states), ham, window, states)
    e, v = lsq.dos_curve(mu, a.eta)
    rho = v / mu.system_size                            # dos_curve returns N * rho(E)

    with open(os.path.join(a.out, "dos_square.dat"), "w") as fh:
        for en, r in zip(e, rho):
            fh.write(f"{en:.8f} {r:.8e}\n")

    integral = float(np.sum(0.5 * (rho[1:] + rho[:-1]) * np.diff(e)))
    exact = exact_square_dos(e)
    sel = (np.abs(e) > 1.0) & (np.abs(e) < 3.0)
    dev = float(np.mean(np.abs(rho[sel] - exact[sel])))
    r0 = float(np.interp(0.0, e, rho))
    r1 = 0.5 * float(np.interp(1.0, e, rho) + np.interp(-1.0, e, rho))
    print(f"N = {n} sites, M = {a.moments}, R = {a.states} random vectors, eta = {a.eta} meV "
          f"(M_eff = {len(e) // 30})")
    print(f"integral of rho over the grid = {integral:.4f}")
    print(f"rho(0) = {r0:.4f} 1/eV, mean rho(+-1) = {r1:.4f} 1/eV (van Hove peak resolved: {r0 > r1})")
    print(f"mean |rho - rho_exact| for 1 < |E| < 3: {dev:.2e} 1/eV "
          f"(rho_exact there: {np.nanmean(exact[sel]):.3f} 1/eV)")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(5, 3.4))
        ax.plot(e, exact, color="0.6", lw=2.5, label="exact (infinite lattice)")
        ax.plot(e, rho, color="C0", lw=1.0, label=f"KPM, M = {a.moments}, R = {a.states}")
        ax.set_xlabel("E (eV)"); ax.set_ylabel(r"$\rho(E)$ (1/eV per site)")
        ax.set_ylim(0, 1.2 * np.nanmax(exact[np.abs(e) > 0.3])); ax.legend(frameon=False)
        fig.tight_layout(); fig.savefig(os.path.join(a.out, "dos_square.png"), dpi=150)
    except ImportError:
        pass
    ok = abs(integral - 1.0) < 0.05 and r0 > r1 and dev < 0.01
    print("checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
