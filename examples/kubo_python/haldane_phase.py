#!/usr/bin/env python3
"""Hall conductivity of the Haldane model on both sides of its topological transition.

Builds the Haldane model on an L x L honeycomb supercell in NumPy (periodic boundaries,
nearest-neighbour hopping t1 = -1 eV, next-nearest-neighbour hopping t2 = 0.15 eV with
phase phi = pi/2, Semenoff mass +m on sublattice A and -m on B), together with the
velocity operators V_a = i[H, X_a] (hbar = 1, lengths in units of the C-C distance).
For each mass it:

  1. computes the rank-2 moments mu_mn = Tr[VX T_m(H~) VY T_n(H~)] and
     Tr[VX T_m VX T_n] with an exact trace over all sites (lsquant.kubo_plan, run_plan);
  2. reconstructs sigma_xy(E) with the Kubo-Bastin formula and sigma_xx(E) with the
     Kubo-Greenwood formula (lsquant.kubo_bastin, lsquant.kubo_greenwood);
  3. extracts the Chern number from the plateau of sigma_xy inside the gap,
     C = 2 pi <sigma_xy>_{|E| < E_avg} / A, with A the supercell area.

The model is a Chern insulator (|C| = 1) for |m| < 3 sqrt(3) t2 sin(phi) = 0.779 eV and a
trivial insulator (C = 0) above. The defaults use m = 0 and m = 1.2 eV.

Usage (repository root, after building the core with Python bindings):
    PYTHONPATH=build/lsquant_core:python python3 examples/kubo_python/haldane_phase.py \
        [--L 16] [--moments 128] [--broadening-meV 10] [--masses 0 1.2] [--out out/kubo_python]
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np

import lsquant as lsq
from lsquant import _core

SQ3 = math.sqrt(3.0)
A1 = np.array([SQ3, 0.0])                  # lattice vectors (C-C distance = 1)
A2 = np.array([SQ3 / 2.0, 1.5])
DELTA_B = np.array([0.0, 1.0])             # position of B in the cell
CELL_AREA = 1.5 * SQ3


def haldane(L, t1=-1.0, t2=0.15, phi=math.pi / 2, mass=0.0):
    """Return (H, VX, VY) as dense complex arrays for the L x L periodic supercell."""
    n = 2 * L * L
    idx = lambda i, j, s: 2 * ((i % L) + L * (j % L)) + s
    pos = np.zeros((n, 2))
    for j in range(L):
        for i in range(L):
            pos[idx(i, j, 0)] = i * A1 + j * A2
            pos[idx(i, j, 1)] = i * A1 + j * A2 + DELTA_B
    h = np.zeros((n, n), dtype=complex)
    vx = np.zeros_like(h)
    vy = np.zeros_like(h)

    def hop(a, b, t, d):
        # H_ab = t for a bond from b to a with displacement d = r_a - r_b;
        # V_ab = i (x_b - x_a) H_ab = -i d H_ab, and the Hermitian conjugate.
        h[a, b] += t
        h[b, a] += np.conj(t)
        vx[a, b] += -1j * d[0] * t
        vx[b, a] += np.conj(-1j * d[0] * t)
        vy[a, b] += -1j * d[1] * t
        vy[b, a] += np.conj(-1j * d[1] * t)

    for j in range(L):
        for i in range(L):
            a = idx(i, j, 0)
            # nearest neighbours of A: B in cells (i, j), (i, j-1), (i+1, j-1)
            for di, dj in ((0, 0), (0, -1), (1, -1)):
                b = idx(i + di, j + dj, 1)
                d = (di * A1 + dj * A2 + DELTA_B)              # r_B - r_A
                hop(b, a, t1, d)
            # next-nearest neighbours along a1, a2 - a1, -a2 (120 degrees apart);
            # phase +phi on A, -phi on B
            for di, dj in ((1, 0), (-1, 1), (0, -1)):
                d = di * A1 + dj * A2
                for s, sign in ((0, 1.0), (1, -1.0)):
                    hop(idx(i + di, j + dj, s), idx(i, j, s), t2 * np.exp(1j * sign * phi), d)
            h[a, a] += mass
            h[idx(i, j, 1), idx(i, j, 1)] -= mass
    return h, vx, vy


def to_core(m):
    n = m.shape[0]
    rows, cols, vals = [0], [], []
    for i in range(n):
        nz = np.nonzero(np.abs(m[i]) > 0)[0]
        cols.extend(int(j) for j in nz)
        vals.extend(complex(m[i, j]) for j in nz)
        rows.append(len(cols))
    return _core.SparseMatrix.from_csr(n, np.asarray(rows, dtype=np.int64),
                                       np.asarray(cols, dtype=np.int64),
                                       np.asarray(vals, dtype=np.complex128))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--L", type=int, default=16)
    ap.add_argument("--moments", type=int, default=128)
    ap.add_argument("--broadening-meV", type=float, default=10.0, dest="eta")
    ap.add_argument("--masses", type=float, nargs="+", default=[0.0, 1.2])
    ap.add_argument("--out", default="out/kubo_python")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    t2, phi = 0.15, math.pi / 2
    m_c = 3.0 * SQ3 * t2 * math.sin(phi)
    ok = True
    print(f"Haldane model, {a.L} x {a.L} cells ({2 * a.L * a.L} sites), M = {a.moments}, "
          f"eta = {a.eta} meV; transition at |m| = {m_c:.3f} eV")
    for mass in a.masses:
        h, vx, vy = haldane(a.L, t2=t2, phi=phi, mass=mass)
        assert np.allclose(h, h.conj().T) and np.allclose(vx, vx.conj().T) and np.allclose(vy, vy.conj().T)
        n = h.shape[0]
        gap_half = abs(abs(mass) - m_c)                          # half-gap at K, K'
        ham = to_core(h)
        window = _core.SpectralWindow(-4.5, 4.5)                 # centred, encloses |E| < 3.5
        ops = [_core.SparseMatrixOperator(to_core(vx)), _core.SparseMatrixOperator(to_core(vy))]
        res = {}
        for name, li, ri, word in (("xy", 0, 1, "VX delta VY delta"), ("xx", 0, 0, "VX delta VX delta")):
            states = _core.LocalStates(n, list(range(n)))
            mu = lsq.run_plan(lsq.kubo_plan(a.moments, left_index=li, right_index=ri, word=word),
                              ham, window, states, operators=ops)
            res[name] = mu
        e, sxy = lsq.kubo_bastin(res["xy"], a.eta)
        e_xx, sxx = lsq.kubo_greenwood(res["xx"], a.eta)      # its own grid of M - 1 points
        area = (n / 2) * CELL_AREA
        e_avg = 0.4 * gap_half
        inside = np.abs(e) < e_avg
        chern = 2.0 * math.pi * float(np.mean(sxy[inside])) / area
        sxx_gap = float(np.mean(sxx[np.abs(e_xx) < e_avg])) / area
        sxx_band = float(np.max(sxx)) / area
        expect = 1 if abs(mass) < m_c else 0
        good = abs(abs(chern) - expect) < 0.05
        ok &= good
        np.savetxt(os.path.join(a.out, f"sigma_xy_m{mass:g}.dat"), np.column_stack([e, sxy / area]),
                   header="E (eV), sigma_xy/A (e = hbar = 1, length unit = C-C distance)")
        np.savetxt(os.path.join(a.out, f"sigma_xx_m{mass:g}.dat"), np.column_stack([e_xx, sxx / area]),
                   header="E (eV), sigma_xx/A (e = hbar = 1, length unit = C-C distance)")
        print(f"m = {mass:4.2f} eV: half-gap {gap_half:.3f} eV, C = {chern:+.3f} (|C| expected {expect}); "
              f"sigma_xx/A in the gap {sxx_gap:.1e}, maximum in the bands {sxx_band:.2f}")
    print("checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
