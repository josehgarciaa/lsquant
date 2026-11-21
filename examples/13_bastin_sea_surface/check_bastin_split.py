#!/usr/bin/env python3
"""Check the Fermi-sea / Fermi-surface split of the Kubo-Bastin Hall conductivity.

One rank-2 moment tensor mu_mn of Tr[VX delta(E-H) VY dG-(E)/dE] is
reconstructed four ways by run_bastin_split.json:

    bastin      full Kubo-Bastin sigma_xy(E)            KuboBastin_...JACKSON.dat
    sea         Fermi-sea part (Green-function term)    KuboBastinSea_...JACKSON.dat
    surf        Fermi-surface part (delta-delta term)   KuboBastinSurf_...JACKSON.dat
    bastin_fft  full Kubo-Bastin on the M Chebyshev     KuboBastin_...JACKSON_FFTgrid.dat
                nodes x_j = cos(pi (j + 1/4) / (M/2))

Checks:
  1. sea + surf reproduces bastin on the common uniform grid.
  2. bastin_fft gives the same Hall plateau as bastin. (Away from the gap the
     two differ: the FFT-node curve integrates dsigma/dE on only M - 1 nodes,
     so it is coarse where sigma_xy changes fast; the deviation is printed.)
  3. In the Haldane gap the surface term vanishes and the sea term carries the
     quantised Hall conductance: C = 2 pi <sigma_xy> / A = +1
     (sigma_xy extensive, e = hbar = 1, A the total area).

Usage:
    python check_bastin_split.py <run_out_dir>
"""
import glob
import os
import sys

import numpy as np

A_CELL = 3.0 * np.sqrt(3.0) / 2.0     # honeycomb cell area, nearest-neighbour distance 1
DIM = 512


def load(pattern):
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise FileNotFoundError(pattern)
    return np.loadtxt(paths[0])


def chern(e, s, half=0.5):
    mask = np.abs(e) < half
    return float(np.mean(s[mask]) / (DIM / 2 * A_CELL) * 2 * np.pi)


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "out/bastin_split"
    bas = load(os.path.join(d, "KuboBastin_NonEqOpVX-VY*JACKSON.dat"))
    sea = load(os.path.join(d, "KuboBastinSea_*JACKSON.dat"))
    srf = load(os.path.join(d, "KuboBastinSurf_*JACKSON.dat"))
    fft = load(os.path.join(d, "KuboBastin_*JACKSON_FFTgrid.dat"))
    ok = True
    scale = float(np.max(np.abs(bas[:, 1])))
    same_grid = np.allclose(bas[:, 0], sea[:, 0]) and np.allclose(bas[:, 0], srf[:, 0])
    diff = float(np.max(np.abs(sea[:, 1] + srf[:, 1] - bas[:, 1]))) / scale
    print(f"grids: bastin/sea/surf {len(bas)} points (same grid: {same_grid}), bastin_fft {len(fft)} points")
    print(f"[1] max |sea + surf - bastin| / max |bastin| = {diff:.3e}")
    good1 = same_grid and diff < 1e-2
    inside = (fft[:, 0] > bas[0, 0]) & (fft[:, 0] < bas[-1, 0])
    dfft = float(np.max(np.abs(np.interp(fft[inside, 0], bas[:, 0], bas[:, 1]) - fft[inside, 1]))) / scale
    print(f"    (information) max |bastin_fft - bastin| / max |bastin| = {dfft:.3e} over "
          f"{int(inside.sum())} FFT nodes; bastin_fft at the top of the band = {fft[-1, 1]:.4f}")
    c_b, c_sea, c_srf = chern(*bas.T), chern(*sea.T), chern(*srf.T)
    c_fft = chern(*fft.T)
    print(f"[2] gap |E| < 0.5 eV: C(bastin_fft) = {c_fft:+.4f} vs C(bastin) = {c_b:+.4f}")
    good2 = abs(c_fft - c_b) < 0.02
    print(f"[3] gap |E| < 0.5 eV: C(bastin) = {c_b:+.4f}, C(sea) = {c_sea:+.4f}, C(surf) = {c_srf:+.4f}")
    good3 = abs(c_b - 1) < 0.02 and abs(c_sea - 1) < 0.02 and abs(c_srf) < 0.02
    for tag, g in (("1", good1), ("2", good2), ("3", good3)):
        print(f"    check {tag}: {'PASS' if g else 'FAIL'}")
    ok = good1 and good2 and good3
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
