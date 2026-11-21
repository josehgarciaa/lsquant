#!/usr/bin/env python3
"""Checks for the spectral projection (shift) of the Bi bilayer.

1. Nothing inside the window is lost: the exact (Gaussian) DOS of the
   projected and the unprojected matrices agree on the oracle grid
   (|E| <= 1.5 eV), as recorded by the run in word_0_oracle.json.
2. The projected spectrum fits in the window: the spectral bounds after
   projection are the window edges, against the full SIESTA range before.
3. Sum rules: the identity word Tr[delta(E-H)] still counts all 72
   spin-orbitals per cell (the states outside the window are moved to the
   nearest window edge, not removed); the projected identity
   Tr[WINDOW delta(E-H)] counts only the kept bands, whose average number per
   k point is in projection.json.
4. Inside the window, away from the edges, the two words give the same DOS.

Usage:
    python check_projection.py <run_out_dir>
"""
import json
import os
import sys

import numpy as np


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else "out/projection_shift"
    proj = json.load(open(os.path.join(d, "projection.json")))
    ok = True
    for w in proj["words"]:
        err = w["oracle_error"]["max_relative"]
        good = err < 1e-10
        ok &= good
        print(f"[1] {w['word']:24s} exact DOS projected vs unprojected on |E| <= 1.5 eV: "
              f"max rel diff = {err:.2e}  {'PASS' if good else 'FAIL'}")
    lo, hi = proj["window"]
    sb, ob = proj["spectral_bounds"], proj["original_bounds"]
    good = abs(sb[0] - lo) < 1e-9 and abs(sb[1] - hi) < 1e-9
    ok &= good
    print(f"[2] spectrum before [{ob[0]:.2f}, {ob[1]:.2f}] eV, after [{sb[0]:.6f}, {sb[1]:.6f}] eV  "
          f"{'PASS' if good else 'FAIL'}")
    dos = np.loadtxt(os.path.join(d, "word_0.dat"))
    win = np.loadtxt(os.path.join(d, "word_1.dat"))
    kept = float(np.mean(proj["kept_bands_per_k"]))
    i0 = float(np.trapezoid(dos[:, 1], dos[:, 0]))
    i1 = float(np.trapezoid(win[:, 1], win[:, 0]))
    good = abs(i0 - 72) < 0.05 and abs(i1 - kept) < 0.01
    ok &= good
    print(f"[3] int Tr[delta] dE = {i0:.3f} (72 spin-orbitals); int Tr[WINDOW delta] dE = {i1:.4f} "
          f"(mean kept bands {kept:.5f})  {'PASS' if good else 'FAIL'}")
    m = np.abs(dos[:, 0]) < 1.5
    diff = float(np.max(np.abs(dos[m, 1] - win[m, 1])))
    good = diff < 0.01 * dos[m, 1].max() + 1e-12
    ok &= good
    print(f"[4] |E| < 1.5 eV: max |DOS - WINDOW DOS| = {diff:.4f} (max DOS {dos[m, 1].max():.3f}) "
          f"states/eV/cell  {'PASS' if good else 'FAIL'}")
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
