#!/usr/bin/env python3
"""Check the curves rebuilt by ``lsquant reconstruct`` from stored moments.

1. Reproducibility: reconstructing the N = 256 moments at 120 meV gives, bit
   for bit, the 120 meV curve that ``lsquant run`` wrote.
2. Analytic limit: the reconstructed curve covers |E| <= 0.95 * 2 eV (the
   reconstruction grid stops at 95 % of the half-width). For the long chain
   (N = 4096) the weight on that interval must approach the bulk value
   (2/pi) arcsin(0.95) = 0.798 per site.
3. Broadening: the 60 meV curve of N = 256 is smoother than the 20 meV one
   (smaller maximum), with the same weight.

Usage:
    python check_reconstruct.py <run_dir_of_example_01> <reconstruct_out_dir>
"""
import filecmp
import os
import sys

import numpy as np


def weight(path, n_sites):
    c = np.loadtxt(path)
    return float(np.trapezoid(c[:, 1], c[:, 0]) / n_sites), float(c[:, 1].max() / n_sites)


def main():
    run_dir, rec_dir = sys.argv[1], sys.argv[2]
    stem = "meanSpectralOp1chain1d_N{n}KPM_M512_stateexact_chain1d_N{n}JACKSON.dat"
    ok = True
    same = filecmp.cmp(os.path.join(run_dir, stem.format(n=256)),
                       os.path.join(rec_dir, "eta120", stem.format(n=256)), shallow=False)
    ok &= same
    print(f"[1] reconstruct at 120 meV == curve written by run: {same}  {'PASS' if same else 'FAIL'}")
    bulk = 2.0 / np.pi * np.arcsin(0.95)
    w4096, _ = weight(os.path.join(rec_dir, "dir_eta60", stem.format(n=4096)), 4096)
    good = abs(w4096 - bulk) < 0.02
    ok &= good
    print(f"[2] N=4096, 60 meV: weight on |E|<1.9 eV = {w4096:.4f} per site "
          f"(bulk {bulk:.4f})  {'PASS' if good else 'FAIL'}")
    w20, p20 = weight(os.path.join(run_dir, "eta20", "meanmoments_chain1d_N256_M512JACKSON.dat"), 256)
    w60, p60 = weight(os.path.join(rec_dir, "eta60", stem.format(n=256)), 256)
    good = p60 < p20 and abs(w60 - w20) < 0.02
    ok &= good
    print(f"[3] N=256: max rho 20 meV = {p20:.3f}, 60 meV = {p60:.3f} 1/eV; "
          f"weights {w20:.4f}, {w60:.4f}  {'PASS' if good else 'FAIL'}")
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
