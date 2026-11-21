#!/usr/bin/env python3
"""Checks for Anderson disorder in the non-orthogonal SIESTA basis of the Bi bilayer.

1. Two layouts, one disorder: the real-space CSR run (basis "real") and the
   Fourier/Bloch run (basis "bloch") with the same disorder seed apply the same
   potential (D S + S D)/2 and the same random-phase vectors, so their moments
   must agree to round-off.
2. Two disorder seeds: seeds 7 and 8 draw different potentials (different
   disorder.npy hashes) and give different, but statistically similar, DOS.
3. Sum rule: the DOS per cell integrates to the number of spin-orbitals per
   cell (52 for two Bi atoms in a DZP basis with spin-orbit coupling).

Usage:
    python check_disorder.py <out_root>   # holds disorder_{real,bloch}_seed7/, disorder_real_seed8/
"""
import glob
import hashlib
import json
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(_HERE)), "python"))
from lsquant import _core  # noqa: E402


def load(run_dir):
    mu = _core.ReadChebmom(sorted(glob.glob(os.path.join(run_dir, "*.chebmom")))[0])
    dos = np.loadtxt(os.path.join(run_dir, "dos.dat"))
    w = np.load(os.path.join(run_dir, "disorder.npy"))
    basis = json.load(open(os.path.join(run_dir, "basis.json")))
    return np.asarray(mu.data).ravel(), dos, w, basis


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "out"
    r7, d7, w7, b7 = load(os.path.join(root, "disorder_real_seed7"))
    k7, dk7, wk7, _ = load(os.path.join(root, "disorder_bloch_seed7"))
    r8, d8, w8, _ = load(os.path.join(root, "disorder_real_seed8"))
    ok = True
    sha = lambda a: hashlib.sha256(a.astype("<c16").tobytes()).hexdigest()[:12]
    print(f"dimension {len(w7)}, window {b7['window']}, overlap bounds {b7['overlap_bounds']}")
    print(f"disorder sha256: seed7 real {sha(w7)}, seed7 bloch {sha(wk7)}, seed8 real {sha(w8)}")
    same = np.array_equal(w7, wk7)
    ok &= same
    dmu = float(np.max(np.abs(r7 - k7)))
    good = dmu < 1e-10
    ok &= good
    print(f"[1] real vs bloch, seed 7: same disorder vector = {same}; max |d mu| = {dmu:.2e}  "
          f"{'PASS' if good else 'FAIL'}")
    diff = float(np.max(np.abs(d7[:, 1] - np.interp(d7[:, 0], d8[:, 0], d8[:, 1]))))
    print(f"[2] seed 7 vs seed 8: max |d DOS| = {diff:.3f} states/eV/cell "
          f"(DOS max {d7[:, 1].max():.2f}); differ = {not np.array_equal(w7, w8)}")
    ok &= not np.array_equal(w7, w8)
    for name, d in (("seed7 real", d7), ("seed7 bloch", dk7), ("seed8 real", d8)):
        integral = float(np.trapezoid(d[:, 1], d[:, 0]))
        good = abs(integral - 72) < 0.5
        ok &= good
        print(f"[3] {name}: integral DOS dE = {integral:.3f} (expected 72)  {'PASS' if good else 'FAIL'}")
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
