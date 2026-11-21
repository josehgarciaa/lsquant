#!/usr/bin/env python3
"""Chern number of the golden Haldane model from the Kubo-Bastin curves of run.sh.

C = 2 pi <sigma_xy>_{|E| < 0.5 eV} / A, with A = (N/2) 3 sqrt(3)/2 the area of the
512-site supercell (C-C distance 1). Reads every KuboBastin_NonEqOpVX-VY*.dat (sigma_xy) under the given output
directories and checks |C - 1| < 0.05 for each (the plateau rounds off at the gap
edges as the broadening grows).
usage: python3 examples/run_files/chern.py <out_dir> [<out_dir> ...]
"""
import glob
import math
import os
import sys

import numpy as np

AREA = (512 / 2) * 1.5 * math.sqrt(3.0)


def main(dirs):
    ok = True
    for d in dirs:
        for f in sorted(glob.glob(os.path.join(d, "**", "KuboBastin_NonEqOpVX-VY*.dat"), recursive=True)):
            e, s = np.loadtxt(f, unpack=True)[:2]
            c = 2.0 * math.pi * float(np.mean(s[np.abs(e) < 0.5])) / AREA
            ok &= abs(c - 1.0) < 0.05
            print(f"{os.path.relpath(f)}: C = {c:+.4f}")
    print("checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
