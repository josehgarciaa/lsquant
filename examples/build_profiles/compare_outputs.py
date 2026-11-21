#!/usr/bin/env python3
"""Compare the files of two run outputs of the same run file (golden vs native build).

Moment files (*.chebmom*) are compared after their two header lines (the native profile
tags the first header line with "# build_profile=native"); curves (*.dat) on every column
after the energy column. Each difference is reported relative to the largest value of the
same kind in the run (all moment files, or all curves), so a curve that is itself noise,
such as the Fermi-surface part of a Hall conductivity, is judged on the scale of the run.
Exits 0 when every file is identical or within 1e-10 of that scale.
usage: python3 examples/build_profiles/compare_outputs.py <dir_golden> <dir_native>
"""
import filecmp
import os
import sys

import numpy as np


def load(path):
    skip = 2 if ".chebmom" in path else 0
    rows = []
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i < skip or line.lstrip().startswith("#"):
                continue
            try:
                rows.append([float(x) for x in line.split()])
            except ValueError:
                continue
    width = max((len(r) for r in rows), default=0)
    data = np.array([r for r in rows if len(r) == width])
    return data[:, 1:] if path.endswith(".dat") and width > 1 else data


def main(a, b):
    names = [n for n in sorted(set(os.listdir(a)) & set(os.listdir(b)))
             if os.path.isfile(os.path.join(a, n)) and (n.endswith(".dat") or ".chebmom" in n)]
    data = {n: (load(os.path.join(a, n)), load(os.path.join(b, n))) for n in names}
    scale = {}
    for n, (va, _) in data.items():
        kind = "curve" if n.endswith(".dat") else "moments"
        scale[kind] = max(scale.get(kind, 0.0), float(np.max(np.abs(va))) if va.size else 0.0)
    ok, worst = True, 0.0
    for n, (va, vb) in data.items():
        kind = "curve" if n.endswith(".dat") else "moments"
        if filecmp.cmp(os.path.join(a, n), os.path.join(b, n), shallow=False):
            print(f"identical  {n}")
            continue
        if va.shape != vb.shape or va.size == 0:
            print(f"DIFFERENT  {n} (layout)")
            ok = False
            continue
        rel = float(np.max(np.abs(va - vb))) / (scale[kind] or 1.0)
        worst = max(worst, rel)
        ok &= rel < 1e-10
        label = "curves" if kind == "curve" else "moments"
        print(f"{'within' if rel < 1e-10 else 'BEYOND'}     {n}: max |golden - native| / max |{label} of the run| = {rel:.2e}")
    print(f"largest relative difference: {worst:.2e}")
    print("checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
