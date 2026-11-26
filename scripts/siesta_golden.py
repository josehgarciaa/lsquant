#!/usr/bin/env python3
"""Reduced goldens for the SIESTA run files (G11) of the stage-7 gate.

Problem: the 24 ``runs/goldens/g11_*.json`` run files (Bi(111) bilayer and PdSe2
monolayer, three velocity modes, three metric routes, density of states and
spin Hall conductivity) were only checked to run. Their outputs cannot be
committed whole: one curve is 1.2 MB and a rank-2 moment file 16 MB. They also
cannot be byte-identical between machines, because the Bloch Hamiltonians are
diagonalised and orthogonalised by LAPACK, whose rounding depends on the BLAS
library in use.

Approach: a golden is a small JSON file per run file with
  * every two-column curve the run wrote (``*.dat``), sampled on at most
    ``MAX_ROWS`` evenly spaced rows, with the number of rows, the largest
    absolute value and the trapezoid integral of the whole curve;
  * the top-level numbers and numeric lists of the exact-oracle files
    (``spin_chern.json``, ``exact_shc.json``), the lists sampled the same way.
``check`` recomputes the same reduction from a fresh run and compares with a
relative tolerance, each quantity relative to the largest absolute value of its
curve (so the near-zero tails of a curve do not dominate). There are two
tolerances, set from the difference measured between two machines with
different BLAS libraries: the reconstructed density of states agreed to 8e-14
and the reconstructed spin Hall curve to 7e-10 (tolerance 1e-7); the exact
oracle, which sums over nearly degenerate bands, agreed only to 7e-7
(tolerance 1e-5).

Usage:
    siesta_golden.py make  <run output dir> <golden.json>
    siesta_golden.py check <run output dir> <golden.json> [rtol_curves [rtol_oracle]]
``check`` prints one line with the largest deviation and exits 0 (within
tolerance) or 1.

Pitfall: integers (Z2, number of occupied bands) are compared exactly; a
quantity present on one side only is a failure, never skipped.
To-do: the moments themselves are not pinned, only what is reconstructed from them.
"""
import glob
import json
import os
import sys

import numpy as np

MAX_ROWS = 256
DEFAULT_RTOL = 1e-7          # reconstructed curves (*.dat)
DEFAULT_ORACLE_RTOL = 1e-5   # numbers and lists of the exact-oracle files
ORACLE_FILES = ("spin_chern.json", "exact_shc.json")


def _sample(values):
    values = np.asarray(values, dtype=float)
    stride = max(1, int(np.ceil(values.shape[0] / MAX_ROWS)))
    return stride, values[::stride]


def reduce_output(out_dir):
    """Return the reduction of one run output directory as a plain dict."""
    curves = {}
    for path in sorted(glob.glob(os.path.join(out_dir, "*.dat"))):
        data = np.loadtxt(path, ndmin=2)
        if data.shape[1] != 2:
            continue
        stride, rows = _sample(data)
        curves[os.path.basename(path)] = {
            "rows": int(data.shape[0]), "stride": stride,
            "max_abs": float(np.max(np.abs(data[:, 1]))),
            "integral": float(np.sum(0.5 * (data[1:, 1] + data[:-1, 1]) * np.diff(data[:, 0]))),
            "x": rows[:, 0].tolist(), "y": rows[:, 1].tolist(),
        }
    scalars, lists = {}, {}
    for name in ORACLE_FILES:
        path = os.path.join(out_dir, name)
        if not os.path.isfile(path):
            continue
        with open(path, encoding="utf-8") as fh:
            content = json.load(fh)
        for key, value in sorted(content.items()):
            label = f"{name}:{key}"
            if isinstance(value, bool):
                continue
            if isinstance(value, (int, float)):
                scalars[label] = value
            elif (isinstance(value, list) and len(value) > 1
                  and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in value)):
                stride, rows = _sample(value)
                lists[label] = {"rows": len(value), "stride": stride,
                                "max_abs": float(np.max(np.abs(value))), "y": rows.tolist()}
    return {"curves": curves, "scalars": scalars, "lists": lists}


def compare(new, gold):
    """Return ({"curves": (dev, where), "oracle": (dev, where)}, structural failures)."""
    worst, problems = {"curves": (0.0, "-"), "oracle": (0.0, "-")}, []
    group = "curves"

    def note(deviation, label):
        if deviation > worst[group][0]:
            worst[group] = (float(deviation), label)

    for kind in ("curves", "lists"):
        if sorted(new[kind]) != sorted(gold[kind]):
            problems.append(f"{kind}: {sorted(new[kind])} != golden {sorted(gold[kind])}")
            continue
        group = "curves" if kind == "curves" else "oracle"
        for name, ref in gold[kind].items():
            got = new[kind][name]
            if got["rows"] != ref["rows"]:
                problems.append(f"{name}: {got['rows']} rows != golden {ref['rows']}")
                continue
            scale = max(ref["max_abs"], 1e-300)
            note(np.max(np.abs(np.array(got["y"]) - np.array(ref["y"]))) / scale, f"{name} values")
            note(abs(got["max_abs"] - ref["max_abs"]) / scale, f"{name} maximum")
            if kind == "curves":
                span = max(abs(ref["x"][-1] - ref["x"][0]), 1e-300)
                note(np.max(np.abs(np.array(got["x"]) - np.array(ref["x"]))) / span, f"{name} grid")
                note(abs(got["integral"] - ref["integral"]) / (scale * span), f"{name} integral")
    group = "oracle"
    if sorted(new["scalars"]) != sorted(gold["scalars"]):
        problems.append(f"scalars: {sorted(new['scalars'])} != golden {sorted(gold['scalars'])}")
    else:
        for name, ref in gold["scalars"].items():
            got = new["scalars"][name]
            if isinstance(ref, int) and isinstance(got, int):
                if got != ref:
                    problems.append(f"{name}: {got} != golden {ref}")
            else:
                note(abs(got - ref) / max(abs(ref), 1.0), name)
    return worst, problems


def main(argv):
    if len(argv) not in (4, 5, 6) or argv[1] not in ("make", "check"):
        print(__doc__)
        return 2
    out_dir, golden = argv[2], argv[3]
    reduction = reduce_output(out_dir)
    if not reduction["curves"]:
        print(f"no two-column curve in {out_dir}")
        return 1
    if argv[1] == "make":
        with open(golden, "w", encoding="utf-8") as fh:
            json.dump(reduction, fh, indent=0, sort_keys=True)
            fh.write("\n")
        print(f"wrote {golden}: {len(reduction['curves'])} curves, {len(reduction['scalars'])} numbers, "
              f"{len(reduction['lists'])} lists")
        return 0
    limits = {"curves": float(argv[4]) if len(argv) >= 5 else DEFAULT_RTOL,
              "oracle": float(argv[5]) if len(argv) == 6 else DEFAULT_ORACLE_RTOL}
    with open(golden, encoding="utf-8") as fh:
        gold = json.load(fh)
    worst, problems = compare(reduction, gold)
    for line in problems:
        print("structure:", line)
    print("; ".join(f"{group} max rel dev {worst[group][0]:.3e} at {worst[group][1]} (rtol {limits[group]:g})"
                    for group in ("curves", "oracle")))
    return 0 if not problems and all(worst[g][0] <= limits[g] for g in limits) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
