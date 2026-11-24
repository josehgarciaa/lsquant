#!/usr/bin/env python3
"""Check basis invariance and analytic bands without recomputing any traces.

Reuse compare_bases.py's loaders and analytic-band oracle. All arrays, including
complex moments and curve coordinates, are compared. The 1e-10 moment budget
is the existing comparison's round-off allowance (reference 1.3e-14).
Curves must remain identical as written. Work is linear in three small rank-2
moment tensors, plus 128 times 512 analytic band values.
Pitfall: all bases must use the same explicit [-9,9] eV spectral window.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.dont_write_bytecode = True
from check_support import one, relative_error, run, same_grid, table
import numpy as np


def check(root, c, args):
    """Check saved traces in three bases and the independent graphene spectrum."""
    from compare_bases import BASES, load, curve, analytic_dos_moments
    for prefix, cv in (("SpectralOp1", "meanSpectralOp1"),
                       ("NonEqOpVX-VX", "KuboGreenWood_NonEqOpVX-VX")):
        mus, cvs = {}, {}
        for basis in BASES:
            one(root / f"bases_{basis}", f"{prefix}*.chebmom")
            table(one(root / f"bases_{basis}", f"{cv}*JACKSON.dat"))
            mu = load(str(root), basis, prefix)
            mus[basis] = np.asarray(mu.data)
            cvs[basis] = curve(str(root), basis, cv)
            c.require(f"{basis} {prefix} window", mu.window.lower == -9 and mu.window.upper == 9)
            expected_dims = [128] if prefix == "SpectralOp1" else [128, 128]
            c.require(f"{basis} {prefix} dimensions", list(mu.dims) == expected_dims)
        same_grid(*cvs.values())
        for basis in ("bloch", "eigen"):
            c.near(f"{prefix} {basis} vs real relative moment error",
                   relative_error(mus[basis], mus["real"]), 0, 1e-10,
                   "exact traces in different bases, floating-point summation")
            c.near(f"{prefix} {basis} vs real max stored curve error",
                   np.max(np.abs(cvs[basis] - cvs["real"])), 0, 0,
                   "README curves identical as written")
    mu = load(str(root), "real", "SpectralOp1")
    ref = analytic_dos_moments(mu.window, mu.dims[0])
    c.near("analytic-band max complex DOS moment error",
           np.max(np.abs(np.asarray(mu.data).ravel() - ref)), 0, 1e-10,
           "analytic spectrum and exact trace round-off")


if __name__ == "__main__":
    raise SystemExit(run(check))
