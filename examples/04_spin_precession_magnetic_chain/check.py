#!/usr/bin/env python3
"""Check saved Larmor precession using the same analytic gate as plot.py.

No model or figure is generated. The three normalized spin components,
unwrapped-phase frequency and transverse amplitude have the existing 1e-10
tolerance. Pitfall: an incomplete or non-finite time series must fail before
normalization. The full stored header and all m=0 rows are validated here.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.dont_write_bytecode = True
from check_support import one, run
import numpy as np


def check(root, c, args):
    """Validate input and call the shared pointwise and frequency comparison."""
    from precession import NTIME, HBAR, J_EX, _m0, check_precession
    for op in ("SX", "SY", "SZ"):
        path = one(root, f"TimeProj{op}-PX*stateexact*chebmomTD")
        with path.open() as stream:
            header = stream.readline().split()
            dims = stream.readline().split()
        c.require(f"{op} published dimensions", int(header[0]) == 128 and
                  list(map(int, dims)) == [32, NTIME])
        data, nt = _m0(path)
        c.require(f"{op} finite m=0 spin", len(data) == NTIME and np.isfinite(data).all())
        if op == "SX":
            c.require("nonzero initial SX", abs(data[0]) > 1e-12)
    if not c.ok:
        return
    print("Expected omega=2*J/hbar, SX=cos(omega*t), SY=-sin(omega*t), SZ=0, amplitude=1;")
    print("tolerance=1e-10 for relative frequency and absolute normalized spin errors (exact-trace round-off).")
    c.require("shared analytic precession gate", check_precession(str(root)))
    c.near("precession period (fs)", 2 * np.pi / (2 * J_EX / HBAR), 20.7, 0.1,
           "README period rounded to one decimal")


if __name__ == "__main__":
    raise SystemExit(run(check, "run_spin"))
