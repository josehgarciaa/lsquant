#!/usr/bin/env python3
"""Check saved high-level API output against the existing low-level API example.

Only level2() is recomputed: two exact Bloch traces of 2x2 graphene blocks.
Level1 files are read, never overwritten. The moment tolerance is the example's
1e-12 round-off allowance; stored text curves require a 1e-5 relative budget
because their energies and values are printed with six significant digits.
The independent analytic-band DOS moments also guard against shared API errors.
Pitfall: the extensive DOS must be divided by 512 orbitals before quoting it.
Work and storage are bounded by 128 moments and the small reconstruction grids.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.dont_write_bytecode = True
from check_support import one, relative_error, run, table
import numpy as np


def check(root, c, args):
    """Compare both API levels and the quoted per-orbital DOS landmarks."""
    from python_api import level2, _core, BROADENING_MEV, dos_curve, kubo_greenwood
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "10_bases_compare"))
    from compare_bases import analytic_dos_moments
    low = level2()
    for key, prefix, cv_prefix, dims, points in (
            ("dos", "SpectralOp1", "meanSpectralOp1", [128], 3840),
            ("kubo", "NonEqOpVX-VX", "KuboGreenWood_NonEqOpVX-VX", [128, 128], 127)):
        saved = _core.ReadChebmom(str(one(root, prefix + "*.chebmom")))
        c.require(f"{key} moments shape", list(saved.dims) == dims)
        c.near(f"{key} API levels relative moment error",
               relative_error(saved.data, low[key][0].data), 0, 1e-12,
               "same exact trace, round-off allowance from existing example")
        reconstruct = dos_curve if key == "dos" else kubo_greenwood
        e1, v1 = reconstruct(saved, BROADENING_MEV)
        e2, v2 = low[key][1]
        c.require(f"{key} reconstruction coordinates", np.array_equal(e1, e2))
        c.near(f"{key} API levels relative curve error", relative_error(v1, v2), 0, 1e-12,
               "reconstruction from identical moments")
        stored = table(one(root, cv_prefix + "*JACKSON.dat"))
        c.near(f"{key} curve point count", len(stored), points, 0, "published grid size")
        c.near(f"{key} stored coordinate relative error", relative_error(stored[:, 0], e2),
               0, 1e-5, "six significant digit text output")
        c.near(f"{key} stored value relative error", relative_error(stored[:, 1], v2),
               0, 1e-5, "six significant digit text output")
        if key == "dos":
            ref = analytic_dos_moments(saved.window, 128)
            c.near("analytic-band DOS moments", np.max(np.abs(np.asarray(saved.data).ravel() - ref)),
                   0, 1e-10, "independent analytic bands, exact trace round-off")
    e, rho = low["dos"][1]
    rho = rho / 512
    c.near("Dirac-point DOS per orbital (1/eV)", np.interp(0, e, rho), 0.0001, 0.0001,
           "four decimal places and Jackson resolution")
    c.near("van Hove DOS per orbital (1/eV)", np.interp(2.7, e, rho), 0.1773, 0.001,
           "finite mesh and 150 meV Jackson resolution")


if __name__ == "__main__":
    raise SystemExit(run(check))
