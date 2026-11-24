#!/usr/bin/env python3
"""Check diffusion, Einstein conductivity and mean free path from saved curves.

The x-projected diffusion convention is MSD/(2*t), with t>0. The DOS is
extensive on 200^2 sites; sigma=rho*D and ell=2*D/v_F use the plotting
script's conventions. These illustrative random-vector results need a 20%
budget for finite time, finite size and trace noise, rather than exact goldens.
Work is linear in the eight curves plus the existing small clean-band mesh.
Pitfall: the weak-disorder maximum is limited by the 150 fs observation time.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.dont_write_bytecode = True
from check_support import one, run, table
import numpy as np


def check(root, c, args):
    """Recompute transport quantities, using independently quoted calibration."""
    from plot import fermi_velocity
    vf = fermi_velocity(2.0)
    # The clean contour is only a reference for a disordered metal.
    c.near("clean contour v_F (eV lattice units)", vf, 2.2, 0.2,
           "README estimate and finite contour shell")
    maxima, conductivities, lengths = [], [], []
    for w, ref_d, ref_s in ((0.5, 21.9, 2.49), (1, 6.9, 0.77),
                            (2, 1.7, 0.19), (4, 0.39, 0.044)):
        label = f"square_L200_W{w:g}"
        paths = list(root.glob(f"meanCorrelation*{label}*EF2.000000*.dat"))
        paths += list(root.glob(f"meanMSD*{label}*EF2.000000*.dat"))
        if len(paths) != 1:
            raise ValueError(f"expected one MSD for {label}, found {len(paths)}")
        msd = table(paths[0])
        t, dx2 = msd.T
        positive = t > 0
        diffusion = dx2[positive] / (2 * t[positive])
        dmax = float(diffusion.max())
        dos = table(one(root, f"meanSpectralOp1{label}*JACKSON.dat"))
        c.require(f"W={w:g} DOS covers E=2", dos[0, 0] <= 2 <= dos[-1, 0])
        rho = float(np.interp(2, dos[:, 0], dos[:, 1])) / 40000
        sigma, ell = rho * dmax, 2 * dmax / vf
        reason = "single random-vector trace, finite size and finite-time maximum"
        c.near(f"W={w:g} D_max", dmax, ref_d, 0.2 * ref_d, reason)
        c.near(f"W={w:g} sigma", sigma, ref_s, 0.2 * ref_s, reason)
        c.near(f"W={w:g} DOS per site (1/eV)", rho, 0.109, 0.025,
               "clean analytic DOS, modified by disorder and broadening")
        print(f"W={w:g} ell = {ell:.6g} lattice spacings")
        maxima.append(dmax)
        conductivities.append(sigma)
        lengths.append(ell)
        if w == 4:
            c.near("W=4 late D / D_max", diffusion[-1] / dmax, 0.8, 0.15,
                   "approximate late-time decay and stochastic trace")
    for label, values in (("D_max", maxima), ("sigma", conductivities), ("ell", lengths)):
        c.require(f"{label} decreases with disorder", np.all(np.diff(values) < 0))
    for index, target, tol in ((0, 20, 5), (1, 6, 2), (3, 0.35, 0.12)):
        c.near(f"ell at sweep index {index} (sites)", lengths[index], target, tol,
               "rounded README estimate and noisy diffusion maximum")
    c.require("strong disorder passes Ioffe-Regel limit", lengths[-1] < 1)


if __name__ == "__main__":
    raise SystemExit(run(check, "run_transport"))
