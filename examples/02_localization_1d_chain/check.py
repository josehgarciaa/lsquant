#!/usr/bin/env python3
"""Check disorder-averaged late-time spreading against the published calibration.

The saved MSD is averaged over 16 disorder samples before taking its square
root. The final 30 percent of times define the plateau, just as in plot.py.
The plotted localization length is calibrated to Born theory at E=1 eV;
raw sqrt(MSD) is not the Lyapunov localization length. All lengths are sites.
Work is linear in 192 short curves, without time evolution or diagonalisation.
Warning: the centre and edge have known departures from Born theory; neither
the Kappus-Wegner value nor the edge asymptote is a precise gate here.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.dont_write_bytecode = True
from check_support import one, run, same_grid, table
import numpy as np


def check(root, c, args):
    """Check all sizes and energies, using raw plateaus before Born calibration."""
    from plot import _load_averaged_msd, xi_from_plateau, xi_born
    energies = [0.0, 1.0, -1.0, 1.9]
    xi = {}
    for n in (256, 1024, 4096):
        for e in energies:
            curves = [table(one(root / f"r{r:02d}",
                      f"mean*chain1d_dis_N{n}_W1*EF{e:.6f}*.dat")) for r in range(1, 17)]
            same_grid(*curves)
            c.near(f"N={n} E={e:g} final time (fs)", curves[0][-1, 0], 600, 0.01,
                   "printed time resolution")
        t, msd = _load_averaged_msd(str(root), n, 1.0, energies, 16)
        xi[n] = {e: xi_from_plateau(t, msd[e]) for e in energies}
        for e in energies:
            c.require(f"N={n} E={e:g} physical MSD", np.isfinite(msd[e]).all() and
                      np.mean(msd[e][t >= 0.7 * t.max()]) > 0)
    for e, ref in ((0.0, 49.0), (1.0, 43.0), (1.9, 13.3)):
        c.near(f"N=4096 E={e:g} sqrt(MSD) (sites)", xi[4096][e], ref, 0.20 * ref,
               "16 disorder samples, one random trace vector and finite-time plateau")
    c.near("E=+/-1 plateau ratio", xi[4096][-1.0] / xi[4096][1.0], 1, 0.25,
           "finite disorder ensemble is symmetric only on average")
    scale = float(xi_born(1.0, 1.0)) / xi[4096][1.0]
    c.near("calibrated centre length (sites)", scale * xi[4096][0.0], 96, 24,
           "Born trend; centre anomaly and disorder scatter")
    c.require("near-edge length below band-interior length", xi[4096][1.9] < xi[4096][1.0])
    for n in (256, 1024):
        for e in energies:
            c.near(f"N={n} E={e:g} length / largest-sample length",
                   xi[n][e] / xi[4096][e], 1, 0.6,
                   "finite-size suppression and disorder scatter; no strict monotonicity assumed")


if __name__ == "__main__":
    raise SystemExit(run(check, "run_localization"))
