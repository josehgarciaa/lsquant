#!/usr/bin/env python3
"""Check quantised Hall transport and the sea/surface decomposition.

Stored conductivities are extensive in e=hbar=1 units: multiply by 2*pi/A
to obtain e^2/h. The gap mask is |E|<0.5 eV. The finite 512-state cell and
128-moment reconstruction shift the plateau to about 1.017, so a 0.03
absolute Chern tolerance is physical; round-off tolerances are inappropriate.
The existing small Bloch calculation independently checks the topological phase.
Pitfall: this does not rerun transport at other sublattice masses.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.dont_write_bytecode = True
from check_support import one, run, same_grid, table
import numpy as np


def check(root, c, args):
    """Compare the gap plateau to the lower-band topological integer."""
    from plot import A_CELL, chern_from_plateau, chern_lower
    curves = [table(one(root, p + "*JACKSON.dat")) for p in (
        "KuboBastin_NonEqOpVX-VY", "KuboBastinSea_NonEqOpVX-VY",
        "KuboBastinSurf_NonEqOpVX-VY", "KuboGreenWood_NonEqOpVX-VX")]
    same_grid(*curves[:3])
    norm = 2 * np.pi / (256 * A_CELL)
    for name, d, target in zip(("total", "sea", "surface", "longitudinal"), curves, (1, 1, 0, 0)):
        mask = np.abs(d[:, 0]) < 0.5
        c.require(f"{name} gap sampled", mask.sum() >= 3)
        c.near(f"{name} mean gap response (e^2/h)",
               chern_from_plateau(512, d[:, 0], d[:, 1]), target, 0.03,
               "finite cell and 128-moment resolution; calibrated total 1.017")
        c.near(f"{name} max gap deviation (e^2/h)",
               np.max(np.abs(d[mask, 1] * norm - target)), 0, 0.05,
               "finite-resolution gap leakage and plateau flatness")
    bas, sea, surf = curves[:3]
    c.near("max |sea+surface-total| / max |total|",
           np.max(np.abs(sea[:, 1] + surf[:, 1] - bas[:, 1])) / np.max(np.abs(bas[:, 1])),
           0, 0.01, "independent finite-grid reconstructions and printed curves")
    for mass, expected in ((0.0, 1), (1.6, 0)):
        c.near(f"Bloch Chern at mass={mass:g} eV", chern_lower(mass=mass), expected,
               1e-10, "integer lattice Berry flux, away from gap closing")
    c.near("critical mass (eV)", 3 * np.sqrt(3) * 0.15, 0.78, 0.005,
           "README value rounded to two decimals")


if __name__ == "__main__":
    raise SystemExit(run(check, "run_haldane"))
