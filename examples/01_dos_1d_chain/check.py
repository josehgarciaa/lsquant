#!/usr/bin/env python3
"""Check intensive chain moments, finite-size revival and the analytic bulk DOS.

The exact trace has only mu_0=1 and mu_m=2 at positive multiples of N.
Curve files are extensive and must be divided by N. Only |E|<1.9 eV is
compared to the bulk DOS, away from its divergent edges. Complexity is linear
in the saved moments and curves; no Hamiltonian is constructed.
Pitfall: the 10 meV request retains at most the 512 available moments.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.dont_write_bytecode = True
from check_support import one, run, table
import numpy as np


def check(root, c, args):
    """Compare the published chain outputs with the closed-form trace and DOS."""
    from plot import read_moments, m_eff
    integrate = getattr(np, "trapezoid", None)
    if integrate is None:
        integrate = np.trapz
    for n in (256, 1024, 4096):
        d = read_moments(one(root, f"moments_chain1d_N{n}_M512.chebmom1D"))
        c.require(f"N={n} header", d["size"] == n and d["M"] == 512 and
                  d["band_width"] == 4.0 and d["band_center"] == 0.0)
        expected = np.zeros(512)
        expected[0] = 1.0
        expected[n::n] = 2.0
        c.near(f"N={n} max moment error", np.max(np.abs(d["mu"] - expected)),
               0, 1e-10, "exact trace round-off, including revival at m=N")
        data = table(one(root / "eta20", f"meanmoments_chain1d_N{n}_M512JACKSON.dat"))
        mask = np.abs(data[:, 0]) < 1.9
        c.require(f"N={n} bulk samples", mask.sum() > 100)
        if n >= 1024:
            ref = 1 / (np.pi * np.sqrt(4 - data[mask, 0] ** 2))
            c.near(f"N={n} max DOS error (1/eV)",
                   np.max(np.abs(data[mask, 1] / n - ref)), 0, 5e-5,
                   "printed energy/DOS resolution; reference error 1.2e-5")
    for eta, expected in ((20, 315), (80, 79), (24, 262), (10, 512)):
        c.near(f"eta={eta} retained moments", min(512, m_eff(eta, 4)), expected,
               0, "Jackson resolution, capped by available moments")
        data = table(one(root / f"eta{eta}", "meanmoments_chain1d_N256_M512JACKSON.dat"))
        c.require(f"eta={eta} DOS nonnegative", data[:, 1].min() >= -1e-8)
        weight = integrate(data[:, 1] / 256, data[:, 0])
        bulk_weight = (2 / np.pi) * np.arcsin(0.95)
        c.near(f"eta={eta} N=256 DOS weight on stored interval", weight, bulk_weight,
               0.025, "finite-chain level counting and Jackson broadening at interval ends")
        if eta == 80:
            mask = np.abs(data[:, 0]) < 1.9
            ref = 1 / (np.pi * np.sqrt(4 - data[mask, 0] ** 2))
            c.near("N=256 eta=80 bulk DOS error (1/eV)",
                   np.max(np.abs(data[mask, 1] / 256 - ref)), 0, 5e-5,
                   "below the finite-size revival; printed energy/DOS resolution")


if __name__ == "__main__":
    raise SystemExit(run(check, "run_dos"))
