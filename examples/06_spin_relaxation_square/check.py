#!/usr/bin/env python3
"""Check the observed spin relaxation, including its non-golden-rule rate law.

The example explicitly reports p=1.26, not the ideal p=2. Fit the saved rates
and the spin curves with the same normalization and fit interval as plot.py.
Rates are fs^-1, lifetimes fs, and spin-flip disorder eV. The budgets below
allow stochastic trace noise, finite size and the broad 140 meV energy filter.
Cost is linear in a few short tables; no spin evolution is performed.
Warning: agreement with this example does not establish the golden-rule regime.
To-do: separate filter, observation-time and finite-size effects on the offset.
"""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.dont_write_bytecode = True
from check_support import one, run, same_grid, table
import numpy as np


def check(root, c, args):
    """Compare observed rates and fits, and test the documented spin anisotropy."""
    from plot import _load_normalised_spin, fit_rate
    rates = table(one(root, "ey_ratelaw.dat"))
    ds, rate = rates.T
    expected_ds = np.array([0.04, 0.06, 0.08, 0.10, 0.14])
    c.require("rate-law disorder grid", ds.shape == expected_ds.shape and
              np.allclose(ds, expected_ds, rtol=0, atol=1e-12))
    if ds.shape != expected_ds.shape:
        return
    c.require("positive rates", np.all(rate > 0))
    for d, value, target in zip(ds, rate, (0.00832, 0.01387, 0.01897, 0.02417, 0.04187)):
        c.near(f"Dsf={d:g} rate (1/fs)", value, target, 0.15 * target,
               "20 random vectors and rounded published rates")
    p = np.polyfit(np.log(ds), np.log(rate), 1)[0]
    b, a = np.polyfit(ds ** 2, rate, 1)
    residual = np.max(np.abs((a + b * ds ** 2 - rate) / rate))
    c.near("log-log exponent", p, 1.26, 0.15, "finite filter/time and stochastic rate fit")
    c.near("quadratic coefficient (1/fs/eV^2)", b, 1.81, 0.3, "observed finite-size fit")
    c.near("rate offset (1/fs)", a, 0.0065, 0.002, "unexplained finite-filter/time offset")
    c.near("max relative quadratic residual", residual, 0.13, 0.08, "rounded observed residual")
    bgr = 4 * np.pi / 0.6582 * 0.1093
    c.near("golden-rule coefficient (1/fs/eV^2)", bgr, 2.09, 0.01, "rounded analytic clean DOS")
    c.near("fitted / golden-rule coefficient", b / bgr, 0.87, 0.15, "finite filter/time and trace noise")
    c.near("weakest-disorder lifetime (fs)", 1 / rate[0], 120, 20, "rounded stochastic rate")

    def spin(label, op, proj):
        curves = [table(one(root, f"meanTimeProj{x}-{proj}{label}*EF2.000000JACKSON.dat"))
                  for x in (op, "1")]
        same_grid(*curves)
        c.require(f"{label} {op} nonzero norm", abs(curves[1][1:, 1].mean()) > 1e-12)
        return _load_normalised_spin(str(root), label, op, proj)

    lifetimes = []
    for d in (0.05, 0.07, 0.1):
        t, s = spin(f"ey_L20_D{d:g}", "SZ", "PZ")
        lifetime = fit_rate(t, s)
        c.require(f"Dsf={d:g} positive EY lifetime (fs)", np.isfinite(lifetime) and lifetime > 0)
        print(f"Dsf={d:g} EY lifetime = {lifetime:.6g} fs")
        lifetimes.append(lifetime)
    c.require("EY lifetime decreases with spin-flip strength", np.all(np.diff(lifetimes) < 0))
    t0, sz0 = spin("rashba_L20_l0.2_W0", "SZ", "PZ")
    tz, sz = spin("rashba_L20_l0.2_W6", "SZ", "PZ")
    tx, sx = spin("rashba_L20_l0.2_W6", "SX", "PX")
    c.require("clean Rashba precession changes sign", sz0.min() < 0)
    # A finite cell need not attain the asymptotic factor of two; check its direction.
    tz_fit, tx_fit = fit_rate(tz, sz), fit_rate(tx, sx)
    print(f"DP fitted lifetimes: Tz={tz_fit:.6g} fs, Tx={tx_fit:.6g} fs (ideal Tx/Tz=2)")
    c.require("DP in-plane spin outlives out-of-plane spin", 0 < tz_fit < tx_fit < np.inf)


if __name__ == "__main__":
    raise SystemExit(run(check, "run_relaxation"))
