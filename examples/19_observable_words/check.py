#!/usr/bin/env python3
"""Check the five saved Kane-Mele observables against the exact-mesh oracle.

The oracle is generated separately by oracle.py, not by this check. Charge
conductivities need 2*pi/(N_cells*A_cell); spin Hall files already use e^2/h.
DOS and spin DOS must be divided by N_cells, leaving four states per cell.
The comparison uses the same interpolation and |E|<3.5 eV window as plot.py.
Work is linear in five curves and the oracle JSON arrays, without any KPM run.
Pitfall: Jackson and Gaussian peaks differ at 25 versus 20 meV resolution;
their documented differences are error budgets, not round-off errors.
Warning: the optional 25 meV oracle tests the quoted conductivity bracket;
the large-cell random-trace illustration is a different run and is not tested.
"""
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.dont_write_bytecode = True
from check_support import one, run, table
import numpy as np


def configure(parser):
    """Allow independent reference files to live alongside any temporary run."""
    parser.add_argument("--oracle", type=Path, help="20 meV oracle JSON")
    parser.add_argument("--oracle-25", type=Path, help="optional 25 meV oracle JSON")


def check(root, c, args):
    """Test oracle landmarks, five full curves, gap values and spectral weights."""
    from plot import curve
    oracle_path = args.oracle or root.parent / "operators" / "km_L16.oracle.json"
    with oracle_path.open() as stream:
        oracle = json.load(stream)
    with one(root, "run.resolved.json").open() as stream:
        resolved = json.load(stream)
    rec = resolved["reconstruction"]
    cells, area = rec["num_cells"], rec["cell_area"]
    c.require("published cell and mesh", cells == 256 and oracle["L"] == 16)
    c.near("cell area (Angstrom^2)", area, 3 * np.sqrt(3) / 2, 1e-12, "honeycomb geometry")
    c.near("oracle broadening (meV)", oracle["eta_meV"], 20, 0, "published reference")
    for key, target in (("gap_eV", 0.922), ("vbm_eV", -0.510), ("cbm_eV", 0.412)):
        c.near(key, oracle[key], target, 0.001, "README band edges rounded to three decimals")
    c.near("spin Chern", oracle["spin_chern"]["spin_chern"], -2, 1e-8, "integer spin-sector topology")
    c.near("Z2", oracle["spin_chern"]["z2"], 1, 0, "spin Chern parity")
    energies = np.asarray(oracle["energies_eV"], dtype=float)
    if not np.isfinite(energies).all() or not np.all(np.diff(energies) > 0):
        raise ValueError("invalid oracle energy grid")
    window = (energies > -3.5) & (energies < 3.5)
    c.require("oracle covers comparison window", energies[0] <= -3.5 and energies[-1] >= 3.5)
    to_e2h = 2 * np.pi / (cells * area)
    specifications = (
        ("dos", "meanSpectralOp1", "dos_per_cell", 1 / cells, 0, 0.45, 0.002),
        ("spin", "meanSpectralOpSZ", "spin_dos_per_cell", 1 / cells, 0, 0.40, 0.002),
        ("sxx", "KuboGreenWood_NonEqOpVX-VX", "sigma_xx", to_e2h, 0, 110, 0.002),
        ("sxy", "KuboBastin_NonEqOpVX-VY", "sigma_xy", to_e2h, 0.0001, 0.20, 0.002),
        ("shc", "SpinHall_NonEqOpJXSZ-VY", "sigma_z_xy", 1, -2.021, 0.06, 0.015),
    )
    integrals = {}
    # np.trapz supports older NumPy installations as well as current ones.
    integrate = getattr(np, "trapezoid", None)
    if integrate is None:
        integrate = np.trapz
    for name, prefix, ref_key, factor, at_zero, max_error, zero_tol in specifications:
        table(one(root, prefix + "*JACKSON.dat"))
        x, values = curve(str(root), prefix + "*JACKSON.dat")
        y = values * factor
        ref = np.asarray(oracle[ref_key], dtype=float)
        if ref.shape != energies.shape or not np.isfinite(ref).all():
            raise ValueError(f"invalid oracle {ref_key}")
        c.require(f"{name} curve covers comparison window", x[0] <= -3.5 and x[-1] >= 3.5)
        sampled = np.interp(energies[window], x, y)
        delta = np.max(np.abs(sampled - ref[window]))
        c.near(f"{name} max |KPM-oracle|", delta, 0, max_error,
               "documented peak differences plus margin for Jackson/Gaussian resolution")
        c.near(f"{name} KPM at EF=0", np.interp(0, x, y), at_zero, zero_tol,
               "gap leakage at finite Jackson resolution")
        ref_zero = -2.018 if name == "shc" else (0.0001 if name == "sxy" else 0)
        c.near(f"{name} oracle at EF=0", np.interp(0, energies, ref), ref_zero,
               0.001 if name == "shc" else 0.0002, "rounded exact-mesh reference")
        integrals[name] = (float(integrate(sampled, energies[window])),
                           float(integrate(ref[window], energies[window])))
        if name in ("dos", "spin", "sxx"):
            for origin, value, target, tolerance in (
                    ("KPM", integrals[name][0], {"dos": 4, "spin": 0, "sxx": 297}[name],
                     {"dos": 0.02, "spin": 0.02, "sxx": 12}[name]),
                    ("oracle", integrals[name][1], {"dos": 4, "spin": 0, "sxx": 313}[name],
                     {"dos": 0.02, "spin": 0.02, "sxx": 3}[name])):
                c.near(f"{name} {origin} integral over |E|<3.5", value, target, tolerance,
                       "finite integration grid, tails and rounded published integral")
        if name in ("dos", "spin", "sxx"):
            peak = {"dos": 4.0, "spin": 1.8, "sxx": 320}[name]
            c.near(f"{name} oracle peak magnitude", np.max(np.abs(ref[window])), peak,
                   0.20 * peak, "approximate peak heights and finite mesh resolution")
        if name == "sxy":
            near_one = (x > 0.7) & (x < 1.3)
            c.near("Hall minimum near +1 eV (e^2/h)", y[near_one].min(), -0.8, 0.25,
                   "approximate band Hall response, finite energy resolution")
    if args.oracle_25:
        with args.oracle_25.open() as stream:
            reference25 = json.load(stream)
        c.require("25 meV reference parameters", reference25["L"] == 16 and reference25["eta_meV"] == 25)
        e25 = np.asarray(reference25["energies_eV"])
        r25 = np.asarray(reference25["sigma_xx"])
        if r25.shape != e25.shape or not np.isfinite(r25).all():
            raise ValueError("invalid 25 meV conductivity")
        mask25 = (e25 > -3.5) & (e25 < 3.5)
        integral25 = float(integrate(r25[mask25], e25[mask25]))
        c.near("25 meV oracle conductivity integral", integral25, 251, 3,
               "rounded Gaussian-broadening reference")
        c.require("KPM conductivity integral bracketed by Gaussian references",
                  integral25 < integrals["sxx"][0] < integrals["sxx"][1])
    else:
        print("25 meV conductivity bracket not checked; supply --oracle-25 for that additional reference.")


if __name__ == "__main__":
    raise SystemExit(run(check, configure=configure))
