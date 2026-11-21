#!/usr/bin/env python3
"""Plot the five observables of example 19 against the exact reference.

Usage (from the repository root, after the run, the shc reconstruction and oracle.py):
    python examples/19_observable_words/plot.py [out_dir] [label]

Reads the curve files of `python -m lsquant run` / `reconstruct` in out_dir
(default examples/19_observable_words/out) and operators/<label>.oracle.json
(default km_L16), writes fig_words.png next to this script and prints the checks.
"""
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
import prl_style as prl


def curve(out, pattern):
    paths = glob.glob(os.path.join(out, pattern))
    if len(paths) != 1:
        raise SystemExit(f"expected one file {pattern} in {out}, found {len(paths)}")
    data = np.loadtxt(paths[0])
    return data[:, 0], data[:, 1]


def main(out=os.path.join(_HERE, "out"), label="km_L16"):
    with open(os.path.join(_HERE, "operators", f"{label}.json")) as f:
        meta = json.load(f)
    with open(os.path.join(_HERE, "operators", f"{label}.oracle.json")) as f:
        oracle = json.load(f)
    n_cells, area = meta["num_cells"], meta["cell_area_A2"]
    to_e2h = 2 * np.pi / (n_cells * area)          # Kubo curves -> e^2/h, as the shc formula
    kpm = {
        "dos": (curve(out, "meanSpectralOp1*JACKSON.dat"), 1 / n_cells),
        "spin": (curve(out, "meanSpectralOpSZ*JACKSON.dat"), 1 / n_cells),
        "sxx": (curve(out, "KuboGreenWood_NonEqOpVX-VX*JACKSON.dat"), to_e2h),
        "sxy": (curve(out, "KuboBastin_NonEqOpVX-VY*JACKSON.dat"), to_e2h),
        "shc": (curve(out, "SpinHall_NonEqOpJXSZ-VY*JACKSON.dat"), 1.0),   # already e^2/h
    }
    ref = {"dos": "dos_per_cell", "spin": "spin_dos_per_cell", "sxx": "sigma_xx",
           "sxy": "sigma_xy", "shc": "sigma_z_xy"}
    E = np.array(oracle["energies_eV"])
    ylabels = {
        "dos": r"$\rho(E)$ (eV$^{-1}$)",
        "spin": r"$\rho_z(E)$ (eV$^{-1}$)",
        "sxx": r"$\sigma_{xx}$ ($e^2/h$)",
        "sxy": r"$\sigma_{xy}$ ($e^2/h$)",
        "shc": r"$\sigma^z_{xy}$ ($e^2/h$)",
    }

    # Checks printed for the README.
    window = (E > -3.5) & (E < 3.5)
    print(f"gap {oracle['gap_eV']:.3f} eV from {oracle['vbm_eV']:.3f} to {oracle['cbm_eV']:.3f} eV; "
          f"spin Chern {oracle['spin_chern']['spin_chern']:.3f}, Z2 {oracle['spin_chern']['z2']}")
    for key, ((x, y), scale) in kpm.items():
        k = np.interp(E[window], x, y * scale)
        r = np.array(oracle[ref[key]])[window]
        line = (f"{key:4s}: KPM(0) = {np.interp(0.0, x, y * scale):+.4f}, oracle(0) = "
                f"{np.interp(0.0, E, oracle[ref[key]]):+.4f}, max |KPM - oracle| = {np.max(np.abs(k - r)):.3g}")
        if key in ("dos", "spin", "sxx"):
            line += f", integral KPM {np.trapezoid(k, E[window]):.4f} oracle {np.trapezoid(r, E[window]):.4f}"
        print(line)

    prl.use("web", height=9.0)
    fig, axes = plt.subplots(5, 1, sharex=True)
    for i, (ax, key) in enumerate(zip(axes, kpm)):
        (x, y), scale = kpm[key]
        ax.axvspan(oracle["vbm_eV"], oracle["cbm_eV"], color="0.9", lw=0)
        ax.plot(E, oracle[ref[key]], label="exact, $k$ mesh", **prl.trace(1))
        ax.plot(x, y * scale, label="KPM", **prl.trace(0))
        ax.set_ylabel(ylabels[key])
        ax.axhline(0.0, color="0.6", lw=0.5)
        prl.panel(ax, f"({'abcde'[i]})")
    axes[0].legend(loc="upper right")
    axes[-1].set_xlabel(r"$E_F$ (eV)")
    axes[-1].set_xlim(-3.5, 3.5)
    prl.save(fig, os.path.join(_HERE, "fig_words"))


if __name__ == "__main__":
    main(*sys.argv[1:3])
