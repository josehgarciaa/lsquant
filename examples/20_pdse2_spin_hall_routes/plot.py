#!/usr/bin/env python3
"""Plot sigma^z_xy(E) of pentagonal PdSe2 from four routes on one energy axis.

Usage (from the repository root):
    python examples/20_pdse2_spin_hall_routes/plot.py [data_dir]

Reads the curve files in data_dir (default: data/ next to this script), writes
fig_routes.png next to this script and prints the plateau of each route in
its own quantum spin Hall gap and the value at E_F.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
import prl_style as prl

# Quantum spin Hall gap of each electronic structure, E - E_F in eV (global gap on
# the k mesh of the exact routes: 24 x 24 for SIESTA, 300 x 300 for Wannier90).
GAP_SIESTA = (1.014, 1.090)
GAP_QE = (1.188, 1.233)

ROUTES = [  # file, label, gap used for the plateau, trace index, line style
    ("siesta_direct_exact.dat", "SIESTA, exact", GAP_SIESTA, 0, "-"),
    ("siesta_projected_exact.dat", "SIESTA, $[-3,3]$ eV, exact", GAP_SIESTA, 1, "--"),
    ("siesta_projected_kpm.dat", "SIESTA, $[-3,3]$ eV, KPM", GAP_SIESTA, 2, ":"),
    ("wannier_covariant_spn.dat", "QE + Wannier90, LSQuant", GAP_QE, 3, "-."),
    ("paoflow.dat", "QE + PAOFLOW", GAP_QE, 4, "--"),
]


def main(data=os.path.join(_HERE, "data")):
    curves = {f: np.loadtxt(os.path.join(data, f)) for f, *_ in ROUTES}
    print(f"{'route':32s} {'gap (eV)':>13s} {'plateau':>8s} {'spread':>7s} {'at E_F':>7s}")
    for f, label, gap, _, _ in ROUTES:
        E, s = curves[f].T
        # The central half of the gap avoids the broadened band edges.
        lo, hi = gap[0] + 0.25 * (gap[1] - gap[0]), gap[1] - 0.25 * (gap[1] - gap[0])
        m = (E >= lo) & (E <= hi)
        print(f"{label:32s} {gap[0]:.3f}-{gap[1]:.3f} {s[m].mean():8.4f} {np.ptp(s[m]):7.4f} "
              f"{np.interp(0.0, E, s):7.4f}")

    prl.use("web", height=6.0)
    fig, (ax1, ax2) = plt.subplots(2, 1)
    for ax, xlim in ((ax1, (-0.5, 2.0)), (ax2, (0.95, 1.32))):
        ax.axvspan(*GAP_SIESTA, color="0.88", lw=0)
        ax.axvspan(*GAP_QE, color="0.94", lw=0)
        for f, label, _, i, ls in ROUTES:
            E, s = curves[f].T
            style = prl.trace(i)
            style["linestyle"] = ls
            ax.plot(E, s, label=label, **style)
        ax.axhline(0.0, color="0.6", lw=0.5)
        ax.set_xlim(*xlim)
        ax.set_ylabel(r"$\sigma^z_{xy}$ ($e^2/h$)")
    ax2.set_ylim(1.2, 2.3)
    ax1.legend(loc="upper left", fontsize="small")
    ax2.set_xlabel(r"$E - E_F$ (eV)")
    prl.panel(ax1, "(a)", loc="lower right")
    prl.panel(ax2, "(b)", loc="lower left")
    prl.save(fig, os.path.join(_HERE, "fig_routes"))


if __name__ == "__main__":
    main(*sys.argv[1:2])
