#!/usr/bin/env python3
"""Plot sigma^z_xy of PdSe2 from five models of different size.

Usage (from the repository root):
    python examples/21_pdse2_downfolding/plot.py

Reads data/ next to this script and the curves of example 20, writes
fig_downfolding.png. Energies are measured from the valence-band maximum.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_EX20 = os.path.join(os.path.dirname(_HERE), "20_pdse2_spin_hall_routes", "data")
sys.path.insert(0, os.path.dirname(_HERE))
import prl_style as prl

# VBM, E - E_F in eV (example 20). The reduced SIESTA model keeps the SIESTA
# bands of its window unchanged, so it has the VBM of SIESTA.
VBM = {"siesta": -0.972, "qe": -0.822, "paoflow": -0.827}
GAP = {"siesta": (1.986, 2.061), "qe": (2.009, 2.055)}

CURVES = [  # path, origin, label, trace index, line style
    (os.path.join(_EX20, "siesta_direct_exact.dat"), "siesta", "SIESTA, all 220 orbitals (exact)", 0, "-"),
    (os.path.join(_EX20, "siesta_projected_kpm.dat"), "siesta", "SIESTA, window $\\pm$3 eV (KPM)", 0, "--"),
    (os.path.join(_HERE, "data", "siesta_reduced12_kpm.dat"), "siesta", "SIESTA, reduced to 12 bands (KPM)", 1, "-"),
    (os.path.join(_EX20, "wannier_kpm.dat"), "qe", "QE + Wannier90, 12 orbitals (KPM)", 3, "-."),
    (os.path.join(_EX20, "paoflow_kpm.dat"), "paoflow", "QE + PAOFLOW, 108 orbitals (KPM)", 4, (0, (1, 1.5))),
]


def main():
    prl.use("web", height=4.2)
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.6, 1]))
    for ax in (ax1, ax2):
        for key in GAP:
            ax.axvspan(*GAP[key], color="0.6", lw=0, alpha=0.18)
        ax.axhline(0.0, color="0.6", lw=0.5)
    for path, origin, label, i, ls in CURVES:
        E, s = np.loadtxt(path).T
        style = dict(prl.trace(i), linestyle=ls)
        if "reduced12" in path:
            # Thin line with circles: apart from the full SIESTA curve it overlaps.
            style.update(linewidth=0.9, marker="o", markevery=12, markersize=3.5, markerfacecolor="none")
        ax1.plot(E - VBM[origin], s, label=label, **style)
        ax2.plot(E - VBM[origin], s, **style)
    ax1.set_xlim(0.0, 2.6)
    ax1.set_ylim(-1.5, 4.2)   # room for the legend above the curves
    ax2.set_xlim(1.94, 2.10)
    ax2.set_ylim(1.7, 2.1)
    ax1.set_ylabel(r"$\sigma^z_{xy}$ ($e^2/h$)")
    for ax in (ax1, ax2):
        ax.set_xlabel(r"$E - E_\mathrm{VBM}$ (eV)")
    ax1.legend(loc="upper left", fontsize="small", frameon=True, framealpha=0.9, edgecolor="0.8")
    prl.panel(ax1, "(a)", loc="lower right")
    prl.panel(ax2, "(b)", loc="lower left")
    prl.save(fig, os.path.join(_HERE, "fig_downfolding"))


if __name__ == "__main__":
    main()
