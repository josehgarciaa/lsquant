#!/usr/bin/env python3
"""Plot sigma^z_xy of pentagonal PdSe2 from every route on one energy axis.

Usage (from the repository root):
    python examples/20_pdse2_spin_hall_routes/plot.py [data_dir]

Reads the curve files in data_dir (default: data/ next to this script) and
writes fig_routes.png next to this script. Energies are measured from the
valence-band maximum (VBM) of each electronic structure.
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

# VBM of each electronic structure, E - E_F in eV. In an insulator each code
# places its Fermi level at a different point of the gap (smearing, electronic
# temperature), so the VBM is the common origin. Computed from the bands of
# each model on a 60 x 60 mesh refined around the maximum.
VBM = {"siesta": -0.972, "qe": -0.822, "paoflow": -0.827, "openmx": -0.606}
# Topological gap of each electronic structure, E - E_VBM in eV.
GAP = {"siesta": (1.986, 2.061), "qe": (2.009, 2.055), "openmx": (1.918, 1.982)}

SMOOTH = [  # file, origin, label, trace index, line style
    ("siesta_projected_kpm.dat", "siesta", "SIESTA, KPM", 0, "-"),
    ("wannier_kpm.dat", "qe", "QE + Wannier90, KPM", 3, "--"),
    ("paoflow_kpm.dat", "paoflow", "QE + PAOFLOW, KPM", 4, (0, (1, 1.5))),
    ("openmx_exact.dat", "openmx", "OpenMX, exact", 1, "-."),
]
# paoflow.dat is PAOFLOW's spin_Hall output times 4 pi: PAOFLOW gives e^2/hbar
# with the spin operator sigma/2 (PAOFLOW.py: sP = 0.5 Pauli; do_Hall.py: sum
# over k divided by N_k and by the cell area), and 2 pi * 2 converts it to e^2/h
# with the Pauli sigma_z.
EXACT = [  # exact curves, shown in the gap panel only: thin line and a marker each
    ("siesta_direct_exact.dat", "siesta", "SIESTA, exact", 0, "o"),
    ("wannier_covariant_spn.dat", "qe", "QE + Wannier90, exact", 3, "s"),
    ("paoflow.dat", "paoflow", "PAOFLOW, own result", 4, "^"),
]


def main(data=os.path.join(_HERE, "data")):
    load = lambda f: np.loadtxt(os.path.join(data, f)).T
    prl.use("web", height=4.2)
    fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw=dict(width_ratios=[1.6, 1]))
    for ax in (ax1, ax2):
        for key in ("siesta", "qe", "openmx"):   # overlaps show darker
            ax.axvspan(*GAP[key], color="0.6", lw=0, alpha=0.18)
        ax.axhline(0.0, color="0.6", lw=0.5)
    for f, origin, label, i, ls in SMOOTH:
        E, s = load(f)
        style = dict(prl.trace(i), linestyle=ls)
        ax1.plot(E - VBM[origin], s, label=label, **style)
        ax2.plot(E - VBM[origin], s, **style)
    for f, origin, label, i, marker in EXACT:
        E, s = load(f)
        ax2.plot(E - VBM[origin], s, label=label, **dict(prl.trace(i), linestyle="-", linewidth=0.8,
                                                           marker=marker, markevery=4, markersize=3.5,
                                                           markerfacecolor="none"))
    ax1.set_xlim(0.0, 2.6)
    ax1.set_ylim(-1.5, 2.6)
    ax2.set_xlim(1.88, 2.12)
    ax2.set_ylim(1.4, 2.15)
    ax1.set_ylabel(r"$\sigma^z_{xy}$ ($e^2/h$)")
    for ax in (ax1, ax2):
        ax.set_xlabel(r"$E - E_\mathrm{VBM}$ (eV)")
    ax1.legend(loc="upper left", fontsize="small")
    ax2.legend(loc="lower center", fontsize="x-small")
    prl.panel(ax1, "(a)", loc="lower right")
    prl.panel(ax2, "(b)", loc="upper left")
    prl.save(fig, os.path.join(_HERE, "fig_routes"))


if __name__ == "__main__":
    main(*sys.argv[1:2])
