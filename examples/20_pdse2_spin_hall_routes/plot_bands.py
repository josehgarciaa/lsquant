#!/usr/bin/env python3
"""Plot the PdSe2 bands of Quantum ESPRESSO and OpenMX from the valence-band maximum.

Usage (from the repository root):
    python examples/20_pdse2_spin_hall_routes/plot_bands.py [data_dir]

Reads data/bands_qe_openmx.npz and writes fig_bands.png next to this script.
The file holds, on the path G-X-S-Y-G-S (120 points per segment), the bands
of the Quantum ESPRESSO + Wannier90 model, the Quantum ESPRESSO eigenvalues at
the points of the 12 x 12 nscf mesh that lie on the path, and the OpenMX bands
(H and S of the .scfout read by LSQuant), each minus its own VBM.
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
from plot import GAP


def main(data=os.path.join(_HERE, "data")):
    d = np.load(os.path.join(data, "bands_qe_openmx.npz"))
    prl.use("single", height=3.0)
    fig, ax = plt.subplots()
    qe, omx = prl.trace(3)["color"], prl.trace(1)["color"]
    ax.plot(d["x"], d["openmx"], color=omx, lw=0.9)
    ax.plot(d["x"], d["wannier"], color=qe, lw=0.9, ls="--")
    ax.plot(d["qe_x"], d["qe"], "o", ms=2.2, mfc="none", mec=qe, mew=0.6)
    ax.plot([], [], color=qe, ls="--", marker="o", ms=2.2, mfc="none", label="QE: points, Wannier90: line")
    ax.plot([], [], color=omx, label="OpenMX")
    for key in ("qe", "openmx"):
        ax.axhspan(*GAP[key], color="0.6", lw=0, alpha=0.25)
    for t in d["ticks"]:
        ax.axvline(t, color="0.85", lw=0.5)
    ax.set_xticks(d["ticks"], [r"$\Gamma$", "X", "S", "Y", r"$\Gamma$", "S"])
    ax.set_xlim(d["ticks"][0], d["ticks"][-1])
    ax.set_ylim(-1.5, 2.6)
    ax.set_ylabel(r"$E-E_\mathrm{VBM}$ (eV)")
    ax.legend(loc="center left", bbox_to_anchor=(0.0, 0.55), fontsize=6.5)
    fig.tight_layout()
    fig.savefig(os.path.join(_HERE, "fig_bands.png"), dpi=300)


if __name__ == "__main__":
    main(*sys.argv[1:])
