#!/usr/bin/env python3
"""Plot the bands and the anomalous Hall conductivity of the YN2 monolayer.

Usage (from the repository root):
    python examples/22_yn2_anomalous_hall/plot.py [data_dir]

Reads the files in data_dir (default: data/ next to this script) and writes
fig_yn2.png next to this script. Energies are measured from the Fermi level
of each electronic structure.
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


def main(data=os.path.join(_HERE, "data")):
    load = lambda f: np.loadtxt(os.path.join(data, f)).T
    prl.use("web", height=3.3)
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, gridspec_kw=dict(width_ratios=[1.1, 1, 1]))

    # (a) Bands: SIESTA (numerical atomic orbitals) and QE + Wannier90 (plane waves).
    b = np.load(os.path.join(data, "bands.npz"))
    s, w = prl.trace(0), prl.trace(3)
    ax1.plot(b["x_siesta"], b["E_siesta"], color=s["color"], lw=1.0)
    ax1.plot(b["x_wannier"], b["E_wannier"], color=w["color"], lw=1.0, ls="--")
    x = b["x_siesta"]
    ticks = [x[0], x[120], x[240], x[-1]]
    for t in ticks:
        ax1.axvline(t, color="0.8", lw=0.5)
    ax1.axhline(0.0, color="0.6", lw=0.5, ls=":")
    ax1.set_xticks(ticks, [r"$\Gamma$", "M", "K", r"$\Gamma$"])
    ax1.set_xlim(x[0], x[-1])
    ax1.set_ylim(-1.2, 1.2)
    ax1.set_ylabel(r"$E-E_F$ (eV)")
    ax1.plot([], [], color=s["color"], lw=1.0, label="SIESTA")
    ax1.plot([], [], color=w["color"], lw=1.0, ls="--", label="QE + Wannier90")
    ax1.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, fontsize=7)

    # (b) The C = 3 plateau: exact sums for both models and the KPM curve.
    for f, i, ls, label in (("siesta_exact.dat", 0, "-", "SIESTA, exact"),
                            ("wannier_exact.dat", 3, "--", "QE + Wannier90, exact"),
                            ("wannier_kpm_ahc.dat", 1, ":", "QE + Wannier90, KPM")):
        E, sig = load(f)
        ax2.plot(E, sig, color=prl.trace(i)["color"], ls=ls, lw=1.2, label=label)
    ax2.axhline(-3.0, color="0.6", lw=0.5)
    ax2.set_xlim(-0.1, 0.1)
    ax2.set_ylim(-3.4, 0.1)
    ax2.set_xlabel(r"$E-E_F$ (eV)")
    ax2.set_ylabel(r"$\sigma_{xy}\,(e^2/h)$")
    ax2.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=6.5)

    # (c) KPM: the Fermi-sea and Fermi-surface parts of the Bastin curve.
    for f, i, ls, label in (("wannier_kpm_ahc.dat", 1, ":", "total"),
                            ("wannier_kpm_sea.dat", 4, "-", "Fermi sea"),
                            ("wannier_kpm_surf.dat", 5, "-.", "Fermi surface")):
        E, sig = load(f)
        ax3.plot(E, sig, color=prl.trace(i)["color"], ls=ls, lw=1.2, label=label)
    ax3.axhline(0.0, color="0.6", lw=0.5)
    ax3.set_xlim(-0.1, 0.1)
    ax3.set_ylim(-3.4, 0.1)
    ax3.set_xlabel(r"$E-E_F$ (eV)")
    ax3.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=6.5)

    for ax, tag in zip((ax1, ax2, ax3), "abc"):
        prl.panel(ax, f"({tag})", loc="lower right",
                  bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
    fig.tight_layout()
    fig.savefig(os.path.join(_HERE, "fig_yn2.png"), dpi=200)


if __name__ == "__main__":
    main(*sys.argv[1:])
