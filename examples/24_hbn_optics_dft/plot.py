#!/usr/bin/env python3
"""Figure of example 24: the optical conductivity of monolayer hBN from seven DFT models.

Reads curves.npz and analysis.json (default: data/ next to this script) and
writes fig_hbn_dft.png.
(a) Re sigma, PBE: SIESTA, PySCF, Quantum ESPRESSO + Wannier90, QE + PAOFLOW.
(b) Re sigma, HSE06: PySCF, QE + Wannier90, QE + PAOFLOW.
(c) Im sigma, PBE, with the Kramers-Kronig transform of Re sigma plus the fitted
    smooth tail of the truncated transform (circles, up to 12 eV).
(d) The 2D polarisability below the gap, every model.
Lines: Chebyshev moments. Markers in (a) and (b): exact sum over states of the same model.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
COLOR = {"siesta": "#000000", "pyscf": "#0072B2", "w90": "#D55E00", "pao": "#009E73"}
CODE = {"siesta": "SIESTA", "pyscf": "PySCF", "w90": "QE + Wannier90", "pao": "QE + PAOFLOW"}
STYLE = {"pbe": "-", "hse": "--"}


def main(data=os.path.join(_HERE, "data")):
    c = np.load(os.path.join(data, "curves.npz"))
    rep = json.load(open(os.path.join(data, "analysis.json")))
    w = c["photon"]
    plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 0.8,
                         "xtick.direction": "in", "ytick.direction": "in",
                         "xtick.top": True, "ytick.right": True})
    fig, ((a, b), (cc, d)) = plt.subplots(2, 2, figsize=(8.4, 6.4))
    marker = dict(ls="none", marker="s", ms=3.2, mfc="none", mew=0.8)
    for ax, functional, title in ((a, "pbe", "PBE"), (b, "hse", "HSE06")):
        for code in COLOR:
            name = f"{code}_{functional}"
            if name not in rep:
                continue
            ax.plot(w, c[f"{name}__re"], color=COLOR[code], lw=1.1)
            ax.plot(w[::25], c[f"{name}__sos_re"][::25], mec=COLOR[code], **marker)
        ax.set_xlim(3, 15)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r"$\hbar\omega$ (eV)")
        ax.set_ylabel(r"Re $\sigma\,(e^2/h)$")
        ax.text(0.03, 0.95, title, transform=ax.transAxes, va="top")
    for code in COLOR:
        name = f"{code}_pbe"
        if name not in rep:
            continue
        cc.plot(w, c[f"{name}__im"], color=COLOR[code], lw=1.1)
        kk = w <= 12.0                                  # the truncated transform fails near 15 eV
        cc.plot(w[kk][::25], c[f"{name}__kk"][kk][::25], ls="none", marker="o", ms=3.4, mfc="none",
                mew=0.8, mec=COLOR[code])
    cc.axhline(0, color="0.8", lw=0.5)
    cc.set_xlim(0, 15)
    cc.set_xlabel(r"$\hbar\omega$ (eV)")
    cc.set_ylabel(r"Im $\sigma\,(e^2/h)$")
    cc.text(0.03, 0.95, "PBE; circles: Kramers–Kronig of Re $\\sigma$", transform=cc.transAxes, va="top",
            fontsize=8.5)
    for name, r in rep.items():
        code, functional = name.split("_")
        below = w < r["direct_gap"] - 0.3
        d.plot(w[below], c[f"{name}__alpha"][below], color=COLOR[code], ls=STYLE[functional], lw=1.1)
    d.set_xlim(0, None)
    d.set_xlabel(r"$\hbar\omega$ (eV)")
    d.set_ylabel(r"$\alpha_{\rm 2D} = -e^2\,{\rm Im}\,\sigma/\hbar\omega$ (Å)")
    d.text(0.03, 0.95, "solid: PBE    dashed: HSE06", transform=d.transAxes, va="top", fontsize=8.5)
    handles = [plt.Line2D([], [], color=COLOR[code], lw=1.4, label=CODE[code]) for code in COLOR
               if any(n.startswith(code + "_") for n in rep)]
    handles.append(plt.Line2D([], [], ls="none", marker="s", mfc="none", mec="0.3", ms=4,
                              label="exact sum over states"))
    fig.legend(handles=handles, loc="upper center", ncol=len(handles), fontsize=8.5, frameon=False,
               handlelength=1.8, bbox_to_anchor=(0.5, 1.0))
    for ax, letter in zip((a, b, cc, d), "abcd"):
        ax.text(0.97, 0.05 if ax is not d else 0.05, f"({letter})", transform=ax.transAxes,
                ha="right", va="bottom", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(_HERE, "fig_hbn_dft.png"), dpi=200)


if __name__ == "__main__":
    main(*sys.argv[1:])
