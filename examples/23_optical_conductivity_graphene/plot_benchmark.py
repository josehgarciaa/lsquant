#!/usr/bin/env python3
"""Plot the optical tests of graphene and hBN written by benchmark.py.

Usage (from this directory):
    python plot_benchmark.py [data_dir]
Reads bench_graphene.npz, bench_hbn.npz and bench_map.npz (default: data/ next
to this script) and writes fig_benchmark.png next to this script.
Top row graphene, bottom row hBN: Re sigma, Im sigma with the Kramers-Kronig
transform of Re sigma, and the (mu, omega) map (graphene) or the 2D
polarisability below the gap (hBN).
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))
sys.path.insert(0, _HERE)
import prl_style as prl                      # noqa: E402
from make_honeycomb import GRAPHENE, HBN     # noqa: E402


def main(data=os.path.join(_HERE, "data")):
    g = np.load(os.path.join(data, "bench_graphene.npz"))
    h = np.load(os.path.join(data, "bench_hbn.npz"))
    mp = np.load(os.path.join(data, "bench_map.npz"))
    prl.use("web", height=6.0)
    fig = plt.figure(figsize=(10.2, 6.0))
    gs = fig.add_gridspec(2, 3)
    a, b, c = (fig.add_subplot(gs[0, i]) for i in range(3))
    d, e, f = (fig.add_subplot(gs[1, i]) for i in range(3))
    k0, k1, k2, k3 = (prl.trace(i)["color"] for i in range(4))
    marker = dict(ls="none", marker="s", ms=3.2, mfc="none", mew=0.8)

    def real_part(ax, m, model, xlim, every, start):
        w = m["photon"][2:]
        ax.plot(w, m["re_xx"][2:], color=k0, lw=1.2, label=r"moments, $\sigma_{xx}$")
        ax.plot(w, m["re_yy"][2:], color=k1, lw=1.0, ls="--", label=r"moments, $\sigma_{yy}$")
        ax.plot(w[start::every], m["exact_re"][2:][start::every], mec=k2, label="exact", **marker)
        ax.axvline(2 * np.hypot(model["delta"], model["t"]), color="0.75", lw=0.6)
        ax.set_xlim(*xlim)
        ax.set_ylim(0, 6.5)
        ax.set_xlabel(r"$\hbar\omega$ (eV)")
        ax.set_ylabel(r"Re $\sigma\,(e^2/h)$")

    def imag_part(ax, m, model, legend):
        w = m["photon"][2:]
        ax.axhline(0, color="0.8", lw=0.5)
        ax.plot(w, m["im_xx"][2:], color=k0, lw=1.2, label=r"Im $\sigma$, moments")
        ax.plot(w[::15], m["im_kk"][2:][::15], mec=k3, label=r"KK of Re $\sigma$", **marker)
        ax.axvline(2 * np.hypot(model["delta"], model["t"]), color="0.75", lw=0.6)
        ax.set_xlim(0, 17)
        ax.set_xlabel(r"$\hbar\omega$ (eV)")
        ax.set_ylabel(r"Im $\sigma\,(e^2/h)$")
        ax.legend(fontsize=8.5, handlelength=1.4, **legend)

    # (a) graphene Re sigma, (b) Im sigma and Kramers-Kronig, (c) the (mu, omega) map.
    real_part(a, g, GRAPHENE, (0, 10), 12, 14)
    a.axhline(np.pi / 2, color="0.5", lw=0.6, ls=":")
    a.legend(loc="upper left", bbox_to_anchor=(0.0, 0.88), fontsize=8.5, handlelength=1.4)
    imag_part(b, g, GRAPHENE, dict(loc="upper right", bbox_to_anchor=(1.0, 0.86)))
    im = c.pcolormesh(mp["photon"], mp["mu"], mp["re"], shading="auto", cmap="Greys",
                      vmin=0, vmax=3.0, rasterized=True)
    c.plot(2 * mp["mu"], mp["mu"], color=k2, lw=1.0, ls="--")
    c.text(1.0, 1.75, r"$\hbar\omega = 2\mu$", color=k2, rotation=45, fontsize=10)
    c.set_xlim(0, 6)
    c.set_ylim(mp["mu"][0], mp["mu"][-1])
    c.set_xlabel(r"$\hbar\omega$ (eV)")
    c.set_ylabel(r"$\mu$ (eV)")
    cb = fig.colorbar(im, ax=c, pad=0.02)
    cb.set_label(r"Re $\sigma\,(e^2/h)$")

    # (d) hBN Re sigma, (e) Im sigma and Kramers-Kronig, (f) polarisability below the gap.
    real_part(d, h, HBN, (6, 17), 15, 0)
    d.axvline(2 * HBN["delta"], color="0.75", lw=0.6)
    d.legend(loc="upper right", fontsize=8.5, handlelength=1.4)
    imag_part(e, h, HBN, dict(loc="lower right"))
    f.plot(h["alpha_photon"][2:], h["alpha_moments"][2:], color=k0, lw=1.2, label="moments")
    f.plot(h["alpha_photon"][2::15], h["alpha_exact"][2::15], mec=k2, label="exact (clean)", **marker)
    f.set_xlim(0, 2 * HBN["delta"])
    f.set_xlabel(r"$\hbar\omega$ (eV)")
    f.set_ylabel(r"$\alpha_\mathrm{2D} = -e^2\,\mathrm{Im}\,\sigma/\hbar\omega$ ($\mathrm{\AA}$)")
    f.legend(loc="upper left", fontsize=8.5)

    tags = zip((a, b, c, d, e, f), "abcdef",
               ("upper left", "upper right", "upper left", "upper left", "upper right", "lower right"))
    for ax, tag, loc in tags:
        prl.panel(ax, f"({tag})", loc=loc, bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
    fig.tight_layout()
    fig.savefig(os.path.join(_HERE, "fig_benchmark.png"), dpi=200)


if __name__ == "__main__":
    main(*sys.argv[1:])
