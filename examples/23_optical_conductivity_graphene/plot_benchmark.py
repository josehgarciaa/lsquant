#!/usr/bin/env python3
"""Plot the optical tests of graphene and hBN written by benchmark.py.

Usage (from this directory):
    python plot_benchmark.py [data_dir]
Reads bench_graphene.npz, bench_hbn.npz and bench_map.npz (default: data/ next
to this script) and writes fig_benchmark.png next to this script.
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
    prl.use("web", height=6.2)
    fig, axes = plt.subplots(2, 2)
    (a, b), (c, d) = axes
    k0, k1, k2, k3 = (prl.trace(i)["color"] for i in range(4))

    # (a) Graphene, mu = 0: field along x and y, exact, universal value, van Hove line.
    w = g["photon"][2:]
    a.plot(w, g["re_xx"][2:], color=k0, lw=1.2, label=r"moments, $\sigma_{xx}$")
    a.plot(w, g["re_yy"][2:], color=k1, lw=1.0, ls="--", label=r"moments, $\sigma_{yy}$")
    a.plot(w[14::12], g["exact_re"][2:][14::12], "s", ms=3.5, mfc="none", mec=k2, mew=0.8, label="exact")
    a.axhline(np.pi / 2, color="0.5", lw=0.6, ls=":")
    a.axvline(2 * abs(GRAPHENE["t"]), color="0.75", lw=0.6)
    a.set_xlim(0, 10)
    a.set_ylim(0, 6)
    a.set_xlabel(r"$\hbar\omega$ (eV)")
    a.set_ylabel(r"Re $\sigma\,(e^2/h)$")
    a.legend(loc="upper left", bbox_to_anchor=(0.0, 0.88), fontsize=8.5, handlelength=1.4)

    # (b) Graphene: Re sigma(mu, omega) map with the Pauli edge 2 mu.
    im = b.pcolormesh(mp["photon"], mp["mu"], mp["re"], shading="auto", cmap="Greys",
                      vmin=0, vmax=3.0, rasterized=True)
    b.plot(2 * mp["mu"], mp["mu"], color=k2, lw=1.0, ls="--")
    b.text(1.1, 1.75, r"$\hbar\omega = 2\mu$", color=k2, rotation=38, fontsize=10)
    b.set_xlim(0, 6)
    b.set_ylim(mp["mu"][0], mp["mu"][-1])
    b.set_xlabel(r"$\hbar\omega$ (eV)")
    b.set_ylabel(r"$\mu$ (eV)")
    cb = fig.colorbar(im, ax=b, pad=0.02)
    cb.set_label(r"Re $\sigma\,(e^2/h)$")

    # (c) hBN, mu = 0: absorption, exact, gap and van Hove lines.
    w = h["photon"][2:]
    c.plot(w, h["re_xx"][2:], color=k0, lw=1.2, label=r"moments, $\sigma_{xx}$")
    c.plot(w, h["re_yy"][2:], color=k1, lw=1.0, ls="--", label=r"moments, $\sigma_{yy}$")
    c.plot(w[::12], h["exact_re"][2:][::12], "s", ms=3.5, mfc="none", mec=k2, mew=0.8, label="exact")
    c.axvline(2 * HBN["delta"], color="0.75", lw=0.6)
    c.axvline(2 * np.hypot(HBN["delta"], HBN["t"]), color="0.75", lw=0.6)
    c.set_xlim(6, 10)
    c.set_ylim(0, None)
    c.set_xlabel(r"$\hbar\omega$ (eV)")
    c.set_ylabel(r"Re $\sigma\,(e^2/h)$")
    c.legend(loc="upper left", bbox_to_anchor=(0.0, 0.93), fontsize=9)

    # (d) hBN below the gap: the 2D polarisability, from Im sigma.
    d.plot(h["alpha_photon"][2:], h["alpha_moments"][2:], color=k0, lw=1.2, label="moments")
    d.plot(h["alpha_photon"][2::12], h["alpha_exact"][2::12], "s", ms=3.5, mfc="none", mec=k2, mew=0.8,
           label="exact (clean)")
    d.set_xlim(0, 5)
    d.set_xlabel(r"$\hbar\omega$ (eV)")
    d.set_ylabel(r"$\alpha_\mathrm{2D} = -e^2\,\mathrm{Im}\,\sigma/\hbar\omega$ ($\mathrm{\AA}$)")
    d.legend(loc="upper left", fontsize=9)

    for ax, tag, loc in zip((a, b, c, d), "abcd", ("upper left", "upper left", "upper right", "lower right")):
        prl.panel(ax, f"({tag})", loc=loc,
                  bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
    fig.tight_layout()
    fig.savefig(os.path.join(_HERE, "fig_benchmark.png"), dpi=200)


if __name__ == "__main__":
    main(*sys.argv[1:])
