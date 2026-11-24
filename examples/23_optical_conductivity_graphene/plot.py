#!/usr/bin/env python3
"""Figure of the optical conductivity of graphene (example 23).

Usage:
    python plot.py <run_out_dir>
Reads sigma_mu0.dat written by optical.py and writes fig_optical.png next to this script.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import prl_style as prl


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    data = np.loadtxt(os.path.join(sys.argv[1], "sigma_mu0.dat"))
    photon = data[:, 0]
    prl.use("web", height=5.6)
    fig, (top, bot) = plt.subplots(2, 1, sharex=True)

    top.plot(photon, data[:, 1], label="moments", **prl.trace(0))
    top.plot(photon[::6], data[::6, 3], label="exact", lw=0, **prl.trace(1, line=False, marker=True))
    top.axhline(np.pi / 2, color="0.4", lw=0.8, ls=":")
    top.set_ylabel(r"Re $\sigma_{xx}$ ($e^2/h$)")
    top.set_ylim(0, 2.6)
    top.legend(loc="lower right")
    prl.panel(top, "(a)")

    for i, (column, label) in enumerate(((1, "0"), (4, "0.5"), (5, "1.0"))):
        bot.plot(photon, data[:, column], label=rf"$\mu = {label}$ eV", **prl.trace(i))
    bot.set_xlabel(r"photon energy $\hbar\omega$ (eV)")
    bot.set_ylabel(r"Re $\sigma_{xx}$ ($e^2/h$)")
    bot.set_ylim(0, 2.6)
    bot.set_xlim(0, photon.max())
    bot.legend(loc="lower right")
    prl.panel(bot, "(b)")
    prl.save(fig, os.path.join(here, "fig_optical"))


if __name__ == "__main__":
    main()
