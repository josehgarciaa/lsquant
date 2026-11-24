#!/usr/bin/env python3
"""Optical tests of graphene and hBN: peaks, onsets, prefactors, polarisation and polarisability.

The tests are those proposed by T. Galvani (private communication, 2026) for
the optical conductivity of the nearest-neighbour honeycomb model, with
on-site energies +Delta (B) and -Delta (N) and hopping t:

1. Graphene, Pauli blocking: at chemical potential mu the absorption starts
   at hbar omega = 2 mu. Measured: the photon energy at which
   Re sigma(mu, omega) reaches half of Re sigma(0, omega).
2. Van Hove peak: the transition between the valence and conduction bands at
   M gives a peak of Re sigma(0, omega) at 2 sqrt(Delta^2 + t^2), for
   graphene (2|t|) and hBN. A peak position does not depend on the
   broadening, so it is tracked precisely (parabola through the maximum).
3. hBN onset: at mu = 0 the absorption starts at the gap 2 Delta. Measured
   as the half-height point of the step, against the same quantity of the
   exact response broadened like the moments.
4. Graphene prefactor: Re sigma(0, omega) at low photon energy is the
   universal pi/2 e^2/h (plus a small band-curvature correction of this
   model), and its limit omega -> 0 is the static Kubo-Bastin value of the
   same moments.
5. hBN polarisability: alpha_2D = -e^2 Im sigma / (hbar omega) below the gap,
   against the clean k-space sum (exact_honeycomb.py). This tests Im sigma,
   its sign and Re epsilon.
6. Polarisation: a hexagonal crystal has sigma_xx = sigma_yy. With the same
   random states the two differ only by the noise of the trace.
7. Kramers-Kronig: Re sigma and Im sigma are not independent. Im sigma from
   the moments is compared with the Kramers-Kronig transform of Re sigma
   from the same moments, over the whole absorption range, so every peak
   of Re sigma must come with its feature in Im sigma.

Usage (from this directory, after the two runs of run_graphene.json and run_hbn.json):
    python benchmark.py ../../out/optical_graphene ../../out/optical_hbn
Writes benchmark.json (every number below) and the curves in bench_*.npz, read by plot_benchmark.py.
"""
import glob
import json
import os
import sys

import numpy as np

from lsquant import _core
from lsquant.reconstruct import kubo_bastin, optical_conductivity

from exact_honeycomb import E2, Honeycomb
from make_honeycomb import GRAPHENE, HBN, lattice

SPIN = 2
CELLS = 400 * 400
PHOTON = np.concatenate(([1e-4, 2e-4], np.arange(0.02, 17.0 + 1e-9, 0.02)))
ENERGY_WINDOW = (-8.6, 8.6)           # eV: contains both spectra (8.10 and 7.87 eV); 8.6 + 17 < 26
KK_RANGE = (0.5, 16.0)               # eV: photon energies where Im sigma is compared with Kramers-Kronig
DOPED = (0.5, 1.0, 1.5)               # eV
MU_MAP = np.arange(0.0, 2.51, 0.05)   # eV: chemical potentials of the (mu, omega) map


_trapezoid = getattr(np, "trapezoid", None) or np.trapz   # NumPy 2 renamed trapz


def kramers_kronig(photon, re_sigma):
    """Im sigma(omega) = -(2 omega/pi) P int_0^W Re sigma(w') / (w'^2 - omega^2) dw'  (exp(-i omega t)).

    Re sigma is given on `photon` (ascending, from near zero to W, beyond the
    largest transition). The principal value is taken by subtracting
    Re sigma(omega): P int_0^W dw' / (w'^2 - omega^2) = ln|(W - omega)/(W + omega)| / (2 omega).
    Trapezoid rule; the subtracted integrand is regular at w' = omega.
    """
    w, f = np.asarray(photon), np.asarray(re_sigma)
    top = w[-1]
    out = np.full(len(w), np.nan)                       # undefined at the upper limit W
    for i, w0 in enumerate(w[:-1]):
        den = w**2 - w0**2
        g = np.where(np.abs(den) > 1e-12, (f - f[i]) / np.where(den == 0, 1, den), 0.0)
        g[i] = np.gradient(f, w)[i] / (2 * w0)            # limit of the subtracted integrand
        principal = _trapezoid(g, w) + f[i] * np.log(abs((top - w0) / (top + w0))) / (2 * w0)
        out[i] = -2 * w0 / np.pi * principal
    return out


def load(out_dir, op):
    path = sorted(glob.glob(os.path.join(out_dir, f"NonEqOp{op}-{op}*.chebmom")))[0]
    return _core.ReadChebmom(path)


def conductivity(moments, model):
    energies, sigma = optical_conductivity(moments, PHOTON, energy_window=ENERGY_WINDOW,
                                           num_cells=CELLS, cell_area=lattice(model["a"])[2],
                                           spin_degeneracy=SPIN)
    return energies, sigma


def at(energies, sigma, mu):
    """sigma(omega) at chemical potential mu, linear in mu between grid points."""
    return np.array([np.interp(mu, energies, row.real) + 1j * np.interp(mu, energies, row.imag)
                     for row in sigma])


def peak(photon, values, near, half=0.6):
    """Photon energy of the maximum within +-half of `near`: parabola through the three top points."""
    sel = np.where(np.abs(photon - near) < half)[0]
    i = sel[np.argmax(values[sel])]
    y0, y1, y2 = values[i - 1:i + 2]
    return photon[i] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) * (photon[i + 1] - photon[i])


def crossing(photon, values, level, start, stop):
    """First photon energy in [start, stop] at which values rises through level (linear)."""
    sel = np.where((photon >= start) & (photon <= stop))[0]
    for i, j in zip(sel[:-1], sel[1:]):
        if values[i] < level <= values[j]:
            return photon[i] + (level - values[i]) / (values[j] - values[i]) * (photon[j] - photon[i])
    return float("nan")


def main():
    out_g, out_h = sys.argv[1], sys.argv[2]
    report, curves, checks = {}, {}, []

    def check(name, value, target, tol, unit):
        ok = abs(value - target) <= tol
        checks.append(dict(test=name, value=float(value), target=float(target), tolerance=float(tol),
                           unit=unit, passed=bool(ok)))
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {value:.4f} {unit} (target {target:.4f}, tolerance {tol:.4f})")

    for name, model, out_dir in (("graphene", GRAPHENE, out_g), ("hbn", HBN, out_h)):
        mx, my = load(out_dir, "VX"), load(out_dir, "VY")
        half_width = 0.5 * (mx.window.upper - mx.window.lower)
        res = np.pi * half_width / mx.dims[0]           # Jackson resolution of one expansion
        width = np.sqrt(2) * res                        # a transition: two expansions (delta and G)
        energies, sx = conductivity(mx, model)
        _, sy = conductivity(my, model)
        exact = Honeycomb(model["delta"], model["t"], model["a"])
        s0x, s0y = at(energies, sx, 0.0), at(energies, sy, 0.0)
        ref0 = exact.re_sigma(PHOTON[2:], 0.0, width)
        curves[name] = dict(photon=PHOTON, re_xx=s0x.real, im_xx=s0x.imag, re_yy=s0y.real, im_yy=s0y.imag,
                            exact_re=np.concatenate(([np.nan, np.nan], ref0)))
        report[name] = dict(moments=int(mx.dims[0]), window=[mx.window.lower, mx.window.upper],
                            resolution_eV=res, transition_width_eV=width)
        vh = 2 * np.hypot(model["delta"], model["t"])
        e_vh = peak(PHOTON, s0x.real, vh)
        e_vh_y = peak(PHOTON, s0y.real, vh)
        e_vh_exact = peak(PHOTON[2:], ref0, vh)
        check(f"{name}: van Hove peak, field along x", e_vh, vh, 0.02, "eV")
        check(f"{name}: van Hove peak, field along y", e_vh_y, vh, 0.02, "eV")
        report[name].update(van_hove_target=vh, van_hove_x=e_vh, van_hove_y=e_vh_y, van_hove_exact=e_vh_exact)
        band = (PHOTON > 0.3) & (PHOTON < 9.5)
        aniso = np.max(np.abs(s0x.real - s0y.real)[band]) / np.max(s0x.real[band])
        report[name]["anisotropy_max_relative"] = aniso
        check(f"{name}: |Re sigma_xx - Re sigma_yy| / max Re sigma, 0.3 to 9.5 eV", aniso, 0.0, 0.03, "")
        # 7. Kramers-Kronig consistency of the two parts from the same moments.
        kk = kramers_kronig(PHOTON, s0x.real)
        sel = (PHOTON > KK_RANGE[0]) & (PHOTON < KK_RANGE[1])
        kk_dev = np.max(np.abs(kk - s0x.imag)[sel]) / np.max(np.abs(s0x)[sel])
        curves[name]["im_kk"] = kk
        report[name]["kramers_kronig_max_relative"] = kk_dev
        check(f"{name}: |Im sigma - KK[Re sigma]| / max |sigma|, {KK_RANGE[0]} to {KK_RANGE[1]} eV",
              kk_dev, 0.0, 0.03, "")

        if name == "graphene":
            # 4. prefactor: universal value, model value and continuity with the static curve.
            low = (PHOTON > 0.5) & (PHOTON < 1.5)
            plateau = np.mean(s0x.real[low])
            plateau_exact = np.mean(ref0[(PHOTON[2:] > 0.5) & (PHOTON[2:] < 1.5)])
            check("graphene: mean Re sigma(0, omega), 0.5 to 1.5 eV, against the exact model",
                  plateau, plateau_exact, 0.02 * plateau_exact, "e^2/h")
            e_dc, kb = kubo_bastin(mx, 1.0)
            dc = SPIN * 2 * np.pi * np.interp(0.0, e_dc, kb) / (CELLS * lattice(model["a"])[2])
            check("graphene: Re sigma(0, 1e-4 eV) against the static Kubo-Bastin value of the same moments",
                  s0x.real[0], dc, 0.01 * dc, "e^2/h")
            report[name].update(plateau=plateau, plateau_exact=plateau_exact, universal=np.pi / 2,
                                optical_at_1em4_eV=s0x.real[0], static_kubo_bastin=dc)
            # 1. Pauli blocking.
            onsets = {}
            for mu in DOPED:
                ratio = at(energies, sx, mu).real / np.where(s0x.real > 0, s0x.real, np.nan)
                w_half = crossing(PHOTON, ratio, 0.5, mu, 2 * mu + 1.0)
                ref_ratio = exact.re_sigma(PHOTON[2:], mu, width) / np.where(ref0 > 0, ref0, np.nan)
                w_ref = crossing(PHOTON[2:], ref_ratio, 0.5, mu, 2 * mu + 1.0)
                onsets[mu] = dict(moments=w_half, exact=w_ref)
                check(f"graphene: half-absorption point at mu = {mu} eV", w_half, 2 * mu, res, "eV")
                curves[name][f"re_mu{mu}"] = at(energies, sx, mu).real
                curves[name][f"exact_mu{mu}"] = np.concatenate(([np.nan, np.nan],
                                                                exact.re_sigma(PHOTON[2:], mu, width)))
            report[name]["pauli_onsets"] = {str(k): v for k, v in onsets.items()}
            sel = PHOTON >= 0.02
            curves["map"] = dict(mu=MU_MAP, photon=PHOTON[sel],
                                 re=np.array([at(energies, sx, mu).real[sel] for mu in MU_MAP]))
        else:
            # 3. onset at the gap 2 Delta.
            gap = 2 * model["delta"]
            above = (PHOTON > gap + 0.2) & (PHOTON < gap + 0.4)
            w_on = crossing(PHOTON, s0x.real, 0.5 * np.mean(s0x.real[above]), gap - 1.0, gap + 0.5)
            ex_full = np.concatenate(([0.0, 0.0], ref0))
            w_on_exact = crossing(PHOTON, ex_full, 0.5 * np.mean(ex_full[above]), gap - 1.0, gap + 0.5)
            check("hbn: half-height point of the absorption step against the exact step", w_on, w_on_exact,
                  0.25 * res, "eV")
            report[name].update(gap=gap, onset_moments=w_on, onset_exact=w_on_exact)
            # 5. polarisability and Im sigma below the gap.
            below = PHOTON < gap - 0.3
            alpha_kpm = -E2 * s0x.imag[below] / (2 * np.pi * PHOTON[below])
            alpha_exact = exact.polarisability(PHOTON[below])
            curves[name].update(alpha_photon=PHOTON[below], alpha_moments=alpha_kpm, alpha_exact=alpha_exact)
            for w0 in (1e-4, 2.0, 4.0):
                i = np.argmin(np.abs(PHOTON[below] - w0))
                check(f"hbn: alpha_2D at {PHOTON[i]:g} eV against the clean k-space sum",
                      alpha_kpm[i], alpha_exact[i], 0.02 * alpha_exact[i], "Angstrom")
            report[name].update(alpha_2d_static_moments=alpha_kpm[0], alpha_2d_static_exact=alpha_exact[0],
                                eps_static_d_3_33=1 + 4 * np.pi * alpha_kpm[0] / 3.33)
            assert np.all(s0x.imag[PHOTON < gap - 1.0] < 0), "Im sigma of the insulator must be negative"
            checks.append(dict(test="hbn: Im sigma < 0 below the gap (exp(-i omega t))", passed=True))
            print("PASS  hbn: Im sigma < 0 at every photon energy below 2 Delta - 1 eV")

    report["checks"] = checks
    with open("benchmark.json", "w") as fh:
        json.dump(report, fh, indent=1, default=float)
    for name, data in curves.items():
        np.savez(f"bench_{name}.npz", **data)
    failed = [c["test"] for c in checks if not c["passed"]]
    print(f"\n{len(checks) - len(failed)} of {len(checks)} checks pass")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
