#!/usr/bin/env python3
"""Optical conductivity of monolayer hBN from seven DFT models: Chebyshev moments against the exact sums.

Usage (from this directory, after the LSQuant runs and sos.py):
    python analyze.py OUT_ROOT SOS_DIR

OUT_ROOT holds the output directories of the run files (out/hbn_siesta_pbe, ...),
SOS_DIR the sum-over-states references written by sos.py. For every model
the script reconstructs sigma_xx and sigma_yy at a chemical potential in the
middle of the gap and checks:

1. Re sigma and Im sigma from the moments against the exact sum over states
   of the same model, mesh, velocity and window (sos.py), up to 14 eV.
2. Kramers-Kronig: Im sigma from the moments against the transform of
   Re sigma from the moments. In these models the absorption continues
   above the largest photon energy the window allows (sigma -> sigma*), so
   the transform, truncated at 15 eV, misses a tail. Below 10 eV that tail
   is smooth and odd in w, -(2w/pi) int_W^inf Re sigma/(w'^2 - w^2) dw' =
   a w + b w^3 + ...; it is fitted and removed, and what remains must be
   small: every peak of Re sigma has its feature in Im sigma.
3. Isotropy: sigma_xx = sigma_yy for a hexagonal crystal.
4. Onset: the photon energy where Re sigma reaches half of its value 0.3 eV
   above the smallest direct gap of the model on the mesh (at K). In 2D the
   absorption starts with a step at that gap, broadened by the resolution.
5. Insulator: Im sigma < 0 below the gap.

and reports the static 2D polarisability alpha_2D = -e^2 Im sigma / (hbar w),
the first absorption peak and the gap of every model. Writes analysis.json
and curves.npz (read by plot.py).
"""
import glob
import json
import os
import sys

import numpy as np

from lsquant import _core
from lsquant.reconstruct import optical_conductivity

E2 = 14.399645                     # e^2 / (4 pi eps0) in eV Angstrom
PHOTON = np.concatenate(([1e-3], np.arange(0.02, 15.0 + 1e-9, 0.02)))
CELL_AREA = 2.504**2 * np.sqrt(3) / 2
# name: (label, spin degeneracy, integration range of the reconstruction in eV). The lower
# limit must lie below every state the moments hold: in HSE06 the N 2s band of PySCF starts
# below -21 eV, where the shift projection piles the states outside its window; a limit
# through that pile cuts the Fermi sea and gives a spurious response at w -> 0.
MODELS = {
    "siesta_pbe": ("SIESTA, PBE, DZP", 1, (-22.0, 3.0)),
    "pyscf_pbe": ("PySCF, PBE, DZVP", 1, (-22.0, 3.0)),
    "pyscf_hse": ("PySCF, HSE06, DZVP", 1, (-22.0, 3.0)),
    "w90_pbe": ("QE + Wannier90, PBE", 2, (-22.0, 3.0)),
    "w90_hse": ("QE + Wannier90, HSE06", 2, (-22.0, 3.0)),
    "pao_pbe": ("QE + PAOFLOW, PBE", 2, (-22.0, 3.0)),
    "pao_hse": ("QE + PAOFLOW, HSE06", 2, (-22.0, 3.0)),
}
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def kramers_kronig(photon, re_sigma):
    """Im sigma(w) = -(2 w / pi) P int_0^W Re sigma(w') / (w'^2 - w^2) dw', subtracted principal value.

    Valid where Re sigma has vanished before the upper limit W; NaN at W.
    """
    w, f = np.asarray(photon), np.asarray(re_sigma)
    top, slope = w[-1], np.gradient(f, w)
    out = np.full(len(w), np.nan)
    for i, w0 in enumerate(w[:-1]):
        den = w**2 - w0**2
        g = np.where(np.abs(den) > 1e-12, (f - f[i]) / np.where(den == 0, 1, den), 0.0)
        g[i] = slope[i] / (2 * w0)
        out[i] = -2 * w0 / np.pi * (_trapezoid(g, w) + f[i] * np.log(abs((top - w0) / (top + w0))) / (2 * w0))
    return out


def moment_file(out_dir, axis):
    """The moment file of Tr[V{axis} delta(E-H) V{axis} dG-(E)/dE] in an output directory.

    The metric executor names it after the operators (...VX-VX...chebmom);
    the projected executor numbers the words (word_<i>.chebmom) in the order
    of the run file, recorded in run.resolved.json.
    """
    named = sorted(glob.glob(os.path.join(out_dir, f"*V{axis}-V{axis}*.chebmom")))
    if named:
        return named[0]
    words = json.load(open(os.path.join(out_dir, "run.resolved.json")))["words"]
    index = next(i for i, w in enumerate(words) if w.replace(" ", "").startswith(f"Tr[V{axis}delta"))
    return os.path.join(out_dir, f"word_{index}.chebmom")


def sigma(out_dir, axis, spin, energy_window, mu):
    path = moment_file(out_dir, axis)
    moments = _core.ReadChebmom(path)
    energies, s = optical_conductivity(moments, PHOTON, energy_window=energy_window,
                                       num_cells=120 * 120, cell_area=CELL_AREA, spin_degeneracy=spin)
    i = np.argmin(np.abs(energies - mu))
    return s[:, i], energies[i]


def onset(w, re, gap):
    """Half-height point of the absorption step at the direct gap: Re sigma = Re sigma(gap + 0.3) / 2."""
    target = 0.5 * np.interp(gap + 0.3, w, re)
    window = (w > gap - 1.0) & (w <= gap + 0.3)
    ww, rr = w[window], re[window]
    j = np.argmax(rr >= target)
    return float(np.interp(target, [rr[j - 1], rr[j]], [ww[j - 1], ww[j]]))


def tail_removed(w, residual, band):
    """Fit a w + b w^3 to Im sigma - KK[Re sigma] in `band` and return what is left there."""
    basis = np.column_stack((w[band], w[band] ** 3))
    coefficients, *_ = np.linalg.lstsq(basis, residual[band], rcond=None)
    return residual[band] - basis @ coefficients, coefficients


def main():
    out_root, sos_dir = sys.argv[1], sys.argv[2]
    report, curves = {}, {}
    print(f"{'model':24s} {'gap':>6s} {'onset':>6s} {'peak':>6s} {'alpha0':>7s} "
          f"{'|KPM-SOS|':>9s} {'KK':>7s} {'xx-yy':>7s}")
    for name, (label, spin, ewin) in MODELS.items():
        out_dir = next((d for d in (os.path.join(out_root, "hbn_" + name), os.path.join(out_root, name))
                        if glob.glob(os.path.join(d, "*.chebmom"))), None)
        if out_dir is None:
            continue
        ref = np.load(os.path.join(sos_dir, name + ".npz"))
        mu = 0.5 * (float(ref["vbm"]) + float(ref["cbm"]))
        sx, mu_used = sigma(out_dir, "X", spin, ewin, mu)
        try:
            sy, _ = sigma(out_dir, "Y", spin, ewin, mu)
        except (FileNotFoundError, IndexError, RuntimeError, StopIteration):
            sy = np.full_like(sx, np.nan)              # sigma_yy not computed (yet)
        exact = np.interp(PHOTON, ref["photon"], ref["sigma_xx"].real) \
            + 1j * np.interp(PHOTON, ref["photon"], ref["sigma_xx"].imag)
        band = (PHOTON > 0.5) & (PHOTON < 14.0)
        scale = np.abs(exact[band]).max()
        dev = float(np.abs(sx - exact)[band].max() / scale)
        kk = kramers_kronig(PHOTON, sx.real)
        kk_band = (PHOTON > 0.5) & (PHOTON < 10.0)
        left, tail = tail_removed(PHOTON, sx.imag - kk, kk_band)
        kk_dev = float(np.abs(left).max() / np.abs(sx[kk_band]).max())
        kk_tail = PHOTON * tail[0] + PHOTON ** 3 * tail[1]
        aniso = float(np.abs(sx - sy)[band].max() / scale)
        gap = float(ref["direct_gap"])
        start = onset(PHOTON, sx.real, gap)
        first = (PHOTON > gap) & (PHOTON < gap + 2.0)
        peak = float(PHOTON[first][np.argmax(sx.real[first])])
        alpha = -E2 * sx.imag / (2 * np.pi * PHOTON)       # e^2/h -> e^2/hbar: divide by 2 pi
        below = PHOTON < gap - 0.5
        insulator = bool(np.all(sx.imag[below] < 0))
        report[name] = dict(label=label, spin_degeneracy=spin, mu=mu_used, direct_gap=gap,
                            direct_gap_k=ref["direct_gap_k"].tolist(), onset_half_height=start,
                            first_peak=peak, alpha_2d_static=float(alpha[0]),
                            kpm_vs_sos_max_relative=dev, kramers_kronig_max_relative=kk_dev,
                            anisotropy_max_relative=aniso, im_sigma_negative_below_gap=insulator)
        curves[name] = dict(re=sx.real, im=sx.imag, re_yy=sy.real, kk=kk + kk_tail, sos_re=exact.real,
                            sos_im=exact.imag, alpha=alpha)
        print(f"{label:24s} {gap:6.3f} {start:6.3f} {peak:6.2f} {alpha[0]:7.4f} "
              f"{dev:9.4f} {kk_dev:7.4f} {aniso:7.4f}  {'Im<0' if insulator else 'Im>0 below gap!'}")
    json.dump(report, open("analysis.json", "w"), indent=1)
    np.savez("curves.npz", photon=PHOTON,
             **{f"{n}__{k}": v for n, c in curves.items() for k, v in c.items()})
    return 0


if __name__ == "__main__":
    sys.exit(main())
