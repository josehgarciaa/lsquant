#!/usr/bin/env python3
"""Exact optical conductivity of one hBN DFT model: sum over states on the k mesh of its run file.

Usage:
    python sos.py RUN.json SPIN_DEGENERACY OUT.npz [--width WIDTH_eV]

The reference that the Chebyshev result of the same run file must reproduce.
It uses the same model, k mesh, velocity operator and projection window as
the LSQuant run, diagonalises H(k) on every k point, and evaluates the
interband Kubo formula of a clean insulator at zero temperature
(time dependence exp(-i omega t)), in units of e^2/hbar,

    sigma_ab(w) = (i g_s / A) sum_k sum_{v,c} v^a_vc v^b_cv / D [L(w - D) + L(w + D)],

with D = E_c - E_v, the velocities hbar v in eV Angstrom, A the area of
the whole sheet (number of k points times the cell area) and L the complex
line shape of a Gaussian of standard deviation `width`,
L(x) = int G(x') / (x - x' + i0) dx' = -i pi w(x / (sqrt(2) width)) / (sqrt(2 pi) width),
with w the Faddeeva function. Im L = -pi G gives Re sigma; Re L, the
principal value, gives Im sigma, so both parts come from the same lines and
are Kramers-Kronig consistent by construction. The default width is
sqrt(2) times the Jackson resolution pi (window width / 2) / M of the run
file, the width of one transition in the Chebyshev reconstruction. The
result is written in e^2/h.

Only states inside projection.window (when the run file has one) take part,
as in the projected LSQuant run: the velocity is projected on those states.
The occupied states are the bands below zero energy at Gamma (the energy
origin of the model sits in the gap or at the valence-band maximum).
"""
import argparse
import json

import numpy as np
from scipy.special import wofz

from lsquant.bases import BlochMetricBasis
from lsquant.models.registry import load_runfile_model

PHOTON = np.arange(0.01, 15.0 + 1e-9, 0.01)     # eV


def line_shape(x, width):
    """Complex Gaussian line L(x) = P int G/(x - x') - i pi G(x)."""
    return -1j * np.pi * wofz(x / (np.sqrt(2) * width)) / (np.sqrt(2 * np.pi) * width)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run")
    ap.add_argument("spin_degeneracy", type=int)
    ap.add_argument("out")
    ap.add_argument("--width", type=float, help="Gaussian standard deviation, eV")
    a = ap.parse_args()
    rf = json.load(open(a.run))
    if a.width is None:
        lo_w, hi_w = rf["expansion"]["window"]
        a.width = np.sqrt(2) * np.pi * (hi_w - lo_w) / 2 / rf["expansion"]["moments"]
    mode = rf["model"].get("params", {}).get("velocity_mode", "berry_connection")
    model = load_runfile_model(rf["model"], mode)
    basis = BlochMetricBasis(model, rf["basis"]["supercell"], velocity_mode=mode)
    E, C = basis.physical_eigenvalues, basis.physical_eigenvectors       # C^+ S C = 1
    lo, hi = rf.get("projection", {}).get("window", (-np.inf, np.inf))
    inside = (E >= lo) & (E <= hi)
    # Occupation by band index: the bands below zero at Gamma (the first k point),
    # so a valence maximum a few meV above the energy origin stays occupied.
    n_occ = int(np.sum(E[0] < 0))
    band = np.arange(E.shape[1])[None, :]
    occ, emp = inside & (band < n_occ), inside & (band >= n_occ)
    area = basis.num_cells * abs(np.linalg.det(np.asarray(model.uc)[:2, :2]))
    vel = {ax: np.conj(C.swapaxes(-1, -2)) @ basis.operator_blocks("V" + ax) @ C for ax in "XY"}
    sigma = {ab: np.zeros(len(PHOTON), complex) for ab in ("xx", "yy")}
    for k in range(basis.num_cells):
        v_idx, c_idx = np.where(occ[k])[0], np.where(emp[k])[0]
        delta = (E[k, c_idx][None, :] - E[k, v_idx][:, None]).ravel()
        shape = line_shape(PHOTON[:, None] - delta, a.width) + line_shape(PHOTON[:, None] + delta, a.width)
        for ab in sigma:
            v = vel[ab[0].upper()][k][np.ix_(v_idx, c_idx)]
            sigma[ab] += shape @ (np.abs(v.ravel()) ** 2 / delta)
    for ab in sigma:
        sigma[ab] *= 1j * a.spin_degeneracy / area * 2 * np.pi          # e^2/hbar -> e^2/h
    # Band edges on the mesh: direct gap at every k and its minimum.
    vbm = np.where(occ, E, -np.inf).max(axis=1)
    cbm = np.where(emp, E, np.inf).min(axis=1)
    kmin = int(np.argmin(cbm - vbm))
    np.savez(a.out, photon=PHOTON, sigma_xx=sigma["xx"], sigma_yy=sigma["yy"], width=a.width,
             spin_degeneracy=a.spin_degeneracy, direct_gap=float((cbm - vbm)[kmin]),
             direct_gap_k=basis.k_points[kmin], vbm=float(vbm.max()), cbm=float(cbm.min()))
    print(f"{a.run}: direct gap {float((cbm - vbm)[kmin]):.3f} eV at k = {basis.k_points[kmin]}, "
          f"max Re sigma_xx {sigma['xx'].real.max():.3f} e^2/h at {PHOTON[np.argmax(sigma['xx'].real)]:.2f} eV")


if __name__ == "__main__":
    main()
