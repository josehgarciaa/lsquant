#!/usr/bin/env python3
"""Optical conductivity of graphene from velocity-velocity Chebyshev moments, and its check.

The run file computes the rank-2 moments of Tr[VX delta(E-H) VX dG-(E)/dE] once.
`lsquant.reconstruct.optical_conductivity` turns them into sigma_xx(mu, omega)
for every chemical potential and photon energy asked. This script evaluates it
at mu = 0 and compares it with the exact interband conductivity of the same
tight-binding model, integrated over the Brillouin zone:

    Re sigma(omega) = g_s pi / (omega A_cell N_k) sum_k |<c k| dH/dk_x |v k>|^2 delta(omega - 2 |f(k)|)

in e^2/hbar, times 2 pi for e^2/h. For photon energies well below the hopping
it tends to the universal value e^2 / (4 hbar) = (pi / 2) e^2/h.

Usage:
    python optical.py <run_out_dir> [cells_per_side]
Writes sigma_mu0.dat: photon energy; Re and Im sigma from the moments and the exact Re sigma
at mu = 0; Re sigma at mu = 0.5 and 1.0 eV.
"""
import glob
import os
import sys

import numpy as np

from lsquant import _core
from lsquant.reconstruct import optical_conductivity

from make_graphene import A1, A2, CELL_AREA, DELTA, T_HOP

SPIN = 2
ENERGY_WINDOW = (-8.5, 8.5)        # eV: contains the spectrum, |E| <= 3 |t| = 8.1 eV
PHOTON = np.linspace(0.05, 3.9, 155)   # eV: 8.5 + 3.9 < 12.5, the window of the moments
DOPED = (0.5, 1.0)                  # eV: chemical potentials of the doped curves


def exact_interband(photon, width, n_k=1200):
    """Exact Re sigma_xx(omega) at mu = 0 in e^2/h; delta broadened by a Gaussian of `width` eV."""
    u = (np.arange(n_k) + 0.5) / n_k
    u1, u2 = np.meshgrid(u, u, indexing="ij")
    recip = 2 * np.pi * np.linalg.inv(np.array([A1, A2])).T       # rows: b1, b2
    kx = u1 * recip[0, 0] + u2 * recip[1, 0]
    ky = u1 * recip[0, 1] + u2 * recip[1, 1]
    f = np.zeros_like(kx, dtype=complex)
    dfx = np.zeros_like(f)
    for d1, d2 in ((0, 0), (-1, 0), (0, -1)):
        bond = DELTA + d1 * A1 + d2 * A2
        phase = T_HOP * np.exp(1j * (kx * bond[0] + ky * bond[1]))
        f += phase
        dfx += 1j * bond[0] * phase
    gap = 2 * np.abs(f).ravel()
    # |<c| dH/dk_x |v>|^2 for H = [[0, f], [f*, 0]]: (Im(f* df/dk_x) / |f|)^2
    weight = (np.imag(np.conj(f) * dfx) / np.abs(f)).ravel() ** 2
    out = np.empty(photon.size)
    for i, w in enumerate(photon):
        delta = np.exp(-0.5 * ((w - gap) / width) ** 2) / (width * np.sqrt(2 * np.pi))
        out[i] = SPIN * np.pi / (w * CELL_AREA * n_k * n_k) * np.sum(weight * delta)
    return 2 * np.pi * out


def main():
    out_dir = sys.argv[1]
    cells = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    path = sorted(glob.glob(os.path.join(out_dir, "*.chebmom")))[0]
    moments = _core.ReadChebmom(path)
    energies, sigma = optical_conductivity(
        moments, PHOTON, energy_window=ENERGY_WINDOW, num_cells=cells * cells,
        cell_area=CELL_AREA, spin_degeneracy=SPIN)
    at_zero = sigma[:, np.argmin(np.abs(energies))]

    half_width = 0.5 * (moments.window.upper - moments.window.lower)
    resolution = np.pi * half_width / moments.dims[0]
    exact = exact_interband(PHOTON, width=resolution)
    doped = [sigma[:, np.argmin(np.abs(energies - mu))].real for mu in DOPED]
    np.savetxt(os.path.join(out_dir, "sigma_mu0.dat"),
               np.column_stack((PHOTON, at_zero.real, at_zero.imag, exact, *doped)),
               header="photon_energy_eV Re_sigma Im_sigma exact_Re_sigma (mu = 0), then Re_sigma at mu = "
                      + ", ".join(f"{mu:g}" for mu in DOPED) + " eV; all in e^2/h")

    mask = (PHOTON > 1.0) & (PHOTON < 3.5)
    deviation = np.max(np.abs(at_zero.real[mask] / exact[mask] - 1))
    plateau = np.mean(at_zero.real[(PHOTON > 0.8) & (PHOTON < 1.2)])
    print(f"moments {moments.dims[0]} x {moments.dims[1]}, resolution {1000 * resolution:.0f} meV")
    print(f"Re sigma(1 eV) = {plateau:.4f} e^2/h   (universal value pi/2 = {np.pi / 2:.4f})")
    print(f"largest deviation from the exact interband curve, 1.0 to 3.5 eV: {100 * deviation:.1f} %")
    ok = deviation < 0.05
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
