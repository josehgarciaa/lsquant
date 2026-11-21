#!/usr/bin/env python3
"""Exact reference for example 19: diagonalise H(k) on the L x L mesh.

An L x L supercell, periodic in both directions, contains exactly the k points
k = (i/L) b1 + (j/L) b2 of an L x L mesh, so the Bloch blocks of
make_kane_mele.py give the exact spectrum of the supercell that the KPM run
uses. We compute:

  - the band edges around E = 0 and the global gap;
  - the charge Hall conductivity sigma_xy(E) and the spin Hall conductivity
    sigma^z_xy(E) (e^2/h, Pauli sigma_z) from the Berry curvature, with the
    occupations smeared by a Gaussian of standard deviation eta;
  - the spin Chern number at E_F = 0 (lsquant.oracle.spin_hall.spin_chern);
  - the Kubo-Greenwood sigma_xx(E) (e^2/h), the DOS and the spin DOS per cell,
    with the delta functions broadened by the same Gaussian.

Usage: python oracle.py 16 [eta_meV]  ->  operators/km_L16.oracle.json
"""
import json
import os
import sys

import numpy as np

from lsquant.oracle import spin_hall
from make_kane_mele import CELL_AREA, bloch_blocks, bond_vector

_HERE = os.path.dirname(os.path.abspath(__file__))


def bloch_matrices(L):
    blocks = bloch_blocks()
    grid = np.arange(L) / L
    k1, k2 = [a.ravel() for a in np.meshgrid(grid, grid, indexing="xy")]   # k1 fastest
    H = np.zeros((L * L, 4, 4), complex)
    VX, VY = np.zeros_like(H), np.zeros_like(H)
    for r, block in blocks.items():
        phase = np.exp(2j * np.pi * (k1 * r[0] + k2 * r[1]))[:, None, None]
        d = np.array([[bond_vector(r, a, b) for b in range(4)] for a in range(4)])
        H += phase * block
        VX += phase * (1j * block * d[..., 0])
        VY += phase * (1j * block * d[..., 1])
    return H, VX, VY


def main(L=16, eta_meV=20.0):
    H, VX, VY = bloch_matrices(L)
    E, C = np.linalg.eigh(H)
    sz = np.diag([1.0, -1.0, 1.0, -1.0]).astype(complex)
    JXSZ = 0.5 * (VX @ sz + sz @ VX)
    eta = eta_meV / 1000.0
    energies = np.linspace(-4.0, 4.0, 1601)
    omega_c = spin_hall.spin_berry_curvature(E, C, VX, VY)
    omega_s = spin_hall.spin_berry_curvature(E, C, JXSZ, VY)
    sigma_c = spin_hall.sigma_xy(E, omega_c, energies, CELL_AREA, broadening=eta)
    sigma_s = spin_hall.sigma_xy(E, omega_s, energies, CELL_AREA, broadening=eta)
    sigma_s0 = float(spin_hall.sigma_xy(E, omega_s, [0.0], CELL_AREA)[0])
    sigma_c0 = float(spin_hall.sigma_xy(E, omega_c, [0.0], CELL_AREA)[0])
    gauss = np.exp(-0.5 * ((energies[:, None, None] - E[None]) / eta) ** 2) / (eta * np.sqrt(2 * np.pi))
    spin_weight = np.einsum("kan,a,kan->kn", C.conj(), np.diag(sz).real, C).real   # <n k|sigma_z|n k>
    dos = gauss.sum(axis=(1, 2)) / (L * L)
    # Kubo-Greenwood: sigma_xx = 2 pi^2 / (N_k A) sum_k sum_nm |VX_nm|^2 g(E-E_n) g(E-E_m),
    # with g the same Gaussian as the KPM reconstruction (e^2/h; V in eV Angstrom).
    vx_eig = np.abs(np.einsum("kan,kab,kbm->knm", C.conj(), VX, C)) ** 2
    sigma_xx = 2 * np.pi ** 2 / (L * L * CELL_AREA) * np.einsum("ekn,knm,ekm->e", gauss, vx_eig, gauss)
    spin_dos = (gauss * spin_weight[None]).sum(axis=(1, 2)) / (L * L)
    vbm, cbm = float(E[:, 1].max()), float(E[:, 2].min())
    chern = spin_hall.spin_chern(E, C, np.broadcast_to(np.eye(4), H.shape), sz, 0.0, (L, L, 1))
    out = dict(L=L, eta_meV=eta_meV, cell_area_A2=CELL_AREA,
               vbm_eV=vbm, cbm_eV=cbm, gap_eV=cbm - vbm,
               sigma_xy_EF0=sigma_c0, sigma_z_xy_EF0=sigma_s0,
               spin_chern={k: (v if isinstance(v, (int, str, type(None))) else float(v))
                           for k, v in chern.items()},
               energies_eV=energies.tolist(), sigma_xy=sigma_c.tolist(), sigma_z_xy=sigma_s.tolist(),
               sigma_xx=sigma_xx.tolist(), dos_per_cell=dos.tolist(), spin_dos_per_cell=spin_dos.tolist(),
               units=dict(sigma="e^2/h (spin: Pauli sigma_z)", dos="states / eV / cell (4 orbitals)"))
    path = os.path.join(_HERE, "operators", f"km_L{L}.oracle.json")
    with open(path, "w") as f:
        json.dump(out, f)
    print(f"gap {cbm - vbm:.4f} eV ({vbm:.4f} to {cbm:.4f}); sigma_xy(0) = {sigma_c0:.4f}; "
          f"sigma^z_xy(0) = {sigma_s0:.4f} e^2/h; spin Chern {chern['spin_chern']:.4f}, Z2 {chern['z2']}")
    print("wrote", path)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 16,
         float(sys.argv[2]) if len(sys.argv) > 2 else 20.0)
