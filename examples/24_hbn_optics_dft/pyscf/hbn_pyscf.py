#!/usr/bin/env python3
"""hBN monolayer with PySCF (Gaussian basis, GTH pseudopotentials), PBE or HSE06.

Usage:
    python hbn_pyscf.py pbe|hse OUT.npz [--mesh 12] [--basis gth-dzvp]

Same cell and positions as ../qe and ../siesta (a = 2.504 Angstrom, 20 Angstrom
along z). Runs a spin-restricted KRKS calculation (PySCF's two-component KGKS has
no range-separated hybrids) on a Gamma-centred mesh x mesh x 1, with
Gaussian density fitting, and writes the LSQuant model with
lsquant.models.pyscf_lcao.write_pyscf_model, which lifts the Kohn-Sham
matrix to two identical spin blocks. The model stores spinors, so optical
runs use spin_degeneracy = 1.
"""
import argparse
import time

import numpy as np
from pyscf.pbc import dft, gto

from lsquant.models.pyscf_lcao import write_pyscf_model

A, C = 2.504, 20.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("xc", choices=["pbe", "hse"])
    ap.add_argument("out")
    ap.add_argument("--mesh", type=int, default=12)
    ap.add_argument("--basis", default="gth-dzvp")
    a = ap.parse_args()
    cell = gto.Cell()
    cell.a = np.array([[A, 0, 0], [-A / 2, A * 3**0.5 / 2, 0], [0, 0, C]])
    frac = np.array([[1 / 3, 2 / 3, 0.5], [2 / 3, 1 / 3, 0.5]])
    pos = frac @ cell.a
    cell.atom = [["B", pos[0]], ["N", pos[1]]]
    cell.unit = "Angstrom"
    cell.basis = a.basis
    cell.pseudo = "gth-pbe"
    cell.verbose = 4
    cell.build()
    kpts = cell.make_kpts([a.mesh, a.mesh, 1])
    mf = dft.KRKS(cell, kpts).density_fit()
    mf.xc = "PBE" if a.xc == "pbe" else "HSE06"
    mf.conv_tol = 1e-9
    t0 = time.time()
    mf.kernel()
    print(f"SCF {a.xc}: converged={mf.converged}, E={mf.e_tot:.8f} Ha, {time.time() - t0:.0f} s")
    info = write_pyscf_model(mf, a.out, mesh=[a.mesh, a.mesh, 1])
    print("model:", info)


if __name__ == "__main__":
    main()
