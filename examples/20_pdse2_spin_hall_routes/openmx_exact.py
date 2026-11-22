"""Exact spin Hall conductivity of an OpenMX model on a dense k mesh.

Usage:
    python openmx_exact.py SCFOUT OUT_PREFIX [--mesh 60] [--eta 0.006] [--broadening 0.006]

Reads the scfout file with the LSQuant OpenMX reader, builds the covariant
velocity and the spin current with BlochMetricBasis, and writes
sigma^z_xy(E - E_F) in e^2/h (Pauli sigma_z) to OUT_PREFIX.dat.

The interband denominator is (E_n - E_m)^2 + eta^2. A small eta removes the
divergent contributions of nearly degenerate pairs between bands, which a
finite k mesh cannot converge; in a gap it changes the plateau by about
(eta / gap)^2, below 1 % here. The occupations are a Gaussian of standard
deviation `broadening`.
"""
import argparse
import os
import sys

import numpy as np
from scipy.special import erf

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "python"))
from lsquant.bases import BlochMetricBasis  # noqa: E402
from lsquant.models.openmx import load_openmx_model  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("scfout")
    ap.add_argument("out_prefix")
    ap.add_argument("--mesh", type=int, default=60, help="N for an N x N x 1 k mesh")
    ap.add_argument("--eta", type=float, default=0.006, help="interband regulator, eV")
    ap.add_argument("--broadening", type=float, default=0.006, help="Gaussian width, eV")
    args = ap.parse_args()

    model = load_openmx_model(args.scfout, dipoles=True)
    area = abs(np.linalg.det(model.uc[:2, :2]))
    basis = BlochMetricBasis(model, [args.mesh, args.mesh, 1], velocity_mode="covariant")
    E, C = basis.physical_eigenvalues, basis.physical_eigenvectors
    Cd = np.conj(np.swapaxes(C, 1, 2))
    J = Cd @ basis.operator_blocks("JXSZ") @ C
    V = Cd @ basis.operator_blocks("VY") @ C
    dE = E[:, :, None] - E[:, None, :]
    pairs = np.abs(dE) > 1e-6
    omega = -2 * np.sum(np.where(pairs, np.imag(J * np.swapaxes(V, 1, 2)) / (dE**2 + args.eta**2), 0), axis=2)
    grid = np.linspace(-0.5, 2.0, 1251)
    occupied = [0.5 * (1 - erf((E - e) / (np.sqrt(2) * args.broadening))) for e in grid]
    sigma = 2 * np.pi * np.array([np.sum(f * omega) for f in occupied]) / (len(E) * area)
    np.savetxt(args.out_prefix + ".dat", np.column_stack([grid, sigma]), fmt="%.6f %.8e",
               header=f"E-E_F (eV)  sigma^z_xy (e^2/h, Pauli sigma_z); OpenMX, exact, covariant, "
                      f"{args.mesh}x{args.mesh}, interband regulator {args.eta * 1e3:g} meV, "
                      f"Gaussian {args.broadening * 1e3:g} meV")


if __name__ == "__main__":
    main()
