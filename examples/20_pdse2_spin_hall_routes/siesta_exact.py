"""Exact spin Hall conductivity of a SIESTA model on a dense k mesh.

Usage:
    python siesta_exact.py SIESTA_DIR OUT_PREFIX --vbm VBM [--label pdse2] [--mesh 60]
        [--eta 0.006] [--broadening 0.006] [--velocity covariant] [--chunk 200]

Reads H, S (.TSHS) and the radial orbitals (.ion.xml) with the LSQuant SIESTA
reader, builds the velocity (covariant: dipoles by quadrature of the radial
orbitals) and the spin current with BlochMetricBasis, and writes
sigma^z_xy(E - E_VBM) in e^2/h (Pauli sigma_z) to OUT_PREFIX_<velocity>.dat.
VBM is the valence-band maximum in the energies of the .TSHS (E_F = 0), from
band_edges.py.

Same formula as openmx_exact.py: interband denominator (E_n - E_m)^2 + eta^2
and Gaussian occupations of standard deviation `broadening`. The k mesh is
processed in chunks to bound the memory.
"""
import argparse
import time

import numpy as np
from scipy.special import erf

from lsquant.bases import BlochMetricBasis
from lsquant.models.siesta import load_siesta_model


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("siesta_dir")
    ap.add_argument("out_prefix")
    ap.add_argument("--vbm", type=float, required=True, help="VBM in eV, in the energies of the .TSHS")
    ap.add_argument("--label", default="pdse2")
    ap.add_argument("--mesh", type=int, default=60)
    ap.add_argument("--eta", type=float, default=0.006)
    ap.add_argument("--broadening", type=float, default=0.006)
    ap.add_argument("--velocity", nargs="+", default=["covariant"])
    ap.add_argument("--chunk", type=int, default=200)
    args = ap.parse_args()

    model = load_siesta_model(args.siesta_dir, args.label, dipoles="covariant" in args.velocity)
    area = abs(np.linalg.det(model.uc[:2, :2]))
    full_mesh = model.k_mesh([args.mesh, args.mesh, 1])
    grid = np.linspace(-0.5, 3.0, 1751)  # E - E_VBM
    acc = {v: np.zeros(len(grid)) for v in args.velocity}
    t0 = time.time()
    for c0 in range(0, len(full_mesh), args.chunk):
        chunk = full_mesh[c0:c0 + args.chunk]
        model.k_mesh = lambda supercell, chunk=chunk: chunk
        for v in args.velocity:
            basis = BlochMetricBasis(model, [args.mesh, args.mesh, 1], velocity_mode=v)
            E, C = basis.physical_eigenvalues - args.vbm, basis.physical_eigenvectors
            Cd = np.conj(np.swapaxes(C, 1, 2))
            J = Cd @ basis.operator_blocks("JXSZ") @ C
            V = Cd @ basis.operator_blocks("VY") @ C
            dE = E[:, :, None] - E[:, None, :]
            pairs = np.abs(dE) > 1e-6
            omega = -2 * np.sum(np.where(pairs, np.imag(J * np.swapaxes(V, 1, 2)) / (dE**2 + args.eta**2), 0),
                                axis=2)
            for i, e in enumerate(grid):
                acc[v][i] += np.sum(0.5 * (1 - erf((E - e) / (np.sqrt(2) * args.broadening))) * omega)
        print(f"{c0 + len(chunk)}/{len(full_mesh)} k, {time.time() - t0:.0f} s", flush=True)
    for v in args.velocity:
        sigma = 2 * np.pi * acc[v] / (len(full_mesh) * area)
        np.savetxt(f"{args.out_prefix}_{v}.dat", np.column_stack([grid, sigma]), fmt="%.6f %.8e",
                   header=f"E-E_VBM (eV)  sigma^z_xy (e^2/h, Pauli sigma_z); SIESTA, exact, {v}, "
                          f"{args.mesh}x{args.mesh}, interband regulator {args.eta * 1e3:g} meV, "
                          f"Gaussian {args.broadening * 1e3:g} meV; VBM = {args.vbm} eV (E_F = 0)")


if __name__ == "__main__":
    main()
