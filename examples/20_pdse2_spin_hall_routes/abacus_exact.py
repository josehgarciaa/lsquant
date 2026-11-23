"""Exact spin Hall conductivity of an ABACUS model on a dense k mesh.

Usage:
    OMP_NUM_THREADS=1 python abacus_exact.py RUN_DIR OUT_PREFIX [--mesh 60] [--eta 0.006]
        [--broadening 0.006] [--velocity covariant] [--chunk 24] [--workers 16]

RUN_DIR holds the INPUT, STRU and OUT.<suffix>/ of one ABACUS run with
``out_mat_hs2 1`` and ``out_mat_r 1``. The script reads H(R), S(R) and r(R)
with the LSQuant ABACUS reader, builds the velocity and the spin current with
BlochMetricBasis, and writes sigma^z_xy(E - E_F) in e^2/h (Pauli sigma_z) to
OUT_PREFIX_<velocity>.dat, on -1.5 to 2.6 eV in 2 meV steps. E_F is the
last E_Fermi of running_scf.log.

Same formula as openmx_exact.py: interband denominator (E_n - E_m)^2 + eta^2
and Gaussian occupations of standard deviation `broadening`. The k mesh is
processed in chunks of `chunk` points by `workers` processes. Run with
OMP_NUM_THREADS=1: for matrices of a few hundred rows, single-threaded
processes are much faster than threaded linear algebra on a shared node.
"""
import argparse
import multiprocessing as mp
import os
import time

import numpy as np
from scipy.special import erf

from lsquant.bases import BlochMetricBasis

from abacus_bands import load


_STATE = {}


def _chunk_sum(c0):
    """sigma accumulators of one chunk of the k mesh, for every velocity mode."""
    model, args, grid, full_mesh = (_STATE[k] for k in ("model", "args", "grid", "mesh"))
    chunk = full_mesh[c0:c0 + args.chunk]
    model.k_mesh = lambda supercell: chunk
    out = {}
    for v in args.velocity:
        basis = BlochMetricBasis(model, [args.mesh, args.mesh, 1], velocity_mode=v)
        E, C = basis.physical_eigenvalues, basis.physical_eigenvectors
        Cd = np.conj(np.swapaxes(C, 1, 2))
        J = Cd @ basis.operator_blocks("JXSZ") @ C
        V = Cd @ basis.operator_blocks("VY") @ C
        dE = E[:, :, None] - E[:, None, :]
        pairs = np.abs(dE) > 1e-6
        omega = -2 * np.sum(np.where(pairs, np.imag(J * np.swapaxes(V, 1, 2)) / (dE**2 + args.eta**2), 0), axis=2)
        out[v] = np.array([np.sum(0.5 * (1 - erf((E - e) / (np.sqrt(2) * args.broadening))) * omega)
                           for e in grid])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("out_prefix")
    ap.add_argument("--mesh", type=int, default=60, help="N for an N x N x 1 k mesh")
    ap.add_argument("--eta", type=float, default=0.006, help="interband regulator, eV")
    ap.add_argument("--broadening", type=float, default=0.006, help="Gaussian width, eV")
    ap.add_argument("--velocity", nargs="+", default=["covariant"])
    ap.add_argument("--chunk", type=int, default=24)
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()

    model, ef, _ = load(args.run_dir)
    area = abs(np.linalg.det(model.uc[:2, :2]))
    full_mesh = model.k_mesh([args.mesh, args.mesh, 1])
    grid = np.linspace(-1.5, 2.6, 2051)
    _STATE.update(model=model, args=args, grid=grid, mesh=full_mesh)
    acc = {v: np.zeros(len(grid)) for v in args.velocity}
    t0 = time.time()
    starts = list(range(0, len(full_mesh), args.chunk))
    with mp.get_context("fork").Pool(args.workers) as pool:
        for n, part in enumerate(pool.imap_unordered(_chunk_sum, starts), 1):
            for v in args.velocity:
                acc[v] += part[v]
            if n % max(1, len(starts) // 20) == 0 or n == len(starts):
                print(f"{n}/{len(starts)} chunks, {time.time() - t0:.0f} s", flush=True)
    for v in args.velocity:
        sigma = 2 * np.pi * acc[v] / (len(full_mesh) * area)
        np.savetxt(f"{args.out_prefix}_{v}.dat", np.column_stack([grid, sigma]), fmt="%.6f %.8e",
                   header=f"E-E_F (eV)  sigma^z_xy (e^2/h, Pauli sigma_z); ABACUS {os.path.basename(os.path.normpath(args.run_dir))}, "
                          f"exact, {v}, {args.mesh}x{args.mesh}, interband regulator {args.eta * 1e3:g} meV, "
                          f"Gaussian {args.broadening * 1e3:g} meV; E_F = {ef} eV")


if __name__ == "__main__":
    main()
