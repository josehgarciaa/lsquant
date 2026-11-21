#!/usr/bin/env python3
"""Compare the Lowdin and generalised routes with exact diagonalisation of H c = E S c.

The moments of Tr[delta(E-H)] in a non-orthogonal basis are traces of T_m of
S^(-1/2) H S^(-1/2) (Lowdin route) or of S^(-1) H (generalised route). Both
matrices have the eigenvalues E_i of the generalised problem H c = E S c, so
both routes must give the exact moments

    mu_m = (2 - delta_m0) / N * sum_i T_m(x_i),   x_i = (E_i - b) / a,

with (a, b) the half-width and centre of the spectral window. The script
solves the generalised eigenproblem with SciPy, checks it against the analytic
chain spectrum 2t cos(2 pi k/N), and compares the moments of every run.

Usage:
    python check_nonorth.py <out_root>   # holds nonorth_{orthogonal,lowdin,generalised,lowdin_chebyshev}/
"""
import glob
import os
import sys

import numpy as np
import scipy.linalg as la

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_REPO, "python"))

from lsquant import _core  # noqa: E402
from lsquant.models.nonorthogonal import read_csr_text  # noqa: E402

OPS = os.path.join(_HERE, "operators")
WINDOW = (-2.2, 2.2)
T = -1.0
RUNS = (("orthogonal", 1e-10), ("lowdin", 1e-10), ("generalised", 1e-10),
        ("lowdin_chebyshev", 1e-8))


def exact_moments(num_moments):
    h = read_csr_text(os.path.join(OPS, "nonorth_chain.HAM.CSR")).toarray()
    s = read_csr_text(os.path.join(OPS, "nonorth_chain.S.CSR")).toarray()
    evals = la.eigh(h, s, eigvals_only=True)
    n = len(evals)
    analytic = np.sort(2 * T * np.cos(2 * np.pi * np.arange(n) / n))
    a = 0.5 * (WINDOW[1] - WINDOW[0])
    b = 0.5 * (WINDOW[1] + WINDOW[0])
    theta = np.arccos(np.clip((evals - b) / a, -1.0, 1.0))
    mu = np.cos(np.outer(np.arange(num_moments), theta)).sum(axis=1) / n
    mu[1:] *= 2.0
    s_eigs = la.eigvalsh(s)
    return mu, evals, analytic, s_eigs


def run_moments(run_dir):
    path = sorted(glob.glob(os.path.join(run_dir, "*.chebmom")))[0]
    return np.asarray(_core.ReadChebmom(path).data).real.ravel()


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "out"
    ok = True
    mu_ref = None
    for name, tol in RUNS:
        mu = run_moments(os.path.join(root, f"nonorth_{name}"))
        if mu_ref is None:
            mu_ref, evals, analytic, s_eigs = exact_moments(len(mu))
            dev = float(np.max(np.abs(evals - analytic)))
            print(f"overlap S: N = {len(evals)}, spectrum [{s_eigs.min():.4f}, {s_eigs.max():.4f}], "
                  f"cond(S) = {s_eigs.max() / s_eigs.min():.3f}")
            print(f"generalised eigenvalues vs 2t cos(2 pi k/N): max deviation {dev:.2e}")
            ok &= dev < 1e-12
        err = float(np.max(np.abs(mu - mu_ref)))
        good = err < tol
        ok &= good
        print(f"{name:18s} max |mu_m - mu_m(exact)| = {err:.3e}  (tol {tol:g})  {'PASS' if good else 'FAIL'}")
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
