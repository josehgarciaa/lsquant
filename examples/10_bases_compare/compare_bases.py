#!/usr/bin/env python3
"""Compare the DOS and sigma_xx of graphene computed in the real, bloch and eigen bases.

All three runs use the same 16 x 16 supercell (k mesh), the same explicit
spectral window and an exact trace, so they evaluate the same traces in three
different representations of the same Hamiltonian:

    real   the 512 x 512 real-space supercell matrix, 512 local trace states
    bloch  the 256 blocks H(k), 2 block-local trace states
    eigen  the diagonal eigenvalues of each H(k), 2 block-local trace states

The moments must agree to round-off. As an independent check, the DOS moments
are also computed from the analytic graphene bands
E(k) = +/- |t| |1 + exp(-i k.a1) + exp(-i k.a2)| on the same mesh.

Usage:
    python compare_bases.py <out_root>    # holds bases_real/, bases_bloch/, bases_eigen/
"""
import glob
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(_HERE)), "python"))
from lsquant import _core  # noqa: E402

BASES = ("real", "bloch", "eigen")
T, N1, N2 = -2.7, 16, 16


def load(root, basis, prefix):
    path = sorted(glob.glob(os.path.join(root, f"bases_{basis}", f"{prefix}*.chebmom")))[0]
    return _core.ReadChebmom(path)


def curve(root, basis, prefix):
    path = sorted(glob.glob(os.path.join(root, f"bases_{basis}", f"{prefix}*JACKSON.dat")))[0]
    return np.loadtxt(path)


def analytic_dos_moments(window, num_moments):
    m1, m2 = np.meshgrid(np.arange(N1) / N1, np.arange(N2) / N2, indexing="ij")
    f = 1 + np.exp(-2j * np.pi * m1) + np.exp(-2j * np.pi * m2)
    e = np.concatenate([np.abs(T * f).ravel(), -np.abs(T * f).ravel()])
    a, b = 0.5 * (window.upper - window.lower), 0.5 * (window.upper + window.lower)
    theta = np.arccos(np.clip((e - b) / a, -1, 1))
    mu = np.cos(np.outer(np.arange(num_moments), theta)).mean(axis=1)
    mu[1:] *= 2
    return mu


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "out"
    ok = True
    for label, prefix, curve_prefix in (("DOS", "SpectralOp1", "meanSpectralOp1"),
                                        ("sigma_xx", "NonEqOpVX-VX", "KuboGreenWood_NonEqOpVX-VX")):
        mus = {b: np.asarray(load(root, b, prefix).data) for b in BASES}
        curves = {b: curve(root, b, curve_prefix) for b in BASES}
        scale = np.max(np.abs(mus["real"]))
        for b in ("bloch", "eigen"):
            dmu = float(np.max(np.abs(mus[b] - mus["real"]))) / scale
            dcv = float(np.max(np.abs(curves[b][:, 1] - curves["real"][:, 1]))
                        / np.max(np.abs(curves["real"][:, 1])))
            good = dmu < 1e-10 and dcv < 1e-10
            ok &= good
            print(f"{label:8s} {b:5s} vs real: max rel |d mu| = {dmu:.2e}, "
                  f"max rel |d curve| = {dcv:.2e}  {'PASS' if good else 'FAIL'}")
    mu_r = load(root, "real", "SpectralOp1")
    ref = analytic_dos_moments(mu_r.window, mu_r.dims[0])
    err = float(np.max(np.abs(np.asarray(mu_r.data).real.ravel() - ref)))
    good = err < 1e-10
    ok &= good
    print(f"DOS real vs analytic bands: max |d mu| = {err:.2e}  {'PASS' if good else 'FAIL'}")
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
