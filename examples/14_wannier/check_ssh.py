#!/usr/bin/env python3
"""Check the Wannier SSH run against the analytic SSH bands.

The SSH chain has two bands E(k) = +/- |t1 + t2 exp(i k a)|. On the
256-point k mesh of the run (k = 2 pi m / (256 a)) the exact DOS moments are

    mu_m = (2 - delta_m0) / (2 N_k) * sum_{k, +/-} T_m((E(k) - b) / a_w),

with (a_w, b) the half-width and centre of the spectral window stored in the
moment file. The script compares them with the moments LSQuant computed from
ssh_hr.dat, and checks that the reconstructed DOS vanishes inside the gap
|E| < |t1 - t2| = 0.4 eV.

Usage:
    python check_ssh.py <run_out_dir>
"""
import glob
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(_HERE)), "python"))
from lsquant import _core  # noqa: E402

T1, T2, NK = -1.0, -0.6, 256


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "out/ssh_dos"
    path = sorted(glob.glob(os.path.join(out, "SpectralOp1*.chebmom")))[0]
    mu = _core.ReadChebmom(path)
    data = np.asarray(mu.data).real.ravel()
    lo, hi = float(mu.window.lower), float(mu.window.upper)
    a_w, b = 0.5 * (hi - lo), 0.5 * (hi + lo)
    k = 2.0 * np.pi * np.arange(NK) / NK
    e = np.abs(T1 + T2 * np.exp(1j * k))
    energies = np.concatenate([-e, e])
    theta = np.arccos(np.clip((energies - b) / a_w, -1, 1))
    ref = np.cos(np.outer(np.arange(len(data)), theta)).mean(axis=1)
    ref[1:] *= 2.0
    err = float(np.max(np.abs(data - ref)))
    print(f"moment file: {os.path.basename(path)}")
    print(f"window [{lo:.6f}, {hi:.6f}] eV; analytic band edges +/-{e.min():.3f}, +/-{e.max():.3f} eV")
    print(f"max |mu_m - mu_m(analytic)| over {len(data)} moments = {err:.3e}")
    curve = sorted(glob.glob(os.path.join(out, "meanSpectralOp1*JACKSON.dat")))[0]
    c = np.loadtxt(curve)
    dos = c[:, 1] / (2 * NK)          # extensive -> per orbital
    gap = np.abs(c[:, 0]) < 0.25
    band = (np.abs(c[:, 0]) > 0.6) & (np.abs(c[:, 0]) < 1.4)
    print(f"DOS per orbital: max inside |E|<0.25 eV = {dos[gap].max():.2e}; "
          f"mean inside 0.6<|E|<1.4 eV = {dos[band].mean():.3f} 1/eV")
    ok = err < 1e-10 and dos[gap].max() < 1e-3 * dos[band].mean()
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
