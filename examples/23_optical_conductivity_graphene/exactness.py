#!/usr/bin/env python3
"""Exact-trace check: LSQuant's moments and reconstruction against the Kubo-Bastin formula.

Problem: the tests of benchmark.py use a random trace, so they agree with the
exact response only within the noise. This script removes the noise: on a
small supercell it computes the moments with an exact trace (every basis
state) and compares sigma(mu, omega) from `optical_conductivity` with the
finite-frequency Kubo-Bastin formula evaluated independently on the
eigenpairs of the same supercell, with the same Jackson-damped expansions of
delta(E - H) and G+-(E):

    sigma_ab = i/(Omega omega) [S_ab(omega) - S_ab(0)],
    S_ab(mu, omega) = int^mu dE Tr[v_a delta_M(E-H) v_b G+_M(E+omega) + v_a G-_M(E-omega) v_b delta_M(E-H)],

with G+ expanded in exp(-i n arccos x) and no complex conjugate anywhere
(derivation.py, step 2). The two must agree to rounding. The script also
checks the normalisation of the moment file, mu_mn = 4 h_m h_n Tr[v_a T_m v_b T_n]/N,
sigma_xx = sigma_yy (the n x n supercell keeps the threefold axis), and, for
the Haldane model, that sigma_xy(omega -> 0) is the static Kubo-Bastin value.

Usage (from the repository root):
    python examples/23_optical_conductivity_graphene/exactness.py [out_dir]
Runs LSQuant three times (graphene, hBN, Haldane; 6 x 6 cells, 72 orbitals,
128 moments; seconds) and prints PASS or FAIL.
"""
import glob
import json
import os
import subprocess
import sys

import numpy as np

from lsquant import _core
from lsquant.reconstruct import optical_conductivity, shc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_honeycomb import GRAPHENE, HBN, lattice, operators, write_csr   # noqa: E402

CELLS, MOMENTS, HALF_WIDTH = 6, 128, 20.0
ENERGY_WINDOW = (-10.0, 10.0)
PHOTON = np.array([1e-4, 0.5, 1.0, 2.0, 3.0, 5.4, 7.25, 8.6, 9.5])
WORDS = {"xx": "Tr[VX delta(E-H) VX dG-(E)/dE]", "yy": "Tr[VY delta(E-H) VY dG-(E)/dE]",
         "xy": "Tr[VX delta(E-H) VY dG-(E)/dE]"}
MODELS = {"graphene": GRAPHENE, "hbn": HBN,
          "haldane": dict(GRAPHENE, t2=-0.4, phi=np.pi / 2)}


def run_moments(name, pars, out):
    ham, vx, vy = operators(CELLS, **pars)
    ops = os.path.join(out, f"ops_{name}")
    os.makedirs(ops, exist_ok=True)
    for op, m in (("HAM", ham), ("VX", vx), ("VY", vy)):
        write_csr(os.path.join(ops, f"{name}.{op}.CSR"), m)
    res = os.path.join(out, f"moments_{name}")
    rf = {"model": {"kind": "csr", "dir": ops, "label": name, "operators": ["VX", "VY"]},
          "basis": {"name": "real"}, "metric": {"policy": "orthogonal"},
          "word": WORDS["xx"], "words": list(WORDS.values()),
          "expansion": {"moments": MOMENTS, "window": [-HALF_WIDTH, HALF_WIDTH], "trace": "exact"},
          "reconstruction": {"kernel": "jackson", "broadening_meV": 100.0, "formulas": ["bastin"]},
          "output": {"dir": res, "legacy_layout": False, "write_chebmom": True}}
    path = os.path.join(out, f"run_{name}.json")
    with open(path, "w") as fh:
        json.dump(rf, fh, indent=1)
    # 72 orbitals: a few threads are fastest; many threads only contend.
    env = dict(os.environ)
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        env.setdefault(var, "2")
    subprocess.run([sys.executable, "-m", "lsquant", "run", path], check=True, capture_output=True, env=env)
    moments = {}
    for key, tag in (("xx", "VX-VX"), ("yy", "VY-VY"), ("xy", "VX-VY")):
        moments[key] = _core.ReadChebmom(glob.glob(os.path.join(res, f"NonEqOp{tag}*.chebmom"))[0])
    return ham.toarray(), vx.toarray(), vy.toarray(), moments


def kubo_bastin_eigen(ham, va, vb, grid, area):
    """sigma_ab(mu, omega) in e^2/h from the eigenpairs, Jackson-damped delta and G (see module doc)."""
    eps, vec = np.linalg.eigh(ham)
    a, b = vec.conj().T @ va @ vec, vec.conj().T @ vb @ vec
    pair = a.T * b                                   # Tr[v_a P_i v_b P_j] = v^a_ji v^b_ij
    n = np.arange(MOMENTS)
    g = ((MOMENTS - n + 1) * np.cos(np.pi * n / (MOMENTS + 1))
         + np.sin(np.pi * n / (MOMENTS + 1)) / np.tan(np.pi / (MOMENTS + 1))) / (MOMENTS + 1)
    coef = (g * np.where(n == 0, 1.0, 2.0))[:, None] * np.cos(n[:, None] * np.arccos(eps / HALF_WIDTH))
    x = grid / HALF_WIDTH
    delta = lambda u: coef.T @ np.cos(n[:, None] * np.arccos(u)) / (np.pi * np.sqrt(1 - u**2))
    green = lambda u, s: -1j * s * coef.T @ np.exp(-1j * s * n[:, None] * np.arccos(u)) / np.sqrt(1 - u**2)

    def big_s(w):
        s = (np.einsum("ix,ij,jx->x", delta(x), pair, green(x + w, 1))
             + np.einsum("ix,ij,jx->x", green(x - w, -1), pair, delta(x)))
        out = np.zeros_like(s)
        out[1:] = np.cumsum(0.5 * (s[1:] + s[:-1]) * np.diff(x))
        return out

    base = big_s(0.0)
    omega_area = CELLS * CELLS * area
    return np.array([2 * np.pi * 1j * (big_s(w / HALF_WIDTH) - base) / (omega_area * HALF_WIDTH * w)
                     for w in PHOTON])


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "out/optical_exactness"
    os.makedirs(out, exist_ok=True)
    failed = 0

    def check(ok, text):
        nonlocal failed
        failed += not ok
        print(("PASS  " if ok else "FAIL  ") + text)

    for name, pars in MODELS.items():
        ham, vx, vy, moments = run_moments(name, pars, out)
        area = lattice(pars["a"])[2]
        dim = ham.shape[0]
        # Normalisation of the file: mu_mn = 4 h_m h_n Tr[v_a T_m v_b T_n] / N.
        eps, vec = np.linalg.eigh(ham)
        t = np.cos(np.arange(MOMENTS)[:, None] * np.arccos(eps / HALF_WIDTH)[None, :])
        a, b = vec.conj().T @ vx @ vec, vec.conj().T @ vy @ vec
        half = np.where(np.arange(MOMENTS) == 0, 0.5, 1.0)
        scale = None
        for key, (va, vb) in (("xx", (a, a)), ("xy", (a, b))):
            trace = np.einsum("ji,ij,mi,nj->mn", va, vb, t, t)      # Tr[v_a T_m v_b T_n]
            want = 4 * np.outer(half, half) * trace / dim
            scale = scale or np.abs(want).max()                     # the xx moments set the scale
            err = np.abs(np.asarray(moments[key].data) - want).max() / scale
            check(err < 1e-10, f"{name}: moment file {key} = 4 h_m h_n Tr[v_a T_m v_b T_n]/N "
                               f"(error {err:.1e} of the largest xx moment)")
        diff_xy = np.abs(np.asarray(moments["xx"].data) - np.asarray(moments["yy"].data)).max()
        check(diff_xy < 1e-9 * np.abs(np.asarray(moments["xx"].data)).max(),
              f"{name}: moments of v_x v_x and v_y v_y equal (max difference {diff_xy:.1e})")
        scale = None
        for key, (va, vb) in (("xx", (vx, vx)), ("xy", (vx, vy))):
            energies, sigma = optical_conductivity(moments[key], PHOTON, energy_window=ENERGY_WINDOW,
                                                   num_cells=CELLS * CELLS, cell_area=area)
            ref = kubo_bastin_eigen(ham, va, vb, energies, area)
            scale = scale or np.abs(ref).max()                      # sigma_xx sets the scale
            err = np.abs(sigma - ref).max() / scale
            check(err < 1e-9, f"{name}: sigma_{key}(mu, omega) = Kubo-Bastin on the eigenpairs "
                              f"(error {err:.1e} of max |sigma_xx|; max |sigma_{key}| "
                              f"{np.abs(ref).max():.3g} e^2/h)")
        if name == "haldane":
            energies, sigma = optical_conductivity(moments["xy"], PHOTON[:1], energy_window=ENERGY_WINDOW,
                                                   num_cells=CELLS * CELLS, cell_area=area)
            grid, static = shc(moments["xy"], 1.0, CELLS * CELLS, area)
            opt0, st0 = np.interp(0.0, energies, sigma[0].real), np.interp(0.0, grid, static)
            check(abs(opt0 - st0) < 0.01 * abs(st0),
                  f"haldane: sigma_xy(0, 1e-4 eV) = {opt0:.4f}, static Kubo-Bastin {st0:.4f} e^2/h")
    print("\nPASS" if not failed else f"\n{failed} FAIL")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
