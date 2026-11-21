#!/usr/bin/env python3
"""Hall conductivity of the Haldane model written in a non-orthogonal basis.

Reads the rotated Haldane fixture of goldens/stage0/nonorth/haldane/complex: the golden
Haldane model (512 sites, C = +1) re-expressed in a local, non-orthogonal basis with
H' = R^+ H R, V' = R^+ V R and the overlap S' = R^+ R (cond(S') = 5.2). The physics is
unchanged, so every correct method must give the Chern number of the orthogonal model.

The script computes the Kubo-Bastin sigma_xy(E) from the VX-VY moments in three ways:

  lowdin       Chebyshev recursion on S^-1/2 H' S^-1/2, velocities S^-1/2 V' S^-1/2;
  generalised  recursion on S^-1 H' (no square roots), velocities S^-1 V', and the
               covariant trace: bra = S T_m(S^-1 H~) S^-1 VX S^-1 |r>;
  naive        recursion on S^-1 H' with S^-1 V' on both sides and no S on the bra, the
               estimator of the legacy non-orthogonal drivers;

each with the exact overlap appliers, and the Chern number from the gap plateau
C = 2 pi <sigma_xy>_{|E| < 0.5 eV} / A. It also runs the generalised route with the
polynomial (Chebyshev) and iterative (Ischo) appliers of S^-1.

Usage (repository root, after building the core with Python bindings):
    PYTHONPATH=build/lsquant_core:python python3 examples/nonorthogonal/haldane_rotated.py
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np

import lsquant as lsq
from lsquant import _core

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
FIX = os.path.join(ROOT, "goldens", "stage0", "nonorth", "haldane", "complex")
M = 128
WINDOW = (-4.5, 4.5)
A_CELL = 1.5 * math.sqrt(3.0)


def read(name):
    return _core.ReadCsrText(os.path.join(FIX, f"haldane.{name}.CSR"))


def s_bounds():
    """Spectrum of S' from the fixture manifest, widened by 10 %: the Chebyshev applier
    needs 0 < lower < upper, and the Gershgorin bound of S' reaches below 0."""
    with open(os.path.join(FIX, "manifest.json")) as fh:
        man = json.load(fh)
    return _core.SpectralBounds(0.9 * man["lambda_min(S')"], 1.1 * man["lambda_max(S')"])


def chern(mu, n):
    e, sxy = lsq.kubo_bastin(mu, 10.0)
    inside = np.abs(e) < 0.5
    return 2.0 * math.pi * float(np.mean(sxy[inside])) / ((n / 2) * A_CELL)


def states(n):
    return _core.LocalStates(n, list(range(n)))


def lowdin(h, vx, vy, s):
    n = s.dimension
    isq = _core.ExactInverseSqrtOverlap(s)
    lo = lambda m: _core.LowdinOperator(isq, _core.SparseMatrixOperator(m))
    return lsq.run_plan(lsq.kubo_plan(M, left_index=0, right_index=1, word="VX delta VY delta"),
                        lo(h), _core.SpectralWindow(*WINDOW), states(n), operators=[lo(vx), lo(vy)])


def generalised(h, vx, vy, s, inv, covariant=True):
    n = s.dimension
    rh = _core.GeneralisedRescaledOperator(h, s, inv, _core.SpectralWindow(*WINDOW))
    table = [_core.MetricInverseProduct(inv, _core.SparseMatrixOperator(vx)),
             _core.MetricInverseProduct(inv, _core.SparseMatrixOperator(vy)),
             inv, _core.SparseMatrixOperator(s)]
    if covariant:   # bra: S^-1, S^-1 VX, Expand, S;  ket: Expand, S^-1 VY
        code = [_core.Expand(M), _core.ApplyOperator(1), _core.MirrorToBra(), _core.ApplyOperator(2),
                _core.ApplyOperator(0), _core.Expand(M), _core.ApplyOperator(3), _core.Trace([1.0 + 0j])]
    else:           # legacy estimator: S^-1 V on both sides, no S on the bra
        code = [_core.Expand(M), _core.ApplyOperator(1), _core.MirrorToBra(), _core.ApplyOperator(0),
                _core.Expand(M), _core.Trace([1.0 + 0j])]
    mu = _core.MomentEngine(rh, table, states(n), M).run(_core.Plan(code, "VX delta VY delta"))
    _core.ApplyKpmWeights(mu, n, True)
    mu.window = _core.SpectralWindow(*WINDOW)
    mu.system_size = n
    return mu


def main():
    h, vx, vy, s = read("HAM"), read("VX"), read("VY"), read("S")
    n = s.dimension
    golden = _core.ReadChebmom(os.path.join(ROOT, "goldens", "stage0", "haldane",
                                            "NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D"))
    ref = np.asarray(golden.data)
    print(f"rotated Haldane model: {n} sites, M = {M}, exact trace, window {WINDOW}")
    print(f"orthogonal golden:          C = {chern(golden, n):+.4f}")
    exact_inv = _core.ExactInverseOverlap(s)
    runs = [
        ("lowdin, exact S^-1/2", lambda: lowdin(h, vx, vy, s)),
        ("generalised, exact S^-1", lambda: generalised(h, vx, vy, s, exact_inv)),
        ("generalised, Chebyshev S^-1 (P = 60)",
         lambda: generalised(h, vx, vy, s, _core.ChebyshevInverseOverlap(s, s_bounds(), 60))),
        ("generalised, Ischo S^-1 (tol 1e-13)",
         lambda: generalised(h, vx, vy, s, _core.IschoInverseOverlap(s, 1e-13, 200))),
        ("naive (legacy estimator)", lambda: generalised(h, vx, vy, s, exact_inv, covariant=False)),
    ]
    ok = True
    for name, fn in runs:
        t0 = time.perf_counter()
        mu = fn()
        dt = time.perf_counter() - t0
        d = float(np.max(np.abs(np.asarray(mu.data) - ref)) / np.max(np.abs(ref)))
        c = chern(mu, n)
        print(f"{name:38s} C = {c:+.4f}   max|mu - mu_golden|/max|mu_golden| = {d:.1e}   {dt:5.1f} s")
        if name.startswith("naive"):
            ok &= d > 1e-2
        else:
            ok &= d < 1e-8 and abs(c - 1.017) < 0.005
    print("checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
