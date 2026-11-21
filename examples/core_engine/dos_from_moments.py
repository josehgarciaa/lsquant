#!/usr/bin/env python3
"""Density of states from a rank-1 .chebmom file (Jackson kernel), with a check
against the analytic DOS of the 1D chain.

The stage-1 "dos" plan stores mu_0 = Tr[T_0(H~)] / N and mu_m = 2 Tr[T_m(H~)] / N for
m >= 1 (the factor 2 of the Chebyshev expansion is already included). The DOS per site
in energy units is then

    rho(E) = sum_m g_m mu_m T_m(x) / (pi a sqrt(1 - x^2)),   x = (E - c) / a,

with a, c the half-width and centre of the spectral window and g_m the Jackson kernel.

usage: dos_from_moments.py file.chebmom [--out dos.dat] [--points 401]
"""
import argparse
import math

import numpy as np


def read_chebmom(path):
    """Return (window, moments) from a version-1 .chebmom file of rank 1."""
    header, data = {}, []
    with open(path) as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        tok = lines[i].split()
        i += 1
        if not tok:
            continue
        if tok[0] == "data":
            break
        header[tok[0]] = tok[1:]
    for line in lines[i:]:
        t = line.split()
        if t:
            data.append(complex(float(t[0]), float(t[1])))
    if header.get("rank", ["1"])[0] != "1":
        raise SystemExit("expected a rank-1 .chebmom file")
    lo, hi = (float(v) for v in header["window"][:2])
    return (lo, hi), np.array(data)


def jackson(M):
    m = np.arange(M)
    q = math.pi / (M + 1)
    return ((M - m + 1) * np.cos(q * m) + np.sin(q * m) / math.tan(q)) / (M + 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("chebmom")
    ap.add_argument("--out", default="dos.dat")
    ap.add_argument("--points", type=int, default=401)
    args = ap.parse_args()

    (lo, hi), mu = read_chebmom(args.chebmom)
    a, c = 0.5 * (hi - lo), 0.5 * (hi + lo)
    M = len(mu)
    g = jackson(M)
    x = np.linspace(-0.999, 0.999, args.points)
    T = np.cos(np.outer(np.arange(M), np.arccos(x)))     # T_m(x)
    rho = (g * mu.real) @ T / (math.pi * a * np.sqrt(1.0 - x * x))
    E = c + a * x
    np.savetxt(args.out, np.column_stack([E, rho]), header="E  rho(E) per site (Jackson, M=%d)" % M)

    # analytic DOS of the chain with |t| = 1, away from the van Hove edges
    inner = np.abs(E) < 1.5
    exact = 1.0 / (math.pi * np.sqrt(4.0 - E[inner] ** 2))
    err = np.max(np.abs(rho[inner] - exact) / exact)
    trapezoid = getattr(np, "trapezoid", None) or np.trapz
    norm = trapezoid(rho, E)
    print("wrote %s: %d points, M = %d" % (args.out, len(E), M))
    print("integral of rho over the band = %.4f (1 expected, minus the weight beyond |x| = 0.999)" % norm)
    print("max relative deviation from 1/(pi sqrt(4 - E^2)) for |E| < 1.5: %.2e" % err)
    if err > 0.05:
        raise SystemExit("deviation above 5 %: check the moments")


if __name__ == "__main__":
    main()
