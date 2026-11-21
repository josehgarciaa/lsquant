#!/usr/bin/env python3
"""Scalar and magnetic Anderson disorder in the Bi bilayer: Kramers pairs.

Loads the SIESTA Bi(111) bilayer of examples/07 (spin-orbit coupling, non-orthogonal basis),
builds a 3 x 3 x 1 supercell and adds Anderson disorder of width W = 0.5 eV with each
basis.disorder.kind, V_ij = (w_i + w_j) S_ij / 2:

  scalar       one value per orbital on both spin components (time-reversal symmetric);
  magnetic     +h on spin up, -h on spin down (a random Zeeman field along z);
  independent  one value per spinor coefficient (scalar and magnetic parts together).

It solves H c = E S c densely and prints the largest splitting of the eigenvalue pairs
(E_2k, E_2k+1). Time-reversal symmetry with spin-orbit coupling makes every level doubly
degenerate (Kramers). The stored SIESTA Hamiltonian is rounded, so the clean pairs are
degenerate only to that rounding; the scalar kind keeps this pairing, and the magnetic and
independent kinds split the pairs by a fraction of W.

usage (repository root): PYTHONPATH=build/golden/lsquant_core:python python3 examples/disorder_kinds/kramers.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.linalg import eigh

from lsquant.disorder_bases import spinor_disorder_values
from lsquant.models.disorder import anderson_potential, supercell_csr
from lsquant.models.siesta import load_siesta_model

BI = os.path.join(os.path.dirname(__file__), "..", "07_spin_hall_bi_bilayer", "siesta")
SUPERCELL = (3, 3, 1)
W, SEED = 0.5, 7


def pair_splitting(h, s, window=5.0):
    """Largest E_2k+1 - E_2k over the whole spectrum and over pairs with |E| < window."""
    e = eigh(h, s, eigvals_only=True)
    d = e[1::2] - e[0::2]
    near = np.abs(0.5 * (e[1::2] + e[0::2])) < window
    return e, float(np.max(d)), float(np.max(d[near]))


def main():
    model = load_siesta_model(os.path.abspath(BI), "bi")
    h = supercell_csr(model, SUPERCELL, model.hoppings)
    s = supercell_csr(model, SUPERCELL, model.overlaps)
    h = 0.5 * (h + h.conj().T)
    s = 0.5 * (s + s.conj().T)
    dim = h.shape[0]
    hd, sd = h.toarray(), s.toarray()
    e0, split0, near0 = pair_splitting(hd, sd)
    print(f"Bi bilayer, supercell {SUPERCELL}: {dim} spinor orbitals, W = {W} eV, seed {SEED}")
    print("               largest Kramers-pair splitting: all levels | levels with |E| < 5 eV")
    print(f"  clean        {split0:.2e} eV | {near0:.2e} eV")
    ok = True
    results = {}
    for kind in ("scalar", "magnetic", "independent"):
        w = spinor_disorder_values(dim, W, SEED, kind)
        v = anderson_potential(s, w).toarray()
        e, split, near = pair_splitting(hd + v, sd)
        results[kind] = near
        print(f"  {kind:12s} {split:.2e} eV | {near:.2e} eV")
    # the scalar kind keeps the clean pairing; the kinds with a magnetic part split it
    ok &= results["scalar"] < 2 * near0 + 1e-9
    ok &= results["magnetic"] > 100 * near0 and results["independent"] > 100 * near0
    print("checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
