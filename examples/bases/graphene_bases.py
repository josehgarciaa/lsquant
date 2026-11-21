#!/usr/bin/env python3
"""Graphene density of states in the real-space and Bloch bases.

Loads the graphene Wannier model of golden G9 (goldens/stage0/bases/graphene/model:
nearest-neighbour hopping t = -1 eV, C-C distance 1 Angstrom) and computes its density of
states on an L x L supercell in two bases:

  * Bloch basis with an exact trace: the Hamiltonian is block diagonal in k, so the
    trace over all L^2 k-points needs only W = 2 vectors (BlockLocalStates);
  * real-space basis with a stochastic trace over R random-phase vectors.

Both reconstruct rho(E) per unit cell with the Jackson kernel. The script compares them
and checks two analytic features of graphene:

  * near the Dirac point rho(E) = 2|E| / (sqrt(3) pi t^2) per cell (two valleys, no spin),
    to leading order in E; the script fits rho/|E| = s0 + s2 E^2 and compares s0;
  * the van Hove peaks at |E| = |t| = 1 eV.

Usage (repository root, after building the core with Python bindings):
    PYTHONPATH=build/lsquant_core:python python3 examples/bases/graphene_bases.py \
        [--L 512] [--moments 2048] [--states 4] [--broadening-meV 15] [--out out/bases]
"""
from __future__ import annotations

import argparse
import math
import os
import sys
import time

import numpy as np

import lsquant as lsq
from lsquant import _core
from lsquant.bases import make_basis
from lsquant.models.wannier import load_wannier_model

MODEL = os.path.join(os.path.dirname(__file__), "..", "..", "goldens", "stage0", "bases", "graphene", "model")


def dos_per_cell(basis, states, window, moments, eta):
    t0 = time.perf_counter()
    plan = lsq.dos_plan(moments, int(states.number_of_states), word="1 delta")
    mu = lsq.run_plan(plan, basis.hamiltonian(), window, states, symmetrise=False)
    e, v = lsq.dos_curve(mu, eta)
    # dos_curve returns system_size * rho(E) with rho normalised to 1 per orbital
    return e, v / basis.num_cells, time.perf_counter() - t0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--L", type=int, default=512)
    ap.add_argument("--moments", type=int, default=2048)
    ap.add_argument("--states", type=int, default=4)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--broadening-meV", type=float, default=15.0, dest="eta")
    ap.add_argument("--out", default="out/bases")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)

    model = load_wannier_model(os.path.abspath(MODEL), seed="graphene")
    sc = (a.L, a.L, 1)
    window = _core.SpectralWindow(-3.3, 3.3)          # graphene band: [-3, 3] eV
    bloch = make_basis("bloch", model, sc)
    real = make_basis("real", model, sc)
    eb, rb, tb = dos_per_cell(bloch, bloch.states("exact"), window, a.moments, a.eta)
    er, rr, tr = dos_per_cell(real, real.states("random", num_states=a.states, seed=a.seed),
                              window, a.moments, a.eta)
    np.savetxt(os.path.join(a.out, "graphene_dos.dat"), np.column_stack([eb, rb, rr]),
               header="E (eV), rho per cell: Bloch exact trace, real space random trace")

    n_orb = bloch.dimension
    integral = float(np.sum(0.5 * (rb[1:] + rb[:-1]) * np.diff(eb)))
    # rho(E)/|E| = s0 + s2 E^2 near the Dirac point (s2: trigonal warping); fit away from
    # E = 0, where the broadening rounds the cusp
    dirac = (np.abs(eb) > 0.1) & (np.abs(eb) < 0.4)
    fit = np.linalg.lstsq(np.column_stack([np.ones(dirac.sum()), eb[dirac] ** 2]),
                          rb[dirac] / np.abs(eb[dirac]), rcond=None)[0]
    slope = float(fit[0])
    slope_exact = 2.0 / (math.sqrt(3.0) * math.pi)
    band = (np.abs(eb) > 0.2) & (np.abs(eb) < 2.8)
    dev = float(np.mean(np.abs(rr[band] - rb[band])) / np.mean(rb[band]))
    e_vh = [float(eb[np.argmax(np.where(eb < 0, rb, 0))]), float(eb[np.argmax(np.where(eb > 0, rb, 0))])]

    print(f"graphene {a.L} x {a.L} cells ({n_orb} orbitals), M = {a.moments}, eta = {a.eta} meV")
    print(f"  Bloch basis, exact trace ({bloch.states('exact').number_of_states} vectors): {tb:6.1f} s")
    print(f"  real basis, {a.states} random vectors:                     {tr:6.1f} s")
    print(f"integral of rho per cell (Bloch) = {integral:.4f} (2 orbitals per cell)")
    print(f"Dirac slope s0 from rho/|E| = s0 + s2 E^2 on 0.1 < |E| < 0.4 eV: {slope:.4f} /eV^2 "
          f"(exact 2/(sqrt(3) pi) = {slope_exact:.4f})")
    print(f"van Hove peaks at E = {e_vh[0]:+.3f} and {e_vh[1]:+.3f} eV (exact -+1)")
    print(f"real (random trace) vs Bloch (exact trace), mean relative deviation in the band: {dev:.2e}")
    ok = (abs(integral - 2.0) < 0.01 and abs(slope / slope_exact - 1.0) < 0.01
          and all(abs(abs(x) - 1.0) < 0.02 for x in e_vh) and dev < 0.01)
    print("checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
