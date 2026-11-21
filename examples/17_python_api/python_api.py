#!/usr/bin/env python3
"""Drive LSQuant from Python, without a run file on disk, at two levels.

Level 1  ``execute_runfile(dict)``: the same engine as ``python -m lsquant run``,
         fed a Python dictionary instead of a JSON file.
Level 2  The building blocks the run file calls: ``load_wannier_model`` ->
         ``make_basis`` -> ``words.compile`` -> ``run_plan`` -> ``dos_curve`` /
         ``kubo_greenwood``. Use this level to change a step (your own
         operator, state set, or post-processing) inside a Python script.

Both levels compute the DOS and sigma_xx of the graphene p_z model of example
10 (16 x 16 k mesh, Bloch basis, exact trace, 128 moments). The script checks
that the two levels give the same moments and curves.

Usage:
    python examples/17_python_api/python_api.py [out_dir]     # default /tmp/lsq_python_api
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
sys.path.insert(0, os.path.join(_REPO, "python"))

from lsquant import _core, words  # noqa: E402
from lsquant.bases import make_basis  # noqa: E402
from lsquant.models.wannier import load_wannier_model  # noqa: E402
from lsquant.moments import run_plan  # noqa: E402
from lsquant.reconstruct import dos_curve, kubo_greenwood  # noqa: E402
from lsquant.runfile import execute_runfile  # noqa: E402

MODEL_DIR = os.path.join(_REPO, "examples", "10_bases_compare", "graphene")
SUPERCELL = (16, 16, 1)
MOMENTS = 128
WINDOW = (-9.0, 9.0)
BROADENING_MEV = 150.0
DOS_WORD = "Tr[delta(E-H)]"
KUBO_WORD = "Tr[VX delta(E-H) VX delta(E-H)]"


def level1(out_dir):
    """Level 1: a run file given as a Python dictionary."""
    run = {
        "model": {"kind": "wannier", "dir": MODEL_DIR, "seed": "graphene"},
        "basis": {"name": "bloch", "supercell": list(SUPERCELL)},
        "metric": {"policy": "orthogonal"},
        "word": DOS_WORD,
        "words": [DOS_WORD, KUBO_WORD],
        "expansion": {"moments": MOMENTS, "window": list(WINDOW), "trace": "exact"},
        "reconstruction": {"kernel": "jackson", "broadening_meV": BROADENING_MEV,
                           "formulas": ["dos", "greenwood"]},
        "output": {"dir": out_dir, "legacy_layout": False, "write_chebmom": True},
    }
    result = execute_runfile(run)
    moments = {}
    for path in result["moment_files"]:
        key = "dos" if os.path.basename(path).startswith("SpectralOp1") else "kubo"
        moments[key] = _core.ReadChebmom(path)
    return result, moments


def level2():
    """Level 2: model -> basis -> compiled word -> moments -> curve, by hand."""
    model = load_wannier_model(MODEL_DIR, seed="graphene")
    basis = make_basis("bloch", model, SUPERCELL)
    window = _core.SpectralWindow(*WINDOW)
    states = basis.states(trace_mode="exact")           # W block-local states
    operators = {"VX": basis.operator("VX")}
    out = {}
    for key, word in (("dos", DOS_WORD), ("kubo", KUBO_WORD)):
        compiled = words.compile(
            word,
            basis={"name": "bloch", "operators": operators},
            params={"expansion": {"moments": MOMENTS, "window": window,
                                  "states": states.number_of_states}},
        )
        mu = run_plan(compiled.plan, basis.hamiltonian(), window, states,
                      operators=compiled.operators, batch_size=MOMENTS,
                      symmetrise=(key == "kubo"))
        mu.number_of_states = basis.dimension               # exact trace over all N_k W states
        curve = (dos_curve if key == "dos" else kubo_greenwood)(mu, BROADENING_MEV)
        out[key] = (mu, curve)
        print(f"level 2: {word:34s} plan = {len(compiled.plan.code)} instructions, "
              f"moments {mu.dims}, curve {len(curve[0])} points")
    return out


def main():
    out_dir = sys.argv[1] if len(sys.argv) > 1 else "/tmp/lsq_python_api"
    result, mom1 = level1(out_dir)
    print(f"level 1: wrote {len(result['written_files'])} files to {result['output_dir']}")
    lev2 = level2()
    ok = True
    for key in ("dos", "kubo"):
        a = np.asarray(mom1[key].data)
        b = np.asarray(lev2[key][0].data)
        dmu = float(np.max(np.abs(a - b)) / np.max(np.abs(a)))
        func = dos_curve if key == "dos" else kubo_greenwood
        e1, v1 = func(mom1[key], BROADENING_MEV)
        e2, v2 = lev2[key][1]
        dcv = float(np.max(np.abs(v1 - v2)) / np.max(np.abs(v1)))
        good = dmu < 1e-12 and dcv < 1e-12 and np.allclose(e1, e2)
        ok &= good
        print(f"{key:5s} level 1 vs level 2: max rel |d mu| = {dmu:.2e}, max rel |d curve| = {dcv:.2e}  "
              f"{'PASS' if good else 'FAIL'}")
    e, rho = lev2["dos"][1]
    rho = rho / (2 * np.prod(SUPERCELL))
    print(f"graphene DOS per orbital at E = 0: {np.interp(0.0, e, rho):.4f} 1/eV; "
          f"at E = |t| = 2.7 eV (van Hove): {np.interp(2.7, e, rho):.4f} 1/eV")
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
