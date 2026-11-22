#!/usr/bin/env python3
"""Developer golden of the OpenMX, Wannier90 and PAOFLOW connectors (PdSe2 monolayer).

Runs the spin Hall run files of an OpenMX model (covariant velocity), a
Wannier90 model (covariant velocity, exact spin) and a PAOFLOW model through
execute_runfile, and compares the KPM curve and the exact projected oracle
with stored reference values. The 12 x 12 mesh and 512 moments keep the run
short: the values are a regression reference, not converged physics.

The models are too large for the repository. Give their directory:
  PYTHONPATH=python:<core> python tests/dev/qe_connectors_check.py DATA_DIR
DATA_DIR holds openmx/pdse2.scfout, w90/ (pdse2_proj_hr.dat, _r.dat, .win,
_centres.xyz, .spn, .eig, _u.mat, _u_dis.mat) and pao/pdse2_paoflow.npz
(write_paoflow_model).
Without DATA_DIR, or with LSQUANT_QE_GOLDEN_DIR unset, the check is skipped.
--write stores new reference values instead of comparing.

The reference values were checked against SIESTA and against PAOFLOW's own
spin Hall curve for the same material (developer validation, Nov 2026).
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
REFERENCE = HERE.parent.parent / "goldens" / "qe_connectors" / "pdse2_reference.json"
ENERGIES = [0.0, 1.21, 1.40]   # E - E_F (eV): trivial gap, near the QE topological gap, metal
COMMON = {
    "basis": {"name": "eigen", "supercell": [12, 12, 1]},
    "metric": {"policy": "orthogonal"},
    "projection": {"kind": "shift", "window": [-3, 3]},
    "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
    "expansion": {"moments": 512, "window": "eigen", "trace": "exact", "seed": 12345},
    "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 6.0,
                       "energy_grid": [-0.5, 2.0, 1251]},
    "output": {"legacy_layout": False, "write_chebmom": True},
}


def models(data):
    return {
        "openmx": {"kind": "openmx", "file": str(data / "openmx" / "pdse2.scfout"),
                   "params": {"velocity_mode": "covariant", "oracle": True}},
        "wannier": {"kind": "wannier", "dir": str(data / "w90"), "seed": "pdse2_proj",
                    "params": {"velocity_mode": "covariant", "spin": "exact",
                               "fermi_energy": -1.3162, "oracle": True}},
        "paoflow": {"kind": "paoflow", "file": str(data / "pao" / "pdse2_paoflow.npz"),
                    "params": {"velocity_mode": "berry_connection", "oracle": True}},
    }


def run(name, model, out):
    from lsquant.runfile import execute_runfile
    rf = dict(COMMON, model=model)
    if model["kind"] == "openmx":
        rf["metric"] = {"policy": "nonorthogonal", "route": "lowdin", "applier": "exact"}
    rf["output"] = dict(COMMON["output"], dir=str(out / name))
    execute_runfile(rf)
    e, s = np.loadtxt(out / name / "word_0.dat").T
    oracle = json.load(open(out / name / "word_0_oracle.json"))
    eo, so = np.array(oracle["energies"]), np.array(oracle["projected"])
    return {"kpm": [float(np.interp(x, e, s)) for x in ENERGIES],
            "exact": [float(np.interp(x, eo, so)) for x in ENERGIES]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("data", nargs="?", default=os.environ.get("LSQUANT_QE_GOLDEN_DIR"))
    parser.add_argument("--write", action="store_true", help="store new reference values")
    parser.add_argument("--only", choices=["openmx", "wannier", "paoflow"])
    parser.add_argument("--tolerance", type=float, default=1e-3, help="absolute, e^2/h")
    args = parser.parse_args(argv)
    if not args.data or not Path(args.data).is_dir():
        print("SKIP: no data directory (argument or LSQUANT_QE_GOLDEN_DIR)")
        return 0
    data = Path(args.data).resolve()
    reference = json.load(open(REFERENCE)) if REFERENCE.is_file() else {}
    results, failures = {}, 0
    with tempfile.TemporaryDirectory() as tmp:
        for name, model in models(data).items():
            if args.only and name != args.only:
                continue
            results[name] = run(name, model, Path(tmp))
            print(f"{name}: KPM {results[name]['kpm']}  exact {results[name]['exact']}")
            if not args.write:
                if name not in reference:
                    print(f"  no reference for {name}")
                    failures += 1
                    continue
                for kind in ("kpm", "exact"):
                    diff = np.abs(np.array(results[name][kind]) - np.array(reference[name][kind]))
                    ok = bool(np.all(diff <= args.tolerance))
                    failures += not ok
                    print(f"  {kind}: max |diff| {diff.max():.2e} {'PASS' if ok else 'FAIL'}")
    if args.write:
        reference.update(results)
        reference["energies_eV_from_E_F"] = ENERGIES
        reference["settings"] = {k: COMMON[k] for k in ("basis", "projection", "expansion", "reconstruction")}
        REFERENCE.parent.mkdir(parents=True, exist_ok=True)
        REFERENCE.write_text(json.dumps(reference, indent=1) + "\n")
        print(f"wrote {REFERENCE}")
        return 0
    print("QE CONNECTORS GOLDEN", "PASS" if failures == 0 else f"FAIL ({failures})")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
