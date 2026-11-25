#!/usr/bin/env python3
"""Write the seven LSQuant run files of example 24 (paths relative to the repository root).

Usage (from the repository root):
    python examples/24_hbn_optics_dft/make_runs.py

Every run computes the moments of Tr[VX delta(E-H) VX dG-(E)/dE] and
Tr[VY delta(E-H) VY dG-(E)/dE] on a 120 x 120 k mesh in the eigenbasis with
an exact trace.

- SIESTA and PySCF (nonorthogonal, spinor-lifted, covariant velocity):
  the local basis has states up to 100-190 eV. A shift projection on
  [-21, 22] eV keeps the states that matter for light up to 15 eV; the
  velocity is projected on them and the others are clipped to the window
  edges, where they do not couple. 1600 moments on [-38, 23] eV.
- Quantum ESPRESSO + Wannier90 and QE + PAOFLOW (orthogonal, spinless):
  eight orbitals, no projection. Wannier90 uses the covariant velocity with
  the position matrix hbn_r.dat; PAOFLOW has no position matrix and uses the
  Berry-connection velocity with the orbital centres. Energies from the
  valence-band maximum of the QE run. 1700 moments on [-38, 26] eV.

Both resolutions are 60 meV. The window holds the spectrum widened by the
largest photon energy (15 eV) below the lowest band and above the
chemical potential.
"""
import json
import os

HERE = "examples/24_hbn_optics_dft"
XX = "Tr[VX delta(E-H) VX dG-(E)/dE]"
YY = "Tr[VY delta(E-H) VY dG-(E)/dE]"
VBM = {"pbe": -4.0569, "hse": -4.8265}      # eV, highest occupied level of the QE scf runs


def run(model, out, *, moments, window, projection=None, nonorthogonal=False):
    r = {"model": model,
         "basis": {"name": "eigen", "supercell": [120, 120, 1]},
         "metric": ({"policy": "nonorthogonal", "route": "lowdin", "applier": "exact"}
                    if nonorthogonal else {"policy": "orthogonal"}),
         "word": XX, "words": [XX, YY],
         "expansion": {"moments": moments, "window": list(window), "trace": "exact", "seed": 2026},
         "reconstruction": {"kernel": "jackson", "broadening_meV": 60.0},
         "output": {"legacy_layout": False, "write_chebmom": True, "dir": f"out/hbn_{out}"}}
    if projection:
        r["projection"] = {"kind": "shift", "window": list(projection)}
    return r


def main():
    cov = {"velocity_mode": "covariant"}
    lcao = dict(moments=1600, window=(-38.0, 23.0), projection=(-21.0, 22.0), nonorthogonal=True)
    files = {
        "run_siesta_pbe.json": run({"kind": "siesta", "dir": f"{HERE}/siesta", "label": "hbn", "params": cov},
                                   "siesta_pbe", **lcao),
        "run_pyscf_pbe.json": run({"kind": "pyscf", "file": f"{HERE}/pyscf/hbn_pyscf_pbe.npz", "params": cov},
                                  "pyscf_pbe", **lcao),
        "run_pyscf_hse.json": run({"kind": "pyscf", "file": f"{HERE}/pyscf/hbn_pyscf_hse.npz", "params": cov},
                                  "pyscf_hse", **lcao),
    }
    for f in ("pbe", "hse"):
        files[f"run_w90_{f}.json"] = run(
            {"kind": "wannier", "dir": f"{HERE}/qe/{f}", "seed": "hbn",
             "params": {"velocity_mode": "covariant", "fermi_energy": VBM[f]}},
            f"w90_{f}", moments=1700, window=(-38.0, 26.0))
        files[f"run_pao_{f}.json"] = run(
            {"kind": "paoflow", "file": f"{HERE}/qe/{f}/hbn_paoflow_{f}.npz",
             "params": {"velocity_mode": "berry_connection"}},
            f"pao_{f}", moments=1700, window=(-38.0, 26.0))
    for name, r in files.items():
        with open(os.path.join(HERE, name), "w") as fh:
            json.dump(r, fh, indent=2)
            fh.write("\n")
    print(" ".join(files))


if __name__ == "__main__":
    main()
