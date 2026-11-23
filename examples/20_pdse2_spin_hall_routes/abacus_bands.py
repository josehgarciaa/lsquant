"""Band quantities of an ABACUS model of PdSe2, and its deviation from a reference.

Usage:
    python abacus_bands.py RUN_DIR OUT.json [--qe PW_OUTPUT] [--shc PREFIX] [--setting NAME]

RUN_DIR holds the INPUT, STRU and OUT.<suffix>/ of one ABACUS run with
``out_mat_hs2 1`` and ``out_mat_r 1``. The script reads the model with the
LSQuant ABACUS reader and writes, with the valence-band maximum (VBM) as the
energy origin:

- the band gap, E_F - VBM (E_F of ABACUS, and the intrinsic level that
  makes the 24 x 24 mesh neutral with the Fermi-Dirac width of the run: in
  a gap of 1.2 eV the bisection of ABACUS can stop anywhere where the
  occupations are 0 or 1 to machine precision) and the quantum spin Hall gap (between bands
  NV + 4 and NV + 5, NV = 100 valence electrons with the PseudoDojo
  pseudopotentials), from the bands on Gamma-X-S-Y-Gamma-S (601 points) and
  on a 24 x 24 mesh;
- with ``--qe``, the RMS and largest deviation from the eigenvalues of a
  Quantum ESPRESSO run (a pw.x output with ``verbosity = 'high'``; its k
  points are read from the output) for the 8 valence and 8 conduction
  bands below and above the gap, each code aligned at its own VBM on those
  k points;
- with ``--shc PREFIX``, the spin Hall plateau (mean over the central half of
  the quantum spin Hall gap) and sigma in the trivial gap (0.6 eV above the
  VBM) of the files ``PREFIX_<N>_covariant.dat`` written by abacus_exact.py,
  and each curve shifted to E - VBM (``PREFIX_<N>_vbm.dat``).
"""
import argparse
import glob
import json
import os
import re

import numpy as np
import scipy.linalg as sla
from scipy.special import expit

from lsquant.models.abacus import load_abacus_model

NV = 100


def load(run_dir):
    inp = open(os.path.join(run_dir, "INPUT")).read()
    suffix = re.search(r"^\s*suffix\s+(\S+)", inp, re.M)
    suffix = suffix.group(1) if suffix else "ABACUS"
    orbital_dir = re.search(r"^\s*orbital_dir\s+(\S+)", inp, re.M)
    orbital_dir = os.path.join(run_dir, orbital_dir.group(1)) if orbital_dir else run_dir
    out = os.path.join(run_dir, f"OUT.{suffix}")
    ef = float(re.findall(r"E_Fermi\s+\S+\s+(\S+)", open(os.path.join(out, "running_scf.log")).read())[-1])
    model = load_abacus_model(os.path.join(out, "data-HR-sparse_SPIN0.csr"),
                              overlap_file=os.path.join(out, "data-SR-sparse_SPIN0.csr"),
                              stru_file=os.path.join(run_dir, "STRU"), orbital_dir=orbital_dir,
                              fermi_energy=ef, position_file=os.path.join(out, "data-rR-sparse.csr"))
    # Large bases with long-range orbitals have many cells; let the Bloch sums
    # hold them (MiB per call).
    model.bloch_memory_mb = 4096.0
    return model, ef, inp


def eigenvalues(model, kf, chunk=200):
    """Eigenvalues (eV from E_F) at fractional k."""
    b = 2 * np.pi * np.linalg.inv(model.uc).T
    out = []
    for c in range(0, len(kf), chunk):
        kc = np.atleast_2d(kf[c:c + chunk]) @ b
        out += [sla.eigvalsh(h, s) for h, s in zip(model.Hk(kc), model.Sk(kc))]
    return np.array(out)


def intrinsic_level(e_mesh, kt):
    """Chemical potential with as many holes below band NV as electrons above it."""
    def excess(mu):
        electrons = expit(-(e_mesh[:, NV:] - mu) / kt).sum()
        holes = expit((e_mesh[:, :NV] - mu) / kt).sum()
        return electrons - holes
    lo, hi = e_mesh[:, NV - 1].max(), e_mesh[:, NV].min()
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if excess(mid) < 0 else (lo, mid)
    return 0.5 * (lo + hi)


def path():
    nodes = np.array([[0, 0, 0], [.5, 0, 0], [.5, .5, 0], [0, .5, 0], [0, 0, 0], [.5, .5, 0]])
    kf = [a + t * (b - a) for a, b in zip(nodes[:-1], nodes[1:]) for t in np.linspace(0, 1, 121)[:-1]]
    return np.array(kf + [nodes[-1]])


def read_pw(path_out):
    """(k_crystal, eigenvalues eV) of a pw.x output written with ``verbosity = 'high'``.

    The k points are the ``cryst. coord.`` list that pw.x prints (any
    K_POINTS form, symmetry-reduced or not).
    """
    text = open(path_out).read()
    cryst = text[re.search(r"cryst\. coord\.\s*\n\s*k\(", text).start():]
    nk = int(re.search(r"number of k points=\s*(\d+)", text).group(1))
    kf = np.array(re.findall(r"k\(\s*\d+\) = \(\s*(\S+)\s+(\S+)\s+(\S+)\)", cryst)[:nk], float)
    blocks = re.findall(r"bands \(ev\):\s*\n(.*?)(?:\n\s*\n\s*(?:occupation|k =|the Fermi|highest|Writing))",
                        text, re.S)
    eig = [np.array(re.findall(r"-?\d+\.\d+", b), float) for b in blocks[-nk:]]
    return kf, np.array(eig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("run_dir")
    ap.add_argument("out_json")
    ap.add_argument("--qe")
    ap.add_argument("--shc")
    ap.add_argument("--setting", default=None)
    args = ap.parse_args()

    model, ef, inp = load(args.run_dir)
    n = 24
    mesh = np.array([[i / n, j / n, 0] for j in range(n) for i in range(n)])
    e_path, e_mesh = eigenvalues(model, path()), eigenvalues(model, mesh)
    both = np.concatenate([e_path, e_mesh])
    vbm = both[:, NV - 1].max()
    result = dict(code="abacus", setting=args.setting or os.path.basename(os.path.normpath(args.run_dir)),
                  pseudo="PseudoDojo NC FR 0.4 PBE standard",
                  basis=re.findall(r"\S+\.orb", open(os.path.join(args.run_dir, "STRU")).read()),
                  kT_eV=round(float(re.search(r"smearing_sigma\s+(\S+)", inp).group(1)) * 13.605693123, 4),
                  smearing=re.search(r"smearing_method\s+(\S+)", inp).group(1),
                  orbitals=int(model.num_wann), E_F_abacus_eV=ef,
                  EF_minus_VBM=float(-vbm), band_gap=float(both[:, NV].min() - vbm),
                  qsh_gap=[float(both[:, NV + 3].max() - vbm), float(both[:, NV + 4].min() - vbm)],
                  velocity="covariant",
                  connection_hermiticity_residual=model.metadata["connection_hermiticity_residual"])
    result["EF_intrinsic_minus_VBM"] = float(intrinsic_level(e_mesh, result["kT_eV"]) - vbm)
    if args.qe:
        kq, eq = read_pw(args.qe)
        ea = eigenvalues(model, kq)
        aq = eq - eq[:, NV - 1].max()
        aa = ea - ea[:, NV - 1].max()
        vb, cb = slice(NV - 8, NV), slice(NV, NV + 8)
        dvb, dcb = aa[:, vb] - aq[:, vb], aa[:, cb] - aq[:, cb]
        result["dev_vs_qe"] = dict(rms_vb=float(np.sqrt((dvb ** 2).mean())), rms_cb=float(np.sqrt((dcb ** 2).mean())),
                                   max=float(max(np.abs(dvb).max(), np.abs(dcb).max())),
                                   reference=os.path.abspath(args.qe), k_points=len(kq))
    if args.shc:
        lo, hi = result["qsh_gap"]
        q = (hi - lo) / 4
        result["plateau"], result["sigma_trivial"] = {}, {}
        for f in sorted(glob.glob(f"{args.shc}_*_covariant.dat")):
            mesh_n = re.search(r"_(\d+)_covariant\.dat$", f).group(1)
            d = np.loadtxt(f)
            e = d[:, 0] - vbm
            sel = (e > lo + q) & (e < hi - q)
            result["plateau"][mesh_n] = float(d[sel, 1].mean())
            result["sigma_trivial"][mesh_n] = float(np.interp(0.6, e, d[:, 1]))
            out = f.replace("_covariant.dat", "_vbm.dat")
            np.savetxt(out, np.column_stack([e, d[:, 1]]), fmt="%.6f %.8e",
                       header="E-VBM (eV)  sigma^z_xy (e^2/h, Pauli sigma_z); " + open(f).readline()[2:].strip())
            result["curve"] = os.path.abspath(out)
    json.dump(result, open(args.out_json, "w"), indent=1)
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
