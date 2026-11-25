#!/usr/bin/env python3
"""Projected PAOFLOW model of hBN from the Quantum ESPRESSO run in this directory.

Usage (in qe/pbe or qe/hse, after scf, open_grid.x and projwfc.x):
    python ../run_paoflow.py OUT.npz [--pthr 0.6] [--outdir tmp]

PAOFLOW projects the Bloch states of the open_grid run (12 x 12 k points,
no symmetry) on the atomic orbitals of the pseudopotentials (B and N: 2s,
2p), keeps the bands whose projectability exceeds pthr, and builds the
eight-orbital H(R). write_paoflow_model stores it for LSQuant. No spin-orbit
coupling: the model is spinless, so optical runs use spin_degeneracy = 2.
"""
import argparse
import os
import re

from PAOFLOW import PAOFLOW

from lsquant.models.orthogonal import write_paoflow_model


def prepare_open_grid_xml(outdir):
    """Give the open_grid.x data file the two scf elements PAOFLOW reads.

    open_grid.x writes no <input> section (PAOFLOW reads the spin-orbit flag
    there) and lists the 144 k points explicitly instead of the
    Monkhorst-Pack grid (PAOFLOW reads nk1, nk2, nk3 there). Both are copied
    from the scf data file; the k points and bands of open_grid are kept.
    """
    target = os.path.join(outdir, "hbn_open.save", "data-file-schema.xml")
    text = open(target).read()
    scf = open(os.path.join(outdir, "hbn.save", "data-file-schema.xml")).read()
    if "<input>" not in text:
        block = re.search(r"\s*<input>.*?</input>", scf, re.S).group(0)
        text = text.replace("\n  <output>", block + "\n  <output>", 1)
    current = re.search(r"<starting_k_points>.*?</starting_k_points>", text, re.S).group(0)
    if "<monkhorst_pack" not in current:
        grid = re.search(r"<starting_k_points>.*?</starting_k_points>", scf, re.S).group(0)
        text = re.sub(r"<starting_k_points>.*?</starting_k_points>", lambda m: grid, text, count=1, flags=re.S)
    open(target, "w").write(text)


def compare_bands(p):
    """Print the PAO bands against the DFT bands at Gamma, M and K of the 12 x 12 mesh.

    PAOFLOW reproduces the bands it kept (projectability above pthr, counted
    by band index) and pushes the others to the shift. Energies relative to
    the DFT Fermi energy.
    """
    import numpy as np
    arrays, attrs = p.data_controller.data_dicts()
    hk = np.fft.fftn(np.asarray(arrays["HRs"])[..., 0], axes=(2, 3, 4))       # PAOFLOW's own transform
    dft = np.asarray(arrays["my_eigsmat"])[:, :, 0] if "my_eigsmat" in arrays else None
    n1, n2, _ = hk.shape[2:5]
    for name, (i, j) in (("Gamma", (0, 0)), ("M", (n1 // 2, 0)), ("K", (n1 // 3, n2 // 3))):
        pao = np.linalg.eigvalsh(hk[:, :, i, j, 0])
        line = f"{name:5s} PAO " + " ".join(f"{e:7.2f}" for e in pao)
        if dft is not None:
            line += "   DFT " + " ".join(f"{e:7.2f}" for e in dft[:10, i * n2 + j])
        print(line)
    print("kept bands (pthr):", attrs.get("bnd"), " shift:", attrs.get("shift"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out")
    ap.add_argument("--pthr", type=float, default=0.6)
    ap.add_argument("--outdir", default="tmp", help="QE outdir of the open_grid run")
    a = ap.parse_args()
    prepare_open_grid_xml(a.outdir)
    p = PAOFLOW.PAOFLOW(savedir=os.path.join(a.outdir, "hbn_open.save"), outputdir="paoflow_out",
                        npool=1, verbose=True)
    p.read_atomic_proj_QE()
    p.projectability(pthr=a.pthr)
    p.pao_hamiltonian()
    compare_bands(p)
    write_paoflow_model(p, a.out)
    print("written", a.out)


if __name__ == "__main__":
    main()
