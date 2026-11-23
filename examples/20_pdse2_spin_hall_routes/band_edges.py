"""Band-edge quantities of a PdSe2 electronic structure, with the VBM as origin.

Usage:
    python band_edges.py siesta DIR OUT_PREFIX [--label pdse2] [--electrons 100]
    python band_edges.py qe BANDS_XML MESH_XML OUT_PREFIX [--electrons 100]

For SIESTA the bands come from the LSQuant SIESTA model (H and S of the .TSHS);
for Quantum ESPRESSO from the data-file-schema.xml of a `bands` run on the
path and of an `nscf` run on the 24 x 24 mesh.

Writes OUT_PREFIX.json and OUT_PREFIX_path.npy (path bands minus the VBM). The
path is G-X-S-Y-G-S with 30 points per segment (151 points); the mesh is
24 x 24 and Gamma-centred. Quantities, in eV:

- VBM, CBM: highest occupied and lowest empty band, given the electron count;
- band_gap = CBM - VBM;
- EF_code_minus_VBM: the Fermi level the code printed, minus the VBM;
- EF_intrinsic_minus_VBM: the chemical potential of the neutral insulator at
  kT = 0.01 eV, from p = n with Boltzmann tails; within a few meV of mid-gap.
  In a gapped system the Fermi level a code prints can be anywhere in the gap,
  because every value in the gap gives the right electron count;
- qsh_gap: the widest global gap between consecutive bands that starts
  1.5 to 2.5 eV above the VBM, [lower edge, upper edge, index above the CBM band].
"""
import argparse
import json
import xml.etree.ElementTree as ET

import numpy as np
from scipy.special import logsumexp

HARTREE = 27.211386245988
KT = 0.01


def path_and_mesh():
    pts = {"G": [0, 0, 0], "X": [0.5, 0, 0], "S": [0.5, 0.5, 0], "Y": [0, 0.5, 0]}
    kp = []
    for a, b in zip("GXSYG", "XSYGS"):
        for t in np.linspace(0, 1, 30, endpoint=False):
            kp.append((1 - t) * np.array(pts[a]) + t * np.array(pts[b]))
    kp.append(pts["S"])
    g = np.arange(24) / 24
    mesh = np.array([[a, b, 0] for a in g for b in g])
    return np.array(kp, float), mesh


def siesta_bands(directory, label):
    import scipy.linalg as sl
    from lsquant.models.bloch import bloch_sum
    from lsquant.models.siesta import load_siesta_model
    m = load_siesta_model(directory, label)
    B = 2 * np.pi * np.linalg.inv(m.uc).T

    def eig(kf):
        kc = kf @ B
        H = bloch_sum(m.r_cart, m.tau, m.hoppings, kc)
        S = bloch_sum(m.r_cart, m.tau, m.overlaps, kc)
        return np.array([sl.eigh(h, s, eigvals_only=True) for h, s in zip(H, S)])

    kp, mesh = path_and_mesh()
    # The TSHS has E_F = 0, so the printed Fermi level is the origin of these energies.
    return eig(kp), eig(mesh), 0.0


def qe_xml(path):
    root = ET.parse(path).getroot()
    bs = root.find("output/band_structure")
    ef = float(bs.find("fermi_energy").text) * HARTREE
    ks = bs.findall("ks_energies")
    E = np.array([[float(x) for x in k.find("eigenvalues").text.split()] for k in ks]) * HARTREE
    return E, ef


def qe_bands(bands_xml, mesh_xml):
    Ep, _ = qe_xml(bands_xml)
    Em, ef = qe_xml(mesh_xml)  # the nscf run keeps the Fermi level of the scf run
    return Ep - ef, Em - ef, 0.0


def edges(Ep, Em, nel, ef_code):
    nv = nel  # spinor bands: one electron per band
    vbm = max(Ep[:, nv - 1].max(), Em[:, nv - 1].max())
    cbm = min(Ep[:, nv].min(), Em[:, nv].min())
    Ev = np.concatenate([Ep[:, :nv].ravel(), Em[:, :nv].ravel()])
    Ec = np.concatenate([Ep[:, nv:].ravel(), Em[:, nv:].ravel()])
    mu = 0.5 * KT * (logsumexp(Ev / KT) - logsumexp(-Ec / KT))
    best = None
    for n in range(nv, Em.shape[1] - 1):
        lo = max(Em[:, n].max(), Ep[:, n].max()) - vbm
        hi = min(Em[:, n + 1].min(), Ep[:, n + 1].min()) - vbm
        if 1.5 < lo < 2.5 and hi > lo and (best is None or hi - lo > best[1] - best[0]):
            best = (lo, hi, n - nv + 1)
    return dict(VBM=vbm, CBM=cbm, band_gap=cbm - vbm, EF_code_minus_VBM=ef_code - vbm,
                EF_intrinsic_minus_VBM=mu - vbm, qsh_gap=best, electrons=nel)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("kind", choices=["siesta", "qe"])
    ap.add_argument("inputs", nargs="+")
    ap.add_argument("--label", default="pdse2")
    ap.add_argument("--electrons", type=int, default=100)
    a = ap.parse_args()
    if a.kind == "siesta":
        Ep, Em, ef = siesta_bands(a.inputs[0], a.label)
    else:
        Ep, Em, ef = qe_bands(a.inputs[0], a.inputs[1])
    out = a.inputs[-1]
    res = edges(Ep, Em, a.electrons, ef)
    np.save(out + "_path.npy", Ep - res["VBM"])
    np.save(out + "_mesh.npy", Em - res["VBM"])
    res["path_bands"] = out + "_path.npy"
    res["mesh_bands"] = out + "_mesh.npy"
    json.dump(res, open(out + ".json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
