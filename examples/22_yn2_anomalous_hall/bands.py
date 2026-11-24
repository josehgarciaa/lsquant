"""YN2 bands along G-M-K-G from the SIESTA model and the Quantum ESPRESSO + Wannier90 model.

Usage:
    python bands.py SIESTA_DIR W90_DIR OUT.npz --fermi-energy EF [--label yn2] [--seed y]

Energies relative to E_F: the E_F of the .TSHS for SIESTA, --fermi-energy
(the Quantum ESPRESSO scf Fermi energy, eV) for the Wannier90 model. Writes
the path coordinate and the bands of both models to OUT.npz.
"""
import argparse

import numpy as np
import scipy.linalg as sla

from lsquant.models.orthogonal import write_wannier_geometry
from lsquant.models.siesta import load_siesta_model
from lsquant.models.wannier import load_wannier_model

NODES = np.array([[0, 0, 0], [0.5, 0, 0], [1 / 3, 1 / 3, 0], [0, 0, 0]])  # G, M, K, G


def path(uc, points=120):
    """Cartesian k along the nodes and the path coordinate, in 1/Angstrom."""
    recip = 2 * np.pi * np.linalg.inv(uc).T
    frac = np.vstack([a + np.linspace(0, 1, points, endpoint=False)[:, None] * (b - a)
                      for a, b in zip(NODES[:-1], NODES[1:])] + [NODES[-1:]])
    k = frac @ recip
    x = np.r_[0, np.cumsum(np.linalg.norm(np.diff(k, axis=0), axis=1))]
    return k, x


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("siesta_dir")
    ap.add_argument("w90_dir")
    ap.add_argument("out")
    ap.add_argument("--fermi-energy", type=float, required=True)
    ap.add_argument("--label", default="yn2")
    ap.add_argument("--seed", default="y")
    a = ap.parse_args()
    siesta = load_siesta_model(a.siesta_dir, a.label, dipoles=False)
    k, x = path(siesta.uc)
    e_siesta = np.array([sla.eigvalsh(h, s) for h, s in zip(siesta.Hk(k), siesta.Sk(k))])
    write_wannier_geometry(a.w90_dir, a.seed)
    wannier = load_wannier_model(a.w90_dir, a.seed)
    kw, xw = path(np.asarray(wannier.uc, float))
    e_wannier = np.linalg.eigvalsh(wannier.Hk(kw)) - a.fermi_energy
    np.savez(a.out, x_siesta=x, E_siesta=e_siesta, x_wannier=xw, E_wannier=e_wannier)


if __name__ == "__main__":
    main()
