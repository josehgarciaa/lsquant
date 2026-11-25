#!/usr/bin/env python3
"""Check a Wannier90 model of hBN against the DFT bands on the mesh and print the decay of H(R).

Usage (in qe/pbe or qe/hse, after wannier90.x):
    python ../check_wannier.py hbn

Reads <seed>_hr.dat, <seed>.win (k list) and <seed>.eig. Prints
- the largest difference between the four valence bands of the model and
  of DFT over the 144 mesh points;
- for the four conduction Wannier bands, the largest and median difference
  to the nearest DFT eigenvalue (the DFT list also holds the vacuum states);
- the bands of both at K and M;
- the largest |H(R)| on each shell of |R|.
Bloch sums use the degeneracies of _hr.dat; H(R) is not folded with
_wsvec.dat, which only matters between mesh points.
"""
import sys

import numpy as np

A1 = np.array([2.504, 0.0, 0.0])
A2 = np.array([-1.252, 2.168528, 0.0])


def read_hr(seed):
    lines = open(seed + "_hr.dat").read().split("\n")
    nw, nr = int(lines[1]), int(lines[2])
    nl = (nr + 14) // 15
    deg = np.array(" ".join(lines[3:3 + nl]).split(), int)
    data = np.array([l.split() for l in lines[3 + nl:] if l.strip()], float)
    r = data[:, :3].astype(int).reshape(nr, nw * nw, 3)[:, 0]
    h = (data[:, 5] + 1j * data[:, 6]).reshape(nr, nw, nw).transpose(0, 2, 1)
    return r, deg, h


def main():
    seed = sys.argv[1]
    r, deg, h = read_hr(seed)
    block = open(seed + ".win").read().split("begin kpoints")[1].split("end kpoints")[0]
    kpts = np.array([l.split() for l in block.strip().split("\n")], float)
    eig = np.loadtxt(seed + ".eig")
    e_dft = eig[:, 2].reshape(len(kpts), -1)

    def bands(k):
        return np.linalg.eigvalsh(np.einsum("r,rij->ij", np.exp(2j * np.pi * r @ k) / deg, h))

    model = np.array([bands(k) for k in kpts])
    valence = np.abs(model[:, :4] - e_dft[:, :4]).max()
    nearest = np.array([[np.abs(e_dft[i] - e).min() for e in model[i, 4:]] for i in range(len(kpts))])
    print(f"valence bands, max |W90 - DFT| on the mesh: {valence:.2e} eV")
    print("conduction Wannier bands 5-8, |W90 - nearest DFT| max:", np.round(nearest.max(axis=0), 3),
          " median:", np.round(np.median(nearest, axis=0), 3))
    for name, k in (("K", [1 / 3, 1 / 3, 0]), ("M", [0.5, 0, 0])):
        i = int(np.argmin(np.linalg.norm(kpts - k, axis=1)))
        print(f"{name}: W90", np.round(model[i], 2), " DFT", np.round(e_dft[i, :6], 2))
    dist = np.linalg.norm(r[:, :1] * A1 + r[:, 1:2] * A2, axis=1)
    size = np.abs(h).max(axis=(1, 2))
    for d0 in sorted(set(np.round(dist, 1)))[:12]:
        print(f"|R| = {d0:5.1f} A  max |H(R)| = {size[np.isclose(np.round(dist, 1), d0)].max():.4f} eV")
    print(f"|R| > {dist.max() - 3:.1f} A: max |H(R)| = {size[dist > dist.max() - 3].max():.2e} eV")


if __name__ == "__main__":
    main()
