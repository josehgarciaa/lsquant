#!/usr/bin/env python3
"""Write a nearest-neighbour graphene p_z model in Wannier90 format.

Writes, inside ``graphene/``:

    graphene_hr.dat   H(R) in the Wannier90 _hr.dat layout (dense W x W per R)
    graphene.uc       the three lattice vectors (rows, Angstrom)
    graphene.xyz      the W = 2 orbital centres (Angstrom)

Lattice: a1 = (a, 0, 0), a2 = (a/2, a sqrt(3)/2, 0), a = 2.46 Angstrom.
Orbital 1 (A) sits at the origin, orbital 2 (B) at (a1 + a2)/3. The three
A -> B bonds connect cell 0 to the cells R = (0,0), (-1,0), (0,-1). The bands
are E(k) = +/- |t| |1 + exp(-i k.a1) + exp(-i k.a2)|, between -3|t| and +3|t|.

Usage:
    python make_graphene.py            # t = -2.7 eV
"""
import os
import sys

A_LATT = 2.46
T_HOP = -2.7


def main(outdir="graphene", t=T_HOP):
    os.makedirs(outdir, exist_ok=True)
    rvecs = [(-1, 0, 0), (0, -1, 0), (0, 0, 0), (0, 1, 0), (1, 0, 0)]
    # H_{ij}(R) = <i, 0 | H | j, R>
    hop = {(0, 0, 0): {(1, 2): t, (2, 1): t},
           (-1, 0, 0): {(1, 2): t},
           (1, 0, 0): {(2, 1): t},
           (0, -1, 0): {(1, 2): t},
           (0, 1, 0): {(2, 1): t}}
    with open(os.path.join(outdir, "graphene_hr.dat"), "w") as f:
        f.write("graphene p_z nearest-neighbour model, written by make_graphene.py\n")
        f.write("2\n")
        f.write(f"{len(rvecs)}\n")
        f.write(" ".join("1" for _ in rvecs) + "\n")
        for r in rvecs:
            for j in (1, 2):          # Wannier90 order: i runs fastest
                for i in (1, 2):
                    v = hop[r].get((i, j), 0.0)
                    f.write(f"{r[0]:5d}{r[1]:5d}{r[2]:5d}{i:5d}{j:5d}{v:12.6f}{0.0:12.6f}\n")
    s3 = 3 ** 0.5
    with open(os.path.join(outdir, "graphene.uc"), "w") as f:
        f.write(f"{A_LATT:.10f} 0.0 0.0\n{A_LATT / 2:.10f} {A_LATT * s3 / 2:.10f} 0.0\n0.0 0.0 20.0\n")
    with open(os.path.join(outdir, "graphene.xyz"), "w") as f:
        f.write("2\n")
        f.write("C_A 0.0 0.0 0.0\n")
        f.write(f"C_B {A_LATT / 2:.10f} {A_LATT * s3 / 6:.10f} 0.0\n")
    print(f"wrote {outdir}/graphene_hr.dat, graphene.uc, graphene.xyz (t = {t} eV)")


if __name__ == "__main__":
    main(t=float(sys.argv[1]) if len(sys.argv) > 1 else T_HOP)
