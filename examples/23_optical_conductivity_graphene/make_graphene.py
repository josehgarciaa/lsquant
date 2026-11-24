#!/usr/bin/env python3
"""Write the Hamiltonian and the velocity of nearest-neighbour graphene as CSR text files.

Problem: the optical conductivity needs H and the velocity operator
v_x = i [H, x] of a periodic supercell, in the sparse text format of the
`csr` model kind. Approach: one p_z orbital per atom, hopping t between
nearest neighbours, n x n cells with periodic boundaries; the element of v_x
between orbitals i and j is i H_ij (x_j - x_i), with the displacement taken
along the bond (minimum image), so it is well defined across the boundary.

Units: energies in eV, lengths in Angstrom, v_x in eV Angstrom (hbar = 1).

Pitfall: the supercell must be larger than two cells in each direction, or a
site reaches the same neighbour through two different bonds.

Usage:
    python make_graphene.py [cells_per_side] [output_dir]
"""
import os
import sys

import numpy as np
from scipy.sparse import coo_matrix

T_HOP = -2.7                       # eV
A_LAT = 2.46                       # Angstrom
A1 = np.array([A_LAT, 0.0])
A2 = np.array([A_LAT / 2, A_LAT * np.sqrt(3) / 2])
CELL_AREA = A1[0] * A2[1]          # Angstrom^2
DELTA = (A1 + A2) / 3              # from sublattice A to sublattice B in the same cell


def operators(n):
    """Return (H, v_x) of an n x n supercell as CSR matrices of dimension 2 n^2."""
    if n < 3:
        raise ValueError("cells_per_side must be at least 3")
    i1, i2 = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    site_a = 2 * (i1 * n + i2).ravel()
    rows, cols, bond_x = [], [], []
    # The three neighbours of an A site: B in the same cell, in cell - a1 and in cell - a2.
    for d1, d2 in ((0, 0), (-1, 0), (0, -1)):
        site_b = 2 * (((i1 + d1) % n) * n + (i2 + d2) % n).ravel() + 1
        bond = DELTA + d1 * A1 + d2 * A2
        rows += [site_a, site_b]
        cols += [site_b, site_a]
        bond_x += [np.full(site_a.size, bond[0]), np.full(site_a.size, -bond[0])]
    rows, cols, bond_x = map(np.concatenate, (rows, cols, bond_x))
    dim = 2 * n * n
    ham = coo_matrix((np.full(rows.size, T_HOP, dtype=complex), (rows, cols)), shape=(dim, dim)).tocsr()
    vel = coo_matrix((1j * T_HOP * bond_x, (rows, cols)), shape=(dim, dim)).tocsr()
    return ham, vel


def write_csr(path, matrix):
    """One matrix in the text format of the `csr` model kind."""
    matrix.sort_indices()
    with open(path, "w") as fh:
        fh.write(f"{matrix.shape[0]} {matrix.nnz}\n")
        np.savetxt(fh, np.column_stack((matrix.data.real, matrix.data.imag)), fmt="%.17g %.17g", newline=" ")
        fh.write("\n")
        np.savetxt(fh, matrix.indices, fmt="%d", newline=" ")
        fh.write("\n")
        np.savetxt(fh, matrix.indptr, fmt="%d", newline=" ")
        fh.write("\n")


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "operators")
    os.makedirs(out, exist_ok=True)
    ham, vel = operators(n)
    write_csr(os.path.join(out, "graphene.HAM.CSR"), ham)
    write_csr(os.path.join(out, "graphene.VX.CSR"), vel)
    print(f"{n} x {n} cells, {ham.shape[0]} orbitals, cell area {CELL_AREA:.6f} A^2 -> {out}")


if __name__ == "__main__":
    main()
