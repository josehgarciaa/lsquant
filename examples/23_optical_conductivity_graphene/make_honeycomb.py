#!/usr/bin/env python3
"""Write H, v_x and v_y of a nearest-neighbour honeycomb model as CSR text files.

Problem: the optical tests need graphene and hBN in the same model, the
velocity along x and along y (the response of a hexagonal crystal must not
depend on the in-plane polarisation), and, for the Hall check, the Haldane
model. Approach: one orbital per site; on-site energy +Delta on sublattice A
(B atom of hBN) and -Delta on sublattice B (N); hopping t between nearest
neighbours; optionally a second-neighbour hopping t2 exp(+-i phi) (Haldane).
n x n cells with periodic boundaries. The element of v_a between orbitals i
and j is i H_ij (a_j - a_i), with the displacement along the bond (minimum
image), so v = i [H, r] is well defined across the boundary.

Units: energies in eV, lengths in Angstrom, velocities in eV Angstrom
(hbar = 1).

Pitfall: the supercell must have at least three cells per side (five with
t2), or a site reaches the same neighbour through two different bonds.

Usage:
    python make_honeycomb.py NAME CELLS OUT_DIR [--delta D] [--t T] [--a A] [--t2 T2 --phi PHI]
Writes OUT_DIR/NAME.HAM.CSR, NAME.VX.CSR and NAME.VY.CSR.
"""
import argparse
import os

import numpy as np
from scipy.sparse import coo_matrix


GRAPHENE = dict(delta=0.0, t=-2.7, a=2.46)
HBN = dict(delta=3.625, t=-2.33, a=2.50)    # T. Galvani's test values: +Delta on B, -Delta on N


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


def lattice(a):
    """Primitive vectors, cell area and the A -> B vector of the same cell."""
    a1 = np.array([a, 0.0])
    a2 = np.array([a / 2, a * np.sqrt(3) / 2])
    return a1, a2, a1[0] * a2[1], (a1 + a2) / 3


def bloch(k, delta, t, a, t2=0.0, phi=0.0):
    """H(k) and dH/dk_x, dH/dk_y (2 x 2 per k) in the convention of `operators`.

    k: array (..., 2). The phase is exp(i k . d) with d the bond vector, so
    dH/dk = i d H_bond: the Bloch form of v = i [H, r].
    """
    a1, a2, _, delta_ab = lattice(a)
    k = np.asarray(k, float)
    shape = k.shape[:-1]
    h = np.zeros(shape + (2, 2), complex)
    dx = np.zeros_like(h)
    dy = np.zeros_like(h)
    h[..., 0, 0], h[..., 1, 1] = delta, -delta
    for d1, d2 in ((0, 0), (-1, 0), (0, -1)):
        bond = delta_ab + d1 * a1 + d2 * a2              # from A to B
        term = t * np.exp(1j * (k @ bond))
        h[..., 0, 1] += term
        dx[..., 0, 1] += 1j * bond[0] * term
        dy[..., 0, 1] += 1j * bond[1] * term
    for sub, sign in ((0, 1), (1, -1)):
        for vec in (a1, a2 - a1, -a2):                   # counter-clockwise second neighbours
            for s in (1, -1):
                term = t2 * np.exp(1j * s * sign * phi) * np.exp(1j * (k @ (s * vec)))
                h[..., sub, sub] += term
                dx[..., sub, sub] += 1j * s * vec[0] * term
                dy[..., sub, sub] += 1j * s * vec[1] * term
    for m in (h, dx, dy):
        m[..., 1, 0] = np.conj(m[..., 0, 1])
    return h, dx, dy


def operators(n, delta, t, a, t2=0.0, phi=0.0):
    """Return (H, v_x, v_y) of an n x n supercell as CSR matrices of dimension 2 n^2."""
    if n < (5 if t2 else 3):
        raise ValueError("too few cells per side for the range of the hopping")
    a1, a2, _, delta_ab = lattice(a)
    i1, i2 = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    site = lambda c1, c2, sub: 2 * ((c1 % n) * n + c2 % n).ravel() + sub
    rows, cols, vals, bonds = [], [], [], []

    def add(src, dst, value, bond):
        rows.append(src)
        cols.append(dst)
        vals.append(np.full(src.size, value, complex))
        bonds.append(np.tile(bond, (src.size, 1)))

    add(site(i1, i2, 0), site(i1, i2, 0), delta, np.zeros(2))
    add(site(i1, i2, 1), site(i1, i2, 1), -delta, np.zeros(2))
    for d1, d2 in ((0, 0), (-1, 0), (0, -1)):
        bond = delta_ab + d1 * a1 + d2 * a2
        add(site(i1, i2, 0), site(i1 + d1, i2 + d2, 1), t, bond)
        add(site(i1 + d1, i2 + d2, 1), site(i1, i2, 0), t, -bond)
    if t2:
        # H_ij = t2 exp(i s phi) for a hop from j to i along +vec (Bloch phase exp(i k.vec) on H_ij with
        # r_j - r_i = vec), the same convention as `bloch`.
        for sub, sign in ((0, 1), (1, -1)):
            for (c1, c2), vec in (((1, 0), a1), ((-1, 1), a2 - a1), ((0, -1), -a2)):
                for s in (1, -1):
                    add(site(i1, i2, sub), site(i1 + s * c1, i2 + s * c2, sub),
                        t2 * np.exp(1j * s * sign * phi), s * vec)
    rows, cols, vals, bonds = map(np.concatenate, (rows, cols, vals, bonds))
    dim = 2 * n * n
    mat = lambda data: coo_matrix((data, (rows, cols)), shape=(dim, dim)).tocsr()
    return mat(vals), mat(1j * vals * bonds[:, 0]), mat(1j * vals * bonds[:, 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("name")
    ap.add_argument("cells", type=int)
    ap.add_argument("out")
    ap.add_argument("--delta", type=float, default=0.0)
    ap.add_argument("--t", type=float, default=-2.7)
    ap.add_argument("--a", type=float, default=2.46)
    ap.add_argument("--t2", type=float, default=0.0)
    ap.add_argument("--phi", type=float, default=0.0)
    p = ap.parse_args()
    os.makedirs(p.out, exist_ok=True)
    ham, vx, vy = operators(p.cells, p.delta, p.t, p.a, p.t2, p.phi)
    for op, m in (("HAM", ham), ("VX", vx), ("VY", vy)):
        write_csr(os.path.join(p.out, f"{p.name}.{op}.CSR"), m)
    area = lattice(p.a)[2]
    print(f"{p.name}: {p.cells} x {p.cells} cells, {ham.shape[0]} orbitals, cell area {area:.6f} A^2 -> {p.out}")


if __name__ == "__main__":
    main()
