#!/usr/bin/env python3
"""Write a Kane-Mele model with Rashba coupling and exchange as csr files.

Honeycomb lattice, bond length a = 1 Angstrom, two sites (A, B) and two spin
components per cell. Orbital order: 4*cell + 2*site + spin (spin 0 = up).

    H = -t   sum_<ij>  c_i^+ c_j                                (nearest neighbours)
        + i lambda_SO sum_<<ij>> nu_ij c_i^+ sigma_z c_j        (intrinsic SOC)
        + i lambda_R  sum_<ij> c_i^+ (sigma x d_ij)_z c_j        (Rashba)
        + M sum_i c_i^+ sigma_z c_i                              (exchange along z)

with nu_ij = +-1 the Kane-Mele sign and d_ij the unit bond vector. The
parameters are set at the top of the file.

For an L x L supercell, periodic in both directions, it writes under
operators/ with the label km_L<L>:

    <label>.HAM.CSR   the Hamiltonian (eV)
    <label>.VX.CSR    VX = i[H, X] (eV Angstrom), from the bond vectors
    <label>.VY.CSR    VY = i[H, Y]
    <label>.SZ.CSR    the Pauli spin operator sigma_z (diagonal +-1)
    <label>.JXSZ.CSR  the spin current {VX, SZ}/2
    <label>.json      the parameters, the number of cells and the cell area

The same Bloch blocks H(k) give the exact reference of oracle.py: an L x L
supercell has exactly the k points of an L x L mesh.

Usage:
    python make_kane_mele.py 16
"""
import json
import os
import sys

import numpy as np
import scipy.sparse as sp

T, LAMBDA_SO, LAMBDA_R, M_EXCHANGE = 1.0, 0.1, 0.05, 0.05     # eV

LATTICE = np.array([[np.sqrt(3) / 2, 1.5, 0.0], [-np.sqrt(3) / 2, 1.5, 0.0], [0.0, 0.0, 10.0]])
SITES = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])          # A, B (Angstrom)
CELL_AREA = float(abs(np.linalg.det(LATTICE[:2, :2])))          # 3 sqrt(3)/2 Angstrom^2
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1.0, -1.0]).astype(complex)


def bloch_blocks():
    """Return {R: 4x4 block} with H = sum_R block_R between cell 0 and cell R."""
    blocks = {}

    def add(r, block):
        r = tuple(r)
        blocks[r] = blocks.get(r, np.zeros((4, 4), complex)) + block

    # Nearest neighbours A(0) -> B(R), with hopping and Rashba.
    for r in [(0, 0), (-1, 0), (0, -1)]:
        d = SITES[1] + np.array(r) @ LATTICE[:2] - SITES[0]
        d = d / np.linalg.norm(d)
        hop = -T * np.eye(2) + 1j * LAMBDA_R * (SX * d[1] - SY * d[0])
        block = np.zeros((4, 4), complex)
        block[0:2, 2:4] = hop
        add(r, block)
        add((-r[0], -r[1]), block.conj().T)
    # Next-nearest neighbours: i lambda nu sigma_z on A, opposite sign on B.
    for r, nu in [((1, 0), 1), ((0, 1), -1), ((-1, 1), 1)]:
        block = 1j * LAMBDA_SO * nu * np.kron(np.diag([1.0, -1.0]), SZ)
        add(r, block)
        add((-r[0], -r[1]), block.conj().T)
    add((0, 0), M_EXCHANGE * np.kron(np.eye(2), SZ))
    return blocks


def bond_vector(r, a, b):
    """Vector from orbital a of cell 0 to orbital b of cell R (Angstrom)."""
    return np.array(r) @ LATTICE[:2] + SITES[b // 2] - SITES[a // 2]


def supercell_operators(L):
    """Return the sparse H, VX, VY of the L x L supercell, orbital 4*cell+o."""
    blocks = bloch_blocks()
    cells = np.arange(L * L)
    ci, cj = cells // L, cells % L
    rows, cols, h, vx, vy = [], [], [], [], []
    for r, block in blocks.items():
        target = ((ci + r[0]) % L) * L + (cj + r[1]) % L
        for a in range(4):
            for b in range(4):
                if block[a, b] == 0:
                    continue
                d = bond_vector(r, a, b)
                rows.append(4 * cells + a)
                cols.append(4 * target + b)
                h.append(np.full(L * L, block[a, b]))
                vx.append(np.full(L * L, 1j * block[a, b] * d[0]))   # (i[H,X])_ab = i H_ab (x_b - x_a)
                vy.append(np.full(L * L, 1j * block[a, b] * d[1]))
    n = 4 * L * L
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    mats = [sp.csr_matrix((np.concatenate(v), (rows, cols)), shape=(n, n)) for v in (h, vx, vy)]
    for m in mats:
        m.sum_duplicates()
        m.sort_indices()
    return mats


def write_csr(path, m):
    """Write the text CSR format read by model.kind = csr."""
    m = sp.csr_matrix(m)
    m.sort_indices()
    vals = np.column_stack((m.data.real, m.data.imag)).ravel()
    with open(path, "w") as f:
        f.write(f"{m.shape[0]} {m.nnz}\n")
        f.write(" ".join(f"{v:.17g}" for v in vals) + "\n")
        f.write(" ".join(map(str, m.indices)) + "\n")
        f.write(" ".join(map(str, m.indptr)) + "\n")


def main(L=16, outdir=None):
    outdir = outdir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "operators")
    os.makedirs(outdir, exist_ok=True)
    label = f"km_L{L}"
    H, VX, VY = supercell_operators(L)
    n = H.shape[0]
    assert abs(H - H.getH()).max() < 1e-14, "H is not Hermitian"
    sz = sp.diags(np.tile([1.0, -1.0], n // 2).astype(complex), format="csr")
    jxsz = 0.5 * (VX @ sz + sz @ VX)
    for name, m in [("HAM", H), ("VX", VX), ("VY", VY), ("SZ", sz), ("JXSZ", jxsz)]:
        write_csr(os.path.join(outdir, f"{label}.{name}.CSR"), m)
    meta = dict(label=label, L=L, orbitals=n, num_cells=L * L, cell_area_A2=CELL_AREA,
                t_eV=T, lambda_SO_eV=LAMBDA_SO, lambda_R_eV=LAMBDA_R, M_eV=M_EXCHANGE,
                bond_length_A=1.0, orbital_order="4*cell + 2*site + spin")
    with open(os.path.join(outdir, f"{label}.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(f"wrote {label}: {n} orbitals, {H.nnz} Hamiltonian entries, cell area {CELL_AREA:.6f} A^2")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 16)
