#!/usr/bin/env python3
"""Build a toy non-orthogonal tight-binding model with a known spectrum.

Start from the periodic 1D chain of N sites, H = t sum_i (|i><i+1| + h.c.),
t = -1 eV, whose eigenvalues are E_k = 2t cos(2 pi k / N). Change basis with a
local, non-unitary transformation R = I + alpha K (K random on the
nearest-neighbour graph, ``lsquant.models.nonorthogonal.make_local_transform``):
the new basis functions overlap, with

    S' = R^+ R,      H' = R^+ H R,      VX' = R^+ VX R.

The generalised problem H' c = E S' c has exactly the chain eigenvalues, so
the non-orthogonal model (H', S') has a spectrum known in closed form, while
S' is a genuine banded overlap matrix with condition number between 5 and 10.

Written under operators/ next to this script:

    chain.HAM.CSR, chain.VX.CSR                                  orthogonal reference
    nonorth_chain.HAM.CSR, nonorth_chain.S.CSR, nonorth_chain.VX.CSR   non-orthogonal model

Usage:
    python make_nonorth_chain.py [N]        # default N = 128
"""
import os
import sys

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(_HERE)), "python"))
from lsquant.models.nonorthogonal import make_local_transform, rotate_basis, write_csr_text  # noqa: E402

T, ALPHA, SEED = -1.0, 0.55, 1   # seed 1 gives cond(S') in [5, 10] at N = 128


def main(n=128):
    out = os.path.join(_HERE, "operators")
    i = np.arange(n)
    j = (i + 1) % n
    h = sp.csr_matrix((np.full(n, T, complex), (i, j)), shape=(n, n))
    h = (h + h.conj().T).tocsr()
    # VX = i[H, X] with x = i (lattice constant 1), minimum image across the boundary
    vx = sp.csr_matrix((np.full(n, 1j * T), (i, j)), shape=(n, n))
    vx = (vx + vx.conj().T).tocsr()
    write_csr_text(os.path.join(out, "chain.HAM.CSR"), h)
    write_csr_text(os.path.join(out, "chain.VX.CSR"), vx)
    r = make_local_transform(h, alpha=ALPHA, complex_entries=True, seed=SEED)
    rot = rotate_basis({"HAM": h, "VX": vx}, r)
    for name in ("HAM", "S", "VX"):
        write_csr_text(os.path.join(out, f"nonorth_chain.{name}.CSR"), rot[name])
    s_eigs = la.eigvalsh(rot["S"].toarray())
    gen = la.eigh(rot["HAM"].toarray(), rot["S"].toarray(), eigvals_only=True)
    exact = np.sort(2 * T * np.cos(2 * np.pi * np.arange(n) / n))
    print(f"N = {n}: spec(S') = [{s_eigs.min():.4f}, {s_eigs.max():.4f}], "
          f"cond(S') = {s_eigs.max() / s_eigs.min():.3f}, nnz(S') = {rot['S'].nnz}")
    print(f"max |eig(H', S') - 2t cos(2 pi k/N)| = {np.max(np.abs(gen - exact)):.2e}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 128)
