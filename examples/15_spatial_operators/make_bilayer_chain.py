#!/usr/bin/env python3
"""Build a biased bilayer chain (two-leg ladder) and its spatial operators.

Two chains ("top" and "bottom" layers) of L sites with lattice constant a,
intra-chain hopping -t, inter-layer hopping -t_perp and a layer bias +/- Delta/2:

    H = -t   sum_{l,i} (|l,i><l,i+1| + h.c.)
        -t_p sum_i     (|top,i><bot,i| + h.c.)
        + Delta/2 sum_i (|top,i><top,i| - |bot,i><bot,i|),   periodic in x.

Its bands are E(k) = -2t cos(ka) +/- s with s = sqrt(t_p^2 + Delta^2/4); both
have the velocity v(k) = 2 t a sin(ka) (eV Angstrom).

Written under operators/ next to this script with label "bilayer_L<L>", orbital
order 2*i + l (l=0 top):

    <label>.HAM.CSR    the Hamiltonian
    <label>.VX.CSR     VX = i[H, X] built from the orbital positions x = i*a
                       (minimum image across the periodic boundary); the word
                       [X,U(t)] of the mean-square displacement uses it
    <label>.A.CSR      layer projector on the top chain   (diagonal 1/0)
    <label>.B.CSR      layer projector on the bottom chain
                       (the run file executor only accepts the catalogue
                       operator names, so the projectors are called A and B)
    <label>.xyz        orbital positions (Angstrom), the input of VX

Usage:
    python make_bilayer_chain.py 256       # layer-resolved DOS (512 orbitals)
    python make_bilayer_chain.py 4096      # mean-square displacement (8192 orbitals)
"""
import os
import sys

T, T_PERP, DELTA, A = 1.0, 0.5, 0.6, 1.0


def write_csr(path, n, rows):
    """rows: list (length n) of dicts {col: complex}."""
    vals, cols, ptr = [], [], [0]
    for r in rows:
        for c in sorted(r):
            v = r[c]
            vals.append(f"{v.real:.17g} {v.imag:.17g}")
            cols.append(str(c))
        ptr.append(len(cols))
    with open(path, "w") as f:
        f.write(f"{n} {len(cols)}\n" + " ".join(vals) + "\n" + " ".join(cols) + "\n"
                + " ".join(map(str, ptr)) + "\n")


_HERE = os.path.dirname(os.path.abspath(__file__))


def main(L=4096, outdir=os.path.join(_HERE, "operators")):
    os.makedirs(outdir, exist_ok=True)
    label = f"bilayer_L{L}"
    n = 2 * L
    idx = lambda i, l: 2 * (i % L) + l
    x = [A * (o // 2) for o in range(n)]
    ham = [dict() for _ in range(n)]
    for i in range(L):
        for l in (0, 1):
            o = idx(i, l)
            ham[o][o] = complex(DELTA / 2 if l == 0 else -DELTA / 2)
            ham[o][idx(i + 1, l)] = complex(-T)
            ham[o][idx(i - 1, l)] = complex(-T)
            ham[o][idx(i, 1 - l)] = complex(-T_PERP)
    vx = [dict() for _ in range(n)]
    for o in range(n):
        for c, h in ham[o].items():
            dx = x[c] - x[o]
            dx -= A * L * round(dx / (A * L))          # minimum image
            if dx != 0.0:
                vx[o][c] = 1j * h * dx                # (i[H,X])_oc = i H_oc (x_c - x_o)
    write_csr(os.path.join(outdir, f"{label}.HAM.CSR"), n, ham)
    write_csr(os.path.join(outdir, f"{label}.VX.CSR"), n, vx)
    write_csr(os.path.join(outdir, f"{label}.A.CSR"), n,
              [{o: 1 + 0j} if o % 2 == 0 else {} for o in range(n)])
    write_csr(os.path.join(outdir, f"{label}.B.CSR"), n,
              [{o: 1 + 0j} if o % 2 == 1 else {} for o in range(n)])
    with open(os.path.join(outdir, f"{label}.xyz"), "w") as f:
        f.write(f"{n}\nbilayer chain, orbital order 2*i + layer\n")
        for o in range(n):
            f.write(f"{'T' if o % 2 == 0 else 'B'} {x[o]:.6f} 0.0 {0.0 if o % 2 == 0 else 3.35:.2f}\n")
    print(f"built {label}: L = {L}, {n} orbitals, t = {T}, t_perp = {T_PERP}, Delta = {DELTA} eV")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 4096)
