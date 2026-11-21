#!/usr/bin/env python3
"""The atom-block operator layout and checkpointing, from run files (native build).

Builds the nearest-neighbour square lattice (t = -1 eV, periodic, one orbital per atom) in
Python, writes its CSR file and the atom positions, and runs the same density-of-states
word through `python -m lsquant run` three ways:

  csr         the default sparse-matrix operator;
  atom_block  numerics.operator = "atom_block": atoms reordered in spatial domains (Morton
              order) and the Hamiltonian stored as a dictionary of atom blocks;
  checkpoint  atom_block with expansion.checkpoint = true (the engine writes the moment
              tensor after every trace state and can resume from it).

Part 1 (L = 64, exact trace over all 4096 sites, M = 256) checks that the three runs give
the same moments: the trace is invariant under the reordering, and the checkpointed run is
byte-identical to the plain atom_block run. Part 2 (L = 1024, 2 random-phase vectors,
M = 512) reports the wall time of the csr and atom_block runs.

The atom_block layout exists in the native build profile only:
    cmake --preset native && cmake --build --preset native
    PYTHONPATH=build/native/lsquant_core:python python3 examples/performance/atom_block.py [out_dir]
"""
from __future__ import annotations

import filecmp
import json
import os
import subprocess
import sys
import time

import numpy as np
import scipy.sparse as sp

from lsquant.models.nonorthogonal import write_csr_text


def square_lattice(L: int, out: str) -> None:
    n = L * L
    rows, cols = [], []
    for y in range(L):
        for x in range(L):
            i = x + L * y
            for j in ((x + 1) % L + L * y, (x - 1) % L + L * y, x + L * ((y + 1) % L), x + L * ((y - 1) % L)):
                rows.append(i)
                cols.append(j)
    h = sp.csr_matrix((np.full(len(rows), -1.0 + 0j), (rows, cols)), shape=(n, n))
    os.makedirs(out, exist_ok=True)
    write_csr_text(os.path.join(out, "square.HAM.CSR"), h)
    with open(os.path.join(out, "positions.xyz"), "w") as fh:
        fh.write(f"{n}\nsquare lattice\n")
        for y in range(L):
            for x in range(L):
                fh.write(f"C {x:.1f} {y:.1f} 0.0\n")


def runfile(model_dir, out, *, moments, trace, states, operator, checkpoint=False):
    rf = {
        "model": {"kind": "csr", "dir": model_dir, "label": "square"},
        "basis": {"name": "real"},
        "metric": {"policy": "orthogonal"},
        "word": "Tr[delta(E-H)]",
        "expansion": {"moments": moments, "window": [-4.2, 4.2], "trace": trace,
                      "states": states, "seed": 12345, "checkpoint": checkpoint},
        "output": {"dir": out, "legacy_names": True},
    }
    if operator == "atom_block":
        rf["numerics"] = {"operator": "atom_block"}
        rf["model"]["positions"] = os.path.join(model_dir, "positions.xyz")
        rf["model"]["orbitals_per_atom"] = 1
    return rf


def run(rf, path):
    with open(path, "w") as fh:
        json.dump(rf, fh, indent=2)
    t0 = time.perf_counter()
    subprocess.run([sys.executable, "-m", "lsquant", "run", path], check=True, stdout=subprocess.DEVNULL)
    return time.perf_counter() - t0


def moments_file(d):
    return os.path.join(d, [f for f in os.listdir(d) if f.endswith(".chebmom1D")][0])


def read_moments(path):
    vals = []
    with open(path) as fh:
        for line in fh.readlines()[2:]:
            p = line.split()
            if len(p) >= 2:
                vals.append(complex(float(p[0]), float(p[1])))
    return np.array(vals)


def main(out="out/performance"):
    ok = True
    # part 1: exactness
    m1 = os.path.abspath(os.path.join(out, "model_L64"))
    square_lattice(64, m1)
    dirs = {}
    for name, op, ck in (("csr", "csr", False), ("atom_block", "atom_block", False), ("checkpoint", "atom_block", True)):
        d = os.path.abspath(os.path.join(out, f"L64_{name}"))
        run(runfile(m1, d, moments=256, trace="exact", states=4096, operator=op, checkpoint=ck), d + ".json")
        dirs[name] = d
    a = read_moments(moments_file(dirs["csr"]))
    b = read_moments(moments_file(dirs["atom_block"]))
    rel = float(np.max(np.abs(a - b)) / np.max(np.abs(a)))
    same = filecmp.cmp(moments_file(dirs["atom_block"]), moments_file(dirs["checkpoint"]), shallow=False)
    print(f"L = 64, exact trace, M = 256: max |mu_csr - mu_atom_block| / max |mu| = {rel:.2e}")
    print(f"checkpointed atom_block run byte-identical to the plain one: {same}")
    ok &= rel < 1e-12 and same
    # part 2: timing
    m2 = os.path.abspath(os.path.join(out, "model_L1024"))
    square_lattice(1024, m2)
    for name in ("csr", "atom_block"):
        d = os.path.abspath(os.path.join(out, f"L1024_{name}"))
        t = run(runfile(m2, d, moments=512, trace="random", states=2, operator=name), d + ".json")
        print(f"L = 1024 ({1024 * 1024} sites), 2 random vectors, M = 512, {name:10s}: {t:6.1f} s "
              f"(OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', 'default')})")
    print("checks:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:2]))
