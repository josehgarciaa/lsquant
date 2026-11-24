"""OpenMX reader: a non-orthogonal spin-orbit model from ``<System.Name>.scfout``.

OpenMX writes the scfout file with ``HS.fileout on``. This module reads the
format of OpenMX 3.9 (scfout version 3, read_scfout.c of the OpenMX source)
and returns the same model class as the SIESTA reader, so every SIESTA
executor, velocity mode and oracle applies.

Conventions and conversions:
- Lengths in the file are Bohr, energies Hartree. The model uses Angstrom and
  eV, with energies measured from the chemical potential: ``H - mu S``, as the
  SIESTA model stores ``E_F = 0``.
- A block of the file is ``<i, home cell | H | j, cell R>`` with atom ``i``
  (``ct_AN``), neighbour ``j = natn[ct_AN][h_AN]`` and lattice vector
  ``R = atv_ijk[ncn[ct_AN][h_AN]]``.
- Non-collinear spin with spin-orbit coupling (``SpinP_switch = 3``):
  ``H_uu = Hks[0] + i iHks[0]``, ``H_dd = Hks[1] + i iHks[1]``,
  ``H_ud = Hks[2] + i (Hks[3] + iHks[2])``, ``H_du(R) = H_ud(-R)^+``. The
  spinor index is reordered from OpenMX's (all up, all down) to LSQuant's
  interleaved (orbital, spin) order.
- ``OLPpo[a]`` is ``<i | (r - tau_i)_a | j, R>`` (row-centred dipole,
  integrated by OpenMX on its real-space grid). The connection of the
  covariant velocity is ``A = D - (tau_j + R - tau_i) S / 2`` (see
  :func:`~lsquant.models.siesta.midpoint_connection`).
"""
from __future__ import annotations

import os

import numpy as np

from .lcao import fold_opposite, opposite_pairs, seed_name_from_file, spinor_lift
from .siesta import SiestaModel
from .units import BOHR_ANGSTROM as BOHR
from .units import HARTREE_EV as HARTREE


class _Reader:
    def __init__(self, path):
        self.data = open(path, "rb").read()
        self.pos = 0

    def take(self, dtype, count):
        dt = np.dtype(dtype)
        out = np.frombuffer(self.data, dt, count, self.pos)
        self.pos += dt.itemsize * count
        return out


def read_scfout(path):
    """Read an OpenMX 3.9 scfout file into a dictionary of numpy arrays.

    Returned blocks are lists over ``(ct_AN, h_AN)`` in file order, with
    atoms numbered from 0. Units are those of the file (Bohr, Hartree).
    """
    r = _Reader(path)
    head = r.take("<i4", 6)
    if head[1] < 0 or head[1] > 4 * 3 + 3:
        raise ValueError(f"{path}: not a little-endian OpenMX scfout file")
    natom, spin, version = int(head[0]), int(head[1]) % 4, int(head[1]) // 4
    if version != 3:
        raise ValueError(f"{path}: scfout version {version}; only version 3 (OpenMX 3.9) is supported")
    ncell = int(head[5])
    order_max = int(r.take("<i4", 1)[0])
    atv = r.take("<f8", 4 * (ncell + 1)).reshape(ncell + 1, 4)[:, 1:]
    atv_ijk = r.take("<i4", 4 * (ncell + 1)).reshape(ncell + 1, 4)[:, 1:]
    norb = r.take("<i4", natom).astype(int)
    fnan = r.take("<i4", natom).astype(int)
    natn = [r.take("<i4", f + 1).astype(int) - 1 for f in fnan]
    ncn = [r.take("<i4", f + 1).astype(int) for f in fnan]
    tv = r.take("<f8", 12).reshape(3, 4)[:, 1:]
    rtv = r.take("<f8", 12).reshape(3, 4)[:, 1:]
    gxyz = r.take("<f8", 4 * natom).reshape(natom, 4)[:, 1:]

    def blocks():
        out = []
        for a in range(natom):
            for h in range(fnan[a] + 1):
                n1, n2 = norb[a], norb[natn[a][h]]
                out.append(r.take("<f8", n1 * n2).reshape(n1, n2))
        return out

    hks = [blocks() for _ in range(spin + 1)]
    ihks = [blocks() for _ in range(3)] if spin == 3 else None
    olp = blocks()
    olppo = [[blocks() for _ in range(order_max)] for _ in range(3)]
    for _ in range(3):
        blocks()  # OLPmo, momentum matrix, not used
    for _ in range(spin + 1):
        blocks()  # DM
    for _ in range(2):
        blocks()  # iDM
    r.take("<i4", 1)  # eigenvalue method flag
    tail = r.take("<f8", 10)
    pairs = [(a, natn[a][h], atv_ijk[ncn[a][h]]) for a in range(natom) for h in range(fnan[a] + 1)]
    return dict(natom=natom, spin=spin, norb=norb, pairs=pairs, tv=tv, rtv=rtv, gxyz=gxyz,
                atv=atv, hks=hks, ihks=ihks, olp=olp, olppo=olppo, chemical_potential=float(tail[0]),
                electronic_temperature=float(tail[1]), valence_electrons=float(tail[8]))


def load_openmx_model(path, *, dipoles=False):
    """Return a :class:`~lsquant.models.siesta.SiestaModel` of an OpenMX scfout file.

    Requires a non-collinear calculation with spin-orbit coupling
    (``scf.SpinPolarization NC``, ``scf.SpinOrbit.Coupling on``).
    ``dipoles=True`` builds the connection of the covariant velocity from
    the position matrix that OpenMX writes in the same file.
    """
    d = read_scfout(path)
    if d["spin"] != 3:
        raise ValueError("load_openmx_model requires a non-collinear spin-orbit calculation "
                         "(scf.SpinPolarization NC, scf.SpinOrbit.Coupling on)")
    natom, norb = d["natom"], d["norb"]
    offset = np.concatenate([[0], np.cumsum(norb)])
    n = int(offset[-1])
    uc = d["tv"] * BOHR
    atom_pos = d["gxyz"] * BOHR
    tau_orb = np.repeat(atom_pos, norb, axis=0)
    cells = sorted({tuple(int(x) for x in R) for _, _, R in d["pairs"]})
    index = {R: i for i, R in enumerate(cells)}
    nr = len(cells)
    uu = np.zeros((nr, n, n), complex)
    dd = np.zeros((nr, n, n), complex)
    ud = np.zeros((nr, n, n), complex)
    S = np.zeros((nr, n, n))
    D = np.zeros((3, nr, n, n))
    hks, ihks, mu = d["hks"], d["ihks"], d["chemical_potential"]
    for b, (a, j, R) in enumerate(d["pairs"]):
        k = index[tuple(int(x) for x in R)]
        rows, cols = slice(offset[a], offset[a + 1]), slice(offset[j], offset[j + 1])
        uu[k, rows, cols] = hks[0][b] + 1j * ihks[0][b]
        dd[k, rows, cols] = hks[1][b] + 1j * ihks[1][b]
        ud[k, rows, cols] = hks[2][b] + 1j * (hks[3][b] + ihks[2][b])
        S[k, rows, cols] = d["olp"][b]
        if dipoles:
            for x in range(3):
                D[x, k, rows, cols] = d["olppo"][x][0][b]
    minus = opposite_pairs(cells)
    du = np.conj(np.swapaxes(ud[minus], 1, 2))
    # Interleaved spinor order: basis function 2*o + s.
    H = np.zeros((nr, 2 * n, 2 * n), complex)
    H[:, 0::2, 0::2], H[:, 1::2, 1::2] = uu, dd
    H[:, 0::2, 1::2], H[:, 1::2, 0::2] = ud, du
    S_spin = spinor_lift(S)
    H = (H - mu * S_spin) * HARTREE
    # OpenMX integrates part of H and S on a real-space grid, so H(R) and
    # H(-R)^+ differ by round-off (about 1e-10 relative). Fold both Hermitian,
    # as the SIESTA reader does, and record the removed part.
    folded = {}
    for name, M in (("hamiltonian", H), ("overlap", S_spin)):
        sym, folded[name] = fold_opposite(M, minus, residual="mean")
        M[...] = sym
    r_int = np.array(cells, int)
    tau = np.repeat(tau_orb, 2, axis=0)
    A = None
    if dipoles:
        r_cart = r_int @ uc
        # A = D - (tau_j + R - tau_i) S / 2, then lifted to spinors.
        shift = (tau_orb[None, None, :, :] + r_cart[:, None, None, :] - tau_orb[None, :, None, :]) / 2
        A_orb = D * BOHR - np.moveaxis(shift, -1, 0) * S[None]
        # Exact integrals give A(-R)^+ = A(R); OpenMX's grid integration meets
        # it to about 1e-4 relative. Fold it and record the removed part.
        A, folded["connection"] = fold_opposite(A_orb, minus, residual="mean")
    seed = seed_name_from_file(path)
    model = SiestaModel(uc, tau, r_int, H, S_spin, dipoles=A, seed_name=seed,
                        metadata=dict(source="openmx", file=os.path.abspath(path),
                                      printed_fermi_level_eV=mu * HARTREE,
                                      valence_electrons=d["valence_electrons"],
                                      hermitian_fold_relative=folded))
    return model
