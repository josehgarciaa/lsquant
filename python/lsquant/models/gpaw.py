"""GPAW LCAO calculations with spin-orbit coupling as LSQuant models.

GPAW computes spin-orbit coupling after the self-consistent cycle
(``gpaw.spinorbit.soc_eigenstates``): it adds, in the basis of the scalar
relativistic Bloch states, the on-site PAW term

    H_soc = sum_a sum_ij |p~_i^a> (dV^a/dr / r) L . sigma <p~_j^a|,

with the matrices ``dVL_vii`` of ``gpaw.spinorbit.soc``. In the LCAO mode
the same term has a matrix in the atomic orbitals,

    A^{ss'}_{MN}(k) = sum_a sum_ij conj(P^a_{Mi}(k)) H^{a,ss'}_{ij} P^a_{Nj}(k),

with ``P^a_{Mi}(k) = <p~_i^a | Phi_{Mk}>`` (``wfs.P_aqMi``) and
``H^{a,ss'} = sum_v dVL^a_v sigma_v[s, s']``. The spinor LCAO Hamiltonian is
``H(k) tensor 1 + A(k)`` with the overlap ``S(k) tensor 1``. Using every
LCAO band, its eigenvalues are those of ``soc_eigenstates``: the two
diagonalise the same operator in two complete bases of the same space.

Two steps:

1. :func:`export_gpaw_lcao` runs in the Python of GPAW after the
   self-consistent cycle. It writes H(k), S(k) and A(k) on the k mesh of the
   calculation, the structure and the Fermi level to an ``.npz`` file.
2. :func:`load_gpaw_model` reads that file in the Python of LSQuant (GPAW is
   not needed) and Fourier transforms the matrices to H(R) and S(R).

Conventions, read in the sources of GPAW 26.7:

- The LCAO matrices of GPAW (``calculate_hamiltonian_matrix``, ``S_qMM``)
  are ``X_{MN}(k) = sum_R exp(-2 pi i k.R) <M, 0|X|N, R>``: the complex
  conjugate of the LSQuant convention ``exp(+2 pi i k.R)``. Checked on the
  overlap, which is short ranged only with this sign; the export stores
  ``conj(H)`` and ``conj(S)``. In that convention the spin-orbit matrix
  is ``A`` above as written, with ``P_aqMi`` as GPAW stores it: then
  ``conj(H) tensor 1 + A`` (overlap ``conj(S) tensor 1``) has the eigenvalues of
  ``soc_eigenstates``, and ``H tensor 1 + A`` does not. The inverse transform on an ``n1 x n2 x n3`` mesh gives
  each matrix element modulo ``n_c`` cells; every element of an atom pair is
  placed at the image ``R`` that minimises ``|R + tau_j - tau_i|`` (images at
  the same distance share the element equally), the atom-gauge range of the
  Bloch sums. The mesh must be larger than twice the range of the orbitals.
- ``calculate_hamiltonian_matrix`` fills one triangle; ``tri2full`` copies
  the other. Hartree is converted with ``ase.units.Ha``.
- Spin blocks ``A^{ss'}`` are stored ``[[uu, ud], [du, dd]]`` per atom pair
  and become interleaved spinors in :func:`~lsquant.models.lcao.assemble_blocks`.
- Orbitals of one atom are contiguous, atoms in the order of the ``Atoms``
  object (``basis_functions.M_a``); every orbital sits on its atom.

Validation (PdSe2 monolayer, GPAW 26.7, dzp, PBE, 12 x 12 mesh, 212
spinor orbitals): ``conj(H) tensor 1 + A`` reproduces ``soc_eigenstates`` to
1e-10 eV at every k of the mesh; the model of :func:`load_gpaw_model`
reproduces them to 6e-7 eV on the mesh and, on Gamma-X-S-Y-Gamma off the
mesh, a non-self-consistent GPAW run with ``soc_eigenstates`` to 6e-10 eV
within 4 eV of E_F. A restart file exported with the script gives the same
matrices as the export inside the run. ``tests/unit/test_gpaw_model.py``
checks the transform on synthetic files.

Warning: the spin-orbit term is not self-consistent; it is the GPAW
approximation. The velocity is ``bare`` or ``berry_connection``: the export
has no position matrix between orbitals, so no ``covariant``.
"""
from __future__ import annotations

import argparse
import itertools
import os
import sys

import numpy as np

FORMAT = "lsquant-gpaw-lcao-1"


def export_gpaw_lcao(calc, path, *, soc=True, check=False):
    """Write the LCAO matrices of a converged GPAW calculation to ``path``.

    Run it in the Python of GPAW, after ``atoms.get_potential_energy()``,
    with the old calculator (``gpaw.calculator.GPAW``) in LCAO mode, without
    k-point symmetry (``symmetry='off'``, so the mesh covers the whole
    Brillouin zone), spin paired, and without domain decomposition
    (``parallel={'domain': 1}``); k points may be distributed over MPI ranks.
    ``soc`` adds A(k); ``check`` also stores the eigenvalues of
    ``soc_eigenstates`` (all LCAO bands are needed: ``nbands='nao'``), which
    :func:`load_gpaw_model` does not use.
    """
    from ase.units import Ha
    from gpaw.mpi import world
    from gpaw.utilities.tools import tri2full

    wfs = calc.wfs
    kd = wfs.kd
    if wfs.mode != "lcao":
        raise ValueError("export_gpaw_lcao needs an LCAO calculation (mode='lcao')")
    if kd.nibzkpts != kd.nbzkpts:
        raise ValueError("the k mesh must cover the whole Brillouin zone: use symmetry='off'")
    if wfs.nspins != 1 or not wfs.collinear:
        raise ValueError("export_gpaw_lcao supports spin-paired calculations only")
    if wfs.gd.comm.size != 1:
        raise ValueError("run without domain decomposition: parallel={'domain': 1}")
    if wfs.S_qMM is None:
        # Restart file: the basis functions and projections are not set up yet
        # (as in gpaw.lcao.tools.get_lcao_hamiltonian).
        wfs.set_positions(calc.spos_ac)
    nao, nk = wfs.setups.nao, kd.nbzkpts
    H = np.zeros((nk, nao, nao), complex)
    S = np.zeros((nk, nao, nao), complex)
    A = np.zeros((nk, 2, 2, nao, nao), complex)
    dvl = {}
    if soc:
        from gpaw.spinorbit import soc as soc_matrices
        dvl = {a: soc_matrices(wfs.setups[a], calc.hamiltonian.xc, D_sp)
               for a, D_sp in calc.density.D_asp.items()}
    for kpt in wfs.kpt_u:
        h = wfs.eigensolver.calculate_hamiltonian_matrix(calc.hamiltonian, wfs, kpt)
        tri2full(h)
        s = wfs.S_qMM[kpt.q].copy()
        tri2full(s)
        # GPAW sums exp(-i k.R); LSQuant exp(+i k.R): store the conjugates.
        H[kpt.k], S[kpt.k] = np.conj(h) * Ha, np.conj(s)
        for a, d in dvl.items():
            pauli = np.array([[d[2], d[0] - 1j * d[1]], [d[0] + 1j * d[1], -d[2]]]) * Ha
            P = wfs.P_aqMi[a][kpt.q]
            for s1, s2 in itertools.product(range(2), range(2)):
                A[kpt.k, s1, s2] += P.conj() @ pauli[s1, s2] @ P.T
    for X in (H, S, A):
        kd.comm.sum(X, 0)
    extra = {}
    if check and soc:
        from gpaw.spinorbit import soc_eigenstates
        extra["soc_eigenvalues"] = soc_eigenstates(calc).eigenvalues()
    atoms, fermi = calc.get_atoms(), calc.get_fermi_level()     # on every rank
    extra.update(_radial_data(wfs.setups))
    if world.rank == 0:
        np.savez_compressed(
            path, format=FORMAT, kpts=kd.bzk_kc, mesh=np.asarray(kd.N_c), H=H, S=S, A=A, soc=bool(soc),
            lattice=atoms.cell.array, positions=atoms.positions, symbols=np.array(atoms.get_chemical_symbols()),
            orbitals_per_atom=np.array([s.nao for s in wfs.setups]), fermi_energy=fermi,
            **extra)


RADIAL_POINTS = 3001


def _radial_data(setups):
    """Tabulate what the position matrix needs, per atom, in Bohr.

    For the basis functions and the projectors of each atom: ``l`` and the
    cut-off of every radial function and its values ``f(r)`` on
    ``RADIAL_POINTS`` uniform points from 0 to the cut-off, where GPAW's
    function is ``f(r) r^l Y_L``, ``r^l Y_L`` the real solid harmonics of
    ``gpaw.sphere.spherical_harmonics.YL`` (also stored, up to l = 3). Plus
    the PAW corrections ``dO_ii`` and ``Delta_iiL`` for L = 1, 2, 3
    (y, z, x).
    """
    from gpaw.sphere.spherical_harmonics import YL

    out = {}
    for a, setup in enumerate(setups):
        for name, splines in (("basis", setup.basis_functions_J), ("proj", setup.pt_j)):
            ls = [sp.get_angular_momentum_number() for sp in splines]
            rcs = [sp.get_cutoff() for sp in splines]
            out[f"{name}_l_{a}"] = np.array(ls, int)
            out[f"{name}_rc_{a}"] = np.array(rcs, float)
            out[f"{name}_f_{a}"] = np.array([sp.map(np.linspace(0.0, rc, RADIAL_POINTS))
                                            for sp, rc in zip(splines, rcs)])
        out[f"dO_{a}"] = np.asarray(setup.dO_ii, float)
        out[f"Delta_{a}"] = np.asarray(setup.Delta_iiL, float)[:, :, 1:4]
    from ase.units import Bohr
    out["bohr"] = float(Bohr)
    rows = [(L, c, *n) for L in range(16) for c, n in YL[L]]
    out["YL"] = np.array(rows, float)
    return out


def read_gpaw_lcao(path):
    """Return the content of a file of :func:`export_gpaw_lcao` as a dict."""
    with np.load(path, allow_pickle=False) as f:
        data = {k: f[k] for k in f.files}
    if str(data.get("format")) != FORMAT:
        raise ValueError(f"{path}: not an LSQuant GPAW export (format {data.get('format')!r})")
    return data


def real_space_blocks(kpts, mesh, matrices, lattice, positions, norb):
    """Inverse Bloch transform of ``matrices[k, ..., M, N]`` to atom-pair blocks.

    Returns ``{(R, i, j): block}`` with integer ``R`` and the trailing axes of
    each block those of the atom pair. Each element is placed at the image of
    ``R`` (modulo ``mesh``) that minimises ``|R + tau_j - tau_i|``; images
    within 1e-6 Angstrom of the minimum share it equally.
    """
    kpts = np.asarray(kpts, float)
    mesh = np.asarray(mesh, int)
    lattice = np.asarray(lattice, float)
    positions = np.asarray(positions, float)
    starts = np.concatenate([[0], np.cumsum(norb)])
    if len(kpts) != int(np.prod(mesh)):
        raise ValueError(f"{len(kpts)} k points do not fill the {tuple(mesh)} mesh")
    shifts = [range(-1, 2) if n > 1 else range(1) for n in mesh]
    blocks = {}
    for R0 in itertools.product(*[range(n) for n in mesh]):
        phase = np.exp(-2j * np.pi * kpts @ np.array(R0)) / len(kpts)
        full = np.tensordot(phase, matrices, axes=(0, 0))
        for i, j in itertools.product(range(len(norb)), repeat=2):
            images = [np.array(R0) + mesh * np.array(m) for m in itertools.product(*shifts)]
            dist = np.array([np.linalg.norm((R @ lattice) + positions[j] - positions[i]) for R in images])
            near = [R for R, d in zip(images, dist) if d < dist.min() + 1e-6]
            block = full[..., starts[i]:starts[i + 1], starts[j]:starts[j + 1]] / len(near)
            for R in near:
                blocks[(tuple(int(x) for x in R), i, j)] = block
    return blocks


class _Functions:
    """The radial functions of one atom times GPAW's solid harmonics (Bohr)."""

    def __init__(self, data, atom, kind, yl):
        self.l = data[f"{kind}_l_{atom}"].astype(int)
        self.rc = data[f"{kind}_rc_{atom}"].astype(float)
        self.f = data[f"{kind}_f_{atom}"].astype(float)
        self.radius = float(self.rc.max())
        self.size = int(np.sum(2 * self.l + 1))
        self.yl = yl

    def groups(self, ratio=1.5):
        """Split the radial functions by cut-off: ``[(rows, functions), ...]``.

        A function joins the current group while its cut-off is within
        ``ratio`` of the group's smallest; each group is integrated on its
        own box, so short functions are not resolved on the grid of long ones.
        """
        order = np.argsort(self.rc, kind="stable")
        bounds = np.concatenate([[0], np.cumsum(2 * self.l + 1)])
        out, current = [], []
        for J in order:
            if current and self.rc[J] > ratio * self.rc[current[0]]:
                out.append(current)
                current = []
            current.append(J)
        out.append(current)
        result = []
        for js in out:
            js = sorted(js)
            rows = np.concatenate([np.arange(bounds[J], bounds[J + 1]) for J in js])
            sub = object.__new__(_Functions)
            sub.l, sub.rc, sub.f, sub.yl = self.l[js], self.rc[js], self.f[js], self.yl
            sub.radius, sub.size = float(sub.rc.max()), int(np.sum(2 * sub.l + 1))
            result.append((rows, sub))
        return result

    def onsite(self, other=None):
        """``<M|N>`` and ``<M|r|N>`` for functions on one centre (the origin), exactly.

        ``other`` (default: these functions) gives the columns.

        Radial integrals on the tabulated grid times angular integrals of the
        harmonics (Gauss-Legendre in cos(theta), uniform in phi: exact for
        these polynomial degrees).
        """
        ct, wt = np.polynomial.legendre.leggauss(24)
        phi = np.linspace(0, 2 * np.pi, 48, endpoint=False)
        st = np.sqrt(1 - ct**2)
        n = np.stack([np.outer(st, np.cos(phi)), np.outer(st, np.sin(phi)), np.outer(ct, np.ones_like(phi))],
                     axis=-1).reshape(-1, 3)
        w = np.outer(wt, np.full(len(phi), 2 * np.pi / len(phi))).ravel()
        other = self if other is None else other

        def harmonics(fs):
            bounds = np.concatenate([[0], np.cumsum(2 * fs.l + 1)])
            ylm = np.zeros((fs.size, len(n)))
            for J, l in enumerate(fs.l):
                for k, L in enumerate(range(l * l, (l + 1) ** 2)):
                    for c, nx, ny, nz in fs.yl[L]:
                        ylm[bounds[J] + k] += c * n[:, 0] ** nx * n[:, 1] ** ny * n[:, 2] ** nz
            return ylm, bounds

        y1, bounds1 = harmonics(self)
        y2, bounds2 = harmonics(other)
        ang0 = (y1 * w) @ y2.T
        ang1 = [(y1 * (w * n[:, v])) @ y2.T for v in range(3)]
        S = np.zeros((self.size, other.size))
        A = np.zeros((3, self.size, other.size))
        for J1, J2 in itertools.product(range(len(self.l)), range(len(other.l))):
            rc = min(self.rc[J1], other.rc[J2])
            r = np.linspace(0.0, rc, 4 * len(self.f[J1]))
            f1 = np.interp(r, np.linspace(0, self.rc[J1], len(self.f[J1])), self.f[J1])
            f2 = np.interp(r, np.linspace(0, other.rc[J2], len(other.f[J2])), other.f[J2])
            base = f1 * f2 * r ** (self.l[J1] + other.l[J2] + 2)
            b1, b2 = slice(bounds1[J1], bounds1[J1 + 1]), slice(bounds2[J2], bounds2[J2 + 1])
            S[b1, b2] = np.trapezoid(base, r) * ang0[b1, b2]
            i1 = np.trapezoid(base * r, r)
            for v in range(3):
                A[v, b1, b2] = i1 * ang1[v][b1, b2]
        return S, A

    def __call__(self, points):
        """Values at ``points`` (P, 3), relative to the atom: shape (size, P)."""
        r = np.linalg.norm(points, axis=1)
        out = np.empty((self.size, len(points)))
        row = 0
        for l, rc, f in zip(self.l, self.rc, self.f):
            radial = np.interp(r, np.linspace(0.0, rc, len(f)), f, right=0.0)
            for L in range(l * l, (l + 1) ** 2):
                harmonic = np.zeros(len(points))
                for c, nx, ny, nz in self.yl[L]:
                    harmonic += c * points[:, 0] ** nx * points[:, 1] ** ny * points[:, 2] ** nz
                out[row] = radial * harmonic
                row += 1
        return out


def _pair_integrals(fi, fj, d, order, dipole=True):
    """``<f_i, 0 | f_j, d>`` and ``<f_i, 0 | (u - d/2) | f_j, d>`` by Gauss-Legendre.

    Tensor nodes on the intersection of the two cut-off boxes; the
    functions vanish beyond their cut-offs. Bohr in, Bohr out.
    """
    lower = np.maximum(-fi.radius, d - fj.radius)
    upper = np.minimum(fi.radius, d + fj.radius)
    S = np.zeros((fi.size, fj.size))
    A = np.zeros((3, fi.size, fj.size))
    if np.any(lower >= upper):
        return S, A
    x, w = np.polynomial.legendre.leggauss(order)
    axes = [0.5 * (upper[c] + lower[c]) + 0.5 * (upper[c] - lower[c]) * x for c in range(3)]
    weights = [0.5 * (upper[c] - lower[c]) * w for c in range(3)]
    yz = np.array(np.meshgrid(axes[1], axes[2], indexing="ij")).reshape(2, -1).T
    wyz = (weights[1][:, None] * weights[2][None, :]).ravel()
    for xv, wx in zip(axes[0], weights[0]):
        points = np.column_stack((np.full(len(yz), xv), yz))
        a, b = fi(points), fj(points - d)
        wa = a * (wx * wyz)
        S += wa @ b.T
        if dipole:
            for c in range(3):
                A[c] += (wa * (points[:, c] - d[c] / 2)) @ b.T
    return S, A


def _grouped_integrals(fi, fj, d, order, dipole=True):
    """:func:`_pair_integrals` over every pair of cut-off groups (see ``_Functions.groups``)."""
    S = np.zeros((fi.size, fj.size))
    A = np.zeros((3, fi.size, fj.size))
    if not np.any(d):
        return fi.onsite(fj)
    for (ri, gi), (rj, gj) in itertools.product(fi.groups(), fj.groups()):
        if np.linalg.norm(d) >= gi.radius + gj.radius:
            continue
        s, a = _pair_integrals(gi, gj, d, order, dipole)
        S[np.ix_(ri, rj)], A[:, ri[:, None], rj[None, :]] = s, a
    return S, A


def gpaw_connection(data, r_int, *, order=48):
    """Midpoint connection ``A_a(R)`` and overlap ``S(R)`` of the PAW LCAO basis.

    For every cell of ``r_int`` and atom pair, in Angstrom:

        A_a,MN(R) = <M, 0| r_a - m_a |N, R>,   m = (tau_M + tau_N + R) / 2,

    with the PAW position operator: the pseudo integral of the basis
    functions (Gauss-Legendre quadrature of :func:`_pair_integrals` on the
    box of each pair of cut-off groups; radial and angular integrals for the
    blocks of one atom) plus

        sum_{a', R'} <M|p_i^{a'}, R'> [(c - m) dO_ij + Delta_r,ij] <p_j^{a'}, R'|N, R>,

    where ``c`` is the centre of the projector atom and
    ``Delta_r = sqrt(4 pi / 3) Delta_iiL`` its dipole correction. The same
    construction without the dipole gives the PAW overlap, which is compared
    with GPAW's S(R) by the caller. Returns ``(S, A)`` with shapes
    ``(N_R, N, N)`` and ``(3, N_R, N, N)`` in the orbital order of GPAW.
    """
    bohr = float(data["bohr"])
    yl = {}
    for L, c, nx, ny, nz in data["YL"]:
        yl.setdefault(int(L), []).append((c, int(nx), int(ny), int(nz)))
    lattice = np.asarray(data["lattice"], float) / bohr
    tau = np.asarray(data["positions"], float) / bohr
    natoms = len(tau)
    basis = [_Functions(data, a, "basis", yl) for a in range(natoms)]
    proj = [_Functions(data, a, "proj", yl) for a in range(natoms)]
    norb = np.array([b.size for b in basis])
    if not np.array_equal(norb, data["orbitals_per_atom"]):
        raise ValueError("the radial functions do not match the orbitals per atom of the export")
    starts = np.concatenate([[0], np.cumsum(norb)])
    dO = [np.asarray(data[f"dO_{a}"]) for a in range(natoms)]
    # Delta_iiL is stored for L = 1, 2, 3 = (y, z, x).
    dr = [np.sqrt(4 * np.pi / 3) * np.asarray(data[f"Delta_{a}"])[:, :, [2, 0, 1]].transpose(2, 0, 1)
          for a in range(natoms)]
    # Projector overlaps <M, 0 | p^a, R'> for every image within reach.
    reach = max(b.radius for b in basis) + max(p.radius for p in proj)
    span = [int(np.ceil(reach / np.linalg.norm(v))) + 1 for v in lattice]
    images = [np.array(m) for m in itertools.product(*[range(-n, n + 1) for n in span])]
    overlap_p = {}
    for i, a in itertools.product(range(natoms), repeat=2):
        for Rp in images:
            d = Rp @ lattice + tau[a] - tau[i]
            if np.linalg.norm(d) < basis[i].radius + proj[a].radius:
                o, _ = _grouped_integrals(basis[i], proj[a], d, order, dipole=False)
                overlap_p[(i, a, tuple(Rp))] = o
    r_int = np.asarray(r_int, int)
    lookup = {tuple(R): n for n, R in enumerate(r_int)}
    nao = int(norb.sum())
    S = np.zeros((len(r_int), nao, nao))
    A = np.zeros((3, len(r_int), nao, nao))
    done = set()
    for n, R in enumerate(r_int):
        for i, j in itertools.product(range(natoms), repeat=2):
            if (n, i, j) in done:
                continue
            d = R @ lattice + tau[j] - tau[i]
            mid = (tau[i] + tau[j] + R @ lattice) / 2
            s_blk, a_blk = _grouped_integrals(basis[i], basis[j], d, order)
            for (ii, a, Rp), oi in overlap_p.items():
                if ii != i:
                    continue
                oj = overlap_p.get((j, a, tuple(np.array(Rp) - R)))
                if oj is None:
                    continue
                c = np.array(Rp) @ lattice + tau[a]
                s_blk += oi @ dO[a] @ oj.T
                for v in range(3):
                    a_blk[v] += oi @ ((c - mid)[v] * dO[a] + dr[a][v]) @ oj.T
            rows, cols = slice(starts[i], starts[i + 1]), slice(starts[j], starts[j + 1])
            S[n, rows, cols], A[:, n, rows, cols] = s_blk, a_blk
            done.add((n, i, j))
            # Translation: <N, 0|r - m|M, -R> = <M, 0|r - m|N, R>^T for real functions.
            partner = lookup.get(tuple(-R))
            if partner is not None and (partner, j, i) not in done:
                S[partner, cols, rows], A[:, partner, cols, rows] = s_blk.T, a_blk.transpose(0, 2, 1)
                done.add((partner, j, i))
    return S, A * bohr


def load_gpaw_model(path, *, soc=True, fermi_energy=None, tolerance=1e-10, dipoles=False,
                    quadrature_order=48, overlap_tolerance=1e-4, dipole_cache=None):
    """Return the LSQuant model of a file written by :func:`export_gpaw_lcao`.

    Non-orthogonal (H and S of the LCAO basis), in eV and Angstrom. With
    ``soc`` the model is spinful: ``H(R) tensor 1 + A(R)`` in interleaved
    spinors; without it, spinless. ``fermi_energy`` (eV) defaults to the Fermi
    level of GPAW and is removed as ``H - E_F S``. Blocks whose largest
    element is below ``tolerance`` (eV, or dimensionless for S) are dropped.

    ``dipoles=True`` adds the midpoint connection of the covariant velocity,
    built by :func:`gpaw_connection` from the radial functions in the file;
    its overlap must match GPAW's S(R) within ``overlap_tolerance`` (the
    quadrature error: 1e-5 at order 48 for PdSe2). ``dipole_cache`` names an
    ``.npz`` that stores the result for the same export, cells and order.
    """
    from .lcao import assemble_blocks, build_lcao_model, fold_opposite, opposite_pairs, seed_name_from_file
    from .siesta import spinor_lift

    data = read_gpaw_lcao(path)
    norb = data["orbitals_per_atom"].astype(int)
    args = (data["kpts"], data["mesh"])
    geometry = (data["lattice"], data["positions"], norb)
    with_soc = bool(soc)
    if with_soc and not bool(data["soc"]):
        raise ValueError(f"{path} was exported without spin-orbit coupling; use soc=False")
    h_blocks = real_space_blocks(*args, data["H"], *geometry)
    s_blocks = real_space_blocks(*args, data["S"], *geometry)
    a_blocks = real_space_blocks(*args, data["A"], *geometry) if with_soc else {}

    def pairs_h():
        for key, h in h_blocks.items():
            if with_soc:
                a = a_blocks[key]
                block = np.block([[h + a[0, 0], a[0, 1]], [a[1, 0], h + a[1, 1]]])
            else:
                block = h
            if np.abs(block).max() > tolerance:
                yield key[0], key[1], key[2], block

    def pairs_s():
        for key, s in s_blocks.items():
            if np.abs(s).max() > tolerance:
                yield key[0], key[1], key[2], s

    r_h, H = assemble_blocks(pairs_h(), norb, spinful=with_soc, name="H")
    r_s, S = assemble_blocks(pairs_s(), norb, spinful=False, name="S")
    ef = float(data["fermi_energy"]) if fermi_energy is None else float(fermi_energy)
    meta = {"source": "gpaw", "file": str(path), "spin_orbit": "non-self-consistent PAW on-site term"
            if with_soc else "none", "k_mesh": [int(n) for n in data["mesh"]],
            "gpaw_fermi_level_eV": float(data["fermi_energy"]),
            "fermi_energy_source": "file" if fermi_energy is None else "run file"}
    model = build_lcao_model(data["lattice"], data["positions"], norb, r_h, H, S, spinful=with_soc,
                             seed_name=seed_name_from_file(path), metadata=meta, fermi_energy=ef,
                             overlap_r_int=r_s)
    if dipoles:
        if "YL" not in data:
            raise ValueError(f"{path} holds no radial functions; export it again to use dipoles=True")
        import hashlib

        key = hashlib.sha256(np.ascontiguousarray(data["S"]).tobytes() + model.r_int.tobytes()
                             + str(quadrature_order).encode()).hexdigest()
        cached = None
        if dipole_cache is not None and os.path.isfile(dipole_cache):
            with np.load(dipole_cache) as f:
                if str(f["key"]) == key:
                    cached = (f["S"], f["A"])
        if cached is None:
            S_q, A_q = gpaw_connection(data, model.r_int, order=quadrature_order)
            if dipole_cache is not None:
                np.savez(dipole_cache, key=key, S=S_q, A=A_q)
        else:
            S_q, A_q = cached
        step = 2 if with_soc else 1
        S_model = model.overlaps[:, ::step, ::step].real
        error = float(np.abs(S_q - S_model).max())
        model.metadata.update(dipole_method=f"PAW position operator, Gauss-Legendre order {quadrature_order}",
                              dipole_overlap_check=error)
        if error > overlap_tolerance:
            raise ValueError(f"the quadrature overlap differs from GPAW's S(R) by {error:.1e} "
                             f"(> {overlap_tolerance:g}); raise quadrature_order")
        # Fold the connection Hermitian, as the SIESTA, OpenMX and ABACUS
        # readers do (audit item B7); the recorded residual is measured before
        # the fold, with its historical normalisation. Caches hold the
        # unfolded A: folding after reading keeps old caches valid.
        A_q, residual = fold_opposite(A_q, opposite_pairs(model.r_int), residual="absolute")
        model.metadata["dipole_hermiticity_residual"] = residual
        lifted = spinor_lift(A_q) if with_soc else A_q
        model.dipoles = np.asarray(lifted, complex)
    return model


def main(argv=None):
    """Export a GPAW LCAO restart file: ``python gpaw.py calc.gpw model.npz``.

    Runs in the Python of GPAW (also under ``mpirun ... gpaw python``);
    LSQuant need not be installed there.
    """
    ap = argparse.ArgumentParser(description="Export the LCAO matrices of a GPAW restart file for LSQuant.")
    ap.add_argument("gpw", help="restart file written by calc.write() after the SCF (LCAO, symmetry='off')")
    ap.add_argument("out", help="output .npz")
    ap.add_argument("--no-soc", action="store_true", help="do not add the spin-orbit matrix")
    ap.add_argument("--check", action="store_true", help="also store the eigenvalues of soc_eigenstates")
    args = ap.parse_args(argv)
    # Run as a script, this file's directory is on sys.path and its name
    # would hide the gpaw package.
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path[:] = [p for p in sys.path if os.path.abspath(p or os.curdir) != here]
    sys.modules.pop("gpaw", None)
    from gpaw.calculator import GPAW

    calc = GPAW(args.gpw, txt=None, parallel={"domain": 1})
    export_gpaw_lcao(calc, args.out, soc=not args.no_soc, check=args.check)
    return 0


if __name__ == "__main__":
    sys.exit(main())
