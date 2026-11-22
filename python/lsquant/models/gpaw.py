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
    if world.rank == 0:
        np.savez_compressed(
            path, format=FORMAT, kpts=kd.bzk_kc, mesh=np.asarray(kd.N_c), H=H, S=S, A=A, soc=bool(soc),
            lattice=atoms.cell.array, positions=atoms.positions, symbols=np.array(atoms.get_chemical_symbols()),
            orbitals_per_atom=np.array([s.nao for s in wfs.setups]), fermi_energy=fermi,
            **extra)


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


def load_gpaw_model(path, *, soc=True, fermi_energy=None, tolerance=1e-10):
    """Return the LSQuant model of a file written by :func:`export_gpaw_lcao`.

    Non-orthogonal (H and S of the LCAO basis), in eV and Angstrom. With
    ``soc`` the model is spinful: ``H(R) tensor 1 + A(R)`` in interleaved
    spinors; without it, spinless. ``fermi_energy`` (eV) defaults to the Fermi
    level of GPAW and is removed as ``H - E_F S``. Blocks whose largest
    element is below ``tolerance`` (eV, or dimensionless for S) are dropped.
    """
    from .lcao import assemble_blocks, build_lcao_model

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
            "gpaw_fermi_level_eV": float(data["fermi_energy"])}
    return build_lcao_model(data["lattice"], data["positions"], norb, r_h, H, S, spinful=with_soc,
                            seed_name="gpaw", metadata=meta, fermi_energy=ef, overlap_r_int=r_s)


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
