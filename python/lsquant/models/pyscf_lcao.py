"""PySCF periodic Gaussian models: export and reader.

A periodic PySCF calculation with spin-orbit coupling (``dft.KGKS(cell,
kpts).x2c1e()``, all-electron or with pseudopotentials) is turned into a
non-orthogonal LSQuant model: H(R), S(R) and the connection A(R) of the
covariant velocity.

Conventions of the file and of the returned model: eV, Angstrom; blocks
``M_ij(R) = <i,0|M|j,R>``; ``M(k) = sum_R exp(ik.(R+tau_j-tau_i)) M(R)``
with every basis function on its atom; interleaved spinors (function up,
function down). PySCF stores two-component matrices as [all alpha, all beta]
blocks; ``write_pyscf_model`` reorders them. Energies are relative to ``fermi_energy``
(H - E_F S).

Where each matrix comes from:

* S(R) and the position matrix Q(R) = <i,0|r|j,R> are analytic Gaussian
  integrals between the home cell and the cell at R. No quadrature and no
  grid: the connection A = Q - (tau_i + tau_j + R) S / 2 is exact in the
  basis (picture change of r in X2C neglected).
* H(R) is the inverse Fourier transform of the Kohn-Sham matrix on a
  Monkhorst-Pack mesh. The Kohn-Sham matrix is built from the converged
  density at the k points of that mesh, which may be denser than the SCF
  mesh (a non-self-consistent step, as in ``mf.get_bands``). Every block is
  assigned to its minimal images in the supercell of the mesh, with equal
  weights for degenerate images (Wannier90's ``use_ws_distance``), so H(k)
  is exact on the mesh and interpolated between.

Only ``write_pyscf_model`` imports PySCF; ``load_pyscf_model`` needs NumPy.
"""
from __future__ import annotations

import itertools
import os
import warnings

import numpy as np

from .units import BOHR_ANGSTROM as BOHR
from .units import HARTREE_EV as HARTREE


def ghf_to_interleaved(nao):
    """Index array that reorders [alpha..., beta...] to interleaved spinors."""
    perm = np.empty(2 * nao, int)
    perm[0::2] = np.arange(nao)
    perm[1::2] = np.arange(nao) + nao
    return perm


def mesh_fractions(mesh):
    """Fractional k points of a Gamma-centred mesh, PySCF's ``make_kpts`` order."""
    axes = [np.arange(n) / n for n in mesh]
    return np.array(list(itertools.product(*axes)))


def minimal_image_blocks(matrices_k, k_frac, mesh, uc, tau, tolerance=1e-5):
    """Real-space blocks of ``M(k)`` given on a full mesh, minimal images.

    ``matrices_k`` has shape (N_k, W, W) in the cell gauge
    ``M(k) = sum_R exp(2 pi i k_frac.R) M(R)``. The inverse transform gives
    the blocks summed over the supercell ``mesh``; each block (i, j, R) is
    placed at the images R + L*mesh, L in {-1, 0, 1}^3, with the shortest
    ``|tau_j + R + L*mesh - tau_i|``, split equally between images that are
    shorter than the shortest plus ``tolerance`` (Angstrom).

    Returns ``r_int`` (N_R, 3) and blocks (N_R, W, W). The blocks reproduce
    ``matrices_k`` exactly at every k of the mesh.
    """
    mesh = np.asarray(mesh, int)
    matrices_k = np.asarray(matrices_k)
    nk, w, _ = matrices_k.shape
    if nk != int(np.prod(mesh)):
        raise ValueError("matrices_k must cover the full mesh")
    base = mesh_fractions(mesh) * mesh  # integer R inside the supercell
    base = np.rint(base).astype(int)
    phase = np.exp(-2j * np.pi * np.asarray(k_frac) @ base.T)  # (N_k, N_R0)
    blocks0 = np.einsum("kr,kij->rij", phase, matrices_k) / nk
    shifts = np.array([L for L in itertools.product((-1, 0, 1), repeat=3)
                       if all(L[a] == 0 or mesh[a] > 1 for a in range(3))])
    tau = np.asarray(tau, float)
    d0 = tau[None, :, :] - tau[:, None, :]  # (W, W, 3): tau_j - tau_i
    out = {}
    for r0, block in zip(base, blocks0):
        cands = r0[None, :] + shifts * mesh[None, :]
        cart = cands @ uc
        dist = np.linalg.norm(d0[None] + cart[:, None, None, :], axis=-1)  # (N_L, W, W)
        dmin = dist.min(axis=0)
        chosen = dist <= dmin[None] + tolerance
        weight = chosen / chosen.sum(axis=0)[None]
        for c, wgt in zip(cands, weight):
            if not wgt.any():
                continue
            key = tuple(int(x) for x in c)
            out.setdefault(key, np.zeros((w, w), matrices_k.dtype))
            out[key] += wgt * block
    r_int = np.array(sorted(out), int)
    return r_int, np.array([out[tuple(r)] for r in r_int])


def edge_ratio(r_int, blocks, mesh):
    """Largest |M(R)| on the outer half of the supercell over the largest with R != 0."""
    mesh = np.asarray(mesh)
    size = np.abs(blocks).max(axis=(1, 2))
    nonzero = np.any(r_int != 0, axis=1)
    reach = np.max(np.abs(r_int) / np.maximum(mesh / 2, 1), axis=1)
    outer = reach >= 0.75
    if not outer.any() or not nonzero.any():
        return 0.0
    return float(size[outer].max() / size[nonzero].max())


def write_pyscf_model(mf, path, *, mesh=None, fermi_energy=None, overlap_cutoff=1e-9,
                      locality_warning=1e-3):
    """Write a converged periodic two-component PySCF calculation to ``path`` (.npz).

    ``mf`` is a converged ``pyscf.pbc.dft.KGKS`` (or ``KGHF``) object, with
    or without ``.x2c1e()``. ``mesh`` is the Monkhorst-Pack mesh on which the
    Kohn-Sham matrix is built (default: the SCF mesh, which must then be a
    full Gamma-centred mesh). ``fermi_energy`` (eV) defaults to the middle of
    the gap at the electron count of the cell; for a metal, pass it.
    """
    from pyscf import gto as mgto

    from .siesta import spinor_lift

    cell = mf.cell
    nao = cell.nao_nr()
    if mf.mo_coeff is None or not getattr(mf, "converged", True):
        raise ValueError("write_pyscf_model needs a converged calculation")
    perm = ghf_to_interleaved(nao)
    uc = np.asarray(cell.lattice_vectors()) * BOHR
    atom_pos = np.asarray(cell.atom_coords()) * BOHR
    ao_atom = np.array([int(lab.split()[0]) for lab in cell.ao_labels()])
    tau = np.repeat(atom_pos[ao_atom], 2, axis=0)
    dim = int(cell.dimension)
    if mesh is None:
        mesh = _infer_mesh(cell, mf.kpts)
    mesh = np.array([int(x) for x in mesh] + [1] * (3 - len(mesh)))
    if dim == 2 or np.allclose(mesh[2], 1):
        mesh[2] = 1
    k_frac = mesh_fractions(mesh)
    kpts = cell.get_abs_kpts(k_frac)
    dm = mf.make_rdm1()
    if len(kpts) == len(mf.kpts) and np.allclose(kpts, mf.kpts):
        fock = mf.get_fock(dm=dm)
    else:
        fock = np.asarray(mf.get_hcore(cell, kpts)) + _band_veff(mf, dm, kpts)
    fock = np.asarray(fock)[:, perm][:, :, perm] * HARTREE
    s_k = spinor_lift(np.asarray(cell.pbc_intor("int1e_ovlp", hermi=1, kpts=kpts)))
    nelec = int(round(cell.nelectron))
    if fermi_energy is None:
        import scipy.linalg as sla
        ev = np.sort(np.concatenate([sla.eigh(f, s, eigvals_only=True) for f, s in zip(fock, s_k)]))
        nocc = nelec * len(kpts)
        fermi_energy = float(0.5 * (ev[nocc - 1] + ev[nocc]))
    r_h, h_blocks = minimal_image_blocks(fock - fermi_energy * s_k, k_frac, mesh, uc, tau)
    ratio = edge_ratio(r_h, h_blocks, mesh)
    if ratio > locality_warning:
        warnings.warn(f"PySCF H(R) has not decayed at the edge of the {mesh.tolist()} supercell "
                      f"(edge/max = {ratio:.2e}); use a denser mesh", RuntimeWarning, stacklevel=2)
    # Exact S(R), Q(R) on every cell within the overlap range. Along a
    # periodic axis with one k point (a slab with vacuum computed as a 3D
    # cell) H(R) is summed over that axis, so S and Q are summed too; the
    # component of Q along that axis is then meaningless and set to zero.
    nimg = np.array([int(n) if a < dim else 0 for a, n in enumerate(cell.nimgs)])
    folded_axes = [a for a in range(3) if mesh[a] == 1]
    acc_s, acc_q = {}, {}
    for r in itertools.product(*[range(-n, n + 1) for n in nimg]):
        # Same basis, every atom moved by R (the coordinates live in _env).
        shifted = cell.copy(deep=True)
        ptr = shifted._atm[:, mgto.PTR_COORD]
        shift = np.asarray(r) @ cell.lattice_vectors()
        for a in range(cell.natm):
            shifted._env[ptr[a]:ptr[a] + 3] += shift
        s = mgto.intor_cross("int1e_ovlp", cell, shifted)
        if np.abs(s).max() <= overlap_cutoff:
            continue
        q = mgto.intor_cross("int1e_r", cell, shifted) * BOHR
        key = tuple(0 if a in folded_axes else int(r[a]) for a in range(3))
        acc_s[key] = acc_s.get(key, 0) + s
        acc_q[key] = acc_q.get(key, 0) + q
    for r in r_h:
        acc_s.setdefault(tuple(r), np.zeros((nao, nao)))
        acc_q.setdefault(tuple(r), np.zeros((3, nao, nao)))
    r_int = np.array(sorted(acc_s), int)
    s_r = np.array([acc_s[tuple(r)] for r in r_int])
    q_r = np.array([acc_q[tuple(r)] for r in r_int])  # (N_R, 3, n, n)
    for a in folded_axes:
        q_r[:, a] = 0.0
    index = {tuple(r): i for i, r in enumerate(r_int)}
    h_r = np.zeros((len(r_int), 2 * nao, 2 * nao), complex)
    for r, b in zip(r_h, h_blocks):
        h_r[index[tuple(r)]] = b
    np.savez(path, uc=uc, tau=tau, r_int=r_int, hoppings=h_r,
             overlap_orb=s_r, position_orb=np.moveaxis(q_r, 1, 0),
             fermi_energy=fermi_energy, mesh=mesh, edge_ratio=ratio,
             nelectron=nelec, basis=str(cell.basis if isinstance(cell.basis, str) else "custom"))
    return dict(fermi_energy=fermi_energy, edge_ratio=ratio, n_cells=len(r_int))


def _band_veff(mf, dm, kpts_band, magnetisation_tolerance=1e-4):
    """Two-component J + V_xc at k points off the SCF mesh, for a pure functional.

    PySCF's two-component XC code does not accept band k points (2.15), so
    J and V_xc are built from the total density with the one-component code
    and placed on both spin blocks. This is exact for a non-magnetic density
    (zero magnetisation everywhere, as in a time-reversal-symmetric crystal);
    a magnetic density is refused.
    """
    from pyscf.pbc.dft import numint as pni
    cell = mf.cell
    n = cell.nao_nr()
    dm = np.asarray(dm)
    rho = dm[:, :n, :n] + dm[:, n:, n:]
    # Magnetisation part of the home-cell density matrix D(R=0) = mean_k D(k);
    # with real basis functions only its symmetric real part enters m(r).
    d0 = dm.mean(axis=0)
    aa, bb, ab, ba = d0[:n, :n], d0[n:, n:], d0[:n, n:], d0[n:, :n]
    sym = lambda m: 0.5 * (m + m.T)
    mag = [sym((aa - bb).real), sym((ab + ba).real), sym((ab - ba).imag)]
    if max(np.abs(x).max() for x in mag) > magnetisation_tolerance:
        raise ValueError("write_pyscf_model(mesh=...) needs a non-magnetic density; "
                         "use the SCF mesh for a magnetic calculation")
    if mf._numint.libxc.is_hybrid_xc(mf.xc):
        raise ValueError("write_pyscf_model(mesh=...) supports pure functionals only")
    vj = mf.with_df.get_jk(rho, hermi=1, kpts=mf.kpts, kpts_band=kpts_band, with_k=False)[0]
    _, _, vxc = pni.KNumInt().nr_rks(cell, mf.grids, mf.xc, rho, kpts=mf.kpts, kpts_band=kpts_band)
    v1 = np.asarray(vj) + np.asarray(vxc)
    return np.array([np.kron(np.eye(2), v) for v in v1])


def _infer_mesh(cell, kpts):
    frac = cell.get_scaled_kpts(kpts)
    frac = np.mod(np.round(frac, 8), 1.0)
    mesh = []
    for a in range(3):
        vals = np.unique(np.round(frac[:, a], 6))
        mesh.append(len(vals))
    if int(np.prod(mesh)) != len(kpts):
        raise ValueError("the SCF k points are not a full Gamma-centred mesh; pass mesh=")
    return mesh


def load_pyscf_model(path, *, connection=True, bloch_memory_mb=4096.0, fermi_energy=None):
    """Return a :class:`~lsquant.models.siesta.SiestaModel` from a file of ``write_pyscf_model``.

    ``connection=True`` adds the midpoint connection of the covariant
    velocity, built from the analytic position matrix. An all-electron
    basis gives large blocks, so the Bloch-sum workspace defaults to 4 GiB.
    ``fermi_energy`` (eV) selects the energy origin: the file stores H with
    its own Fermi energy already subtracted, so an explicit value shifts H by
    ``-(x - E_F(file)) S``; the default keeps the file's origin.
    """
    from .siesta import SiestaModel, midpoint_connection, spinor_lift
    from .lcao import seed_name_from_file

    d = np.load(path)
    uc, tau, r_int = d["uc"], d["tau"], d["r_int"]
    s_orb = d["overlap_orb"]
    S = spinor_lift(s_orb)
    ef_file = float(d["fermi_energy"])
    if fermi_energy is None:
        origin, source = ef_file, "file"
        hoppings = d["hoppings"]
    else:
        origin, source = float(fermi_energy), "run file"
        hoppings = np.asarray(d["hoppings"], complex)
        shift = origin - ef_file
        if shift != 0.0:
            hoppings = hoppings - shift * S
    A = None
    if connection:
        q = d["position_orb"]  # (3, N_R, n, n)
        tau_orb = tau[0::2]
        a_orb = midpoint_connection(q, s_orb, tau_orb, r_int @ uc)
        A = np.asarray(spinor_lift(a_orb), complex)
    seed = seed_name_from_file(path)
    return SiestaModel(uc, tau, r_int, hoppings, S, dipoles=A, seed_name=seed,
                       bloch_memory_mb=float(bloch_memory_mb),
                       metadata=dict(source="pyscf", file=os.path.abspath(str(path)),
                                     printed_fermi_level_eV=ef_file,
                                     fermi_energy_subtracted_eV=origin,
                                     fermi_energy_source=source,
                                     mesh=d["mesh"].tolist(), edge_ratio=float(d["edge_ratio"]),
                                     valence_electrons=int(d["nelectron"])))
