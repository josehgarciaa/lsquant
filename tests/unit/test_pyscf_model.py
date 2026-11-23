"""PySCF periodic Gaussian connector.

Without PySCF: the reorder of two-component matrices, the minimal-image
inverse Fourier transform and the reader on a synthetic file. With PySCF (the
test is skipped otherwise): a two-component calculation of a hydrogen chain
written and read back, whose eigenvalues must match PySCF on the mesh.
"""
import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

from lsquant.models.pyscf_lcao import (ghf_to_interleaved, load_pyscf_model,  # noqa: E402
                                  mesh_fractions, minimal_image_blocks)


def _bloch(r_int, blocks, k_frac):
    return np.einsum("kr,rij->kij", np.exp(2j * np.pi * np.asarray(k_frac) @ np.asarray(r_int).T), blocks)


def test_ghf_to_interleaved():
    perm = ghf_to_interleaved(3)
    assert perm.tolist() == [0, 3, 1, 4, 2, 5]
    m = np.arange(36).reshape(6, 6)
    out = m[perm][:, perm]
    assert out[1, 1] == m[3, 3] and out[0, 1] == m[0, 3]


def test_minimal_images_recover_short_range_blocks():
    rng = np.random.default_rng(1)
    uc = np.diag([3.0, 3.5, 20.0])
    tau = np.array([[0.2, 0.1, 0.0], [1.6, 1.9, 0.3]])
    w = len(tau)
    cells = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (-1, 2, 0)]
    blocks = {}
    for r in cells:
        a = rng.normal(size=(w, w)) + 1j * rng.normal(size=(w, w))
        blocks[r] = a if r != (0, 0, 0) else a + a.conj().T
    for r in list(blocks):
        if r != (0, 0, 0):
            blocks[tuple(-x for x in r)] = blocks[r].conj().T
    r_ref = np.array(sorted(blocks))
    m_ref = np.array([blocks[tuple(r)] for r in r_ref])
    mesh = (6, 6, 1)
    k = mesh_fractions(mesh)
    r_out, m_out = minimal_image_blocks(_bloch(r_ref, m_ref, k), k, mesh, uc, tau)
    # Exact on the mesh ...
    assert np.allclose(_bloch(r_out, m_out, k), _bloch(r_ref, m_ref, k))
    # ... and, since the blocks fit in the supercell, exact everywhere.
    k_off = rng.random((7, 3)) * [1, 1, 0]
    assert np.allclose(_bloch(r_out, m_out, k_off), _bloch(r_ref, m_ref, k_off))


def test_reader_builds_spinor_overlap_and_midpoint_connection(tmp_path):
    rng = np.random.default_rng(2)
    uc = np.diag([3.0, 3.0, 15.0])
    tau_orb = np.array([[0.0, 0.0, 0.0], [1.5, 1.0, 0.2]])
    r_int = np.array([[-1, 0, 0], [0, 0, 0], [1, 0, 0]])
    n = len(tau_orb)
    s1 = rng.normal(scale=0.1, size=(n, n))
    s_orb = np.array([s1.T, np.eye(n) + 0.05 * (np.ones((n, n)) - np.eye(n)), s1])
    # A Hermitian connection B (B(-R)^T = B(R)) and the position matrix it implies.
    b1 = rng.normal(size=(3, n, n))
    b0 = rng.normal(size=(3, n, n))
    B = np.array([np.swapaxes(b1, 1, 2), b0 + np.swapaxes(b0, 1, 2), b1]).transpose(1, 0, 2, 3)
    mid = (tau_orb[None, :, None, :] + tau_orb[None, None, :, :] + (r_int @ uc)[:, None, None, :]) / 2
    Q = B + np.moveaxis(mid, -1, 0) * s_orb[None]
    h = rng.normal(size=(3, 2 * n, 2 * n)) + 1j * rng.normal(size=(3, 2 * n, 2 * n))
    h = np.array([h[2].conj().T, h[1] + h[1].conj().T, h[2]])
    path = tmp_path / "chain.npz"
    np.savez(path, uc=uc, tau=np.repeat(tau_orb, 2, axis=0), r_int=r_int, hoppings=h,
             overlap_orb=s_orb, position_orb=Q, fermi_energy=-1.0, mesh=np.array([4, 1, 1]),
             edge_ratio=0.0, nelectron=2, basis="custom")
    model = load_pyscf_model(path)
    assert model.overlaps.shape == (3, 2 * n, 2 * n)
    assert np.allclose(model.overlaps[:, 0::2, 0::2], s_orb)
    assert np.allclose(model.overlaps[:, 0::2, 1::2], 0)
    assert np.allclose(model.dipoles[:, :, 0::2, 0::2], B)
    assert np.allclose(model.dipoles[:, :, 1::2, 1::2], B)
    assert model.metadata["source"] == "pyscf"


def test_hydrogen_chain_round_trip(tmp_path):
    pytest.importorskip("pyscf")
    from pyscf.pbc import gto, scf
    from lsquant.models.pyscf_lcao import BOHR, HARTREE, write_pyscf_model

    cell = gto.Cell()
    cell.a = np.diag([1.6, 8.0, 8.0])
    cell.atom = [("H", (0.0, 0.0, 0.0)), ("H", (0.75, 0.1, 0.0))]
    cell.basis = "sto-3g"
    cell.verbose = 0
    cell.build()
    kpts = cell.make_kpts([6, 1, 1])
    mf = scf.KGHF(cell, kpts).density_fit()
    mf.kernel()
    path = tmp_path / "h2.npz"
    info = write_pyscf_model(mf, path)
    model = load_pyscf_model(path)
    k_frac = mesh_fractions([6, 1, 1])
    hk = _bloch(model.r_int, model.hoppings, k_frac)
    sk = _bloch(model.r_int, model.overlaps, k_frac)
    import scipy.linalg as sla
    ours = np.sort([sla.eigh(a, b, eigvals_only=True) for a, b in zip(hk, sk)], axis=1)
    ref = np.sort(np.array(mf.mo_energy) * HARTREE - info["fermi_energy"], axis=1)
    assert np.allclose(ours, ref, atol=1e-6)
    assert np.allclose(model.uc, cell.lattice_vectors() * BOHR)


def test_band_veff_equals_two_component_veff_on_the_scf_mesh():
    """The off-mesh J + V_xc route reproduces PySCF's two-component veff at the SCF k points."""
    pytest.importorskip("pyscf")
    from pyscf.pbc import gto, dft
    from lsquant.models.pyscf_lcao import _band_veff

    cell = gto.Cell()
    cell.a = np.diag([1.6, 8.0, 8.0])
    cell.atom = [("H", (0.0, 0.0, 0.0)), ("H", (0.75, 0.1, 0.0))]
    cell.basis = "sto-3g"
    cell.verbose = 0
    cell.build()
    mf = dft.KGKS(cell, cell.make_kpts([4, 1, 1])).density_fit()
    mf.xc = "lda"
    mf.kernel()
    dm = mf.make_rdm1()
    assert np.allclose(_band_veff(mf, dm, mf.kpts), mf.get_veff(cell, dm), atol=1e-8)
