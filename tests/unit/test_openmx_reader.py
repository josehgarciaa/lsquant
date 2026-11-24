"""OpenMX scfout reader on a synthetic version-3 file written by the test.

The file has two atoms with one orbital each, a non-collinear spin-orbit
Hamiltonian split over Hks[0..3] and iHks[0..2] as OpenMX stores it, the
overlap, the row-centred position matrix and the chemical potential. The
reader must give H(k) - mu S(k) in eV with interleaved spinors, S(k) lifted
to spinors, and the covariant connection.
"""
import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant.models.openmx import BOHR, HARTREE, load_openmx_model  # noqa: E402
from lsquant.runfile import RunFile  # noqa: E402

CELLS = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]


def _fixture(rng):
    """Random Hermitian spinor blocks per lattice vector, orbital basis of 2."""
    def herm_pair():
        a = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        return {(0, 0, 0): 0.5 * (a + a.conj().T), (1, 0, 0): a * 0.3, (-1, 0, 0): (a * 0.3).conj().T}
    uu, dd = herm_pair(), herm_pair()
    b = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    ud = {(0, 0, 0): 0.2 * b, (1, 0, 0): 0.1 * b.T, (-1, 0, 0): 0.15 * b}
    s0 = np.array([[1.0, 0.1], [0.1, 1.0]])
    s1 = rng.normal(scale=0.05, size=(2, 2))
    S = {(0, 0, 0): s0, (1, 0, 0): s1, (-1, 0, 0): s1.T}
    D = {R: rng.normal(scale=0.2, size=(3, 2, 2)) for R in CELLS}
    return uu, dd, ud, S, D


def _write_scfout(path, uu, dd, ud, S, D, mu, gxyz, tv):
    natom, ncell = 2, 2
    # Neighbour list per atom: itself at R = 0 first (as OpenMX), then every
    # other (atom, cell) pair of the three cells.
    nb = {a: [(a, 0)] + [(j, c) for c in range(3) for j in range(natom) if (j, c) != (a, 0)]
          for a in range(natom)}
    with open(path, "wb") as fh:
        w = lambda arr, dt: fh.write(np.asarray(arr, dt).tobytes())
        w([natom, 3 + 4 * 3, natom, 0, 0, ncell], "<i4")
        w([1], "<i4")                                    # order_max
        atv = np.array([[0, 0, 0, 0], [0, tv[0, 0], tv[0, 1], tv[0, 2]], [0, -tv[0, 0], -tv[0, 1], -tv[0, 2]]])
        w(atv, "<f8")
        w([[0, 0, 0, 0], [0, 1, 0, 0], [0, -1, 0, 0]], "<i4")
        w([1, 1], "<i4")                                 # Total_NumOrbs
        w([len(nb[a]) - 1 for a in range(natom)], "<i4")  # FNAN
        for a in range(natom):
            w([j + 1 for j, _ in nb[a]], "<i4")          # natn (1-based)
        for a in range(natom):
            w([c for _, c in nb[a]], "<i4")              # ncn
        w(np.c_[np.zeros(3), tv], "<f8")
        w(np.c_[np.zeros(3), 2 * np.pi * np.linalg.inv(tv).T], "<f8")
        w(np.c_[np.zeros(2), gxyz], "<f8")
        def element(M, a, j, c):
            return M[CELLS[c]][a, j]
        def write_all(fn):
            for a in range(natom):
                for j, c in nb[a]:
                    w([fn(a, j, c)], "<f8")
        im_split = 0.4  # part of Im H_ud stored in Hks[3], the rest in iHks[2]
        write_all(lambda a, j, c: element(uu, a, j, c).real)
        write_all(lambda a, j, c: element(dd, a, j, c).real)
        write_all(lambda a, j, c: element(ud, a, j, c).real)
        write_all(lambda a, j, c: im_split * element(ud, a, j, c).imag)
        write_all(lambda a, j, c: element(uu, a, j, c).imag)
        write_all(lambda a, j, c: element(dd, a, j, c).imag)
        write_all(lambda a, j, c: (1 - im_split) * element(ud, a, j, c).imag)
        write_all(lambda a, j, c: element(S, a, j, c))
        for x in range(3):
            write_all(lambda a, j, c, x=x: D[CELLS[c]][x, a, j])
        for _ in range(3 + 4 + 2):                        # OLPmo, DM, iDM
            write_all(lambda a, j, c: 0.0)
        w([0], "<i4")                                    # eigenvalue method flag
        w([mu, 300.0, 0, 0, 0, 0, 0, 0, 2.0, 0.0], "<f8")
        w([0], "<i4")                                    # no input-file lines


def test_openmx_reader_hamiltonian_overlap_and_connection(tmp_path):
    rng = np.random.default_rng(4)
    uu, dd, ud, S, D = _fixture(rng)
    # H_du(R) = H_ud(-R)^+ makes the spinor Hamiltonian Hermitian.
    mu = -0.12
    tv = np.diag([7.0, 30.0, 30.0])
    gxyz = np.array([[0.5, 1.0, 0.0], [3.0, 2.0, 0.5]])
    path = str(tmp_path / "toy.scfout")
    _write_scfout(path, uu, dd, ud, S, D, mu, gxyz, tv)
    model = load_openmx_model(path, dipoles=True)
    assert model.num_wann == 4
    uc = tv * BOHR
    tau = np.repeat(gxyz * BOHR, 2, axis=0)
    k = rng.uniform(-1, 1, size=(3, 3)) @ (2 * np.pi * np.linalg.inv(uc).T)
    k[:, 1:] = 0

    def spinor(R):
        H = np.zeros((4, 4), complex)
        H[0::2, 0::2], H[1::2, 1::2] = uu[R], dd[R]
        H[0::2, 1::2] = ud[R]
        H[1::2, 0::2] = ud[tuple(-x for x in R)].conj().T
        return H
    expected_h, expected_s = [], []
    for kk in k:
        h = np.zeros((4, 4), complex)
        s = np.zeros((4, 4), complex)
        for R in CELLS:
            d = np.array(R) @ uc + tau[None, :, :] - tau[:, None, :]
            phase = np.exp(1j * d @ kk)
            sR = np.kron(S[R], np.eye(2))
            h += phase * (spinor(R) - mu * sR) * HARTREE
            s += phase * sR
        expected_h.append(h)
        expected_s.append(s)
    np.testing.assert_allclose(model.Hk(k), expected_h, atol=1e-10)
    np.testing.assert_allclose(model.Sk(k), expected_s, atol=1e-12)
    assert model.metadata["printed_fermi_level_eV"] == pytest.approx(mu * HARTREE)
    # A(R) = D - (tau_j + R - tau_i) S / 2, folded Hermitian, on the up-up block.
    i0 = [tuple(r) for r in model.r_int].index((1, 0, 0))
    im = [tuple(r) for r in model.r_int].index((-1, 0, 0))
    t = gxyz * BOHR
    raw = lambda R, x: D[R][x] * BOHR - ((t[None, :, x] + np.array(R) @ uc[:, x] - t[:, None, x]) / 2) * S[R]
    for x in range(3):
        want = 0.5 * (raw((1, 0, 0), x) + raw((-1, 0, 0), x).T)
        np.testing.assert_allclose(model.dipoles[x, i0][0::2, 0::2].real, want, atol=1e-12)
        np.testing.assert_allclose(model.dipoles[x, im], np.conj(model.dipoles[x, i0]).T, atol=1e-12)


def test_runfile_validates_openmx():
    base = dict(word="Tr[JXSZ delta(E-H) VY dG-(E)/dE]", expansion={"moments": 8},
                basis={"name": "eigen", "supercell": [2, 2, 1]})
    with pytest.raises(ValueError, match="requires 'file'"):
        RunFile.from_dict(dict(base, model={"kind": "openmx"}))
    # fermi_energy selects the energy origin (H - x S); spin stays refused.
    rf = RunFile.from_dict(dict(base, model={"kind": "openmx", "file": "x.scfout",
                                             "params": {"fermi_energy": 0.1}}))
    assert rf.model["params"]["fermi_energy"] == 0.1
    with pytest.raises(ValueError, match="spin"):
        RunFile.from_dict(dict(base, model={"kind": "openmx", "file": "x.scfout",
                                            "params": {"spin": "labels"}}))
