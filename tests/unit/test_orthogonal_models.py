"""Orthogonal models in the metric executor: Wannier90 exact spin and PAOFLOW.

Three checks:
- ``read_wannier_spin`` rebuilds ``V^+ S V`` when the rows of ``_u_dis.mat``
  list only the in-window bands, packed at the top, with a window that moves
  between k points;
- the Wannier adapter gives :class:`BlochMetricBasis` the same velocities as
  :meth:`WannierModel.velocity_blocks` (up to the Wannier90 sign) for a
  position matrix with off-diagonal and R != 0 elements;
- a PAOFLOW model written by ``write_paoflow_model`` reproduces PAOFLOW's
  ``H(k) = sum_R exp(-i k.R) H(R)`` on the grid's k points and carries its
  on-site spin.
"""
import os
import shutil
import sys
from types import SimpleNamespace

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GRAPHENE = os.path.join(REPO, "goldens", "stage0", "bases", "graphene", "model")
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")
FortranFile = pytest.importorskip("scipy.io").FortranFile
from lsquant.bases import BlochMetricBasis  # noqa: E402
from lsquant.models.orthogonal import (OrthogonalModel, WannierSpin, load_paoflow_model,  # noqa: E402
                                       load_wannier_metric_model, read_wannier_spin,
                                       write_paoflow_model, write_wannier_geometry)
from lsquant.models.wannier import load_wannier_model  # noqa: E402
from lsquant.runfile import RunFile, uses_metric_executor  # noqa: E402


def _random_hermitian(rng, n):
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return a + a.conj().T


def _write_umat(path, kpts, mats, header_dims):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("written by test_orthogonal_models\n")
        fh.write(" ".join(str(d) for d in header_dims) + "\n")
        for k, m in zip(kpts, mats):
            fh.write("\n" + " ".join(f"{x:.10f}" for x in k) + "\n")
            for z in m.ravel(order="F"):
                fh.write(f"{z.real:.15e} {z.imag:.15e}\n")


def test_exact_spin_with_packed_disentanglement_rows(tmp_path):
    rng = np.random.default_rng(3)
    nk, nb, nw, seed = 4, 7, 4, "toy"
    kpts = np.array([[i / 2, j / 2, 0.0] for j in range(2) for i in range(2)])
    win = (-1.0, 1.0)
    eig = np.zeros((nk, nb))
    U, Ud, V = [], [], []
    for k in range(nk):
        # Five bands in the window, starting at a different band index per k.
        first = k % 2
        e = np.full(nb, 5.0)
        e[:first] = -5.0
        e[first:first + 5] = np.linspace(-0.9, 0.9, 5)
        eig[k] = e
        q, _ = np.linalg.qr(rng.normal(size=(5, nw)) + 1j * rng.normal(size=(5, nw)))
        u, _ = np.linalg.qr(rng.normal(size=(nw, nw)) + 1j * rng.normal(size=(nw, nw)))
        packed = np.zeros((nb, nw), complex)
        packed[:5] = q
        Ud.append(packed)
        U.append(u)
        full = np.zeros((nb, nw), complex)
        full[first:first + 5] = q @ u
        V.append(full)
    S_b = np.array([[_random_hermitian(rng, nb) for _ in range(nk)] for _ in range(3)])
    d = str(tmp_path)
    _write_umat(os.path.join(d, seed + "_u.mat"), kpts, U, (nk, nw, nw))
    _write_umat(os.path.join(d, seed + "_u_dis.mat"), kpts, Ud, (nk, nw, nb))
    with open(os.path.join(d, seed + ".eig"), "w", encoding="utf-8") as fh:
        for k in range(nk):
            for n in range(nb):
                fh.write(f"{n + 1} {k + 1} {eig[k, n]:.10f}\n")
    with open(os.path.join(d, seed + ".win"), "w", encoding="utf-8") as fh:
        fh.write(f"num_wann = {nw}\nnum_bands = {nb}\ndis_win_min = {win[0]}\ndis_win_max = {win[1]}\n")
    with FortranFile(os.path.join(d, seed + ".spn"), "w") as f:
        f.write_record(np.array("written by test".ljust(60), dtype="S60"))
        f.write_record(np.array([nb, nk], dtype=np.int32))
        for k in range(nk):
            rec = [[S_b[a, k, m, n] for a in range(3)] for n in range(nb) for m in range(n + 1)]
            f.write_record(np.array(rec, dtype=np.complex128).view(np.float64).ravel())
    k_read, S_w = read_wannier_spin(d, seed)
    np.testing.assert_allclose(k_read, kpts)
    expected = np.einsum("kba,xkbc,kcd->xkad", np.conj(V), S_b, np.array(V))
    np.testing.assert_allclose(S_w, expected, atol=1e-10)
    # The interpolated spin returns S_W on the mesh it came from (one k point
    # along z must give the single image R_z = 0, not two).
    uc = np.diag([3.0, 4.0, 20.0])
    spin = WannierSpin(d, seed, uc, np.zeros((nw, 3)))
    assert spin.mesh == (2, 2, 1)
    k_cart = kpts @ (2 * np.pi * np.linalg.inv(uc).T)
    for a, axis in enumerate("xyz"):
        herm = 0.5 * (expected[a] + np.conj(np.swapaxes(expected[a], 1, 2)))
        np.testing.assert_allclose(spin.blocks(k_cart, axis), herm, atol=1e-10)


def _write_r_dat(model_dir, seed, rng):
    """Write a ``_r.dat`` with random elements on every R of the ``_hr.dat``."""
    with open(os.path.join(model_dir, seed + "_hr.dat"), encoding="utf-8") as fh:
        lines = [ln.split() for ln in fh.read().splitlines()[1:] if ln.strip()]
    width, count = int(lines[0][0]), int(lines[1][0])
    index, seen = 2, 0
    while seen < count:
        seen += len(lines[index])
        index += 1
    r_list = list(dict.fromkeys(tuple(map(int, ln[:3])) for ln in lines[index:]))
    out = ["written by test_orthogonal_models", str(width), str(count)]
    for r in r_list:
        for n in range(width):
            for m in range(width):
                vals = rng.normal(scale=0.3, size=6)
                out.append(f"{r[0]} {r[1]} {r[2]} {m + 1} {n + 1} " + " ".join(f"{x:.10f}" for x in vals))
    with open(os.path.join(model_dir, seed + "_r.dat"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")


@pytest.mark.parametrize("mode", ["bare", "berry_connection", "covariant"])
def test_wannier_adapter_matches_velocity_blocks(tmp_path, mode):
    target = tmp_path / "graphene"
    shutil.copytree(GRAPHENE, target)
    _write_r_dat(str(target), "graphene", np.random.default_rng(11))
    adapter = load_wannier_metric_model(str(target), "graphene", dipoles=True)
    basis = BlochMetricBasis(adapter, [5, 4, 1], velocity_mode=mode)
    wannier = load_wannier_model(str(target), seed="graphene")
    for axis in "xy":
        reference = wannier.velocity_blocks(basis.k_points, axis, mode=mode)
        np.testing.assert_allclose(basis.operator_blocks("V" + axis.upper()), -reference, atol=1e-10)
    np.testing.assert_allclose(basis.hk_array, wannier.Hk(basis.k_points), atol=1e-12)


def test_metric_basis_uses_model_spin_in_spin_current():
    rng = np.random.default_rng(5)
    uc = np.diag([2.0, 2.5, 10.0])
    tau = rng.uniform(0, 1, size=(4, 3))
    tau[:, 2] = 0
    r_int = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0]])
    h0 = _random_hermitian(rng, 4)
    t = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    t2 = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    hoppings = np.array([h0, t, t.conj().T, t2, t2.conj().T])
    spin = np.array([_random_hermitian(rng, 4) for _ in range(3)])
    model = OrthogonalModel(uc, tau, r_int, hoppings, spin=spin)
    basis = BlochMetricBasis(model, [3, 3, 1], velocity_mode="berry_connection")
    sz = basis.operator_blocks("SZ")
    np.testing.assert_allclose(sz, np.broadcast_to(spin[2], sz.shape))
    vx = basis.operator_blocks("VX")
    np.testing.assert_allclose(basis.operator_blocks("JXSZ"), 0.5 * (vx @ sz + sz @ vx), atol=1e-12)


def test_paoflow_round_trip(tmp_path):
    rng = np.random.default_rng(9)
    n, alat_bohr = 4, 7.0
    a_vectors = np.diag([1.0, 1.2, 3.0])
    hrs = np.zeros((n, n, 3, 3, 1, 1), complex)
    for i in range(3):
        for j in range(3):
            hrs[:, :, i, j, 0, 0] = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    # Hermiticity of H(k): H(-R) = H(R)^+ on the 3 x 3 FFT grid.
    for i in range(3):
        for j in range(3):
            hrs[:, :, (-i) % 3, (-j) % 3, 0, 0] = hrs[:, :, i, j, 0, 0].conj().T
    hrs[:, :, 0, 0, 0, 0] = 0.5 * (hrs[:, :, 0, 0, 0, 0] + hrs[:, :, 0, 0, 0, 0].conj().T)
    tau_bohr = rng.uniform(0, 3, size=(n, 3))
    dnm = tau_bohr[:, None, :] - tau_bohr[None, :, :]
    sj = np.array([0.5 * _random_hermitian(rng, n) for _ in range(3)])
    arrays = dict(HRs=hrs, a_vectors=a_vectors, Dnm=dnm, Sj=sj)
    attrs = dict(alat=alat_bohr)
    fake = SimpleNamespace(data_controller=SimpleNamespace(data_dicts=lambda: (arrays, attrs)))
    path = str(tmp_path / "toy_paoflow.npz")
    # Random H(R) does not decay, so the locality check must warn.
    with pytest.warns(RuntimeWarning, match="has not decayed"):
        write_paoflow_model(fake, path)
    assert float(np.load(path)["edge_ratio"]) > 1e-2
    model = load_paoflow_model(path)
    bohr = 0.529177210903
    uc = a_vectors * alat_bohr * bohr
    # On the grid's own k points H(k) does not depend on which lattice image
    # carries each element (the writer picks the Wigner-Seitz one).
    frac = np.array([[i / 3, j / 3, 0] for i in range(3) for j in range(3)])
    k = frac @ (2 * np.pi * np.linalg.inv(uc).T)
    # PAOFLOW convention in the cell gauge, then the atom-gauge phase.
    R = np.array([[i if i <= 1 else i - 3, j if j <= 1 else j - 3, 0] for i in range(3) for j in range(3)]) @ uc
    blocks = np.array([hrs[:, :, i, j, 0, 0] for i in range(3) for j in range(3)])
    h_cell = np.einsum("kr,rab->kab", np.exp(-1j * k @ R.T), blocks)
    tau = (tau_bohr - tau_bohr[0]) * bohr
    p = np.exp(1j * k @ tau.T)
    np.testing.assert_allclose(model.Hk(k), p.conj()[:, :, None] * h_cell * p[:, None, :], atol=1e-10)
    np.testing.assert_allclose(model.spin_blocks(k, "z")[0], 2 * sj[2])


def test_runfile_routes_orthogonal_models():
    assert uses_metric_executor({"kind": "paoflow", "file": "x.npz"})
    assert uses_metric_executor({"kind": "wannier", "dir": ".", "seed": "s", "params": {"spin": "exact"}})
    assert not uses_metric_executor({"kind": "wannier", "dir": ".", "seed": "s"})
    base = dict(word="Tr[JXSZ delta(E-H) VY dG-(E)/dE]", expansion={"moments": 8},
                basis={"name": "eigen", "supercell": [2, 2, 1]})
    with pytest.raises(ValueError, match="velocity_mode"):
        RunFile.from_dict(dict(base, model={"kind": "paoflow", "file": "x.npz",
                                            "params": {"velocity_mode": "covariant"}}))
    with pytest.raises(ValueError, match="requires 'file'"):
        RunFile.from_dict(dict(base, model={"kind": "paoflow"}))
    with pytest.raises(ValueError, match="spin"):
        RunFile.from_dict(dict(base, model={"kind": "wannier", "dir": ".", "seed": "s",
                                            "params": {"spin": "diagonal"}}))


def test_write_wannier_geometry_from_win_and_centres(tmp_path):
    d, seed = str(tmp_path), "w"
    with open(os.path.join(d, seed + ".win"), "w", encoding="utf-8") as fh:
        fh.write("num_wann = 2\nbegin unit_cell_cart\nbohr\n10.0 0 0\n0 12.0 0\n0 0 30.0\nend unit_cell_cart\n")
    with open(os.path.join(d, seed + "_centres.xyz"), "w", encoding="utf-8") as fh:
        fh.write("4\nWannier centres\nX 0.1 0.2 0.3\nX 1.0 2.0 3.0\nPd 0 0 0\nSe 1 1 1\n")
    uc_path, xyz_path = write_wannier_geometry(d, seed)
    np.testing.assert_allclose(np.loadtxt(uc_path), np.diag([10.0, 12.0, 30.0]) * 0.529177210903)
    lines = open(xyz_path, encoding="utf-8").read().split("\n")
    assert lines[0] == "2" and lines[2].split()[1:] == ["1.00000000", "2.00000000", "3.00000000"]
