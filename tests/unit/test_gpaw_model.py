"""GPAW LCAO exports in the metric executors.

The fixtures are files in the layout of ``export_gpaw_lcao``, written here
from random short-ranged H(R), S(R) and spin-orbit A(R) of a two-atom
crystal (3 and 2 orbitals): the k matrices on a 6 x 5 x 1 mesh are their
Bloch sums ``sum_R exp(+2 pi i k.R) X(R)``. The checks:
- the inverse transform puts every block at its nearest image and the model
  reproduces ``H(k) tensor 1 + A(k)`` and ``S(k) tensor 1`` at k points off
  the mesh (the transform is exact for a range below half the mesh);
- without spin-orbit coupling the model is spinless; the Fermi level of the
  file, or the one given, is removed as H - E_F S;
- a file that is not an export, a mesh that the k points do not fill, and
  ``soc=True`` on a file without A are refused;
- the run file routes ``model.kind = "gpaw"`` and runs the DOS and the spin
  Hall word.
The exporter itself needs GPAW; it is checked against ``soc_eigenstates``
when GPAW is installed.
"""
import json
import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant.models.gpaw import FORMAT, load_gpaw_model, real_space_blocks  # noqa: E402
from lsquant.models.lcao import cell_gauge_hk  # noqa: E402
from lsquant.models.orthogonal import load_runfile_model  # noqa: E402
from lsquant.runfile import RunFile, execute_runfile, uses_metric_executor  # noqa: E402

LATTICE = np.array([[3.1, 0.0, 0.0], [0.4, 3.3, 0.0], [0.0, 0.0, 15.0]])
POSITIONS = np.array([[0.2, 0.3, 7.0], [1.9, 1.8, 7.6]])
NORB = np.array([3, 2])
CELLS = [(a, b, 0) for a in range(-1, 2) for b in range(-1, 2)]
MESH = (6, 5, 1)
EF = -1.7


def _hermitian_set(rng, n, scale, onsite=0.0, extra=()):
    mats = {}
    for R in CELLS:
        neg = tuple(-x for x in R)
        if neg in mats:
            mats[R] = mats[neg].conj().swapaxes(-1, -2)
            continue
        m = scale * (rng.normal(size=extra + (n, n)) + 1j * rng.normal(size=extra + (n, n)))
        if R == neg:
            m = 0.5 * (m + m.conj().swapaxes(-1, -2)) + onsite * np.eye(n)
        mats[R] = m
    return mats


def _spin_orbit(rng, n, scale):
    """``A^{ss'}(R)`` with ``A(-R) = A(R)^+`` as a 2n x 2n matrix, stored as (2, 2, n, n)."""
    full = _hermitian_set(rng, 2 * n, scale)
    return {R: np.array([[m[:n, :n], m[:n, n:]], [m[n:, :n], m[n:, n:]]]) for R, m in full.items()}


def _bloch(mats, k_frac):
    return np.array([sum(np.exp(2j * np.pi * np.dot(k, R)) * m for R, m in mats.items()) for k in k_frac])


def _mesh():
    return np.array([[i / MESH[0], j / MESH[1], 0.0] for i in range(MESH[0]) for j in range(MESH[1])])


def _fixture(seed=3):
    rng = np.random.default_rng(seed)
    n = int(NORB.sum())
    H = _hermitian_set(rng, n, 0.3, onsite=0.0)
    S = _hermitian_set(rng, n, 0.02, onsite=1.0)
    A = _spin_orbit(rng, n, 0.05)
    return H, S, A


def _write(path, H, S, A, *, soc=True, fmt=FORMAT, kpts=None, mesh=MESH):
    k = _mesh() if kpts is None else kpts
    np.savez(path, format=fmt, kpts=k, mesh=np.array(mesh), H=_bloch(H, k), S=_bloch(S, k),
             A=_bloch(A, k), soc=soc, lattice=LATTICE, positions=POSITIONS, symbols=np.array(["Pd", "Se"]),
             orbitals_per_atom=NORB, fermi_energy=EF)


def _spinor(h, a):
    return np.block([[h + a[0, 0], a[0, 1]], [a[1, 0], h + a[1, 1]]])


def _interleave(n):
    """Permutation from the [up block, down block] layout to interleaved spinors."""
    return np.array([s * n + i for i in range(n) for s in range(2)])


def test_nearest_image_and_exact_transform():
    H, _, _ = _fixture()
    k = _mesh()
    blocks = real_space_blocks(k, MESH, _bloch(H, k), LATTICE, POSITIONS, NORB)
    starts = np.concatenate([[0], np.cumsum(NORB)])
    for (R, i, j), block in blocks.items():
        if np.abs(block).max() > 1e-12:
            assert R in CELLS
            np.testing.assert_allclose(block, H[R][starts[i]:starts[i + 1], starts[j]:starts[j + 1]], atol=1e-12)


def test_model_matches_spinor_bloch_sum(tmp_path):
    H, S, A = _fixture()
    path = tmp_path / "pdse2_export.npz"
    _write(path, H, S, A)
    model = load_gpaw_model(path)
    assert model.spinful and model.metadata["source"] == "gpaw"
    assert model.metadata["fermi_energy_subtracted_eV"] == EF
    k = np.array([[0.13, 0.27, 0.0], [-0.31, 0.08, 0.0], [0.45, -0.22, 0.0]])
    perm = _interleave(int(NORB.sum()))
    Hk, Sk, Ak = _bloch(H, k), _bloch(S, k), _bloch(A, k)
    got_h = cell_gauge_hk(model, k)
    k_cart = k @ (2 * np.pi * np.linalg.inv(LATTICE).T)
    from lsquant.models.bloch import bloch_sum
    got_s = bloch_sum(model.r_cart, model.tau, model.overlaps, k_cart, cell=True)
    for n in range(len(k)):
        want_s = np.kron(np.eye(2), Sk[n])[np.ix_(perm, perm)]
        want_h = _spinor(Hk[n], Ak[n])[np.ix_(perm, perm)] - EF * want_s
        np.testing.assert_allclose(got_h[n], want_h, atol=1e-10)
        np.testing.assert_allclose(got_s[n], want_s, atol=1e-12)


def test_spinless_and_fermi_override(tmp_path):
    H, S, A = _fixture()
    path = tmp_path / "export.npz"
    _write(path, H, S, A, soc=False)
    model = load_gpaw_model(path, soc=False, fermi_energy=0.4)
    assert not model.spinful and model.hoppings.shape[1] == NORB.sum()
    k = np.array([[0.21, -0.17, 0.0]])
    np.testing.assert_allclose(cell_gauge_hk(model, k)[0], _bloch(H, k)[0] - 0.4 * _bloch(S, k)[0], atol=1e-10)
    with pytest.raises(ValueError, match="without spin-orbit"):
        load_gpaw_model(path)


def test_refusals(tmp_path):
    H, S, A = _fixture()
    bad = tmp_path / "other.npz"
    _write(bad, H, S, A, fmt="something-else")
    with pytest.raises(ValueError, match="not an LSQuant GPAW export"):
        load_gpaw_model(bad)
    short = tmp_path / "short.npz"
    _write(short, H, S, A, kpts=_mesh()[:-1])
    with pytest.raises(ValueError, match="do not fill"):
        load_gpaw_model(short)


def test_runfile_routes_and_runs_gpaw(tmp_path):
    H, S, A = _fixture()
    path = tmp_path / "export.npz"
    _write(path, H, S, A)
    model = {"kind": "gpaw", "file": str(path), "params": {"fermi_energy": 0.1}}
    assert uses_metric_executor(model)
    assert load_runfile_model(model, "berry_connection").metadata["fermi_energy_subtracted_eV"] == 0.1
    base = dict(word="Tr[delta(E-H)]", expansion={"moments": 8}, basis={"name": "eigen", "supercell": [2, 2, 1]})
    assert RunFile.from_dict(dict(base, model=model)).model["kind"] == "gpaw"
    with pytest.raises(ValueError, match="requires 'file'"):
        RunFile.from_dict(dict(base, model={"kind": "gpaw"}))
    RunFile.from_dict(dict(base, model=dict(model, params={"velocity_mode": "covariant"})))
    with pytest.raises(ValueError, match="velocity_mode for gpaw"):
        RunFile.from_dict(dict(base, model=dict(model, params={"velocity_mode": "position"})))
    with pytest.raises(ValueError, match="no radial functions"):
        load_runfile_model(model, "covariant")
    with pytest.raises(ValueError, match="true or false"):
        RunFile.from_dict(dict(base, model=dict(model, params={"soc": "yes"})))
    for word, formula, name in [("Tr[delta(E-H)]", "dos", "dos"),
                                ("Tr[JXSZ delta(E-H) VY dG-(E)/dE]", "shc", "shc")]:
        out = tmp_path / name
        execute_runfile({
            "model": {"kind": "gpaw", "file": str(path)},
            "basis": {"name": "eigen", "supercell": [3, 3, 1]},
            "word": word,
            "expansion": {"moments": 16, "window": "eigen", "trace": "exact"},
            "reconstruction": {"kernel": "jackson", "formulas": [formula], "broadening_meV": 500.0,
                               "energy_grid": [-3.0, 3.0, 13]},
            "output": {"dir": str(out)},
        })
        assert json.loads((out / "run.resolved.json").read_text())["model"]["kind"] == "gpaw"


def test_exporter_matches_soc_eigenstates(tmp_path):
    """With GPAW installed: export a tiny LCAO calculation and compare with ``soc_eigenstates``."""
    pytest.importorskip("gpaw")
    import scipy.linalg as sl
    from ase.build import bulk
    from gpaw.calculator import GPAW

    from lsquant.models.gpaw import export_gpaw_lcao, read_gpaw_lcao

    atoms = bulk("Si", "diamond", a=5.43)
    atoms.calc = GPAW(mode="lcao", basis="dzp", xc="PBE", h=0.25, kpts={"size": (3, 3, 3), "gamma": True},
                      symmetry="off", nbands="nao", txt=None, parallel={"domain": 1})
    atoms.get_potential_energy()
    path = tmp_path / "si.npz"
    export_gpaw_lcao(atoms.calc, path, check=True)
    d = read_gpaw_lcao(path)
    for n in range(len(d["kpts"])):
        h = np.block([[d["H"][n] + d["A"][n, 0, 0], d["A"][n, 0, 1]], [d["A"][n, 1, 0], d["H"][n] + d["A"][n, 1, 1]]])
        e = sl.eigh(h, np.kron(np.eye(2), d["S"][n]), eigvals_only=True)
        np.testing.assert_allclose(e, d["soc_eigenvalues"][n], atol=1e-6)


C0, C1 = 0.28209479177387814, 0.4886025119029199     # GPAW's r^l Y_L for l = 0, 1
YL_TABLE = np.array([(0, C0, 0, 0, 0), (1, C1, 0, 1, 0), (2, C1, 0, 0, 1), (3, C1, 1, 0, 0)], float)


def _radial(alpha, rc, l, norm):
    from lsquant.models.gpaw import RADIAL_POINTS
    r = np.linspace(0.0, rc, RADIAL_POINTS)
    return norm * np.exp(-alpha * r**2) * np.ones_like(r) if l >= 0 else None


def _gaussian_atoms(positions, lattice, dO=0.0):
    """Atoms with a normalised s and p Gaussian (alpha = 1) and one s projector."""
    alpha, rc = 1.0, 6.0
    ns = 1 / (C0 * (np.pi / (2 * alpha)) ** 0.75)
    npx = 1 / (C1 * np.sqrt((np.pi / (2 * alpha)) ** 1.5 / (4 * alpha)))
    data = {"bohr": 1.0, "YL": YL_TABLE, "lattice": np.asarray(lattice, float),
            "positions": np.asarray(positions, float), "orbitals_per_atom": np.array([4] * len(positions))}
    for a in range(len(positions)):
        data[f"basis_l_{a}"] = np.array([0, 1])
        data[f"basis_rc_{a}"] = np.array([rc, rc])
        data[f"basis_f_{a}"] = np.array([_radial(alpha, rc, 0, ns), _radial(alpha, rc, 1, npx)])
        data[f"proj_l_{a}"] = np.array([0])
        data[f"proj_rc_{a}"] = np.array([3.0])
        data[f"proj_f_{a}"] = np.array([_radial(2.0, 3.0, 0, 1.0)])
        data[f"dO_{a}"] = np.array([[dO]])
        data[f"Delta_{a}"] = np.zeros((1, 1, 3))
    return data, alpha, ns, npx


def test_connection_on_site_analytic():
    from lsquant.models.gpaw import gpaw_connection
    data, alpha, ns, npx = _gaussian_atoms([[0.0, 0.0, 0.0]], np.diag([30.0, 30.0, 30.0]))
    S, A = gpaw_connection(data, np.array([[0, 0, 0]]), order=48)
    # linear interpolation of the tabulated radial functions: about 2e-6
    np.testing.assert_allclose(S[0], np.eye(4), atol=1e-5)
    want = ns * C0 * npx * C1 * (np.pi / (2 * alpha)) ** 1.5 / (4 * alpha)
    # <s|x|p_x>, <s|y|p_y>, <s|z|p_z>: p order y, z, x (GPAW's L = 1, 2, 3)
    for axis, orbital in ((0, 3), (1, 1), (2, 2)):
        assert abs(A[axis, 0, 0, orbital] - want) < 1e-5
        assert abs(A[axis, 0, orbital, 0] - want) < 1e-5
    assert np.abs(A[:, 0, 1:, 1:]).max() < 1e-10


def test_connection_two_atoms_translation_and_paw_overlap():
    from lsquant.models.gpaw import _Functions, _pair_integrals, gpaw_connection
    lattice = np.diag([4.0, 30.0, 30.0])
    positions = [[0.3, 0.2, 0.0], [1.9, -0.4, 0.7]]
    cells = np.array([[n, 0, 0] for n in range(-3, 4)])
    data0, *_ = _gaussian_atoms(positions, lattice, dO=0.0)
    data1, *_ = _gaussian_atoms(positions, lattice, dO=0.3)
    S0, A0 = gpaw_connection(data0, cells, order=40)
    S1, A1 = gpaw_connection(data1, cells, order=40)
    # A(-R)^T = A(R), here checked against an independent pseudo integral of the partner block.
    yl = {}
    for L, c, nx, ny, nz in YL_TABLE:
        yl.setdefault(int(L), []).append((c, int(nx), int(ny), int(nz)))
    f = [_Functions(data0, a, "basis", yl) for a in range(2)]
    R = np.array([1, 0, 0])
    d = R @ lattice + np.array(positions[0]) - np.array(positions[1])
    s, a = _pair_integrals(f[1], f[0], d, 40)
    n = int(np.where((cells == R).all(axis=1))[0][0])
    np.testing.assert_allclose(A0[:, n, 4:, :4], a, atol=1e-8)
    np.testing.assert_allclose(S0[n, 4:, :4], s, atol=1e-8)
    minus = [int(np.where((cells == -r).all(axis=1))[0][0]) for r in cells]
    np.testing.assert_allclose(A1, A1[:, minus].transpose(0, 1, 3, 2), atol=1e-12)
    # The PAW term adds sum <M|p> dO <p|N> to the overlap and shifts nothing else for a symmetric pair.
    assert np.abs(S1 - S0).max() > 1e-3
    assert np.all(np.linalg.eigvalsh(S1[cells.tolist().index([0, 0, 0])]) > 0)
