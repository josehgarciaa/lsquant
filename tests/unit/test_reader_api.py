"""Reader API alignments: the dipoles keyword, seed names, the memory bound,
the dipole cache path and the GPAW Hermitian fold.

The readers that build the connection of the covariant velocity all take
``dipoles`` (default False): with ``dipoles=True`` the model carries
``model.dipoles`` and with ``dipoles=False`` it stays None. ``connection``
(PySCF) and ``positions`` (Wannier90) remain as deprecated aliases that win
with a DeprecationWarning. gpaw and abacus name their models after the file
stem; ``model.params.dipole_cache`` resolves like ``model.file``;
``bloch_memory_mb`` bounds the workspace of the orthogonal models and the
run-file value wins over the per-kind default (4096 MiB for PySCF, else 256).
The GPAW connection is folded Hermitian after loading (old caches hold the
unfolded A and stay valid).
"""
import hashlib
import os
import sys
import warnings

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")

from lsquant.models.registry import CONNECTORS, load_runfile_model  # noqa: E402

K = np.array([0.13, 0.27, 0.0])
CELLS = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]


def _rng_blocks(rng, cells, n, scale, *, identity=False):
    mats = {}
    for R in cells:
        neg = tuple(-x for x in R)
        if neg in mats:
            mats[R] = mats[neg].conj().T
            continue
        m = scale * (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
        if R == neg:
            m = 0.5 * (m + m.conj().T) + (np.eye(n) if identity else 0)
        mats[R] = m
    return mats


def _wannier_dir(tmp_path, rng):
    from lsquant.models.orthogonal import load_wannier_metric_model
    n = 3
    mats = _rng_blocks(rng, CELLS, n, 0.4)
    np.savetxt(tmp_path / "w.uc", np.diag([3.0, 3.2, 20.0]))
    (tmp_path / "w.xyz").write_text(
        "3\ncomment\nWF1 0.0 0.0 0.0\nWF2 1.0 0.2 0.0\nWF3 0.3 1.1 0.0\n")
    with open(tmp_path / "w_hr.dat", "w", encoding="utf-8") as fh:
        fh.write("written by test_reader_api\n3\n3\n1 1 1\n")
        for R in CELLS:
            for i in range(3):
                for j in range(3):
                    fh.write(f"{R[0]} {R[1]} {R[2]} {i + 1} {j + 1} "
                             f"{mats[R][i, j].real:.16e} {mats[R][i, j].imag:.16e}\n")
    pos = 0.1 * rng.normal(size=(len(CELLS), n, n))
    with open(tmp_path / "w_r.dat", "w", encoding="utf-8") as fh:
        fh.write("test position matrix\n")
        for R, m in zip(CELLS, pos):
            for i in range(n):
                for j in range(n):
                    fh.write(f"{R[0]} {R[1]} {R[2]} {i + 1} {j + 1} "
                             f"{m[i, j].real:.16e} {m[i, j].imag:.16e}\n")
    return lambda **kw: load_wannier_metric_model(str(tmp_path), "w", **kw)


def _paoflow_file(tmp_path, rng):
    from lsquant.models.orthogonal import load_paoflow_model
    n = 4
    mats = _rng_blocks(rng, CELLS, n, 0.4)
    spin = np.zeros((3, n, n))
    spin[2] = np.eye(n)
    np.savez(tmp_path / "p.npz", uc=np.diag([3.0, 3.1, 18.0]), tau=np.zeros((n, 3)),
             r_int=np.array(CELLS, int), hoppings=np.array([mats[R] for R in CELLS]),
             spin=spin, fermi_energy=0.0, mesh=np.array([1, 1, 1]), edge_ratio=0.0)
    return lambda **kw: load_paoflow_model(str(tmp_path / "p.npz"), **kw)


def _pyscf_file(tmp_path, rng):
    from lsquant.models.pyscf_lcao import load_pyscf_model
    nao = 2
    h = _rng_blocks(rng, CELLS, nao, 0.4)
    s = _rng_blocks(rng, CELLS, nao, 0.05, identity=True)
    hoppings = np.zeros((len(CELLS), 2 * nao, 2 * nao), complex)
    for i, R in enumerate(CELLS):
        hoppings[i] = np.kron(h[R], np.eye(2))
    np.savez(tmp_path / "pyscf.npz", uc=np.diag([3.0, 3.1, 18.0]),
             tau=np.zeros((2 * nao, 3)), r_int=np.array(CELLS, int), hoppings=hoppings,
             overlap_orb=np.array([s[R] for R in CELLS]),
             position_orb=np.zeros((3, len(CELLS), nao, nao)), fermi_energy=0.0,
             mesh=np.array([1, 1, 1]), edge_ratio=0.0, nelectron=2, basis="test")
    return lambda **kw: load_pyscf_model(str(tmp_path / "pyscf.npz"), **kw)


def _gpaw_file(tmp_path, rng):
    from lsquant.models.gpaw import FORMAT, load_gpaw_model
    lattice = np.diag([3.0, 3.2, 15.0])
    positions = np.array([[0.1, 0.2, 7.0]])
    norb = np.array([2])
    mesh = (4, 3, 1)
    kpts = np.array([[i / mesh[0], j / mesh[1], 0.0]
                     for i in range(mesh[0]) for j in range(mesh[1])])
    cells = [(a, b, 0) for a in (-1, 0, 1) for b in (-1, 0, 1)]

    def bloch(mats):
        return np.array([sum(np.exp(2j * np.pi * np.dot(k, R)) * m for R, m in mats.items())
                         for k in kpts])

    h = _rng_blocks(rng, cells, 2, 0.3)
    s = _rng_blocks(rng, cells, 2, 0.02, identity=True)
    np.savez(tmp_path / "export.npz", format=np.array(FORMAT), kpts=kpts, mesh=np.array(mesh),
             H=np.conj(bloch(h)), S=np.conj(bloch(s)), soc=False, lattice=lattice,
             positions=positions, symbols=np.array(["Si"]), orbitals_per_atom=norb,
             fermi_energy=0.0)
    return lambda **kw: load_gpaw_model(str(tmp_path / "export.npz"), soc=False, **kw)


def _abacus_dir(tmp_path, rng, *, with_position=True):
    from lsquant.models.abacus import load_abacus_model
    n = 2
    h = _rng_blocks(rng, CELLS, n, 0.4)
    s = _rng_blocks(rng, CELLS, n, 0.05, identity=True)
    (tmp_path / "STRU").write_text("""ATOMIC_SPECIES
Bi 208.98 Bi.upf

LATTICE_CONSTANT
7.5

LATTICE_VECTORS
1.0 0.0 0.0
0.5 0.8 0.0
0.0 0.0 4.0

ATOMIC_POSITIONS
Cartesian_angstrom
Bi
0.0
1
0.0 0.0 5.0 1 1 1
""")

    def write_csr(name, mats):
        out = ["STEP: 0", f"Matrix Dimension of H(R): {n}", f"Matrix number of H(R): {len(mats)}"]
        for R, m in zip(CELLS, mats):
            m = np.real(m)
            rows, cols = np.nonzero(m)
            out.append(f"{R[0]} {R[1]} {R[2]} {len(rows)}")
            out.append(" ".join(repr(float(v)) for v in m[rows, cols]))
            out.append(" ".join(map(str, cols)))
            out.append(" ".join(map(str, np.concatenate(
                [[0], np.cumsum(np.bincount(rows, minlength=n))]))))
        (tmp_path / name).write_text("\n".join(out) + "\n")

    write_csr("hr_tot", np.array([h[R] for R in CELLS]))
    write_csr("srs1_nao.csr", np.array([s[R] for R in CELLS]))
    if with_position:
        q = np.zeros((3, len(CELLS), n, n))
        q[0] = 0.2
        out = ["STEP: 0", f"Matrix Dimension of r(R): {n}", f"Matrix number of r(R): {len(CELLS)}"]
        for r, R in enumerate(CELLS):
            out.append(f"{R[0]} {R[1]} {R[2]}")
            for a in range(3):
                rows, cols = np.nonzero(q[a, r])
                out.append(str(len(rows)))
                out.append(" ".join(repr(float(v)) for v in q[a, r][rows, cols]))
                out.append(" ".join(map(str, cols)))
                out.append(" ".join(map(str, np.concatenate(
                    [[0], np.cumsum(np.bincount(rows, minlength=n))]))))
        (tmp_path / "rr.csr").write_text("\n".join(out) + "\n")

    def load(**kw):
        if with_position:
            kw.setdefault("position_file", str(tmp_path / "rr.csr"))
        return load_abacus_model(str(tmp_path / "hr_tot"),
                                 overlap_file=str(tmp_path / "srs1_nao.csr"),
                                 stru_file=str(tmp_path / "STRU"),
                                 orbitals_per_element={"Bi": n}, **kw)
    return load


@pytest.mark.parametrize("kind,builder", [
    ("wannier", _wannier_dir),
    ("pyscf", _pyscf_file),
    ("abacus", _abacus_dir),
])
def test_dipoles_false_leaves_connection_none(tmp_path, kind, builder):
    load = builder(tmp_path, np.random.default_rng(41))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert load(dipoles=False).dipoles is None


@pytest.mark.parametrize("kind,builder", [
    ("wannier", _wannier_dir),
    ("pyscf", _pyscf_file),
    ("abacus", _abacus_dir),
])
def test_dipoles_true_builds_connection(tmp_path, kind, builder):
    assert builder(tmp_path, np.random.default_rng(42))(dipoles=True).dipoles is not None


def test_wannier_positions_alias_warns_and_wins(tmp_path):
    load = _wannier_dir(tmp_path, np.random.default_rng(31))
    assert load(dipoles=True).dipoles is not None
    assert load(dipoles=False).dipoles is None
    with pytest.warns(DeprecationWarning, match="positions"):
        assert load(positions=True).dipoles is not None
    with pytest.warns(DeprecationWarning, match="positions"):
        assert load(dipoles=False, positions=True).dipoles is not None   # the alias wins


def test_pyscf_connection_alias_warns_and_wins(tmp_path):
    load = _pyscf_file(tmp_path, np.random.default_rng(43))
    assert load().dipoles is None            # new default: no connection
    with pytest.warns(DeprecationWarning, match="connection"):
        assert load(connection=True).dipoles is not None
    with pytest.warns(DeprecationWarning, match="connection"):
        assert load(dipoles=True, connection=False).dipoles is None      # the alias wins


@pytest.mark.parametrize("kind,builder,stem", [
    ("gpaw", _gpaw_file, "export"),
    ("abacus", _abacus_dir, "hr_tot"),
])
def test_seed_name_is_the_file_stem(tmp_path, kind, builder, stem):
    load = builder(tmp_path, np.random.default_rng(44))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        assert load(dipoles=False).seed_name == stem


def test_runfile_applies_memory_bound(tmp_path):
    load = _pyscf_file(tmp_path, np.random.default_rng(45))
    model = {"kind": "pyscf", "file": str(tmp_path / "pyscf.npz"), "params": {}}
    assert load_runfile_model(model, "bare").bloch_memory_mb == 4096.0   # pyscf default stays
    model["params"] = {"bloch_memory_mb": 64.0}
    assert load_runfile_model(model, "bare").bloch_memory_mb == 64.0     # run file wins
    assert CONNECTORS["pyscf"].default_bloch_memory_mb == 4096.0
    for connector in CONNECTORS.values():
        if connector.kind != "pyscf":
            assert connector.default_bloch_memory_mb == 256.0


def test_orthogonal_model_memory_bound_bit_identical(tmp_path):
    from lsquant.models.orthogonal import load_paoflow_model
    rng = np.random.default_rng(46)
    n = 6
    cells = [(a, b, 0) for a in range(-2, 3) for b in range(-2, 3)]
    mats = _rng_blocks(rng, cells, n, 0.4)
    spin = np.zeros((3, n, n))
    spin[2] = np.eye(n)
    np.savez(tmp_path / "m.npz", uc=np.diag([3.0, 3.1, 18.0]), tau=np.zeros((n, 3)),
             r_int=np.array(cells, int), hoppings=np.array([mats[R] for R in cells]),
             spin=spin, fermi_energy=0.0, mesh=np.array([1, 1, 1]), edge_ratio=0.0)
    model = load_paoflow_model(str(tmp_path / "m.npz"))
    k = np.tile(K, (5, 1))
    model.bloch_memory_mb = 1e-4     # far below one k point: one k at a time
    tiny_hk = model.Hk(k)
    tiny_dh, tiny_dh_cell = model.derivatives(k, 0)[0], model.derivatives(k, 0, cell=True)[0]
    model.bloch_memory_mb = 1e6
    np.testing.assert_array_equal(model.Hk(k), tiny_hk)
    np.testing.assert_array_equal(model.derivatives(k, 0)[0], tiny_dh)
    np.testing.assert_array_equal(model.derivatives(k, 0, cell=True)[0], tiny_dh_cell)


def test_wannier_spin_blocks_memory_bound():
    from lsquant.models.orthogonal import WannierSpin
    rng = np.random.default_rng(47)
    spin = WannierSpin.__new__(WannierSpin)   # blocks() only needs the arrays
    spin.uc = np.eye(3) * 2.0
    spin.tau = rng.normal(size=(4, 3))
    spin.r_int = np.array([[a, b, 0] for a in (-1, 0, 1) for b in (-1, 0, 1)], int)
    spin.spin_r = 0.01 * rng.normal(size=(3, len(spin.r_int), 4, 4))
    k = np.tile(K, (7, 1))
    spin.bloch_memory_mb = 1e-4
    tiny = spin.blocks(k, "z").copy()
    spin.bloch_memory_mb = 1e6
    np.testing.assert_array_equal(spin.blocks(k, "z"), tiny)


def test_gpaw_dipole_cache_path_resolves_like_model_file(tmp_path):
    from lsquant.runfile import RunFile
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "export.npz").write_bytes(b"placeholder")
    run = {"kind": "gpaw", "file": "sub/export.npz",
           "params": {"dipole_cache": "cache/dipoles.npz"}}
    rf = RunFile.from_dict(dict(word="Tr[delta(E-H)]", expansion={"moments": 8},
                                basis={"name": "eigen", "supercell": [2, 2, 1]}, model=run),
                           repo_root=str(tmp_path))
    assert rf.model["params"]["dipole_cache"] == str(tmp_path / "cache" / "dipoles.npz")
    assert CONNECTORS["gpaw"].path_params == ("dipole_cache",)


def test_gpaw_connection_is_folded(tmp_path):
    from lsquant.models.gpaw import FORMAT, load_gpaw_model
    from lsquant.models.lcao import opposite_pairs
    rng = np.random.default_rng(49)
    cells = [(a, b, 0) for a in (-1, 0, 1) for b in (-1, 0, 1)]
    mesh = (4, 3, 1)
    kpts = np.array([[i / mesh[0], j / mesh[1], 0.0]
                     for i in range(mesh[0]) for j in range(mesh[1])])

    def bloch(mats):
        return np.array([sum(np.exp(2j * np.pi * np.dot(k, R)) * m for R, m in mats.items())
                         for k in kpts])

    h = _rng_blocks(rng, cells, 2, 0.3)
    s = _rng_blocks(rng, cells, 2, 0.02, identity=True)
    np.savez(tmp_path / "export.npz", format=np.array(FORMAT), kpts=kpts, mesh=np.array(mesh),
             H=np.conj(bloch(h)), S=np.conj(bloch(s)), soc=False,
             lattice=np.diag([3.0, 3.2, 15.0]), positions=np.array([[0.1, 0.2, 7.0]]),
             symbols=np.array(["Si"]), orbitals_per_atom=np.array([2]), fermi_energy=0.0)
    order = 4
    base = load_gpaw_model(str(tmp_path / "export.npz"), soc=False)
    key = hashlib.sha256(np.ascontiguousarray(np.load(tmp_path / "export.npz")["S"]).tobytes()
                         + base.r_int.tobytes() + str(order).encode()).hexdigest()
    unfolded = rng.normal(size=(3, len(base.r_int), 2, 2))   # not Hermitian-paired
    cache = tmp_path / "cache.npz"
    np.savez(cache, key=key, S=np.asarray(base.overlaps.real), A=unfolded)
    model = load_gpaw_model(str(tmp_path / "export.npz"), soc=False, dipoles=True,
                            quadrature_order=order, dipole_cache=str(cache))
    assert model.dipoles is not None
    partner = opposite_pairs(model.r_int)
    mirrored = model.dipoles[:, partner, :, :].conj().transpose(0, 1, 3, 2)
    np.testing.assert_allclose(model.dipoles, mirrored, atol=1e-14)
    # the residual is measured before the fold, with its historical normalisation
    assert model.metadata["dipole_hermiticity_residual"] == pytest.approx(
        float(np.abs(unfolded - unfolded[:, partner, :, :].transpose(0, 1, 3, 2)).max()))
