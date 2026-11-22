"""HamGNN Hamiltonians in the metric executors.

The fixtures are written here: a ``graph_data.npz`` with two structures
(plain objects in place of ``torch_geometric`` ``Data``; the export only
reads their attributes) and a ``prediction_hamiltonian.npy`` in the layout
of HamGNN's output, for the three spin modes. The checks:
- H(k) and S(k) of the LSQuant model equal the sum of HamGNN's
  ``band_cal.py`` (padded ``nao_max`` blocks, masked to the real orbitals,
  spin-major global order) after the spinor permutation, the unit
  conversion and the Fermi-energy shift;
- the second structure of the file is cut from the right rows; a wrong
  ``nao_max`` or spin mode is refused;
- ``prediction=None`` exports the DFT blocks of ``graph_data.npz``;
- the spin operator of the metric basis is ``S (I tensor sigma)``;
- the run file routes ``model.kind = "hamgnn"`` and runs DOS and spin Hall.
"""
import json
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant.bases import BlochMetricBasis  # noqa: E402
from lsquant.models.hamgnn import (BOHR_ANGSTROM, HARTREE_EV, load_hamgnn_model,  # noqa: E402
                                   main, write_hamgnn_model)
from lsquant.models.lcao import cell_gauge_hk  # noqa: E402
from lsquant.models.orthogonal import load_runfile_model  # noqa: E402
from lsquant.runfile import RunFile, execute_runfile, uses_metric_executor  # noqa: E402

NAO = 6
Z = np.array([34, 1])
MASK = {34: np.array([0, 3, 4, 5]), 1: np.array([0])}   # Se: s + p, H: s, inside nao_max = 6
NORB = np.array([4, 1])
LATTICE = np.array([[6.0, 0.0, 0.0], [2.6, 5.0, 0.0], [0.0, 0.0, 28.0]])      # Bohr
POSITIONS = np.array([[0.2, 0.3, 13.0], [3.1, 2.0, 14.0]])                    # Bohr
SHIFTS = [(1, 0, 0), (0, 1, 0), (1, -1, 0)]
K = np.array([[0.13, 0.27, 0.0], [-0.31, 0.08, 0.0], [0.45, -0.22, 0.0]])


def _edges():
    """Both directions of every bond, as HamGNN stores them (no self-edge at R = 0)."""
    edges = [(0, 1, (0, 0, 0)), (1, 0, (0, 0, 0))]
    for R in SHIFTS:
        neg = tuple(-x for x in R)
        for i, j in ((0, 0), (1, 1), (0, 1), (1, 0)):
            edges += [(i, j, R), (j, i, neg)]
    return edges


def _padded(rng, i, j, spinful, scale=1.0, real=False):
    """A random padded block, non-zero only on the real orbitals of atoms ``i`` and ``j``."""
    w = 2 * NAO if spinful else NAO
    out = np.zeros((w, w), complex)
    rows = np.concatenate([MASK[Z[i]], MASK[Z[i]] + NAO]) if spinful else MASK[Z[i]]
    cols = np.concatenate([MASK[Z[j]], MASK[Z[j]] + NAO]) if spinful else MASK[Z[j]]
    m = rng.normal(size=(len(rows), len(cols))) + (0 if real else 1j * rng.normal(size=(len(rows), len(cols))))
    out[np.ix_(rows, cols)] = scale * m
    return out


def _structure(rng, spinful, collinear=False):
    """Hermitian blocks ``(onsite, offsite)`` and the overlap of one structure."""
    edges = _edges()
    index = {(i, j, R): e for e, (i, j, R) in enumerate(edges)}
    on, off, son, soff = [], [None] * len(edges), [], [None] * len(edges)
    real = not spinful or collinear
    for a in range(2):
        m = _padded(rng, a, a, spinful, real=real)
        on.append(0.5 * (m + m.conj().T))
        s = 0.05 * _padded(rng, a, a, False, real=True)
        s = 0.5 * (s + s.T)
        s[MASK[Z[a]], MASK[Z[a]]] += 1.0
        son.append(s.real)
    for e, (i, j, R) in enumerate(edges):
        partner = index[(j, i, tuple(-x for x in R))]
        if off[partner] is not None:
            off[e], soff[e] = off[partner].conj().T, soff[partner].T
        else:
            off[e] = 0.3 * _padded(rng, i, j, spinful, real=real)
            soff[e] = 0.05 * _padded(rng, i, j, False, real=True).real
    if collinear:   # no spin mixing: zero up-down blocks
        for blocks in (on, off):
            for b in blocks:
                b[:NAO, NAO:] = 0
                b[NAO:, :NAO] = 0
    return edges, np.array(on), np.array(off), np.array(son), np.array(soff)


def _graph(edges, son, soff):
    edge_index = np.array([[e[0] for e in edges], [e[1] for e in edges]])
    cell_shift = np.array([e[2] for e in edges])
    return SimpleNamespace(z=Z.copy(), cell=LATTICE[None].copy(), pos=POSITIONS.copy(), edge_index=edge_index,
                           cell_shift=cell_shift, nbr_shift=cell_shift @ LATTICE,
                           Son=son.reshape(2, -1), Soff=soff.reshape(len(edges), -1))


def _hamgnn_convention(blocks):
    """HamGNN stores the up-down block as -ud (so du = -ud^+), see models/hamgnn.py."""
    out = np.array(blocks)
    out[:, :NAO, NAO:] *= -1
    return out


def _flat(on, off, spin):
    """The prediction rows of one structure in HamGNN's output layout."""
    blocks = np.concatenate([on, off])
    if spin == "soc":
        blocks = _hamgnn_convention(blocks)
    if spin == "none":
        return blocks.real.reshape(len(blocks), -1)
    if spin == "collinear":
        both = np.stack([blocks[:, :NAO, :NAO], blocks[:, NAO:, NAO:]], axis=1)
        return both.real.reshape(len(blocks), -1)
    return np.concatenate([blocks.real, blocks.imag]).reshape(2 * len(blocks), -1)


def _write_fixture(tmp_path, spin, seed=0):
    """Two structures in graph_data.npz and the prediction file; return the blocks of the second."""
    rng = np.random.default_rng(seed)
    spinful = spin != "none"
    graphs, rows = {}, []
    for n in range(2):
        edges, on, off, son, soff = _structure(rng, spinful, collinear=spin == "collinear")
        graphs[n] = _graph(edges, son, soff)
        rows.append(_flat(on, off, spin))
    np.savez(tmp_path / "graph_data.npz", graph=np.array(graphs, dtype=object))
    np.save(tmp_path / "prediction_hamiltonian.npy", np.concatenate(rows).ravel())
    return edges, on, off, son, soff


def _band_cal_sum(edges, on, off, son, soff, spinful):
    """HamGNN's band_cal.py: padded blocks, exp(+2 pi i k.R), masked, then [[uu, ud], [du, dd]]."""
    keep = np.concatenate([MASK[z] + a * NAO for a, z in enumerate(Z)])

    def bloch(onsite, offsite, k):
        out = np.zeros((2 * NAO, 2 * NAO), complex)
        for a in range(2):
            out[a * NAO:(a + 1) * NAO, a * NAO:(a + 1) * NAO] += onsite[a]
        for e, (i, j, R) in enumerate(edges):
            out[i * NAO:(i + 1) * NAO, j * NAO:(j + 1) * NAO] += np.exp(2j * np.pi * k @ np.array(R)) * offsite[e]
        return out[np.ix_(keep, keep)]

    hk, sk = [], []
    for k in K:
        s = bloch(son, soff, k)
        if spinful:
            parts = [[bloch(on[:, a * NAO:(a + 1) * NAO, b * NAO:(b + 1) * NAO],
                            off[:, a * NAO:(a + 1) * NAO, b * NAO:(b + 1) * NAO], k) for b in range(2)]
                     for a in range(2)]
            hk.append(np.block(parts))
            sk.append(np.kron(np.eye(2), s))
        else:
            hk.append(bloch(on, off, k))
            sk.append(s)
    return np.array(hk), np.array(sk)


def _spin_major(matrix):
    n = NORB.sum()
    perm = np.concatenate([2 * np.arange(n), 2 * np.arange(n) + 1])
    return matrix[:, perm][:, :, perm]


@pytest.mark.parametrize("spin", ["none", "collinear", "soc"])
def test_export_matches_band_cal(tmp_path, spin):
    edges, on, off, son, soff = _write_fixture(tmp_path, spin)
    out = tmp_path / "model.npz"
    write_hamgnn_model(tmp_path / "graph_data.npz", tmp_path / "prediction_hamiltonian.npy", out,
                       nao_max=NAO, spin=spin, structure=1)
    model = load_hamgnn_model(out, fermi_energy=-1.5)
    spinful = spin != "none"
    assert model.spinful is spinful and model.metadata["orthogonal"] is False
    np.testing.assert_allclose(model.uc, LATTICE * BOHR_ANGSTROM, atol=1e-12)
    h_ref, s_ref = _band_cal_sum(edges, on, off, son, soff, spinful)
    h_ref = h_ref * HARTREE_EV + 1.5 * s_ref
    hk = cell_gauge_hk(model, K)
    s_cart = K @ (2 * np.pi * np.linalg.inv(model.uc).T)
    p = np.exp(1j * s_cart @ model.tau.T)
    sk = p[:, :, None] * model.Sk(s_cart) * p.conj()[:, None, :]   # atom gauge back to the cell gauge
    if spinful:
        hk, sk = _spin_major(hk), _spin_major(sk)
    np.testing.assert_allclose(hk, h_ref, atol=1e-10)
    np.testing.assert_allclose(sk, s_ref, atol=1e-12)
    assert model.metadata["hamiltonian_hermiticity_residual"] < 1e-12


def test_command_line_and_dft_blocks(tmp_path):
    edges, on, off, son, soff = _write_fixture(tmp_path, "soc", seed=1)
    graphs = np.load(tmp_path / "graph_data.npz", allow_pickle=True)["graph"].item()
    flat = _flat(on, off, "soc")
    half = len(flat) // 2
    g = graphs[1]
    g.Hon, g.Hoff = flat[:2], flat[2:half]
    g.iHon, g.iHoff = flat[half:half + 2], flat[half + 2:]
    np.savez(tmp_path / "graph_dft.npz", graph=np.array(graphs, dtype=object))
    main([str(tmp_path / "graph_dft.npz"), "dft", str(tmp_path / "dft.npz"), "--nao-max", str(NAO),
          "--spin", "soc", "--structure", "1"])
    main([str(tmp_path / "graph_data.npz"), str(tmp_path / "prediction_hamiltonian.npy"),
          str(tmp_path / "pred.npz"), "--nao-max", str(NAO), "--spin", "soc", "--structure", "1"])
    dft, pred = load_hamgnn_model(tmp_path / "dft.npz"), load_hamgnn_model(tmp_path / "pred.npz")
    assert dft.metadata["hamiltonian_source"] == "dft" and pred.metadata["hamiltonian_source"] == "prediction"
    np.testing.assert_allclose(cell_gauge_hk(dft, K), cell_gauge_hk(pred, K), atol=1e-12)


def test_wrong_layout_is_refused(tmp_path):
    _write_fixture(tmp_path, "soc")
    args = (tmp_path / "graph_data.npz", tmp_path / "prediction_hamiltonian.npy", tmp_path / "m.npz")
    with pytest.raises(ValueError, match="nao_max"):
        write_hamgnn_model(*args, nao_max=NAO + 1, spin="soc")
    with pytest.raises(ValueError, match="spin mode"):
        write_hamgnn_model(*args, nao_max=NAO, spin="none")
    with pytest.raises(ValueError, match="out of range"):
        write_hamgnn_model(*args, nao_max=NAO, spin="soc", structure=2)
    np.savez(tmp_path / "other.npz", a=np.zeros(3))
    with pytest.raises(ValueError, match="not a HamGNN export"):
        load_hamgnn_model(tmp_path / "other.npz")


def test_basis_def_overrides_overlap_mask(tmp_path):
    _write_fixture(tmp_path, "none")
    args = (tmp_path / "graph_data.npz", tmp_path / "prediction_hamiltonian.npy")
    write_hamgnn_model(*args, tmp_path / "a.npz", nao_max=NAO, structure=1)
    write_hamgnn_model(*args, tmp_path / "b.npz", nao_max=NAO, structure=1, basis_def=MASK)
    a, b = np.load(tmp_path / "a.npz"), np.load(tmp_path / "b.npz")
    np.testing.assert_array_equal(a["orbital_index"], b["orbital_index"])
    np.testing.assert_array_equal(a["orbitals_per_atom"], NORB)
    assert str(a["dft_code"]) == "unknown" and (a["orbital_l"] == -1).all()


def test_dft_code_records_the_basis(tmp_path, monkeypatch):
    _write_fixture(tmp_path, "none")
    from lsquant.models.hamgnn import PADDED_SHELLS, padded_l
    for (code, nao), n in (("openmx", 14), 14), (("openmx", 26), 26), (("abacus", 27), 27), (("abacus", 40), 40):
        assert len(padded_l(code, nao)) == n
    ell = padded_l("openmx", 26)
    np.testing.assert_array_equal(ell[[0, 2, 3, 8, 9, 18, 19, 25]], [0, 0, 1, 1, 2, 2, 3, 3])
    # No HamGNN table has nao_max = 6: l is unknown (-1) ...
    args = (tmp_path / "graph_data.npz", tmp_path / "prediction_hamiltonian.npy")
    write_hamgnn_model(*args, tmp_path / "a.npz", nao_max=NAO, structure=1, dft_code="openmx")
    model = load_hamgnn_model(tmp_path / "a.npz")
    assert model.metadata["dft_code"] == "openmx" and set(model.metadata["orbital_l"]) == {-1}
    # ... and with a table s, s, s, p for the 6 padded orbitals, Se (s + p) and H (s) get their l.
    monkeypatch.setitem(PADDED_SHELLS, ("openmx", NAO), "sssp")
    write_hamgnn_model(*args, tmp_path / "c.npz", nao_max=NAO, structure=1, dft_code="openmx")
    assert load_hamgnn_model(tmp_path / "c.npz").metadata["orbital_l"] == [0, 1, 1, 1, 0]
    with pytest.raises(ValueError, match="dft_code"):
        write_hamgnn_model(tmp_path / "graph_data.npz", None, tmp_path / "b.npz", nao_max=NAO, dft_code="vasp")


def test_graph_without_overlap_gives_orthogonal_model(tmp_path):
    edges, on, off, son, soff = _write_fixture(tmp_path, "none")
    graphs = np.load(tmp_path / "graph_data.npz", allow_pickle=True)["graph"].item()
    for g in graphs.values():
        del g.Son, g.Soff
    np.savez(tmp_path / "graph_noS.npz", graph=np.array(graphs, dtype=object))
    args = (tmp_path / "graph_noS.npz", tmp_path / "prediction_hamiltonian.npy", tmp_path / "m.npz")
    with pytest.raises(ValueError, match="basis_def"):
        write_hamgnn_model(*args, nao_max=NAO, structure=1)
    write_hamgnn_model(*args, nao_max=NAO, structure=1, basis_def=MASK)
    model = load_hamgnn_model(tmp_path / "m.npz", fermi_energy=0.4)
    assert type(model).__name__ == "OrthogonalModel" and model.spinful is False
    h_ref, _ = _band_cal_sum(edges, on, off, son, soff, False)
    np.testing.assert_allclose(cell_gauge_hk(model, K), h_ref * HARTREE_EV - 0.4 * np.eye(NORB.sum()), atol=1e-10)
    basis = BlochMetricBasis(model, [3, 3, 1], velocity_mode="bare")
    assert basis.operator_blocks("VX").shape[-1] == NORB.sum()


def test_spin_operator_and_spinless_refusal(tmp_path):
    _write_fixture(tmp_path, "soc")
    write_hamgnn_model(tmp_path / "graph_data.npz", tmp_path / "prediction_hamiltonian.npy",
                       tmp_path / "soc.npz", nao_max=NAO, spin="soc")
    basis = BlochMetricBasis(load_hamgnn_model(tmp_path / "soc.npz"), [3, 3, 1], velocity_mode="berry_connection")
    sigma_z = np.kron(np.eye(NORB.sum()), np.diag([1.0, -1.0]))
    np.testing.assert_allclose(basis.operator_blocks("SZ"), basis.sk_array @ sigma_z, atol=1e-12)
    spinless_dir = tmp_path / "nosoc"
    spinless_dir.mkdir()
    _write_fixture(spinless_dir, "none")
    write_hamgnn_model(spinless_dir / "graph_data.npz", spinless_dir / "prediction_hamiltonian.npy",
                       spinless_dir / "m.npz", nao_max=NAO)
    spinless = BlochMetricBasis(load_hamgnn_model(spinless_dir / "m.npz"), [3, 3, 1], velocity_mode="bare")
    with pytest.raises(ValueError, match="spinless"):
        spinless.operator_blocks("SZ")


def test_runfile_routes_and_runs_hamgnn(tmp_path):
    _write_fixture(tmp_path, "soc")
    path = tmp_path / "model.npz"
    write_hamgnn_model(tmp_path / "graph_data.npz", tmp_path / "prediction_hamiltonian.npy", path,
                       nao_max=NAO, spin="soc")
    model = {"kind": "hamgnn", "file": str(path), "params": {"fermi_energy": 0.2}}
    assert uses_metric_executor(model)
    assert load_runfile_model(model, "berry_connection").metadata["fermi_energy_subtracted_eV"] == 0.2
    base = dict(word="Tr[delta(E-H)]", expansion={"moments": 8}, basis={"name": "eigen", "supercell": [2, 2, 1]})
    assert RunFile.from_dict(dict(base, model=model)).model["kind"] == "hamgnn"
    with pytest.raises(ValueError, match="requires 'file'"):
        RunFile.from_dict(dict(base, model={"kind": "hamgnn"}))
    with pytest.raises(ValueError, match="covariant"):
        RunFile.from_dict(dict(base, model=dict(model, params={"velocity_mode": "covariant"})))
    with pytest.raises(ValueError):
        RunFile.from_dict(dict(base, model=dict(model, params={"nao_max": 26})))
    runs = (("Tr[delta(E-H)]", "dos", None), ("Tr[JXSZ delta(E-H) VY dG-(E)/dE]", "shc", None),
            ("Tr[JXSZ delta(E-H) VY dG-(E)/dE]", "shc", {"kind": "shift", "window": [-20.0, 20.0]}))
    for n, (word, formula, projection) in enumerate(runs):
        out = tmp_path / f"{formula}{n}"
        execute_runfile({
            "model": {"kind": "hamgnn", "file": str(path), "params": {"velocity_mode": "berry_connection"}},
            "basis": {"name": "eigen", "supercell": [3, 3, 1]},
            "word": word,
            "expansion": {"moments": 32, "window": "eigen", "trace": "exact"},
            "reconstruction": {"kernel": "jackson", "formulas": [formula], "broadening_meV": 2000.0,
                               "energy_grid": [-20.0, 20.0, 21]},
            "output": {"dir": str(out)},
            **({"projection": projection} if projection else {}),
        })
        resolved = json.loads((out / "run.resolved.json").read_text())
        assert resolved["model"]["kind"] == "hamgnn"
        if projection is None:
            assert json.loads((out / "basis.json").read_text())["source"] == "hamgnn"
