"""DeepH-format Hamiltonians (DeepH-pack, DeepH-E3, xDeepH, MACE-H) in the metric executors.

The fixtures are written here, in the two layouts of the DeepH codes, from
a random two-atom model with spin-orbit coupling. The checks:
- H(k) and S(k) of the LSQuant model equal the sums of the DeepH tools
  (``dense_calc.py`` of DeepH-pack for the processed-structure layout,
  ``matrix_obj.py`` of DeepH-dock for the DeepH-pack layout), after the
  spinor permutation and the Fermi-energy shift;
- the split storage ``H = I_2 tensor H_0 + H_KB`` is rebuilt;
- the spin operator of the metric basis is ``S (I tensor sigma)``;
- a spinless model refuses spin words; an asymmetric prediction is
  symmetrized with a warning;
- the run file routes ``model.kind = "deeph"`` and checks its params.
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
h5py = pytest.importorskip("h5py")
from lsquant.bases import BlochMetricBasis  # noqa: E402
from lsquant.models.deeph import load_deeph_model  # noqa: E402
from lsquant.models.lcao import cell_gauge_hk, orbital_slices  # noqa: E402
from lsquant.models.orthogonal import load_runfile_model  # noqa: E402
from lsquant.runfile import RunFile, uses_metric_executor  # noqa: E402

LATTICE = np.array([[3.0, 0.0, 0.0], [1.4, 2.6, 0.0], [0.0, 0.0, 15.0]])
POSITIONS = np.array([[0.1, 0.2, 7.0], [1.6, 1.1, 7.4]])
ORBITAL_TYPES = [[0, 1], [0]]          # atom 0: s + p (4 orbitals), atom 1: s (1 orbital)
NORB = np.array([4, 1])
CELLS = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, -1, 0), (-1, 1, 0)]


def _pair_blocks(rng, spin, scale=1.0, hermitian=True, complex_=True):
    """Random blocks ``M_ij(R)`` with ``M_ji(-R) = M_ij(R)^+`` (``spin`` = 1 or 2)."""
    blocks = {}
    for R in CELLS:
        for i in range(2):
            for j in range(2):
                key = (R, i, j)
                partner = (tuple(-x for x in R), j, i)
                if partner in blocks and hermitian:
                    blocks[key] = blocks[partner].conj().T
                    continue
                shape = (spin * NORB[i], spin * NORB[j])
                m = rng.normal(size=shape) + (1j * rng.normal(size=shape) if complex_ else 0)
                if key == partner:
                    m = 0.5 * (m + m.conj().T)
                blocks[key] = scale * m
    return blocks


def _overlap_blocks(rng):
    blocks = _pair_blocks(rng, 1, scale=0.05, complex_=False)
    for i in range(2):
        blocks[((0, 0, 0), i, i)] = np.eye(NORB[i]) + blocks[((0, 0, 0), i, i)]
    return blocks


def _deeph_sum(blocks, k_frac, spin):
    """``sum_R exp(+2 pi i k.R) M(R)`` in DeepH order: atom by atom, each atom up then down."""
    starts = np.concatenate([[0], np.cumsum(spin * NORB)])
    out = np.zeros((len(k_frac), starts[-1], starts[-1]), complex)
    for (R, i, j), m in blocks.items():
        phase = np.exp(2j * np.pi * (k_frac @ np.array(R)))
        out[:, starts[i]:starts[i + 1], starts[j]:starts[j + 1]] += phase[:, None, None] * m
    return out


def _write_processed(path, h_blocks, s_blocks, spinful, fermi):
    """The processed-structure layout of DeepH-pack, DeepH-E3, xDeepH and MACE-H."""
    np.savetxt(path / "lat.dat", LATTICE.T)
    np.savetxt(path / "site_positions.dat", POSITIONS.T)
    (path / "orbital_types.dat").write_text("\n".join(" ".join(map(str, o)) for o in ORBITAL_TYPES) + "\n")
    (path / "info.json").write_text(json.dumps({"isspinful": spinful, "fermi_level": fermi}))
    for name, blocks in (("hamiltonians_pred.h5", h_blocks), ("overlaps.h5", s_blocks)):
        with h5py.File(path / name, "w") as f:
            for (R, i, j), m in blocks.items():
                f[str([R[0], R[1], R[2], i + 1, j + 1])] = m


def _write_chunked(path, blocks):
    keys = list(blocks)
    entries = np.concatenate([blocks[k].ravel() for k in keys])
    with h5py.File(path, "w") as f:
        f["atom_pairs"] = np.array([[*R, i, j] for R, i, j in keys])
        f["chunk_boundaries"] = np.concatenate([[0], np.cumsum([blocks[k].size for k in keys])])
        f["chunk_shapes"] = np.array([blocks[k].shape for k in keys])
        f["entries"] = entries


def _write_pack(path, fermi, spinful, split=False):
    """The DeepH-pack layout: POSCAR, info.json, chunked HDF5 files."""
    frac = POSITIONS @ np.linalg.inv(LATTICE)
    lines = ["toy", "1.0"] + [" ".join(f"{x:.12f}" for x in v) for v in LATTICE]
    lines += ["Pd Se", "1 1", "Direct"] + [" ".join(f"{x:.12f}" for x in p) for p in frac]
    (path / "POSCAR").write_text("\n".join(lines) + "\n")
    info = {"atoms_quantity": 2, "orbits_quantity": 5, "orthogonal_basis": False, "spinful": spinful,
            "fermi_energy_eV": fermi, "elements_orbital_map": {"Pd": [0, 1], "Se": [0]}}
    if split:
        info["split_hkb"] = True
    (path / "info.json").write_text(json.dumps(info))


def _to_deeph_order(model_matrix, spinful):
    rows = np.concatenate(orbital_slices(NORB, spinful))
    return model_matrix[:, rows][:, :, rows]


K = np.array([[0.13, 0.27, 0.0], [-0.31, 0.08, 0.0], [0.45, -0.22, 0.0]])


def test_processed_structure_layout_matches_dense_calc(tmp_path):
    rng = np.random.default_rng(1)
    h, s = _pair_blocks(rng, 2), _overlap_blocks(rng)
    _write_processed(tmp_path, h, s, True, 1.3)
    model = load_deeph_model(tmp_path)
    assert model.metadata["layout"] == "processed-structure" and model.spinful
    s_ref = _deeph_sum({k: np.kron(np.eye(2), m) for k, m in s.items()}, K, 2)
    h_ref = _deeph_sum(h, K, 2) - 1.3 * s_ref
    np.testing.assert_allclose(_to_deeph_order(cell_gauge_hk(model, K), True), h_ref, atol=1e-12)
    k_cart = K @ (2 * np.pi * np.linalg.inv(LATTICE).T)
    hk = model.Hk(k_cart)
    np.testing.assert_allclose(hk, np.conj(np.swapaxes(hk, 1, 2)), atol=1e-12)
    # Atom gauge: the same matrix up to the orbital phases exp(i k.tau).
    p = np.exp(1j * k_cart @ model.tau.T)
    np.testing.assert_allclose(hk, p.conj()[:, :, None] * cell_gauge_hk(model, K) * p[:, None, :], atol=1e-12)
    # An inference directory may hold rlat.dat (reciprocal vectors as columns) instead of lat.dat.
    np.savetxt(tmp_path / "rlat.dat", 2 * np.pi * np.linalg.inv(LATTICE))
    (tmp_path / "lat.dat").unlink()
    np.testing.assert_allclose(load_deeph_model(tmp_path).uc, LATTICE, atol=1e-12)


def test_deeph_pack_layout_and_split_kb_storage(tmp_path):
    rng = np.random.default_rng(2)
    h0, kb, s = _pair_blocks(rng, 1, complex_=False), _pair_blocks(rng, 2, scale=0.2), _overlap_blocks(rng)
    _write_pack(tmp_path, -2.0, True, split=True)
    _write_chunked(tmp_path / "hamiltonian_pred.h5", h0)
    _write_chunked(tmp_path / "hkb.h5", kb)
    _write_chunked(tmp_path / "overlap.h5", s)
    model = load_deeph_model(tmp_path)
    assert model.metadata["layout"] == "deeph-pack"
    # DeepH-dock assembles all up components, then all down ones.
    n = NORB.sum()
    starts = np.concatenate([[0], np.cumsum(NORB)])

    def global_blocked(blocks, spinful):
        out = np.zeros((len(K), 2 * n, 2 * n), complex)
        for (R, i, j), m in blocks.items():
            phase = np.exp(2j * np.pi * (K @ np.array(R)))[:, None, None]
            full = m if spinful else np.kron(np.eye(2), m)
            for si in range(2):
                for sj in range(2):
                    sub = full[si * NORB[i]:(si + 1) * NORB[i], sj * NORB[j]:(sj + 1) * NORB[j]]
                    rows = slice(starts[i] + si * n, starts[i + 1] + si * n)
                    cols = slice(starts[j] + sj * n, starts[j + 1] + sj * n)
                    out[:, rows, cols] += phase * sub
        return out

    s_ref = global_blocked(s, False)
    h_ref = global_blocked(h0, False) + global_blocked(kb, True) + 2.0 * s_ref
    perm = np.concatenate([2 * np.arange(n), 2 * np.arange(n) + 1])
    model_h = cell_gauge_hk(model, K)[:, perm][:, :, perm]
    np.testing.assert_allclose(model_h, h_ref, atol=1e-12)


def test_spin_operator_and_spinless_refusal(tmp_path):
    rng = np.random.default_rng(3)
    spinful_dir, spinless_dir = tmp_path / "soc", tmp_path / "nosoc"
    spinful_dir.mkdir(), spinless_dir.mkdir()
    s = _overlap_blocks(rng)
    _write_processed(spinful_dir, _pair_blocks(rng, 2), s, True, 0.0)
    _write_processed(spinless_dir, _pair_blocks(rng, 1, complex_=False), s, False, 0.0)
    basis = BlochMetricBasis(load_deeph_model(spinful_dir), [3, 3, 1], velocity_mode="berry_connection")
    sigma_z = np.kron(np.eye(NORB.sum()), np.diag([1.0, -1.0]))
    np.testing.assert_allclose(basis.operator_blocks("SZ"), basis.sk_array @ sigma_z, atol=1e-12)
    spinless = BlochMetricBasis(load_deeph_model(spinless_dir), [3, 3, 1], velocity_mode="bare")
    assert spinless.operator_blocks("VX").shape[-1] == NORB.sum()
    with pytest.raises(ValueError, match="spinless"):
        spinless.operator_blocks("SZ")
    with pytest.raises(ValueError, match="spinless"):
        spinless.operator_blocks("JXSZ")


def test_asymmetric_prediction_is_symmetrized(tmp_path):
    rng = np.random.default_rng(4)
    h = _pair_blocks(rng, 1, complex_=False)
    h[((1, 0, 0), 0, 1)] = h[((1, 0, 0), 0, 1)] + 0.1
    _write_processed(tmp_path, h, _overlap_blocks(rng), False, 0.0)
    with pytest.warns(RuntimeWarning, match="differs from H"):
        model = load_deeph_model(tmp_path)
    assert model.metadata["hamiltonian_hermiticity_residual"] > 1e-3
    hk = model.Hk(K @ (2 * np.pi * np.linalg.inv(LATTICE).T))
    np.testing.assert_allclose(hk, np.conj(np.swapaxes(hk, 1, 2)), atol=1e-12)


def test_runfile_routes_deeph(tmp_path):
    rng = np.random.default_rng(5)
    _write_processed(tmp_path, _pair_blocks(rng, 2), _overlap_blocks(rng), True, 0.7)
    model = {"kind": "deeph", "dir": str(tmp_path), "params": {"fermi_energy": 0.2}}
    assert uses_metric_executor(model)
    loaded = load_runfile_model(model, "berry_connection")
    assert loaded.metadata["fermi_energy_subtracted_eV"] == 0.2
    base = dict(word="Tr[JXSZ delta(E-H) VY dG-(E)/dE]", expansion={"moments": 8},
                basis={"name": "eigen", "supercell": [2, 2, 1]})
    rf = RunFile.from_dict(dict(base, model=model))
    assert rf.model["kind"] == "deeph"
    with pytest.raises(ValueError, match="requires 'dir'"):
        RunFile.from_dict(dict(base, model={"kind": "deeph"}))
    with pytest.raises(ValueError, match="covariant"):
        RunFile.from_dict(dict(base, model=dict(model, params={"velocity_mode": "covariant"})))
    with pytest.raises(ValueError, match="spin"):
        RunFile.from_dict(dict(base, model=dict(model, params={"spin": "exact"})))


def test_runfile_executes_deeph_dos_and_spin_hall(tmp_path):
    rng = np.random.default_rng(6)
    model_dir = tmp_path / "model"
    model_dir.mkdir()
    _write_processed(model_dir, _pair_blocks(rng, 2, scale=0.3), _overlap_blocks(rng), True, 0.0)
    from lsquant.runfile import execute_runfile
    runs = (("Tr[delta(E-H)]", "dos", None), ("Tr[JXSZ delta(E-H) VY dG-(E)/dE]", "shc", None),
            ("Tr[JXSZ delta(E-H) VY dG-(E)/dE]", "shc", {"kind": "shift", "window": [-1.0, 1.0]}))
    for n, (word, formula, projection) in enumerate(runs):
        out = tmp_path / f"{formula}{n}"
        execute_runfile({
            "model": {"kind": "deeph", "dir": str(model_dir), "params": {"velocity_mode": "berry_connection"}},
            "basis": {"name": "eigen", "supercell": [3, 3, 1]},
            "word": word,
            "expansion": {"moments": 32, "window": "eigen", "trace": "exact"},
            "reconstruction": {"kernel": "jackson", "formulas": [formula], "broadening_meV": 200.0,
                               "energy_grid": [-1.0, 1.0, 21]},
            "output": {"dir": str(out)},
            **({"projection": projection} if projection else {}),
        })
        resolved = json.loads((out / "run.resolved.json").read_text())
        assert resolved["model"]["kind"] == "deeph"
        if projection is None:
            assert json.loads((out / "basis.json").read_text())["source"] == "deeph"
