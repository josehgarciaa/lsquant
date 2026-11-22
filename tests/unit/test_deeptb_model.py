"""DeePTB Hamiltonians (DeePTB-SK and DeePTB-E3) from DeePTB's sparse export.

The export directories are written here, as ``TBSystem.to_pardiso_json``
and ``TBSystem.to_pardiso`` write them, from a random spinless model. The
checks:
- H(k) and S(k) equal the sum of DeePTB's Julia backend
  (``exp(+2 pi i k.R)``, block ``"i_j_R"`` in rows of atom ``i``) with the
  Fermi energy removed; both structure formats give the same model;
- without ``predicted_overlaps.h5`` the model is orthogonal;
- a spinful export (spin-orbit term missing from the files) is refused;
- the run file routes ``model.kind = "deeptb"`` and runs a DOS.
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
from lsquant.models.deeptb import load_deeptb_model  # noqa: E402
from lsquant.models.lcao import cell_gauge_hk  # noqa: E402
from lsquant.runfile import RunFile, execute_runfile, uses_metric_executor  # noqa: E402

CELL = np.array([[3.8, 0.0, 0.0], [0.0, 3.8, 0.0], [0.0, 0.0, 20.0]])
POSITIONS = np.array([[0.0, 0.0, 10.0], [1.9, 1.9, 10.3]])
SYMBOLS = ["Mo", "S"]
COUNTS = {"Mo": 9, "S": 4}            # Mo: 1s1p1d, S: 1s1p
NORB = np.array([9, 4])
CELLS = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (-1, 1, 0)]
K = np.array([[0.1, 0.2, 0.0], [0.37, -0.15, 0.0], [-0.42, 0.33, 0.0]])


def _blocks(rng, scale, onsite_identity=False):
    """Canonical blocks (i <= j) of a real model, completed with partners as DeePTB writes them."""
    blocks = {}
    for R in CELLS:
        for i in range(2):
            for j in range(2):
                key, partner = (R, i, j), (tuple(-x for x in R), j, i)
                if key in blocks or partner in blocks:
                    continue
                m = scale * rng.normal(size=(NORB[i], NORB[j]))
                if key == partner:
                    m = 0.5 * (m + m.T) + (np.eye(NORB[i]) if onsite_identity else 0)
                blocks[key] = m
                blocks[partner] = m.T
    return blocks


def _write_h5(path, blocks):
    with h5py.File(path, "w") as f:
        group = f.create_group("0")
        for (R, i, j), m in blocks.items():
            group.create_dataset(f"{i}_{j}_{R[0]}_{R[1]}_{R[2]}", data=m)


def _write_json_structure(path, spinful=False):
    data = {"basis_info": {"basis": {"Mo": ["1s", "1p", "1d"], "S": ["1s", "1p"]}, "orbital_counts": COUNTS,
                           "site_norbits": [int(n) * (2 if spinful else 1) for n in NORB],
                           "total_orbitals": int(NORB.sum()) * (2 if spinful else 1), "spinful": spinful},
            "structure": {"cell": CELL.tolist(), "pbc": [True, True, False], "nsites": 2,
                          "symbols": SYMBOLS, "positions": POSITIONS.tolist()},
            "meta": {"version": "1.0", "generator": "test"}}
    (path / "structure.json").write_text(json.dumps(data))


def _julia_sum(blocks, k_frac):
    starts = np.concatenate([[0], np.cumsum(NORB)])
    out = np.zeros((len(k_frac), starts[-1], starts[-1]), complex)
    for (R, i, j), m in blocks.items():
        out[:, starts[i]:starts[i + 1], starts[j]:starts[j + 1]] += \
            np.exp(2j * np.pi * (k_frac @ np.array(R)))[:, None, None] * m
    return out


def test_export_matches_julia_backend_and_both_structure_formats(tmp_path):
    rng = np.random.default_rng(11)
    h, s = _blocks(rng, 1.0), _blocks(rng, 0.05, onsite_identity=True)
    new, old = tmp_path / "json", tmp_path / "dat"
    for d in (new, old):
        d.mkdir()
        _write_h5(d / "predicted_hamiltonians.h5", h)
        _write_h5(d / "predicted_overlaps.h5", s)
    _write_json_structure(new)
    np.savetxt(old / "cell.dat", CELL)
    np.savetxt(old / "positions.dat", POSITIONS)
    np.savetxt(old / "atomic_numbers.dat", [42, 16], fmt="%d")
    (old / "basis.dat").write_text(str({"Mo": "1s1p1d", "S": "1s1p"}))
    model = load_deeptb_model(new, fermi_energy=-1.5)
    assert model.metadata["source"] == "deeptb" and model.metadata["orthogonal"] is False
    reference = _julia_sum(h, K) + 1.5 * _julia_sum(s, K)
    np.testing.assert_allclose(cell_gauge_hk(model, K), reference, atol=1e-12)
    assert model.Sk(K @ (2 * np.pi * np.linalg.inv(CELL).T)).shape == (3, 13, 13)
    legacy = load_deeptb_model(old, fermi_energy=-1.5)
    np.testing.assert_allclose(cell_gauge_hk(legacy, K), reference, atol=1e-12)
    np.testing.assert_allclose(legacy.tau, model.tau)


def test_orthogonal_export_and_spinful_refusal(tmp_path):
    rng = np.random.default_rng(12)
    h = _blocks(rng, 1.0)
    _write_h5(tmp_path / "predicted_hamiltonians.h5", h)
    _write_json_structure(tmp_path)
    model = load_deeptb_model(tmp_path)
    assert model.metadata["orthogonal"] is True and model.spinful is False
    np.testing.assert_allclose(cell_gauge_hk(model, K), _julia_sum(h, K), atol=1e-12)
    _write_json_structure(tmp_path, spinful=True)
    with pytest.raises(ValueError, match="spin-orbit term is not in the files"):
        load_deeptb_model(tmp_path)


def test_runfile_routes_deeptb(tmp_path):
    rng = np.random.default_rng(13)
    model_dir = tmp_path / "model"
    model_dir.mkdir()
    _write_h5(model_dir / "predicted_hamiltonians.h5", _blocks(rng, 0.5))
    _write_h5(model_dir / "predicted_overlaps.h5", _blocks(rng, 0.05, onsite_identity=True))
    _write_json_structure(model_dir)
    model = {"kind": "deeptb", "dir": str(model_dir), "params": {"fermi_energy": 0.3}}
    assert uses_metric_executor(model)
    base = dict(word="Tr[delta(E-H)]", expansion={"moments": 32, "window": "eigen", "trace": "exact"},
                basis={"name": "eigen", "supercell": [3, 3, 1]})
    with pytest.raises(ValueError, match="requires 'dir'"):
        RunFile.from_dict(dict(base, model={"kind": "deeptb"}))
    with pytest.raises(ValueError, match="hamiltonian_file"):
        RunFile.from_dict(dict(base, model=dict(model, params={"hamiltonian_file": "x.h5"})))
    out = tmp_path / "dos"
    execute_runfile(dict(base, model=model, output={"dir": str(out)},
                         reconstruction={"kernel": "jackson", "formulas": ["dos"], "broadening_meV": 200.0,
                                         "energy_grid": [-1.0, 1.0, 21]}))
    assert json.loads((out / "basis.json").read_text())["fermi_energy_subtracted_eV"] == 0.3
