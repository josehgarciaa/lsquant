"""ThreeBodyTB (TB3) models in the metric executors.

The fixtures are XML files written here in the layout of ThreeBodyTB's
``write_tb_crys`` (a two-atom C-H crystal with s+p and s orbitals). The
checks:
- H(k) and S(k) of the LSQuant model equal TB3's ``Hk``: the sum
  ``conj(H(r)) exp(-2 pi i k.R_r)``, its Hermitian part, plus the charge
  term ``S(k) .* h1``, in eV, with the Fermi energy of the file removed;
- a gzipped file reads the same; an orthogonal file gives an orthogonal model;
- two spin channels give a spinful model without spin mixing;
- ``scfspin``, a wrong orbital count and an unknown element are refused;
- the run file routes ``model.kind = "threebodytb"`` and runs the DOS.
"""
import gzip
import json
import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant.models.lcao import cell_gauge_hk  # noqa: E402
from lsquant.models.orthogonal import load_runfile_model  # noqa: E402
from lsquant.models.threebodytb import BOHR_ANGSTROM, RYDBERG_EV, load_threebodytb_model  # noqa: E402
from lsquant.runfile import RunFile, execute_runfile, uses_metric_executor  # noqa: E402

A = np.array([[4.6, 0.0, 0.0], [2.3, 4.0, 0.0], [0.0, 0.0, 30.0]])     # Bohr, rows
COORDS = np.array([[0.0, 0.0, 0.5], [1 / 3, 1 / 3, 0.52]])
TYPES = ["C", "H"]
NWAN = 5
CELLS = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (1, -1, 0), (-1, 1, 0)]
K = np.array([[0.13, 0.27, 0.0], [-0.31, 0.08, 0.0], [0.45, -0.22, 0.0]])
EFERMI = -0.21      # Rydberg


def _hermitian_set(rng, scale, complex_=True, identity=False):
    """``M[r]`` on CELLS with ``M(-R) = M(R)^+`` (TB3 stores both)."""
    mats = {}
    for R in CELLS:
        neg = tuple(-x for x in R)
        if neg in mats:
            mats[R] = mats[neg].conj().T
            continue
        m = scale * (rng.normal(size=(NWAN, NWAN)) + (1j * rng.normal(size=(NWAN, NWAN)) if complex_ else 0))
        if R == neg:
            m = 0.5 * (m + m.conj().T) + (np.eye(NWAN) if identity else 0)
        mats[R] = m
    return np.array([mats[R] for R in CELLS])


def _arr2str(a):
    return "\n".join(" ".join(f"{x: .16e}" for x in row) for row in np.atleast_2d(a)) + "\n"


def _write_xml(path, H, S, *, nonorth=True, h1=None, scfspin=False, gz=False):
    """A ``tb_crys`` XML in the layout of ThreeBodyTB's ``write_tb_crys``."""
    nspin = H.shape[0]
    lines = []
    for r in range(len(CELLS)):
        for m in range(NWAN):
            for n in range(NWAN):
                hs = " ".join(f"{H[s, r, m, n].real: .10f} {H[s, r, m, n].imag: .10f}" for s in range(nspin))
                sr = S[r, m, n] if nonorth else (1.0 if m == n and CELLS[r] == (0, 0, 0) else 0.0)
                lines.append(f"{r + 1:3d} {m + 1:3d} {n + 1:3d} {hs} {np.real(sr): .10f} {np.imag(sr): .10f} ")
    ind = "\n".join(" ".join(f"{x:7d}" for x in R) for R in CELLS)
    h1_text = _arr2str(h1) if h1 is not None else _arr2str(np.zeros((NWAN, NWAN)))
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<root>
  <crystal>
    <A>{_arr2str(A)}</A>
    <nat>2</nat>
    <coords>{_arr2str(COORDS)}</coords>
    <types>{"  ".join(TYPES)}</types>
  </crystal>
  <nelec>5.0</nelec>
  <nspin>{nspin}</nspin>
  <dftenergy>0.0</dftenergy>
  <scf>{str(h1 is not None).lower()}</scf>
  <efermi>{EFERMI}</efermi>
  <energy>-1.0</energy>
  <tightbinding>
    <nr>{len(CELLS)}</nr>
    <nwan>{NWAN}</nwan>
    <nspin>{nspin}</nspin>
    <nonorth>{str(nonorth).lower()}</nonorth>
    <ind_arr>{ind}</ind_arr>
    <r_dict>Dict()</r_dict>
    <scf>{str(h1 is not None).lower()}</scf>
    <h1>{h1_text}</h1>
    <scfspin>{str(scfspin).lower()}</scfspin>
    <H>{chr(10).join(lines)}
</H>
  </tightbinding>
</root>
"""
    data = xml.encode()
    path.write_bytes(gzip.compress(data) if gz else data)


def _tb3_hk(H, S, h1, k, spin=0, nonorth=True):
    """ThreeBodyTB's Hk: Hermitian part of sum conj(H) exp(-2 pi i k.R), plus S(k) .* h1."""
    phase = np.exp(-2j * np.pi * np.array(CELLS) @ k)
    hk = np.einsum("r,rmn->mn", phase, np.conj(H[spin]))
    sk = np.einsum("r,rmn->mn", phase, np.conj(S)) if nonorth else np.zeros_like(hk)
    hk = 0.5 * (hk + hk.conj().T)
    if h1 is not None:
        hk = hk + sk * h1
    return 0.5 * (hk + hk.conj().T), 0.5 * (sk + sk.conj().T)


def _fixture(seed=0, nspin=1):
    rng = np.random.default_rng(seed)
    H = np.array([_hermitian_set(rng, 0.1) for _ in range(nspin)])
    S = _hermitian_set(rng, 0.02, complex_=False, identity=True)
    h1 = rng.normal(size=(NWAN, NWAN)) * 0.05
    return H, S, 0.5 * (h1 + h1.T)


@pytest.mark.parametrize("gz", [False, True])
def test_matches_threebodytb_hk(tmp_path, gz):
    H, S, h1 = _fixture()
    path = tmp_path / ("model.xml.gz" if gz else "model.xml")
    _write_xml(path, H, S, h1=h1, gz=gz)
    model = load_threebodytb_model(path)
    assert model.spinful is False and model.metadata["orthogonal"] is False
    np.testing.assert_allclose(model.uc, A * BOHR_ANGSTROM, atol=1e-12)
    np.testing.assert_allclose(model.tau[:4], np.repeat([COORDS[0] @ A * BOHR_ANGSTROM], 4, axis=0), atol=1e-12)
    assert model.metadata["fermi_energy_subtracted_eV"] == pytest.approx(EFERMI * RYDBERG_EV)
    assert model.metadata["orbital_labels"] == [(0, "s"), (0, "pz"), (0, "px"), (0, "py"), (1, "s")]
    hk = cell_gauge_hk(model, K)
    for n, k in enumerate(K):
        h_ref, s_ref = _tb3_hk(H, S, h1, k)
        np.testing.assert_allclose(hk[n], (h_ref - EFERMI * s_ref) * RYDBERG_EV, atol=1e-8)
    # Generalized eigenvalues agree with TB3's (k in Cartesian, atom gauge).
    k_cart = K @ (2 * np.pi * np.linalg.inv(model.uc).T)
    from scipy.linalg import eigh
    h_ref, s_ref = _tb3_hk(H, S, h1, K[0])
    ours = eigh(model.Hk(k_cart)[0], model.Sk(k_cart)[0], eigvals_only=True)
    np.testing.assert_allclose(ours, eigh(h_ref, s_ref, eigvals_only=True) * RYDBERG_EV - EFERMI * RYDBERG_EV,
                               atol=1e-8)


def test_band_recipe_of_the_guide(tmp_path):
    """The bands of docs/how-to/use-threebodytb.md equal TB3's generalized eigenvalues on a path."""
    from scipy.linalg import eigh
    H, S, h1 = _fixture(5)
    _write_xml(tmp_path / "m.xml", H, S, h1=h1)
    model = load_threebodytb_model(tmp_path / "m.xml")
    nodes = np.array([[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0], [0, 0, 0]])
    path = np.concatenate([np.linspace(a, b, 5, endpoint=False) for a, b in zip(nodes[:-1], nodes[1:])]
                          + [nodes[-1:]])
    k_cart = path @ (2 * np.pi * np.linalg.inv(model.uc).T)
    bands = np.array([eigh(h, s, eigvals_only=True) for h, s in zip(model.Hk(k_cart), model.Sk(k_cart))])
    for k, ours in zip(path, bands):
        h_ref, s_ref = _tb3_hk(H, S, h1, k)
        np.testing.assert_allclose(ours, (eigh(h_ref, s_ref, eigvals_only=True) - EFERMI) * RYDBERG_EV, atol=1e-8)


def test_orthogonal_and_fermi_override(tmp_path):
    H, S, _ = _fixture(1)
    _write_xml(tmp_path / "m.xml", H, S, nonorth=False)
    model = load_threebodytb_model(tmp_path / "m.xml", fermi_energy=0.5)
    assert type(model).__name__ == "OrthogonalModel"
    hk = cell_gauge_hk(model, K)
    for n, k in enumerate(K):
        h_ref, _ = _tb3_hk(H, S, None, k, nonorth=False)
        np.testing.assert_allclose(hk[n], h_ref * RYDBERG_EV - 0.5 * np.eye(NWAN), atol=1e-8)


def test_two_spin_channels(tmp_path):
    H, S, h1 = _fixture(2, nspin=2)
    _write_xml(tmp_path / "m.xml", H, S, h1=h1)
    model = load_threebodytb_model(tmp_path / "m.xml", fermi_energy=0.0)
    assert model.spinful is True
    hk = cell_gauge_hk(model, K)
    for n, k in enumerate(K):
        for s in range(2):
            h_ref, _ = _tb3_hk(H, S, h1, k, spin=s)
            np.testing.assert_allclose(hk[n][s::2, s::2], h_ref * RYDBERG_EV, atol=1e-8)
        np.testing.assert_allclose(hk[n][0::2, 1::2], 0, atol=1e-12)


def test_refusals(tmp_path):
    H, S, h1 = _fixture(3)
    _write_xml(tmp_path / "spin.xml", H, S, h1=h1, scfspin=True)
    with pytest.raises(ValueError, match="h1spin"):
        load_threebodytb_model(tmp_path / "spin.xml")
    _write_xml(tmp_path / "ok.xml", H, S, h1=h1)
    with pytest.raises(ValueError, match="nwan"):
        load_threebodytb_model(tmp_path / "ok.xml", orbitals_per_element={"H": 4})
    text = (tmp_path / "ok.xml").read_text().replace("<types>C  H</types>", "<types>C  Xq</types>")
    (tmp_path / "unknown.xml").write_text(text)
    with pytest.raises(ValueError, match="orbitals_per_element"):
        load_threebodytb_model(tmp_path / "unknown.xml")
    unknown = load_threebodytb_model(tmp_path / "unknown.xml", orbitals_per_element={"Xq": 1})
    assert unknown.num_wann == NWAN and unknown.metadata["orbital_labels"] is None


def test_runfile_routes_and_runs_threebodytb(tmp_path):
    H, S, h1 = _fixture(4)
    path = tmp_path / "model.xml"
    _write_xml(path, H, S, h1=h1)
    model = {"kind": "threebodytb", "file": str(path), "params": {"fermi_energy": 0.1}}
    assert uses_metric_executor(model)
    assert load_runfile_model(model, "bare").metadata["fermi_energy_subtracted_eV"] == 0.1
    base = dict(word="Tr[delta(E-H)]", expansion={"moments": 8}, basis={"name": "eigen", "supercell": [2, 2, 1]})
    assert RunFile.from_dict(dict(base, model=model)).model["kind"] == "threebodytb"
    with pytest.raises(ValueError, match="requires 'file'"):
        RunFile.from_dict(dict(base, model={"kind": "threebodytb"}))
    with pytest.raises(ValueError, match="covariant"):
        RunFile.from_dict(dict(base, model=dict(model, params={"velocity_mode": "covariant"})))
    RunFile.from_dict(dict(base, model=dict(model, params={"orbitals_per_element": {"C": 4, "H": 1}})))
    with pytest.raises(ValueError, match="positive integers"):
        RunFile.from_dict(dict(base, model=dict(model, params={"orbitals_per_element": {"C": 0}})))
    out = tmp_path / "dos"
    execute_runfile({
        "model": {"kind": "threebodytb", "file": str(path)},
        "basis": {"name": "eigen", "supercell": [3, 3, 1]},
        "word": "Tr[delta(E-H)]",
        "expansion": {"moments": 32, "window": "eigen", "trace": "exact"},
        "reconstruction": {"kernel": "jackson", "formulas": ["dos"], "broadening_meV": 500.0,
                           "energy_grid": [-5.0, 5.0, 21]},
        "output": {"dir": str(out)},
    })
    assert json.loads((out / "run.resolved.json").read_text())["model"]["kind"] == "threebodytb"
    assert json.loads((out / "basis.json").read_text())["source"] == "threebodytb"
    with pytest.raises(ValueError, match="spinless"):
        execute_runfile({
            "model": {"kind": "threebodytb", "file": str(path)},
            "basis": {"name": "eigen", "supercell": [2, 2, 1]},
            "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
            "expansion": {"moments": 8, "window": "eigen", "trace": "exact"},
            "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 500.0,
                               "energy_grid": [-1.0, 1.0, 5]},
            "output": {"dir": str(tmp_path / "shc")},
        })
