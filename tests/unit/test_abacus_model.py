"""ABACUS-format Hamiltonians (NextHAM ``hr_tot``) in the metric executors.

The fixtures are written here: an ABACUS ``STRU`` (two species, one with
two atoms), numerical orbital file headers, and sparse H(R) and S(R) files
in the text and the binary layouts that NextHAM's ``add_element.py`` reads
and writes. The spinful H is built as NextHAM builds it: per atom pair a
``[[uu, ud], [du, dd]]`` block, transposed to orbital-major with spin
fastest. The checks:
- H(k) and S(k) of the LSQuant model equal pyatb's sum
  ``sum_R exp(+2 pi i k.R) X(R)``, in eV, with the Fermi energy removed;
- the binary file reads like the text file; the orbital counts come from
  the orbital files, or from ``orbitals_per_element``;
- the spin operator of the metric basis is +1 on the up components of
  NextHAM's blocks;
- a dimension that does not match the basis is refused;
- the run file routes ``model.kind = "abacus"``, resolves the param paths
  and runs DOS and spin Hall.
"""
import json
import os
import struct
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant.bases import BlochMetricBasis  # noqa: E402
from lsquant.models.abacus import ANGSTROM_BOHR, RYDBERG_EV, load_abacus_model, read_abacus_csr  # noqa: E402
from lsquant.models.lcao import cell_gauge_hk  # noqa: E402
from lsquant.models.orthogonal import load_runfile_model  # noqa: E402
from lsquant.runfile import RunFile, execute_runfile, uses_metric_executor  # noqa: E402

ALAT = 7.5                                                   # Bohr
VECTORS = np.array([[1.0, 0.0, 0.0], [0.5, 0.8, 0.0], [0.0, 0.0, 4.0]])
DIRECT = np.array([[0.0, 0.0, 0.5], [0.5, 0.0, 0.5], [0.3, 0.4, 0.55]])   # Bi, Bi, Se
SHELLS = {"Bi": [1, 1], "Se": [1]}                          # Bi: 1s 1p (4 orbitals), Se: 1s
NORB = np.array([4, 4, 1])
CELLS = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
K = np.array([[0.13, 0.27, 0.0], [-0.31, 0.08, 0.0], [0.45, -0.22, 0.0]])


def _write_basis(path):
    (path / "STRU").write_text(f"""ATOMIC_SPECIES
Bi 208.98 Bi.upf
Se 78.97 Se.upf   // comment

NUMERICAL_ORBITAL
Bi_gga_7au_100Ry_1s1p.orb
Se_gga_8au_100Ry_1s.orb

LATTICE_CONSTANT
{ALAT}

LATTICE_VECTORS
{' '.join(map(str, VECTORS[0]))}
{' '.join(map(str, VECTORS[1]))}
{' '.join(map(str, VECTORS[2]))}

ATOMIC_POSITIONS
Direct
Bi
0.0
2
{' '.join(map(str, DIRECT[0]))} 1 1 1
{' '.join(map(str, DIRECT[1]))} m 1 1 1
Se
0.0
1
{' '.join(map(str, DIRECT[2]))} 1 1 1
""")
    for name, shells in (("Bi_gga_7au_100Ry_1s1p.orb", SHELLS["Bi"]), ("Se_gga_8au_100Ry_1s.orb", SHELLS["Se"])):
        labels = "SPDF"
        header = ["-" * 40, "Element                     X", "Energy Cutoff(Ry)          100",
                  "Radius Cutoff(a.u.)         7", f"Lmax                        {len(shells) - 1}"]
        header += [f"Number of {labels[l]}orbital-->       {n}" for l, n in enumerate(shells)]
        (path / name).write_text("\n".join(header) + "\n" + "-" * 40 + "\nMesh 5\ndr 0.01\n")


def _hermitian_set(rng, n, scale, complex_, identity=False):
    mats = {}
    for R in CELLS:
        neg = tuple(-x for x in R)
        if neg in mats:
            mats[R] = mats[neg].conj().T
            continue
        m = scale * (rng.normal(size=(n, n)) + (1j * rng.normal(size=(n, n)) if complex_ else 0))
        if R == neg:
            m = 0.5 * (m + m.conj().T) + (np.eye(n) if identity else 0)
        mats[R] = m
    return np.array([mats[R] for R in CELLS])


def _nextham_interleave(spin_major):
    """NextHAM's ``reshape((2, n, 2, n)).transpose((1, 0, 3, 2))`` on a full spin-major matrix."""
    n = spin_major.shape[-1] // 2
    return spin_major.reshape(2, n, 2, n).transpose(1, 0, 3, 2).reshape(2 * n, 2 * n)


def _write_csr_text(path, mats):
    complex_ = np.iscomplexobj(mats)
    out = ["STEP: 0", f"Matrix Dimension of H(R): {mats.shape[1]}", f"Matrix number of H(R): {len(mats)}"]
    for R, m in zip(CELLS, mats):
        rows, cols = np.nonzero(m)
        out.append(f"{R[0]} {R[1]} {R[2]} {len(rows)}")
        if len(rows):
            data = m[rows, cols]
            out.append(" ".join(f"({v.real:.15e},{v.imag:.15e})" for v in data) if complex_
                       else " ".join(repr(float(v.real)) for v in data))
            out.append(" ".join(map(str, cols)))
            out.append(" ".join(map(str, np.concatenate([[0], np.cumsum(np.bincount(rows, minlength=m.shape[0]))]))))
    path.write_text("\n".join(out) + "\n")


def _write_csr_binary(path, mats):
    n = mats.shape[1]
    complex_ = np.iscomplexobj(mats)
    buf = [struct.pack("iii", 0, n, len(mats))]
    for R, m in zip(CELLS, mats):
        rows, cols = np.nonzero(m)
        buf.append(struct.pack("iiii", *R, len(rows)))
        if len(rows):
            data = m[rows, cols]
            buf.append(np.asarray(data, complex if complex_ else float).tobytes())
            buf.append(cols.astype(np.int32).tobytes())
        buf.append(np.concatenate([[0], np.cumsum(np.bincount(rows, minlength=n))]).astype(np.int32).tobytes())
    path.write_bytes(b"".join(buf))


def _fixture(path, spinful, seed=0, binary=False):
    rng = np.random.default_rng(seed)
    n = NORB.sum()
    _write_basis(path)
    s = _hermitian_set(rng, n, 0.03, False, identity=True).real
    if spinful:
        h_spin_major = _hermitian_set(rng, 2 * n, 0.01, True)
        h = np.array([_nextham_interleave(m) for m in h_spin_major])
        s_file = np.array([np.kron(m, np.eye(2)) for m in s]).astype(complex)
    else:
        h_spin_major = None
        h = _hermitian_set(rng, n, 0.01, False).real
        s_file = s
    write = _write_csr_binary if binary else _write_csr_text
    write(path / "hr_tot", h.astype(complex) if spinful else h)
    write(path / "srs1_nao.csr", s_file)
    return h, s_file, h_spin_major


def _pyatb_sum(mats, k):
    return np.einsum("r,rmn->mn", np.exp(2j * np.pi * np.array(CELLS) @ k), mats)


@pytest.mark.parametrize("spinful,binary", [(False, False), (True, False), (True, True)])
def test_matches_pyatb_sum(tmp_path, spinful, binary):
    h, s, _ = _fixture(tmp_path, spinful, binary=binary)
    model = load_abacus_model(tmp_path / "hr_tot", overlap_file=tmp_path / "srs1_nao.csr",
                              stru_file=tmp_path / "STRU", fermi_energy=-2.0)
    assert model.spinful is spinful and model.metadata["orthogonal"] is False
    lattice = VECTORS * ALAT / ANGSTROM_BOHR
    np.testing.assert_allclose(model.uc, lattice, atol=1e-12)
    np.testing.assert_allclose(model.tau[0], DIRECT[0] @ lattice, atol=1e-12)
    assert model.metadata["orbital_l"] == [0, 1, 1, 1] * 2 + [0]
    hk = cell_gauge_hk(model, K)
    for n, k in enumerate(K):
        s_ref = _pyatb_sum(s, k)
        np.testing.assert_allclose(hk[n], _pyatb_sum(h, k) * RYDBERG_EV + 2.0 * s_ref, atol=1e-10)


def test_text_and_binary_files_agree(tmp_path):
    text_dir, bin_dir = tmp_path / "t", tmp_path / "b"
    text_dir.mkdir(), bin_dir.mkdir()
    _fixture(text_dir, True, seed=3)
    _fixture(bin_dir, True, seed=3, binary=True)
    for name in ("hr_tot", "srs1_nao.csr"):
        rt, mt = read_abacus_csr(text_dir / name)
        rb, mb = read_abacus_csr(bin_dir / name)
        np.testing.assert_array_equal(rt, rb)
        np.testing.assert_allclose(mt, mb, atol=1e-14)


def test_spin_operator_follows_nextham_blocks(tmp_path):
    _fixture(tmp_path, True, seed=1)
    model = load_abacus_model(tmp_path / "hr_tot", overlap_file=tmp_path / "srs1_nao.csr", stru_file=tmp_path / "STRU")
    # On-site H = sigma_z on every orbital, built in NextHAM's spin-major block and interleaved.
    n = NORB.sum()
    sz = _nextham_interleave(np.kron(np.diag([1.0, -1.0]), np.eye(n)))
    np.testing.assert_allclose(sz, np.kron(np.eye(n), np.diag([1.0, -1.0])))
    basis = BlochMetricBasis(model, [2, 2, 1], velocity_mode="berry_connection")
    np.testing.assert_allclose(basis.operator_blocks("SZ"), basis.sk_array @ sz, atol=1e-12)


def test_orbital_counts_and_refusals(tmp_path):
    _fixture(tmp_path, False, seed=2)
    args = dict(overlap_file=tmp_path / "srs1_nao.csr", stru_file=tmp_path / "STRU")
    moved = tmp_path / "orbitals"
    moved.mkdir()
    for f in tmp_path.glob("*.orb"):
        f.rename(moved / f.name)
    with pytest.raises(ValueError, match="orbital_dir"):
        load_abacus_model(tmp_path / "hr_tot", **args)
    a = load_abacus_model(tmp_path / "hr_tot", orbital_dir=moved, **args)
    b = load_abacus_model(tmp_path / "hr_tot", orbitals_per_element={"Bi": 4, "Se": 1}, **args)
    np.testing.assert_allclose(cell_gauge_hk(a, K), cell_gauge_hk(b, K), atol=1e-12)
    assert set(b.metadata["orbital_l"]) == {-1}
    with pytest.raises(ValueError, match="dimension"):
        load_abacus_model(tmp_path / "hr_tot", orbitals_per_element={"Bi": 9, "Se": 1}, **args)


def test_runfile_routes_and_runs_abacus(tmp_path):
    _fixture(tmp_path, True, seed=4)
    model = {"kind": "abacus", "file": "hr_tot",
             "params": {"overlap_file": "srs1_nao.csr", "stru_file": "STRU", "fermi_energy": 0.1}}
    base = dict(word="Tr[delta(E-H)]", expansion={"moments": 8}, basis={"name": "eigen", "supercell": [2, 2, 1]})
    runfile = tmp_path / "run.json"
    runfile.write_text(json.dumps(dict(base, model=model)))
    rf = RunFile.from_file(str(runfile))
    assert rf.model["params"]["stru_file"] == str(tmp_path / "STRU")
    assert uses_metric_executor(rf.model)
    assert load_runfile_model(rf.model, "bare").metadata["fermi_energy_subtracted_eV"] == 0.1
    with pytest.raises(ValueError, match="stru_file"):
        RunFile.from_dict(dict(base, model={"kind": "abacus", "file": str(tmp_path / "hr_tot"),
                                            "params": {"overlap_file": "x"}}))
    with pytest.raises(ValueError, match="covariant"):
        RunFile.from_dict(dict(base, model=dict(model, params=dict(model["params"], velocity_mode="covariant"))))
    full = {"kind": "abacus", "file": str(tmp_path / "hr_tot"),
            "params": {"overlap_file": str(tmp_path / "srs1_nao.csr"), "stru_file": str(tmp_path / "STRU")}}
    for n, (word, formula) in enumerate((("Tr[delta(E-H)]", "dos"), ("Tr[JXSZ delta(E-H) VY dG-(E)/dE]", "shc"))):
        out = tmp_path / f"{formula}{n}"
        execute_runfile({
            "model": full,
            "basis": {"name": "eigen", "supercell": [3, 3, 1]},
            "word": word,
            "expansion": {"moments": 32, "window": "eigen", "trace": "exact"},
            "reconstruction": {"kernel": "jackson", "formulas": [formula], "broadening_meV": 500.0,
                               "energy_grid": [-2.0, 2.0, 21]},
            "output": {"dir": str(out)},
        })
        assert json.loads((out / "basis.json").read_text())["source"] == "abacus"
