"""Small regressions for SHC run metadata, directory filtering and curve labels."""

import json
from pathlib import Path
import sys

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))
lsq = pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant import _core
from lsquant.cli import main
from lsquant.runfile import execute_runfile, reconstruct_target, validate_runfile_dict


def config(kind="csr"):
    model = {"kind": "csr", "dir": "goldens/stage0/chain1d/operators", "label": "chain1d"}
    basis = {"name": "real"}
    if kind == "wannier":
        model = {"kind": "wannier", "dir": "examples/14_wannier/ssh", "seed": "ssh"}
        basis = {"name": "bloch", "supercell": [2, 1, 1]}
    return {
        "model": model, "basis": basis,
        "word": "Tr[VX delta(E-H) VX dG-(E)/dE]",
        "expansion": {"moments": 8, "window": [-4.0, 4.0], "trace": "exact", "seed": 12345},
        "reconstruction": {"formulas": ["shc"], "broadening_meV": 500.0,
                           "num_cells": 2, "cell_area": 3.0},
        "output": {"write_chebmom": True, "legacy_layout": False},
    }


def small_csr(raw, tmp_path):
    """Use a two-site Hermitian model so exact traces stay inexpensive."""
    ops = tmp_path / "ops"
    ops.mkdir()
    (ops / "toy.HAM.CSR").write_text("2 2\n-1 0 -1 0\n1 0\n0 1 2\n")
    (ops / "toy.VX.CSR").write_text("2 2\n0 1 0 -1\n1 0\n0 1 2\n")
    raw["model"] = {"kind": "csr", "dir": str(ops), "label": "toy"}


@pytest.mark.parametrize("kind", ["csr", "wannier"])
@pytest.mark.parametrize("missing", [("num_cells",), ("cell_area",), ("num_cells", "cell_area")])
@pytest.mark.parametrize("target", [False, True])
def test_shc_missing_geometry_refused_at_validation(kind, missing, target):
    raw = config(kind)
    if target:
        raw["reconstruction"]["formulas"] = ["bastin"]
        raw["reconstruction"]["targets"] = [{"formulas": ["shc"], "broadening_meV": 500.0}]
    for key in missing:
        del raw["reconstruction"][key]
    with pytest.raises(ValueError, match="num_cells.*cell_area"):
        validate_runfile_dict(raw)


def test_shc_missing_geometry_creates_no_output(tmp_path):
    raw = config()
    small_csr(raw, tmp_path)
    del raw["model"]["kind"]
    del raw["reconstruction"]["num_cells"]
    out = tmp_path / "refused"
    with pytest.raises(ValueError, match="num_cells.*cell_area"):
        execute_runfile(raw, out_dir_override=str(out), repo_root=str(REPO))
    assert not out.exists()


@pytest.mark.parametrize("kind", ["csr", "wannier"])
def test_shc_run_matches_explicit_reconstruction(tmp_path, kind):
    raw = config(kind)
    if kind == "csr":
        small_csr(raw, tmp_path)
    result = execute_runfile(raw, out_dir_override=str(tmp_path / "run"), repo_root=str(REPO))
    explicit = reconstruct_target(result["moment_files"][0], formula="shc", broadening_meV=500.0,
                                  out_dir=str(tmp_path / "explicit"))
    assert len(result["curve_files"]) == len(explicit["curve_files"]) == 1
    assert Path(result["curve_files"][0]).read_bytes() == Path(explicit["curve_files"][0]).read_bytes()


def test_shc_reconstruction_target_matches_explicit(tmp_path):
    raw = config()
    small_csr(raw, tmp_path)
    raw["reconstruction"]["formulas"] = ["bastin"]
    raw["reconstruction"]["targets"] = [{"formulas": ["shc"], "broadening_meV": 500.0,
                                           "subdir": "target", "also_top_dir": True}]
    result = execute_runfile(raw, out_dir_override=str(tmp_path / "run"), repo_root=str(REPO))
    explicit = reconstruct_target(result["moment_files"][0], formula="shc", broadening_meV=500.0,
                                  out_dir=str(tmp_path / "explicit"))
    target = next(Path(p) for p in result["curve_files"] if Path(p).parent.name == "target")
    assert target.read_bytes() == Path(explicit["curve_files"][0]).read_bytes()
    assert (tmp_path / "run" / target.name).read_bytes() == target.read_bytes()


def moments(tmp_path, name, word, dims=(8, 8), legacy=False):
    mu = _core.ChebyshevMoments(list(dims))
    mu.word = word
    mu.system_size = 2
    mu.window = _core.SpectralWindow(-2.0, 2.0)
    mu.data[:] = 0.0
    mu.data.flat[0] = 1.0
    path = tmp_path / name
    if legacy:
        _core.WriteLegacyChebmom1D(str(path), mu)
    else:
        _core.WriteChebmom(str(path), mu)
    return path


def test_directory_skips_incompatible_rank_and_shape(tmp_path, capsys):
    moments(tmp_path, "a.chebmom", "VX delta VY delta")
    dos = moments(tmp_path, "b.chebmom", "delta", dims=(8,))
    greenwood = moments(tmp_path, "c.chebmom", "Tr[VX delta(E-H) VX delta(E-H)]")
    time_word = moments(tmp_path, "d.chebmom", "VX delta U(t)^+ VY U(t)")
    (tmp_path / "run.resolved.json").write_text(json.dumps({"reconstruction": {"num_cells": 2, "cell_area": 3.0}}))
    result = reconstruct_target(str(tmp_path), formula="shc", broadening_meV=500.0)
    assert len(result["curve_files"]) == 1
    assert {entry["moment_file"] for entry in result["skipped_files"]} == {
        str(dos), str(greenwood), str(time_word)}
    assert all(entry["reason"] for entry in result["skipped_files"])
    assert main(["reconstruct", str(tmp_path), "--formula", "shc", "--broadening-meV", "500"]) == 0
    output = capsys.readouterr()
    assert "Skipped" in output.err and dos.name in output.err and greenwood.name in output.err


def test_directory_no_compatible_files_fails_before_writing(tmp_path):
    moments(tmp_path, "dos.chebmom1D", "delta", dims=(8,), legacy=True)
    with pytest.raises(ValueError, match="No moment files.*shc"):
        reconstruct_target(str(tmp_path), formula="shc")
    assert not list(tmp_path.glob("*.dat"))


def test_legacy_directory_skips_dos_and_prints_json(tmp_path, capsys):
    source = REPO / "goldens" / "stage0" / "chain1d"
    names = ["SpectralOp1chain1dKPM_M64_stateexact.chebmom1D",
             "NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D"]
    for name in names:
        (tmp_path / name).write_bytes((source / name).read_bytes())
    assert main(["reconstruct", str(tmp_path), "--formula", "bastin", "--broadening-meV", "500", "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert len(result["curve_files"]) == 1
    assert [entry["moment_file"] for entry in result["skipped_files"]] == [str(tmp_path / names[0])]


@pytest.mark.parametrize("word, expected", [
    ("JXSZ delta VY delta", "energy_eV sigma_z_xy [e^2/h, sigma_z Pauli]"),
    ("JYSX delta VZ delta", "energy_eV sigma_x_yz [e^2/h, sigma_x Pauli]"),
    ("VX delta VY delta", "energy_eV normalised_bastin[VX,VY] [2*pi/(num_cells*cell_area)]"),
])
def test_shc_header_follows_stored_word(tmp_path, word, expected):
    path = moments(tmp_path, "tensor.chebmom", word)
    (tmp_path / "run.resolved.json").write_text(json.dumps({"reconstruction": {"num_cells": 2, "cell_area": 3.0}}))
    result = reconstruct_target(str(path), formula="shc", broadening_meV=500.0)
    curve = Path(result["curve_files"][0])
    assert curve.read_text().splitlines()[0] == "# " + expected
    mu = _core.ReadChebmom(str(path))
    energy, values = lsq.shc(mu, 500.0, 2, 3.0)
    np.testing.assert_array_equal(np.loadtxt(curve), np.column_stack((energy, values)))
