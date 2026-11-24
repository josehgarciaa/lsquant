"""Unit tests of stage 7: observable words, run files, the CLI and the MCP tools.

They need the compiled module ``lsquant._core`` (see test_python_stage2.py); without it
every test is skipped. The MCP tests call the tool functions directly and need no ``mcp``
package.

What is checked:
  * words: parsing, the catalogue, compilation to a plan with one Expand per spectral
    slot, refusal of unsupported words, ``describe``;
  * run files: g1.json and g5.json of runs/goldens reproduce the goldens G1 and the four G5
    curves byte for byte; strict validation refuses unknown keys, a missing word and
    ``window = "eigen"`` on a csr model;
  * the CLI: ``inspect`` prints JSON and exits 0, a refused run file exits non-zero and
    writes no output directory;
  * the MCP tools: ``words_describe``; ``run`` and ``status`` on g1.json inside the
    configured root; explicit input/output and tool-target escapes are refused,
    including connector paths and symbolic links, without disclosing private targets.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time

import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GOLD = os.path.join(REPO, "goldens", "stage0")
RUNS = os.path.join(REPO, "runs", "goldens")
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

lsq = pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant import words  # noqa: E402
from lsquant.runfile import RunFile, execute_runfile, validate_runfile_dict  # noqa: E402


def load(name):
    with open(os.path.join(RUNS, name)) as fh:
        return json.load(fh)


# ------------------------------------------------------------------ words
def test_parse_and_compile_two_slot_word():
    w = words.parse("Tr[VX delta(E-H) VY delta(E-H)]")
    assert w is not None
    cw = words.compile("Tr[VX delta(E-H) VY delta(E-H)]", basis="real",
                       params={"moments": 16, "states": 4, "window": [-2.0, 2.0]})
    kinds = [type(i).__name__ for i in cw.plan.code]
    assert kinds.count("Expand") == 2


def test_catalogue_lists_entries():
    cat = words.list_catalogue()
    assert len(cat) == 7
    assert {"dos_identity", "kubo_greenwood", "kubo_bastin", "msd"} <= {e["shape_id"] for e in cat}
    assert all("canonical_example" in e and "plan_builder" in e for e in cat)


def test_describe_returns_mapping():
    info = words.describe("Tr[delta(E-H)]")
    assert isinstance(info, dict) and info


@pytest.mark.parametrize("bad", ["Tr[VX G+(E) VY delta(E-H)]", "Tr[VX delta(E-H) VY", "VX delta(E-H)"])
def test_unsupported_or_malformed_words_are_refused(bad):
    with pytest.raises(Exception):
        words.compile(bad, basis="real", params={"moments": 8, "states": 1, "window": [-1.0, 1.0]})


# ------------------------------------------------------------------ run files
def test_g1_runfile_reproduces_golden(tmp_path):
    res = execute_runfile(load("g1.json"), out_dir_override=str(tmp_path), repo_root=REPO)
    name = "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"
    assert any(p.endswith(name) for p in res["written_files"])
    assert (tmp_path / name).read_bytes() == open(os.path.join(GOLD, "chain1d", name), "rb").read()


def test_g5_runfile_reproduces_golden_curves(tmp_path):
    execute_runfile(load("g5.json"), out_dir_override=str(tmp_path), repo_root=REPO)
    for g in ("KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
              "KuboBastinSea_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
              "KuboBastinSurf_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat",
              "KuboGreenWood_NonEqOpVX-VXhaldaneKPM_M128x128_stateexactJACKSON.dat"):
        assert (tmp_path / g).read_bytes() == open(os.path.join(GOLD, "haldane", g), "rb").read(), g


def test_validation_refuses_unknown_key():
    raw = load("g1.json")
    raw["expansion"]["momentz"] = 64
    with pytest.raises(ValueError):
        validate_runfile_dict(raw)


def test_validation_refuses_missing_word():
    raw = load("g1.json")
    del raw["word"]
    with pytest.raises(ValueError):
        RunFile.from_dict(raw, repo_root=REPO)


def test_eigen_window_refused_on_csr_model():
    raw = load("g1.json")
    raw["expansion"]["window"] = "eigen"
    with pytest.raises(ValueError):
        RunFile.from_dict(raw, repo_root=REPO)


# ------------------------------------------------------------------ CLI
def run_cli(*args, cwd=REPO):
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(os.path.abspath(p) for p in sys.path if p)
    return subprocess.run([sys.executable, "-m", "lsquant", *args], cwd=cwd, env=env,
                          capture_output=True, text=True, timeout=600)


def test_cli_inspect_prints_json():
    r = run_cli("inspect", os.path.join(RUNS, "g1.json"))
    assert r.returncode == 0, r.stderr
    assert isinstance(json.loads(r.stdout), dict)


def test_cli_refused_runfile_exits_nonzero(tmp_path):
    raw = load("g1.json")
    raw["unknown_section"] = {}
    rf = tmp_path / "bad.json"
    rf.write_text(json.dumps(raw))
    out = tmp_path / "out"
    r = run_cli("run", str(rf), "--out", str(out))
    assert r.returncode != 0
    assert r.stderr.strip()
    assert not out.exists()


# ------------------------------------------------------------------ MCP tools
def test_mcp_words_describe():
    from lsquant import mcp_server
    assert json.loads(mcp_server.tool_words_describe("Tr[delta(E-H)]"))


MCP_PATH_FIELDS = (
    "model.dir", "model.csr_dir", "model.file", "model.hr_file",
    "model.positions", "model.atom_map", "model.params.overlap_file",
    "model.params.stru_file", "model.params.orbital_dir",
    "model.params.position_file", "model.params.dipole_cache",
    "model.params.hamiltonian_file", "output.dir",
)


def _skip_unless_field_is_available(field):
    """``model.hr_file`` needs ``numerics.operator = "atom_block"``, a native-profile operator;
    the golden profile refuses the run file before the path guard is reached."""
    from lsquant import _core
    if field == "model.hr_file" and _core.build_profile != "native":
        pytest.skip("model.hr_file needs the native build profile (atom_block operator)")


def _mcp_path_run(field, value):
    raw = load("g1.json")
    raw["model"]["dir"] = "."
    if field.startswith("model.params."):
        kind = ("gpaw" if field.endswith("dipole_cache") else
                "deeph" if field.endswith("hamiltonian_file") else "abacus")
        raw["model"] = {"kind": kind, "dir": ".", "file": "input",
                        "params": {"velocity_mode": "bare"}}
        if kind == "abacus":
            raw["model"]["params"].update(overlap_file="input", stru_file="input")
    if field == "model.hr_file":
        raw["numerics"] = {"operator": "atom_block"}
    target = raw
    parts = field.split(".")
    for part in parts[:-1]:
        target = target[part]
    target[parts[-1]] = value
    return raw


@pytest.mark.parametrize("field", MCP_PATH_FIELDS)
@pytest.mark.parametrize("escape", ["absolute", "traversal", "symlink"])
def test_mcp_run_refuses_all_path_escapes(tmp_path, monkeypatch, field, escape):
    from lsquant import mcp_server
    _skip_unless_field_is_available(field)
    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "private"
    outside.mkdir()
    (root / "link").symlink_to(outside, target_is_directory=True)
    value = {"absolute": str(outside), "traversal": "../private",
             "symlink": "link"}[escape]
    monkeypatch.setattr(mcp_server, "_CONFIGURED_ROOT", str(root))
    calls = []
    monkeypatch.setattr(mcp_server.subprocess, "Popen", lambda *a, **k: calls.append(a))
    with pytest.raises(ValueError) as error:
        mcp_server.tool_run(_mcp_path_run(field, value))
    message = str(error.value)
    assert field in message
    assert "outside the server root" in message
    if escape != "absolute":
        assert str(outside) not in message
    assert not calls
    assert not list(root.glob(".mcp_run_*"))
    assert not (root / "out").exists()


@pytest.mark.parametrize("field", MCP_PATH_FIELDS)
def test_mcp_run_accepts_inside_paths(tmp_path, monkeypatch, field):
    from lsquant import mcp_server
    _skip_unless_field_is_available(field)
    (tmp_path / "input").touch()
    monkeypatch.setattr(mcp_server, "_CONFIGURED_ROOT", str(tmp_path))
    monkeypatch.setattr(mcp_server, "_ACTIVE_RUNS", {})
    calls = []
    monkeypatch.setattr(mcp_server.subprocess, "Popen", lambda *a, **k: calls.append(k))
    raw = _mcp_path_run(field, "input")
    if field == "output.dir":
        raw["output"]["dir"] = "."
    result = json.loads(mcp_server.tool_run(raw))
    assert len(calls) == 1
    assert calls[0]["cwd"] == str(tmp_path)
    if field == "output.dir":
        assert result["output_dir"] == str(tmp_path)
    mcp_server._ACTIVE_RUNS[result["run_id"]].log_handle.close()


@pytest.mark.parametrize("tool,key", [("tool_inspect", "path"),
                                      ("tool_reconstruct", "output_dir"),
                                      ("tool_run", "runfile_json")])
@pytest.mark.parametrize("escape", ["absolute", "traversal", "symlink"])
def test_mcp_tool_targets_refuse_escapes(tmp_path, monkeypatch, tool, key, escape):
    from lsquant import mcp_server
    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "private"
    outside.touch()
    (root / "link").symlink_to(outside)
    value = {"absolute": str(outside), "traversal": "../private",
             "symlink": "link"}[escape]
    monkeypatch.setattr(mcp_server, "_CONFIGURED_ROOT", str(root))
    with pytest.raises(ValueError) as error:
        getattr(mcp_server, tool)(value)
    assert key in str(error.value)
    assert "outside the server root" in str(error.value)
    if escape != "absolute":
        assert str(outside) not in str(error.value)


@pytest.fixture
def mcp_no_dispatch(monkeypatch):
    from lsquant import mcp_server

    def refuse_dispatch(*args, **kwargs):
        pytest.fail("path validation must refuse the run before dispatch")

    monkeypatch.setattr(mcp_server.subprocess, "Popen", refuse_dispatch)


def test_mcp_run_checks_source_relative_path(tmp_path, monkeypatch, mcp_no_dispatch):
    from lsquant import mcp_server
    root = tmp_path / "root"
    nested = root / "configs"
    nested.mkdir(parents=True)
    outside = tmp_path / "private"
    outside.touch()
    (nested / "input").symlink_to(outside)
    source = nested / "run.json"
    source.write_text(json.dumps(_mcp_path_run("model.positions", "input")))
    monkeypatch.setattr(mcp_server, "_CONFIGURED_ROOT", str(root))
    with pytest.raises(ValueError, match="model.positions.*outside the server root") as error:
        mcp_server.tool_run("configs/run.json")
    assert str(outside) not in str(error.value)


@pytest.mark.parametrize("at_root", [False, True])
def test_mcp_run_preserves_input_resolution(tmp_path, monkeypatch, at_root):
    from lsquant import mcp_server
    nested = tmp_path / "configs"
    nested.mkdir()
    (nested / "input").touch()
    if at_root:
        (tmp_path / "input").touch()
    source = nested / "run.json"
    source.write_text(json.dumps(_mcp_path_run("model.positions", "input")))
    monkeypatch.setattr(mcp_server, "_CONFIGURED_ROOT", str(tmp_path))
    monkeypatch.setattr(mcp_server, "_ACTIVE_RUNS", {})
    calls = []
    monkeypatch.setattr(mcp_server.subprocess, "Popen", lambda *a, **k: calls.append(a))
    result = json.loads(mcp_server.tool_run("configs/run.json"))
    with open(calls[0][0][4], encoding="utf-8") as fh:
        payload = json.load(fh)
    expected = (tmp_path if at_root else nested) / "input"
    assert payload["model"]["positions"] == str(expected)
    mcp_server._ACTIVE_RUNS[result["run_id"]].log_handle.close()


def test_mcp_run_checks_deeph_path_relative_to_model_dir(tmp_path, monkeypatch, mcp_no_dispatch):
    from lsquant import mcp_server
    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "private"
    outside.mkdir()
    model_dir = root / "model"
    model_dir.mkdir()
    (model_dir / "link").symlink_to(outside, target_is_directory=True)
    # Only the path relative to model.dir encounters the escaping symlink.
    raw = _mcp_path_run("model.params.hamiltonian_file", "link/input.h5")
    raw["model"]["dir"] = str(model_dir)
    monkeypatch.setattr(mcp_server, "_CONFIGURED_ROOT", str(root))
    with pytest.raises(ValueError, match="model.params.hamiltonian_file.*outside the server root"):
        mcp_server.tool_run(raw)


def test_mcp_run_checks_gpaw_cache_write_suffix(tmp_path, monkeypatch, mcp_no_dispatch):
    from lsquant import mcp_server
    root = tmp_path / "root"
    root.mkdir()
    outside = tmp_path / "private.npz"
    outside.touch()
    (root / "cache.npz").symlink_to(outside)
    monkeypatch.setattr(mcp_server, "_CONFIGURED_ROOT", str(root))
    with pytest.raises(ValueError, match="model.params.dipole_cache.*outside the server root") as error:
        mcp_server.tool_run(_mcp_path_run("model.params.dipole_cache", "cache"))
    assert str(outside) not in str(error.value)


def test_mcp_run_checks_generator_symlink(tmp_path, monkeypatch, mcp_no_dispatch):
    from lsquant import mcp_server
    from lsquant.runfile import GENERATOR_CATALOGUE
    name, relative = next(iter(GENERATOR_CATALOGUE.items()))
    root = tmp_path / "root"
    script = root / relative
    script.parent.mkdir(parents=True)
    outside = tmp_path / "private.py"
    outside.touch()
    script.symlink_to(outside)
    raw = load("g1.json")
    raw["model"] = {"kind": "generator", "name": name}
    monkeypatch.setattr(mcp_server, "_CONFIGURED_ROOT", str(root))
    with pytest.raises(ValueError, match="model.name.*outside the server root") as error:
        mcp_server.tool_run(raw)
    assert str(outside) not in str(error.value)


def test_mcp_default_output_and_inline_model_data(tmp_path, monkeypatch):
    from lsquant import mcp_server
    raw = _mcp_path_run("model.dir", ".")
    raw.pop("output")
    raw["model"].update(positions=[[0.0, 0.0, 0.0]], atom_map=[0])
    monkeypatch.setattr(mcp_server, "_CONFIGURED_ROOT", str(tmp_path))
    monkeypatch.setattr(mcp_server, "_ACTIVE_RUNS", {})
    monkeypatch.setattr(mcp_server.subprocess, "Popen", lambda *a, **k: None)
    result = json.loads(mcp_server.tool_run(raw))
    assert result["output_dir"] == str(tmp_path)
    mcp_server._ACTIVE_RUNS[result["run_id"]].log_handle.close()


def test_mcp_run_status_and_root(tmp_path, monkeypatch):
    import shutil
    from lsquant import mcp_server
    # the run is a subprocess started in the root: give it absolute import paths
    monkeypatch.setenv("PYTHONPATH", os.pathsep.join(os.path.abspath(p) for p in sys.path if p))
    shutil.copytree(os.path.join(GOLD, "chain1d", "operators"), tmp_path / "ops")
    mcp_server.set_configured_root(str(tmp_path))
    raw = load("g1.json")
    raw["model"]["dir"] = "ops"
    raw["output"]["dir"] = "out"
    started = json.loads(mcp_server.tool_run(raw))
    for _ in range(1200):
        st = json.loads(mcp_server.tool_status(started["run_id"]))
        if st["status"] != "running":
            break
        time.sleep(0.5)
    assert st["status"] == "succeeded", st
    name = "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"
    assert (tmp_path / "out" / name).read_bytes() == open(os.path.join(GOLD, "chain1d", name), "rb").read()
    outside = load("g1.json")
    outside["model"]["dir"] = os.path.join(GOLD, "chain1d", "operators")     # outside the root
    with pytest.raises(Exception):
        mcp_server.tool_run(outside)
