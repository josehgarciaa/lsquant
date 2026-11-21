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
    configured root; a model directory outside the root is refused.
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
