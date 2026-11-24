"""Calculation identity gates for run-file checkpoints; requires the compiled core."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import sys

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))
pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant import _core  # noqa: E402
from lsquant.runfile import execute_runfile  # noqa: E402


@pytest.fixture
def calculation(tmp_path):
    # Off-diagonal entries make random-phase moments depend on the state seed.
    matrix = tmp_path / "toy.HAM.CSR"
    matrix.write_text("3 7\n0.2 0  0.4 0  0.4 0  -0.3 0  0.5 0  0.5 0  0.1 0\n"
                      "0 1 0 1 2 1 2\n0 2 5 7\n")
    raw = {
        "model": {"kind": "csr", "dir": str(tmp_path), "label": "toy", "operators": []},
        "basis": {"name": "real"},
        "word": "Tr[delta(E-H)]",
        "expansion": {"moments": 8, "window": [-2.0, 2.0], "trace": "random",
                      "states": 4, "seed": 17, "checkpoint": True, "state_block": 1},
        "output": {"legacy_layout": True, "write_chebmom": False},
    }
    yield raw, matrix
    _core.ConfigureEngineDiagnostics(False, False, -1)


def run(raw, out):
    return execute_runfile(raw, out_dir_override=str(out), repo_root=str(REPO))


def interrupt(raw, out):
    # The public diagnostics hook is the target of LSQ_CHECKPOINT_ABORT_AFTER.
    _core.ConfigureEngineDiagnostics(False, False, 1)
    try:
        with pytest.raises(RuntimeError, match="LSQ_CHECKPOINT_ABORT_AFTER"):
            run(raw, out)
    finally:
        _core.ConfigureEngineDiagnostics(False, False, -1)
    paths = list(out.glob("checkpoint_*.lsqckpt"))
    assert len(paths) == 1
    return paths[0]


def moments(result):
    return [Path(p).read_bytes() for p in result["moment_files"]]


def test_same_calculation_resumes_byte_identically(calculation, tmp_path):
    raw, _ = calculation
    out = tmp_path / "resume"
    checkpoint = interrupt(raw, out)
    first = checkpoint.read_bytes()
    assert interrupt(raw, out) == checkpoint
    # A second interruption advances to the next state, proving an actual resume.
    assert checkpoint.read_bytes() != first
    assert moments(run(raw, out)) == moments(run(raw, tmp_path / "fresh"))
    assert not checkpoint.exists()
    # Pins the new namespace as well as retaining the existing resume guarantee.
    assert checkpoint.name.startswith("checkpoint_v2_")


@pytest.mark.parametrize("change", ["window", "seed", "hamiltonian"])
def test_changed_calculation_ignores_checkpoint(calculation, tmp_path, change):
    raw, matrix = calculation
    out = tmp_path / "resume"
    checkpoint = interrupt(raw, out)
    original = checkpoint.read_bytes()
    if change == "window":
        raw["expansion"]["window"] = [-3.0, 3.0]
    elif change == "seed":
        raw["expansion"]["seed"] = 23
    else:
        stat = matrix.stat()
        matrix.write_text(matrix.read_text().replace("0.2 0", "0.7 0", 1))
        os.utime(matrix, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        assert matrix.stat().st_size == stat.st_size
    with pytest.warns(RuntimeWarning, match="Ignoring checkpoint.*fingerprint"):
        resumed = run(raw, out)
    assert checkpoint.read_bytes() == original
    assert moments(resumed) == moments(run(raw, tmp_path / "fresh"))


def test_legacy_checkpoint_is_ignored(calculation, tmp_path):
    raw, _ = calculation
    out = tmp_path / "resume"
    checkpoint = interrupt(raw, out)
    tag = hashlib.sha256(raw["word"].encode()).hexdigest()[:12]
    legacy = out / f"checkpoint_{tag}.lsqckpt"
    checkpoint.rename(legacy)
    with pytest.warns(RuntimeWarning, match="Ignoring checkpoint.*legacy"):
        resumed = run(raw, out)
    assert legacy.exists()
    assert moments(resumed) == moments(run(raw, tmp_path / "fresh"))
