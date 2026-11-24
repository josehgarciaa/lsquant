"""Regressions for accepted run files that silently omit requested work.

Tiny diagonal CSR inputs keep execution checks independent of external models.
The generator round trip uses its smallest ring and the fixed default seed.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))
pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant.runfile import RunFile, execute_runfile, validate_runfile_dict


@pytest.fixture
def raw(tmp_path):
    ops = tmp_path / "ops"
    ops.mkdir()
    # CSR text uses complex pairs, column indices, then row pointers (zero based).
    (ops / "tiny.HAM.CSR").write_text("2 2\n-0.5 0 0.5 0\n0 1\n0 1 2\n")
    (ops / "tiny.PTOP.CSR").write_text("2 2\n1 0 0 0\n0 1\n0 1 2\n")
    return dict(model=dict(kind="csr", dir=str(ops), label="tiny"),
                word="Tr[delta(E-H)]",
                expansion=dict(moments=8, window=[-1, 1], trace="exact"))


def test_csr_nonzero_disorder_refused(raw):
    raw["basis"] = dict(disorder=dict(W=0.5))
    with pytest.raises(ValueError, match=r"basis.disorder.W.*csr"):
        validate_runfile_dict(raw)
    raw["basis"]["disorder"]["W"] = 0
    validate_runfile_dict(raw)


def test_missing_explicit_unused_operator_refused(raw, tmp_path):
    raw["model"]["operators"] = ["PTOP", "ABSENT"]
    with pytest.raises(FileNotFoundError, match=r"tiny\.ABSENT\.CSR"):
        execute_runfile(raw, out_dir_override=str(tmp_path / "out"))


@pytest.mark.parametrize("word,formula", [
    (word, formula)
    for word in ("Tr[A delta(E-H) U(t)^+ B U(t)]",
                 "Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]",
                 "Tr[P U(t)^+ SX U(t) delta(E-H) P]")
    for formula in (None, "time_correlation", "time_correlation_curve")
    if not ("Tr[P " in word and formula is None)
])
def test_time_curve_requires_fermi_energies(raw, word, formula):
    raw["word"] = word
    raw["time"] = dict(steps=2, dt=0.1)
    raw["reconstruction"] = dict(broadening_meV=100, fermi_energies=[])
    if formula is not None:
        raw["reconstruction"]["formulas"] = [formula]
    if "Tr[P " in word:
        for row_formula in (None, ["row_m0"], ["row0"]):
            raw["reconstruction"]["formulas"] = row_formula
            validate_runfile_dict(raw)
        raw["reconstruction"]["formulas"] = [formula]
    with pytest.raises(ValueError, match="nonempty reconstruction.fermi_energies"):
        validate_runfile_dict(raw)
    raw["reconstruction"]["fermi_energies"] = [0.0]
    validate_runfile_dict(raw)
    raw["reconstruction"]["broadening_meV"] = None
    raw["reconstruction"]["fermi_energies"] = []
    validate_runfile_dict(raw)


def test_custom_operator_runs(raw, tmp_path):
    raw["model"]["operators"] = ["PTOP"]
    raw["word"] = "Tr[PTOP delta(E-H)]"
    rf = RunFile.from_dict(raw)
    result = execute_runfile(rf, out_dir_override=str(tmp_path / "out"))
    assert result["moment_files"]
    assert all(Path(p).is_file() for p in result["moment_files"])


def test_generator_resolved_file_runs_again(tmp_path):
    raw = dict(model=dict(kind="generator", name="make_chain", tutorial=1,
                          params=dict(sizes=[4], hopping_t=-1.0)),
               word="Tr[delta(E-H)]",
               expansion=dict(moments=8, window=[-2, 2], trace="exact"),
               reconstruction=dict(broadening_meV=100, formulas=["dos"]))
    first = execute_runfile(raw, out_dir_override=str(tmp_path / "first"), repo_root=str(REPO))
    receipt = json.loads(Path(first["resolved_runfile"]).read_text())
    assert "disorder" not in receipt["basis"]
    second = execute_runfile(receipt, out_dir_override=str(tmp_path / "second"), repo_root=str(REPO))
    assert second["moment_files"]
    assert [Path(p).read_bytes() for p in first["moment_files"]] == [
        Path(p).read_bytes() for p in second["moment_files"]]


def test_projection_auto_refused_before_model_load():
    raw = dict(model=dict(kind="siesta", dir="unused", label="tiny"),
               basis=dict(name="bloch"), metric=dict(policy="nonorthogonal"),
               projection=dict(kind="hard", window=[-1, 1]),
               word="Tr[delta(E-H)]", expansion=dict(moments=8, state_block="auto"))
    with pytest.raises(ValueError, match="Projected runs require an integer expansion.state_block"):
        validate_runfile_dict(raw)


def test_auto_block_receipt_matches_execution(raw, tmp_path, monkeypatch):
    from lsquant import runfile
    actual = []
    resolve = runfile._resolve_state_block

    def record(expansion, dim, moments, states):
        block = resolve(expansion, dim, moments, states)
        actual.append((block, expansion["state_block_rule"]))
        return block

    monkeypatch.setattr(runfile, "_resolve_state_block", record)
    raw["expansion"]["state_block"] = "auto"
    result = execute_runfile(raw, out_dir_override=str(tmp_path / "out"))
    receipt = json.loads(Path(result["resolved_runfile"]).read_text())
    assert actual
    assert receipt["expansion"]["state_block"] == actual[-1][0]
    assert receipt["expansion"]["state_block_rule"] == actual[-1][1]


def test_time_curve_target_requires_own_fermi_energies(raw):
    raw["word"] = "Tr[A delta(E-H) U(t)^+ B U(t)]"
    raw["time"] = dict(steps=2, dt=0.1)
    target = dict(broadening_meV=100, formulas=["time_correlation_curve"])
    raw["reconstruction"] = dict(targets=[target])
    with pytest.raises(ValueError, match="nonempty reconstruction.targets.fermi_energies"):
        validate_runfile_dict(raw)
    target["fermi_energies"] = [0.0]
    validate_runfile_dict(raw)
