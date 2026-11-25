"""The reduced-golden tool of the SIESTA run files (scripts/siesta_golden.py).

The tool is the comparison the stage-7 gate uses for the G11 run files, so its
own behaviour is pinned here on synthetic output directories: a reduction
compares equal to itself, a change above the tolerance is found wherever it is
(sampled value, maximum, integral, oracle number, integer), and a missing or
extra quantity is a failure rather than something skipped.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("siesta_golden", REPO / "scripts" / "siesta_golden.py")
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)


def write_output(directory, scale=1.0, bump_row=None, z2=1, chern=2.0):
    directory.mkdir()
    x = np.linspace(-3.0, 3.0, 1000)
    y = scale * np.exp(-x * x)
    if bump_row is not None:
        y[bump_row] += 1e-3
    np.savetxt(directory / "dos.dat", np.column_stack((x, y)), fmt="%.17g", header="energy_eV DOS")
    (directory / "spin_chern.json").write_text(json.dumps(
        {"spin_chern": chern, "z2": z2, "convention": "text is ignored", "list": [0.1, 0.2, 0.3]}))
    return directory


def run(argv):
    return tool.main(["siesta_golden.py", *map(str, argv)])


def test_a_run_matches_its_own_golden(tmp_path, capsys):
    out = write_output(tmp_path / "run")
    assert run(["make", out, tmp_path / "g.json"]) == 0
    assert run(["check", out, tmp_path / "g.json"]) == 0
    text = capsys.readouterr().out
    assert "curves max rel dev 0.000e+00" in text and "oracle max rel dev 0.000e+00" in text
    golden = json.loads((tmp_path / "g.json").read_text())
    assert golden["curves"]["dos.dat"]["rows"] == 1000
    assert len(golden["curves"]["dos.dat"]["y"]) <= tool.MAX_ROWS
    assert golden["scalars"] == {"spin_chern.json:spin_chern": 2.0, "spin_chern.json:z2": 1}


def test_rounding_noise_passes_and_a_real_change_fails(tmp_path):
    run(["make", write_output(tmp_path / "ref"), tmp_path / "g.json"])
    assert run(["check", write_output(tmp_path / "noise", scale=1 + 1e-10), tmp_path / "g.json"]) == 0
    assert run(["check", write_output(tmp_path / "moved", scale=1 + 1e-5), tmp_path / "g.json"]) == 1
    assert run(["check", write_output(tmp_path / "moved2", scale=1 + 1e-5), tmp_path / "g.json", 1e-4]) == 0


def test_the_oracle_has_its_own_looser_tolerance(tmp_path):
    run(["make", write_output(tmp_path / "ref"), tmp_path / "g.json"])
    assert run(["check", write_output(tmp_path / "a", chern=2.0 + 2e-6), tmp_path / "g.json"]) == 0
    assert run(["check", write_output(tmp_path / "b", chern=2.0 + 2e-6), tmp_path / "g.json", 1e-9, 1e-7]) == 1


def test_a_change_between_the_sampled_rows_is_seen_through_the_integral(tmp_path, capsys):
    run(["make", write_output(tmp_path / "ref"), tmp_path / "g.json"])
    stride = json.loads((tmp_path / "g.json").read_text())["curves"]["dos.dat"]["stride"]
    assert stride > 1
    capsys.readouterr()
    assert run(["check", write_output(tmp_path / "bump", bump_row=stride * 100 + 1), tmp_path / "g.json"]) == 1
    assert "dos.dat integral" in capsys.readouterr().out


@pytest.mark.parametrize("change", [dict(z2=0), dict(chern=2.001)])
def test_oracle_numbers_are_compared(tmp_path, change):
    run(["make", write_output(tmp_path / "ref"), tmp_path / "g.json"])
    assert run(["check", write_output(tmp_path / "other", **change), tmp_path / "g.json"]) == 1


def test_a_missing_curve_or_oracle_file_is_a_failure(tmp_path, capsys):
    run(["make", write_output(tmp_path / "ref"), tmp_path / "g.json"])
    out = write_output(tmp_path / "partial")
    (out / "spin_chern.json").unlink()
    assert run(["check", out, tmp_path / "g.json"]) == 1
    assert "structure:" in capsys.readouterr().out
    empty = tmp_path / "empty"
    empty.mkdir()
    assert run(["check", empty, tmp_path / "g.json"]) == 1
