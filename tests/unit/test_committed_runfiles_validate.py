"""Keep committed run inputs compatible with inspection without loading models.

The JSON branch of inspect_target constructs RunFile.from_dict and then
serializes defaults. Calling that branch also covers resolution and defaults,
without executing any model or writing outputs. The manifest comes from Git
so local scratch JSON cannot accidentally become a required input.
"""

from pathlib import Path
import subprocess
import sys

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))
pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant.runfile import inspect_target

RUNFILES = [p for p in subprocess.check_output(
    ["git", "ls-files", "-z", "--", "runs", "examples"], cwd=REPO
).decode().split("\0") if p and (
    (p.startswith("runs/") and p.endswith(".json")) or
    (p.startswith("examples/") and Path(p).match("run*.json")))]


def test_runfile_inventory():
    """Guard against a discovery that finds nothing: 116 inputs were committed when this was written.

    Pitfall: an exact count would fail on every new run file; the floor only catches a broken
    discovery (wrong directory, no Git metadata).
    """
    assert len(RUNFILES) >= 116


@pytest.mark.parametrize("relative", RUNFILES, ids=RUNFILES)
def test_committed_runfile_inspection(relative):
    """Accept each input through the CLI's read-only JSON inspection path."""
    info = inspect_target(str(REPO / relative))
    assert info["kind"] == "runfile"
    assert info["resolved_runfile"]["word"]
