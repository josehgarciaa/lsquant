"""Thread binding at import and the workspace bound of the Bloch sum.

- Importing the package must not bind the process to one core: core binding
  (``OMP_PROC_BIND``/``OMP_PLACES``) is set only when the process is the
  ``lsquant`` command, and caller values always win.
- ``bloch_sum`` with a workspace bound too small for one k point gives the
  same result as with a generous bound (one k point at a time), instead of
  stopping.
"""
import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "python"))

from lsquant import runtime  # noqa: E402
from lsquant.models.bloch import bloch_sum  # noqa: E402


@pytest.mark.parametrize("argv, expected", [
    (["python", "-m", "lsquant", "run", "x.json"], True),
    (["python", "-m", "lsquant.cli", "run"], True),
    (["/usr/bin/python3", "/opt/env/bin/lsquant", "run", "x.json"], True),
    (["python", "my_script.py"], False),
    (["python", "-c", "import lsquant"], False),
    (["python", "-m", "pytest", "tests/unit"], False),
    (["python", "examples/20/openmx_exact.py", "lsquant"], False),
])
def test_cli_detection(argv, expected):
    assert runtime.launched_as_cli(argv) is expected


@pytest.mark.parametrize("cli", [False, True])
def test_binding_only_for_the_command(monkeypatch, cli):
    for key in ("OMP_PROC_BIND", "OMP_PLACES", "OMP_NUM_THREADS"):
        monkeypatch.delenv(key, raising=False)
    argv = ["python", "-m", "lsquant", "run"] if cli else ["python", "script.py"]
    runtime.configure_threads(argv)
    assert (os.environ.get("OMP_PROC_BIND") == "spread") is cli
    assert (os.environ.get("OMP_PLACES") == "cores") is cli
    assert int(os.environ["OMP_NUM_THREADS"]) >= 1


def test_caller_values_win(monkeypatch):
    monkeypatch.setenv("OMP_PROC_BIND", "false")
    monkeypatch.setenv("OMP_NUM_THREADS", "3")
    runtime.configure_threads(["python", "-m", "lsquant", "run"])
    assert os.environ["OMP_PROC_BIND"] == "false"
    assert os.environ["OMP_NUM_THREADS"] == "3"


@pytest.mark.parametrize("derivative", [None, 0])
def test_bloch_sum_tiny_workspace_matches(derivative):
    rng = np.random.default_rng(1)
    nr, width = 7, 40
    r_cart = rng.normal(size=(nr, 3))
    centres = rng.normal(size=(width, 3))
    blocks = rng.normal(size=(nr, width, width)) + 1j * rng.normal(size=(nr, width, width))
    k = rng.normal(size=(5, 3))
    big = bloch_sum(r_cart, centres, blocks, k, derivative=derivative, memory_mb=256.0)
    tiny = bloch_sum(r_cart, centres, blocks, k, derivative=derivative, memory_mb=1e-3)
    assert np.allclose(big, tiny, atol=1e-12)
