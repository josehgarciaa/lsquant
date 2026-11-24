"""Thread binding at import and the workspace bound of the Bloch sum.

- Importing the package must not bind the process to one core: core binding
  (``OMP_PROC_BIND``/``OMP_PLACES``) is set only when the process is the
  ``lsquant`` command and has its own cores (a scheduler CPU set) or
  ``LSQUANT_BIND_THREADS=1``; caller values always win.
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


@pytest.mark.parametrize("cli, confined, choice, bound", [
    (False, True, None, False),    # a plain import never binds
    (False, True, "1", False),
    (True, True, None, True),      # the command with its own cores binds
    (True, False, None, False),    # the command on a shared node does not
    (True, False, "1", True),      # explicit request
    (True, True, "0", False),      # explicit refusal
])
def test_binding_only_for_a_confined_command(monkeypatch, cli, confined, choice, bound):
    for key in ("OMP_PROC_BIND", "OMP_PLACES", "OMP_NUM_THREADS", "LSQUANT_BIND_THREADS"):
        monkeypatch.delenv(key, raising=False)
    if choice is not None:
        monkeypatch.setenv("LSQUANT_BIND_THREADS", choice)
    monkeypatch.setattr(runtime, "owns_its_cores", lambda: confined)
    argv = ["python", "-m", "lsquant", "run"] if cli else ["python", "script.py"]
    runtime.configure_threads(argv)
    assert (os.environ.get("OMP_PROC_BIND") == "spread") is bound
    assert (os.environ.get("OMP_PLACES") == "cores") is bound
    assert int(os.environ["OMP_NUM_THREADS"]) >= 1


def test_owns_its_cores_compares_affinity_with_the_machine(monkeypatch):
    monkeypatch.setattr(os, "cpu_count", lambda: 8)
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: set(range(8)), raising=False)
    assert runtime.owns_its_cores() is False
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: {2, 3}, raising=False)
    assert runtime.owns_its_cores() is True


def test_bad_binding_choice_is_refused(monkeypatch):
    monkeypatch.setenv("LSQUANT_BIND_THREADS", "maybe")
    with pytest.raises(ValueError, match="LSQUANT_BIND_THREADS"):
        runtime.binding_requested(["python", "-m", "lsquant", "run"])


def test_caller_values_win(monkeypatch):
    monkeypatch.setenv("OMP_PROC_BIND", "false")
    monkeypatch.setenv("OMP_NUM_THREADS", "3")
    monkeypatch.setenv("LSQUANT_BIND_THREADS", "1")
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
