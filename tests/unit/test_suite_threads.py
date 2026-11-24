"""The unit suite must not start a machine-wide OpenMP team (see tests/conftest.py)."""
from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

CONFTEST = Path(__file__).resolve().parents[1] / "conftest.py"


def load_conftest():
    spec = importlib.util.spec_from_file_location("suite_conftest", CONFTEST)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_default_thread_count_is_bounded():
    conf = load_conftest()
    assert conf.default_thread_count(1) == 1
    assert conf.default_thread_count(2) == 2
    assert conf.default_thread_count(64) == conf.MAX_DEFAULT_THREADS
    assert conf.default_thread_count(0) == 1


def run_conftest(env):
    code = ("import importlib.util, os; "
            f"s = importlib.util.spec_from_file_location('c', {str(CONFTEST)!r}); "
            "m = importlib.util.module_from_spec(s); s.loader.exec_module(m); "
            "print(os.environ['OMP_NUM_THREADS'])")
    return subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True,
                          check=True).stdout.strip()


def test_unset_variable_gets_the_bounded_default():
    env = {k: v for k, v in os.environ.items() if k != "OMP_NUM_THREADS"}
    assert 1 <= int(run_conftest(env)) <= load_conftest().MAX_DEFAULT_THREADS


def test_caller_value_wins():
    assert run_conftest(dict(os.environ, OMP_NUM_THREADS="7")) == "7"


def test_suite_runs_with_a_set_thread_count():
    # conftest.py is loaded before this module, so the variable is always set here.
    assert int(os.environ["OMP_NUM_THREADS"]) >= 1
