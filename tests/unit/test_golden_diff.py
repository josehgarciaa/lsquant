"""Unit tests for scripts/golden_diff.py, the numeric localiser used by the golden gate.

golden_diff.py is called by scripts/check_goldens.sh only after `cmp` has found a byte
difference. Its job is to say *where* the first difference is (in the units of the file
layout: moment index m, (m, n) pair, (m, t) pair, or row/column) and *how large* the
largest difference is, and to classify it as formatting only, rounding level
(< 1e-12 relative) or defect level. These tests build small files of each layout in a
temporary directory and check the printed line and the exit code.

Run with:  python3 -m pytest tests/unit
"""
import importlib.util
import math
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "golden_diff.py"

# Import the script as a module so its helpers can be tested directly.
_spec = importlib.util.spec_from_file_location("golden_diff", SCRIPT)
golden_diff = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(golden_diff)


def run(new_text, gold_text, tmp_path, ext=".dat"):
    """Write the two files, run golden_diff.py on them, return (exit code, stdout)."""
    new = tmp_path / ("new" + ext)
    gold = tmp_path / ("gold" + ext)
    new.write_text(new_text)
    gold.write_text(gold_text)
    proc = subprocess.run([sys.executable, str(SCRIPT), str(new), str(gold)],
                          capture_output=True, text=True)
    return proc.returncode, proc.stdout.strip()


# ----------------------------------------------------------------- helper functions
def test_rel_basic():
    d, r = golden_diff._rel(1.0, 1.0 + 1e-13)
    assert d == pytest.approx(1e-13, rel=1e-3)
    assert r == pytest.approx(1e-13, rel=1e-3)


def test_rel_both_zero_is_zero():
    assert golden_diff._rel(0.0, 0.0) == (0.0, 0.0)


def test_rel_nan_equals_nan():
    assert golden_diff._rel(math.nan, math.nan) == (0.0, 0.0)


def test_rel_nan_against_finite_is_infinite():
    d, r = golden_diff._rel(math.nan, 1.0)
    assert math.isinf(d) and math.isinf(r)


def test_rel_same_infinity_is_zero():
    assert golden_diff._rel(math.inf, math.inf) == (0.0, 0.0)


def test_num_parses_and_rejects():
    assert golden_diff._num("1.5e-3") == pytest.approx(1.5e-3)
    assert golden_diff._num("abc") is None


# ----------------------------------------------------------------- classification
def test_identical_values_formatting_only(tmp_path):
    code, out = run("1.0 2.0\n", "1.00 2.0\n", tmp_path)
    assert code == 1
    assert "identical values (formatting only)" in out


def test_rounding_level(tmp_path):
    code, out = run("1.0000000000000002 2.0\n", "1.0 2.0\n", tmp_path)
    assert code == 1
    assert "rounding level (< 1e-12)" in out


def test_defect_level(tmp_path):
    code, out = run("1.001 2.0\n", "1.0 2.0\n", tmp_path)
    assert code == 1
    assert "DEFECT level (>= 1e-12)" in out
    assert "first-diff line 1 col 1" in out


def test_only_line_endings_differ(tmp_path):
    code, out = run("1.0 2.0\n", "1.0 2.0", tmp_path)
    assert code == 0
    assert "no line differs" in out


# ----------------------------------------------------------------- structure errors
def test_line_count_differs(tmp_path):
    code, out = run("1\n2\n", "1\n2\n3\n", tmp_path)
    assert code == 1
    assert "STRUCTURE differs: 2 lines new vs 3 golden" in out


def test_token_count_differs(tmp_path):
    code, out = run("1 2\n", "1 2 3\n", tmp_path)
    assert code == 1
    assert "STRUCTURE differs at line 1" in out


def test_non_numeric_token(tmp_path):
    code, out = run("abc 2\n", "abd 2\n", tmp_path)
    assert code == 1
    assert "NON-NUMERIC token differs" in out


# ----------------------------------------------------------------- moment layouts
HEADER = "2048 4.0 0.0\n"


def test_chebmom1d_index_is_m(tmp_path):
    gold = HEADER + "3\n" + "1.0 0.0\n0.5 0.0\n0.25 0.0\n"
    new = HEADER + "3\n" + "1.0 0.0\n0.5 0.0\n0.30 0.0\n"
    code, out = run(new, gold, tmp_path, ext=".chebmom1D")
    assert code == 1
    assert "m=2 re" in out


def test_chebmom2d_index_is_m_n(tmp_path):
    # M0 = 2, M1 = 3: flat index k = m * 3 + n; change k = 4 -> (m, n) = (1, 1), imag part.
    rows = ["%d.0 0.0" % k for k in range(6)]
    gold = HEADER + "2 3\n" + "\n".join(rows) + "\n"
    rows[4] = "4.0 1.0"
    new = HEADER + "2 3\n" + "\n".join(rows) + "\n"
    code, out = run(new, gold, tmp_path, ext=".chebmom2D")
    assert code == 1
    assert "(m,n)=(1,1) im" in out


def test_chebmomtd_index_is_m_t(tmp_path):
    # M = 2, Nt = 4: flat index k = m * 4 + t; change k = 6 -> (m, t) = (1, 2).
    rows = ["%d.0 0.0" % k for k in range(8)]
    gold = "2048 4.0 0.0 4 0.1\n2 4\n" + "\n".join(rows) + "\n"
    rows[6] = "7.0 0.0"
    new = "2048 4.0 0.0 4 0.1\n2 4\n" + "\n".join(rows) + "\n"
    code, out = run(new, gold, tmp_path, ext=".chebmomTD")
    assert code == 1
    assert "(m,t)=(1,2) re" in out


def test_header_difference_reported_by_line_and_column(tmp_path):
    gold = HEADER + "1\n1.0 0.0\n"
    new = "2048 4.5 0.0\n1\n1.0 0.0\n"
    code, out = run(new, gold, tmp_path, ext=".chebmom1D")
    assert code == 1
    assert "first-diff line 1 col 2" in out


def test_usage_without_arguments():
    proc = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True)
    assert proc.returncode == 2
    assert "Usage" in proc.stdout
