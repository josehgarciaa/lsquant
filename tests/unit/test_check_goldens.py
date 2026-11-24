"""Unit tests for the golden gate scripts/check_goldens.sh itself.

The gate's stage 0 runs six legacy drivers and compares their output files with
goldens/stage0. Here the drivers are replaced by small fake executables that write the
golden files (or a corrupted copy) into the working directory, so the gate's own logic is
tested in seconds and without a build:

  * all outputs identical            -> exit 0, "GATE PASS (stage 0)"
  * one output corrupted             -> exit 1, "GATE FAIL", a DIFFERS line with the
                                        golden_diff.py localisation
  * one output not produced          -> exit 1, "(not produced: ...)"
  * a driver missing from build_dir  -> exit 2 before any run
  * unknown stage, missing build dir -> exit 2
  * the gate never writes into goldens/
  * tolerance mode (GOLDEN_TOL=1): a curve off by one unit of its last printed digit,
    a curve value at the noise level (below 1e-12 of the curve maximum) and moments off
    by 1e-13 relative pass as TOLPASS; five units, or 1e-9 relative on the moments,
    still fail; without GOLDEN_TOL all of them are DIFFERS

Run with:  python3 -m pytest tests/unit
"""
import hashlib
import os
import pathlib
import shlex
import stat
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parents[2]
GATE = REPO / "scripts" / "check_goldens.sh"
GOLD = REPO / "goldens" / "stage0"

DRIVERS = [
    "inline_compute-kpm-nonEqOp",
    "inline_compute-kpm-spectralOp",
    "inline_kuboBastinFromChebmom",
    "inline_kuboBastinSeaFromChebmom",
    "inline_kuboBastinSurfFromChebmom",
    "inline_kuboGreenwoodFromChebmom",
]

# Fake drivers. They reproduce the legacy output file names from the arguments and copy
# the matching golden into the current directory. FAKE_MODE=corrupt changes one digit of
# every produced chain1d moment file; FAKE_MODE=skip makes spectralOp produce nothing.
FAKE_NONEQ = r"""#!/usr/bin/env bash
# args: <label> <op1> <op2> <M> <state>
f="NonEqOp$2-$3$1KPM_M$4x$4_state$5.chebmom2D"
cp "$FAKE_GOLD/$1/$f" "$f"
if [ "${FAKE_MODE:-}" = corrupt ] && [ "$1" = chain1d ]; then
  awk 'NR == 3 { $1 = $1 + 0.5 } { print }' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
fi
"""
FAKE_SPECTRAL = r"""#!/usr/bin/env bash
# args: <label> <op> <M> <state>
[ "${FAKE_MODE:-}" = skip ] && exit 0
f="SpectralOp$2$1KPM_M$3_state$4.chebmom1D"
cp "$FAKE_GOLD/$1/$f" "$f"
"""
# Every Kubo driver copies all Haldane .dat goldens; the gate only checks their presence
# and content, not which driver wrote which file. FAKE_MODE=dat1 / dat5 shifts one
# conductivity value of the Kubo-Greenwood curve by 1 or 5 units of its last printed digit;
# FAKE_MODE=noise moves its smallest value by 3e-13 of the curve's maximum.
FAKE_KUBO = r"""#!/usr/bin/env bash
cp "$FAKE_GOLD"/haldane/*.dat .
case "${FAKE_MODE:-}" in
  dat1|dat5) python3 "$FAKE_PERTURB" curve KuboGreenWood_*.dat "${FAKE_MODE#dat}" ;;
  noise)     python3 "$FAKE_PERTURB" noise KuboGreenWood_*.dat 0 ;;
esac
"""
# FAKE_MODE=mom13 / mom9 scales every chain1d moment by (1 + 1e-13) or (1 + 1e-9).
FAKE_NONEQ_TOL = r"""#!/usr/bin/env bash
f="NonEqOp$2-$3$1KPM_M$4x$4_state$5.chebmom2D"
cp "$FAKE_GOLD/$1/$f" "$f"
case "${FAKE_MODE:-}" in
  mom13) [ "$1" = chain1d ] && python3 "$FAKE_PERTURB" moments "$f" 1e-13 ;;
  mom9)  [ "$1" = chain1d ] && python3 "$FAKE_PERTURB" moments "$f" 1e-9 ;;
esac
exit 0
"""
# Helper used by the fake drivers to perturb a produced file in place.
PERTURB = r"""
import math, sys
kind, path, amount = sys.argv[1], sys.argv[2], sys.argv[3]
lines = open(path).read().split("\n")
def column(lines):
    # (line index, value) of every row with a numeric second column
    out = []
    for i, line in enumerate(lines):
        t = line.split()
        if len(t) >= 2:
            out.append((i, float(t[1])))
    return out
if kind == "curve":
    # the largest |value| of the curve: shift it by `amount` units of its last printed
    # (6th) significant digit, printed again with 6 significant digits
    i, x = max(column(lines), key=lambda p: abs(p[1]))
    unit = 10.0 ** (math.floor(math.log10(abs(x))) - 5)
    t = lines[i].split()
    t[1] = "%g" % (x + int(amount) * unit)
    lines[i] = " ".join(t)
elif kind == "noise":
    # the smallest nonzero |value|: add 3e-13 of the curve's largest |value| (below the
    # 1e-12 noise floor of the tolerance rule, far above one unit of that small value)
    col = column(lines)
    gmax = max(abs(v) for _, v in col)
    i, x = min((p for p in col if p[1] != 0.0), key=lambda p: abs(p[1]))
    t = lines[i].split()
    t[1] = "%g" % (x + 3e-13 * gmax)
    lines[i] = " ".join(t)
else:
    # scale every moment value (lines 3..) by (1 + amount), 15 significant digits
    f = 1.0 + float(amount)
    for i in range(2, len(lines)):
        t = lines[i].split()
        if t:
            lines[i] = " ".join("%.15g" % (float(v) * f) for v in t)
open(path, "w").write("\n".join(lines))
"""


def _write_exe(path, text):
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@pytest.fixture
def fake_build(tmp_path):
    build = tmp_path / "build"
    build.mkdir()
    _write_exe(build / "inline_compute-kpm-nonEqOp", FAKE_NONEQ)
    _write_exe(build / "nonEq_tol", FAKE_NONEQ_TOL)
    _write_exe(build / "inline_compute-kpm-spectralOp", FAKE_SPECTRAL)
    for name in DRIVERS[2:]:
        _write_exe(build / name, FAKE_KUBO)
    return build


def run_gate(args, mode="", tmp_path=None, tol=False):
    env = dict(os.environ, FAKE_GOLD=str(GOLD), FAKE_MODE=mode,
               GOLDEN_TOL="1" if tol else "0")
    if tmp_path is not None:
        perturb = tmp_path / "perturb.py"
        perturb.write_text(PERTURB)
        env["FAKE_PERTURB"] = str(perturb)
    if tmp_path is not None:
        env["TMPDIR"] = str(tmp_path)
    return subprocess.run(["bash", str(GATE)] + [str(a) for a in args],
                          capture_output=True, text=True, env=env)


def _tree_digest(root):
    h = hashlib.sha256()
    for p in sorted(root.rglob("*")):
        if p.is_file():
            h.update(str(p.relative_to(root)).encode())
            h.update(p.read_bytes())
    return h.hexdigest()


def test_pass_when_all_outputs_identical(fake_build, tmp_path):
    proc = run_gate(["0", fake_build], tmp_path=tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert "GATE PASS (stage 0)" in proc.stdout
    assert "DIFFERS" not in proc.stdout
    # chain1d at T = 16, 1, 64 (2 files each) and Haldane at T = 16 (2 moments + .dat files)
    assert proc.stdout.count("IDENTICAL") >= 8


def test_fail_and_localise_when_output_corrupted(fake_build, tmp_path):
    proc = run_gate(["0", fake_build], mode="corrupt", tmp_path=tmp_path)
    assert proc.returncode == 1
    assert "GATE FAIL (stage 0)" in proc.stdout
    assert "DIFFERS   [T=16] chain1d/NonEqOpVX-VX" in proc.stdout
    # the numeric localisation from golden_diff.py is printed under the DIFFERS line
    assert "first-diff line 3" in proc.stdout


@pytest.mark.parametrize("tol", [False, True])
def test_python_override_is_used_for_gate_helpers(fake_build, tmp_path, monkeypatch, tol):
    """Exercise both the numeric diagnostic and tolerance-classifier calls.

    Pitfall: a passing byte comparison never invokes Python in stage 0, so
    corrupt a moment file to ensure the override is actually exercised.
    """
    marker = tmp_path / "interpreter-calls"
    wrapper = tmp_path / "python wrapper"
    _write_exe(wrapper, '#!/usr/bin/env bash\n'
               'printf "%s\\n" "$1" >> "$PYTHON_CALLS"\n'
               f'exec {shlex.quote(sys.executable)} "$@"\n')
    monkeypatch.setenv("PYTHON", str(wrapper))
    monkeypatch.setenv("PYTHON_CALLS", str(marker))
    proc = run_gate(["0", fake_build], mode="corrupt", tmp_path=tmp_path, tol=tol)
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "first-diff line 3" in proc.stdout
    assert marker.exists(), "gate ignored PYTHON"
    calls = marker.read_text().splitlines()
    assert str(REPO / "scripts/golden_diff.py") in calls
    if tol:
        assert any(path.endswith("/tolcmp.py") for path in calls)


def test_fail_when_output_not_produced(fake_build, tmp_path):
    proc = run_gate(["0", fake_build], mode="skip", tmp_path=tmp_path)
    assert proc.returncode == 1
    assert "(not produced: SpectralOp1chain1dKPM_M64_stateexact.chebmom1D)" in proc.stdout


def test_missing_driver_stops_before_running(fake_build, tmp_path):
    (fake_build / "inline_kuboGreenwoodFromChebmom").unlink()
    proc = run_gate(["0", fake_build], tmp_path=tmp_path)
    assert proc.returncode == 2
    assert "legacy driver" in proc.stdout and "missing" in proc.stdout
    assert "stage 0 needs the legacy inline_* drivers" in proc.stdout
    assert "cmake --preset golden" in proc.stdout
    assert "IDENTICAL" not in proc.stdout


def test_unknown_stage(fake_build, tmp_path):
    proc = run_gate(["99", fake_build], tmp_path=tmp_path)
    assert proc.returncode == 2
    assert "no gate defined" in proc.stdout


def test_missing_build_dir(tmp_path):
    proc = run_gate(["0", tmp_path / "does-not-exist"], tmp_path=tmp_path)
    assert proc.returncode == 2
    assert "build dir" in proc.stdout


def test_no_arguments_prints_usage():
    proc = run_gate([])
    assert proc.returncode == 2
    assert "usage" in proc.stdout


def test_gate_never_writes_into_goldens(fake_build, tmp_path):
    before = _tree_digest(GOLD)
    run_gate(["0", fake_build], mode="corrupt", tmp_path=tmp_path)
    assert _tree_digest(GOLD) == before


# ----------------------------------------------------------------- tolerance mode
@pytest.fixture
def tol_build(fake_build):
    """Fake build whose nonEqOp driver can perturb the moments (FAKE_MODE=mom13|mom9)."""
    (fake_build / "inline_compute-kpm-nonEqOp").unlink()
    (fake_build / "nonEq_tol").rename(fake_build / "inline_compute-kpm-nonEqOp")
    return fake_build


def test_curve_one_unit_is_tolpass_only_in_tolerance_mode(tol_build, tmp_path):
    strict = run_gate(["0", tol_build], mode="dat1", tmp_path=tmp_path)
    assert strict.returncode == 1
    assert "DIFFERS   [T=16] haldane/KuboGreenWood" in strict.stdout
    tol = run_gate(["0", tol_build], mode="dat1", tmp_path=tmp_path, tol=True)
    assert tol.returncode == 0, tol.stdout
    assert "TOLPASS   [T=16] haldane/KuboGreenWood" in tol.stdout
    assert "tolerance-PASS" in tol.stdout and "GATE PASS (stage 0)" in tol.stdout


def test_curve_five_units_fails_in_tolerance_mode(tol_build, tmp_path):
    proc = run_gate(["0", tol_build], mode="dat5", tmp_path=tmp_path, tol=True)
    assert proc.returncode == 1
    assert "DIFFERS   [T=16] haldane/KuboGreenWood" in proc.stdout
    assert "beyond one unit" in proc.stdout


def test_moments_1e13_are_tolpass(tol_build, tmp_path):
    proc = run_gate(["0", tol_build], mode="mom13", tmp_path=tmp_path, tol=True)
    assert proc.returncode == 0, proc.stdout
    assert "TOLPASS   [T=16] chain1d/NonEqOpVX-VX" in proc.stdout


def test_moments_1e9_fail_in_tolerance_mode(tol_build, tmp_path):
    proc = run_gate(["0", tol_build], mode="mom9", tmp_path=tmp_path, tol=True)
    assert proc.returncode == 1
    assert "DIFFERS   [T=16] chain1d/NonEqOpVX-VX" in proc.stdout
    assert "<= 1e-12 required" in proc.stdout


def test_curve_noise_level_value_is_tolpass(tol_build, tmp_path):
    strict = run_gate(["0", tol_build], mode="noise", tmp_path=tmp_path)
    assert strict.returncode == 1
    tol = run_gate(["0", tol_build], mode="noise", tmp_path=tmp_path, tol=True)
    assert tol.returncode == 0, tol.stdout
    assert "TOLPASS   [T=16] haldane/KuboGreenWood" in tol.stdout
