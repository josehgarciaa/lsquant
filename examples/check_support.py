"""Small numerical gates for saved example outputs, without running simulations.

Each gate reports a measured value, an expectation and an absolute tolerance.
Missing, ambiguous, non-finite or malformed input is a failure, never a reason
to regenerate data. Work is linear in the small curve/moment tables loaded.

Pitfall: importing an existing plotting helper can initialise Matplotlib. The
CLI confines its cache to temporary scratch and disables Python bytecode writes.
Set TMPDIR to the test's tmp_path to confine all scratch to that directory.
Warning: these checks assume the published run configurations, not arbitrary
parameter scans. To-do: validate runtime and memory limits with a built package.
"""
import argparse
import os
from pathlib import Path
import sys
import tempfile

# These gates have a one-core contract even when the surrounding shell was
# configured for production runs. Set this before NumPy or the engine loads.
for _variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                  "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_variable] = "1"

import numpy as np


class Checks:
    """Accumulate failures so all available numerical discrepancies are printed."""

    def __init__(self):
        self.ok = True

    def near(self, label, value, expected, tolerance, reason):
        """Use absolute error; callers supply units and the physical error budget."""
        good = bool(np.isfinite(value) and abs(value - expected) <= tolerance)
        self.ok &= good
        print(f"{label}: {value:.12g}; expected {expected:.12g} +/- {tolerance:g} "
              f"({reason}) {'PASS' if good else 'FAIL'}")

    def require(self, label, condition):
        """Report structural conditions which have no continuous tolerance."""
        self.near(label, float(bool(condition)), 1.0, 0.0, "required condition")


def one(root, pattern):
    """Reject ambiguous files instead of silently testing an unrelated run."""
    paths = sorted(Path(root).glob(pattern))
    if len(paths) != 1:
        raise ValueError(f"expected one {pattern} in {root}, found {len(paths)}")
    return paths[0]


def table(path):
    """Read a finite, increasing two-column energy/time curve, in stored units."""
    data = np.loadtxt(path, ndmin=2)
    if (data.shape[0] < 2 or data.shape[1] != 2 or
            not np.isfinite(data).all() or not np.all(np.diff(data[:, 0]) > 0)):
        raise ValueError(f"invalid two-column curve: {path}")
    return data


def same_grid(*curves):
    """Require identical stored coordinates before subtracting sampled curves."""
    if not all(np.array_equal(curves[0][:, 0], d[:, 0]) for d in curves[1:]):
        raise ValueError("curve grids differ")


def relative_error(actual, expected):
    """Compare entire complex arrays with a global reference scale, not per entry."""
    actual, expected = np.asarray(actual), np.asarray(expected)
    if actual.shape != expected.shape or not actual.size:
        raise ValueError("empty or differently shaped arrays")
    if not np.isfinite(actual).all() or not np.isfinite(expected).all():
        raise ValueError("non-finite array")
    scale = float(np.max(np.abs(expected)))
    return float(np.max(np.abs(actual - expected))) / (scale if scale else 1.0)


def run(check, subdir=None, configure=None):
    """Run a read-only check with failure status 1, including invalid input.

    Temporary cache writes are isolated and removed; no output is regenerated.
    Imported numerical helpers retain their original module paths and behavior.
    """
    parser = argparse.ArgumentParser(description=check.__doc__)
    parser.add_argument("out_dir", type=Path)
    if configure:
        configure(parser)
    args = parser.parse_args()
    root = args.out_dir
    if subdir and (root / subdir).is_dir():
        root = root / subdir
    checks = Checks()
    old_cache = os.environ.get("MPLCONFIGDIR")
    old_bytecode = sys.dont_write_bytecode
    try:
        with tempfile.TemporaryDirectory(prefix="example-check-") as scratch:
            os.environ["MPLCONFIGDIR"] = scratch
            sys.dont_write_bytecode = True
            if not root.is_dir():
                raise ValueError(f"not an output directory: {root}")
            check(root, checks, args)
    except (Exception, SystemExit) as exc:
        checks.ok = False
        print(f"input/check error: {exc}")
    finally:
        sys.dont_write_bytecode = old_bytecode
        if old_cache is None:
            os.environ.pop("MPLCONFIGDIR", None)
        else:
            os.environ["MPLCONFIGDIR"] = old_cache
    print("CHECK PASS" if checks.ok else "CHECK FAIL")
    return 0 if checks.ok else 1
