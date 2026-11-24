"""Guard dependency floors required by source APIs and NumPy compatibility.

Read metadata as text so this check also runs on supported Python 3.10,
which lacks tomllib. No installed dependency version substitutes for the
published requirement. Pitfall: SciPy's NumPy support matters even when its
individual linear algebra and special-function APIs predate NumPy 2.0.
"""

from pathlib import Path
import re

REPO = Path(__file__).resolve().parents[2]


def _minimum(name):
    metadata = (REPO / "python/pyproject.toml").read_text()
    dependencies = re.search(r"dependencies\s*=\s*\[(.*?)\]", metadata,
                             re.DOTALL).group(1)
    match = re.search(rf'"{name}\s*>=\s*(\d+)\.(\d+)(?:\.\d+)?"',
                      dependencies)
    assert match, f"{name} needs an explicit minimum version"
    return tuple(map(int, match.groups()))


def test_numpy_api_minimum():
    """Require NumPy 2.0 for the committed trapezoid integrations."""
    source = (REPO / "python/lsquant/models/gpaw.py").read_text()
    assert "np.trapezoid(" in source
    assert _minimum("numpy") >= (2, 0)


def test_scipy_numpy_compatibility_minimum():
    """Require the SciPy 1.13 NumPy 2.0 compatibility floor."""
    assert _minimum("numpy") >= (2, 0)
    assert _minimum("scipy") >= (1, 13)
