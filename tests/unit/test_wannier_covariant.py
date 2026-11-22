"""Covariant Wannier velocity reads ``_r.dat`` and caches it on a frozen model.

The fixture is the two-orbital graphene model of the stage-0 goldens. Its
``_r.dat`` is written here with diagonal positions equal to the Wannier
centres at R = 0 and no off-diagonal elements. For that position matrix the
commutator term only converts the cell derivative into the atom gauge, so the
covariant velocity must equal the ``berry_connection`` velocity.
"""
import os
import shutil
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MODEL = os.path.join(REPO, "goldens", "stage0", "bases", "graphene", "model")
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant.models.wannier import load_wannier_model  # noqa: E402


def _write_r_dat(model_dir, seed, centres):
    """Write a ``_r.dat`` with ``<m|r|n> = delta_mn centre_n`` at R = 0."""
    with open(os.path.join(model_dir, seed + "_hr.dat"), encoding="utf-8") as fh:
        lines = [ln.split() for ln in fh.read().splitlines()[1:] if ln.strip()]
    width, count = int(lines[0][0]), int(lines[1][0])
    index, seen = 2, 0
    while seen < count:
        seen += len(lines[index])
        index += 1
    r_list = list(dict.fromkeys(tuple(map(int, ln[:3])) for ln in lines[index:]))
    out = ["written by test_wannier_covariant", str(width), str(count)]
    for r in r_list:
        for n in range(width):
            for m in range(width):
                pos = centres[m] if (r == (0, 0, 0) and m == n) else np.zeros(3)
                vals = " ".join(f"{x:.10f} 0.0" for x in pos)
                out.append(f"{r[0]} {r[1]} {r[2]} {m + 1} {n + 1} {vals}")
    with open(os.path.join(model_dir, seed + "_r.dat"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")


@pytest.fixture()
def model_dir(tmp_path):
    target = tmp_path / "graphene"
    shutil.copytree(MODEL, target)
    centres = load_wannier_model(str(target), seed="graphene").tau
    _write_r_dat(str(target), "graphene", centres)
    return str(target)


@pytest.mark.parametrize("direction", ["x", "y"])
def test_covariant_runs_on_frozen_model_and_caches(model_dir, direction):
    model = load_wannier_model(model_dir, seed="graphene")
    k = np.array([[0.3, 0.1, 0.0], [-0.7, 0.4, 0.0]])
    first = model.velocity_blocks(k, direction, mode="covariant")
    assert hasattr(model, "_position_cache")
    second = model.velocity_blocks(k, direction, mode="covariant")
    np.testing.assert_array_equal(first, second)


@pytest.mark.parametrize("direction", ["x", "y"])
def test_covariant_with_centre_positions_equals_berry_connection(model_dir, direction):
    model = load_wannier_model(model_dir, seed="graphene")
    rng = np.random.default_rng(7)
    k = rng.uniform(-2.0, 2.0, size=(6, 3))
    k[:, 2] = 0.0
    cov = model.velocity_blocks(k, direction, mode="covariant")
    bc = model.velocity_blocks(k, direction, mode="berry_connection")
    np.testing.assert_allclose(cov, bc, atol=1e-10)
    np.testing.assert_allclose(cov, np.conj(np.swapaxes(cov, -1, -2)), atol=1e-12)
