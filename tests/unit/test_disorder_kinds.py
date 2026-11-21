"""Unit tests of the selectable Anderson disorder (basis.disorder.kind).

They need the compiled module ``lsquant._core`` (see test_python_stage2.py); without it
every test is skipped.

What is checked:
  * spinor_disorder_values: ``scalar`` puts one value per orbital on both spin components,
    ``magnetic`` puts +h on spin up and -h on spin down, ``independent`` is the previous
    one-draw-per-coefficient layout; every value lies in [-W/2, W/2]; an odd dimension is
    refused for scalar and magnetic;
  * anderson_potential with an orthonormal spinor overlap: the scalar kind gives u (x) I_2,
    the magnetic kind h (x) sigma_z;
  * on the rotated Kane-Mele fixture of the stage-6 checks (interleaved spinors), the
    disordered Bloch and real-space bases give the same DOS moments for the magnetic kind,
    and the scalar and magnetic kinds give different moments;
  * the run file: kind defaults to scalar, unknown kinds are refused, and magnetic or
    independent kinds are refused for models without the SIESTA spinor layout.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
for p in (os.path.join(REPO, "python"), os.path.join(REPO, "tests", "dev")):
    if p not in sys.path:
        sys.path.insert(0, p)

pytest.importorskip("lsquant", reason="lsquant._core not built")
sparse = pytest.importorskip("scipy.sparse")
from lsquant.bases import make_anderson_disorder  # noqa: E402
from lsquant.disorder_bases import spinor_disorder_values  # noqa: E402
from lsquant.models.disorder import anderson_potential  # noqa: E402

N, W, SEED = 40, 1.0, 7


def test_scalar_values():
    w = spinor_disorder_values(N, W, SEED, "scalar").real
    assert np.array_equal(w[0::2], w[1::2])
    assert np.all(np.abs(w) <= W / 2)


def test_magnetic_values():
    w = spinor_disorder_values(N, W, SEED, "magnetic").real
    assert np.array_equal(w[0::2], -w[1::2])
    assert np.all(np.abs(w) <= W / 2) and np.any(w != 0)


def test_independent_values_are_the_previous_layout():
    assert np.array_equal(spinor_disorder_values(N, W, SEED, "independent"), make_anderson_disorder(N, W, SEED))


def test_scalar_and_magnetic_share_the_orbital_draw():
    s = spinor_disorder_values(N, W, SEED, "scalar").real
    m = spinor_disorder_values(N, W, SEED, "magnetic").real
    assert np.array_equal(s[0::2], m[0::2])


def test_odd_dimension_and_unknown_kind_refused():
    with pytest.raises(ValueError):
        spinor_disorder_values(N + 1, W, SEED, "magnetic")
    with pytest.raises(ValueError):
        spinor_disorder_values(N, W, SEED, "random_field")


def test_potentials_in_an_orthonormal_spinor_basis():
    s = sparse.identity(N, dtype=complex, format="csr")
    ws = spinor_disorder_values(N, W, SEED, "scalar")
    wm = spinor_disorder_values(N, W, SEED, "magnetic")
    u = ws.real[0::2]
    sigma_z = np.diag([1.0, -1.0])
    assert np.allclose(anderson_potential(s, ws).toarray(), np.kron(np.diag(u), np.eye(2)))
    assert np.allclose(anderson_potential(s, wm).toarray(), np.kron(np.diag(u), sigma_z))


# ------------------------------------------------------------------ bases on the KM fixture
@pytest.fixture(scope="module")
def km():
    pytest.importorskip("stage_6c_check")
    from stage_6c_check import fixtures
    return fixtures()[1]


def dos_moments(basis, M=16):
    from lsquant.apps.stage6b_siesta import compute_moments, DOS_WORD
    mu, _ = compute_moments(basis, DOS_WORD, moments=M, trace="random", states=1, seed=12345)
    return np.asarray(mu.data)


def test_magnetic_disorder_bloch_equals_real(km):
    from lsquant.disorder_bases import DisorderedBlochMetricBasis, RealSpaceMetricBasis
    mesh = [4, 4, 1]
    bloch = DisorderedBlochMetricBasis(km, mesh, disorder_w=1.0, disorder_seed=7, disorder_kind="magnetic",
                                       velocity_mode="covariant", route="lowdin")
    real = RealSpaceMetricBasis(km, mesh, disorder_w=1.0, disorder_seed=7, disorder_kind="magnetic",
                                velocity_mode="covariant", route="lowdin", applier="exact")
    assert np.array_equal(bloch.disorder_values, real.disorder_values)
    a, b = dos_moments(bloch), dos_moments(real)
    assert np.max(np.abs(a - b)) < 1e-10
    scalar = DisorderedBlochMetricBasis(km, mesh, disorder_w=1.0, disorder_seed=7, disorder_kind="scalar",
                                        velocity_mode="covariant", route="lowdin")
    assert np.max(np.abs(dos_moments(scalar) - a)) > 1e-6


# ------------------------------------------------------------------ run file
def g1():
    with open(os.path.join(REPO, "runs", "goldens", "g1.json")) as fh:
        return json.load(fh)


def test_runfile_default_and_refusals():
    from lsquant.runfile import RunFile, validate_runfile_dict
    rf = RunFile.from_dict(g1(), repo_root=REPO)
    assert rf.basis["disorder"]["kind"] == "scalar"
    for kind in ("magnetic", "independent", "unknown"):
        raw = g1()
        raw["basis"]["disorder"] = {"W": 1.0, "seed": 7, "kind": kind}
        with pytest.raises(ValueError):
            validate_runfile_dict(raw)
