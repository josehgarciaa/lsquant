"""Unit tests of the stage-2 Python layer (lsquant package over the pybind11 module).

They need the compiled module ``lsquant._core``: put ``<build>/lsquant_core`` and
``python`` on PYTHONPATH, for example

    PYTHONPATH=build/lsquant_core:python python3 -m pytest tests/unit/test_python_stage2.py

Without the module every test is skipped, so ``pytest tests/unit`` still runs on a machine
without a build.

What is checked:
  * jackson_damping against the closed form, g_0 = 1, input validation;
  * DeltaFunction.table (SymPy) against legacy_delta_table (legacy operation order);
  * dos_plan / kubo_plan: instruction lists, words, weights, input validation;
  * run_plan with the dos plan on a random Hermitian matrix against dense diagonalisation;
  * run_plan with both plans on the golden chain: identical to G1 and G2;
  * dos_curve on G1: byte-identical text to G3, weight inside the grid, input validation.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GOLD = os.path.join(REPO, "goldens", "stage0", "chain1d")
G1 = os.path.join(GOLD, "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D")
G2 = os.path.join(GOLD, "NonEqOpVX-VXchain1dKPM_M64x64_stateexact.chebmom2D")
G3 = os.path.join(GOLD, "meanSpectralOp1chain1dKPM_M64_stateexactJACKSON.dat")
OPS = os.path.join(GOLD, "operators")

if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

lsq = pytest.importorskip("lsquant", reason="lsquant._core not built (see module docstring)")
_core = lsq._core


# ------------------------------------------------------------------ helpers
def dense_to_core(h: np.ndarray):
    """Build a _core.SparseMatrix (CSR) from a dense complex matrix."""
    n = h.shape[0]
    rows, cols, vals = [0], [], []
    for i in range(n):
        nz = np.nonzero(h[i])[0]
        cols.extend(int(j) for j in nz)
        vals.extend(complex(h[i, j]) for j in nz)
        rows.append(len(cols))
    return _core.SparseMatrix.from_csr(
        n,
        np.asarray(rows, dtype=np.int64),
        np.asarray(cols, dtype=np.int64),
        np.asarray(vals, dtype=np.complex128),
    )


def random_hermitian(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return 0.5 * (a + a.conj().T)


# ------------------------------------------------------------------ kernels
def test_jackson_closed_form():
    m_tot = 37
    g = lsq.jackson_damping(m_tot)
    phi = math.pi / (m_tot + 1)
    m = np.arange(m_tot)
    ref = ((m_tot - m + 1) * np.cos(m * phi) + np.sin(m * phi) / np.tan(phi)) * phi / math.pi
    assert g.shape == (m_tot,)
    assert g[0] == pytest.approx(1.0, abs=1e-15)
    assert np.max(np.abs(g - ref)) < 1e-14
    assert np.all(g > 0) and np.all(g <= 1.0 + 1e-15)


def test_jackson_rejects_nonpositive():
    with pytest.raises(ValueError):
        lsq.jackson_damping(0)


def test_delta_table_sympy_matches_legacy():
    m_eff, alpha = 21, 0.95
    n = 30 * m_eff
    x = [-alpha + i * (2.0 * alpha) / (n - 1) for i in range(n)]
    sym = lsq.DeltaFunction.table(x, m_eff)
    leg = lsq.legacy_delta_table(x, m_eff)
    assert sym.shape == leg.shape
    assert float(np.max(np.abs(sym - leg))) < 1e-14


# ------------------------------------------------------------------ plans
def test_dos_plan_structure():
    plan = lsq.dos_plan(8, 4, word="1 delta")
    assert plan.word == "1 delta"
    assert len(plan.code) == 2


def test_kubo_plan_structure():
    plan = lsq.kubo_plan(8, left_index=0, right_index=0, word="VX delta VX delta")
    assert plan.word == "VX delta VX delta"
    assert len(plan.code) == 6


@pytest.mark.parametrize("args", [(0, 4), (8, 0), (-1, 1)])
def test_dos_plan_rejects_bad_input(args):
    with pytest.raises(ValueError):
        lsq.dos_plan(*args, word="1 delta")


def test_kubo_plan_rejects_bad_input():
    with pytest.raises(ValueError):
        lsq.kubo_plan(0, left_index=0, right_index=0, word="VX delta VX delta")


# ------------------------------------------------------------------ run_plan
def test_dos_moments_against_dense_diagonalisation():
    n, m_tot = 24, 16
    h = random_hermitian(n, seed=3)
    eps = np.linalg.eigvalsh(h)
    lo, hi = eps.min() - 0.5, eps.max() + 0.5
    window = _core.SpectralWindow(lo, hi)
    states = _core.LocalStates(n, list(range(n)))
    mu = lsq.run_plan(lsq.dos_plan(m_tot, n, word="1 delta"), dense_to_core(h), window, states)
    x = eps * window.scale + window.shift           # H~ = (H - c)/a: scale = 1/a, shift = -c/a
    assert np.all(np.abs(x) < 1.0)
    ref = np.array([np.sum(np.cos(m * np.arccos(x))) for m in range(m_tot)])
    ref *= np.where(np.arange(m_tot) == 0, 1.0, 2.0) / n
    assert mu.rank == 1 and list(mu.dims) == [m_tot]
    assert np.max(np.abs(mu.data - ref)) < 1e-12


def test_golden_chain_plans_are_identical_to_g1_g2():
    ham = _core.ReadCsrText(os.path.join(OPS, "chain1d.HAM.CSR"))
    vx = _core.SparseMatrixOperator(_core.ReadCsrText(os.path.join(OPS, "chain1d.VX.CSR")))
    dim = ham.dimension
    window = _core.SpectralWindow(-2.0, 2.0)
    states = _core.LocalStates(dim, list(range(dim)))
    mu_dos = lsq.run_plan(lsq.dos_plan(64, dim, word="1 delta"), ham, window, states, batch_size=64)
    mu_kubo = lsq.run_plan(lsq.kubo_plan(64, left_index=0, right_index=0, word="VX delta VX delta"), ham, window, states, operators=[vx],
                           batch_size=64, symmetrise=True)
    assert np.array_equal(mu_dos.data, _core.ReadChebmom(G1).data)
    assert np.array_equal(mu_kubo.data, _core.ReadChebmom(G2).data)


# ------------------------------------------------------------------ dos_curve
def test_dos_curve_reproduces_g3_text(tmp_path):
    e, v = lsq.dos_curve(_core.ReadChebmom(G1), 120.0)
    out = tmp_path / "g3.dat"
    _core.WriteCurveText(str(out), e, v)
    assert out.read_bytes() == open(G3, "rb").read()


def test_dos_curve_sum_rule():
    mu = _core.ReadChebmom(G1)
    e, v = lsq.dos_curve(mu, 20.0)
    # values are N * rho(E). G1 uses the window [-2, 2], equal to the band, and the grid
    # stops at |x| = 0.95: the exact weight inside is (2/pi) asin(0.95) = 0.798. The
    # Jackson broadening at 20 meV changes it by much less than 0.01.
    integral = float(np.sum(0.5 * (v[1:] + v[:-1]) * np.diff(e))) / mu.system_size
    assert abs(integral - 2.0 / math.pi * math.asin(0.95)) < 0.01


def test_dos_curve_energy_axis_for_centred_and_offset_windows():
    """dos_curve maps x to E = x / scale + shift = a x - c/a.

    For a window centred on 0 this is E = a x, the exact inverse of the rescaling, and the
    output equals the legacy one. For a window with centre c it gives a x - c/a instead of
    a x + c (the legacy drivers give a x + c/a); the time-evolution stage corrects the axis
    (golden G4c). The test pins both cases.
    """
    for lo, hi in ((-2.0, 2.0), (-1.0, 3.0)):
        mu = _core.ReadChebmom(G1)
        mu.window = _core.SpectralWindow(lo, hi)
        e, _ = lsq.dos_curve(mu, 120.0)
        a, c = 0.5 * (hi - lo), 0.5 * (hi + lo)
        m_eff = len(e) // 30
        x = np.array([-0.95 + i * 1.9 / (30 * m_eff - 1) for i in range(30 * m_eff)])
        assert np.max(np.abs(e - (a * x - c / a))) < 1e-12
        if c == 0.0:
            assert np.max(np.abs(e - (a * x + c))) < 1e-12


def test_dos_curve_rejects_bad_input():
    mu = _core.ReadChebmom(G1)
    with pytest.raises(ValueError):
        lsq.dos_curve(mu, 0.0)
    with pytest.raises(ValueError):
        lsq.dos_curve(mu, 10.0, alpha=1.0)
    with pytest.raises(ValueError):
        lsq.dos_curve(_core.ReadChebmom(G2), 10.0)
