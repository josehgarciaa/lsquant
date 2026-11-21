"""Unit tests of the stage-3 Python layer: rank-2 moments and the Kubo reconstructions.

They need the compiled module ``lsquant._core`` (see test_python_stage2.py); without it
every test is skipped.

What is checked:
  * kubo_plan with two different Hermitian operators on a random model against
    Tr[A T_m(H~) B T_n(H~)] from dense diagonalisation, with the legacy weights
    (2 - delta_m0)(2 - delta_n0)/N;
  * the symmetrised tensor satisfies mu(m, n) = conj(mu(n, m));
  * window_from_gershgorin: symmetric about 0, encloses the spectrum, padding 10 %;
  * jackson_damping_2d is the outer product of jackson_damping;
  * GreenRetarded / GreenRetardedDerivative (SymPy) against legacy_green_tables;
  * the Haldane curves (Greenwood, Bastin, sea, surface) rebuilt from the golden moments are
    byte-identical to the golden curves, and the graphene FFT-grid curves too;
  * Bastin = sea + surface on the Haldane golden, up to the difference of the two
    integration rules (< 0.3 % of the maximum);
  * Kubo-Greenwood sigma_xx >= 0 on the Haldane golden (it is a sum of squares);
  * input validation of kubo_plan and of the reconstructions.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
GOLD = os.path.join(REPO, "goldens", "stage0")
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

lsq = pytest.importorskip("lsquant", reason="lsquant._core not built")
_core = lsq._core

G5XY = "NonEqOpVX-VYhaldaneKPM_M128x128_stateexact"
G5XX = "NonEqOpVX-VXhaldaneKPM_M128x128_stateexact"
G4XX = "NonEqOpVX-VXgraphene_goldenKPM_M64x64_statedefault"


def dense_to_core(h):
    n = h.shape[0]
    rows, cols, vals = [0], [], []
    for i in range(n):
        nz = np.nonzero(h[i])[0]
        cols.extend(int(j) for j in nz)
        vals.extend(complex(h[i, j]) for j in nz)
        rows.append(len(cols))
    return _core.SparseMatrix.from_csr(n, np.asarray(rows, dtype=np.int64),
                                       np.asarray(cols, dtype=np.int64),
                                       np.asarray(vals, dtype=np.complex128))


def random_hermitian(n, seed):
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return 0.5 * (a + a.conj().T)


def chebyshev_matrices(h, window, m_tot):
    """T_m(H~) for m < m_tot, with H~ = H * scale + shift, by diagonalisation."""
    eps, u = np.linalg.eigh(h)
    x = eps * window.scale + window.shift
    assert np.all(np.abs(x) < 1.0)
    return [u @ np.diag(np.cos(m * np.arccos(x))) @ u.conj().T for m in range(m_tot)]


# ------------------------------------------------------------------ rank-2 moments
@pytest.fixture(scope="module")
def dense_model():
    n, m_tot = 20, 10
    h = random_hermitian(n, seed=11)
    a = random_hermitian(n, seed=12)
    b = random_hermitian(n, seed=13)
    eps = np.linalg.eigvalsh(h)
    window = _core.SpectralWindow(eps.min() - 0.3, eps.max() + 0.3)
    t = chebyshev_matrices(h, window, m_tot)
    w = np.where(np.arange(m_tot) == 0, 1.0, 2.0)
    ref = np.array([[np.trace(a @ t[m] @ b @ t[k]) for k in range(m_tot)] for m in range(m_tot)])
    ref *= np.outer(w, w) / n
    return dict(n=n, m=m_tot, h=h, a=a, b=b, window=window, ref=ref)


def run_kubo(d, symmetrise):
    states = _core.LocalStates(d["n"], list(range(d["n"])))
    ops = [_core.SparseMatrixOperator(dense_to_core(d["a"])),
           _core.SparseMatrixOperator(dense_to_core(d["b"]))]
    plan = lsq.kubo_plan(d["m"], left_index=0, right_index=1, word="A delta B delta")
    return lsq.run_plan(plan, dense_to_core(d["h"]), d["window"], states, operators=ops,
                        symmetrise=symmetrise)


def test_kubo_plan_against_dense_trace(dense_model):
    mu = run_kubo(dense_model, symmetrise=False)
    data = np.asarray(mu.data).reshape(dense_model["m"], dense_model["m"])
    assert mu.rank == 2
    scale = np.max(np.abs(dense_model["ref"]))
    assert np.max(np.abs(data - dense_model["ref"])) < 1e-12 * scale


def test_symmetrised_tensor_is_hermitian(dense_model):
    mu = run_kubo(dense_model, symmetrise=True)
    data = np.asarray(mu.data).reshape(dense_model["m"], dense_model["m"])
    assert np.max(np.abs(data - data.conj().T)) == 0.0
    sym = 0.5 * (dense_model["ref"] + dense_model["ref"].conj().T)
    assert np.max(np.abs(data - sym)) < 1e-12 * np.max(np.abs(sym))


def test_kubo_plan_rejects_bad_input():
    with pytest.raises(ValueError):
        lsq.kubo_plan(0, left_index=0, right_index=1, word="A delta B delta")


# ------------------------------------------------------------------ window, kernels, tables
def test_window_from_gershgorin_encloses_spectrum():
    h = random_hermitian(16, seed=5)
    w = lsq.window_from_gershgorin(dense_to_core(h))
    eps = np.linalg.eigvalsh(h)
    radius = np.max(np.sum(np.abs(h), axis=1))     # Gershgorin bound of a Hermitian matrix
    assert w.lower == pytest.approx(-w.upper)
    assert w.lower < eps.min() and w.upper > eps.max()
    assert w.upper == pytest.approx(1.10 * radius, rel=1e-12)


def test_jackson_2d_is_outer_product():
    g2 = lsq.jackson_damping_2d(7, 5)
    assert g2.shape == (7, 5)
    assert np.max(np.abs(g2 - np.outer(lsq.jackson_damping(7), lsq.jackson_damping(5)))) < 1e-15


def test_green_tables_sympy_matches_legacy():
    m_eff, n = 17, 200
    x = [-0.95 + i * 1.9 / (n - 1) for i in range(n)]
    leg_g, leg_dg = lsq.legacy_green_tables(x, m_eff)
    g = lsq.GreenRetarded.table(x, m_eff)
    dg = lsq.GreenRetardedDerivative.table(x, m_eff)
    assert np.max(np.abs(g - leg_g)) < 1e-12
    assert np.max(np.abs(dg - leg_dg)) < 1e-12 * max(1.0, np.max(np.abs(leg_dg)))


# ------------------------------------------------------------------ reconstructions
def read_dat(path):
    return np.loadtxt(path)


@pytest.mark.parametrize("func,moments,golden", [
    ("kubo_greenwood", G5XX, "KuboGreenWood_" + G5XX),
    ("kubo_bastin", G5XY, "KuboBastin_" + G5XY),
    ("kubo_bastin_sea", G5XY, "KuboBastinSea_" + G5XY),
    ("kubo_bastin_surf", G5XY, "KuboBastinSurf_" + G5XY),
])
def test_haldane_curves_identical_to_goldens(tmp_path, func, moments, golden):
    mu = _core.ReadChebmom(os.path.join(GOLD, "haldane", moments + ".chebmom2D"))
    e, v = getattr(lsq, func)(mu, 10.0)
    out = tmp_path / "curve.dat"
    _core.WriteCurveText(str(out), e, v)
    assert out.read_bytes() == open(os.path.join(GOLD, "haldane", golden + "JACKSON.dat"), "rb").read()


@pytest.mark.parametrize("formula,golden", [
    ("greenwood", "KuboGreenWood_" + G4XX + "JACKSON_FFTgrid.dat"),
    ("bastin", "KuboBastin_" + G4XX + "JACKSON_FFTgrid.dat"),
])
def test_graphene_fft_grid_identical_to_goldens(tmp_path, formula, golden):
    mu = _core.ReadChebmom(os.path.join(GOLD, "graphene", G4XX + ".chebmom2D"))
    e, v = lsq.kubo_fft_grid(mu, 10.0, formula=formula)
    out = tmp_path / "curve.dat"
    _core.WriteCurveText(str(out), e, v)
    assert out.read_bytes() == open(os.path.join(GOLD, "graphene", "fft_grid", golden), "rb").read()


def test_bastin_equals_sea_plus_surface():
    mu = _core.ReadChebmom(os.path.join(GOLD, "haldane", G5XY + ".chebmom2D"))
    e1, total = lsq.kubo_bastin(mu, 10.0)
    e2, sea = lsq.kubo_bastin_sea(mu, 10.0)
    e3, surf = lsq.kubo_bastin_surf(mu, 10.0)
    assert np.array_equal(e1, e2) and np.array_equal(e1, e3)
    # kubo_bastin integrates with the forward rectangle rule, kubo_bastin_sea with the
    # trapezoid rule, both on 30 M points: they agree up to that quadrature error,
    # 0.23 % of the maximum on this golden (largest near |E| = 1 eV, 8e-8 at E = 0).
    d = np.abs(total - (sea + surf))
    assert np.max(d) < 3e-3 * np.max(np.abs(total))
    assert d[np.argmin(np.abs(e1))] < 1e-6 * np.max(np.abs(total))


def test_greenwood_is_non_negative():
    mu = _core.ReadChebmom(os.path.join(GOLD, "haldane", G5XX + ".chebmom2D"))
    _, sxx = lsq.kubo_greenwood(mu, 10.0)
    assert np.min(sxx) > -1e-9 * np.max(sxx)


def test_reconstructions_reject_bad_input():
    mu2 = _core.ReadChebmom(os.path.join(GOLD, "haldane", G5XY + ".chebmom2D"))
    mu1 = _core.ReadChebmom(os.path.join(GOLD, "chain1d", "SpectralOp1chain1dKPM_M64_stateexact.chebmom1D"))
    for f in (lsq.kubo_greenwood, lsq.kubo_bastin, lsq.kubo_bastin_sea, lsq.kubo_bastin_surf):
        with pytest.raises(ValueError):
            f(mu2, 0.0)
        with pytest.raises(ValueError):
            f(mu1, 10.0)
