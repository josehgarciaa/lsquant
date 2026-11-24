"""Unit tests of the stage-4 Python layer: time evolution with Chebyshev series.

They need the compiled module ``lsquant._core`` (see test_python_stage2.py); without it
every test is skipped. Every check compares with dense matrix exponentials,
U(t) = exp(-i H t / hbar) with hbar = 0.6582119624 eV fs, on a random 12 x 12 model:

  * time_evolution_series: phase * sum_m c_m T_m(x) = exp(-i E dt / hbar) at every
    eigenvalue, for a centred and an off-centre window;
  * projected_evolution_plan: mu[m, n] = w_m Tr[P U(t_n)^+ O U(t_n) T_m(H~) P];
  * time_correlation_plan:    mu[m, n] = w_m Tr[A T_m(H~) U(t_n)^+ B U(t_n)];
  * msd_plan with shift-free V/a agrees with the dense MSD for both windows;
  * time_step_from_tmax, the corrected energy axis E = a x + c of the uniform-grid Kubo
    reconstructions, and time_correlation_curve against its definition.
w_m = (2 - delta_m0) / N, t_n = n dt, dt = tmax / (steps - 1).
"""

from __future__ import annotations

import math
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
sl = pytest.importorskip("scipy.linalg")
from lsquant.moments import HBAR, velocity_window  # noqa: E402

N, M, STEPS, TMAX = 12, 8, 5, 2.0


def herm(n, seed):
    r = np.random.default_rng(seed)
    a = r.normal(size=(n, n)) + 1j * r.normal(size=(n, n))
    return 0.5 * (a + a.conj().T)


def core(m):
    n = m.shape[0]
    rows, cols, vals = [0], [], []
    for i in range(n):
        nz = np.nonzero(np.abs(m[i]) > 0)[0]
        cols.extend(int(j) for j in nz)
        vals.extend(complex(m[i, j]) for j in nz)
        rows.append(len(cols))
    return _core.SparseMatrix.from_csr(n, np.asarray(rows, dtype=np.int64),
                                       np.asarray(cols, dtype=np.int64),
                                       np.asarray(vals, dtype=np.complex128))


@pytest.fixture(scope="module", params=[(-3.0, 3.0), (-8.0, 6.0)], ids=["centred", "offcentre"])
def model(request):
    lo, hi = request.param
    h = 0.3 * herm(N, 3)                       # spectrum inside [-2, 2.2]
    eps, u = np.linalg.eigh(h)
    w = _core.SpectralWindow(lo, hi)
    x = eps * w.scale + w.shift
    assert np.all(np.abs(x) < 1.0)
    t = [u @ np.diag(np.cos(m * np.arccos(x))) @ u.conj().T for m in range(M)]
    dt = lsq.time_step_from_tmax(STEPS, TMAX)
    coeffs, phase = lsq.time_evolution_series(w, dt)
    return dict(h=h, eps=eps, x=x, w=w, t=t, dt=dt, coeffs=coeffs, phase=phase,
                centred=(lo + hi == 0.0),
                wts=np.where(np.arange(M) == 0, 1.0, 2.0) / N)


def ut(d, time):
    return sl.expm(-1j * d["h"] * time / HBAR)


def states():
    return _core.LocalStates(N, list(range(N)))


def as_array(mu):
    return np.asarray(mu.data).reshape(mu.dims)


def test_time_step_from_tmax():
    assert lsq.time_step_from_tmax(5, 2.0) == 0.5


def test_bessel_series_reproduces_propagator(model):
    d = model
    series = np.array([sum(d["coeffs"][m] * math.cos(m * math.acos(xx)) for m in range(len(d["coeffs"])))
                       for xx in d["x"]])
    assert np.max(np.abs(d["phase"] * series - np.exp(-1j * d["eps"] * d["dt"] / HBAR))) < 1e-13


def test_projected_evolution_against_dense(model):
    d = model
    o, p = herm(N, 4), np.diag((np.arange(N) % 2).astype(float))
    plan = lsq.projected_evolution_plan(M, STEPS, N, d["coeffs"], operator_index=0, projector_index=1,
                                        phase=d["phase"])
    mu = lsq.run_plan(plan, core(d["h"]), d["w"], states(),
                      operators=[_core.SparseMatrixOperator(core(o)), _core.SparseMatrixOperator(core(p))])
    ref = np.array([[d["wts"][m] * np.trace(p @ ut(d, k * d["dt"]).conj().T @ o @ ut(d, k * d["dt"]) @ d["t"][m] @ p)
                     for k in range(STEPS)] for m in range(M)])
    assert list(mu.dims) == [M, STEPS]
    assert np.max(np.abs(as_array(mu) - ref)) < 1e-12 * np.max(np.abs(ref))


def test_time_correlation_against_dense(model):
    d = model
    a, b = herm(N, 5), herm(N, 6)
    plan = lsq.time_correlation_plan(M, STEPS, N, d["coeffs"], left_index=0, right_index=1, phase=d["phase"])
    mu = lsq.run_plan(plan, core(d["h"]), d["w"], states(),
                      operators=[_core.SparseMatrixOperator(core(a)), _core.SparseMatrixOperator(core(b))])
    ref = np.array([[d["wts"][m] * np.trace(a @ d["t"][m] @ ut(d, k * d["dt"]).conj().T @ b @ ut(d, k * d["dt"]))
                     for k in range(STEPS)] for m in range(M)])
    assert np.max(np.abs(as_array(mu) - ref)) < 1e-12 * np.max(np.abs(ref))


def test_msd_against_dense(model):
    d = model
    x_op = np.diag(np.random.default_rng(1).normal(size=N))
    v = 1j * (d["h"] @ x_op - x_op @ d["h"])
    plan = lsq.msd_plan(M, STEPS, N, d["coeffs"], velocity_index=0, phase=d["phase"])
    mu = lsq.run_plan(plan, core(d["h"]), d["w"], states(), operators=[_core.RescaledOperator(core(v), velocity_window(d["w"]))])

    def comm(time):
        return x_op @ ut(d, time) - ut(d, time) @ x_op

    ref = np.array([[d["wts"][m] * np.trace(comm(k * d["dt"]).conj().T @ d["t"][m] @ comm(k * d["dt"]))
                     for k in range(STEPS)] for m in range(M)])
    err = np.max(np.abs(as_array(mu) - ref)) / np.max(np.abs(ref))
    assert err < 1e-10
    if not d["centred"]:
        shifted = lsq.run_plan(
            plan, core(d["h"]), d["w"], states(),
            operators=[_core.RescaledOperator(core(v), d["w"])],
        )
        # A Hamiltonian shift is not part of the commutator velocity.
        assert np.max(np.abs(as_array(shifted) - ref)) > 0.1 * np.max(np.abs(ref))


def test_kubo_uniform_axis_is_exact_offcentre():
    g4c = os.path.join(GOLD, "graphene", "offcentre_enclosing")
    mu = _core.ReadChebmom(os.path.join(g4c, "NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexact.chebmom2D"))
    e, _ = lsq.kubo_bastin(mu, 10.0)
    gold = np.loadtxt(os.path.join(g4c, "KuboBastin_NonEqOpVX-VYgraphene_goldenKPM_M128x128_stateexactJACKSON.dat"))
    assert np.max(np.abs(e - gold[:, 0])) < 1e-5            # golden energies printed to 6 digits
    lo, hi = mu.window.lower, mu.window.upper
    assert lo + hi != 0.0


def test_time_correlation_curve_definition():
    g6 = os.path.join(GOLD, "chain1d_ballistic", "CorrelationVXchain1dKPM_M256_state20exact.chebmomTD")
    mu = _core.ReadChebmom(g6)
    t, val = lsq.time_correlation_curve(mu, 10.0, 0.5)
    m_tot, n_t = mu.dims
    a = 0.5 * (mu.window.upper - mu.window.lower)
    c = 0.5 * (mu.window.upper + mu.window.lower)
    eta = 2.0 * 10.0 / 1000.0 / (2.0 * a)
    m_eff = min(int(math.ceil(math.pi / eta)), m_tot)
    g = lsq.jackson_damping(m_eff)
    xf = (0.5 - c) / a
    delta = np.array([math.cos(m * math.acos(xf)) / (math.pi * math.sqrt(1 - xf * xf)) for m in range(m_eff)])
    data = np.asarray(mu.data).reshape(m_tot, n_t)
    ref = np.array([np.sum(delta * g * data[:m_eff, n].real) / a for n in range(n_t)])
    assert len(val) == n_t
    assert np.max(np.abs(val - ref)) < 1e-10 * max(1.0, np.max(np.abs(ref)))


@pytest.mark.parametrize("radius", [2.0, 2.2, 3.7])
def test_velocity_window_preserves_centred_bits(radius):
    import struct

    original = _core.SpectralWindow(-radius, radius)
    symmetric = velocity_window(original)
    for field in ("scale", "shift"):
        assert struct.pack("d", getattr(symmetric, field)) == struct.pack("d", getattr(original, field))
    assert math.copysign(1.0, symmetric.shift) == -1.0


@pytest.mark.parametrize("bounds", [(-3.0, 3.0), (-4.0, 2.5)], ids=["centred", "offcentre"])
def test_msd_runfile_against_dense(tmp_path, bounds):
    from scipy.sparse import csr_matrix
    from lsquant.models.nonorthogonal import write_csr_text
    from lsquant.runfile import execute_runfile

    n, moments, steps, dt = 4, 8, 4, 0.2
    h = np.diag(np.full(n - 1, -1.0), 1) + np.diag(np.full(n - 1, -1.0), -1)
    x = np.diag(np.arange(n, dtype=float))
    v = 1j * (h @ x - x @ h)
    ops = tmp_path / "ops"
    write_csr_text(ops / "chain.HAM.CSR", csr_matrix(h))
    write_csr_text(ops / "chain.VX.CSR", csr_matrix(v))
    raw = {
        "model": {"kind": "csr", "dir": str(ops), "label": "chain", "operators": ["VX"]},
        "basis": {"name": "real"},
        "metric": {"policy": "orthogonal"},
        "word": "Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]",
        "expansion": {"moments": moments, "window": list(bounds), "trace": "exact", "states": n},
        "time": {"steps": steps, "dt": dt},
        "reconstruction": {"broadening_meV": None},
        "output": {"dir": str(tmp_path / "out"), "legacy_layout": False},
    }
    result = execute_runfile(raw, repo_root=REPO)
    paths = [p for p in result["written_files"] if ".chebmom" in p]
    assert len(paths) == 1
    mu = _core.ReadChebmom(paths[0])
    eps, u = np.linalg.eigh(h)
    w = _core.SpectralWindow(*bounds)
    assert np.all(np.abs(eps * w.scale + w.shift) < 1.0)
    polys = [u @ np.diag(np.cos(m * np.arccos(eps * w.scale + w.shift))) @ u.conj().T
             for m in range(moments)]
    ref = np.empty((moments, steps), dtype=complex)
    for k in range(steps):
        propagator = sl.expm(-1j * h * k * dt / HBAR)
        comm = x @ propagator - propagator @ x
        for m in range(moments):
            ref[m, k] = (1.0 if m == 0 else 2.0) / n * np.trace(comm.conj().T @ polys[m] @ comm)
    assert np.max(np.abs(as_array(mu) - ref)) < 1e-10
