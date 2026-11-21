"""Unit tests of the stage-6 non-orthogonal metric: overlap appliers and the two routes.

They need the compiled module ``lsquant._core`` (see test_python_stage2.py); without it
every test is skipped. The tests build a small rotated model, so they run in seconds. The full-size
identities against the goldens G1 and G5 run in tests/dev/stage_6_check.py and in
examples/nonorthogonal/haldane_rotated.py.

What is checked:
  * the exact appliers: S (S^-1 v) = v and (S^-1/2)^2 S v = v;
  * Ischo (iterative) and Chebyshev (polynomial, bounds from spec(S)) S^-1 against the
    exact applier;
  * on a small random model rewritten with R = I + alpha K (N = 96): the Loewdin route
    (DOS, Kubo A-B) and the generalised route (DOS, Kubo A-B with the covariant plan
    bra: S^-1, S^-1 A, Expand, S) equal the orthogonal moments;
  * the naive generalised estimator (S^-1 A on both sides, no S on the bra) differs from
    the orthogonal moments by O(1);
  * OrthogonalMetric and NonOrthogonalMetric report what they are.
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



def op(m):
    return _core.SparseMatrixOperator(m)


def apply(o, v):
    return np.asarray(o.multiply(np.asarray(v, dtype=np.complex128)))


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


# ---------------------------------------------------------------- routes on a small rotated model
# A random banded Hermitian H and two Hermitian operators A, B (N = 96), rewritten in the
# non-orthogonal basis R = I + alpha K of lsquant.models.nonorthogonal. Every trace in the
# coefficient representation must equal the orthogonal one, computed here with the same plans.
N_SMALL, M_SMALL = 96, 24
WIN = (-6.0, 6.0)


@pytest.fixture(scope="module")
def small():
    sp = pytest.importorskip("scipy.sparse")
    from lsquant.models.nonorthogonal import make_local_transform, rotate_basis
    rng = np.random.default_rng(7)

    def banded(width, scale):
        m = np.zeros((N_SMALL, N_SMALL), dtype=complex)
        for k in range(1, width + 1):
            d = scale * (rng.normal(size=N_SMALL - k) + 1j * rng.normal(size=N_SMALL - k))
            m += np.diag(d, k)
        m = m + m.conj().T + np.diag(scale * rng.normal(size=N_SMALL))
        return sp.csr_matrix(m)

    orth = {"HAM": banded(3, 0.6), "A": banded(2, 1.0), "B": banded(2, 1.0)}
    rot = rotate_basis(orth, make_local_transform(orth["HAM"], alpha=0.55, complex_entries=True, seed=12345))

    def core(m):
        m = sp.csr_matrix(m)
        m.sort_indices()
        return _core.SparseMatrix.from_csr(N_SMALL, m.indptr.astype(np.int64), m.indices.astype(np.int64),
                                           m.data.astype(np.complex128))

    eps = np.linalg.eigvalsh(orth["HAM"].toarray())
    assert WIN[0] < eps.min() and eps.max() < WIN[1]
    s_eig = np.linalg.eigvalsh(rot["S"].toarray())
    return {"orth": {k: core(v) for k, v in orth.items()}, "rot": {k: core(v) for k, v in rot.items()},
            "s_eig": (float(s_eig.min()), float(s_eig.max()))}


def test_exact_appliers(small):
    s = small["rot"]["S"]
    inv, isq = _core.ExactInverseOverlap(s), _core.ExactInverseSqrtOverlap(s)
    v = np.random.default_rng(1).normal(size=N_SMALL) + 0j
    assert rel(apply(op(s), apply(inv, v)), v) < 1e-12
    assert rel(apply(isq, apply(isq, apply(op(s), v))), v) < 1e-12


def test_approximate_appliers(small):
    s = small["rot"]["S"]
    lo, hi = small["s_eig"]
    exact = _core.ExactInverseOverlap(s)
    ischo = _core.IschoInverseOverlap(s, 1e-13, 200)
    cheb = _core.ChebyshevInverseOverlap(s, _core.SpectralBounds(0.9 * lo, 1.1 * hi), 60)
    rng = np.random.default_rng(2)
    v = rng.normal(size=N_SMALL) + 1j * rng.normal(size=N_SMALL)
    ref = apply(exact, v)
    assert rel(apply(ischo, v), ref) < 1e-10
    # PCG with the incomplete-Cholesky preconditioner; on this banded S the preconditioner
    # is already exact, so last_iterations may be 0. The reported residual meets the tolerance.
    assert ischo.last_error <= 1e-13
    assert rel(apply(cheb, v), ref) < 1e-10
    assert cheb.effective_order <= 60


def exact_states():
    return _core.LocalStates(N_SMALL, list(range(N_SMALL)))


def orth_dos(m):
    return np.asarray(lsq.run_plan(lsq.dos_plan(M_SMALL, N_SMALL, word="1 delta"), m["HAM"],
                                   _core.SpectralWindow(*WIN), exact_states()).data)


def orth_kubo(m):
    mu = lsq.run_plan(lsq.kubo_plan(M_SMALL, left_index=0, right_index=1, word="A delta B delta"), m["HAM"],
                      _core.SpectralWindow(*WIN), exact_states(), operators=[op(m["A"]), op(m["B"])])
    return np.asarray(mu.data)


def test_lowdin_dos_equals_orthogonal(small):
    r = small["rot"]
    isq = _core.ExactInverseSqrtOverlap(r["S"])
    mu = lsq.run_plan(lsq.dos_plan(M_SMALL, N_SMALL, word="1 delta"), _core.LowdinOperator(isq, op(r["HAM"])),
                      _core.SpectralWindow(*WIN), exact_states())
    assert np.max(np.abs(np.asarray(mu.data) - orth_dos(small["orth"]))) < 1e-12


def test_generalised_dos_equals_orthogonal(small):
    r = small["rot"]
    rh = _core.GeneralisedRescaledOperator(r["HAM"], r["S"], _core.ExactInverseOverlap(r["S"]),
                                           _core.SpectralWindow(*WIN))
    mu = _core.MomentEngine(rh, [], exact_states(), M_SMALL).run(lsq.dos_plan(M_SMALL, N_SMALL, word="1 delta"))
    assert np.max(np.abs(np.asarray(mu.data) - orth_dos(small["orth"]))) < 1e-12


def test_lowdin_kubo_equals_orthogonal(small):
    r = small["rot"]
    isq = _core.ExactInverseSqrtOverlap(r["S"])
    lo = lambda m: _core.LowdinOperator(isq, op(m))
    mu = lsq.run_plan(lsq.kubo_plan(M_SMALL, left_index=0, right_index=1, word="A delta B delta"), lo(r["HAM"]),
                      _core.SpectralWindow(*WIN), exact_states(), operators=[lo(r["A"]), lo(r["B"])])
    ref = orth_kubo(small["orth"])
    assert np.max(np.abs(np.asarray(mu.data) - ref)) < 1e-11 * np.max(np.abs(ref))


def generalised_kubo(r, covariant):
    s = r["S"]
    inv = _core.ExactInverseOverlap(s)
    rh = _core.GeneralisedRescaledOperator(r["HAM"], s, inv, _core.SpectralWindow(*WIN))
    # operator table: 0 S^-1 A, 1 S^-1 B, 2 S^-1, 3 S
    table = [_core.MetricInverseProduct(inv, op(r["A"])), _core.MetricInverseProduct(inv, op(r["B"])), inv, op(s)]
    if covariant:
        # ket: Expand, S^-1 B; bra: S^-1, S^-1 A, Expand, S  -> Tr[S^-1 A T_m S^-1 B T_n]
        code = [_core.Expand(M_SMALL), _core.ApplyOperator(1), _core.MirrorToBra(),
                _core.ApplyOperator(2), _core.ApplyOperator(0), _core.Expand(M_SMALL), _core.ApplyOperator(3),
                _core.Trace([1.0 + 0j])]
    else:
        code = [_core.Expand(M_SMALL), _core.ApplyOperator(1), _core.MirrorToBra(),
                _core.ApplyOperator(0), _core.Expand(M_SMALL), _core.Trace([1.0 + 0j])]
    mu = _core.MomentEngine(rh, table, exact_states(), M_SMALL).run(_core.Plan(code, "A delta B delta"))
    _core.ApplyKpmWeights(mu, N_SMALL, True)
    return np.asarray(mu.data)


def test_generalised_covariant_kubo_equals_orthogonal(small):
    ref = orth_kubo(small["orth"])
    assert np.max(np.abs(generalised_kubo(small["rot"], True) - ref)) < 1e-11 * np.max(np.abs(ref))


def test_naive_generalised_kubo_is_not_covariant(small):
    ref = orth_kubo(small["orth"])
    assert np.max(np.abs(generalised_kubo(small["rot"], False) - ref)) > 1e-2 * np.max(np.abs(ref))


def test_metric_policies(small):
    assert _core.OrthogonalMetric(4).is_orthogonal
    s = small["rot"]["S"]
    m = _core.NonOrthogonalMetric(s, _core.ExactInverseOverlap(s), _core.ExactInverseSqrtOverlap(s))
    assert not m.is_orthogonal
