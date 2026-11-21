"""Unit tests of the materials layer (stages 6b-8): SIESTA reader, metric forms, disorder.

They need the compiled module ``lsquant._core`` (see test_python_stage2.py); the SIESTA
tests also need ``sisl``. Without them the tests are skipped. The heavy identities
(spin Hall against the exact-diagonalisation oracle, disorder in the Bloch and real
bases, projection) run in tests/dev/stage_6b_check.py, stage_6c_check.py and
stage_6d_check.py.

What is checked:
  * metric_forms: the symbolic proof that the Loewdin and generalised routes give the same
    trace for every catalogue word shape, and that the covariant velocity is Hermitian;
  * the SIESTA reader on the committed Bi bilayer inputs (examples/07): spinor blocks,
    S(R = 0) Hermitian (to the 1e-7 relative TSHS rounding) and positive definite, H(k) and S(k) Hermitian at a generic k,
    and the generalised eigenvalues real and inside the eigen window;
  * anderson_potential: V = (w_i + w_j) S_ij / 2, equal to w S for a uniform w, Hermitian,
    and refused for complex or mis-shaped values.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
BI = os.path.join(REPO, "examples", "07_spin_hall_bi_bilayer", "siesta")
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")
sparse = pytest.importorskip("scipy.sparse")


# ------------------------------------------------------------------ symbolic forms
def test_routes_have_equal_traces_for_every_word_shape():
    from lsquant import metric_forms
    result = metric_forms.trace_equivalence()
    assert result["hamiltonian_similarity"]
    shapes = [k for k, v in result.items() if isinstance(v, dict) and "equal" in v]
    assert {"dos", "greenwood", "bastin", "time_correlation"} <= set(shapes)
    assert all(result[k]["equal"] for k in shapes)


def test_covariant_velocity_is_hermitian():
    from lsquant import metric_forms
    verdict = metric_forms.hermiticity(metric_forms.velocity_forms("covariant"))
    # physical V, the Loewdin S^-1/2 V S^-1/2 and the anticommutator are Hermitian; the
    # generalised S^-1 V is S-Hermitian (S^-1 V)^+ S = S (S^-1 V)
    for name in ("physical", "lowdin", "anticommutator"):
        assert verdict[name]["hermitian"], name
    assert verdict["generalised"]["s_hermitian"]


# ------------------------------------------------------------------ SIESTA reader
@pytest.fixture(scope="module")
def bi_model():
    pytest.importorskip("sisl")
    from lsquant.models.siesta import load_siesta_model
    return load_siesta_model(BI, "bi")


def bloch(model, matrices, k):
    phases = np.exp(2j * np.pi * (model.r_int @ k))
    return np.tensordot(phases, matrices, axes=1)


def test_bi_model_blocks(bi_model):
    m = bi_model
    assert m.num_wann % 2 == 0                       # spinor blocks: W = 2 N_orb
    r0 = int(np.flatnonzero(np.all(m.r_int == 0, axis=1))[0])
    s0 = m.overlaps[r0]
    # TSHS stores rounded values: Hermitian to the 1e-7 relative bound the reader documents
    assert np.max(np.abs(s0 - s0.conj().T)) < 1e-7 * np.max(np.abs(s0))
    assert np.linalg.eigvalsh(0.5 * (s0 + s0.conj().T)).min() > 0.0


def test_bi_bloch_matrices_hermitian_and_spectrum_real(bi_model):
    m = bi_model
    k = np.array([0.137, 0.291, 0.0])
    hk, sk = bloch(m, m.hoppings, k), bloch(m, m.overlaps, k)
    assert np.max(np.abs(hk - hk.conj().T)) < 1e-7 * np.max(np.abs(hk))
    assert np.max(np.abs(sk - sk.conj().T)) < 1e-7 * np.max(np.abs(sk))
    from scipy.linalg import eigh
    e = eigh(hk, sk, eigvals_only=True)
    assert np.all(np.isfinite(e))
    assert e.min() > -40.0 and e.max() < 70.0          # inside the SIESTA eigen window [-32.8, 64.5]


# ------------------------------------------------------------------ disorder potential
def test_anderson_potential_formula():
    from lsquant.models.disorder import anderson_potential
    rng = np.random.default_rng(4)
    a = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
    s = a @ a.conj().T + 8 * np.eye(8)                  # Hermitian positive definite overlap
    w = rng.uniform(-0.5, 0.5, size=8)
    v = anderson_potential(sparse.csr_matrix(s), w).toarray()
    assert np.max(np.abs(v - 0.5 * (np.diag(w) @ s + s @ np.diag(w)))) < 1e-12
    assert np.max(np.abs(v - v.conj().T)) < 1e-12
    vu = anderson_potential(sparse.csr_matrix(s), np.full(8, 0.3)).toarray()
    assert np.max(np.abs(vu - 0.3 * s)) < 1e-12


def test_anderson_potential_refuses_bad_values():
    from lsquant.models.disorder import anderson_potential
    s = sparse.identity(4, dtype=complex, format="csr")
    with pytest.raises(ValueError):
        anderson_potential(s, np.ones(3))
    with pytest.raises(ValueError):
        anderson_potential(s, np.ones(4) * (1 + 1j))
