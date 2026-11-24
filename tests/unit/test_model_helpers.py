"""Shared helpers of the model readers: units, spinor lift, Hermitian pairing.

The readers used to carry their own copy of the Bohr/Hartree/Rydberg
constants, of the spinor lift ``S tensor I_2`` and of the opposite-image
Hermitian fold. This module pins the one implementation of each: the
constants as literals (the Rydberg values are the codes' own, kept distinct),
the lift of a single matrix and of a stack against ``np.kron``, and the
pairing (a missing ``-R`` added as ``M(R)^+``, a present pair folded to its
mean) with each residual normalisation the readers write to ``basis.json``.

Pure NumPy checks: no model files, no optional dependency (sisl, h5py, ...).
"""
import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")

from lsquant.models import units  # noqa: E402
from lsquant.models.lcao import (  # noqa: E402
    fold_opposite,
    hermitian_pairs,
    opposite_pairs,
    seed_name_from_file,
    spinor_overlap,
)
from lsquant.models.siesta import spinor_lift, spinor_scalar  # noqa: E402


def test_unit_constants_pinned():
    # The Rydberg values are the codes' own conversion constants and differ
    # from each other and from HARTREE_EV / 2 on purpose.
    assert units.BOHR_ANGSTROM == 0.529177210903
    assert units.HARTREE_EV == 27.211386245988
    assert units.RYDBERG_EV_ABACUS == 13.605698066
    assert units.RYDBERG_EV_THREEBODYTB == 13.6056980659
    assert units.ANGSTROM_BOHR == 1.8897259886
    # The readers keep their historical names, now aliases of units.py.
    from lsquant.models import abacus, hamgnn, openmx, pyscf_lcao, threebodytb
    assert abacus.RYDBERG_EV is units.RYDBERG_EV_ABACUS
    assert abacus.ANGSTROM_BOHR is units.ANGSTROM_BOHR
    assert threebodytb.RYDBERG_EV is units.RYDBERG_EV_THREEBODYTB
    assert threebodytb.BOHR_ANGSTROM is units.BOHR_ANGSTROM
    assert hamgnn.BOHR_ANGSTROM is units.BOHR_ANGSTROM
    assert hamgnn.HARTREE_EV is units.HARTREE_EV
    assert openmx.BOHR is units.BOHR_ANGSTROM
    assert openmx.HARTREE is units.HARTREE_EV
    assert pyscf_lcao.BOHR is units.BOHR_ANGSTROM
    assert pyscf_lcao.HARTREE is units.HARTREE_EV


def test_spinor_lift_matrix_and_stack():
    rng = np.random.default_rng(7)
    a = rng.normal(size=(3, 3))
    b = a + 1j * rng.normal(size=(3, 3))
    for x in (a, b):
        lift = spinor_lift(x)
        reference = np.kron(x, np.eye(2))
        np.testing.assert_array_equal(lift, reference)
        assert lift.shape == (6, 6)
    # The dtype follows the input (the model constructors convert to complex).
    assert spinor_lift(a).dtype == np.float64
    assert spinor_lift(b).dtype == np.complex128
    # Leading axes are kept: an (N_R, n, n) stack lifts block by block.
    stack = rng.normal(size=(4, 3, 3)) + 1j * rng.normal(size=(4, 3, 3))
    lifted = spinor_lift(stack)
    assert lifted.shape == (4, 6, 6)
    for n in range(4):
        np.testing.assert_array_equal(lifted[n], np.kron(stack[n], np.eye(2)))
    # The historical names are aliases of the one function.
    assert spinor_scalar is spinor_lift
    assert spinor_overlap is spinor_lift


def test_opposite_pairs():
    r_int = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 2, -3]])
    partner = opposite_pairs(r_int)
    for n, p in enumerate(partner):
        np.testing.assert_array_equal(-r_int[n], r_int[p])
    with pytest.raises(ValueError, match="not closed under R -> -R"):
        opposite_pairs(np.array([[0, 0, 0], [1, 0, 0]]))


def _random_stack(rng, r_int, n):
    mats = np.zeros((len(r_int), n, n), complex)
    for i, r in enumerate(r_int):
        if tuple(r) == (0, 0, 0):
            h = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            mats[i] = h + h.conj().T
        else:
            mats[i] = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return mats


def test_hermitian_pairs_adds_missing_opposite():
    rng = np.random.default_rng(11)
    r_int = np.array([[0, 0, 0], [1, 0, 0]])
    mats = _random_stack(rng, r_int, 3)
    r_out, folded, residual = hermitian_pairs(r_int, mats)
    assert [tuple(r) for r in r_out] == [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]
    # The added -R block is M(R)^+ and every pair is exactly Hermitian.
    np.testing.assert_array_equal(folded[2], folded[1].conj().T)
    # The fold is the mean with the mirrored stack, and the residual is the
    # pre-fold mismatch over (2 max|M|), rebuilt here from the same inputs.
    extended = np.concatenate([mats, [mats[1].conj().T]])
    mirrored = extended[[0, 2, 1]]
    np.testing.assert_array_equal(folded, 0.5 * (extended + mirrored))
    scale = max(float(np.abs(extended).max()), 1e-300)
    assert residual == pytest.approx(float(np.abs(extended - mirrored).max()) / (2 * scale))


def test_hermitian_pairs_present_pair_and_residual_modes():
    rng = np.random.default_rng(13)
    r_int = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
    mats = _random_stack(rng, r_int, 3)
    r_out, folded, residual = hermitian_pairs(r_int, mats)
    assert [tuple(r) for r in r_out] == [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]
    partner = opposite_pairs(r_out)
    mirrored = mats[partner].conj().swapaxes(1, 2)
    np.testing.assert_array_equal(folded, 0.5 * (mats + mirrored))
    scale = max(float(np.abs(mats).max()), 1e-300)
    # "mean2": max|M - M(-R)^+| / (2 max|M|) (the H(R)/S(R) readers).
    assert residual == pytest.approx(float(np.abs(mats - mirrored).max()) / (2 * scale))
    # "mean": max|M - fold| / max|M| (OpenMX, ABACUS): same value, other form.
    _, mean = fold_opposite(mats, partner, residual="mean")
    assert mean == pytest.approx(float(np.abs(mats - folded).max()) / scale)
    # "absolute": max|M - M(-R)^+| (SIESTA's certificate, GPAW's annotation).
    _, absolute = fold_opposite(mats, partner, residual="absolute")
    assert absolute == pytest.approx(float(np.abs(mats - mirrored).max()))


def test_fold_opposite_leading_axis_and_unknown_mode():
    # The connection stacks keep a leading axis: (3, N_R, n, n).
    rng = np.random.default_rng(17)
    r_int = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
    mats = rng.normal(size=(3, 3, 3, 3)) + 1j * rng.normal(size=(3, 3, 3, 3))
    partner = opposite_pairs(r_int)
    fold, residual = fold_opposite(mats, partner, residual="mean")
    assert fold.shape == mats.shape
    for v in range(3):
        np.testing.assert_array_equal(fold[v], 0.5 * (mats[v] + mats[v, partner].conj().swapaxes(1, 2)))
    assert residual == pytest.approx(
        float(np.abs(mats - fold).max()) / max(float(np.abs(mats).max()), 1e-300))
    with pytest.raises(ValueError, match="unknown residual normalisation"):
        fold_opposite(mats, partner, residual="bogus")


def test_seed_name_from_file():
    assert seed_name_from_file("/data/models/silicon.scfout") == "silicon"
    assert seed_name_from_file("hamiltonian.npz") == "hamiltonian"
    assert seed_name_from_file(os.path.join("dir", "with.dot", "name")) == "name"
