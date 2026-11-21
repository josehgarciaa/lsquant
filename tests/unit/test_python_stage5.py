"""Unit tests of the stage-5 bases: real space, Bloch (k space) and per-k eigenbasis.

They need the compiled module ``lsquant._core`` (see test_python_stage2.py); without it
every test is skipped. The model is the graphene Wannier model of golden G9
(goldens/stage0/bases/graphene/model), on a 6 x 4 x 1 supercell (48 orbitals). No
wannier2sparse binary is needed: RealBasis builds the supercell CSR itself.

What is checked:
  * BlochTransform is unitary (B B^+ x = B^+ B x = x);
  * exact-trace DOS and VX-VX moments are identical in the real, Bloch and eigen bases;
  * random-phase moments are identical in the real and Bloch bases (same physical vectors);
  * with Anderson disorder (W = 1 eV, seed 7) the real and Bloch moments agree;
  * EigenBasis.exact_dos integrates to the number of orbitals per cell (2);
  * EigenBasis rejects disorder, make_basis rejects an unknown basis name.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
MODEL = os.path.join(REPO, "goldens", "stage0", "bases", "graphene", "model")
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

lsq = pytest.importorskip("lsquant", reason="lsquant._core not built")
_core = lsq._core
from lsquant.bases import make_basis  # noqa: E402
from lsquant.models.wannier import load_wannier_model  # noqa: E402

SC = (6, 4, 1)
M = 32
WINDOW = (-3.6, 3.6)


@pytest.fixture(scope="module")
def model():
    return load_wannier_model(MODEL, seed="graphene")


def dos(basis, states):
    plan = lsq.dos_plan(M, int(states.number_of_states), word="1 delta")
    mu = lsq.run_plan(plan, basis.hamiltonian(), _core.SpectralWindow(*WINDOW), states, symmetrise=False)
    return np.asarray(mu.data)


def kubo_xx(basis, states):
    plan = lsq.kubo_plan(M, left_index=0, right_index=0, word="VX delta VX delta")
    mu = lsq.run_plan(plan, basis.hamiltonian(), _core.SpectralWindow(*WINDOW), states,
                      operators=[basis.operator("VX")])
    return np.asarray(mu.data)


def test_bloch_transform_is_unitary(model):
    b = make_basis("bloch", model, SC)
    fwd = b.transform
    inv = fwd.adjoint()
    rng = np.random.default_rng(3)
    x = rng.normal(size=b.dimension) + 1j * rng.normal(size=b.dimension)
    y = np.asarray(inv.multiply(np.asarray(fwd.multiply(x))))
    z = np.asarray(fwd.multiply(np.asarray(inv.multiply(x))))
    assert np.max(np.abs(y - x)) < 1e-13 * np.max(np.abs(x))
    assert np.max(np.abs(z - x)) < 1e-13 * np.max(np.abs(x))


@pytest.mark.parametrize("name", ["bloch", "eigen"])
def test_exact_moments_identical_across_bases(model, name):
    real = make_basis("real", model, SC)
    other = make_basis(name, model, SC)
    d_real, d_other = dos(real, real.states("exact")), dos(other, other.states("exact"))
    k_real, k_other = kubo_xx(real, real.states("exact")), kubo_xx(other, other.states("exact"))
    assert np.max(np.abs(d_real - d_other)) < 1e-12
    assert np.max(np.abs(k_real - k_other)) < 1e-12 * np.max(np.abs(k_real))
    assert d_real[0] == pytest.approx(1.0, abs=1e-14)


def test_random_phase_moments_identical_real_bloch(model):
    real = make_basis("real", model, SC)
    bloch = make_basis("bloch", model, SC)
    d1 = dos(real, real.states("random", num_states=4, seed=12345))
    d2 = dos(bloch, bloch.states("random", num_states=4, seed=12345))
    assert np.max(np.abs(d1 - d2)) < 1e-12


def test_disorder_real_and_bloch_agree(model):
    real = make_basis("real", model, SC, disorder_w=1.0, disorder_seed=7)
    bloch = make_basis("bloch", model, SC, disorder_w=1.0, disorder_seed=7)
    d1 = dos(real, real.states("random", num_states=4, seed=12345))
    d2 = dos(bloch, bloch.states("random", num_states=4, seed=12345))
    clean = dos(make_basis("real", model, SC), real.states("random", num_states=4, seed=12345))
    assert np.max(np.abs(d1 - d2)) < 1e-12
    assert np.max(np.abs(d1 - clean)) > 1e-3          # the disorder changes the moments


def test_exact_dos_sum_rule(model):
    e = np.linspace(-3.6, 3.6, 2881)
    d = make_basis("eigen", model, (24, 24, 1)).exact_dos(e, 0.04)
    assert np.sum(d) * (e[1] - e[0]) == pytest.approx(2.0, abs=1e-9)


def test_eigenbasis_rejects_disorder(model):
    with pytest.raises(ValueError):
        make_basis("eigen", model, SC, disorder_w=1.0)


def test_make_basis_rejects_unknown_name(model):
    with pytest.raises(ValueError):
        make_basis("plane-wave", model, SC)
