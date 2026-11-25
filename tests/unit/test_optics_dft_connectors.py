"""The DFT connectors used by the optical run files (example 24).

- ``hermitize_blocks`` makes M(-R) = M(R)^dagger exact, so the Bloch sum is
  Hermitian at every k, and refuses a block without its partner.
- ``write_paoflow_model`` writes a spinless model (no spin-orbit in the DFT
  run) with a zero spin operator instead of failing.
- The metric executor accepts the longitudinal words
  ``Tr[VX delta(E-H) VX dG-(E)/dE]`` and ``Tr[VY ...]`` together, writes one
  rank-2 moment file per word and the curves ``sigma_xx``/``sigma_yy``, and
  refuses ``oracle: true`` for a run with only longitudinal words.
"""
import json
import os
import sys
from types import SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from lsquant.models.orthogonal import load_paoflow_model, write_paoflow_model  # noqa: E402
from lsquant.models.pyscf_lcao import hermitize_blocks  # noqa: E402

XX = "Tr[VX delta(E-H) VX dG-(E)/dE]"
YY = "Tr[VY delta(E-H) VY dG-(E)/dE]"


def _bloch(r_int, blocks, k_frac):
    return np.einsum("kr,rij->kij", np.exp(2j * np.pi * np.asarray(k_frac) @ np.asarray(r_int).T), blocks)


def test_hermitize_blocks_gives_hermitian_bloch_sums():
    rng = np.random.default_rng(3)
    r_int = np.array([(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)])
    exact = {}
    for r in ((0, 0, 0), (1, 0, 0), (0, 1, 0)):
        a = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
        exact[r] = a + a.conj().T if r == (0, 0, 0) else a
        exact[tuple(-x for x in r)] = exact[r].conj().T
    clean = np.array([exact[tuple(r)] for r in r_int])
    noisy = clean + 1e-9 * (rng.normal(size=clean.shape) + 1j * rng.normal(size=clean.shape))
    fixed = hermitize_blocks(r_int, noisy)
    for i, r in enumerate(r_int):
        j = next(j for j, s in enumerate(r_int) if np.all(s == -r))
        np.testing.assert_array_equal(fixed[i], fixed[j].conj().T)
    hk = _bloch(r_int, fixed, rng.random((6, 3)) * [1, 1, 0])
    np.testing.assert_allclose(hk, hk.conj().swapaxes(-1, -2), atol=1e-15)
    assert np.abs(fixed - clean).max() < 1e-8                       # only the noise is touched
    with pytest.raises(ValueError, match="no partner"):
        hermitize_blocks(r_int[:2], noisy[:2])


def test_spinless_paoflow_model_has_zero_spin(tmp_path):
    n = 2
    hrs = np.zeros((n, n, 3, 3, 1, 1), complex)
    hrs[:, :, 0, 0, 0, 0] = np.diag([-1.0, 1.0])
    hrs[0, 1, 1, 0, 0, 0] = hrs[1, 0, 2, 0, 0, 0] = 0.3            # H(-R) = H(R)^+
    arrays = dict(HRs=hrs, a_vectors=np.diag([1.0, 1.0, 4.0]), Dnm=np.zeros((n, n, 3)))
    attrs = dict(alat=4.0, dftSO=False)
    calls = []
    fake = SimpleNamespace(data_controller=SimpleNamespace(data_dicts=lambda: (arrays, attrs)),
                           spin_operator=lambda **_: calls.append(1))
    path = str(tmp_path / "spinless.npz")
    write_paoflow_model(fake, path)
    assert not calls                                                  # no spin operator requested
    np.testing.assert_array_equal(np.load(path)["spin"], 0)
    model = load_paoflow_model(path)
    assert model.Hk(np.zeros((1, 3))).shape == (1, n, n)


def _two_band_insulator(tmp_path):
    """Gapped honeycomb-like toy: on-site +-2 eV, hoppings 1 eV along a1 and a2, spinless."""
    cells = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
    hop = np.array([[0.0, 1.0], [0.0, 0.0]])
    blocks = {(0, 0, 0): np.array([[2.0, 1.0], [1.0, -2.0]], complex)}
    for r in ((1, 0, 0), (0, 1, 0)):
        blocks[r] = hop.astype(complex)
        blocks[tuple(-x for x in r)] = hop.T.astype(complex)
    path = tmp_path / "toy.npz"
    np.savez(path, uc=np.array([[2.5, 0, 0], [-1.25, 2.165, 0], [0, 0, 20.0]]),
             tau=np.array([[0.0, 1.443, 0.0], [1.25, 0.722, 0.0]]), r_int=np.array(cells),
             hoppings=np.array([blocks[c] for c in cells]), spin=np.zeros((3, 2, 2)),
             fermi_energy=0.0, mesh=np.array([3, 3, 1]), edge_ratio=0.0)
    return str(path)


def _runfile(tmp_path, words, oracle=False):
    return {
        "model": {"kind": "paoflow", "file": _two_band_insulator(tmp_path),
                  "params": {"velocity_mode": "berry_connection", "oracle": oracle}},
        "basis": {"name": "eigen", "supercell": [6, 6, 1]},
        "metric": {"policy": "orthogonal"},
        "word": words[0], "words": words,
        "expansion": {"moments": 32, "window": [-12.0, 12.0], "trace": "exact", "seed": 1},
        "reconstruction": {"kernel": "jackson", "broadening_meV": 500.0},
        "output": {"legacy_layout": False, "write_chebmom": True, "dir": str(tmp_path / "out")},
    }


def test_longitudinal_words_run_on_an_orthogonal_dft_model(tmp_path):
    from lsquant.runfile import execute_runfile
    execute_runfile(_runfile(tmp_path, [XX, YY]))
    out = tmp_path / "out"
    names = sorted(os.listdir(out))
    assert any("VX-VX" in n and n.endswith(".chebmom") for n in names), names
    assert any("VY-VY" in n and n.endswith(".chebmom") for n in names), names
    curves = json.dumps(names) + "".join(open(out / n).read() for n in names if n.endswith(".json"))
    assert "sigma_xx" in curves and "sigma_yy" in curves


def test_oracle_is_refused_for_longitudinal_words_only(tmp_path):
    from lsquant.runfile import execute_runfile
    with pytest.raises(ValueError, match="longitudinal"):
        execute_runfile(_runfile(tmp_path, [XX, YY], oracle=True))


def test_paoflow_writer_keeps_the_hexagonal_symmetry_between_mesh_points(tmp_path):
    """A hexagonal model with hoppings on the Nyquist shell of a 4 x 4 grid.

    The writer must place every grid element on its Wigner-Seitz image: the
    model is then recovered exactly at every k, and E(k) = E(C3 k). Folding
    the grid indices to a fixed parallelogram puts part of the hoppings of
    length 2a on lattice vectors of length 2 sqrt(3) a.
    """
    bohr = 0.529177210903
    a = 2.5
    uc = np.array([[a, 0, 0], [-a / 2, a * np.sqrt(3) / 2, 0], [0, 0, 15.0]])
    tau = np.array([[0.0, a / np.sqrt(3), 0.0], [a / 2, a / (2 * np.sqrt(3)), 0.0]])
    hop = [(0, 0, (0, 0), 1.5), (1, 1, (0, 0), -1.5)]                   # (m, n, R', t)
    for r in ((0, 0), (-1, 0), (0, 1)):                                   # B-N nearest neighbours
        bond = np.array(r + (0,)) @ uc + tau[1] - tau[0]
        assert abs(np.linalg.norm(bond) - a / np.sqrt(3)) < 1e-9
        hop.append((0, 1, r, -2.3))
    for r in ((1, 0), (0, 1), (1, 1), (2, 0), (0, 2), (2, 2)):            # same sublattice, a and 2a
        t = 0.2 if max(r) == 1 else 0.05
        for m in (0, 1):
            hop.append((m, m, r, t))
    # Hermitian partners H_nm(-R') = H_mn(R')^*; the on-site energies are their own partners.
    hop += [(n, m, tuple(-x for x in r), np.conj(t)) for m, n, r, t in hop if r != (0, 0) or m != n]

    def h_true(k):
        h = np.zeros((2, 2), complex)
        for m, n, r, t in hop:
            h[m, n] += t * np.exp(1j * k @ (np.array(r + (0,)) @ uc + tau[n] - tau[m]))
        return h

    n = 4
    hrs = np.zeros((2, 2, n, n, 1, 1), complex)
    for m, nn, r, t in hop:                                               # PAOFLOW grid index of R = -R'
        hrs[m, nn, (-r[0]) % n, (-r[1]) % n, 0, 0] += t
    alat = a / bohr
    arrays = dict(HRs=hrs, a_vectors=uc / a, Dnm=np.repeat((tau / bohr)[:, None, :], 2, axis=1))
    fake = SimpleNamespace(data_controller=SimpleNamespace(data_dicts=lambda: (arrays, dict(alat=alat))))
    path = str(tmp_path / "hex.npz")
    with pytest.warns(RuntimeWarning):                                   # 2a hoppings sit on the grid edge
        write_paoflow_model(fake, path)
    model = load_paoflow_model(path)
    b = 2 * np.pi * np.linalg.inv(uc).T
    rng = np.random.default_rng(5)
    k = (rng.random((6, 3)) * [1, 1, 0]) @ b
    np.testing.assert_allclose(model.Hk(k), np.array([h_true(q) for q in k]), atol=1e-12)
    c3 = np.array([[-0.5, -np.sqrt(3) / 2, 0], [np.sqrt(3) / 2, -0.5, 0], [0, 0, 1]])
    for q in k:
        np.testing.assert_allclose(np.linalg.eigvalsh(model.Hk(q[None])[0]),
                                   np.linalg.eigvalsh(model.Hk((c3 @ q)[None])[0]), atol=1e-12)
