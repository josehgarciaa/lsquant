"""Anomalous Hall word in the run file of a DFT-model executor.

A Haldane Chern insulator (C = +-1), written as a PAOFLOW model file, runs
through ``execute_runfile`` with ``Tr[VX delta(E-H) VY dG-(E)/dE]``:

- the ``ahc`` curve is the quantised plateau in the gap, and the exact
  occupied-band sum (``exact_ahc.json``) gives the same value;
- the ``sea`` and ``surf`` parts add up to the full curve, and only the sea
  part survives in the gap.
"""
import json
import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant.runfile import execute_runfile  # noqa: E402


def _haldane_npz(path, t1=-1.0, t2=0.15, phi=np.pi / 2):
    """Honeycomb lattice, A at 0 and B at (a1+a2)/3; H_ij(R) = <i,0|H|j,R>."""
    uc = np.array([[1.0, 0.0, 0.0], [0.5, np.sqrt(3) / 2, 0.0], [0.0, 0.0, 10.0]])
    tau = np.array([[0.0, 0.0, 0.0], (uc[0] + uc[1]) / 3])
    blocks = {}

    def add(R, i, j, value):
        # Each bond once, with its Hermitian partner H_ji(-R) = conj(H_ij(R)).
        minus_r = tuple(-x for x in R)
        blocks.setdefault(R, np.zeros((2, 2), complex))[i, j] += value
        blocks.setdefault(minus_r, np.zeros((2, 2), complex))[j, i] += np.conj(value)

    for R in [(0, 0, 0), (-1, 0, 0), (0, -1, 0)]:
        add(R, 0, 1, t1)
    for R in [(1, 0, 0), (-1, 1, 0), (0, -1, 0)]:
        add(R, 0, 0, t2 * np.exp(1j * phi))
        add(R, 1, 1, t2 * np.exp(-1j * phi))
    r_int = np.array(list(blocks))
    hoppings = np.array([blocks[tuple(r)] for r in r_int])
    np.savez(path, uc=uc, tau=tau, r_int=r_int, hoppings=hoppings,
             spin=np.zeros((3, 2, 2)), fermi_energy=0.0, mesh=np.array([3, 3, 1]), edge_ratio=0.0)


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("ahc")
    model = str(tmp / "haldane.npz")
    _haldane_npz(model)
    rf = {
        "model": {"kind": "paoflow", "file": model, "params": {"oracle": True}},
        "basis": {"name": "eigen", "supercell": [36, 36, 1]},
        "metric": {"policy": "orthogonal"},
        "word": "Tr[VX delta(E-H) VY dG-(E)/dE]",
        "expansion": {"moments": 512, "window": "eigen", "trace": "exact", "seed": 1},
        "reconstruction": {"kernel": "jackson", "formulas": ["ahc", "sea", "surf"],
                           "broadening_meV": 20.0, "energy_grid": [-0.3, 0.3, 61],
                           "fermi_energies": [0.0]},
        "output": {"legacy_layout": False, "dir": str(tmp / "out")},
    }
    execute_runfile(rf)
    return tmp / "out"


def _curve(path, energy):
    e, v = np.loadtxt(path).T
    return float(np.interp(energy, e, v))


def test_plateau_is_quantised_and_matches_exact_sum(run):
    plateau = _curve(run / "ahc.dat", 0.0)
    assert abs(abs(plateau) - 1.0) < 0.02
    exact = json.load(open(run / "exact_ahc.json"))
    assert exact["units"] == "e^2/h"
    assert abs(np.interp(0.0, exact["energies"], exact["values"]) - plateau) < 0.02


def test_sea_and_surface_add_up(run):
    for energy in (-2.0, -1.0, 0.0, 1.0):
        total = _curve(run / "ahc.dat", energy)
        sea, surf = _curve(run / "ahc_sea.dat", energy), _curve(run / "ahc_surf.dat", energy)
        # The three curves sit on slightly different energy grids: compare after interpolation.
        assert abs(sea + surf - total) < 1e-3
    assert abs(_curve(run / "ahc_surf.dat", 0.0)) < 1e-3
