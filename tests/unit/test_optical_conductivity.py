"""Optical conductivity reconstruction (finite-frequency Kubo-Bastin).

The reference is a direct NumPy transcription of the formulas of the
post-processing script this reconstruction was ported from (T. Galvani,
private communication, 2026): the Gamma matrix as an explicit
outer product, a term-by-term double sum, the cumulative trapezoid, the
subtraction of the omega = 0 row and the division by omega. The C++ kernel
uses a different algorithm (half-contractions and Horner evaluation), so the
agreement checks the algebra and not a copy of the code.
"""
import numpy as np
import pytest

from lsquant import _core
from lsquant.kernels import jackson_damping_2d
from lsquant.reconstruct import dielectric_function, optical_conductivity


def gamma_matrix(energy, omega, m0, m1):
    a0 = np.cos(np.arange(m0) * np.arccos(energy)) / np.sqrt(1 - energy**2)
    a1 = np.cos(np.arange(m1) * np.arccos(energy)) / np.sqrt(1 - energy**2)
    plus = np.exp(1j * np.arange(m1) * np.arccos(energy + omega)) / np.sqrt(1 - (energy + omega) ** 2)
    minus = -np.exp(-1j * np.arange(m0) * np.arccos(energy - omega)) / np.sqrt(1 - (energy - omega) ** 2)
    return np.outer(a0, plus) + np.outer(minus, a1)


def random_moments(m0, m1, seed=3):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(m0, m1)) + 1j * rng.normal(size=(m0, m1))


def test_kernel_matches_the_direct_double_sum():
    mu = random_moments(9, 6)
    energies = np.linspace(-0.6, 0.55, 7)
    omegas = np.array([0.0, 1e-5, 0.1, 0.35])
    got = _core.ContractRank2Optical(mu, energies, omegas)
    assert got.shape == (4, 7)
    want = np.array([[np.sum(gamma_matrix(e, w, 9, 6) * mu) for e in energies] for w in omegas])
    assert np.allclose(got, want, rtol=1e-12, atol=1e-12)


def test_kernel_refuses_energies_outside_the_window():
    mu = random_moments(4, 4)
    with pytest.raises(ValueError, match="spectral window"):
        _core.ContractRank2Optical(mu, np.array([0.8]), np.array([0.3]))
    with pytest.raises(ValueError, match="spectral window"):
        _core.ContractRank2Optical(mu, np.array([-0.8]), np.array([0.3]))


def make_moments(size=12, half_width=10.0):
    mu = _core.ChebyshevMoments([size, size])
    mu.system_size = 50
    mu.window = _core.SpectralWindow(-half_width, half_width)
    mu.data[:] = random_moments(size, size, seed=5)
    return mu


def test_conductivity_matches_the_reference_script_formulas():
    mu = make_moments()
    size, half_width, orbitals = 12, 10.0, 50
    num_cells, cell_area, points = 25, 5.24, 64
    photon = np.array([1e-4, 0.5, 2.0, 4.0])
    energies, sigma = optical_conductivity(
        mu, photon, energy_window=(-5.0, 5.0), num_cells=num_cells, cell_area=cell_area,
        num_energies=points, spin_degeneracy=2)

    # The script: moments * N / 4, reduced energies, omega grid that starts at zero.
    scaled = np.asarray(mu.data) * orbitals / 4 * jackson_damping_2d(size, size)
    grid = np.linspace(-5.0 / half_width, 5.0 / half_width, points)
    omegas = np.concatenate(([0.0], photon / half_width))
    small = np.array([[np.sum(gamma_matrix(e, w, size, size) * scaled) for e in grid] for w in omegas])
    big = np.zeros_like(small)
    big[:, 1:] = np.cumsum(0.5 * (small[:, 1:] + small[:, :-1]) * np.diff(grid), axis=1)
    prefactor = 2 * (4 / np.pi) / (num_cells * cell_area) / half_width**2
    script = np.conj((big[1:] - big[0]) / omegas[1:, None] * prefactor)   # in e^2/hbar

    assert np.allclose(energies, grid * half_width)
    assert sigma.shape == (4, points)
    assert np.allclose(sigma, 2 * np.pi * script, rtol=1e-10, atol=1e-10)   # e^2/h = 2 pi e^2/hbar


def test_fewer_moments_and_default_grid():
    mu = make_moments()
    energies, sigma = optical_conductivity(
        mu, [1.0], energy_window=(-5.0, 5.0), num_cells=1, cell_area=1.0, num_moments=8)
    assert energies.size == 16 and sigma.shape == (1, 16)


@pytest.mark.parametrize("kwargs, text", [
    (dict(photon_energies=[0.0]), "positive"),
    (dict(photon_energies=[6.0]), "wider expansion.window"),
    (dict(energy_window=(2.0, -2.0)), "lower < upper"),
    (dict(num_moments=13), "num_moments"),
    (dict(num_cells=0), "num_cells"),
])
def test_bad_requests_are_refused(kwargs, text):
    args = dict(photon_energies=[1.0], energy_window=(-5.0, 5.0), num_cells=1, cell_area=1.0)
    args.update(kwargs)
    photon = args.pop("photon_energies")
    with pytest.raises(ValueError, match=text):
        optical_conductivity(make_moments(), photon, **args)


def test_rank_one_moments_are_refused():
    mu = _core.ChebyshevMoments([8])
    mu.window = _core.SpectralWindow(-1.0, 1.0)
    with pytest.raises(ValueError, match="rank-2"):
        optical_conductivity(mu, [0.1], energy_window=(-0.5, 0.5), num_cells=1, cell_area=1.0)


def test_dielectric_function_of_a_sheet():
    sigma = np.array([[1.0 + 0.5j, 2.0], [0.25j, 1.0]])
    eps = dielectric_function(sigma, [1.0, 2.0], thickness=3.35)
    assert np.allclose(eps[0], 1 + 2j * 14.399645 * sigma[0] / (1.0 * 3.35))
    assert np.allclose(eps[1], 1 + 2j * 14.399645 * sigma[1] / (2.0 * 3.35))
    with pytest.raises(ValueError):
        dielectric_function(sigma, [1.0, 2.0], thickness=0.0)
