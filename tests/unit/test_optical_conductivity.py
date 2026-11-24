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
from lsquant.reconstruct import dielectric_function, optical_conductivity, shc


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


def script_sigma(data, size, half_width, orbitals, num_cells, cell_area, points, photon):
    """The reference script, transcribed: moments * N / 4, reduced energies, omega grid from zero."""
    scaled = data * orbitals / 4 * jackson_damping_2d(size, size)
    grid = np.linspace(-5.0 / half_width, 5.0 / half_width, points)
    omegas = np.concatenate(([0.0], photon / half_width))
    small = np.array([[np.sum(gamma_matrix(e, w, size, size) * scaled) for e in grid] for w in omegas])
    big = np.zeros_like(small)
    big[:, 1:] = np.cumsum(0.5 * (small[:, 1:] + small[:, :-1]) * np.diff(grid), axis=1)
    prefactor = 2 * (4 / np.pi) / (num_cells * cell_area) / half_width**2
    return grid, np.conj((big[1:] - big[0]) / omegas[1:, None] * prefactor)   # in e^2/hbar


def test_conductivity_matches_the_reference_script_formulas():
    """Real moments (one velocity): the script, value for value. Complex moments: the script applied
    to conj(mu), because the script's final conj returns sigma_ba for the moments of (v_a, v_b)."""
    size, half_width, orbitals = 12, 10.0, 50
    num_cells, cell_area, points = 25, 5.24, 64
    photon = np.array([1e-4, 0.5, 2.0, 4.0])
    for real in (True, False):
        mu = make_moments()
        if real:
            sym = np.asarray(mu.data).real
            mu.data[:] = sym + sym.T
        energies, sigma = optical_conductivity(
            mu, photon, energy_window=(-5.0, 5.0), num_cells=num_cells, cell_area=cell_area,
            num_energies=points, spin_degeneracy=2)
        data = np.asarray(mu.data) if real else np.conj(np.asarray(mu.data))
        grid, script = script_sigma(data, size, half_width, orbitals, num_cells, cell_area, points, photon)
        assert np.allclose(energies, grid * half_width)
        assert sigma.shape == (4, points)
        assert np.allclose(sigma, 2 * np.pi * script, rtol=1e-10, atol=1e-10)   # e^2/h = 2 pi e^2/hbar


def exact_trace_moments(h, va, vb, size, half_width):
    """mu_mn = 4 h_m h_n Tr[v_a T_m(H) v_b T_n(H)] / N, the normalisation LSQuant writes."""
    dim = h.shape[0]
    hs = h / half_width
    cheb = [np.eye(dim), hs]
    for _ in range(size - 2):
        cheb.append(2 * hs @ cheb[-1] - cheb[-2])
    half = np.where(np.arange(size) == 0, 0.5, 1.0)
    tr = np.array([[np.trace(va @ cheb[i] @ vb @ cheb[j]) for j in range(size)] for i in range(size)])
    mu = _core.ChebyshevMoments([size, size])
    mu.system_size = dim
    mu.window = _core.SpectralWindow(-half_width, half_width)
    mu.data[:] = 4 * np.outer(half, half) * tr / dim
    return mu


def bastin_from_eigenpairs(h, va, vb, size, half_width, grid, photon, num_cells, cell_area):
    """sigma_ab = i/(Omega omega)[S(omega) - S(0)] with the Jackson-damped expansions of delta and G+-
    evaluated on the eigenpairs: G+ with exp(-i n arccos x), no conjugation anywhere (e^2/h)."""
    eps, vec = np.linalg.eigh(h)
    a = vec.conj().T @ va @ vec
    b = vec.conj().T @ vb @ vec
    pair = a.T * b                                      # Tr[v_a P_i v_b P_j] = v^a_ji v^b_ij
    n = np.arange(size)
    jackson = ((size - n + 1) * np.cos(np.pi * n / (size + 1))
               + np.sin(np.pi * n / (size + 1)) / np.tan(np.pi / (size + 1))) / (size + 1)
    coef = (jackson * np.where(n == 0, 1.0, 2.0))[:, None] * np.cos(n[:, None] * np.arccos(eps / half_width))
    x = grid / half_width

    def delta(u):
        return coef.T @ np.cos(n[:, None] * np.arccos(u)) / (np.pi * np.sqrt(1 - u**2))

    def green(u, s):
        return -1j * s * coef.T @ np.exp(-1j * s * n[:, None] * np.arccos(u)) / np.sqrt(1 - u**2)

    def area(w):
        s = (np.einsum("ix,ij,jx->x", delta(x), pair, green(x + w, 1))
             + np.einsum("ix,ij,jx->x", green(x - w, -1), pair, delta(x)))
        out = np.zeros_like(s)
        out[1:] = np.cumsum(0.5 * (s[1:] + s[:-1]) * np.diff(x))
        return out

    base = area(0.0)
    return np.array([2 * np.pi * 1j * (area(w / half_width) - base)
                     / (num_cells * cell_area * half_width * w) for w in photon])


def test_exact_trace_moments_give_the_kubo_bastin_formula_for_every_pair():
    """The moments of an exact trace through optical_conductivity equal the Kubo-Bastin formula
    evaluated independently on the eigenpairs, for sigma_xx and for the Hall pair sigma_xy."""
    rng = np.random.default_rng(11)
    dim, size, half_width = 14, 24, 4.0
    herm = lambda m: (m + m.conj().T) / 2
    h = herm(rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim)))
    h *= 1.8 / np.max(np.abs(np.linalg.eigvalsh(h)))      # spectrum inside the energy window
    vx = herm(rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim)))
    vy = herm(rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim)))
    photon = np.array([1e-3, 0.3, 1.0, 1.7])
    for va, vb in ((vx, vx), (vx, vy), (vy, vx)):
        mu = exact_trace_moments(h, va, vb, size, half_width)
        energies, sigma = optical_conductivity(mu, photon, energy_window=(-2.2, 2.2), num_cells=3,
                                               cell_area=2.0, num_energies=96)
        want = bastin_from_eigenpairs(h, va, vb, size, half_width, energies, photon, 3, 2.0)
        assert np.allclose(sigma, want, rtol=1e-9, atol=1e-9 * np.abs(want).max())


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


def test_hall_limit_is_the_static_kubo_bastin_value():
    """sigma_xy(mu, omega -> 0) from optical_conductivity has the sign and size of the static curve."""
    rng = np.random.default_rng(4)
    dim, size, half_width = 16, 32, 4.0
    herm = lambda m: (m + m.conj().T) / 2
    h = herm(rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim)))
    h *= 1.8 / np.max(np.abs(np.linalg.eigvalsh(h)))      # spectrum inside the energy window
    vx = herm(rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim)))
    vy = herm(rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim)))
    mu = exact_trace_moments(h, vx, vy, size, half_width)
    energies, sigma = optical_conductivity(mu, [1e-4], energy_window=(-2.2, 2.2), num_cells=1,
                                           cell_area=1.0, num_energies=512)
    grid, static = shc(mu, 1e-3, 1, 1.0)            # 2 pi KuboBastin/(N_k A), e^2/h, full moments
    for e in (-0.8, -0.2, 0.4, 1.0):
        want = np.interp(e, grid, static)
        got = np.interp(e, energies, sigma[0].real)
        assert abs(got - want) < 0.03 * np.abs(static).max(), (e, got, want)
