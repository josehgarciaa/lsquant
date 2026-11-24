"""Check dimensionless spectral kernels against independent closed forms.

Small, seeded grids avoid endpoint singularities and keep every table tiny.
Pitfall: these are coefficient kernels; contraction supplies the factor two
for nonzero Chebyshev orders, rather than the table itself.
"""

import math
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant import spectral_functions as sf


def _grid():
    return np.random.default_rng(12345).uniform(-0.95, 0.95, 13)


@pytest.mark.parametrize("kernel", [sf.DeltaFunction, sf.legacy_delta_table])
def test_delta_closed_form(kernel):
    r"""Check :math:`c_m=T_m(x)/(\pi\sqrt{1-x^2})`, with half the weight.

    Here :math:`T_m(x)=\cos(m\arccos x)` and the full Chebyshev weight is
    :math:`2/(\pi\sqrt{1-x^2})`. Scalar math is independent of lambdify.
    """
    x = _grid()
    expected = [[math.cos(m * math.acos(t)) / (math.pi * math.sqrt(1-t*t))
                 for t in x] for m in range(12)]
    actual = kernel.table(x, 12) if hasattr(kernel, "table") else kernel(x, 12)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    assert actual.shape == (12, 13)
    assert actual.flags.c_contiguous


@pytest.mark.parametrize("legacy", [False, True])
def test_green_closed_forms(legacy):
    r"""Check retarded coefficients and their analytic energy derivative.

    With :math:`s=\sqrt{1-x^2}` and :math:`q=e^{-im\arccos x}`,
    :math:`g_m=-iq/s` and :math:`g'_m=q(m/s^2-ix/s^3)`.
    Unequal table lengths guard the two independent moment axes.
    """
    x = _grid()
    s = np.sqrt(1-x*x)
    def expected(n, derivative=False):
        m = np.arange(n)[:, None]
        # Euler's formula avoids the production complex power/exponential.
        angle = m * np.arccos(x)
        q = np.cos(angle) - 1j * np.sin(angle)
        return q * (m / s**2 - 1j*x / s**3) if derivative else -1j*q/s
    if legacy:
        green, derivative = sf.legacy_green_tables(x, 9, 12)
    else:
        green = sf.GreenRetarded.table(x, 9)
        derivative = sf.GreenRetardedDerivative.table(x, 12)
    np.testing.assert_allclose(green, expected(9), rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(derivative, expected(12, True), rtol=1e-12, atol=1e-12)
    default = sf.legacy_green_tables(x, 9)[1] if legacy else derivative[:9]
    np.testing.assert_allclose(default, expected(9, True), rtol=1e-12, atol=1e-12)


def test_gaussian_closed_form():
    r"""Check :math:`e^{-x^2/(2\eta^2)}/(\sqrt{2\pi}\eta)`.

    GaussianDelta uses its parameter axis for width, not moment order.
    """
    x = _grid()
    eta = 0.3
    expected = np.exp(-x*x/(2*eta*eta)) / (math.sqrt(2*math.pi)*eta)
    np.testing.assert_allclose(sf.GaussianDelta.evaluate(x, eta), expected,
                               rtol=1e-12, atol=1e-12)


def test_time_evolution_series():
    r"""Check :math:`\sum_m c_m(t)T_m(x)=e^{-ixt}` at small times.

    A scalar factorial series independently evaluates :math:`J_m(t)`;
    24 terms suffice at these dimensionless times without large arrays.
    """
    times = np.array([0.0, 0.2, 0.7])
    expected = np.array([
        [(1 if m == 0 else 2) * (-1j)**m * sum(
            (-1)**k * (t/2)**(2*k+m) /
            (math.factorial(k)*math.factorial(k+m)) for k in range(24))
         for t in times] for m in range(18)])
    np.testing.assert_allclose(sf.TimeEvolution.table(times, 18), expected,
                               rtol=1e-12, atol=1e-12)
    x = _grid()
    polynomials = np.cos(np.arange(18)[:, None]*np.arccos(x))
    np.testing.assert_allclose(expected.T @ polynomials,
                               np.exp(-1j*times[:, None]*x),
                               rtol=1e-12, atol=1e-12)
    for t in times:
        coeff = sf.bessel_series_coefficients(t)
        np.testing.assert_allclose(coeff, expected[:len(coeff), list(times).index(t)],
                                   rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("count", [0, -1])
def test_nonpositive_moment_count(count):
    """Reject empty expansions before allocating a coefficient table."""
    for kernel in (sf.DeltaFunction, sf.GreenRetarded,
                   sf.GreenRetardedDerivative, sf.TimeEvolution):
        with pytest.raises(ValueError, match="positive"):
            kernel.table([0.2], count)
    with pytest.raises(ValueError, match="positive"):
        sf.legacy_delta_table([0.2], count)
    with pytest.raises(ValueError, match="positive"):
        sf.legacy_green_tables([0.2], 2, count)


def test_symbolic_constant_broadcast():
    """Broadcast a constant kernel over both axes without aliasing output.

    Pitfall: lambdify returns a scalar for expressions with neither symbol;
    the public table still promises a contiguous two-dimensional array.
    """
    import sympy as sp
    kernel = sf.SpectralFunction(sp.Integer(3))
    actual = kernel.table([0.1, 0.2], 4)
    np.testing.assert_array_equal(actual, np.full((4, 2), 3.0))
    assert actual.dtype == np.float64
    assert actual.flags.c_contiguous


@pytest.mark.parametrize("time", [0.0, 0.2, 0.7])
def test_legacy_bessel_closed_form(time):
    r"""Check :math:`c_m=(2-\delta_{m0})(-i)^m J_m(t)` via factorials.

    The independent convergent series is inexpensive at these small times.
    The first sub-tolerance Bessel term is excluded from the returned series.
    """
    def bessel(m):
        return sum((-1)**k * (time/2)**(2*k+m) /
                   (math.factorial(k)*math.factorial(k+m))
                   for k in range(24))
    expected = [bessel(0)]
    for m in range(1, 18):
        if abs(bessel(m)) <= 1e-15:
            break
        expected.append(2*(-1j)**m*bessel(m))
    np.testing.assert_allclose(sf.legacy_bessel_coefficients(time), expected,
                               rtol=1e-12, atol=1e-12)
