#!/usr/bin/env python3
"""Derive the optical-conductivity formula of the reference script with SymPy, step by step.

The reference is the post-processing script of T. Galvani (private
communication, 2026), which `lsquant.reconstruct.optical_conductivity`
reproduces. The script evaluates

    sigma(mu, omega) = P conj[ (G(mu, w) - G(mu, 0)) / w ],
    G(mu, w)         = int_{x_min}^{mu/D} dx  sum_mn  mu~_mn g_m g_n Gamma_mn(x, w),
    Gamma_mn(x, w)   = T_m(x)/sqrt(1-x^2) * exp(+i n arccos(x+w))/sqrt(1-(x+w)^2)
                     - exp(-i m arccos(x-w))/sqrt(1-(x-w)^2) * T_n(x)/sqrt(1-x^2),
    P                = g_s (4/pi) / (A D^2),        mu~_mn = (N/4) mu_mn(file),

with x = E/D, w = hbar omega/D, D the half-width of the spectral window, A the
area of the supercell, g_m the Jackson kernel and mu_mn the moments that
LSQuant writes. This script proves each step from the finite-frequency
Kubo-Bastin formula,

    sigma_ab(mu, omega) = i / (Omega omega) [S_ab(mu, omega) - S_ab(mu, 0)],
    S_ab(mu, omega)     = int_{-inf}^{mu} dE Tr[ v_a delta(E-H) v_b G+(E+omega)
                                                + v_a G-(E-omega) v_b delta(E-H) ],

(hbar = 1, e^2/hbar units, Omega the area of the sample), and states what it
implies for the conventions: the moment normalisation, the complex
conjugate, the order of the two velocities, the units and the dielectric
function. Each step ends in an assertion; the script prints PASS when all
hold. It runs in a few seconds and needs only SymPy and NumPy.

Usage:
    python derivation.py
"""
import numpy as np
import sympy as sp

I, pi = sp.I, sp.pi


def step(title):
    print(f"\n== {title}")


def check(expr, what):
    """Assert that a SymPy expression simplifies to zero."""
    assert sp.simplify(expr) == 0, what
    print("   ok:", what)


# ---------------------------------------------------------------------------
step("1. Chebyshev expansion of delta(x - y): exact for every polynomial of degree < M")
x, y = sp.symbols("x y", real=True)
M = 6


def c(n):
    return 1 if n == 0 else 2          # (2 - delta_n0) = 4 h_n with h_0 = 1/2, h_n = 1


# delta_M(x - y) = sum_{n<M} c_n T_n(y) T_n(x) / (pi sqrt(1-x^2)); integrate against p(x).
p = 3 * x**5 - x**4 + 2 * x**2 - 7 * x + 1
theta = sp.symbols("theta", real=True)
moments = [sp.integrate((p * sp.chebyshevt(n, x)).subs(x, sp.cos(theta)), (theta, 0, pi)) for n in range(M)]
projected = sum(c(n) * sp.chebyshevt(n, y) * moments[n] / pi for n in range(M))
check(sp.expand(projected - p.subs(x, y)), "int p(x) delta_M(x - y) dx = p(y) for deg p < M")

# ---------------------------------------------------------------------------
step("2. Chebyshev expansion of G+(x) = 1/(x - y + i0): coefficients -i c_n exp(-i n arccos x)/sqrt(1-x^2)")
r, phi, eps = sp.symbols("r phi epsilon", positive=True)
th = sp.symbols("theta", positive=True)
z = r * sp.exp(-I * th)
# sum_n c_n r^n exp(-i n theta) cos(n phi) = -1 + sum_n z^n e^{i n phi} + sum_n z^n e^{-i n phi}  (|z| < 1)
series = -1 + 1 / (1 - z * sp.exp(I * phi)) + 1 / (1 - z * sp.exp(-I * phi))
at_one = series.subs(r, 1)
check(sp.simplify((-I / sp.sin(th) * at_one).rewrite(sp.exp) - (1 / (sp.cos(th) - sp.cos(phi))).rewrite(sp.exp)),
      "-i/sin(theta) * sum_n c_n e^{-i n theta} T_n(y) = 1/(x - y), x = cos theta, y = cos phi")
# The regulator r = e^{-eps} < 1 is theta -> theta - i eps, i.e. x -> cos(theta - i eps):
check(series.subs(r, sp.exp(-eps)) - at_one.subs(th, th - I * eps), "r = exp(-eps) is theta -> theta - i eps")
shift = sp.expand_trig(sp.cos(th - I * eps))
im_shift = sp.im(sp.expand_complex(shift))
check(im_shift - sp.sin(th) * sp.sinh(eps), "Im cos(theta - i eps) = sin(theta) sinh(eps) > 0 on 0 < theta < pi")
print("   so the series with exp(-i n arccos x) converges to 1/(x - y + i0) = G+ ;"
      " its complex conjugate, exp(+i n arccos x), gives G- = 1/(x - y - i0).")

# ---------------------------------------------------------------------------
step("3. Kubo-Bastin S(mu, omega) in the eigenbasis is the Kubo-Greenwood conductivity")
E, mu, w = sp.symbols("E mu omega", real=True)
ea, eb, eta, Om = sp.symbols("epsilon_a epsilon_b eta Omega", real=True)
K = sp.symbols("K")       # v^a_{ab} v^b_{ba} for one pair of eigenstates (a, b)
# Tr[v_a delta(E-H) v_b G+(E+omega) + v_a G-(E-omega) v_b delta(E-H)] for the pair (a, b):
integrand = K * (sp.DiracDelta(E - ea) / (E + w - eb + I * eta) + sp.DiracDelta(E - eb) / (E - w - ea - I * eta))
S = sp.integrate(integrand, (E, -sp.oo, mu))
fa, fb = sp.Heaviside(mu - ea), sp.Heaviside(mu - eb)
check(S - K * (fa - fb) / (w + ea - eb + I * eta), "S = K (f_a - f_b)/(omega + e_a - e_b + i0)")
sigma_bastin = I / (Om * w) * (S - S.subs(w, 0))
sigma_kg = I / Om * K * (fa - fb) / ((eb - ea) * (w + ea - eb + I * eta))
check((sigma_bastin - sigma_kg).subs(eta, 0), "i [S(omega) - S(0)]/(Omega omega) = Kubo-Greenwood, e^{-i omega t}")
print("   Re sigma_xx = pi/(Omega omega) sum (f_a - f_b)|v_ab|^2 delta(omega - e_b + e_a) >= 0: absorption is positive.")
# Below the gap (one valence a, one conduction b, |v|^2 = K, e = e_b - e_a > omega):
e = sp.symbols("e", positive=True)
kg = lambda f1, f2, e1, e2_: I / Om * K * (f1 - f2) / ((e2_ - e1) * (w + e1 - e2_))
both = kg(1, 0, ea, ea + e) + kg(0, 1, ea + e, ea)          # pairs (v, c) and (c, v), f_v = 1, f_c = 0
check(both - (-I * w) * (2 * K / (Om * e * (e**2 - w**2))),
      "insulator: sigma = -i omega alpha, alpha = 2|v|^2/(Omega e (e^2 - omega^2)) > 0, so Im sigma < 0")
# The static limit is Bastin's DC formula: d/domega G+(E+omega) = dG+/dE, d/domega G-(E-omega) = -dG-/dE.
check(sp.limit(sigma_bastin.subs(eta, 0), w, 0) - I / Om * K * (fa - fb) * (-1) / (ea - eb)**2,
      "sigma(omega -> 0) = (i/Omega) K (f_a - f_b)(-1)/(e_a - e_b)^2 = the static Kubo-Bastin term")

# ---------------------------------------------------------------------------
step("4. Insert the expansions: the Gamma matrix, the complex conjugate and the prefactor")
m, n = sp.symbols("m n", integer=True, nonnegative=True)
xx, ww = sp.symbols("x w", real=True)
# With |x| < 1 and |x +- w| < 1 the angles arccos(x), arccos(x + w), arccos(x - w) and the square
# roots sqrt(1 - u^2) are real; name them a, b, c and s0, sp_, sm.
a, b, cc = sp.symbols("a b c", real=True)
s0, spl, smi = sp.symbols("s_0 s_+ s_-", positive=True)
Tm, Tn = sp.cos(m * a), sp.cos(n * a)
delta_m = lambda T: T / (pi * s0)                       # delta(x - H) term, times c_m T_m(H)
g_plus = -I * sp.exp(-I * n * b) / spl                  # G+(x + w - H), times c_n T_n(H)   (step 2)
g_minus = I * sp.exp(I * m * cc) / smi                  # G-(x - w - H), times c_m T_m(H)
# Coefficient of c_m c_n Tr[v_a T_m(H) v_b T_n(H)] in the integrand s(x, w):
coeff = delta_m(Tm) * g_plus + g_minus * delta_m(Tn)
gamma = Tm / s0 * sp.exp(I * n * b) / spl - sp.exp(-I * m * cc) / smi * Tn / s0     # the script's Gamma_mn
check(sp.expand(coeff - (-I / pi) * sp.conjugate(gamma)), "coefficient = (-i/pi) conj(Gamma_mn(x, w))")
print("   With c_m c_n Tr[v_a T_m v_b T_n] = N mu_mn(file) = 4 mu~_mn (step 5), the integrand is\n"
      "   s(x, w) = -(4i/pi) sum_mn mu~_mn conj(Gamma_mn(x, w)).")
# sigma = g_s i/(Omega omega)[S(omega) - S(0)], S = (1/D) int dx s, omega = D w:
D, A, gs, Gw, G0 = sp.symbols("D A g_s G_w G_0")
sigma = gs * I / (A * D * ww) * (1 / D) * (-4 * I / pi) * (Gw - G0)       # G = int sum mu~ conj(Gamma)
P = gs * 4 / (pi * A * D**2)
check(sigma - P * (Gw - G0) / ww, "sigma = P [G*(w) - G*(0)]/w,  P = g_s (4/pi)/(A D^2),  A = supercell area")
print("   G*(w) = int sum mu~_mn conj(Gamma_mn) = conj( int sum conj(mu~_mn) Gamma_mn ).\n"
      "   The script takes conj( int sum mu~_mn Gamma_mn ): the same when mu~ is real.")

# ---------------------------------------------------------------------------
step("5. The moments: normalisation, reality and the order of the velocities (numerical, exact trace)")
rng = np.random.default_rng(7)
dim, nm = 9, 6
h = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
h = (h + h.conj().T) / 8
va = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
va = va + va.conj().T
vb = rng.normal(size=(dim, dim)) + 1j * rng.normal(size=(dim, dim))
vb = vb + vb.conj().T
hs = h / (1.01 * np.max(np.abs(np.linalg.eigvalsh(h))))
cheb = [np.eye(dim), hs]
for _ in range(nm - 2):
    cheb.append(2 * hs @ cheb[-1] - cheb[-2])
tr = lambda a, b: np.array([[np.trace(a @ cheb[i] @ b @ cheb[j]) for j in range(nm)] for i in range(nm)])
mu_ab, mu_ba, mu_aa = tr(va, vb), tr(vb, va), tr(va, va)
assert np.allclose(mu_aa.imag, 0) and np.allclose(mu_aa, mu_aa.T)
print("   ok: Tr[v T_m v T_n] is real and symmetric (one velocity): the conj of the script is exact for sigma_xx")
assert np.allclose(mu_ab.conj(), mu_ba)
print("   ok: conj(Tr[v_a T_m v_b T_n]) = Tr[v_b T_m v_a T_n]: with the moments of (v_a, v_b) the script's\n"
      "       conj returns sigma_ba, i.e. the Hall part with the opposite sign. LSQuant conjugates the moments\n"
      "       first, so optical_conductivity returns sigma_ab for every pair.")
print("   LSQuant writes mu_mn(file) = 4 h_m h_n Tr[v_a T_m(H) v_b T_n(H)] / N (h_0 = 1/2, h_n = 1),\n"
      "   checked to 1e-14 on an exact trace in benchmark.py; so mu~ = (N/4) mu(file) = h_m h_n Tr[...],\n"
      "   and c_m c_n Tr[...] = 4 h_m h_n Tr[...] = 4 mu~_mn, as used in step 4.")

# ---------------------------------------------------------------------------
step("6. Units and the dielectric function")
e2, hw, d, sig_hbar = sp.symbols("e2 hbar_omega d sigma_hat", positive=True)   # e2 = e^2/(4 pi eps0) = 14.4 eV A
# Velocities in eV A, energies in eV, areas in A^2 with hbar = 1: P (D w)^-1 ... is dimensionless, in e^2/hbar.
sig_h = 2 * pi * sig_hbar                                    # e^2/h = (e^2/hbar)/(2 pi)
eps_gauss = 1 + 4 * pi * I * e2 * sig_hbar / (hw * d)        # 1 + 4 pi i sigma_3D/omega, sigma_3D = sigma/d
eps_script = 1 + I * sig_hbar * 4 * pi * e2 / (d * D) * D / hw   # i conj[...]/w^2 P 4 pi e2/(H D), w = hbar_omega/D
eps_lsquant = 1 + 2 * I * e2 * sig_h / (hw * d)              # dielectric_function(sigma in e^2/h)
check(eps_script - eps_gauss, "script: epsilon = 1 + 4 pi i e^2 sigma/(hbar omega d), sigma in e^2/hbar")
check(eps_lsquant - eps_gauss, "LSQuant: epsilon = 1 + 2 i e^2 sigma/(hbar omega d), sigma in e^2/h")
alpha = sp.symbols("alpha_2D", positive=True)
check(sp.im(sp.expand(eps_gauss.subs(sig_hbar, -I * hw * alpha / e2))), "sigma = -i omega alpha: epsilon real")
print("   2D polarisability (a length): alpha_2D = -e^2 Im(sigma)/(hbar omega), sigma in e^2/hbar;\n"
      "   epsilon = 1 + 4 pi alpha_2D/d for a sheet of thickness d.")

print("\nPASS: the script's formula is the finite-frequency Kubo-Bastin formula, term by term.")
