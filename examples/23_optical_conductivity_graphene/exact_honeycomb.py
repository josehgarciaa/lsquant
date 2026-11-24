#!/usr/bin/env python3
"""Exact optical response of the nearest-neighbour honeycomb model on a dense k mesh.

Problem: the tests of benchmark.py need references that do not come from the
Chebyshev moments: the absorption Re sigma_aa(omega) of graphene and hBN at
any chemical potential, and the clean polarisability of hBN below its gap.
Approach: the two-band Bloch Hamiltonian of make_honeycomb.bloch, whose
derivative is the Bloch form of v = i[H, r], diagonalised on an n x n mesh
(midpoints, so no k point sits on a Dirac point). With hbar = 1, velocities
in eV Angstrom and areas in Angstrom^2:

    Re sigma_aa(omega) = g_s pi / (omega N_k A) sum_k |<c|v_a|v>|^2 (f_v - f_c) delta(omega - e_cv)   [e^2/hbar]
    alpha(omega)       = 2 g_s / (N_k A) sum_k |<c|v_a|v>|^2 / (e_cv (e_cv^2 - omega^2)),  Im sigma = -omega alpha

(an insulator at mu in the gap; derivation.py, step 3). delta is a Gaussian
of the given width. The functions return e^2/h (times 2 pi) and, for alpha,
the 2D polarisability as a length: alpha_2D = e^2 alpha / (4 pi eps_0) with
e^2/(4 pi eps_0) = 14.399645 eV Angstrom, so that epsilon = 1 + 4 pi alpha_2D/d
for a sheet of thickness d.
"""
import numpy as np

from make_honeycomb import bloch, lattice

E2 = 14.399645          # e^2 / (4 pi eps_0) in eV Angstrom


class Honeycomb:
    """Bands and interband velocity matrix elements on an n x n midpoint mesh."""

    def __init__(self, delta, t, a, n_k=1200, chunk=250_000):
        self.area = lattice(a)[2]
        a1, a2, _, _ = lattice(a)
        recip = 2 * np.pi * np.linalg.inv(np.array([a1, a2])).T
        u = (np.arange(n_k) + 0.5) / n_k
        u1, u2 = np.meshgrid(u, u, indexing="ij")
        k = np.column_stack((u1.ravel(), u2.ravel())) @ recip
        self.n_k = k.shape[0]
        ev, ec, vx2, vy2 = [], [], [], []
        for s in range(0, self.n_k, chunk):
            h, dx, dy = bloch(k[s:s + chunk], delta, t, a)
            e, c = np.linalg.eigh(h)
            cd = np.conj(np.swapaxes(c, 1, 2))
            ev.append(e[:, 0])
            ec.append(e[:, 1])
            vx2.append(np.abs((cd @ dx @ c)[:, 1, 0]) ** 2)
            vy2.append(np.abs((cd @ dy @ c)[:, 1, 0]) ** 2)
        self.ev, self.ec = np.concatenate(ev), np.concatenate(ec)
        self.v2 = {"x": np.concatenate(vx2), "y": np.concatenate(vy2)}

    def re_sigma(self, photon, mu, width, axis="x", spin=2):
        """Re sigma_aa(omega) in e^2/h at chemical potential mu (zero temperature), Gaussian delta."""
        weight = self.v2[axis] * ((self.ev < mu).astype(float) - (self.ec < mu).astype(float))
        keep = weight != 0
        gap, weight = (self.ec - self.ev)[keep], weight[keep]
        out = np.empty(len(photon))
        for i, w in enumerate(photon):
            near = np.abs(gap - w) < 8 * width
            g = np.exp(-0.5 * ((w - gap[near]) / width) ** 2) / (width * np.sqrt(2 * np.pi))
            out[i] = spin * np.pi / (w * self.n_k * self.area) * np.sum(weight[near] * g)
        return 2 * np.pi * out

    def polarisability(self, photon, axis="x", spin=2):
        """Clean alpha_2D(omega) in Angstrom for mu in the gap; valid below the gap."""
        gap = self.ec - self.ev
        return np.array([E2 * 2 * spin / (self.n_k * self.area)
                         * np.sum(self.v2[axis] / (gap * (gap**2 - w**2))) for w in photon])

    def im_sigma_insulator(self, photon, axis="x", spin=2):
        """Im sigma_aa(omega) in e^2/h of the clean insulator below the gap: -omega alpha."""
        return -2 * np.pi * np.asarray(photon) * self.polarisability(photon, axis, spin) / E2
