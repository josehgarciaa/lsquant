"""Atom-pair blocks of a machine-learning Hamiltonian to an LSQuant model.

Machine-learning electronic-structure codes predict H(R), and often take
S(R) from the DFT code, as one block per atom pair: ``M_{ij}(R)`` couples
the orbitals of atom ``i`` in the home cell with those of atom ``j`` in
cell ``R``. The readers of each code (``models/deeph.py`` and the others)
convert their file layout to the convention used here, and this module
builds the model that :class:`~lsquant.bases.BlochMetricBasis` reads.

Conventions of this module:

- ``M_{ij}(R) = <i, 0| M | j, R>``, so that
  ``M(k) = sum_R exp(+i k.(R + tau_j - tau_i)) M(R)`` (atom gauge of
  :func:`~lsquant.models.bloch.bloch_sum`). Readers state how they reach it.
- Energies in eV, lengths in Angstrom, lattice vectors as rows.
- Every orbital sits on its atom: ``tau`` repeats the atomic position.
  The order of the orbitals inside an atom does not enter LSQuant, because
  no operator acts on the angular part of an orbital.
- Spinful blocks of one atom pair are ``[[uu, ud], [du, dd]]`` with
  ``n_i x n_j`` sub-blocks (the layout of DeepH and HamGNN); they are
  stored with interleaved spinors (orbital up, orbital down), the layout of
  the SIESTA model, so that the Pauli matrices on the spinor labels are the
  spin operator.
"""
from __future__ import annotations

import numpy as np

from .bloch import bloch_sum

_SPIN = (0, 1)


def orbital_slices(norb, spinful):
    """Return ``rows[a]``, the interleaved global indices of atom ``a``.

    For a spinful atom the indices follow the pair-block order of its rows:
    all up components, then all down components.
    """
    norb = np.asarray(norb, int)
    starts = np.concatenate([[0], np.cumsum(norb)])
    rows = []
    for a, n in enumerate(norb):
        spinless = starts[a] + np.arange(n)
        rows.append(np.concatenate([2 * spinless + s for s in _SPIN]) if spinful else spinless)
    return rows


def assemble_blocks(pairs, norb, *, spinful, name="H"):
    """Sum atom-pair blocks into ``(r_int, M)`` with ``M`` of shape ``(N_R, W, W)``.

    ``pairs`` yields ``(R, i, j, block)`` with integer ``R`` (three lattice
    indices), 0-based atoms and ``block`` of shape ``(n_i, n_j)`` or, for a
    spinful model, ``(2 n_i, 2 n_j)`` in the ``[[uu, ud], [du, dd]]`` layout.
    A repeated ``(R, i, j)`` is an error, not a sum.
    """
    rows = orbital_slices(norb, spinful)
    width = sum(len(r) for r in rows)
    index, mats, seen = {}, [], set()
    for R, i, j, block in pairs:
        key = (tuple(int(x) for x in R), int(i), int(j))
        if key in seen:
            raise ValueError(f"{name}: block R={key[0]} i={key[1]} j={key[2]} appears twice")
        seen.add(key)
        block = np.asarray(block)
        if block.shape != (len(rows[key[1]]), len(rows[key[2]])):
            raise ValueError(f"{name}: block R={key[0]} i={key[1]} j={key[2]} has shape {block.shape}, "
                             f"expected {(len(rows[key[1]]), len(rows[key[2]]))}")
        if key[0] not in index:
            index[key[0]] = len(mats)
            mats.append(np.zeros((width, width), complex))
        mats[index[key[0]]][np.ix_(rows[key[1]], rows[key[2]])] = block
    if not mats:
        raise ValueError(f"{name}: no blocks")
    return np.array(list(index), int), np.array(mats)


def spinor_overlap(overlap):
    """Lift a spinless ``S(R)`` to interleaved spinors, ``S tensor I_2``."""
    return np.array([np.kron(s, np.eye(2)) for s in overlap])


def hermitian_pairs(r_int, mats):
    """Return ``(r_int, M, residual)`` with ``M(-R) = M(R)^+`` for every R.

    A missing ``-R`` is added as ``M(R)^+`` (codes that store one direction
    of each bond). Present pairs are replaced by their mean
    ``(M(R) + M(-R)^+)/2``; ``residual`` is the largest
    ``|M(R) - M(-R)^+| / (2 max|M|)`` before the mean, which measures how far
    a predicted Hamiltonian is from Hermitian. The result is exactly Hermitian
    at every k, as the Bloch sums of LSQuant require.
    """
    r_int = np.asarray(r_int, int)
    mats = np.asarray(mats, complex)
    lookup = {tuple(r): n for n, r in enumerate(r_int)}
    r_out, m_out = [tuple(r) for r in r_int], list(mats)
    for r in list(lookup):
        neg = tuple(-x for x in r)
        if neg not in lookup:
            lookup[neg] = len(r_out)
            r_out.append(neg)
            m_out.append(mats[lookup[r]].conj().T)
    m_out = np.array(m_out)
    partner = np.array([lookup[tuple(-x for x in r)] for r in r_out])
    mirrored = m_out[partner].conj().swapaxes(1, 2)
    scale = max(float(np.abs(m_out).max()), 1e-300)
    residual = float(np.abs(m_out - mirrored).max()) / (2 * scale)
    return np.array(r_out, int), 0.5 * (m_out + mirrored), residual


def check_spin_scalar(r_int, overlap, name="S"):
    """Raise unless an interleaved spinful ``S(R)`` equals ``S_0 tensor I_2``."""
    up, down = overlap[:, 0::2, 0::2], overlap[:, 1::2, 1::2]
    scale = max(float(np.abs(overlap).max()), 1e-300)
    error = max(float(np.abs(up - down).max()), float(np.abs(overlap[:, 0::2, 1::2]).max()),
                float(np.abs(overlap[:, 1::2, 0::2]).max())) / scale
    if error > 1e-8:
        raise ValueError(f"{name}: the spinful overlap is not S_0 tensor I_2 (relative error {error:.1e})")


def build_lcao_model(uc, positions, norb, r_int, hamiltonian, overlap=None, *, spinful,
                     seed_name, metadata, fermi_energy=0.0, overlap_r_int=None,
                     asymmetry_warning=1e-3):
    """Return the LSQuant model of an atom-pair Hamiltonian.

    ``hamiltonian`` and ``overlap`` are the outputs of
    :func:`assemble_blocks` (interleaved spinors when ``spinful``);
    ``overlap`` may be spinless for a spinful model, in which case it is
    lifted to ``S tensor I_2``. Both are made Hermitian with
    :func:`hermitian_pairs`; the residuals are stored in ``metadata``, and a
    Hamiltonian residual above ``asymmetry_warning`` raises a
    ``RuntimeWarning``. ``fermi_energy`` (eV) is removed as
    ``H(R) - E_F S(R)``, so run-file energies are measured from it.

    Without ``overlap`` the model is an
    :class:`~lsquant.models.orthogonal.OrthogonalModel`; with it, a
    :class:`~lsquant.models.siesta.SiestaModel` (the non-orthogonal model of
    the metric executors). The attribute ``spinful`` tells
    :class:`~lsquant.bases.BlochMetricBasis` whether spin operators exist.
    """
    import warnings

    from .orthogonal import OrthogonalModel
    from .siesta import SiestaModel

    uc = np.asarray(uc, float)
    norb = np.asarray(norb, int)
    positions = np.asarray(positions, float)
    if uc.shape != (3, 3) or positions.shape != (len(norb), 3):
        raise ValueError("lattice must be 3x3 and positions one row per atom")
    tau = np.repeat(positions, norb * (2 if spinful else 1), axis=0)
    r_int, hamiltonian, h_residual = hermitian_pairs(r_int, hamiltonian)
    meta = dict(metadata)
    meta.update(hamiltonian_hermiticity_residual=h_residual, fermi_energy_subtracted_eV=float(fermi_energy),
                spinful=bool(spinful), atoms=int(len(norb)), orbitals_per_atom=[int(n) for n in norb])
    if h_residual > asymmetry_warning:
        warnings.warn(f"H(R) differs from H(-R)^+ by {h_residual:.1e} of max|H|; LSQuant uses the mean "
                      "of the two (check the prediction)", RuntimeWarning, stacklevel=3)
    if overlap is None:
        hamiltonian = hamiltonian.copy()
        origin = np.where(~np.any(r_int, axis=1))[0]
        if len(origin) != 1:
            raise ValueError("the Hamiltonian has no R = 0 block")
        hamiltonian[origin[0]] -= float(fermi_energy) * np.eye(hamiltonian.shape[1])
        model = OrthogonalModel(uc, tau, r_int, hamiltonian, seed_name=seed_name, metadata=meta)
        if not spinful:
            model.metadata["spin_operator"] = "none"
        model.spinful = bool(spinful)
        return model
    s_r_int = r_int if overlap_r_int is None else np.asarray(overlap_r_int, int)
    s_r_int, overlap, s_residual = hermitian_pairs(s_r_int, overlap)
    if spinful and overlap.shape[1] * 2 == hamiltonian.shape[1]:
        overlap = spinor_overlap(overlap)
    if overlap.shape[1:] != hamiltonian.shape[1:]:
        raise ValueError(f"S(R) blocks {overlap.shape[1:]} do not match H(R) blocks {hamiltonian.shape[1:]}")
    if spinful:
        check_spin_scalar(s_r_int, overlap)
    # One R list for H and S: H on S-only cells is zero and S on H-only cells is zero.
    cells = {tuple(r): n for n, r in enumerate(r_int)}
    for r in s_r_int:
        if tuple(r) not in cells:
            cells[tuple(r)] = len(cells)
    width = hamiltonian.shape[1]
    H = np.zeros((len(cells), width, width), complex)
    S = np.zeros((len(cells), width, width), complex)
    H[[cells[tuple(r)] for r in r_int]] = hamiltonian
    S[[cells[tuple(r)] for r in s_r_int]] = overlap
    H -= float(fermi_energy) * S
    meta["overlap_hermiticity_residual"] = s_residual
    model = SiestaModel(uc, tau, np.array(list(cells), int), H, S, seed_name=seed_name,
                        metadata=meta, input_tolerance=1e-8)
    # SiestaModel records the SIESTA energy reference; this model records its own.
    model.metadata.pop("stored_fermi_level", None)
    model.metadata["orthogonal"] = False
    model.spinful = bool(spinful)
    return model


def cell_gauge_hk(model, k_frac):
    """``sum_R exp(+2 pi i k.R) H(R)`` at fractional k, for tests and checks."""
    k_cart = np.atleast_2d(k_frac) @ (2 * np.pi * np.linalg.inv(model.uc).T)
    return bloch_sum(model.r_cart, model.tau, model.hoppings, k_cart, cell=True)
