"""ABACUS-format Hamiltonians: NextHAM predictions (and ABACUS's own matrices).

NextHAM (Yin, Dai, Pan and He, ICLR 2026) is an equivariant network that
predicts the Hamiltonian of ABACUS in its numerical atomic orbitals, with
spin-orbit coupling. Its post-processing (``NextHAM2Bands``:
``output_hs_data.py`` writes ``hr_cut``, ``add_element.py`` completes it
beyond the cut-off with the ABACUS matrix of the same structure and writes
``hr_tot``) stores the predicted H(R) in the sparse format of ABACUS,
the format of ABACUS's own ``hrs1_nao.csr`` / ``data-HR-sparse_SPIN0.csr``.
The overlap is ABACUS's ``srs1_nao.csr`` and the structure its ``STRU``.
LSQuant reads these three files; ABACUS's own H(R) file can be read the
same way, which compares the prediction with DFT.

Files and conventions, read in the sources of NextHAM
(``NextHAM2Bands/add_element.py``, ``output_hs_data.py``) and of pyatb,
the band code NextHAM uses (``src/pyatb/io/abacus_read_xr.py``,
``abacus_read_stru.py``, ``src/cpp/core/base_data.cpp``):

- Sparse matrix file (text): ``STEP: 0``, ``Matrix Dimension of H(R): N``,
  ``Matrix number of H(R): M``; then per R a line ``R1 R2 R3 nnz`` followed,
  if ``nnz > 0``, by the CSR ``data`` (real numbers, or ``(re,im)`` pairs
  for ``nspin = 4``), the 0-based column ``indices`` and the ``indptr``
  (N + 1 numbers). The binary variant (``STEP`` missing) holds the header
  ``int32 step, N, M``, per R ``int32 R1, R2, R3, nnz``, then ``nnz`` float64
  or complex128 values, ``nnz`` int32 indices and N + 1 int32 ``indptr``.
- Rows are orbitals of the home cell, columns orbitals of cell R, and
  ``X(k) = sum_R exp(+2 pi i k.R) X(R)`` (pyatb ``get_exp_ikR``): the
  LSQuant convention ``<i, 0|X|j, R>``, no sign change.
- Orbital order: atoms species by species in the order of
  ``ATOMIC_POSITIONS``, the orbitals of one atom contiguous (by l, then
  zeta, then m). With spin-orbit coupling (``nspin = 4``) the spin index runs
  fastest (orbital up, orbital down), which is the interleaved layout of
  LSQuant (NextHAM transposes its ``[[uu, ud], [du, dd]]`` blocks into it).
- Orbitals per atom: from the numerical orbital file of each species named
  in ``STRU`` (``NUMERICAL_ORBITAL``): line 5 ``Lmax``, then one line per l
  ending with the number of radial functions; ``n_l (2l + 1)`` orbitals.
  These files are the basis of the training data.
- Units: H in Rydberg (ABACUS; NextHAM divides its eV by 13.605698066), S
  dimensionless, ``LATTICE_CONSTANT`` in Bohr. Converted to eV and Angstrom
  with pyatb's constants.
- Energy reference: the absolute scale of ABACUS; give the Fermi energy of
  the ABACUS run (``E_Fermi`` in its log, eV) as ``fermi_energy``.

Pitfall: NextHAM predicts H only within its cut-off (7 Angstrom in
``add_element.py``); ``hr_tot`` takes the matrix elements beyond it from the
ABACUS file given as ``weak_file`` in its pipeline. The model is therefore
partly predicted and partly DFT.

Position matrix (``out_mat_r 1``): ``data-rR-sparse.csr`` holds
``Q_a,ij(R) = <i, 0| r_a |j, R>`` in Bohr, the absolute position (ABACUS adds
``tau_i S_ij(R)`` to the integral of ``r - tau_i``; ``cal_r_overlap_R.cpp``).
Layout: the same three header lines (``r(R)``); per R a line ``R1 R2 R3``,
then for x, y and z a line ``nnz`` followed, if ``nnz > 0``, by the CSR
``data``, ``indices`` and ``indptr`` (real numbers). For ``nspin = 4`` the
dimension is that of H, with zero spin off-diagonal elements. With
``position_file`` the model gets the connection
``A = Q - (tau_i + tau_j + R) S / 2`` of the covariant velocity (see
:func:`~lsquant.models.siesta.midpoint_connection`); exact two-centre
integrals give ``A(-R)^+ = A(R)``. Without it the orbitals sit on their
atoms with no intra-atomic dipoles (``bare`` and ``berry_connection`` only).
Checked on PdSe2 (ABACUS 3.10.1, PseudoDojo NC FR, Dojo-NC-FR DZP
orbitals, 300 spinor orbitals): the eigenvalues of the model equal those
ABACUS prints to 4e-4 meV at seven k points of the SCF mesh; the
anti-Hermitian part of the connection is 1.6e-8 of its largest element; the
exact covariant spin Hall plateau is 1.948 e^2/h (48x48 to 72x72), against
1.979 (SIESTA) and 1.971 (OpenMX), and it is 0.005 in the trivial gap, where
``berry_connection`` gives -0.23.
The orbital order and the sign convention of ABACUS's real spherical
harmonics do not matter here: H, S and r are in the same basis, and the spin
operator is the Pauli matrix on the interleaved spin index.
"""
from __future__ import annotations

import os
import re
import struct

import numpy as np

#: Rydberg in eV and Angstrom in Bohr, as in pyatb's constants.py; kept in
#: lsquant.models.units under their own names (the Rydberg differs from the
#: ThreeBodyTB value and from HARTREE_EV / 2 on purpose).
from .units import ANGSTROM_BOHR
from .units import RYDBERG_EV_ABACUS as RYDBERG_EV
_POSITION_KINDS = ("Direct", "Cartesian", "Cartesian_au", "Cartesian_angstrom")


def read_abacus_csr(path):
    """Return ``(r_int, mats)`` of an ABACUS sparse matrix file, text or binary.

    ``mats`` has shape ``(M, N, N)``: complex for ``nspin = 4`` files, real
    otherwise.
    """
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:4] == b"STEP":
        return _read_text(raw.decode(), path)
    return _read_binary(raw, path)


def _read_text(text, path):
    lines = iter(text.splitlines())
    next(lines)
    n = int(next(lines).split()[-1])
    m = int(next(lines).split()[-1])
    r_int = np.zeros((m, 3), int)
    mats = None
    for r in range(m):
        head = next(lines).split()
        r_int[r] = [int(float(x)) for x in head[:3]]
        nnz = int(head[3])
        if mats is None:
            mats = np.zeros((m, n, n), complex)
        if nnz == 0:
            continue
        data_line = next(lines)
        if "(" in data_line:
            pairs = re.findall(r"\(([^)]*)\)", data_line)
            data = np.array([complex(float(a), float(b)) for a, b in (p.split(",") for p in pairs)])
        else:
            data = np.array(data_line.split(), float)
        indices = np.array(next(lines).split(), int)
        indptr = np.array(next(lines).split(), int)
        if len(data) != nnz or len(indices) != nnz or len(indptr) != n + 1:
            raise ValueError(f"{path}: R = {r_int[r]} has inconsistent CSR sizes")
        rows = np.repeat(np.arange(n), np.diff(indptr))
        mats[r, rows, indices] = data
    if mats is None:
        mats = np.zeros((0, n, n), complex)
    return r_int, mats


def _read_binary(raw, path):
    step, n, m = struct.unpack_from("iii", raw, 0)

    def parse(width):
        offset, r_int, mats = 12, np.zeros((m, 3), int), np.zeros((m, n, n), complex)
        for r in range(m):
            r1, r2, r3, nnz = struct.unpack_from("iiii", raw, offset)
            offset += 16
            r_int[r] = (r1, r2, r3)
            if nnz:
                dtype = complex if width == 16 else float
                data = np.frombuffer(raw, dtype=dtype, count=nnz, offset=offset)
                offset += nnz * width
                indices = np.frombuffer(raw, dtype=np.int32, count=nnz, offset=offset)
                offset += nnz * 4
            indptr = np.frombuffer(raw, dtype=np.int32, count=n + 1, offset=offset)
            offset += (n + 1) * 4
            if nnz:
                rows = np.repeat(np.arange(n), np.diff(indptr))
                mats[r, rows, indices] = data
        return offset, r_int, mats

    for width in (8, 16):
        try:
            offset, r_int, mats = parse(width)
        except (struct.error, ValueError):
            continue
        if offset == len(raw):
            return r_int, mats
    raise ValueError(f"{path} is neither an ABACUS text nor binary sparse matrix file")


def read_abacus_rr(path):
    """Return ``(r_int, Q)`` of an ABACUS ``data-rR-sparse.csr``, text or binary.

    ``Q`` has shape ``(3, M, N, N)``, the position matrix in Bohr.
    """
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:4] == b"STEP":
        lines = iter(raw.decode().splitlines())
        next(lines)
        n = int(next(lines).split()[-1])
        m = int(next(lines).split()[-1])
        r_int, q = np.zeros((m, 3), int), np.zeros((3, m, n, n))
        for r in range(m):
            r_int[r] = [int(float(x)) for x in next(lines).split()[:3]]
            for a in range(3):
                nnz = int(next(lines).split()[0])
                if nnz == 0:
                    continue
                data = np.array(next(lines).split(), float)
                indices = np.array(next(lines).split(), int)
                indptr = np.array(next(lines).split(), int)
                if len(data) != nnz or len(indices) != nnz or len(indptr) != n + 1:
                    raise ValueError(f"{path}: R = {r_int[r]}, direction {a} has inconsistent CSR sizes")
                q[a, r, np.repeat(np.arange(n), np.diff(indptr)), indices] = data
        return r_int, q
    _, n, m = struct.unpack_from("iii", raw, 0)
    offset, r_int, q = 12, np.zeros((m, 3), int), np.zeros((3, m, n, n))
    for r in range(m):
        r_int[r] = struct.unpack_from("iii", raw, offset)
        offset += 12
        for a in range(3):
            (nnz,) = struct.unpack_from("i", raw, offset)
            offset += 4
            if nnz == 0:
                continue
            data = np.frombuffer(raw, float, nnz, offset)
            offset += 8 * nnz
            indices = np.frombuffer(raw, np.int32, nnz, offset)
            offset += 4 * nnz
            indptr = np.frombuffer(raw, np.int32, n + 1, offset)
            offset += 4 * (n + 1)
            q[a, r, np.repeat(np.arange(n), np.diff(indptr)), indices] = data
    if offset != len(raw):
        raise ValueError(f"{path} is neither a text nor a binary ABACUS r(R) file")
    return r_int, q


def _strip(line):
    return re.sub(r"(#|//).*", "", line).strip()


def read_stru(path):
    """Return ``(lattice_angstrom, species, positions_angstrom)`` of an ABACUS STRU file.

    ``species`` lists ``(label, orbital_file, n_atoms)`` in file order;
    ``positions`` has one row per atom in that order.
    """
    lines = [x for x in (_strip(line) for line in open(path, encoding="utf-8")) if x]
    keys = {k: lines.index(k) for k in ("ATOMIC_SPECIES", "NUMERICAL_ORBITAL", "LATTICE_CONSTANT",
                                         "LATTICE_VECTORS", "ATOMIC_POSITIONS") if k in lines}
    for k in ("ATOMIC_SPECIES", "LATTICE_CONSTANT", "LATTICE_VECTORS", "ATOMIC_POSITIONS"):
        if k not in keys:
            raise ValueError(f"{path}: no {k} block (a STRU with latname is not supported)")
    labels = []
    i = keys["ATOMIC_SPECIES"] + 1
    while i < len(lines) and lines[i] not in keys:
        labels.append(lines[i].split()[0])
        i += 1
    orbital_files = [None] * len(labels)
    if "NUMERICAL_ORBITAL" in keys:
        start = keys["NUMERICAL_ORBITAL"] + 1
        orbital_files = [lines[start + s].split()[0] for s in range(len(labels))]
    alat = float(lines[keys["LATTICE_CONSTANT"] + 1].split()[0]) / ANGSTROM_BOHR
    vectors = np.array([lines[keys["LATTICE_VECTORS"] + 1 + a].split()[:3] for a in range(3)], float)
    lattice = vectors * alat
    k = keys["ATOMIC_POSITIONS"] + 1
    kind = lines[k].split()[0]
    if kind not in _POSITION_KINDS:
        raise ValueError(f"{path}: ATOMIC_POSITIONS {kind!r} is not supported; use one of {_POSITION_KINDS}")
    k += 1
    species, positions = [], []
    for s, label in enumerate(labels):
        if lines[k].split()[0] != label:
            raise ValueError(f"{path}: ATOMIC_POSITIONS lists {lines[k].split()[0]!r} where {label!r} is expected")
        count = int(lines[k + 2].split()[0])
        coords = np.array([lines[k + 3 + a].split()[:3] for a in range(count)], float)
        k += 3 + count
        if kind == "Direct":
            coords = coords @ lattice
        elif kind == "Cartesian":
            coords = coords * alat
        elif kind == "Cartesian_au":
            coords = coords / ANGSTROM_BOHR
        species.append((label, orbital_files[s], count))
        positions.append(coords)
    return lattice, species, np.concatenate(positions)


def read_orbital_shells(path):
    """Number of radial functions per l (``[n_s, n_p, ...]``) of an ABACUS numerical orbital file."""
    with open(path, encoding="utf-8") as f:
        lines = [f.readline() for _ in range(5)]
        if "lmax" not in lines[4].lower():
            raise ValueError(f"{path}: line 5 is not 'Lmax' (not an ABACUS orbital file)")
        lmax = int(lines[4].split()[-1])
        return [int(f.readline().split()[-1]) for _ in range(lmax + 1)]


def _element(label):
    return re.match(r"[A-Z][a-z]?", label).group(0)


def _shells_of(label, orbital_file, search_dirs, orbitals_per_element):
    """``(n_orbitals, l of every orbital or None)`` of one atom of species ``label``."""
    if orbitals_per_element:
        for key in (label, _element(label)):
            if key in orbitals_per_element:
                return int(orbitals_per_element[key]), None
    if orbital_file:
        for d in search_dirs:
            p = os.path.join(d, orbital_file)
            if os.path.isfile(p):
                shells = read_orbital_shells(p)
                ell = [l for l, n in enumerate(shells) for _ in range(n * (2 * l + 1))]
                return len(ell), ell
    raise ValueError(f"no orbital file for species {label!r} ({orbital_file}) in {search_dirs}; set "
                     "params.orbital_dir or params.orbitals_per_element")


def load_abacus_model(hamiltonian_file, *, overlap_file, stru_file, orbital_dir=None,
                      orbitals_per_element=None, fermi_energy=0.0, position_file=None):
    """Read an ABACUS-format H(R) (NextHAM ``hr_tot``, or ABACUS's own) with its S(R) and STRU.

    Returns a non-orthogonal model (:class:`~lsquant.models.siesta.SiestaModel`)
    in eV and Angstrom, spinful when the matrix dimension is twice the
    number of orbitals (``nspin = 4``). ``fermi_energy`` (eV) is removed as
    ``H - E_F S``. The orbital files are searched in ``orbital_dir``, then
    next to ``stru_file``. ``position_file`` (ABACUS's ``data-rR-sparse.csr``)
    adds the connection of the covariant velocity.
    """
    from .lcao import build_lcao_model

    paths = [os.path.abspath(os.path.expanduser(str(p))) for p in (hamiltonian_file, overlap_file, stru_file)]
    h_path, s_path, stru_path = paths
    lattice, species, positions = read_stru(stru_path)
    search = [os.path.abspath(os.path.expanduser(str(orbital_dir)))] if orbital_dir else []
    search.append(os.path.dirname(stru_path))
    norb, ell, elements = [], [], []
    for label, orbital_file, count in species:
        n, l_list = _shells_of(label, orbital_file, search, orbitals_per_element)
        norb += [n] * count
        ell += (l_list or [-1] * n) * count
        elements += [label] * count
    norb = np.array(norb, int)
    h_r, hamiltonian = read_abacus_csr(h_path)
    s_r, overlap = read_abacus_csr(s_path)
    width = hamiltonian.shape[1]
    if width == norb.sum():
        spinful = False
    elif width == 2 * norb.sum():
        spinful = True
    else:
        raise ValueError(f"H(R) has dimension {width}, but STRU and the orbital files give {norb.sum()} "
                         "orbitals; check the orbital files (the basis of the training data)")
    meta = dict(source="abacus", hamiltonian_file=h_path, overlap_file=s_path, stru_file=stru_path,
                elements=elements, orbital_l=[int(x) for x in ell],
                basis="ABACUS numerical atomic orbitals of STRU (NUMERICAL_ORBITAL)",
                orbital_positions="atom-centred: every orbital on its atom")
    model = build_lcao_model(lattice, positions, norb, h_r, hamiltonian * RYDBERG_EV, overlap, spinful=spinful,
                             seed_name=os.path.basename(h_path), metadata=meta,
                             fermi_energy=float(fermi_energy), overlap_r_int=s_r)
    if position_file is not None:
        q_path = os.path.abspath(os.path.expanduser(str(position_file)))
        _attach_connection(model, q_path)
    return model


def _attach_connection(model, q_path):
    """Give ``model`` the midpoint connection built from ABACUS's r(R) file."""
    from .siesta import fold_opposite, midpoint_connection, opposite_pairs, spinor_lift

    q_r, q = read_abacus_rr(q_path)
    width = model.num_wann
    if q.shape[2] * 2 == width:
        q = spinor_lift(q)
    elif q.shape[2] != width:
        raise ValueError(f"r(R) has dimension {q.shape[2]}, H(R) {width}; use the files of one ABACUS run")
    cells = {tuple(r): n for n, r in enumerate(model.r_int)}
    full = np.zeros((3, len(cells), width, width))
    missing = 0.0
    for n, r in enumerate(q_r):
        if tuple(r) in cells:
            full[:, cells[tuple(r)]] = q[:, n]
        else:
            missing = max(missing, float(np.abs(q[:, n]).max()))
    if missing > 1e-6:
        raise ValueError(f"r(R) has elements up to {missing:.1e} Bohr on cells where S(R) is zero; "
                         "use the files of one ABACUS run")
    a = midpoint_connection(full / ANGSTROM_BOHR, model.overlaps, model.tau, model.r_cart)
    sym, residual = fold_opposite(a, opposite_pairs(model.r_int), residual="mean")
    model.metadata["connection_hermiticity_residual"] = residual
    model.metadata["position_file"] = q_path
    model.dipoles = sym
