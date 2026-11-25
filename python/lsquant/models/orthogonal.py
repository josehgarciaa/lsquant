"""Orthogonal models for the metric executor: Wannier90 with exact spin, and PAOFLOW.

The run-file executor of SIESTA builds every operator through
:class:`~lsquant.bases.BlochMetricBasis`, whose model protocol is ``Hk``,
``Sk``, ``derivatives``, ``Ak``, ``k_mesh``, ``tau`` and ``uc``. The classes
here give an orthogonal model (``S = 1``) that protocol, so a Wannier90 or a
PAOFLOW Hamiltonian reaches the same velocity modes, physical spin operators,
spectral projection, exact oracle and KPM moments as a SIESTA model.

Gauge. ``H(k)_ab = sum_R exp(+i k.(R + tau_b - tau_a)) H_ab(R)`` (atom gauge,
:func:`~lsquant.models.wannier.hk_blocks`). Its derivative is
``dH_atom = dH_cell + i[H, tau]``, so the covariant velocity of Wannier90,
``dH_cell + i[H, A_full]``, is ``dH_atom + i[H, A_full - tau]``: the
connection handed to :func:`~lsquant.models.nonorthogonal.bloch_velocity` is
``A = A_full - diag(tau)``. Velocities use the ``+dH`` sign of SIESTA; the
Wannier90 method :meth:`WannierModel.velocity_blocks` returns ``-dH``. A Hall
or spin Hall word holds two velocities, so the two signs give the same result.

Spin. With spin-orbit coupling the Wannier functions mix spin, so the
diagonal Pauli matrix on the spinor labels is not the spin operator.
:func:`read_wannier_spin` builds ``S_W(k) = V^+ S_B V`` from the DFT spin
matrices (``.spn``) with ``V = U_dis U``; the rows of ``_u_dis.mat`` list
only the bands inside the outer window, packed at the top.
"""
from __future__ import annotations

import os
import re

import numpy as np

from .bloch import bloch_sum, k_mesh
from .registry import load_runfile_model  # noqa: F401  (re-export; the dispatcher lives in the connector registry)
from .units import BOHR_ANGSTROM
from .wannier import hk_blocks, load_wannier_model, read_position_matrix

_AXES = {"x": 0, "y": 1, "z": 2}


def _k_chunk_count(n_points, fixed_bytes, per_k_bytes, memory_mb):
    """k points per Bloch-sum chunk under a MiB workspace bound.

    The sums are independent per k, so chunking the k points is
    bit-identical to one call; the bound only sets the chunk size.
    """
    chunk = max(1, int((float(memory_mb) * 1024**2 - fixed_bytes) // per_k_bytes))
    return min(chunk, n_points)


def _hk_blocks_bounded(r_cart, tau, hoppings, k_points, memory_mb):
    """``hk_blocks`` over k chunks bounded by ``memory_mb`` (bit-identical).

    Pitfall: ``hk_blocks`` and ``bloch_sum`` dress the atom gauge in a
    different order, so the bound is applied by chunking the k points here
    instead of switching routines.
    """
    points = np.atleast_2d(np.asarray(k_points, dtype=float))
    nr, width = hoppings.shape[0], hoppings.shape[1]
    chunk = _k_chunk_count(len(points), 2 * hoppings.nbytes,
                           16 * nr + 32 * width + 16 * width * width, memory_mb)
    if chunk >= len(points):
        return hk_blocks(r_cart, tau, hoppings, points)
    out = np.empty((len(points), width, width), dtype=np.complex128)
    for start in range(0, len(points), chunk):
        out[start:start + chunk] = hk_blocks(r_cart, tau, hoppings, points[start:start + chunk])
    return out


def _read_win_window(path):
    """Return the outer window (eV) of a ``.win`` file; missing keys mean all bands."""
    lo, hi = -np.inf, np.inf
    if not os.path.isfile(path):
        return lo, hi
    for line in open(path, encoding="utf-8"):
        body = re.split(r"[!#]", line)[0].strip().lower().replace(":", "=")
        if "=" not in body:
            continue
        key, value = (s.strip() for s in body.split("=", 1))
        if key == "dis_win_min":
            lo = float(value.split()[0])
        elif key == "dis_win_max":
            hi = float(value.split()[0])
    return lo, hi


def _read_spn(path):
    """Read ``<psi_mk|sigma_a|psi_nk>`` (Pauli) from a formatted or unformatted ``.spn``.

    Returns an array of shape ``(3, N_k, N_b, N_b)``. Both files store the
    upper triangle ``m <= n`` with ``n`` in the outer loop and three
    components per element.
    """
    try:
        from scipy.io import FortranFile
        with FortranFile(path, "r") as f:
            f.read_record("S60")
            nb, nk = (int(x) for x in f.read_ints(np.int32))
            records = [f.read_reals(np.float64).view(np.complex128).reshape(-1, 3)
                       for _ in range(nk)]
    except (ValueError, OSError, TypeError):
        lines = open(path, encoding="utf-8").read().split("\n")
        nb, nk = (int(x) for x in lines[1].split()[:2])
        values = np.array([[float(x) for x in ln.split()[:2]] for ln in lines[2:] if ln.strip()])
        flat = values[:, 0] + 1j * values[:, 1]
        per_k = 3 * nb * (nb + 1) // 2
        records = [flat[k * per_k:(k + 1) * per_k].reshape(-1, 3) for k in range(nk)]
    S = np.zeros((3, nk, nb, nb), complex)
    pairs = np.array([(m, n) for n in range(nb) for m in range(n + 1)])
    m_idx, n_idx = pairs[:, 0], pairs[:, 1]
    for k, rec in enumerate(records):
        for a in range(3):
            S[a, k, m_idx, n_idx] = rec[:, a]
            S[a, k, n_idx, m_idx] = np.conj(rec[:, a])
    return S


def _read_umat(path, nk):
    """Read ``_u.mat`` or ``_u_dis.mat``: fractional k points and ``(N_k, rows, cols)``."""
    lines = open(path, encoding="utf-8").read().splitlines()
    nk_file, n1, n2 = (int(x) for x in lines[1].split()[:3])
    if nk_file != nk:
        raise ValueError(f"{path}: {nk_file} k points, expected {nk}")
    # _u.mat stores num_wann x num_wann; _u_dis.mat stores num_bands rows, num_wann columns.
    rows, cols = (n2, n1) if n2 != n1 else (n1, n2)
    i, kpts, mats = 2, [], []
    for _ in range(nk):
        while not lines[i].strip():
            i += 1
        kpts.append([float(x) for x in lines[i].split()[:3]])
        i += 1
        v = np.array([[float(x) for x in lines[i + j].split()[:2]] for j in range(rows * cols)])
        i += rows * cols
        mats.append((v[:, 0] + 1j * v[:, 1]).reshape(rows, cols, order="F"))
    return np.array(kpts), np.array(mats)


def read_wannier_spin(model_dir, seed, *, unitarity_tolerance=1e-6):
    """Return ``(k_frac, S_W)`` with ``S_W`` of shape ``(3, N_k, W, W)`` (Pauli).

    Reads ``<seed>.spn``, ``<seed>_u.mat`` and, for a disentangled model,
    ``<seed>_u_dis.mat``, ``<seed>.eig`` and the outer window of
    ``<seed>.win``. Raises ``ValueError`` when ``V^+ V`` is not the identity,
    which signals a wrong window or band count.
    """
    base = os.path.join(model_dir, seed)
    S_b = _read_spn(base + ".spn")
    _, nk, nb, _ = S_b.shape
    probe = open(base + "_u.mat", encoding="utf-8").read().splitlines()[1].split()
    nw = int(probe[1])
    kpts, U = _read_umat(base + "_u.mat", nk)
    if os.path.isfile(base + "_u_dis.mat"):
        _, Ud = _read_umat(base + "_u_dis.mat", nk)
        eig = np.loadtxt(base + ".eig")
        if eig.shape[0] != nk * nb:
            raise ValueError(f"{base}.eig has {eig.shape[0]} rows, expected {nk * nb}")
        E = eig[:, 2].reshape(nk, nb)
        lo, hi = _read_win_window(base + ".win")
        V = np.zeros((nk, nb, nw), complex)
        for k in range(nk):
            inside = np.where((E[k] >= lo) & (E[k] <= hi))[0]
            V[k, inside, :] = Ud[k, :len(inside), :] @ U[k]
    else:
        if nb != nw:
            raise ValueError(f"{base}: {nb} bands and {nw} Wannier functions but no _u_dis.mat")
        V = U
    error = float(np.abs(np.einsum("kba,kbc->kac", V.conj(), V) - np.eye(nw)).max())
    if error > unitarity_tolerance:
        raise ValueError(f"V^+V differs from the identity by {error:.1e}; check the outer window")
    S_w = np.einsum("kba,xkbc,kcd->xkad", V.conj(), S_b, V)
    return kpts, S_w


class WannierSpin:
    """Exact Wannier spin interpolated from the Wannier90 k mesh.

    ``S(R) = (1/N_k) sum_k exp(-2 pi i k.R) S_W(k)`` on the R range of the
    ``n1 x n2 x n3`` mesh, centred on the origin. :meth:`blocks` returns the
    atom-gauge matrix at Cartesian k points. ``bloch_memory_mb`` bounds the
    workspace of the R sum, chunking the k points (bit-identical per k).
    """

    bloch_memory_mb = 256.0

    def __init__(self, model_dir, seed, uc, tau):
        kpts, S_w = read_wannier_spin(model_dir, seed)
        n = [len(np.unique(np.round(kpts[:, a] % 1.0, 8))) for a in range(3)]
        # Centred range: -(m//2) .. m - m//2 - 1, so a single point gives R = 0 only.
        ranges = [range(-(m // 2), m - m // 2) for m in n]
        self.r_int = np.array([[a, b, c] for a in ranges[0] for b in ranges[1] for c in ranges[2]])
        phase = np.exp(-2j * np.pi * (kpts @ self.r_int.T))
        self.spin_r = np.einsum("kr,xkab->xrab", phase, S_w) / len(kpts)
        self.uc, self.tau = np.asarray(uc, float), np.asarray(tau, float)
        self.mesh = tuple(n)

    def blocks(self, k_points, component):
        a = _AXES[component.lower()]
        k = np.atleast_2d(k_points)
        kf = k @ self.uc.T / (2 * np.pi)
        spin = self.spin_r[a]
        nr, width = len(self.r_int), len(self.tau)
        chunk = _k_chunk_count(len(k), 2 * spin.nbytes,
                               16 * nr + 16 * width * width, self.bloch_memory_mb)
        cell = np.empty((len(k), width, width), dtype=np.complex128)
        phases = np.exp(2j * np.pi * (kf @ self.r_int.T))
        for start in range(0, len(k), chunk):
            stop = start + chunk
            cell[start:stop] = np.einsum("kr,rab->kab", phases[start:stop], spin, optimize=True)
        p = np.exp(1j * (k @ self.tau.T))
        atom = p.conj()[:, :, None] * cell * p[:, None, :]
        return 0.5 * (atom + np.conj(np.swapaxes(atom, 1, 2)))


class OrthogonalModel:
    """An orthogonal tight-binding model in the protocol of BlochMetricBasis.

    ``hoppings`` follow the atom-gauge convention of
    :func:`~lsquant.models.wannier.hk_blocks`. ``positions`` is the optional
    full position matrix ``(r_int, A)`` of :func:`read_position_matrix`
    (enables ``covariant``). ``spin`` is ``None`` (diagonal Pauli matrices on
    interleaved spinor labels), an object with ``blocks(k, component)``, or a
    ``(3, W, W)`` array of a k-independent on-site spin operator.
    """

    input_tolerance = 1e-10
    bloch_memory_mb = 256.0
    physical_spin = None

    def __init__(self, uc, tau, r_int, hoppings, *, seed_name="model", positions=None,
                 spin=None, metadata=None):
        self.uc = np.asarray(uc, float)
        self.tau = np.asarray(tau, float)
        self.r_int = np.asarray(r_int, int)
        self.r_cart = self.r_int @ self.uc
        self.hoppings = np.asarray(hoppings, complex)
        self.num_wann = len(self.tau)
        if self.hoppings.shape != (len(self.r_int), self.num_wann, self.num_wann):
            raise ValueError("hoppings must have shape (N_R, W, W)")
        self.seed_name = seed_name
        self.positions = positions
        self.dipoles = None if positions is None else positions[1]
        self.spin = spin
        self.metadata = dict(metadata or {})
        self.metadata.update(spin_ordering="interleaved", gauge="atom", energy_unit="eV",
                             length_unit="Angstrom", orthogonal=True,
                             spin_operator=("labels" if spin is None else
                                            "on-site" if isinstance(spin, np.ndarray) else "exact"))

    def k_mesh(self, supercell):
        return k_mesh(self.uc, supercell)

    def Hk(self, k):
        return _hk_blocks_bounded(self.r_cart, self.tau, self.hoppings,
                                  np.atleast_2d(k), self.bloch_memory_mb)

    def Sk(self, k):
        n = len(np.atleast_2d(k))
        return np.broadcast_to(np.eye(self.num_wann, dtype=complex), (n, self.num_wann, self.num_wann)).copy()

    def derivatives(self, k, axis, *, cell=False):
        """Return ``(dH_a, dS_a=0)``: atom derivative, or cell derivative in the atom gauge."""
        k = np.atleast_2d(k)
        if cell:
            dH = bloch_sum(self.r_cart, np.zeros_like(self.tau), self.hoppings, k,
                           derivative=axis, cell=True, memory_mb=self.bloch_memory_mb)
            p = np.exp(1j * (k @ self.tau.T))
            dH = p.conj()[:, :, None] * dH * p[:, None, :]
        else:
            # +dH_atom/dk. dhk_blocks keeps the wannier2sparse sign, -dH_atom/dk.
            dH = bloch_sum(self.r_cart, self.tau, self.hoppings, k, derivative=axis,
                           memory_mb=self.bloch_memory_mb)
        return dH, np.zeros_like(dH)

    def Ak(self, k, axis):
        """Return ``A_full(k) - diag(tau_a)`` in the atom gauge (Angstrom)."""
        if self.positions is None:
            raise ValueError("covariant velocity requires the position matrix (_r.dat)")
        r_int, A = self.positions
        k = np.atleast_2d(k)
        full = bloch_sum(r_int @ self.uc, self.tau, A[axis], k,
                         memory_mb=self.bloch_memory_mb)
        return full - np.diag(self.tau[:, axis])[None, :, :]

    def spin_blocks(self, k, component):
        """Return the physical Pauli spin matrix, or None for spinor labels."""
        if self.spin is None:
            return None
        k = np.atleast_2d(k)
        if isinstance(self.spin, np.ndarray):
            block = self.spin[_AXES[component.lower()]]
            return np.broadcast_to(block, (len(k),) + block.shape).copy()
        return self.spin.blocks(k, component)


def load_wannier_metric_model(model_dir, seed, *, spin="labels", dipoles=None, positions=None,
                              fermi_energy=None):
    """Wannier90 model for the metric executor.

    ``spin`` is ``"labels"`` (diagonal Pauli matrices on the spinor labels,
    as Wannier90 orders each spinor projection up then down) or ``"exact"``
    (:class:`WannierSpin`). ``dipoles`` selects the position matrix read from
    ``<seed>_r.dat``: ``None`` (the default) reads it when the file is
    present, ``True`` reads it, ``False`` skips it. The deprecated
    ``positions`` keyword has the same three values: when given it wins and a
    ``DeprecationWarning`` says so. ``fermi_energy`` (eV) selects the energy
    origin and is removed as ``H(R) - E_F 1``; a Wannier90 file carries none,
    so the default is 0.0. Energies, windows and Fermi energies of the run
    file are measured from it.
    """
    if positions is not None:
        import warnings
        warnings.warn("load_wannier_metric_model(positions=...) is deprecated; use dipoles=... "
                      "instead", DeprecationWarning, stacklevel=2)
        dipoles = positions
    base = os.path.join(model_dir, seed)
    if not (os.path.isfile(base + ".uc") and os.path.isfile(base + ".xyz")):
        write_wannier_geometry(model_dir, seed)
    w = load_wannier_model(model_dir, seed)
    hoppings = np.array(w.hoppings, complex)
    origin = 0.0 if fermi_energy is None else float(fermi_energy)
    origin_index = np.where(~np.any(w.r_int, axis=1))[0]
    if len(origin_index) != 1:
        raise ValueError("the Wannier model has no unique R = 0 block")
    hoppings[origin_index[0]] -= origin * np.eye(w.num_wann)
    if dipoles is None:
        dipoles = os.path.isfile(base + "_r.dat")
    pos = read_position_matrix(base + "_r.dat", base + "_hr.dat") if dipoles else None
    if spin not in ("labels", "exact"):
        raise ValueError("Wannier spin must be 'labels' or 'exact'")
    spin_obj = WannierSpin(model_dir, seed, w.uc, w.tau) if spin == "exact" else None
    return OrthogonalModel(w.uc, w.tau, w.r_int, hoppings, seed_name=seed, positions=pos,
                           spin=spin_obj, metadata=dict(source="wannier90", model_dir=os.path.abspath(model_dir),
                                                        fermi_energy_subtracted_eV=origin,
                                                        fermi_energy_source="none" if fermi_energy is None else "run file"))


def write_wannier_geometry(model_dir, seed):
    """Write ``<seed>.uc`` and ``<seed>.xyz`` from ``<seed>.win`` and ``<seed>_centres.xyz``.

    ``.uc`` holds the three lattice vectors of ``unit_cell_cart`` (Angstrom;
    a ``bohr`` unit line is converted). ``.xyz`` holds the Wannier centres
    that Wannier90 writes with ``write_xyz = true``: the lines labelled ``X``
    of ``<seed>_centres.xyz``, in Wannier-function order. Returns the paths.
    """
    base = os.path.join(model_dir, seed)
    lines = open(base + ".win", encoding="utf-8").read().splitlines()
    cell = None
    for i, line in enumerate(lines):
        if line.split("!")[0].strip().lower().startswith("begin unit_cell_cart"):
            rows, j = [], i + 1
            while not lines[j].strip().lower().startswith("end"):
                rows.append(lines[j].split("!")[0].split())
                j += 1
            unit = rows[0][0].lower() if rows and len(rows[0]) == 1 else "ang"
            vectors = np.array([r for r in rows if len(r) == 3], float)
            cell = vectors * (BOHR_ANGSTROM if unit.startswith("bohr") else 1.0)
            break
    if cell is None or cell.shape != (3, 3):
        raise ValueError(f"{base}.win: no unit_cell_cart block with three vectors")
    if not os.path.isfile(base + "_centres.xyz"):
        raise ValueError(f"{base}_centres.xyz not found; run Wannier90 with write_xyz = true")
    centres = [list(map(float, ln.split()[1:4])) for ln in open(base + "_centres.xyz", encoding="utf-8")
               if len(ln.split()) == 4 and ln.split()[0] == "X"]
    np.savetxt(base + ".uc", cell, fmt="%.14f")
    with open(base + ".xyz", "w", encoding="utf-8") as fh:
        fh.write(f"{len(centres)}\n")
        for n, c in enumerate(centres):
            fh.write(f"WF{n + 1} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}\n")
    return base + ".uc", base + ".xyz"


def wigner_seitz_blocks(hr, uc, tau, *, tolerance=1e-5):
    """Place the H(R) of a periodic FFT grid on Wigner-Seitz images, orbital pair by orbital pair.

    ``hr[m, n, i, j, l]`` is PAOFLOW's block on grid index (i, j, l) of an
    n1 x n2 x n3 grid, with ``H(k) = sum_R exp(-i k.R) H(R)``. Each index
    stands for every lattice vector R = (i, j, l) + (s1 n1, s2 n2, s3 n3).
    The element (m, n) is put on the image whose bond, in the stored
    convention R' = -R (``exp(+i k.R')``, atom gauge), is shortest:
    |R'.uc + tau_n - tau_m|. Images whose lengths agree within
    ``tolerance`` (relative) share the element equally. On the DFT mesh this
    changes nothing; between the mesh points it keeps the point symmetry of
    the lattice, which a fixed parallelogram of R (the grid indices folded to
    [-n/2, n/2]) breaks on a hexagonal lattice. Wannier90 does the same with
    ``use_ws_distance``. The pair (n, m) at -R' has the same bond length and
    gets the same weights, so H(k) stays Hermitian. Returns ``(r_int, blocks)``.
    """
    hr = np.asarray(hr)
    nawf = hr.shape[0]
    mesh = np.array(hr.shape[2:5])
    uc, tau = np.asarray(uc, float), np.asarray(tau, float)
    shifts = np.array([(a, b, c) for a in (range(-2, 3) if mesh[0] > 1 else [0])
                       for b in (range(-2, 3) if mesh[1] > 1 else [0])
                       for c in (range(-2, 3) if mesh[2] > 1 else [0])])
    bond0 = tau[None, :, :] - tau[:, None, :]                          # tau_n - tau_m, (m, n, 3)
    accumulated = {}
    for index in np.ndindex(*mesh):
        images = -(np.asarray(index)[None, :] + shifts * mesh[None, :])  # stored R' for each shift
        lengths = np.linalg.norm((images @ uc)[:, None, None, :] + bond0[None], axis=-1)  # (s, m, n)
        shortest = lengths.min(axis=0)
        chosen = lengths <= shortest[None] * (1 + tolerance) + tolerance
        weights = chosen / chosen.sum(axis=0)[None]
        block = hr[:, :, index[0], index[1], index[2]]
        for s in np.where(chosen.any(axis=(1, 2)))[0]:
            key = tuple(int(x) for x in images[s])
            if key not in accumulated:
                accumulated[key] = np.zeros((nawf, nawf), complex)
            accumulated[key] += weights[s] * block
    keys = sorted(accumulated)
    return np.array(keys, int), np.array([accumulated[k] for k in keys])


def write_paoflow_model(paoflow, path, *, sh_l=None, sh_j=None, locality_warning=1e-2):
    """Write the PAO Hamiltonian, orbital positions and spin of a PAOFLOW run to ``path`` (.npz).

    Call after ``pao_hamiltonian()``. PAOFLOW transforms with
    ``H(k) = sum_R exp(-i k.R) H(R)``; the file stores ``H(-R)`` so that the
    LSQuant convention ``exp(+i k.R)`` gives the same ``H(k)``. Energies are
    relative to the DFT Fermi energy, as in PAOFLOW. The spin is PAOFLOW's
    on-site operator times two (Pauli units). Each element of the grid is
    placed on its Wigner-Seitz image (:func:`wigner_seitz_blocks`), so H(k)
    between the DFT k points keeps the symmetry of the lattice.

    ``sh_l``, ``sh_j`` are passed to PAOFLOW's ``spin_operator``: one entry
    per l shell of every atom in basis order, with ``sh_j`` the j of the
    first radial channel of that shell in the pseudopotential. PAOFLOW's own
    default assumes ``j = l - 1/2`` first; fully relativistic UPF files that
    store ``j = l + 1/2`` first (PseudoDojo, for example) need them.

    The file also stores ``edge_ratio``, the largest ``|H(R)|`` on the
    outermost R shell of the PAOFLOW grid over the largest ``|H(R)|`` with
    R != 0. Above ``locality_warning`` a ``RuntimeWarning`` is raised: H(R)
    has not decayed, and H(k) between the DFT k points can hold spurious
    states (a denser nscf mesh or a richer projection basis is needed).
    """
    import warnings
    arrays, attrs = paoflow.data_controller.data_dicts()
    if sh_l is not None or sh_j is not None:
        if sh_l is None or sh_j is None:
            raise ValueError("give both sh_l and sh_j")
        arrays.pop("sh_l", None)
        arrays.pop("sh_j", None)
        paoflow.spin_operator(sh_l=np.asarray(sh_l), sh_j=np.asarray(sh_j))
        arrays, attrs = paoflow.data_controller.data_dicts()
    elif "Sj" not in arrays and attrs.get("dftSO", False):
        paoflow.spin_operator()
        arrays, attrs = paoflow.data_controller.data_dicts()
    if "Sj" not in arrays:
        # Spinless model (no spin-orbit in the DFT run): no spin operator; the
        # optical and charge words then need spin_degeneracy = 2.
        nawf = np.asarray(arrays["HRs"]).shape[0]
        arrays["Sj"] = np.zeros((3, nawf, nawf), complex)
    alat = float(attrs["alat"]) * BOHR_ANGSTROM
    uc = np.asarray(arrays["a_vectors"], float) * alat
    hr = np.asarray(arrays["HRs"])[:, :, :, :, :, 0]
    n1, n2, n3 = hr.shape[2:5]

    tau = (np.asarray(arrays["Dnm"])[:, 0, :] - np.asarray(arrays["Dnm"])[0, 0, :]) * BOHR_ANGSTROM
    r_arr, blocks = wigner_seitz_blocks(hr, uc, tau)
    size = np.abs(blocks).max(axis=(1, 2))
    reach = np.max(np.abs(r_arr) / np.maximum(np.array([n1, n2, n3]) // 2, 1), axis=1)
    inner = size[np.any(r_arr != 0, axis=1)].max() if len(r_arr) > 1 else 1.0
    edge_ratio = float(size[reach >= 1].max() / inner) if np.any(reach >= 1) and inner > 0 else 0.0
    if edge_ratio > locality_warning:
        warnings.warn(f"PAOFLOW H(R) has not decayed at the edge of the {n1}x{n2}x{n3} grid "
                      f"(edge/max = {edge_ratio:.2e}); H(k) off the DFT mesh may hold spurious states",
                      RuntimeWarning, stacklevel=2)
    np.savez(path, uc=uc, tau=tau, r_int=r_arr, hoppings=blocks,
             spin=2.0 * np.asarray(arrays["Sj"]), fermi_energy=float(attrs.get("Efermi", 0.0)),
             mesh=np.array([n1, n2, n3]), edge_ratio=edge_ratio)


def load_paoflow_model(path, *, fermi_energy=None):
    """Read a file of :func:`write_paoflow_model` as an orthogonal model with on-site spin.

    ``fermi_energy`` (eV) selects the energy origin. The file stores H on
    PAOFLOW's scale (its Fermi energy already the origin) and records that
    energy; an explicit value shifts ``H(R=0)`` by ``-(x - E_F(file)) 1``,
    the default keeps the file's origin.
    """
    from .lcao import seed_name_from_file

    d = np.load(path)
    seed = seed_name_from_file(path)
    ef_file = float(d["fermi_energy"])
    if fermi_energy is None:
        origin, source = ef_file, "file"
        hoppings = d["hoppings"]
    else:
        origin, source = float(fermi_energy), "run file"
        hoppings = np.array(d["hoppings"], complex)
        shift = origin - ef_file
        if shift != 0.0:
            r0 = np.where(~np.any(np.asarray(d["r_int"], int), axis=1))[0]
            if len(r0) != 1:
                raise ValueError("the PAOFLOW model has no unique R = 0 block")
            hoppings[r0[0]] -= shift * np.eye(hoppings.shape[1])
    return OrthogonalModel(d["uc"], d["tau"], d["r_int"], hoppings, seed_name=seed,
                           spin=np.asarray(d["spin"], complex),
                           metadata=dict(source="paoflow", file=os.path.abspath(path),
                                         fermi_energy_subtracted_eV=origin,
                                         fermi_energy_source=source))
