"""Spatial domains, orbital permutation and atom-block array assembly.

Pure array code backing ``numerics.operator = "atom_block"`` (native profile): builds the
orbital permutation (atoms contiguous, atoms grouped in spatial domains, Morton order inside
and across domains, interior atoms before boundary atoms) and the plain arrays consumed by
``_core.AtomBlockOperator`` (atom offsets, neighbour lists, block ids, the deduplicated block
dictionary, the on-site diagonal). No physics lives here: callers supply positions, atom maps,
neighbour graphs and block values.
"""
from __future__ import annotations

import os
from typing import Any

import numpy as np


def read_hr_file(path: str) -> np.ndarray:
    """Read an ``H(R)`` hopping list: header ``W nhops``, then ``r1 r2 i j re im`` rows."""
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().split()
        if len(header) != 2:
            raise ValueError(f"H(R) file '{path}' must start with 'W nhops'")
        w, nhops = int(header[0]), int(header[1])
        hops = np.loadtxt(fh, ndmin=2)
    if hops.shape != (nhops, 6):
        raise ValueError(f"H(R) file '{path}': expected {nhops} hopping rows, got {hops.shape}")
    if int(hops[:, 2].max()) >= w or int(hops[:, 3].max()) >= w:
        raise ValueError(f"H(R) file '{path}': orbital index out of range for W={w}")
    return hops


def read_positions(source: Any) -> np.ndarray:
    """Read orbital positions from a path or return an array as ``(N, 3)`` float64.

    Accepted file formats: raw text with ``x y z`` per line, or an ``.xyz``-like file
    (first line the count, second a comment, then ``tag x y z`` per line). Blank lines and
    ``#`` comments are skipped.
    """
    if isinstance(source, (list, tuple, np.ndarray)):
        arr = np.asarray(source, dtype=np.float64)
        if arr.ndim != 2 or arr.shape[1] not in (2, 3):
            raise ValueError("model.positions array must have shape (N, 2) or (N, 3)")
        if arr.shape[1] == 2:
            arr = np.column_stack([arr, np.zeros(len(arr))])
        return np.ascontiguousarray(arr)
    path = os.fspath(source)
    with open(path, encoding="utf-8") as fh:
        first = fh.readline()
        while first and (not first.strip() or first.lstrip().startswith("#")):
            first = fh.readline()
        if not first:
            raise ValueError(f"positions file '{path}' is empty")
        fields = first.split()
        if len(fields) == 1:
            count = int(fields[0])
            fh.readline()  # XYZ comment (may be empty).
            arr = np.loadtxt(fh, usecols=(-3, -2, -1), ndmin=2)
            if len(arr) != count:
                raise ValueError(f"positions file '{path}': expected {count} rows, got {len(arr)}")
        else:
            fh.seek(0)
            arr = np.loadtxt(fh, usecols=(-3, -2, -1), ndmin=2)
    return np.ascontiguousarray(arr)


def atom_map_from_spec(
    dim: int,
    orbitals_per_atom: int | None = None,
    atom_map: Any = None,
) -> np.ndarray:
    """Return the atom index of every orbital (int64, consecutive grouping).

    Exactly one of ``orbitals_per_atom`` (uniform width) or ``atom_map`` (path to one integer
    per line, or an array of length ``dim``) may be given.
    """
    if orbitals_per_atom is not None and atom_map is not None:
        raise ValueError("model gives both orbitals_per_atom and atom_map; expected one")
    if atom_map is not None:
        if isinstance(atom_map, (list, tuple, np.ndarray)):
            amap = np.asarray(atom_map, dtype=np.int64)
        else:
            amap = np.loadtxt(os.fspath(atom_map), dtype=np.int64)
        if amap.ndim != 1 or amap.size != dim:
            raise ValueError(f"atom_map must hold one entry per orbital ({dim})")
        if amap.min() < 0 or np.any(np.diff(amap) < 0):
            raise ValueError("atom_map must be non-decreasing (atoms are consecutive orbitals)")
        return amap
    if orbitals_per_atom is None:
        raise ValueError(
            "numerics.operator='atom_block' requires model.orbitals_per_atom or model.atom_map"
        )
    w = int(orbitals_per_atom)
    if w <= 0 or dim % w != 0:
        raise ValueError(f"orbitals_per_atom={w} does not divide the dimension {dim}")
    return np.repeat(np.arange(dim // w, dtype=np.int64), w)


def _spread_bits(v: np.ndarray) -> np.ndarray:
    """Spread the low 21 bits of ``v`` (uint64) with two zero bits between each."""
    v = v.astype(np.uint64) & np.uint64(0x1FFFFF)
    v = (v | (v << np.uint64(32))) & np.uint64(0x1F00000000FFFF)
    v = (v | (v << np.uint64(16))) & np.uint64(0x1F0000FF0000FF)
    v = (v | (v << np.uint64(8))) & np.uint64(0x100F00F00F00F00F)
    v = (v | (v << np.uint64(4))) & np.uint64(0x10C30C30C30C30C3)
    v = (v | (v << np.uint64(2))) & np.uint64(0x1249249249249249)
    return v


def morton_codes(coords: np.ndarray) -> np.ndarray:
    """Morton (Z-order) code of integer coordinates ``(M, d<=3)``; low 21 bits per axis."""
    coords = np.asarray(coords)
    codes = _spread_bits(coords[:, 0])
    for d in range(1, coords.shape[1]):
        codes |= _spread_bits(coords[:, d]) << np.uint64(d)
    return codes


def default_domain_atoms(max_width: int, state_block: int = 1, vectors_in_flight: int = 3) -> int:
    """Default domain size: one domain's working set within half of a 1 MiB L2 cache.

    The working set per atom is ``W * vectors_in_flight * min(S, 8)`` complex scalars;
    clamped to [64, 16384] atoms so domains never degenerate to single atoms or single boxes.
    """
    per_atom = max(1, max_width) * vectors_in_flight * min(max(state_block, 1), 8) * 16
    return int(np.clip((1 << 19) // per_atom, 64, 16384))


def spatial_permutation(
    atom_pos: np.ndarray,
    pair_a: np.ndarray,
    pair_b: np.ndarray,
    domain_atoms: int,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Order atoms into spatial domains; return ``(atom_order, info)``.

    ``atom_order[new] = old`` groups atoms into boxes of about ``domain_atoms`` atoms, orders
    boxes by the Morton code of their coordinates, atoms inside a box by the Morton code of
    refined sub-box coordinates, interior atoms before atoms with a neighbour in another box.
    ``pair_a``/``pair_b`` are the (old) atom indices of the neighbour graph (used only for the
    interior/boundary split).
    """
    pos = np.asarray(atom_pos, dtype=np.float64)
    num_atoms = pos.shape[0]
    lo = pos.min(axis=0)
    hi = pos.max(axis=0)
    extent = np.maximum(hi - lo, 0.0)
    keep = extent > 0.05 * extent.max()
    if not np.any(keep):
        keep[:] = True
    dims = int(keep.sum())
    sub = (pos[:, keep] - lo[keep]) / np.maximum(extent[keep], 1e-300)
    geo = float(np.exp(np.log(np.maximum(extent[keep], 1e-300)).mean()))
    side = (num_atoms / max(domain_atoms, 1)) ** (1.0 / dims)
    grid = np.maximum(1, np.rint(side * extent[keep] / geo).astype(np.int64))
    box = np.minimum((sub * grid).astype(np.int64), grid - 1)
    box_code = morton_codes(box[:, :3] if dims <= 3 else box)
    refined = np.minimum((sub * grid * 16).astype(np.int64), grid * 16 - 1)
    sub_code = morton_codes(refined[:, :3] if dims <= 3 else refined)
    domain_id = np.ravel_multi_index(box.T, tuple(int(g) for g in grid)) if dims > 1 else box[:, 0]
    pa, pb = np.asarray(pair_a), np.asarray(pair_b)
    boundary = np.zeros(num_atoms, dtype=bool)
    cross = domain_id[pa] != domain_id[pb]
    boundary[pa[cross]] = True
    boundary[pb[cross]] = True
    order = np.lexsort((sub_code, boundary.astype(np.int8), box_code))
    info = {
        "num_atoms": int(num_atoms),
        "num_domains": int(len(np.unique(domain_id))),
        "domain_atoms": int(domain_atoms),
        "box_grid": [int(g) for g in grid],
        "boundary_atoms": int(boundary.sum()),
    }
    return order.astype(np.int64), info


def orbital_permutation(atom_order: np.ndarray, atom_of_orbital: np.ndarray) -> np.ndarray:
    """Expand an atom order into an orbital permutation ``perm[new] = old``.

    Orbitals of one atom stay consecutive in their original relative order.
    ``atom_of_orbital`` must be non-decreasing (atoms are consecutive orbital ranges in
    ascending atom order, as produced by :func:`atom_map_from_spec`).
    """
    atom_of_orbital = np.asarray(atom_of_orbital)
    if atom_of_orbital.size and np.any(np.diff(atom_of_orbital) < 0):
        raise ValueError(
            "orbital_permutation: atom_of_orbital must be non-decreasing "
            "(atoms are consecutive orbital ranges in ascending atom order)"
        )
    counts = np.bincount(atom_of_orbital)
    offsets = np.concatenate(([0], np.cumsum(counts)))
    widths = counts[atom_order]
    starts = np.repeat(offsets[atom_order], widths)
    new_starts = np.repeat(np.cumsum(widths) - widths, widths)
    return starts + np.arange(len(atom_of_orbital), dtype=np.int64) - new_starts


def _pad4(n: int) -> int:
    # Tight packing: dictionary blocks are stored back to back with no padding (the C++
    # operator validates dictOffset[d+1] - dictOffset[d] == rows*cols). W >= 2 block sizes
    # are natural multiples of 64 B (2x2 = 4 and 12x12 = 144 complex scalars), so only the
    # scalar W = 1 case is unaligned, where SIMD alignment is irrelevant.
    return n


def blocks_from_csr(
    rowptr: np.ndarray,
    col: np.ndarray,
    val: np.ndarray,
    atom_of_orbital: np.ndarray,
    perm: np.ndarray,
    diagonal: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Convert a CSR operator into atom-block arrays under the orbital permutation ``perm``.

    ``perm[new] = old``. The diagonal (entries whose row and column are the same orbital) is
    split off and returned separately when ``diagonal`` is not given; pass an explicit array to
    supply the on-site values (e.g. Anderson disorder) instead of taking them from the CSR.
    """
    from . import _core

    arrays = _core.AtomBlocksFromCsr(
        np.ascontiguousarray(rowptr, dtype=np.int64),
        np.ascontiguousarray(col, dtype=np.int64),
        np.ascontiguousarray(val, dtype=np.complex128),
        np.ascontiguousarray(atom_of_orbital, dtype=np.int64),
        np.ascontiguousarray(perm, dtype=np.int64),
    )
    if diagonal is not None:
        arrays["diagonal"] = np.ascontiguousarray(diagonal, dtype=np.complex128)[perm]
    return arrays


def blocks_from_hoppings(
    hops: np.ndarray,
    shape: tuple[int, int],
    positions: np.ndarray,
    orbitals_per_atom: int,
    domain_atoms: int,
) -> tuple[dict[str, np.ndarray], np.ndarray, dict[str, Any]]:
    """Build atom-block arrays directly from an ``H(R)`` list on a periodic supercell.

    ``hops`` rows are ``(r1, r2, i, j, re, im)`` (orbital ``i`` of cell 0 to orbital ``j`` of
    cell ``R``), ``shape`` is ``(L1, L2)`` and ``positions`` the ``(N, 3)`` orbital positions
    of the supercell in natural cell-major order. One atom = one cell of
    ``orbitals_per_atom`` orbitals. The dictionary is assembled per hopping vector, so a
    crystal never materialises a CSR. Returns ``(arrays, perm, info)``.
    """
    l1, l2 = int(shape[0]), int(shape[1])
    w = int(orbitals_per_atom)
    num_atoms = l1 * l2
    n = num_atoms * w
    hops = np.asarray(hops, dtype=np.float64)
    atom_of_orbital = np.repeat(np.arange(num_atoms, dtype=np.int64), w)
    atom_pos = positions.reshape(num_atoms, w, 3).mean(axis=1)
    r_set = np.unique(hops[:, :2].astype(np.int64), axis=0)
    rmax = int(np.abs(r_set).max()) if len(r_set) else 0
    cell1 = np.arange(num_atoms, dtype=np.int64) // l2
    cell2 = np.arange(num_atoms, dtype=np.int64) % l2
    pair_a_l: list[np.ndarray] = []
    pair_b_l: list[np.ndarray] = []
    for r1, r2 in r_set:
        pair_a_l.append(np.arange(num_atoms, dtype=np.int64))
        pair_b_l.append((((cell1 + r1) % l1) * l2 + (cell2 + r2) % l2))
    pair_a = np.concatenate(pair_a_l)
    pair_b = np.concatenate(pair_b_l)
    atom_order, info = spatial_permutation(atom_pos, pair_a, pair_b, domain_atoms)
    info["hopping_range"] = rmax
    new_of_old_atom = np.empty(num_atoms, dtype=np.int64)
    new_of_old_atom[atom_order] = np.arange(num_atoms, dtype=np.int64)
    cells_of_new = np.stack([cell1[atom_order], cell2[atom_order]], axis=1)

    k_num = len(r_set)
    nbr = np.empty((num_atoms, k_num), dtype=np.int32)
    for k, (r1, r2) in enumerate(r_set):
        nb1 = (cells_of_new[:, 0] + r1) % l1
        nb2 = (cells_of_new[:, 1] + r2) % l2
        nbr[:, k] = new_of_old_atom[nb1 * l2 + nb2].astype(np.int32)

    blocks = hoppings_dictionary(hops, w, positions[:w], np.zeros((3, 2)), None)
    dict_rows = np.full(k_num, w, dtype=np.int32)
    dict_cols = np.full(k_num, w, dtype=np.int32)
    dict_offset = np.zeros(k_num + 1, dtype=np.int64)
    cur = 0
    for k in range(k_num):
        dict_offset[k] = cur
        cur += _pad4(w * w)
    dict_offset[-1] = cur
    dict_values = np.zeros(cur, dtype=np.complex128)
    for k in range(k_num):
        dict_values[dict_offset[k] : dict_offset[k] + w * w] = blocks[k].ravel()

    atom_offset = np.arange(0, n + 1, w, dtype=np.int64)
    pair_offset = np.arange(0, (num_atoms + 1) * k_num, k_num, dtype=np.int64)
    perm = atom_order.repeat(w) * w + np.tile(np.arange(w, dtype=np.int64), num_atoms)
    arrays = {
        "atom_offset": atom_offset,
        "pair_offset": pair_offset,
        "nbr": nbr.ravel(),
        "block_id": np.tile(np.arange(k_num, dtype=np.int32), num_atoms),
        "dict_offset": dict_offset,
        "dict_rows": dict_rows,
        "dict_cols": dict_cols,
        "dict_values": dict_values,
        "diagonal": np.zeros(0, dtype=np.complex128),
    }
    return arrays, perm, info


def velocity_arrays_from_hoppings(
    hops: np.ndarray,
    base_arrays: dict[str, np.ndarray],
    positions_cell: np.ndarray,
    lattice: np.ndarray,
    axis: int,
) -> dict[str, np.ndarray]:
    """Velocity arrays ``i [H, X_axis]`` sharing the pair structure of ``base_arrays``.

    ``base_arrays`` must come from :func:`blocks_from_hoppings` (one dictionary block per
    hopping vector, ``np.unique`` order); the velocity dictionary is
    :func:`hoppings_dictionary` on the same hopping list, so deduplication is identical to the
    Hamiltonian's on a crystal. ``positions_cell`` are the ``(W, 3)`` orbital positions inside
    one cell and ``lattice`` the ``(3, 2)`` matrix of primitive vectors (columns). The bond
    vector is the unwrapped ``R . a + pos_j - pos_i``, identical to the minimum-image
    convention whenever the hopping range is below half the supercell.
    """
    blocks = hoppings_dictionary(
        hops,
        int(base_arrays["dict_rows"][0]),
        positions_cell,
        lattice,
        axis,
    )
    w = blocks.shape[1]
    dict_values = np.zeros_like(base_arrays["dict_values"])
    for k in range(len(blocks)):
        start = int(base_arrays["dict_offset"][k])
        dict_values[start : start + w * w] = blocks[k].ravel()
    out = dict(base_arrays)
    out["dict_values"] = dict_values
    return out


def hoppings_dictionary(
    hops: np.ndarray,
    orbitals_per_atom: int,
    positions_cell: np.ndarray,
    lattice: np.ndarray,
    axis: int | None = None,
) -> np.ndarray:
    """Per-hopping-vector ``W x W`` blocks of ``H`` (axis None) or of ``i [H, X_axis]``.

    ``positions_cell`` are the ``(W, 3)`` orbital positions inside one cell and ``lattice``
    the ``(3, 2)`` matrix whose columns are the two primitive vectors. Returns blocks in the
    order of ``np.unique(hops[:, :2])`` — the dictionary order of ``blocks_from_hoppings``.
    """
    hops = np.asarray(hops, dtype=np.float64)
    w = int(orbitals_per_atom)
    r_set = np.unique(hops[:, :2].astype(np.int64), axis=0)
    r_index = {(int(r1), int(r2)): k for k, (r1, r2) in enumerate(r_set)}
    blocks = np.zeros((len(r_set), w, w), dtype=np.complex128)
    lat = np.asarray(lattice, dtype=np.float64)
    pos = np.asarray(positions_cell, dtype=np.float64)
    for r1, r2, i, j, re, im in hops:
        k = r_index[(int(r1), int(r2))]
        value = complex(re, im)
        if axis is not None:
            dr = r1 * lat[:, 0] + r2 * lat[:, 1] + pos[int(j)] - pos[int(i)]
            value = 1j * dr[int(axis)] * value
        blocks[k, int(i), int(j)] += value
    return blocks


def to_core(arrays: dict[str, Any], hermitian: bool = True):
    """Instantiate ``_core.AtomBlockOperator`` from an arrays dict.

    The dict keys match the binding's argument names exactly (``atom_offset``,
    ``pair_offset``, ``nbr``, ``block_id``, ``dict_offset``, ``dict_rows``, ``dict_cols``,
    ``dict_values``, ``diagonal``).
    """
    from . import _core

    return _core.AtomBlockOperator(
        atom_offset=np.ascontiguousarray(arrays["atom_offset"], dtype=np.int64),
        pair_offset=np.ascontiguousarray(arrays["pair_offset"], dtype=np.int64),
        nbr=np.ascontiguousarray(arrays["nbr"], dtype=np.int32),
        block_id=np.ascontiguousarray(arrays["block_id"], dtype=np.int32),
        dict_offset=np.ascontiguousarray(arrays["dict_offset"], dtype=np.int64),
        dict_rows=np.ascontiguousarray(arrays["dict_rows"], dtype=np.int32),
        dict_cols=np.ascontiguousarray(arrays["dict_cols"], dtype=np.int32),
        dict_values=np.ascontiguousarray(arrays["dict_values"], dtype=np.complex128),
        diagonal=np.ascontiguousarray(arrays["diagonal"], dtype=np.complex128),
        hermitian=hermitian,
    )


def to_core_cuda(arrays: dict[str, Any], hermitian: bool = True):
    """Instantiate ``_core.CudaAtomBlockOperator`` from an arrays dict (CUDA build only)."""
    from . import _core

    return _core.CudaAtomBlockOperator(
        atom_offset=np.ascontiguousarray(arrays["atom_offset"], dtype=np.int64),
        pair_offset=np.ascontiguousarray(arrays["pair_offset"], dtype=np.int64),
        nbr=np.ascontiguousarray(arrays["nbr"], dtype=np.int32),
        block_id=np.ascontiguousarray(arrays["block_id"], dtype=np.int32),
        dict_offset=np.ascontiguousarray(arrays["dict_offset"], dtype=np.int64),
        dict_rows=np.ascontiguousarray(arrays["dict_rows"], dtype=np.int32),
        dict_cols=np.ascontiguousarray(arrays["dict_cols"], dtype=np.int32),
        dict_values=np.ascontiguousarray(arrays["dict_values"], dtype=np.complex128),
        diagonal=np.ascontiguousarray(arrays["diagonal"], dtype=np.complex128),
        hermitian=hermitian,
    )
