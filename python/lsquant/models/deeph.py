"""DeepH-format Hamiltonians: DeepH-pack, DeepH-E3, xDeepH, MACE-H and DeepH-2.

These codes write a predicted H(R) per atom pair next to the DFT overlap
S(R). Two directory layouts exist; :func:`load_deeph_model` detects which.

**Processed-structure layout** (DeepH-pack before its JAX version,
DeepH-E3, xDeepH and MACE-H, which is built on DeepH-E3):

- ``lat.dat``: lattice vectors as *columns*, Angstrom (or ``rlat.dat``,
  the reciprocal vectors as columns, when an inference directory has no
  ``lat.dat``).
- ``site_positions.dat``: Cartesian atomic positions as columns, Angstrom.
- ``orbital_types.dat``: one line per atom with the angular momentum l of
  each radial orbital; atom ``a`` has ``sum(2 l + 1)`` orbitals.
- ``info.json``: ``isspinful``; ``fermi_level`` (eV) when the DFT
  preprocessing wrote it.
- ``hamiltonians_pred.h5`` (prediction) or ``hamiltonians.h5`` (DFT) and
  ``overlaps.h5`` (also read as ``overlap.h5``): one dataset per atom
  pair named ``"[R1, R2, R3, i, j]"``
  with 1-based atoms; the dataset read by h5py has the orbitals of atom
  ``i`` as rows (the Julia tools of DeepH transpose on read, the Python
  tools do not). Energies in eV.

**DeepH-pack layout** (the JAX DeepH-pack and DeepH-dock tools):

- ``POSCAR`` (VASP 5 format, Angstrom) and ``info.json`` with ``spinful``,
  ``orthogonal_basis``, ``fermi_energy_eV`` and ``elements_orbital_map``
  (l of each radial orbital, per element).
- ``hamiltonian_pred.h5`` (prediction) or ``hamiltonian.h5`` (DFT) and
  ``overlap.h5``, each with ``atom_pairs`` (rows ``[R1, R2, R3, i, j]``,
  0-based atoms), ``chunk_boundaries``, ``chunk_shapes`` and ``entries``.
- With ``split_hkb: true`` in ``info.json``, ``hamiltonian*.h5`` holds the
  spinless ``H_0`` and ``hkb.h5`` the spinful Kleinman-Bylander part;
  ``H = I_2 tensor H_0 + H_KB``.

Conventions converted here (both layouts use the same ones):

- Fourier sign: DeepH sums ``H(k) = sum_R exp(+2 pi i k.R) H_ij(R)`` with
  ``H_ij(R) = <i, 0|H|j, R>`` (``dense_calc.py`` of DeepH-pack,
  ``matrix_obj.py`` of DeepH-dock). This is the LSQuant convention; only the
  orbital phases of the atom gauge are added.
- Spinors: a spinful block of an atom pair is ``[[uu, ud], [du, dd]]``
  (``openmx_get_data.jl`` of DeepH-pack; the ABACUS interface permutes to
  this layout). LSQuant interleaves the spinors, see
  :mod:`lsquant.models.lcao`.
- The overlap of a spinful model is stored spinless and lifted to
  ``S tensor I_2``; spin-expanded overlap blocks are accepted when they are
  spin-diagonal and equal.
- Energies: H is stored on the absolute scale of the DFT code.
  ``fermi_energy`` (eV) is subtracted as ``H - E_F S``; by default it is the
  value in ``info.json``. A predicted Hamiltonian can have its Fermi level
  elsewhere: pass the value from the prediction when it differs.

The order of the m components inside a shell differs between DFT codes and
is not converted: no LSQuant operator acts on it.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from .lcao import assemble_blocks, build_lcao_model
from .siesta import spinor_lift

LEGACY_FILES = ("hamiltonians_pred.h5", "hamiltonians.h5")
PACK_FILES = ("hamiltonian_pred.h5", "hamiltonian.h5")


def _h5py():
    try:
        import h5py
    except ImportError as error:  # pragma: no cover - depends on the environment
        raise ImportError("reading DeepH files needs h5py (pip install h5py)") from error
    return h5py


def _array(dataset):
    """Return a dataset as a numpy array; Julia writes complex as fields r and i."""
    data = np.asarray(dataset[()])
    if data.dtype.names and {"r", "i"} <= set(data.dtype.names):
        data = data["r"] + 1j * data["i"]
    return data


def read_poscar(path):
    """Return ``(lattice, elements, cartesian)`` of a VASP 5 POSCAR (Angstrom)."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    scale = float(lines[1].split()[0])
    if scale <= 0:
        raise ValueError(f"{path}: a negative (volume) scale factor is not supported")
    lattice = np.array([[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]) * scale
    symbols = lines[5].split()
    if all(s.lstrip("-").isdigit() for s in symbols):
        raise ValueError(f"{path}: the POSCAR has no element line (VASP 4 format)")
    counts = [int(x) for x in lines[6].split()]
    line = 7
    if lines[line].strip()[:1] in ("s", "S"):
        line += 1
    direct = lines[line].strip()[:1] in ("d", "D")
    n = sum(counts)
    coords = np.array([[float(x) for x in lines[line + 1 + a].split()[:3]] for a in range(n)])
    cartesian = coords @ lattice if direct else coords * scale
    elements = [s for s, c in zip(symbols, counts) for _ in range(c)]
    return lattice, elements, cartesian


def _first_file(root, names, explicit):
    if explicit is not None:
        path = root / explicit
        if not path.is_file():
            raise FileNotFoundError(path)
        return path
    for name in names:
        if (root / name).is_file():
            return root / name
    raise FileNotFoundError(f"no {' or '.join(names)} in {root}")


def _legacy_pairs(path):
    """Yield ``(R, i, j, block)`` from a ``"[R1, R2, R3, i, j]"`` file (1-based atoms)."""
    h5py = _h5py()
    with h5py.File(path, "r") as f:
        for key in f.keys():
            values = [int(x) for x in key.strip()[1:-1].split(",")]
            if len(values) != 5:
                raise ValueError(f"{path}: dataset name {key!r} is not [R1, R2, R3, i, j]")
            yield values[:3], values[3] - 1, values[4] - 1, _array(f[key])


def _pack_pairs(path):
    """Yield ``(R, i, j, block)`` from a chunked DeepH-pack file (0-based atoms)."""
    h5py = _h5py()
    with h5py.File(path, "r") as f:
        pairs = np.asarray(f["atom_pairs"][()], int)
        bounds = np.asarray(f["chunk_boundaries"][()], int)
        shapes = np.asarray(f["chunk_shapes"][()], int)
        entries = _array(f["entries"])
    for n, row in enumerate(pairs):
        yield row[:3], row[3], row[4], entries[bounds[n]:bounds[n + 1]].reshape(shapes[n])


def _overlap(path_pairs, norb, spinful):
    """Assemble S(R), spinless or spin-expanded as the block shapes say."""
    blocks = list(path_pairs)
    expanded = spinful and any(np.asarray(b).shape == (2 * norb[i], 2 * norb[j]) for _, i, j, b in blocks)
    return assemble_blocks(blocks, norb, spinful=expanded, name="S")


def load_deeph_model(directory, *, hamiltonian_file=None, fermi_energy=None):
    """Read a DeepH-format directory as an LSQuant model (non-orthogonal unless stated).

    ``hamiltonian_file`` is a file name inside ``directory``; by default the
    prediction (``hamiltonians_pred.h5`` or ``hamiltonian_pred.h5``) and, if
    absent, the DFT file. ``fermi_energy`` (eV) overrides the value of
    ``info.json``. See the module docstring for the conventions.
    """
    root = Path(directory).expanduser().resolve()
    info = json.loads((root / "info.json").read_text(encoding="utf-8")) if (root / "info.json").is_file() else {}
    meta = dict(source="deeph", directory=str(root))
    hkb = None
    if (root / "POSCAR").is_file() and "elements_orbital_map" in info:
        uc, elements, positions = read_poscar(root / "POSCAR")
        orbital_map = info["elements_orbital_map"]
        norb = np.array([sum(2 * int(l) + 1 for l in orbital_map[e]) for e in elements], int)
        spinful = bool(info.get("spinful", False))
        orthogonal = bool(info.get("orthogonal_basis", False))
        split = bool(info.get("split_hkb", info.get("hamiltonian_storage") not in (None, "full")))
        h_path = _first_file(root, PACK_FILES, hamiltonian_file)
        info_fermi = info.get("fermi_energy_eV")
        pairs, s_path = _pack_pairs, root / "overlap.h5"
        if split:
            if not spinful:
                raise ValueError("info.json: split_hkb needs spinful = true")
            hkb = root / "hkb.h5"
        meta.update(layout="deeph-pack", elements=elements)
    elif (root / "orbital_types.dat").is_file() and any((root / f).is_file() for f in ("lat.dat", "rlat.dat")):
        if (root / "lat.dat").is_file():
            uc = np.loadtxt(root / "lat.dat").T
        else:
            # rlat.dat holds the reciprocal vectors b_i (2 pi / Angstrom) as columns.
            uc = 2 * np.pi * np.linalg.inv(np.loadtxt(root / "rlat.dat").T).T
        positions = np.loadtxt(root / "site_positions.dat").reshape(3, -1).T
        lines = [ln.split() for ln in (root / "orbital_types.dat").read_text(encoding="utf-8").splitlines()
                 if ln.strip()]
        norb = np.array([sum(2 * int(l) + 1 for l in ln) for ln in lines], int)
        spinful = bool(info.get("isspinful", False))
        orthogonal = False
        h_path = _first_file(root, LEGACY_FILES, hamiltonian_file)
        info_fermi = info.get("fermi_level")
        pairs = _legacy_pairs
        s_path = root / "overlaps.h5" if (root / "overlaps.h5").is_file() else root / "overlap.h5"
        meta.update(layout="processed-structure")
    else:
        raise FileNotFoundError(f"{root}: neither POSCAR + info.json (DeepH-pack) nor lat.dat + "
                                "orbital_types.dat (DeepH-E3, MACE-H, older DeepH-pack)")
    if len(norb) != len(positions):
        raise ValueError(f"{root}: {len(norb)} atoms in the orbital list, {len(positions)} positions")
    if hkb is not None:
        base_r, base = assemble_blocks(pairs(h_path), norb, spinful=False, name="H_0")
        kb_r, kb = assemble_blocks(pairs(hkb), norb, spinful=True, name="H_KB")
        cells = {tuple(r): n for n, r in enumerate(kb_r)}
        hamiltonian = list(kb)
        r_int = [tuple(r) for r in kb_r]
        for r, block in zip(base_r, base):
            if tuple(r) not in cells:
                cells[tuple(r)] = len(r_int)
                r_int.append(tuple(r))
                hamiltonian.append(np.zeros_like(kb[0]))
            hamiltonian[cells[tuple(r)]] = hamiltonian[cells[tuple(r)]] + spinor_lift(block)
        r_int, hamiltonian = np.array(r_int, int), np.array(hamiltonian)
        meta["hamiltonian_file"] = [str(h_path), str(hkb)]
    else:
        r_int, hamiltonian = assemble_blocks(pairs(h_path), norb, spinful=spinful, name="H")
        meta["hamiltonian_file"] = str(h_path)
    overlap, s_r_int = None, None
    if not orthogonal:
        if not s_path.is_file():
            raise FileNotFoundError(f"{s_path} not found; DeepH needs the DFT overlap next to the prediction")
        s_r_int, overlap = _overlap(pairs(s_path), norb, spinful)
        meta["overlap_file"] = str(s_path)
    if fermi_energy is None:
        fermi_energy = float(info_fermi) if info_fermi is not None else 0.0
        meta["fermi_energy_source"] = "file" if info_fermi is not None else "none"
    else:
        meta["fermi_energy_source"] = "run file"
    seed = os.path.basename(str(root))
    return build_lcao_model(uc, positions, norb, r_int, hamiltonian, overlap, spinful=spinful,
                            seed_name=seed, metadata=meta, fermi_energy=float(fermi_energy),
                            overlap_r_int=s_r_int)
