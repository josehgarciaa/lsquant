"""DeePTB Hamiltonians (DeePTB-SK and DeePTB-E3) from DeePTB's sparse export.

DeePTB writes its real-space Hamiltonian for its Julia/Pardiso band code
with ``TBSystem.to_pardiso_json(output_dir)`` (also run by
``dptb pdso ... -i model.pth -stu structure``) or the older
``TBSystem.to_pardiso(output_dir)``. The directory holds:

- ``predicted_hamiltonians.h5``: group ``"0"`` with one dataset per atom
  pair named ``"i_j_R1_R2_R3"`` (0-based atoms). DeePTB completes every
  block with its partner ``"j_i_-R1_-R2_-R3"`` = block^+ before writing.
- ``predicted_overlaps.h5``: the same for S(R), only for models with an
  overlap (non-orthogonal DeePTB-SK and DeePTB-E3 models trained on LCAO
  data); without it the model is orthogonal.
- ``structure.json`` (``to_pardiso_json``): ``structure.cell`` (rows,
  Angstrom), ``structure.positions`` (Cartesian, Angstrom),
  ``structure.symbols`` and ``basis_info.orbital_counts`` (spinless
  orbitals per element), ``basis_info.spinful``; or, from ``to_pardiso``,
  ``cell.dat``, ``positions.dat``, ``atomic_numbers.dat`` and ``basis.dat``
  (a Python dict such as ``{'Si': '1s1p1d'}``).

Conventions converted here:

- Fourier sign: DeePTB's Julia backend (``pardiso/utils/hamiltonian.jl``)
  sums ``H(k) = sum_R exp(+2 pi i k.R) H(R)`` with the block of
  ``"i_j_R"`` in rows of atom ``i`` and columns of atom ``j``, which is the
  LSQuant convention ``<i, 0|H|j, R>``. DeePTB's Python ``HR2HK`` uses
  ``exp(-i k.R)``: for its real spinless blocks this is H(-k) = H(k)*, with
  the same bands.
- Units: eV and Angstrom (ASE).
- Energies are on DeePTB's scale; the files hold no Fermi energy. Pass the
  one of ``TBSystem.get_efermi`` as ``fermi_energy``.
- Spin: the export holds spinless blocks only. DeePTB-SK adds its on-site
  spin-orbit coupling when it builds H(k), and this term is not in the
  files, so a model exported with ``spinful: true`` is refused.
"""
from __future__ import annotations

import ast
import json
import os
import re
from pathlib import Path

import numpy as np

from .lcao import assemble_blocks, build_lcao_model

_L_SIZE = {"s": 1, "p": 3, "d": 5, "f": 7, "g": 9}
_SYMBOLS = ("X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se "
            "Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy "
            "Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn").split()


def _count_orbitals(compressed):
    """Orbitals of a DeePTB basis string such as ``'2s1p1d'``."""
    return sum(int(n) * _L_SIZE[l] for n, l in re.findall(r"(\d+)([spdfg])", compressed))


def _read_structure(root):
    """Return ``(cell, positions, symbols, orbital_counts, spinful)`` of an export directory."""
    if (root / "structure.json").is_file():
        data = json.loads((root / "structure.json").read_text(encoding="utf-8"))
        st, basis = data["structure"], data["basis_info"]
        return (np.array(st["cell"], float), np.array(st["positions"], float), list(st["symbols"]),
                {k: int(v) for k, v in basis["orbital_counts"].items()}, bool(basis.get("spinful", False)))
    if (root / "cell.dat").is_file():
        cell = np.loadtxt(root / "cell.dat").reshape(3, 3)
        positions = np.loadtxt(root / "positions.dat").reshape(-1, 3)
        numbers = np.loadtxt(root / "atomic_numbers.dat", dtype=int).reshape(-1)
        basis = ast.literal_eval((root / "basis.dat").read_text(encoding="utf-8"))
        return (cell, positions, [_SYMBOLS[z] for z in numbers],
                {k: _count_orbitals(v) for k, v in basis.items()}, False)
    raise FileNotFoundError(f"{root}: no structure.json (to_pardiso_json) or cell.dat (to_pardiso)")


def _pairs(path):
    """Yield ``(R, i, j, block)`` from a DeePTB export file."""
    import h5py
    with h5py.File(path, "r") as f:
        group = f["0"] if "0" in f else f
        for key in group.keys():
            parts = [int(x) for x in key.split("_")]
            if len(parts) != 5:
                raise ValueError(f"{path}: dataset name {key!r} is not i_j_R1_R2_R3")
            yield parts[2:], parts[0], parts[1], np.asarray(group[key][()])


def load_deeptb_model(directory, *, fermi_energy=None):
    """Read a DeePTB export directory (see the module docstring) as an LSQuant model.

    ``fermi_energy`` (eV) selects the energy origin; the files carry none, so
    the default is 0.0.
    """
    root = Path(directory).expanduser().resolve()
    cell, positions, symbols, counts, spinful = _read_structure(root)
    norb = np.array([counts[s] for s in symbols], int)
    h_path = root / "predicted_hamiltonians.h5"
    if not h_path.is_file():
        raise FileNotFoundError(f"{h_path} not found; export with TBSystem.to_pardiso_json")
    if spinful:
        raise ValueError("structure.json marks a spinful (spin-orbit) DeePTB model, but DeePTB exports "
                         "only its spinless blocks; the spin-orbit term is not in the files")
    r_int, hamiltonian = assemble_blocks(_pairs(h_path), norb, spinful=False, name="H")
    if fermi_energy is None:
        origin, source = 0.0, "none"
    else:
        origin, source = float(fermi_energy), "run file"
    meta = dict(source="deeptb", directory=str(root), hamiltonian_file=str(h_path), elements=symbols,
                fermi_energy_source=source)
    overlap = s_r_int = None
    s_path = root / "predicted_overlaps.h5"
    if s_path.is_file():
        s_r_int, overlap = assemble_blocks(_pairs(s_path), norb, spinful=False, name="S")
        meta["overlap_file"] = str(s_path)
    return build_lcao_model(cell, positions, norb, r_int, hamiltonian, overlap, spinful=False,
                            seed_name=os.path.basename(str(root)), metadata=meta,
                            fermi_energy=origin, overlap_r_int=s_r_int)
