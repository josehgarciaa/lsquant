"""ThreeBodyTB (TB3) tight-binding models from their XML file.

ThreeBodyTB.jl (Garrity and Choudhary) builds a non-orthogonal
tight-binding model of a crystal from two- and three-body terms fitted to
DFT. In Julia, ``tbc = calc_tb_fast(crys)`` (or ``scf_energy``) gives a
``tb_crys`` object, and ``write_tb_crys("model.xml", tbc)`` writes it as
XML (also read as ``.xml.gz``). LSQuant reads that file; no export code is
needed. The format was read in the sources of ThreeBodyTB.jl
(``src/TB.jl``: ``write_tb_crys``, ``read_tb_crys``, ``Hk``;
``src/Crystal.jl``: ``orbital_index``; ``src/Atomdata.jl``;
``src/GlobalUnits.jl``).

File content used here:

- ``crystal/A``: lattice vectors as rows, Bohr; ``crystal/coords``:
  fractional positions, one row per atom; ``crystal/types``: elements.
- ``tightbinding/ind_arr``: one lattice vector R (integers) per row;
  ``tightbinding/H``: lines ``r m n Re(H) Im(H) [Re(H_dn) Im(H_dn)] Re(S)
  Im(S)`` with 1-based ``r``, ``m``, ``n``; two H columns when the root
  ``nspin`` is 2 (spin-polarized hoppings).
- ``tightbinding/nonorth``: whether S(R) is used; ``tightbinding/scf`` and
  ``h1``: the self-consistent charge term; ``efermi``: Fermi energy.
- Energies in Rydberg.

Conventions converted here:

- Fourier sign. TB3's ``Hk`` sums ``H_mn(k) = sum_r conj(H_mn(r))
  exp(-2 pi i k.R_r)`` (Julia's ``H[m,n,:]' * exp_ikr``), then takes the
  Hermitian part. Written in the LSQuant form ``sum_R exp(+2 pi i k.R)
  M_mn(R)``, this is ``M(R) = conj(H(-R))``; the same holds for S. The
  Hermitian part is the mean of M(R) and M(-R)^+ of
  :func:`~lsquant.models.lcao.hermitian_pairs`.
- Self-consistent charges. With ``scf = true`` TB3 adds ``S(k) .* h1``
  (element by element) to H(k). Since ``h1`` does not depend on R this is
  ``M(R) += S(R) .* h1`` in real space. ``h1`` is read as Julia reads it
  (``reshape(h1, nwan, nwan)``, column-major).
- Orbitals. TB3 orders orbitals atom by atom, with the shells of
  ``Atomdata.jl`` for each element (``s``: 1, ``p``: 3, ``d``: 5, ``f``: 7
  orbitals). The table :data:`ORBITALS` is copied from that file; pass
  ``orbitals_per_element`` for an element it lacks or a modified database.
  Every orbital sits on its atom.
- Units: Rydberg to eV and Bohr to Angstrom with TB3's own constants
  (``GlobalUnits.jl``), so that energies match TB3's printed eV.
- Fermi energy: ``efermi`` of the file (Rydberg) unless the run file gives
  ``fermi_energy`` (eV); it is removed as ``H - E_F S``.
- Spin: ThreeBodyTB has no spin-orbit coupling. A file with two spin
  channels gives a spinful model without spin mixing (interleaved spinors,
  ``[[H_up, 0], [0, H_dn]]`` per orbital pair).

Warning: ``write_tb_crys`` does not write ``h1spin``, the self-consistent
magnetic term. A file with ``scfspin = true`` is refused, because the model
in the file is not the one TB3 diagonalizes.

Pitfall: the file holds the model at the charges of its last SCF step. A
``tb_crys`` written before ``scf_energy`` has ``scf = false`` and no charge
term, even if the model needs self-consistency (ionic systems).

To-do: check against TB3's own ``Hk`` on a real model (for example
``calc_tb_fast`` of rock-salt NaCl after ``scf_energy``) once Julia and
ThreeBodyTB are available on a test machine.
"""
from __future__ import annotations

import gzip
import os
import xml.etree.ElementTree as ET

import numpy as np

#: Rydberg in eV and Bohr in Angstrom, as in ThreeBodyTB's GlobalUnits.jl
#: (kept in lsquant.models.units under their own names; the Rydberg differs
#: from the ABACUS value and from HARTREE_EV / 2 on purpose).
from .units import BOHR_ANGSTROM
from .units import RYDBERG_EV_THREEBODYTB as RYDBERG_EV
from .siesta import spinor_lift
_SHELL = {"s": 1, "p": 3, "d": 5, "f": 7}
#: Orbital order inside a shell, from ``orbital_index`` in ThreeBodyTB's Crystal.jl.
_LABELS = {"s": ["s"], "p": ["pz", "px", "py"], "d": ["dz2", "dxz", "dyz", "dx2_y2", "dxy"],
           "f": ["fz3", "fxz2", "fyz2", "fz_x2y2", "fxyz", "fx_x2_3y2", "fy_3x2_y2"]}

#: Shells of each element in ThreeBodyTB's Atomdata.jl (usnistgov/ThreeBodyTB.jl 0752192).
ORBITALS = {
    "H": "s", "Li": "sp", "Be": "sp", "B": "sp", "C": "sp", "N": "sp", "O": "sp", "F": "sp",
    "Na": "sp", "Mg": "sp", "Al": "sp", "Si": "sp", "P": "sp", "S": "sp", "Cl": "sp",
    "K": "sdp", "Ca": "sdp", "Sc": "sdp", "Ti": "sdp", "V": "sdp", "Cr": "sdp", "Mn": "sdp",
    "Fe": "sdp", "Co": "sdp", "Ni": "sdp", "Cu": "sdp", "Zn": "sdp", "Ga": "sdp", "Ge": "sp",
    "As": "sp", "Se": "sp", "Br": "sp", "Rb": "sdp", "Sr": "sdp", "Y": "sdp", "Zr": "sdp",
    "Nb": "sdp", "Mo": "sdp", "Tc": "sdp", "Ru": "sdp", "Rh": "sdp", "Pd": "sdp", "Ag": "sdp",
    "Cd": "sdp", "In": "sdp", "Sn": "sp", "Sb": "sp", "Te": "sp", "I": "sp", "Cs": "sdp",
    "Ba": "sdp", "La": "sd", "Hf": "sdp", "Ta": "sdp", "W": "sdp", "Re": "sdp", "Os": "sdp",
    "Ir": "sdp", "Pt": "sdp", "Au": "sdp", "Hg": "sdp", "Tl": "sdp", "Pb": "sp", "Bi": "sp",
    "Hx": "s", "X": "s", "Xa": "s",
}


def orbital_count(element, orbitals_per_element=None):
    """Orbitals of one atom of ``element`` in a ThreeBodyTB model."""
    if orbitals_per_element and element in orbitals_per_element:
        return int(orbitals_per_element[element])
    if element not in ORBITALS:
        raise ValueError(f"element {element!r} is not in ThreeBodyTB's table; give "
                         "params.orbitals_per_element")
    return sum(_SHELL[s] for s in ORBITALS[element])


def orbital_labels(elements, orbitals_per_element=None):
    """``[(atom, label), ...]`` in ThreeBodyTB's order, or ``None`` if an element has no table entry.

    The shells follow the order of ``Atomdata.jl`` (``s, d, p`` for the
    transition metals), as ``orbital_index`` does.
    """
    out = []
    for a, element in enumerate(elements):
        if (orbitals_per_element and element in orbitals_per_element) or element not in ORBITALS:
            return None
        out += [(a, label) for shell in ORBITALS[element] for label in _LABELS[shell]]
    return out


def _text(node, path, default=None):
    found = node.find(path)
    if found is None or found.text is None:
        if default is not None:
            return default
        raise ValueError(f"ThreeBodyTB file: missing <{path}>")
    return found.text.strip()


def _floats(text):
    return np.array(text.split(), float)


def _bool(text):
    return text.strip().lower() == "true"


def read_threebodytb_xml(path):
    """Return the raw content of a ``write_tb_crys`` file as a dict (Rydberg, Bohr, TB3 indices)."""
    with open(path, "rb") as f:
        raw = f.read()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    root = ET.fromstring(raw)
    crystal, tb = root.find("crystal"), root.find("tightbinding")
    if crystal is None or tb is None:
        raise ValueError(f"{path} is not a ThreeBodyTB tb_crys file (no <crystal> or <tightbinding>)")
    nat = int(_text(crystal, "nat"))
    nr, nwan = int(_text(tb, "nr")), int(_text(tb, "nwan"))
    nspin = int(_text(root, "nspin", "1"))
    out = dict(
        lattice=_floats(_text(crystal, "A")).reshape(3, 3),
        coords=_floats(_text(crystal, "coords")).reshape(nat, 3),
        types=_text(crystal, "types").split(),
        nspin=nspin, nr=nr, nwan=nwan,
        nonorth=_bool(_text(tb, "nonorth")),
        ind_arr=np.rint(_floats(_text(tb, "ind_arr"))).astype(int).reshape(nr, 3),
        scf=_bool(_text(tb, "scf", "false")),
        scfspin=_bool(_text(tb, "scfspin", "false")),
        efermi=float(_text(root, "efermi", "0.0")),
        nelec=float(_text(root, "nelec", "0.0")),
    )
    if len(out["types"]) != nat:
        raise ValueError(f"{path}: {len(out['types'])} types for {nat} atoms")
    h1 = tb.find("h1")
    out["h1"] = None
    if out["scf"] and h1 is not None and h1.text and h1.text.split():
        out["h1"] = _floats(h1.text).reshape(nwan, nwan, order="F")   # Julia's reshape
    H = np.zeros((nspin, nr, nwan, nwan), complex)
    S = np.zeros((nr, nwan, nwan), complex)
    for line in _text(tb, "H").splitlines():
        sp = line.split()
        if len(sp) not in (7, 9):
            continue
        if len(sp) != 5 + 2 * nspin:
            raise ValueError(f"{path}: H line with {len(sp)} columns, expected {5 + 2 * nspin} for nspin={nspin}")
        r, m, n = int(sp[0]) - 1, int(sp[1]) - 1, int(sp[2]) - 1
        for s in range(nspin):
            H[s, r, m, n] = float(sp[3 + 2 * s]) + 1j * float(sp[4 + 2 * s])
        S[r, m, n] = float(sp[3 + 2 * nspin]) + 1j * float(sp[4 + 2 * nspin])
    out["H"], out["S"] = H, S
    return out


def load_threebodytb_model(path, *, fermi_energy=None, orbitals_per_element=None):
    """Read a ThreeBodyTB ``write_tb_crys`` file (see the module docstring) as an LSQuant model."""
    from .lcao import build_lcao_model

    path = os.path.abspath(os.path.expanduser(str(path)))
    d = read_threebodytb_xml(path)
    if d["scfspin"]:
        raise ValueError("the ThreeBodyTB file has scfspin = true, but write_tb_crys does not write the "
                         "magnetic term h1spin; the model in the file is incomplete")
    norb = np.array([orbital_count(t, orbitals_per_element) for t in d["types"]], int)
    if norb.sum() != d["nwan"]:
        raise ValueError(f"the elements give {norb.sum()} orbitals but the file has nwan = {d['nwan']}; "
                         "give params.orbitals_per_element as in your ThreeBodyTB database")
    if d["scf"] and not d["nonorth"] and d["h1"] is not None:
        raise ValueError("an orthogonal ThreeBodyTB model with a charge term is not supported")
    # TB3 sums conj(H(r)) exp(-2 pi i k.R_r): M(R) = conj(H(-R)), M(-R_r) = conj(H(r)).
    r_int = -d["ind_arr"]
    hamiltonian = np.conj(d["H"])
    overlap = np.conj(d["S"]) if d["nonorth"] else None
    if d["h1"] is not None:
        hamiltonian = hamiltonian + overlap[None] * d["h1"][None, None]
    hamiltonian = hamiltonian * RYDBERG_EV
    spinful = d["nspin"] == 2
    if spinful:
        up, dn = np.diag([1.0, 0.0]), np.diag([0.0, 1.0])
        hamiltonian = np.array([np.kron(hu, up) + np.kron(hd, dn) for hu, hd in zip(*hamiltonian)])
        if overlap is not None:
            overlap = spinor_lift(overlap)
    else:
        hamiltonian = hamiltonian[0]
    uc = d["lattice"] * BOHR_ANGSTROM
    positions = d["coords"] @ uc
    e_f = d["efermi"] * RYDBERG_EV if fermi_energy is None else float(fermi_energy)
    meta = dict(source="threebodytb", file=path, elements=list(d["types"]), scf=bool(d["scf"]),
                nonorthogonal=bool(d["nonorth"]), nspin=int(d["nspin"]),
                fermi_energy_from_file_eV=d["efermi"] * RYDBERG_EV, nelec=d["nelec"],
                basis="ThreeBodyTB minimal atomic orbitals (fit to DFT, see Atomdata.jl)",
                orbital_positions="atom-centred: every orbital on its atom",
                orbital_labels=orbital_labels(d["types"], orbitals_per_element))
    model = build_lcao_model(uc, positions, norb, r_int, hamiltonian, overlap, spinful=spinful,
                             seed_name=os.path.basename(path).split(".")[0], metadata=meta,
                             fermi_energy=e_f)
    if spinful:
        model.metadata["spin_operator"] = "collinear (no spin-orbit coupling)"
    return model
