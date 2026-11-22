"""HamGNN Hamiltonians: export in the HamGNN environment, read in LSQuant.

HamGNN writes no stand-alone Hamiltonian file. A prediction run
(``stage: test`` in its config) saves ``prediction_hamiltonian.npy`` in the
log directory, and the structures live in the ``graph_data.npz`` of the
run: a pickled dictionary of ``torch_geometric`` ``Data`` objects, which
only opens where ``torch`` and ``torch_geometric`` are installed. The export
is therefore split in two steps:

1. :func:`write_hamgnn_model` runs in the Python environment of HamGNN and
   writes one structure to a plain ``.npz`` (numpy arrays only). This file
   imports nothing from LSQuant at module level, so it also runs as a
   script without LSQuant installed:
   ``python hamgnn.py graph_data.npz prediction_hamiltonian.npy model.npz --nao-max 26``.
2. :func:`load_hamgnn_model` reads that file in LSQuant (``model.kind =
   "hamgnn"`` in a run file).

HamGNN conventions, read in its sources (``DFT_interfaces/openmx/band_cal.py``,
``graph_data_gen.py``, ``openmx_postprocess/read_openmx.c``,
``abacus/band_cal.py``, ``siesta/read_siesta.py``), and their conversion:

- Blocks: ``Hon[a]`` is the on-site block of atom ``a``; ``Hoff[e]`` couples
  atom ``edge_index[0, e]`` in the home cell (rows) with atom
  ``edge_index[1, e]`` in cell ``cell_shift[e]`` (columns), and both
  directions of every bond are stored. HamGNN sums
  ``H(k) = sum_e exp(+i k.nbr_shift[e]) Hoff[e] + Hon`` with
  ``nbr_shift = cell_shift @ cell``. This is the LSQuant convention
  ``<i, 0|H|j, R>`` with ``exp(+i k.R)``: no sign change.
- Prediction layout: per structure, in the order of ``graph_data.npz``,
  ``Hon`` rows then ``Hoff`` rows, each a flattened ``nao_max x nao_max``
  block (spinless); ``(2, nao_max, nao_max)`` with spin up then down
  (collinear spin, ``spin_colinear``); or, with spin-orbit coupling
  (``soc_switch``), the real parts of all ``Hon`` and ``Hoff`` blocks
  followed by their imaginary parts, each ``2 nao_max x 2 nao_max`` in the
  layout ``[[uu, ud], [du, dd]]``.
- Orbitals: every atom is padded to ``nao_max`` orbitals. The real orbitals
  are those with a non-zero diagonal of the on-site overlap ``Son``
  (normalized atomic orbitals have ``S_ii = 1``; padding has 0); HamGNN's
  own tables (``basis_def``) can be passed instead.
- Spin-flip blocks with spin-orbit coupling: HamGNN does not store the
  Hermitian completion. Its graphs (``graph_data_gen.py``: real ``ud`` and
  ``du`` both ``Hks[2]``, imaginary ``du = -(Hks[3] + iHks[2])``) and its
  Uni-HamGNN output (``hamgnn_output.py``, ``so3`` basis: ``ud = xi L_y +
  i xi L_x``, ``du = xi L_y - i xi L_x`` with antisymmetric ``xi L``) both
  satisfy ``du = -ud^+`` on every bond, so ``[[uu, ud], [du, dd]]`` is not
  Hermitian. HamGNN's ``band_cal.py`` never sees this: ``torch.linalg.eigh``
  reads the lower triangle, which (S being spin-diagonal) holds only ``du``.
  The reader does the same: it keeps ``du`` and rebuilds
  ``ud(R) = du(-R)^+``. Flipping the sign of both spin-flip blocks would be a
  spin rotation (bands and ``SZ`` unchanged, ``SX`` and ``SY`` reversed), so
  the choice follows HamGNN's own band structure.
- Overlap: ``Son`` and ``Soff`` are spinless; a spinful model uses
  ``S tensor I_2``. A graph without ``Son`` gives an orthogonal model; its
  real orbitals must then come from ``basis_def``.
- Units: Hartree and Bohr in ``graph_data.npz`` and in the prediction for
  the OpenMX, ABACUS and SIESTA/HONPAS interfaces (HamGNN converts SIESTA's
  Rydberg). The reader converts to eV and Angstrom.
- Energy reference: the absolute scale of the DFT code. The files hold no
  Fermi energy; give it as ``fermi_energy`` (eV) in the run file.

Pitfall: the prediction file holds no structure index. The rows of each
structure are found from the block counts of ``graph_data.npz``, in its
order, as HamGNN's ``band_cal.py`` does; a different ``graph_data.npz``
with the same total size would be cut wrongly. The size check catches a
wrong ``nao_max`` or spin mode, not a reordered file.

Warning: the orbital order inside an atom is HamGNN's (OpenMX or ABACUS
order). LSQuant needs no angular information, so it is kept as it is; do
not use the file for orbital-resolved projections without the DFT code's
table.

To-do: check the export on a real HamGNN prediction with spin-orbit
coupling against the DFT Hamiltonian of the same structure (``prediction
= None``), for the OpenMX and the ABACUS interfaces.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

FORMAT = "lsquant-hamgnn-1"
HARTREE_EV = 27.211386245988
BOHR_ANGSTROM = 0.529177210903
SPIN_MODES = ("none", "collinear", "soc")
DFT_CODES = ("openmx", "abacus", "siesta")

#: Shells of HamGNN's ``nao_max`` padding, one letter per shell (s: 1, p: 3,
#: d: 5, f: 7 orbitals), from ``DFT_interfaces/openmx/utils.py`` (``basis_def_14``,
#: ``_19``, ``_26``; SIESTA uses the 19 layout with another m order) and
#: ``DFT_interfaces/abacus/utils.py`` (``basis_def_13_abacus`` ... ``_40_abacus``).
PADDED_SHELLS = {
    ("openmx", 14): "sssppd", ("openmx", 19): "sssppdd", ("openmx", 26): "sssppddf",
    ("siesta", 19): "sssppdd",
    ("abacus", 13): "ssppd", ("abacus", 15): "ssssppd", ("abacus", 27): "ssssppddf",
    ("abacus", 40): "sssspppp" + "ddff",
}
_SHELL_L = {"s": 0, "p": 1, "d": 2, "f": 3}


def padded_l(dft_code, nao_max):
    """Angular momentum of every index of the ``nao_max`` padding, or ``None`` if unknown."""
    shells = PADDED_SHELLS.get((dft_code, int(nao_max)))
    if shells is None:
        return None
    out = np.concatenate([[_SHELL_L[x]] * (2 * _SHELL_L[x] + 1) for x in shells]).astype(int)
    assert len(out) == nao_max
    return out


def _array(value):
    """A numpy array from a torch tensor or an array-like."""
    if hasattr(value, "detach"):
        value = value.detach().cpu().numpy()
    elif hasattr(value, "numpy"):
        value = value.numpy()
    return np.asarray(value)


def _graphs(graph_data):
    """Structures of ``graph_data.npz`` in HamGNN's order."""
    with np.load(graph_data, allow_pickle=True) as f:
        graph = f["graph"].item()
    return list(graph.values()) if isinstance(graph, dict) else list(graph)


def _block_counts(graph):
    """``(atoms, edges)`` of one structure."""
    return len(_array(graph.z).reshape(-1)), int(_array(graph.edge_index).shape[1])


def _rows_per_structure(graph, spin):
    natoms, nedges = _block_counts(graph)
    return 2 * (natoms + nedges) if spin == "soc" else natoms + nedges


def _split_prediction(prediction, graphs, structure, nao_max, spin):
    """Return the ``(n_blocks, w, w)`` blocks of one structure, complex, padded.

    ``w`` is ``nao_max`` (spinless) or ``2 nao_max`` in the
    ``[[uu, ud], [du, dd]]`` layout (collinear spin and spin-orbit).
    """
    width = {"none": nao_max ** 2, "collinear": 2 * nao_max ** 2, "soc": 4 * nao_max ** 2}[spin]
    prediction = np.asarray(prediction)
    rows = [_rows_per_structure(g, spin) for g in graphs]
    expected = sum(rows) * width
    if prediction.size != expected:
        raise ValueError(
            f"the prediction has {prediction.size} numbers, but {len(graphs)} structures with "
            f"nao_max={nao_max} and spin={spin!r} need {expected}; check nao_max, the spin mode "
            "(soc_switch, spin_colinear) and that the prediction was run on this graph_data.npz")
    prediction = prediction.reshape(-1, width)
    start = sum(rows[:structure])
    mine = prediction[start:start + rows[structure]]
    n = nao_max
    if spin == "none":
        return mine.reshape(-1, n, n).astype(complex)
    if spin == "collinear":
        both = mine.reshape(-1, 2, n, n)
        out = np.zeros((len(mine), 2 * n, 2 * n), complex)
        out[:, :n, :n] = both[:, 0]
        out[:, n:, n:] = both[:, 1]
        return out
    real, imag = np.split(mine.reshape(-1, 2 * n, 2 * n), 2, axis=0)
    return real + 1j * imag


def _target_blocks(graph, nao_max, spin):
    """The DFT blocks stored in ``graph_data.npz`` (same layout as a prediction)."""
    if spin == "soc":
        flat = np.concatenate([_array(graph.Hon), _array(graph.Hoff), _array(graph.iHon), _array(graph.iHoff)])
    else:
        flat = np.concatenate([_array(graph.Hon), _array(graph.Hoff)])
    return _split_prediction(flat, [graph], 0, nao_max, spin)


def _orbital_masks(graph, nao_max, basis_def):
    """Indices of the real orbitals of every atom inside its ``nao_max`` padding."""
    z = _array(graph.z).reshape(-1).astype(int)
    if basis_def is not None:
        return [np.asarray(basis_def[int(a)], int) for a in z]
    if getattr(graph, "Son", None) is None:
        raise ValueError("graph_data.npz holds no overlap (Son), so the real orbitals cannot be found "
                         "from S; pass basis_def (HamGNN's table for your nao_max)")
    son = _array(graph.Son).reshape(len(z), nao_max, nao_max)
    masks = [np.flatnonzero(np.abs(np.diagonal(s)) > 0.5) for s in son]
    for a, m in enumerate(masks):
        if len(m) == 0:
            raise ValueError(f"atom {a} (Z={z[a]}) has no orbital with S_ii > 0.5; pass basis_def")
    return masks


def write_hamgnn_model(graph_data, prediction, path, *, nao_max, spin="none", structure=0,
                       basis_def=None, dft_code=None):
    """Write one HamGNN structure and its Hamiltonian to ``path`` (.npz) for LSQuant.

    Run in the Python environment of HamGNN (``torch`` and
    ``torch_geometric`` are needed to open ``graph_data.npz``).

    - ``graph_data``: the ``graph_data.npz`` the prediction was run on.
    - ``prediction``: ``prediction_hamiltonian.npy`` of that run, or
      ``None`` for the DFT Hamiltonian stored in ``graph_data.npz`` (to
      check the pipeline against DFT).
    - ``nao_max``: as in the HamGNN config (14, 19 or 26 for OpenMX; 27 or
      40 for ABACUS).
    - ``spin``: ``"none"``, ``"collinear"`` (HamGNN ``spin_colinear``) or
      ``"soc"`` (HamGNN ``soc_switch``).
    - ``structure``: index of the structure in ``graph_data.npz``.
    - ``basis_def``: optional ``{Z: orbital indices}``, HamGNN's table of
      real orbitals; by default the orbitals with ``S_ii > 0.5``.
    - ``dft_code``: ``"openmx"``, ``"abacus"`` or ``"siesta"``, the code of
      the training data. It fixes the atomic-orbital basis; with it the file
      records the angular momentum ``l`` of every orbital (from HamGNN's
      padding tables), which identifies the basis in ``basis.json``.

    The file stores the real orbitals only, in Hartree and Bohr, as HamGNN
    holds them; :func:`load_hamgnn_model` converts.
    """
    if spin not in SPIN_MODES:
        raise ValueError(f"spin must be one of {SPIN_MODES}, got {spin!r}")
    if dft_code is not None and dft_code not in DFT_CODES:
        raise ValueError(f"dft_code must be one of {DFT_CODES}, got {dft_code!r}")
    graphs = _graphs(graph_data)
    if not 0 <= structure < len(graphs):
        raise ValueError(f"structure {structure} is out of range ({len(graphs)} structures)")
    graph = graphs[structure]
    if prediction is None:
        blocks = _target_blocks(graph, nao_max, spin)
        source = "dft"
    else:
        blocks = _split_prediction(np.load(prediction), graphs, structure, nao_max, spin)
        source = "prediction"
    natoms, nedges = _block_counts(graph)
    edge_index = _array(graph.edge_index).astype(int)
    cells = np.concatenate([np.zeros((natoms, 3), int), _array(graph.cell_shift).reshape(-1, 3).astype(int)])
    atom_i = np.concatenate([np.arange(natoms), edge_index[0]])
    atom_j = np.concatenate([np.arange(natoms), edge_index[1]])
    has_overlap = getattr(graph, "Son", None) is not None and getattr(graph, "Soff", None) is not None
    overlap = (np.concatenate([_array(graph.Son), _array(graph.Soff)]).reshape(-1, nao_max, nao_max)
               if has_overlap else None)
    masks = _orbital_masks(graph, nao_max, basis_def)

    def spin_index(m):
        return np.concatenate([m, m + nao_max]) if spin != "none" else m

    h_entries, s_entries = [], []
    for b in range(natoms + nedges):
        mi, mj = masks[atom_i[b]], masks[atom_j[b]]
        h_entries.append(blocks[b][np.ix_(spin_index(mi), spin_index(mj))].ravel())
        if has_overlap:
            s_entries.append(overlap[b][np.ix_(mi, mj)].ravel())
    ell = padded_l(dft_code, nao_max)
    orbital_l = (np.concatenate([ell[m] for m in masks]) if ell is not None
                 else np.full(sum(len(m) for m in masks), -1)).astype(int)
    np.savez(
        path, dft_code=np.array(dft_code or "unknown"), orbital_l=orbital_l, format=np.array(FORMAT), source=np.array(source), spin=np.array(spin),
        energy_unit=np.array("Hartree"), length_unit=np.array("Bohr"),
        lattice=_array(graph.cell).reshape(3, 3).astype(float),
        positions=_array(graph.pos).reshape(-1, 3).astype(float),
        atomic_numbers=_array(graph.z).reshape(-1).astype(int),
        orbitals_per_atom=np.array([len(m) for m in masks], int),
        orbital_index=np.concatenate(masks).astype(int), nao_max=np.array(int(nao_max)),
        cells=cells, atom_i=atom_i, atom_j=atom_j,
        hamiltonian=np.concatenate(h_entries),
        **({"overlap": np.concatenate(s_entries)} if has_overlap else {}),
    )


def _blocks(entries, cells, atom_i, atom_j, sizes):
    """Yield ``(R, i, j, block)`` from the flat entries of the export file."""
    offset = 0
    for R, i, j in zip(cells, atom_i, atom_j):
        shape = (sizes[i], sizes[j])
        n = shape[0] * shape[1]
        if offset + n > len(entries):
            raise ValueError("the export file is truncated (fewer entries than blocks)")
        yield R, i, j, entries[offset:offset + n].reshape(shape)
        offset += n
    if offset != len(entries):
        raise ValueError("the export file has more entries than its blocks")


def _complete_spin_flip(r_int, hamiltonian):
    """Rebuild the up-down blocks as ``ud(R) = du(-R)^+`` (interleaved spinors).

    See the module docstring: HamGNN's spin-flip blocks obey ``du = -ud^+``,
    and its band code uses ``du``. Every ``R`` needs its ``-R`` (HamGNN stores
    both directions of every bond).
    """
    lookup = {tuple(r): n for n, r in enumerate(np.asarray(r_int, int))}
    out = np.array(hamiltonian, complex)
    for n, r in enumerate(np.asarray(r_int, int)):
        partner = lookup.get(tuple(-r))
        if partner is None:
            raise ValueError(f"H(R) has R = {tuple(r)} but not -R; the spin-flip blocks cannot be completed")
        out[n, 0::2, 1::2] = hamiltonian[partner][1::2, 0::2].conj().T
    return out


def load_hamgnn_model(path, *, fermi_energy=0.0):
    """Read a file of :func:`write_hamgnn_model` as an LSQuant model.

    Energies become eV and lengths Angstrom; ``fermi_energy`` (eV, on the
    scale of the DFT code that HamGNN learned) is removed as
    ``H - E_F S``. The model is non-orthogonal
    (:class:`~lsquant.models.siesta.SiestaModel`), or orthogonal when the
    file has no overlap; with collinear spin or spin-orbit coupling it is
    spinful, with interleaved spinors. The file gives no velocity and no spin
    operator: the metric basis builds the velocity from dH/dk and dS/dk
    with the orbitals on their atoms, and the spin from the Pauli matrices
    on the spinor labels.
    """
    from .lcao import assemble_blocks, build_lcao_model

    path = os.path.abspath(os.path.expanduser(str(path)))
    with np.load(path, allow_pickle=False) as f:
        data = {k: f[k] for k in f.files}
    if str(data.get("format", "")) != FORMAT:
        raise ValueError(f"{path} is not a HamGNN export of LSQuant (format {FORMAT}); "
                         "write it with write_hamgnn_model")
    spin = str(data["spin"])
    if spin not in SPIN_MODES:
        raise ValueError(f"{path}: unknown spin mode {spin!r}")
    spinful = spin != "none"
    norb = data["orbitals_per_atom"].astype(int)
    cells, atom_i, atom_j = data["cells"].astype(int), data["atom_i"].astype(int), data["atom_j"].astype(int)
    h_sizes = norb * (2 if spinful else 1)
    h_pairs = ((R, i, j, b * HARTREE_EV)
               for R, i, j, b in _blocks(data["hamiltonian"], cells, atom_i, atom_j, h_sizes))
    r_int, hamiltonian = assemble_blocks(h_pairs, norb, spinful=spinful, name="H")
    if spin == "soc":
        hamiltonian = _complete_spin_flip(r_int, hamiltonian)
    s_r_int = overlap = None
    if "overlap" in data:
        s_r_int, overlap = assemble_blocks(_blocks(data["overlap"], cells, atom_i, atom_j, norb), norb,
                                           spinful=False, name="S")
    meta = dict(source="hamgnn", file=path, hamiltonian_source=str(data["source"]), spin=spin,
                atomic_numbers=[int(z) for z in data["atomic_numbers"]], nao_max=int(data["nao_max"]),
                dft_code=str(data.get("dft_code", "unknown")),
                orbital_l=[int(x) for x in data.get("orbital_l", [])],
                orbital_positions="atom-centred: every orbital on its atom")
    return build_lcao_model(data["lattice"] * BOHR_ANGSTROM, data["positions"] * BOHR_ANGSTROM, norb,
                            r_int, hamiltonian, overlap, spinful=spinful,
                            seed_name=os.path.splitext(os.path.basename(path))[0], metadata=meta,
                            fermi_energy=float(fermi_energy), overlap_r_int=s_r_int)


def main(argv=None):
    """Command line of :func:`write_hamgnn_model` (run in the HamGNN environment)."""
    parser = argparse.ArgumentParser(description="Export a HamGNN Hamiltonian for LSQuant")
    parser.add_argument("graph_data", help="graph_data.npz of the prediction run")
    parser.add_argument("prediction", help="prediction_hamiltonian.npy, or 'dft' for the DFT Hamiltonian")
    parser.add_argument("output", help="output .npz for LSQuant")
    parser.add_argument("--nao-max", type=int, required=True)
    parser.add_argument("--spin", choices=SPIN_MODES, default="none")
    parser.add_argument("--structure", type=int, default=0)
    parser.add_argument("--dft-code", choices=DFT_CODES, default=None,
                        help="DFT code of the training data (fixes the orbital basis)")
    args = parser.parse_args(argv)
    write_hamgnn_model(args.graph_data, None if args.prediction == "dft" else args.prediction, args.output,
                       nao_max=args.nao_max, spin=args.spin, structure=args.structure, dft_code=args.dft_code)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    sys.exit(main())
