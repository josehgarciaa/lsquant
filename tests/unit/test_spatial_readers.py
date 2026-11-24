"""Exercise small spatial inputs and periodic block assembly contracts.

All file fixtures live under tmp_path. Arrays use at most four atoms and
one orbital per atom, so expected permutations and dictionaries are explicit.
Pitfall: hopping displacements are unwrapped primitive-cell displacements,
whereas neighbour indices wrap around the periodic supercell.
"""

from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant import spatial as s


def test_hr_valid(tmp_path):
    """Read six-column hopping records without discarding imaginary parts."""
    path = tmp_path / "h.dat"
    path.write_text("2 2\n0 0 0 1 2 3\n1 0 1 0 2 -3\n")
    np.testing.assert_array_equal(s.read_hr_file(path),
                                 [[0, 0, 0, 1, 2, 3], [1, 0, 1, 0, 2, -3]])


@pytest.mark.parametrize("text, message", [
    ("", "must start"), ("2\n", "must start"),
    ("two 1\n0 0 0 0 1 0\n", "invalid literal"),
    ("2 2\n0 0 0 1 1 0\n", "expected 2"),
    ("2 1\n0 0 0 1 1\n", "expected 1"),
    ("2 1\n0 0 0 2 1 0\n", "out of range"),
    ("2 1\n0 0 0 1 bad 0\n", "convert"),
])
def test_hr_invalid(tmp_path, text, message):
    """Require actionable ValueErrors for broken headers and hopping rows."""
    path = tmp_path / "h.dat"
    path.write_text(text)
    with pytest.raises(ValueError, match=message):
        s.read_hr_file(path)


@pytest.mark.parametrize("text", [
    "# positions\n\n0 1 2\n3 4 5\n",
    "2\nunit cell\nA 0 1 2\nB 3 4 5\n",
])
def test_positions_valid(tmp_path, text):
    """Accept raw coordinates and tagged XYZ records in the same order."""
    path = tmp_path / "positions.dat"
    path.write_text(text)
    actual = s.read_positions(path)
    np.testing.assert_array_equal(actual, [[0, 1, 2], [3, 4, 5]])
    assert actual.dtype == np.float64
    assert actual.flags.c_contiguous


@pytest.mark.parametrize("text, message", [
    ("# nothing\n", "empty"),
    ("2\ncomment\nA 0 1 2\n", "expected 2 rows"),
    ("0 1 bad\n", "convert"),
    ("two\ncomment\n", "invalid literal"),
])
def test_positions_invalid(tmp_path, text, message):
    """Expose missing records and nonnumeric coordinates as ValueErrors."""
    path = tmp_path / "positions.dat"
    path.write_text(text)
    with pytest.raises(ValueError, match=message):
        s.read_positions(path)


def test_positions_short_row(tmp_path):
    """Require a clear exception instead of an internal column IndexError."""
    path = tmp_path / "positions.dat"
    path.write_text("0 1\n")
    with pytest.raises(ValueError, match="column|coordinate|shape"):
        s.read_positions(path)


def test_positions_arrays():
    """Pad two-dimensional input with a zero third coordinate."""
    np.testing.assert_array_equal(s.read_positions([[1, 2], [3, 4]]),
                                 [[1, 2, 0], [3, 4, 0]])
    for value in ([1, 2], [[1]], [[1, 2, 3, 4]]):
        with pytest.raises(ValueError, match="shape"):
            s.read_positions(value)


def test_atom_maps(tmp_path):
    """Keep uniform and explicit consecutive orbital grouping equivalent."""
    expected = [0, 0, 1, 1]
    path = tmp_path / "map.dat"
    path.write_text("0\n0\n1\n1\n")
    for kwargs in ({"orbitals_per_atom": 2}, {"atom_map": expected},
                   {"atom_map": path}):
        np.testing.assert_array_equal(s.atom_map_from_spec(4, **kwargs), expected)
    for kwargs in ({}, {"orbitals_per_atom": 3}, {"orbitals_per_atom": 0},
                   {"atom_map": [0, 1]}, {"atom_map": [0, 1, 0, 1]},
                   {"atom_map": [-1, 0, 0, 1]},
                   {"orbitals_per_atom": 2, "atom_map": expected}):
        with pytest.raises(ValueError):
            s.atom_map_from_spec(4, **kwargs)
    np.testing.assert_array_equal(s.orbital_permutation(
        np.array([1, 0]), np.array(expected)), [2, 3, 0, 1])
    with pytest.raises(ValueError, match="non-decreasing"):
        s.orbital_permutation(np.array([0, 1]), np.array([1, 0]))


def test_morton_and_domain_limits():
    r"""Check interleaved axis bits and the :math:`[64,16384]` cache clamp."""
    coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0],
                       [0, 0, 1], [2, 0, 0], [1, 1, 1]])
    np.testing.assert_array_equal(s.morton_codes(coords), [0, 1, 2, 4, 8, 7])
    assert s.default_domain_atoms(1) == 10922
    assert s.default_domain_atoms(1, vectors_in_flight=1) == 16384
    assert s.default_domain_atoms(1000, state_block=8) == 64
    assert s.default_domain_atoms(2, state_block=99) == s.default_domain_atoms(2, 8)


def test_spatial_order():
    """Place interior atoms before boundary atoms within each domain."""
    pos = np.array([[0., 0, 0], [1., 0, 0], [2., 0, 0], [3., 0, 0]])
    order, info = s.spatial_permutation(pos, np.array([0, 1, 2]),
                                       np.array([1, 2, 3]), 2)
    np.testing.assert_array_equal(order, [0, 1, 3, 2])
    assert info["num_domains"] == 2
    assert info["boundary_atoms"] == 2


def test_hopping_blocks_and_velocity():
    r"""Check periodic neighbours and :math:`v_x=i[H,X]` dictionaries.

    The two-cell scalar crystal has one positive-cell hop of amplitude two.
    Both cells wrap to the other cell, while the bond length remains three.
    """
    hops = np.array([[1, 0, 0, 0, 2., 0.]])
    pos = np.array([[0., 0, 0], [3., 0, 0]])
    arrays, perm, info = s.blocks_from_hoppings(hops, (2, 1), pos, 1, 64)
    np.testing.assert_array_equal(perm, [0, 1])
    np.testing.assert_array_equal(arrays["nbr"], [1, 0])
    np.testing.assert_array_equal(arrays["dict_values"], [2])
    np.testing.assert_array_equal(arrays["atom_offset"], [0, 1, 2])
    np.testing.assert_array_equal(arrays["pair_offset"], [0, 1, 2])
    assert info["hopping_range"] == 1
    lattice = np.array([[3., 0], [0, 1], [0, 0]])
    velocity = s.velocity_arrays_from_hoppings(hops, arrays, pos[:1], lattice, 0)
    np.testing.assert_array_equal(velocity["dict_values"], [6j])
    assert velocity["nbr"] is arrays["nbr"]
    np.testing.assert_array_equal(arrays["dict_values"], [2])
    doubled = s.hoppings_dictionary(np.concatenate([hops, hops]), 1,
                                    pos[:1], lattice)
    np.testing.assert_array_equal(doubled, [[[4]]])
    operator = s.to_core(arrays, hermitian=False)
    assert operator is not None


def test_csr_blocks_diagonal_override():
    """Preserve explicit onsite values under a reversed orbital permutation."""
    arrays = s.blocks_from_csr(np.array([0, 2, 4]),
                               np.array([0, 1, 0, 1]),
                               np.array([1., 2., 2., 3.]),
                               np.array([0, 1]), np.array([1, 0]))
    np.testing.assert_array_equal(arrays["diagonal"], [3, 1])
    overridden = s.blocks_from_csr(np.array([0, 2, 4]),
                                   np.array([0, 1, 0, 1]),
                                   np.array([1., 2., 2., 3.]),
                                   np.array([0, 1]), np.array([1, 0]),
                                   diagonal=np.array([7., 8.]))
    np.testing.assert_array_equal(overridden["diagonal"], [8, 7])


def test_cuda_adapter_arguments(monkeypatch):
    """Check array conversion without requiring an optional device backend.

    Only the binding constructor is substituted; all adapter conversions run.
    """
    from lsquant import _core
    monkeypatch.setattr(_core, "CudaAtomBlockOperator", lambda **kw: kw,
                        raising=False)
    integer_keys = ("atom_offset", "pair_offset", "nbr", "block_id",
                    "dict_offset", "dict_rows", "dict_cols")
    arrays = {key: [0, 1] for key in integer_keys}
    arrays.update(dict_values=[2.], diagonal=[3.])
    out = s.to_core_cuda(arrays, hermitian=False)
    assert out["hermitian"] is False
    assert out["nbr"].dtype == np.int32
    assert out["atom_offset"].dtype == np.int64
    assert out["dict_values"].dtype == np.complex128
    assert all(value.flags.c_contiguous for key, value in out.items()
               if key != "hermitian")
