"""Reader regression gate using small synthetic inputs and portable reductions.

Only imported input builders create files, all beneath tmp_path. Numerical
references are generated separately and must never be inferred by a test.
The mutation gate uses an orthogonal model so a deliberately inconsistent
translation cannot be rejected by a Hermiticity check before comparison.
"""
from copy import deepcopy

import numpy as np
import pytest

from reader_fingerprints import (
    CASES, REGENERATE, assert_fingerprints_equal, build_case, fingerprint,
    read_fixture,
)


@pytest.mark.parametrize("case", CASES)
def test_reader_fingerprints(case, tmp_path):
    expected = read_fixture()
    assert case in expected, f"Missing fingerprint for {case}; run: {REGENERATE}"
    assert_fingerprints_equal(fingerprint(build_case(case, tmp_path)), expected[case])


def test_fingerprint_detects_hopping_and_translation(tmp_path):
    """Require each independent perturbation to fail the actual comparator.

    A diagonal hopping increment is Hermitian. A one-cell R translation updates
    both integer and Cartesian vectors, so it represents a lattice translation
    rather than merely a stale derived array. Deep copies isolate the changes;
    the unchanged copy must pass before either negative control is meaningful.
    """
    model = build_case("wannier-default", tmp_path)
    expected = fingerprint(model)
    assert_fingerprints_equal(fingerprint(deepcopy(model)), expected)

    changed = deepcopy(model)
    changed.hoppings[0, 0, 0] += 1e-8
    with pytest.raises(AssertionError):
        assert_fingerprints_equal(fingerprint(changed), expected)

    changed = deepcopy(model)
    changed.r_int[0] += np.array([1, 0, 0])
    changed.r_cart = changed.r_int @ changed.uc
    with pytest.raises(AssertionError):
        assert_fingerprints_equal(fingerprint(changed), expected)
