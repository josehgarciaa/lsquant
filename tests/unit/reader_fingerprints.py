"""Portable numerical summaries of the small, file-backed reader examples.

The fixture records values in the model's native units (eV and Angstrom),
without paths, timings, byte hashes, or spectral decompositions. Reductions
retain complex real and imaginary parts separately. Work and temporary memory
are linear in the number of array elements, apart from the two Bloch sums.

Warning: regenerate only after reviewing an intended numerical change.
Pitfall: these summaries detect changes, but cannot uniquely identify an array.
SIESTA has no self-generated reader input in the unit suite; its existing
example is deliberately outside this small synthetic-input gate.
"""
from __future__ import annotations

import argparse
import importlib
import json
import os
from pathlib import Path
import sys
import tempfile

if __name__ == "__main__":
    for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[_name] = "1"

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = Path(__file__).with_name("data") / "reader_fingerprints.json"
REGENERATE = "python tests/unit/reader_fingerprints.py --write"
K_FRACTIONAL = np.array([[0.137, 0.219, 0.071], [-0.4, 0.31, 0.0]])
# Each builder is imported, rather than duplicated, from the existing tests.
BUILDERS = {
    "wannier": ("_wannier_files", 21),
    "paoflow": ("_paoflow_file", 22),
    "pyscf": ("_pyscf_file", 23),
    "hamgnn": ("_hamgnn_file", 24),
    "deeptb": ("_deeptb_dir", 25),
    "abacus": ("_abacus_dir", 26),
    "gpaw": ("_gpaw_file", 27),
    "deeph": ("_deeph_dir", 28),
    "threebodytb": ("_threebodytb_file", 29),
    "openmx": ("_openmx_file", 30),
}
CASES = tuple(f"{kind}-{origin}" for kind in BUILDERS
              for origin in ("default", "explicit")) + (
                  "wannier-positions", "pyscf-dipoles", "openmx-dipoles")


def _summary(array):
    """Reduce one array in C order with reproducible seed-zero real weights.

    A real linear functional of a complex array is complex; keeping both parts
    catches imaginary hopping changes. Norms use flattened Euclidean length,
    which is the Frobenius norm for matrix stacks as well as matrices.
    """
    array = np.asarray(array)
    flat = array.astype(np.complex128, copy=False).ravel(order="C")
    total = np.sum(flat)
    projection = np.sum(np.random.default_rng(0).standard_normal(flat.size) * flat)
    return {"shape": list(array.shape),
            "sum": [float(total.real), float(total.imag)],
            "sum_abs": float(np.sum(np.abs(flat))),
            "norm": float(np.sqrt(np.sum(np.abs(flat) ** 2))),
            "linear": [float(projection.real), float(projection.imag)]}


def _arrays(value):
    """Preserve absent fields and nested position/spin arrays explicitly."""
    if value is None:
        return None
    if isinstance(value, dict):
        return {str(key): _arrays(item) for key, item in sorted(value.items())}
    if isinstance(value, (tuple, list)):
        return [_arrays(item) for item in value]
    if isinstance(value, np.ndarray):
        return _summary(value)
    if isinstance(value, (int, float, complex, np.number)):
        return _summary(np.asarray(value))
    # The exact Wannier spin object owns its real-space arrays and mesh.
    return {key: _arrays(item) for key, item in sorted(vars(value).items())
            if isinstance(item, (np.ndarray, tuple, list, dict))}


def fingerprint(model):
    """Return JSON-compatible invariants of a reader's numerical model.

    Both model classes store uc, tau, r_int, r_cart and hoppings. SiestaModel
    additionally stores overlaps, dipoles and optional physical_spin;
    OrthogonalModel stores optional (position R vectors, position matrices)
    and on-site or real-space spin arrays. Optional fields remain explicit.
    Fermi metadata is selected by name so new bookkeeping fields are included
    without admitting transient input paths. Reciprocal points are fractional;
    Hk/Sk receive Cartesian inverse-Angstrom vectors in the atom gauge.
    """
    arrays = {name: _arrays(getattr(model, name, None)) for name in (
        "uc", "tau", "r_int", "r_cart", "hoppings", "overlaps", "dipoles",
        "positions", "spin", "physical_spin")}
    k = K_FRACTIONAL @ (2 * np.pi * np.linalg.inv(model.uc).T)
    arrays["Hk"] = _summary(model.Hk(k))
    arrays["Sk"] = _summary(model.Sk(k))
    metadata = {key: (value.item() if isinstance(value, np.generic) else value)
                for key, value in model.metadata.items() if "fermi" in key.lower()}
    return {"arrays": arrays, "fermi": metadata,
            "num_wann": int(model.num_wann),
            "spinful": getattr(model, "spinful", None)}


def assert_fingerprints_equal(actual, expected, path="model"):
    """Compare structure exactly and reductions at rtol 1e-10.

    Cancellation can make sums and projections vanish even for large arrays.
    Their absolute floor is 1e-12 times the reference Frobenius norm, with a
    unit floor for empty or zero arrays. Never scale by the changed array:
    a large regression must not inflate its own acceptance threshold.
    """
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and actual.keys() == expected.keys(), path
        if "shape" in expected and "norm" in expected:
            assert actual["shape"] == expected["shape"], f"{path}.shape"
            atol = 1e-12 * max(1.0, expected["norm"])
            for key in ("sum", "sum_abs", "norm", "linear"):
                np.testing.assert_allclose(actual[key], expected[key], rtol=1e-10,
                                           atol=atol, equal_nan=False,
                                           err_msg=f"{path}.{key}")
        else:
            for key in expected:
                assert_fingerprints_equal(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(actual) == len(expected), path
        for index, (left, right) in enumerate(zip(actual, expected)):
            assert_fingerprints_equal(left, right, f"{path}[{index}]")
    elif isinstance(expected, float):
        np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-12,
                                   equal_nan=False, err_msg=path)
    else:
        assert actual == expected, path


def build_case(case, tmp_path):
    """Load one existing synthetic input, with isolated files and a fixed seed.

    Optional dependencies retain the skips of the original tests. The explicit
    origin exercises H - E_F S as well as metadata, without any extra input.
    """
    unit_dir = str(Path(__file__).resolve().parent)
    if unit_dir not in sys.path:
        sys.path.insert(0, unit_dir)
    module = importlib.import_module("test_connector_validation")
    kind, origin = case.rsplit("-", 1)
    name, seed = BUILDERS[kind]
    result = getattr(module, name)(tmp_path, np.random.default_rng(seed))
    load = result[0] if isinstance(result, tuple) else result
    if origin == "dipoles":
        if kind == "pyscf":
            from lsquant.models.pyscf_lcao import load_pyscf_model
            return load_pyscf_model(tmp_path / "pyscf.npz", dipoles=True)
        from lsquant.models.openmx import load_openmx_model
        return load_openmx_model(tmp_path / "toy.scfout", dipoles=True)
    if origin == "positions":
        examples = importlib.import_module("test_orthogonal_models")
        examples._write_r_dat(tmp_path, "w", np.random.default_rng(0))
        from lsquant.models.orthogonal import load_wannier_metric_model
        return load_wannier_metric_model(str(tmp_path), "w", dipoles=True)
    return load() if origin == "default" else load(0.85)


def read_fixture():
    """Fail explicitly until numerical references have been generated."""
    assert FIXTURE.is_file(), f"fixture not generated; run: {REGENERATE}"
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


def main():
    """Generate all cases atomically only when every input can be loaded.

    To-do: review and add the resulting JSON to version control after checking
    these tests. Missing optional packages must not create a partial fixture.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", required=True)
    parser.parse_args()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[name] = "1"
    sys.path.insert(0, str(ROOT / "python"))
    import pytest

    results = {}
    with tempfile.TemporaryDirectory(prefix="reader-fingerprints-") as directory:
        for case in CASES:
            tmp_path = Path(directory) / case
            tmp_path.mkdir()
            try:
                results[case] = fingerprint(build_case(case, tmp_path))
            except pytest.skip.Exception as error:
                raise RuntimeError(f"Cannot generate {case}: {error}") from error
    content = json.dumps(results, indent=2, sort_keys=True, allow_nan=False) + "\n"
    FIXTURE.parent.mkdir(parents=True, exist_ok=True)
    temporary = FIXTURE.with_suffix(".json.tmp")
    temporary.write_text(content, encoding="utf-8")
    temporary.replace(FIXTURE)
    print(f"Generated {len(results)} cases in {FIXTURE.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
