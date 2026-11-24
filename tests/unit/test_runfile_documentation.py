"""Keep the maintained schema and compatibility-key diagnostics discoverable.

Validation uses no model files here: dummy paths exercise only the schema,
so these checks stay small and need only the normal compiled Python binding.
"""

import pathlib
import re
import sys
import warnings

import pytest

REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python"))
pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant import runfile  # noqa: E402


def test_reference_documents_every_allowed_key():
    documented = set(re.findall(r"`([^`\n]+)`", (
        REPO / "docs/reference/run-file.md"
    ).read_text()))
    for name, value in vars(runfile).items():
        if not name.startswith("ALLOWED_"):
            continue
        # Generator parameters are a mapping of names to allowed-key sets.
        groups = value.values() if isinstance(value, dict) else (value,)
        for keys in groups:
            assert isinstance(keys, (set, frozenset)), name
            assert keys <= documented, f"{name}: undocumented keys {keys - documented}"


def test_readme_links_to_schema_and_explains_time_and_routes():
    readme = (REPO / "python/lsquant/README.md").read_text()
    assert "../../docs/reference/run-file.md" in readme
    assert "`anticommutator`" in readme
    assert "There is no default `tmax`" in readme
    assert "takes precedence when both are given" in readme
    assert "It does not set KPM curve grids" in readme
    assert "`reconstruction.fft_grid` does not select a grid" in readme
    assert "| `model.kind` |" not in readme


def _raw(reconstruction, model=None):
    return {
        "model": model or {"kind": "csr", "dir": "unused", "label": "tiny"},
        "word": "Tr[VX delta(E-H) VX delta(E-H)]",
        "expansion": {"moments": 8, "window": [-2.0, 2.0]},
        "reconstruction": reconstruction,
    }


@pytest.mark.parametrize("rec, keys", [
    ({"fft_grid": True}, ["fft_grid"]),
    ({"num_points": 17}, ["num_points"]),
    ({"fft_grid": True, "num_points": 17}, ["fft_grid", "num_points"]),
])
def test_inert_grid_keys_warn_once_and_remain_resolved(rec, keys):
    with pytest.warns(RuntimeWarning, match="has no effect") as caught:
        resolved = runfile.RunFile.from_dict(_raw(rec))
    assert len(caught) == 1
    for key in keys:
        assert f"reconstruction.{key}" in str(caught[0].message)
        assert resolved.reconstruction[key] == rec[key]


@pytest.mark.parametrize("rec", [
    {},
    {"fft_grid": False, "num_points": None},
    {"fft_grid": True, "formulas": ["greenwood_fft"]},
    {"fft_grid": True, "formulas": ["KUBO_GREENWOOD_FFT"]},
    {"fft_grid": True, "targets": [{"formulas": ["greenwood_fft"],
                                   "broadening_meV": 10.0}]},
])
def test_effective_or_unset_grid_controls_do_not_warn(rec):
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        runfile.validate_runfile_dict(_raw(rec))


@pytest.mark.parametrize("model", [
    {"kind": "siesta", "dir": "unused", "label": "tiny"},
    {"kind": "wannier", "dir": "unused", "seed": "tiny",
     "params": {"velocity_mode": "bare"}},
    {"kind": "openmx", "file": "unused.scfout"},
])
def test_oracle_grid_count_is_not_reported_inert_for_metric_models(model):
    raw = _raw({"num_points": 17}, model)
    raw["word"] = "Tr[delta(E-H)]"
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        runfile.validate_runfile_dict(raw)


def test_original_wannier_executor_warns_for_grid_count():
    raw = _raw({"num_points": 17},
               {"kind": "wannier", "dir": "unused", "seed": "tiny"})
    with pytest.warns(RuntimeWarning, match="num_points has no effect"):
        runfile.validate_runfile_dict(raw)
