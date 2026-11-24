"""Pin diagnostics without changing legacy moments, curves, or metric gates.

Small CSR fixtures and synthetic spectral tensors avoid external readers.
The metric failure is injected so advice does not depend on BLAS roundoff.
"""

import json
import math
import os
import sys
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "python"))
pytest.importorskip("lsquant", reason="lsquant._core not built")
from lsquant import _core
from lsquant.runfile import execute_runfile, RunFile
from lsquant.reconstruct import dos_curve, _truncate_and_damp_2d, time_correlation_curve


def run_spec(tmp_path, states=1, legacy=True):
    operators = tmp_path / "operators"
    operators.mkdir(exist_ok=True)
    (operators / "toy.HAM.CSR").write_text("2 2\n-0.5 0 0.5 0\n0 1\n0 1 2\n")
    return dict(model=dict(kind="csr", dir=str(operators), label="toy"),
                basis=dict(name="real"), metric=dict(policy="orthogonal"),
                word="Tr[delta(E-H)]",
                expansion=dict(moments=4, window=[-1, 1], trace="exact", states=states),
                reconstruction=dict(broadening_meV=1, formulas=["dos"]),
                output=dict(legacy_names=legacy, legacy_layout=True))


@pytest.mark.parametrize("legacy", [True, False])
def test_partial_exact_warns_and_records_coverage(tmp_path, legacy):
    spec = run_spec(tmp_path, legacy=legacy)
    with pytest.warns(RuntimeWarning, match="first 1 of 2 basis functions, not the trace"):
        result = execute_runfile(spec, out_dir_override=str(tmp_path / "out"), repo_root=REPO)
    receipt = json.loads((tmp_path / "out" / "run.resolved.json").read_text())
    assert receipt["expansion"]["trace_coverage"] == dict(states=1, dimension=2)
    RunFile.from_dict(receipt)  # output-only provenance must remain replayable
    entry = receipt["receipt"]["curve_resolution"][0]
    assert entry["requested_broadening_meV"] == 1
    assert entry["effective_broadening_meV"] == pytest.approx([math.pi * 250])
    if legacy:
        assert any("stateexact" in os.path.basename(p) for p in result["moment_files"])
    else:
        curves = json.loads((tmp_path / "out" / "curves.json").read_text())
        assert next(iter(curves.values()))["effective_broadening_meV"] == pytest.approx([math.pi * 250])


def synthetic_moments(rank, td=False):
    dims = [4] if rank == 1 else [4, 3]
    data = np.ones(dims, dtype=complex)
    return SimpleNamespace(rank=rank, dims=dims, data=data, system_size=2,
                           window=_core.SpectralWindow(-1, 1),
                           number_of_time_steps=3 if td else 0, time_step=0.1)


@pytest.mark.parametrize("route", ["dos", "rank2", "time"])
def test_resolution_cap_warns_without_changing_arrays(route):
    mu = synthetic_moments(1 if route == "dos" else 2, td=route == "time")
    def reconstruct(broadening=1):
        if route == "dos":
            return dos_curve(mu, broadening)
        if route == "rank2":
            return _truncate_and_damp_2d(mu, broadening, broadening)
        return time_correlation_curve(mu, broadening, 0)
    with pytest.warns(RuntimeWarning, match="Requested broadening 1 meV.*effective broadening"):
        actual = reconstruct()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        reference = reconstruct(2)
    for a, b in zip(actual, reference):
        assert np.array_equal(a, b)
    if route == "dos":
        assert len(actual[0]) == 30 * 4
    elif route == "rank2":
        assert actual[1:] == (4, 3)
    else:
        assert len(actual[0]) == 3


def test_no_cap_or_partial_warning_for_full_trace(tmp_path):
    spec = run_spec(tmp_path, states=2)
    spec["reconstruction"]["broadening_meV"] = 1000
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        execute_runfile(spec, out_dir_override=str(tmp_path / "out"), repo_root=REPO)
    assert not [w for w in caught if w.category is RuntimeWarning]
    receipt = json.loads((tmp_path / "out" / "run.resolved.json").read_text())
    assert "trace_coverage" not in receipt["expansion"]


@pytest.mark.parametrize("applier", ["exact", "chebyshev"])
def test_metric_failure_advice_preserves_gate(monkeypatch, applier):
    pytest.importorskip("scipy.sparse")
    from lsquant.models.siesta import SiestaModel
    from lsquant.disorder_bases import RealSpaceMetricBasis
    overlap = np.array([[1, 0.999], [0.999, 1]], dtype=complex)
    model = SiestaModel(np.eye(3), np.zeros((2, 3)), np.zeros((1, 3), dtype=int),
                        np.diag([-0.001, 0.001])[None], overlap[None])
    monkeypatch.setattr(RealSpaceMetricBasis, "_metric_diagnostics",
                        lambda self: dict(inverse_probe=2e-10, inverse_sqrt_probe=1e-10))
    with pytest.raises(ValueError) as exc:
        RealSpaceMetricBasis(model, [1, 1, 1], applier=applier, velocity_mode="bare", tol=1e-13)
    message = str(exc.value)
    if applier == "exact":
        assert "rounding floor of an ill-conditioned overlap" in message
        assert "estimated cond(S)=1999" in message
        assert "raise metric.tol" in message and "2e-10" in message
        assert "increase Chebyshev" not in message
    else:
        assert "increase Chebyshev order or reduce coeff_tol" in message
    RealSpaceMetricBasis(model, [1, 1, 1], applier=applier, velocity_mode="bare", tol=2e-10)


@pytest.mark.parametrize("basis_class", ["RealBasis", "BlochBasis", "EigenBasis"])
def test_basis_local_exact_warning_distinguishes_block_trace(basis_class):
    from lsquant import bases
    basis = SimpleNamespace(dimension=4, num_wann=2, _disorder_vec=None, seed=12345)
    cls = getattr(bases, basis_class)
    with pytest.warns(RuntimeWarning, match="first 1 of 4 basis functions, not the trace"):
        states = cls.states(basis, "exact", 1)
    assert states.number_of_states == 1
    if basis_class != "RealBasis":
        basis.num_cells = 2
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            cls.states(basis, "exact", 2)
        assert not caught


def test_time_receipt_excludes_time_axis(tmp_path):
    from lsquant.reconstruct import record_resolution
    path = tmp_path / "run.resolved.json"
    path.write_text('{}')
    record_resolution(path, synthetic_moments(2, td=True), 1)
    entry = json.loads(path.read_text())["receipt"]["curve_resolution"][0]
    assert entry["effective_broadening_meV"] == pytest.approx([math.pi * 250])
