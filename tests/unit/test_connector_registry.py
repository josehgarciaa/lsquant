"""The connector registry drives every hand-written model.kind list.

One registry entry per metric-executor connector (models/registry.py) states
the path keys, the allowed params, the velocity modes and the loader; the
run-file validation chain, ``uses_metric_executor``, ``ML_MODEL_KINDS`` and
the executor dispatch are derived from it. These tests pin the derivation
against the run-file validator without touching any model data file.

Checks:
- every registered kind appears in the unknown-kind error (``known_kinds``);
- a model missing its path keys raises the entry's exact error text;
- the entry's allowed params are exactly what the validator accepts;
- ``uses_metric_executor`` agrees with the registry for every kind, and the
  Wannier90 trigger params are the ones the registry declares;
- the ABACUS covariant rule admits ``covariant`` only with ``position_file``;
- the velocity-mode refusal texts stay the ones the readers and the docs
  quote;
- ``runfile.ML_MODEL_KINDS`` keeps its original content, built from the
  registry;
- ``lsquant.models.orthogonal.load_runfile_model`` is still the same
  function as the registry dispatcher.
"""
import os
import sys

import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")

from lsquant.models import registry  # noqa: E402
from lsquant.models.orthogonal import load_runfile_model as orthogonal_load_runfile_model  # noqa: E402
from lsquant.models.registry import CONNECTORS, effective_velocity_modes  # noqa: E402
from lsquant.runfile import ML_MODEL_KINDS, RunFile, uses_metric_executor  # noqa: E402

BASE = dict(word="Tr[delta(E-H)]", expansion={"moments": 8},
            basis={"name": "eigen", "supercell": [2, 2, 1]})

_BAD_KEY = "definitely_not_a_param"


def _model_with_paths(connector, params=None):
    model = {"kind": connector.kind}
    for key in connector.path_keys:
        model[key] = "x"
    if params is not None:
        model["params"] = params
    return model


def _validated_params(connector):
    # Wannier90 validates params only with a trigger param; give it one.
    if connector.metric_trigger is not None:
        return {"velocity_mode": "bare", "oracle": True}
    return {"oracle": True}


def test_every_registered_kind_is_known():
    with pytest.raises(ValueError) as exc:
        RunFile.from_dict(dict(BASE, model={"kind": "no_such_kind"}))
    message = str(exc.value)
    assert message.startswith("Unknown model.kind 'no_such_kind'; expected one of ")
    for kind in CONNECTORS:
        assert repr(kind) in message
    for kind in ("csr", "wannier", "siesta", "generator"):
        assert repr(kind) in message


def test_missing_path_error_texts():
    for connector in CONNECTORS.values():
        with pytest.raises(ValueError) as exc:
            RunFile.from_dict(dict(BASE, model={"kind": connector.kind}))
        assert str(exc.value) == connector.missing_path_error


def test_allowed_params_accepted_and_refused():
    for connector in CONNECTORS.values():
        accepted = _model_with_paths(connector, _validated_params(connector))
        rf = RunFile.from_dict(dict(BASE, model=accepted))
        assert rf.model["kind"] == connector.kind
        params = dict(_validated_params(connector))
        params[_BAD_KEY] = 1
        refused = _model_with_paths(connector, params)
        with pytest.raises(ValueError) as exc:
            RunFile.from_dict(dict(BASE, model=refused))
        assert str(exc.value) == (
            f"Unknown key '{_BAD_KEY}' in run file (model.params ({connector.kind})); "
            f"allowed keys are {sorted(connector.allowed_params)}"
        )


def test_uses_metric_executor_per_kind():
    for connector in CONNECTORS.values():
        model = _model_with_paths(connector)
        expected = connector.metric_trigger is None
        assert uses_metric_executor(model) is expected, connector.kind
    assert uses_metric_executor({"kind": "csr", "dir": "d", "label": "l"}) is False
    assert uses_metric_executor({"kind": "generator", "name": "make_chain"}) is False
    assert uses_metric_executor({"kind": "unknown"}) is False


def test_wannier_trigger_params():
    wannier = CONNECTORS["wannier"]
    assert wannier.metric_trigger == ("velocity_mode", "spin", "fermi_energy")
    model = {"kind": "wannier", "dir": ".", "seed": "s"}
    assert uses_metric_executor(model) is False
    assert uses_metric_executor(dict(model, params={})) is False
    for key in wannier.metric_trigger:
        assert uses_metric_executor(dict(model, params={key: "labels"})) is True
    # A non-Mapping params cannot hold a trigger, as before.
    assert uses_metric_executor(dict(model, params=["velocity_mode"])) is False


def test_wannier_params_unvalidated_without_trigger():
    # Preserved behaviour (a known gap, fixed by the validation unification):
    # a Wannier90 model without a trigger param keeps its params unchecked.
    model = {"kind": "wannier", "dir": ".", "seed": "s", "params": {_BAD_KEY: 1}}
    rf = RunFile.from_dict(dict(BASE, model=model))
    assert rf.model["params"] == {_BAD_KEY: 1}


def test_abacus_covariant_rule():
    abacus = CONNECTORS["abacus"]
    assert abacus.covariant_rule == "param"
    assert abacus.covariant_param == "position_file"
    assert abacus.path_params == ("overlap_file", "stru_file", "orbital_dir", "position_file")
    assert effective_velocity_modes(abacus, {}) == ("bare", "berry_connection")
    assert effective_velocity_modes(abacus, {"position_file": "p"}) == (
        "bare", "berry_connection", "covariant")
    model = {"kind": "abacus", "file": "hr",
             "params": {"overlap_file": "s", "stru_file": "STRU", "velocity_mode": "covariant"}}
    with pytest.raises(ValueError) as exc:
        RunFile.from_dict(dict(BASE, model=model))
    assert str(exc.value) == (
        "model.params.velocity_mode for abacus must be 'bare' or 'berry_connection'; "
        "'covariant' needs model.params.position_file (ABACUS's data-rR-sparse.csr)"
    )
    accepted = dict(model, params=dict(model["params"], position_file="pos"))
    rf = RunFile.from_dict(dict(BASE, model=accepted))
    assert rf.model["params"]["position_file"] == "pos"
    # GPAW admits covariant through its radial quadrature, without any param.
    gpaw = CONNECTORS["gpaw"]
    assert effective_velocity_modes(gpaw, {}) == ("bare", "berry_connection", "covariant")
    # PAOFLOW has no position matrix at all.
    paoflow = CONNECTORS["paoflow"]
    assert effective_velocity_modes(paoflow, {"position_file": "p"}) == (
        "bare", "berry_connection")


def test_velocity_mode_error_texts():
    for connector in CONNECTORS.values():
        params = {"velocity_mode": "bogus"}
        model = _model_with_paths(connector, params)
        with pytest.raises(ValueError) as exc:
            RunFile.from_dict(dict(BASE, model=model))
        expected = connector.velocity_error.format(
            kind=connector.kind,
            modes=list(effective_velocity_modes(connector, params)))
        assert str(exc.value) == expected, connector.kind


def test_spin_fermi_refusals():
    for kind in ("openmx", "pyscf", "paoflow"):
        connector = CONNECTORS[kind]
        for bad in ("spin", "fermi_energy"):
            model = _model_with_paths(connector, {bad: 0.1})
            with pytest.raises(ValueError) as exc:
                RunFile.from_dict(dict(BASE, model=model))
            assert str(exc.value) == connector.refusal
    # Wannier90 accepts both; SIESTA and the H(R)/S(R) kinds refuse spin by
    # key validation (it is not in their allowed params).
    wannier = _model_with_paths(CONNECTORS["wannier"],
                                {"velocity_mode": "bare", "spin": "labels", "fermi_energy": 0.0})
    assert RunFile.from_dict(dict(BASE, model=wannier)).model["kind"] == "wannier"


def test_ml_model_kinds_matches_registry():
    assert ML_MODEL_KINDS == {
        "deeph": ("dir", frozenset({"hamiltonian_file"})),
        "deeptb": ("dir", frozenset()),
        "hamgnn": ("file", frozenset()),
        "threebodytb": ("file", frozenset({"orbitals_per_element"})),
        "abacus": ("file", frozenset({"overlap_file", "stru_file", "orbital_dir",
                                      "orbitals_per_element", "position_file"})),
        "gpaw": ("file", frozenset({"soc", "quadrature_order", "dipole_cache"})),
    }
    assert list(ML_MODEL_KINDS) == [k for k, c in CONNECTORS.items() if c.ml_style]


def test_metric_executor_names_message():
    expected_names = ", ".join(
        k for k, c in CONNECTORS.items() if c.metric_trigger is None)
    with pytest.raises(ValueError) as exc:
        RunFile.from_dict(dict(BASE, model={"kind": "csr", "dir": "d", "label": "l"},
                               reference=True))
    assert str(exc.value) == (
        f"Expected a metric-executor model ({expected_names}, or wannier with "
        "params.velocity_mode or params.spin) when 'reference' is set, got 'csr'"
    )


def test_load_runfile_model_is_reexported():
    assert orthogonal_load_runfile_model is registry.load_runfile_model
    for connector in CONNECTORS.values():
        assert connector.loader is not None, connector.kind
    with pytest.raises(ValueError, match="has no metric executor"):
        orthogonal_load_runfile_model({"kind": "csr", "dir": "d", "label": "l"}, "bare")
