"""Registry of the connector kinds that feed the metric executors.

One entry per ``model.kind`` handled by the BlochMetricBasis executors
(``apps/stage6b_siesta.py`` and ``projection/runfile.py``): SIESTA, OpenMX,
PySCF, Wannier90, PAOFLOW, and the H(R)/S(R) Hamiltonians read by
``models/<kind>.py`` (DeepH, DeepTB, HamGNN, ThreeBodyTB, ABACUS, GPAW).

The run file used to hand-write the list of model kinds in about eight places
(a validation chain per kind, the dispatch, ``ML_MODEL_KINDS``, the missing
-path error texts). Every such list is now derived from this registry, so
adding a connector means adding one entry here plus one reader module.

Each entry states, as data:

- ``path_keys``: the ``model`` keys that must hold a path (``dir``+``label``,
  ``dir``+``seed``, ``dir``, or ``file``) and ``path_hint``: the parenthetical
  hint of the missing-path error, where there is one;
- ``model_keys``: the accepted keys of ``model`` itself (``kind`` and
  ``params`` plus the path keys and any key the executors read; a key of
  another kind is refused, audit item B8);
- ``allowed_params``: the accepted keys of ``model.params``;
- ``velocity_modes`` and ``covariant_rule``: the accepted velocity modes,
  where ``always``/``never``/``param`` says when ``covariant`` joins the base
  modes (ABACUS admits it only with its ``position_file`` param);
- ``required_params``: params keys that must hold a value (the ABACUS S(R)
  file and STRU), with ``required_param_hint`` as the parenthetical why;
- ``path_params``: params keys that hold paths resolved like ``model.file``
  (today the four ABACUS keys);
- ``refusal``: text refusing ``model.params.spin`` when the reader stores its
  own spinors (OpenMX, PySCF, PAOFLOW);
- ``metric_trigger``: params keys that route the kind to the metric executor
  (only Wannier90: ``velocity_mode``/``spin``/``fermi_energy``; every other
  connector always runs there). Without a trigger key a Wannier90 model keeps
  the plain CSR executor and its params may not use the executor-only keys
  (``oracle``, ``oracle_degeneracy_tolerance``, ``bloch_memory_mb``);
- ``keep_metric_policy``: whether the executor keeps ``metric.policy`` as
  given (OpenMX, PySCF) or is switched to the exact Lowdin route;
- ``loader``: a callable ``(model, velocity_mode, phase) -> model`` that
  imports the reader module lazily, so optional dependencies (sisl, h5py,
  ...) never load at package import.

Validation. :func:`validate_connector_params` applies ONE rule set to every
connector (audit items B3, B4, B5): the params must be a JSON object or null;
unknown keys are refused; every numeric param is validated for every kind
that accepts it, with a ``ValueError`` naming the key on a wrong type (never
a ``TypeError``); one wording per message. ``fermi_energy`` (eV, on the scale
of the DFT code's output) selects the energy origin for every connector:
omitted or null means the reference energy of the file when the file has one,
else 0.0; an explicit value replaces it (see each reader). The readers record
``fermi_energy_subtracted_eV`` (the origin used) and ``fermi_energy_source``
(``"file"``, ``"run file"`` or ``"none"``) in the model metadata.

Pitfall: the Fermi-level metadata keys of the readers predate this policy
(``printed_fermi_level_eV``, ``gpaw_fermi_level_eV``,
``fermi_energy_from_file_eV``, ...). They stay, unchanged: the two new keys
are added next to them, not renamed over them.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

import numpy as np

#: Base ``model.params`` keys accepted by the four orthogonal models
#: (wannier, openmx, pyscf, paoflow): ``fermi_energy`` selects the energy
#: origin (see the module docstring); ``spin`` is Wannier90 only.
METRIC_BASE_PARAMS = frozenset(
    {"velocity_mode", "spin", "fermi_energy", "oracle",
     "oracle_degeneracy_tolerance", "bloch_memory_mb"}
)

#: Base ``model.params`` keys accepted by the H(R)/S(R) Hamiltonians
#: (deeph, deeptb, hamgnn, threebodytb, abacus, gpaw).
ML_BASE_PARAMS = frozenset(
    {"velocity_mode", "fermi_energy", "oracle",
     "oracle_degeneracy_tolerance", "bloch_memory_mb"}
)

#: Params that only make sense with the metric executor: refused on a
#: trigger kind (wannier) when no trigger param selects the executor.
_TRIGGER_ONLY_PARAMS = ("oracle", "oracle_degeneracy_tolerance", "bloch_memory_mb")


@dataclass(frozen=True)
class ModelConnector:
    """One ``model.kind`` of the metric executors, as data (see module docstring)."""

    kind: str
    path_keys: tuple[str, ...]
    path_hint: str
    allowed_params: frozenset[str]
    velocity_modes: tuple[str, ...]
    model_keys: frozenset[str]
    loader: Callable[[Mapping[str, Any], str, Any], Any] | None
    covariant_rule: str = "never"           # "always" | "never" | "param"
    covariant_param: str = ""               # params key that admits covariant
    path_params: tuple[str, ...] = ()       # params keys resolved like model.file
    refusal: str | None = None              # text refusing model.params.spin
    metric_trigger: tuple[str, ...] | None = None   # None: always metric executor
    keep_metric_policy: bool = False        # True: keep metric.policy as given
    ml_style: bool = False                  # True: H(R)/S(R) input like ML kinds
    spin_values: tuple[str, ...] = ()
    required_params: tuple[str, ...] = ()
    required_param_hint: str = ""

    @property
    def missing_path_error(self) -> str:
        """The exact error when a ``path_keys`` value is missing (one wording)."""
        keys = " and ".join(repr(k) for k in self.path_keys)
        message = f"model.kind='{self.kind}' requires {keys}"
        return f"{message} {self.path_hint}" if self.path_hint else message


def _load_siesta(model, velocity_mode, phase):
    from .siesta import load_siesta_model
    params = model.get("params") or {}
    return load_siesta_model(
        model["dir"], model["label"], dipoles=velocity_mode == "covariant",
        quadrature_order=params.get("quadrature_order", 24),
        quadrature_tolerance=params.get("quadrature_tolerance", 1e-5),
        quadrature_max_order=params.get("quadrature_max_order", 384),
        fermi_energy=params.get("fermi_energy"),
        phase=phase)


def _load_wannier(model, velocity_mode, phase):
    from .orthogonal import load_wannier_metric_model
    params = model.get("params") or {}
    return load_wannier_metric_model(
        model["dir"], model["seed"], spin=params.get("spin", "labels"),
        positions=True if velocity_mode == "covariant" else False,
        fermi_energy=params.get("fermi_energy"))


def _load_openmx(model, velocity_mode, phase):
    from .openmx import load_openmx_model
    params = model.get("params") or {}
    return load_openmx_model(model["file"], dipoles=velocity_mode == "covariant",
                             fermi_energy=params.get("fermi_energy"))


def _load_pyscf(model, velocity_mode, phase):
    from .pyscf_lcao import load_pyscf_model
    params = model.get("params") or {}
    return load_pyscf_model(model["file"], connection=velocity_mode == "covariant",
                            fermi_energy=params.get("fermi_energy"))


def _load_paoflow(model, velocity_mode, phase):
    if velocity_mode == "covariant":
        raise ValueError("a PAOFLOW model has no position matrix; use bare or berry_connection")
    from .orthogonal import load_paoflow_model
    params = model.get("params") or {}
    return load_paoflow_model(model["file"], fermi_energy=params.get("fermi_energy"))


def _load_deeph(model, velocity_mode, phase):
    from .deeph import load_deeph_model
    params = model.get("params") or {}
    return load_deeph_model(model["dir"], hamiltonian_file=params.get("hamiltonian_file"),
                            fermi_energy=_optional_float(params.get("fermi_energy")))


def _load_deeptb(model, velocity_mode, phase):
    from .deeptb import load_deeptb_model
    params = model.get("params") or {}
    return load_deeptb_model(model["dir"], fermi_energy=params.get("fermi_energy"))


def _load_hamgnn(model, velocity_mode, phase):
    from .hamgnn import load_hamgnn_model
    params = model.get("params") or {}
    return load_hamgnn_model(model["file"], fermi_energy=params.get("fermi_energy"))


def _load_threebodytb(model, velocity_mode, phase):
    from .threebodytb import load_threebodytb_model
    params = model.get("params") or {}
    return load_threebodytb_model(
        model["file"], fermi_energy=_optional_float(params.get("fermi_energy")),
        orbitals_per_element=params.get("orbitals_per_element"))


def _load_abacus(model, velocity_mode, phase):
    from .abacus import load_abacus_model
    params = model.get("params") or {}
    return load_abacus_model(
        model["file"], overlap_file=params["overlap_file"], stru_file=params["stru_file"],
        orbital_dir=params.get("orbital_dir"),
        orbitals_per_element=params.get("orbitals_per_element"),
        fermi_energy=params.get("fermi_energy"),
        position_file=params.get("position_file") if velocity_mode == "covariant" else None)


def _load_gpaw(model, velocity_mode, phase):
    from .gpaw import load_gpaw_model
    params = model.get("params") or {}
    return load_gpaw_model(
        model["file"], soc=bool(params.get("soc", True)),
        fermi_energy=_optional_float(params.get("fermi_energy")),
        dipoles=velocity_mode == "covariant",
        quadrature_order=params.get("quadrature_order", 48),
        dipole_cache=params.get("dipole_cache"))


#: Every kind accepts ``kind``, ``params`` and ``operators`` (the word
#: compiler reads ``model.operators`` for any kind).
_BASE_MODEL_KEYS = frozenset({"kind", "params", "operators"})

#: ``{kind: connector}`` in the order in which run-file messages list the
#: kinds; ``ML_MODEL_KINDS`` and ``known_kinds`` are slices of this order.
CONNECTORS: dict[str, ModelConnector] = {
    c.kind: c for c in (
        ModelConnector(
            "siesta",
            ("dir", "label"),
            "",
            frozenset({"velocity_mode", "fermi_energy", "oracle", "oracle_degeneracy_tolerance",
                       "bloch_memory_mb", "quadrature_order", "quadrature_max_order",
                       "quadrature_tolerance", "nonorthogonal_disorder"}),
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"dir", "label"},
            _load_siesta,
            covariant_rule="always",
        ),
        ModelConnector(
            "openmx",
            ("file",),
            "(the .scfout file)",
            METRIC_BASE_PARAMS,
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"file"},
            _load_openmx,
            covariant_rule="always",
            refusal="an OpenMX model carries its spinors; remove model.params.spin",
            keep_metric_policy=True,
        ),
        ModelConnector(
            "pyscf",
            ("file",),
            "(see write_pyscf_model)",
            METRIC_BASE_PARAMS,
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"file"},
            _load_pyscf,
            covariant_rule="always",
            refusal="a PySCF model carries its spinors; remove model.params.spin",
            keep_metric_policy=True,
        ),
        ModelConnector(
            "wannier",
            ("dir", "seed"),
            "",
            METRIC_BASE_PARAMS,
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"dir", "seed", "csr_dir"},
            _load_wannier,
            covariant_rule="always",
            metric_trigger=("velocity_mode", "spin", "fermi_energy"),
            spin_values=("labels", "exact"),
        ),
        ModelConnector(
            "paoflow",
            ("file",),
            "(see write_paoflow_model)",
            METRIC_BASE_PARAMS,
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"file"},
            _load_paoflow,
            refusal="a PAOFLOW model carries its own on-site spin; remove model.params.spin",
        ),
        ModelConnector(
            "deeph",
            ("dir",),
            "",
            ML_BASE_PARAMS | {"hamiltonian_file"},
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"dir"},
            _load_deeph,
            ml_style=True,
        ),
        ModelConnector(
            "deeptb",
            ("dir",),
            "",
            ML_BASE_PARAMS,
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"dir"},
            _load_deeptb,
            ml_style=True,
        ),
        ModelConnector(
            "hamgnn",
            ("file",),
            "",
            ML_BASE_PARAMS,
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"file"},
            _load_hamgnn,
            ml_style=True,
        ),
        ModelConnector(
            "threebodytb",
            ("file",),
            "",
            ML_BASE_PARAMS | {"orbitals_per_element"},
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"file"},
            _load_threebodytb,
            ml_style=True,
        ),
        ModelConnector(
            "abacus",
            ("file",),
            "",
            ML_BASE_PARAMS | {"overlap_file", "stru_file", "orbital_dir",
                               "orbitals_per_element", "position_file"},
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"file"},
            _load_abacus,
            covariant_rule="param",
            covariant_param="position_file",
            path_params=("overlap_file", "stru_file", "orbital_dir", "position_file"),
            required_params=("overlap_file", "stru_file"),
            required_param_hint="(the S(R) file and the STRU of the same structure)",
            ml_style=True,
        ),
        ModelConnector(
            "gpaw",
            ("file",),
            "",
            ML_BASE_PARAMS | {"soc", "quadrature_order", "dipole_cache"},
            ("bare", "berry_connection"),
            _BASE_MODEL_KEYS | {"file"},
            _load_gpaw,
            covariant_rule="always",
            ml_style=True,
        ),
    )
}

#: ``model`` keys of the kinds without a connector entry, besides ``kind``
#: and ``params`` (see :data:`_BASE_MODEL_KEYS` for ``operators``).
NON_CONNECTOR_MODEL_KEYS: dict[str, frozenset[str]] = {
    "csr": frozenset({"dir", "label", "hr_file", "positions",
                      "orbitals_per_atom", "atom_map", "supercell", "operators"}),
    "generator": frozenset({"name", "tutorial", "operators"}),
}


def effective_velocity_modes(connector: ModelConnector,
                             params: Mapping[str, Any]) -> tuple[str, ...]:
    """Return the velocity modes accepted by ``connector`` given ``params``.

    ``covariant`` joins the base ``velocity_modes`` when the entry's rule says
    ``always`` (SIESTA, OpenMX, PySCF, Wannier90, GPAW) or ``param`` and the
    named params key holds a value (ABACUS with ``position_file``).
    """
    modes = list(connector.velocity_modes)
    if connector.covariant_rule == "always" or (
            connector.covariant_rule == "param" and params.get(connector.covariant_param)):
        modes.append("covariant")
    return tuple(modes)


def velocity_mode_error(connector: ModelConnector,
                        params: Mapping[str, Any]) -> str:
    """The one velocity-mode error wording, with why ``covariant`` is missing."""
    modes = list(effective_velocity_modes(connector, params))
    message = f"model.params.velocity_mode for {connector.kind} must be one of {modes}"
    if "covariant" not in modes:
        if connector.covariant_rule == "param":
            message += (f"; 'covariant' needs model.params.{connector.covariant_param} "
                        "(ABACUS's data-rR-sparse.csr)")
        else:
            message += "; no position matrix, so no 'covariant'"
    return message


def metric_trigger_clause() -> str:
    """The run-file wording of the Wannier90 trigger params, from the registry."""
    connector = next(c for c in CONNECTORS.values() if c.metric_trigger is not None)
    names = [f"params.{key}" for key in connector.metric_trigger]
    return ", ".join(names[:-1]) + " or " + names[-1]


def _params_object(kind: str, raw_params: Any) -> Mapping[str, Any]:
    if raw_params is None:
        return {}
    if not isinstance(raw_params, Mapping):
        raise ValueError(f"model.params ({kind}) must be a JSON object")
    return raw_params


def _refuse_unknown_param_keys(kind: str, allowed_params: frozenset[str],
                               params: Mapping[str, Any]) -> None:
    for key in params.keys():
        if key not in allowed_params:
            raise ValueError(
                f"Unknown key '{key}' in run file (model.params ({kind})); "
                f"allowed keys are {sorted(allowed_params)}"
            )


def _refuse_trigger_only_params(connector: ModelConnector,
                                params: Mapping[str, Any]) -> None:
    for key in _TRIGGER_ONLY_PARAMS:
        if key in params:
            raise ValueError(
                f"model.params.{key} for kind='{connector.kind}' needs the metric executor; "
                "set params.velocity_mode, params.spin or params.fermi_energy to use it"
            )


def _validate_numeric_params(params: Mapping[str, Any]) -> None:
    """Validate the numeric and boolean params of the one rule set.

    ``fermi_energy`` admits null (the default); every other value must be a
    real JSON number of the right range, so a value of the wrong type is a
    ``ValueError`` naming the key, never a bare ``float()`` ``TypeError``.
    """
    if "fermi_energy" in params:
        value = params["fermi_energy"]
        if value is not None and (
                isinstance(value, bool) or not isinstance(value, (int, float))
                or not np.isfinite(value)):
            raise ValueError("model.params.fermi_energy must be a finite number or null (eV)")
    for key, message in (
            ("oracle_degeneracy_tolerance",
             "oracle_degeneracy_tolerance must be positive and finite (eV)"),
            ("bloch_memory_mb", "bloch_memory_mb must be positive and finite (MiB)"),
            ("quadrature_tolerance", "model.params.quadrature_tolerance must be positive and finite")):
        if key in params:
            value = params[key]
            if (isinstance(value, bool) or not isinstance(value, (int, float))
                    or not np.isfinite(value) or value <= 0):
                raise ValueError(message)
    if "oracle" in params and not isinstance(params["oracle"], bool):
        raise ValueError("model.params.oracle must be true or false")
    for key in ("quadrature_order", "quadrature_max_order"):
        if key in params:
            value = params[key]
            if isinstance(value, bool) or not isinstance(value, int) or value < 2:
                raise ValueError(f"model.params.{key} must be an integer >= 2")


def validate_connector_params(connector: ModelConnector,
                              model: Mapping[str, Any],
                              *, for_metric_executor: bool = True) -> None:
    """Validate the ``params`` of a connector model with the one rule set.

    Every connector: the params must be a JSON object or null; unknown keys
    are refused; ``velocity_mode`` uses the one wording; the numeric params
    are validated for every kind that accepts them. A trigger kind validated
    without its trigger param (a plain Wannier90 model) may not use the
    executor-only params and gets no further checks. The path keys of
    ``model`` itself are checked by the caller
    (``runfile.validate_runfile_dict``), as is the per-kind ``model`` key set.
    """
    params = _params_object(connector.kind, model.get("params"))
    _refuse_unknown_param_keys(connector.kind, connector.allowed_params, params)
    if connector.metric_trigger is not None and not for_metric_executor:
        _refuse_trigger_only_params(connector, params)
        return
    if connector.refusal is not None and "spin" in params:
        raise ValueError(connector.refusal)
    if params.get("velocity_mode", "berry_connection") not in effective_velocity_modes(connector, params):
        raise ValueError(velocity_mode_error(connector, params))
    _validate_numeric_params(params)
    if "spin" in params and params["spin"] not in connector.spin_values:
        raise ValueError("model.params.spin must be 'labels' or 'exact'")
    if "nonorthogonal_disorder" in params and not isinstance(
            params["nonorthogonal_disorder"], bool):
        raise ValueError("model.params.nonorthogonal_disorder must be a boolean")
    if "orbitals_per_element" in params:
        per_element = params["orbitals_per_element"]
        if per_element is not None and (not isinstance(per_element, Mapping) or not all(
                isinstance(v, int) and not isinstance(v, bool) and v > 0
                for v in per_element.values())):
            raise ValueError(
                "model.params.orbitals_per_element must map element symbols to positive integers"
            )
    for name in connector.required_params:
        if not params.get(name):
            message = f"model.params.{name} is required for model.kind='{connector.kind}'"
            if connector.required_param_hint:
                message += f" {connector.required_param_hint}"
            raise ValueError(message)
    if "soc" in params and not isinstance(params["soc"], bool):
        raise ValueError("model.params.soc must be true or false")


def validate_model_keys(kind: str, model_keys: frozenset[str],
                        model: Mapping[str, Any]) -> None:
    """Refuse ``model`` keys the kind does not read (audit item B8).

    A null value means "absent" in this schema, and ``run.resolved.json``
    round trips carry every known key set to null, so only a non-null value
    of another kind's key is refused.
    """
    for key in model.keys():
        if key not in model_keys and model[key] is not None:
            raise ValueError(
                f"Unknown key '{key}' in run file (model.kind='{kind}'); "
                f"allowed keys are {sorted(model_keys)}"
            )


def validate_bare_params(kind: str, model: Mapping[str, Any]) -> None:
    """Validate the ``params`` of a kind that uses none (today only ``csr``)."""
    params = _params_object(kind, model.get("params"))
    _refuse_unknown_param_keys(kind, frozenset(), params)


def load_runfile_model(model, velocity_mode, *, phase=None):
    """Load the model of a run file for the metric executors (SIESTA, OpenMX, Wannier90,
    PySCF, PAOFLOW, and the machine-learning Hamiltonians of ``runfile.ML_MODEL_KINDS``).

    ``model`` is the run-file ``model`` mapping. SIESTA keeps its own reader;
    a Wannier90 model reads ``_r.dat`` for the covariant velocity and the
    spin set by ``params.spin``; a PAOFLOW model has no position matrix, so
    ``covariant`` is refused. A machine-learning model is read by its module
    (``models/<kind>.py``) and has no position matrix either.
    """
    kind = model["kind"]
    connector = CONNECTORS.get(kind)
    if connector is None or connector.loader is None:
        raise ValueError(f"model.kind {kind!r} has no metric executor")
    return connector.loader(model, velocity_mode, phase)


def _optional_float(value):
    return None if value is None else float(value)
