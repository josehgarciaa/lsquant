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
  ``dir``+``seed``, ``dir``, or ``file``), and ``missing_path_error``: the
  exact text raised when they do not;
- ``allowed_params``: the accepted keys of ``model.params``;
- ``velocity_modes`` and ``covariant_rule``: the accepted velocity modes,
  where ``always``/``never``/``param`` says when ``covariant`` joins the base
  modes (ABACUS admits it only with its ``position_file`` param);
- ``path_params``: params keys that hold paths resolved like ``model.file``
  (today the four ABACUS keys);
- ``refusal``: text refusing ``spin``/``fermi_energy`` when the reader stores
  them itself (OpenMX, PySCF, PAOFLOW);
- ``metric_trigger``: params keys that route the kind to the metric executor
  (only Wannier90: ``velocity_mode``/``spin``/``fermi_energy``; every other
  connector always runs there);
- ``keep_metric_policy``: whether the executor keeps ``metric.policy`` as
  given (OpenMX, PySCF) or is switched to the exact Lowdin route;
- ``loader``: a callable ``(model, velocity_mode, phase) -> model`` that
  imports its reader module lazily, so optional dependencies (sisl, h5py,
  ...) never load at package import.

Checks. Every entry names its params checks, in order, in ``checks``; the
order is part of the behaviour because a run file that violates two checks at
once must raise the same error it raises today (SIESTA validates its velocity
mode last, after the numeric params; the H(R)/S(R) kinds validate it first).

Pitfall: this module is a behaviour-preserving refactor. The four orthogonal
kinds (wannier, openmx, pyscf, paoflow) do not validate ``fermi_energy``
finiteness or ``bloch_memory_mb`` positivity, a Wannier90 model without a
trigger param keeps its params unvalidated, and ``fermi_energy: null`` raises
``TypeError`` in some loaders but means "the file value" in others. Those
discrepancies are catalogued bugs, kept byte-identical here on purpose.

To-do: unify the three validation flavours (orthogonal / SIESTA / H(R)S(R))
into one validator with a single wording per message, and validate ``params``
of every kind including plain Wannier90 and ``csr``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

import numpy as np

#: Base ``model.params`` keys accepted by the four orthogonal models
#: (wannier, openmx, pyscf, paoflow).
_METRIC_BASE_PARAMS = frozenset(
    {"velocity_mode", "spin", "fermi_energy", "oracle",
     "oracle_degeneracy_tolerance", "bloch_memory_mb"}
)

#: Base ``model.params`` keys accepted by the H(R)/S(R) Hamiltonians
#: (deeph, deeptb, hamgnn, threebodytb, abacus, gpaw).
_ML_BASE_PARAMS = frozenset(
    {"velocity_mode", "fermi_energy", "oracle",
     "oracle_degeneracy_tolerance", "bloch_memory_mb"}
)


@dataclass(frozen=True)
class ModelConnector:
    """One ``model.kind`` of the metric executors, as data (see module docstring)."""

    kind: str
    path_keys: tuple[str, ...]
    missing_path_error: str
    allowed_params: frozenset[str]
    velocity_modes: tuple[str, ...]
    velocity_error: str
    loader: Callable[[Mapping[str, Any], str, Any], Any] | None
    covariant_rule: str = "never"           # "always" | "never" | "param"
    covariant_param: str = ""               # params key that admits covariant
    path_params: tuple[str, ...] = ()       # params keys resolved like model.file
    refusal: str | None = None              # text refusing spin/fermi_energy
    metric_trigger: tuple[str, ...] | None = None   # None: always metric executor
    keep_metric_policy: bool = False        # True: keep metric.policy as given
    ml_style: bool = False                  # True: H(R)/S(R) input like ML kinds
    checks: tuple[str, ...] = ("velocity", "tolerance")
    spin_values: tuple[str, ...] = ()
    params_object_error: str = ""
    # True: refuse any non-Mapping params, even a falsy one (the original
    # SIESTA block checked ``params is not None``; the other validators
    # applied ``params or {}`` first and so accept e.g. an empty list).
    params_object_strict: bool = False

    def __post_init__(self):
        if not self.params_object_error:
            object.__setattr__(
                self, "params_object_error",
                f"model.params ({self.kind}) must be a JSON object",
            )


def _load_siesta(model, velocity_mode, phase):
    from .siesta import load_siesta_model
    params = model.get("params") or {}
    return load_siesta_model(
        model["dir"], model["label"], dipoles=velocity_mode == "covariant",
        quadrature_order=int(params.get("quadrature_order", 24)),
        quadrature_tolerance=float(params.get("quadrature_tolerance", 1e-5)),
        quadrature_max_order=int(params.get("quadrature_max_order", 384)),
        phase=phase)


def _load_wannier(model, velocity_mode, phase):
    from .orthogonal import load_wannier_metric_model
    params = model.get("params") or {}
    return load_wannier_metric_model(
        model["dir"], model["seed"], spin=params.get("spin", "labels"),
        positions=True if velocity_mode == "covariant" else False,
        fermi_energy=float(params.get("fermi_energy", 0.0)))


def _load_openmx(model, velocity_mode, phase):
    from .openmx import load_openmx_model
    return load_openmx_model(model["file"], dipoles=velocity_mode == "covariant")


def _load_pyscf(model, velocity_mode, phase):
    from .pyscf_lcao import load_pyscf_model
    return load_pyscf_model(model["file"], connection=velocity_mode == "covariant")


def _load_paoflow(model, velocity_mode, phase):
    if velocity_mode == "covariant":
        raise ValueError("a PAOFLOW model has no position matrix; use bare or berry_connection")
    from .orthogonal import load_paoflow_model
    return load_paoflow_model(model["file"])


def _load_deeph(model, velocity_mode, phase):
    from .deeph import load_deeph_model
    params = model.get("params") or {}
    return load_deeph_model(model["dir"], hamiltonian_file=params.get("hamiltonian_file"),
                            fermi_energy=_optional_float(params.get("fermi_energy")))


def _load_deeptb(model, velocity_mode, phase):
    from .deeptb import load_deeptb_model
    params = model.get("params") or {}
    return load_deeptb_model(model["dir"], fermi_energy=float(params.get("fermi_energy", 0.0)))


def _load_hamgnn(model, velocity_mode, phase):
    from .hamgnn import load_hamgnn_model
    params = model.get("params") or {}
    return load_hamgnn_model(model["file"], fermi_energy=float(params.get("fermi_energy", 0.0)))


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
        fermi_energy=float(params.get("fermi_energy", 0.0)),
        position_file=params.get("position_file") if velocity_mode == "covariant" else None)


def _load_gpaw(model, velocity_mode, phase):
    from .gpaw import load_gpaw_model
    params = model.get("params") or {}
    return load_gpaw_model(
        model["file"], soc=bool(params.get("soc", True)),
        fermi_energy=_optional_float(params.get("fermi_energy")),
        dipoles=velocity_mode == "covariant",
        quadrature_order=int(params.get("quadrature_order", 48)),
        dipole_cache=params.get("dipole_cache"))


#: ``{kind: connector}`` in the order in which run-file messages list the
#: kinds; ``ML_MODEL_KINDS`` and ``known_kinds`` are slices of this order.
CONNECTORS: dict[str, ModelConnector] = {
    c.kind: c for c in (
        ModelConnector(
            "siesta",
            ("dir", "label"),
            "SIESTA model requires dir and label",
            frozenset({"velocity_mode", "oracle", "oracle_degeneracy_tolerance",
                       "bloch_memory_mb", "quadrature_order", "quadrature_max_order",
                       "quadrature_tolerance", "nonorthogonal_disorder"}),
            ("bare", "berry_connection"),
            "unknown SIESTA velocity_mode",
            _load_siesta,
            covariant_rule="always",
            # SIESTA validates its velocity mode after the numeric params, as
            # the original inline block did.
            checks=("nonorthogonal_disorder", "tolerance", "bloch_memory", "velocity"),
            params_object_error="SIESTA model.params must be a JSON object",
            params_object_strict=True,
        ),
        ModelConnector(
            "openmx",
            ("file",),
            "Run file with model.kind='openmx' requires 'file' (the .scfout file)",
            _METRIC_BASE_PARAMS,
            ("bare", "berry_connection"),
            "model.params.velocity_mode for {kind} must be one of {modes}",
            _load_openmx,
            covariant_rule="always",
            refusal="an OpenMX model carries its spinors and stores energies from the "
                    "chemical potential; remove model.params.spin and model.params.fermi_energy",
            checks=("spin_fermi_refusal", "velocity", "tolerance"),
            keep_metric_policy=True,
        ),
        ModelConnector(
            "pyscf",
            ("file",),
            "Run file with model.kind='pyscf' requires 'file' (see write_pyscf_model)",
            _METRIC_BASE_PARAMS,
            ("bare", "berry_connection"),
            "model.params.velocity_mode for {kind} must be one of {modes}",
            _load_pyscf,
            covariant_rule="always",
            refusal="a PySCF model carries its spinors and stores energies from its Fermi "
                    "energy; remove model.params.spin and model.params.fermi_energy",
            checks=("spin_fermi_refusal", "velocity", "tolerance"),
            keep_metric_policy=True,
        ),
        ModelConnector(
            "wannier",
            ("dir", "seed"),
            "Run file with model.kind='wannier' requires 'dir' and 'seed'",
            _METRIC_BASE_PARAMS,
            ("bare", "berry_connection"),
            "model.params.velocity_mode for {kind} must be one of {modes}",
            _load_wannier,
            covariant_rule="always",
            metric_trigger=("velocity_mode", "spin", "fermi_energy"),
            checks=("velocity", "tolerance", "spin_values"),
            spin_values=("labels", "exact"),
        ),
        ModelConnector(
            "paoflow",
            ("file",),
            "Run file with model.kind='paoflow' requires 'file' (see write_paoflow_model)",
            _METRIC_BASE_PARAMS,
            ("bare", "berry_connection"),
            "model.params.velocity_mode for {kind} must be one of {modes}",
            _load_paoflow,
            refusal="a PAOFLOW model carries its own on-site spin and energies from E_F; "
                    "remove model.params.spin and model.params.fermi_energy",
            checks=("spin_fermi_refusal", "velocity", "tolerance"),
        ),
        ModelConnector(
            "deeph",
            ("dir",),
            "Run file with model.kind='deeph' requires 'dir'",
            _ML_BASE_PARAMS | {"hamiltonian_file"},
            ("bare", "berry_connection"),
            "model.params.velocity_mode for {kind} must be 'bare' or 'berry_connection' "
            "(no position matrix, so no 'covariant')",
            _load_deeph,
            checks=("velocity", "fermi_finite", "tolerance", "bloch_memory"),
            ml_style=True,
        ),
        ModelConnector(
            "deeptb",
            ("dir",),
            "Run file with model.kind='deeptb' requires 'dir'",
            _ML_BASE_PARAMS,
            ("bare", "berry_connection"),
            "model.params.velocity_mode for {kind} must be 'bare' or 'berry_connection' "
            "(no position matrix, so no 'covariant')",
            _load_deeptb,
            checks=("velocity", "fermi_finite", "tolerance", "bloch_memory"),
            ml_style=True,
        ),
        ModelConnector(
            "hamgnn",
            ("file",),
            "Run file with model.kind='hamgnn' requires 'file'",
            _ML_BASE_PARAMS,
            ("bare", "berry_connection"),
            "model.params.velocity_mode for {kind} must be 'bare' or 'berry_connection' "
            "(no position matrix, so no 'covariant')",
            _load_hamgnn,
            checks=("velocity", "fermi_finite", "tolerance", "bloch_memory"),
            ml_style=True,
        ),
        ModelConnector(
            "threebodytb",
            ("file",),
            "Run file with model.kind='threebodytb' requires 'file'",
            _ML_BASE_PARAMS | {"orbitals_per_element"},
            ("bare", "berry_connection"),
            "model.params.velocity_mode for {kind} must be 'bare' or 'berry_connection' "
            "(no position matrix, so no 'covariant')",
            _load_threebodytb,
            checks=("velocity", "fermi_finite", "tolerance", "bloch_memory",
                    "orbitals_per_element"),
            ml_style=True,
        ),
        ModelConnector(
            "abacus",
            ("file",),
            "Run file with model.kind='abacus' requires 'file'",
            _ML_BASE_PARAMS | {"overlap_file", "stru_file", "orbital_dir",
                               "orbitals_per_element", "position_file"},
            ("bare", "berry_connection"),
            "model.params.velocity_mode for abacus must be 'bare' or 'berry_connection'; "
            "'covariant' needs model.params.position_file (ABACUS's data-rR-sparse.csr)",
            _load_abacus,
            covariant_rule="param",
            covariant_param="position_file",
            path_params=("overlap_file", "stru_file", "orbital_dir", "position_file"),
            checks=("velocity", "fermi_finite", "tolerance", "bloch_memory",
                    "orbitals_per_element", "abacus_files"),
            ml_style=True,
        ),
        ModelConnector(
            "gpaw",
            ("file",),
            "Run file with model.kind='gpaw' requires 'file'",
            _ML_BASE_PARAMS | {"soc", "quadrature_order", "dipole_cache"},
            ("bare", "berry_connection"),
            "model.params.velocity_mode for {kind} must be 'bare' or 'berry_connection' "
            "(no position matrix, so no 'covariant')",
            _load_gpaw,
            covariant_rule="always",
            checks=("velocity", "fermi_finite", "tolerance", "bloch_memory", "gpaw_soc"),
            ml_style=True,
        ),
    )
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


def _refuse_unknown_param_keys(connector: ModelConnector,
                               params: Mapping[str, Any]) -> None:
    for key in params.keys():
        if key not in connector.allowed_params:
            raise ValueError(
                f"Unknown key '{key}' in run file (model.params ({connector.kind})); "
                f"allowed keys are {sorted(connector.allowed_params)}"
            )


def _check_spin_fermi_refusal(connector: ModelConnector, params: Mapping[str, Any]) -> None:
    if connector.refusal is not None and ("spin" in params or "fermi_energy" in params):
        raise ValueError(connector.refusal)


def _check_velocity(connector: ModelConnector, params: Mapping[str, Any]) -> None:
    modes = effective_velocity_modes(connector, params)
    if params.get("velocity_mode", "berry_connection") not in modes:
        raise ValueError(connector.velocity_error.format(
            kind=connector.kind, modes=list(modes)))


def _check_fermi_finite(connector: ModelConnector, params: Mapping[str, Any]) -> None:
    if "fermi_energy" in params and not np.isfinite(float(params["fermi_energy"])):
        raise ValueError("model.params.fermi_energy must be a finite number (eV)")


def _check_tolerance(connector: ModelConnector, params: Mapping[str, Any]) -> None:
    tolerance = float(params.get("oracle_degeneracy_tolerance", 1e-6))
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("oracle_degeneracy_tolerance must be positive and finite (eV)")


def _check_bloch_memory(connector: ModelConnector, params: Mapping[str, Any]) -> None:
    memory = float(params.get("bloch_memory_mb", 256.0))
    if not np.isfinite(memory) or memory <= 0:
        raise ValueError("bloch_memory_mb must be positive and finite (MiB)")


def _check_spin_values(connector: ModelConnector, params: Mapping[str, Any]) -> None:
    if "spin" in params and params["spin"] not in connector.spin_values:
        raise ValueError("model.params.spin must be 'labels' or 'exact'")


def _check_nonorthogonal_disorder(connector: ModelConnector,
                                  params: Mapping[str, Any]) -> None:
    if "nonorthogonal_disorder" in params and not isinstance(
            params["nonorthogonal_disorder"], bool):
        raise ValueError("model.params.nonorthogonal_disorder must be a boolean")


def _check_orbitals_per_element(connector: ModelConnector,
                                params: Mapping[str, Any]) -> None:
    per_element = params.get("orbitals_per_element")
    if per_element is not None and (not isinstance(per_element, Mapping) or not all(
            isinstance(v, int) and not isinstance(v, bool) and v > 0
            for v in per_element.values())):
        raise ValueError(
            "model.params.orbitals_per_element must map element symbols to positive integers"
        )


def _check_abacus_files(connector: ModelConnector, params: Mapping[str, Any]) -> None:
    for name in ("overlap_file", "stru_file"):
        if not params.get(name):
            raise ValueError(f"model.params.{name} is required for model.kind='abacus' "
                             "(the S(R) file and the STRU of the same structure)")


def _check_gpaw_soc(connector: ModelConnector, params: Mapping[str, Any]) -> None:
    if not isinstance(params.get("soc", True), bool):
        raise ValueError("model.params.soc must be true or false")


_PARAM_CHECKS = {
    "spin_fermi_refusal": _check_spin_fermi_refusal,
    "velocity": _check_velocity,
    "fermi_finite": _check_fermi_finite,
    "tolerance": _check_tolerance,
    "bloch_memory": _check_bloch_memory,
    "spin_values": _check_spin_values,
    "nonorthogonal_disorder": _check_nonorthogonal_disorder,
    "orbitals_per_element": _check_orbitals_per_element,
    "abacus_files": _check_abacus_files,
    "gpaw_soc": _check_gpaw_soc,
}


def validate_connector_params(connector: ModelConnector,
                              model: Mapping[str, Any]) -> None:
    """Validate the ``params`` of a connector model against its registry entry.

    Runs the entry's checks in their declared order, so a run file that
    violates two checks raises the same error it raised under the original
    hand-written per-kind chain. The path keys of ``model`` itself are
    checked by the caller (``runfile.validate_runfile_dict``), which also
    decides whether a Wannier90 model without a trigger param is validated
    at all (it is not, matching the original behaviour).
    """
    raw_params = model.get("params")
    if raw_params is None:
        params = {}
    elif isinstance(raw_params, Mapping):
        params = raw_params
    elif connector.params_object_strict or raw_params:
        raise ValueError(connector.params_object_error)
    else:
        # Falsy non-Mapping (e.g. []): the original validators saw `or {}`.
        params = {}
    _refuse_unknown_param_keys(connector, params)
    for name in connector.checks:
        _PARAM_CHECKS[name](connector, params)


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
