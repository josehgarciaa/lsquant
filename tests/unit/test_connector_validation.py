"""One validator and one Fermi-energy policy for every connector kind.

Part A pins the rule set of ``lsquant.models.registry.validate_connector_params``
for every registered kind (parametrised over the registry, so a new connector
is covered automatically) and for ``csr``: the params-object rule with its one
message, the unknown-key rule, the numeric rules shared by every kind that
accepts them (a wrong type is a ``ValueError`` naming the key), the one
velocity-mode wording, the per-kind ``model`` key set (audit B8) and the
missing-path texts.

Part B pins the Fermi-energy policy of the readers on tiny fixture files,
built here the way the reader tests build them: the default is unchanged,
``fermi_energy: null`` is the default, passing the file's own value gives a
bit-identical model, and an explicit value shifts the generalised eigenvalues
of ``H(k), S(k)`` at one k point by exactly ``-(x - x_file)`` eV. A reader
whose fixture needs an optional package that is absent is skipped, as the
reader tests skip them. Run files are not needed: the loaders are called
directly, as the metric executors call them.
"""
import gzip
import json
import os
import sys

import numpy as np
import pytest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if os.path.join(REPO, "python") not in sys.path:
    sys.path.insert(0, os.path.join(REPO, "python"))

pytest.importorskip("lsquant", reason="lsquant._core not built")

from lsquant.models.registry import (  # noqa: E402
    CONNECTORS,
    NON_CONNECTOR_MODEL_KEYS,
    effective_velocity_modes,
    metric_trigger_clause,
    velocity_mode_error,
)
from lsquant.runfile import RunFile, execute_runfile  # noqa: E402

BASE = dict(word="Tr[delta(E-H)]", expansion={"moments": 8},
            basis={"name": "eigen", "supercell": [2, 2, 1]})
K = np.array([0.13, 0.27, 0.0])
X_FILE = -0.35
X_NEW = 0.85
CONNECTOR_LIST = list(CONNECTORS.values())
KINDS = [c.kind for c in CONNECTOR_LIST]


def _paths(connector):
    return {key: "x" for key in connector.path_keys}


def _accepted_params(connector):
    """Valid params of the kind (with the ABACUS required files)."""
    params = ({"velocity_mode": "bare", "oracle": True}
              if connector.metric_trigger is not None else {"oracle": True})
    if connector.kind == "abacus":
        params = dict(params, overlap_file="s", stru_file="STRU")
    return params


def _model(connector, params="__absent__", **extra):
    model = {"kind": connector.kind, **_paths(connector)}
    if params != "__absent__":
        model["params"] = params
    model.update(extra)
    return model


def _validated(connector, **kwargs):
    return RunFile.from_dict(dict(BASE, model=_model(connector, **kwargs)))


# --------------------------------------------------------------------- A.1

@pytest.mark.parametrize("connector", CONNECTOR_LIST, ids=KINDS)
@pytest.mark.parametrize("bad", [[], "x", 5, True])
def test_params_must_be_json_object(connector, bad):
    _validated(connector)                      # absent
    _validated(connector, params=None)         # null
    _validated(connector, params={})           # empty
    with pytest.raises(ValueError) as exc:
        _validated(connector, params=bad)
    assert str(exc.value) == f"model.params ({connector.kind}) must be a JSON object"


def test_csr_params_object_rule():
    model = {"kind": "csr", "dir": "d", "label": "l", "params": {"unused": 1}}
    with pytest.raises(ValueError) as exc:
        RunFile.from_dict(dict(BASE, model=model))
    assert str(exc.value) == (
        "Unknown key 'unused' in run file (model.params (csr)); allowed keys are []"
    )
    for bad in ([], "x"):
        with pytest.raises(ValueError) as exc:
            RunFile.from_dict(dict(BASE, model={"kind": "csr", "dir": "d", "label": "l",
                                                "params": bad}))
        assert str(exc.value) == "model.params (csr) must be a JSON object"


# --------------------------------------------------------------------- A.2

@pytest.mark.parametrize("connector", CONNECTOR_LIST, ids=KINDS)
def test_unknown_params_key_refused(connector):
    _validated(connector, params=_accepted_params(connector))
    params = dict(_accepted_params(connector), zz_unknown=1)
    with pytest.raises(ValueError) as exc:
        _validated(connector, params=params)
    assert str(exc.value) == (
        f"Unknown key 'zz_unknown' in run file (model.params ({connector.kind})); "
        f"allowed keys are {sorted(connector.allowed_params)}"
    )


# --------------------------------------------------------------------- A.3

@pytest.mark.parametrize("connector", CONNECTOR_LIST, ids=KINDS)
@pytest.mark.parametrize("bad", ["abc", True, [1], float("nan"), float("inf")])
def test_fermi_energy_finite_number_or_null(connector, bad):
    params = dict(_accepted_params(connector), fermi_energy=None)
    _validated(connector, params=params)
    params = dict(_accepted_params(connector), fermi_energy=0.4)
    _validated(connector, params=params)
    params = dict(_accepted_params(connector), fermi_energy=bad)
    with pytest.raises(ValueError) as exc:
        _validated(connector, params=params)
    assert str(exc.value) == "model.params.fermi_energy must be a finite number or null (eV)"


_POSITIVE_FINITE = (
    ("oracle_degeneracy_tolerance",
     "oracle_degeneracy_tolerance must be positive and finite (eV)"),
    ("bloch_memory_mb", "bloch_memory_mb must be positive and finite (MiB)"),
    ("quadrature_tolerance", "model.params.quadrature_tolerance must be positive and finite"),
)


def _numeric_cases():
    for connector in CONNECTOR_LIST:
        for key, message in _POSITIVE_FINITE:
            if key in connector.allowed_params:
                for bad in (0, -1.0, "x", True):
                    yield connector, key, bad, message
        for key in ("quadrature_order", "quadrature_max_order"):
            if key in connector.allowed_params:
                for bad in (1, 2.5, "24", True):
                    yield connector, key, bad, f"model.params.{key} must be an integer >= 2"
        if "oracle" in connector.allowed_params:
            yield connector, "oracle", "yes", "model.params.oracle must be true or false"
            yield connector, "oracle", 1, "model.params.oracle must be true or false"


@pytest.mark.parametrize(
    "connector,key,bad,message",
    list(_numeric_cases()),
    ids=[f"{c.kind}-{k}-{b!r}" for c, k, b, _ in _numeric_cases()],
)
def test_numeric_params_wrong_type_names_the_key(connector, key, bad, message):
    params = dict(_accepted_params(connector), **{key: bad})
    with pytest.raises(ValueError) as exc:
        _validated(connector, params=params)
    assert str(exc.value) == message


# --------------------------------------------------------------------- A.4

@pytest.mark.parametrize("connector", CONNECTOR_LIST, ids=KINDS)
def test_velocity_mode_one_wording(connector):
    _validated(connector)                                   # default
    _validated(connector, params=dict(_accepted_params(connector),
                                      velocity_mode="bare"))
    params = dict(_accepted_params(connector), velocity_mode="bogus")
    with pytest.raises(ValueError) as exc:
        _validated(connector, params=params)
    assert str(exc.value) == velocity_mode_error(connector, params)
    assert str(exc.value).startswith(
        f"model.params.velocity_mode for {connector.kind} must be one of ")
    modes = effective_velocity_modes(connector, params)
    if "covariant" in modes:
        assert ";" not in str(exc.value)
    elif connector.covariant_rule == "param":
        assert str(exc.value).endswith(
            f"; 'covariant' needs model.params.{connector.covariant_param} "
            "(ABACUS's data-rR-sparse.csr)")
    else:
        assert str(exc.value).endswith("; no position matrix, so no 'covariant'")
    accepted = _accepted_params(connector)
    if "covariant" in effective_velocity_modes(connector, accepted):
        params = dict(accepted, velocity_mode="covariant")
        if connector.covariant_rule == "param":
            params[connector.covariant_param] = "p"
        _validated(connector, params=params)


@pytest.mark.parametrize("connector", CONNECTOR_LIST, ids=KINDS)
def test_missing_path_error_text(connector):
    with pytest.raises(ValueError) as exc:
        RunFile.from_dict(dict(BASE, model={"kind": connector.kind}))
    assert str(exc.value) == connector.missing_path_error
    assert str(exc.value).startswith(f"model.kind='{connector.kind}' requires ")


# --------------------------------------------------------------------- A.5

def _foreign_key(connector):
    return "seed" if "seed" not in connector.model_keys else "label"


@pytest.mark.parametrize("connector", CONNECTOR_LIST, ids=KINDS)
def test_model_keys_per_kind(connector):
    foreign = _foreign_key(connector)
    with pytest.raises(ValueError) as exc:
        _validated(connector, **{foreign: "x"})
    assert str(exc.value) == (
        f"Unknown key '{foreign}' in run file (model.kind='{connector.kind}'); "
        f"allowed keys are {sorted(connector.model_keys)}"
    )
    # A null value means absent (run.resolved.json round trips carry nulls).
    _validated(connector, **{foreign: None})
    _validated(connector, operators=["VX"])   # read by the word compiler for any kind


def test_model_keys_csr_and_generator():
    with pytest.raises(ValueError) as exc:
        RunFile.from_dict(dict(BASE, model={"kind": "csr", "dir": "d", "label": "l",
                                            "seed": "s"}))
    assert "model.kind='csr'" in str(exc.value) and "'seed'" in str(exc.value)
    allowed = NON_CONNECTOR_MODEL_KEYS["csr"] | {"kind", "params"}
    with pytest.raises(ValueError) as exc:
        RunFile.from_dict(dict(BASE, model={"kind": "csr", "dir": "d", "label": "l",
                                            "file": "x"}))
    assert str(exc.value) == (
        f"Unknown key 'file' in run file (model.kind='csr'); allowed keys are {sorted(allowed)}"
    )
    csr = {"kind": "csr", "dir": "d", "label": "l", "hr_file": "h",
           "positions": "p", "orbitals_per_atom": 2, "atom_map": "a",
           "supercell": [2, 2], "operators": None}
    RunFile.from_dict(dict(BASE, model=csr))
    with pytest.raises(ValueError) as exc:
        RunFile.from_dict(dict(BASE, model={"kind": "generator", "name": "make_chain",
                                            "seed": 1}))
    assert "model.kind='generator'" in str(exc.value)
    gen = {"kind": "generator", "name": "make_chain", "tutorial": 1, "params": {}}
    RunFile.from_dict(dict(BASE, model=gen))


# --------------------------------------------------------------------- A.6

def test_metric_executor_messages_name_all_trigger_params():
    assert metric_trigger_clause() == "params.velocity_mode, params.spin or params.fermi_energy"
    with pytest.raises(ValueError) as exc:
        RunFile.from_dict(dict(BASE, model={"kind": "csr", "dir": "d", "label": "l"},
                               reference=True))
    assert "or wannier with params.velocity_mode, params.spin or params.fermi_energy" in str(exc.value)


def test_dispatch_message_without_wrong_reason():
    raw = dict(BASE, basis={"name": "real", "supercell": [2, 2, 1]},
               model={"kind": "paoflow", "file": "x.npz"})
    with pytest.raises(ValueError) as exc:
        execute_runfile(raw)
    assert str(exc.value) == (
        "model.kind='paoflow' runs in the clean bloch or eigen basis "
        "(real basis and disorder are available for siesta)"
    )


# ================================================================= B

def _generalized_eigs(model, k):
    import scipy.linalg as sla
    H = model.Hk(np.atleast_2d(k))[0]
    S = model.Sk(np.atleast_2d(k))[0]
    return sla.eigh(H, S)[0]


def _assert_fermi_policy(load, *, x_file):
    """Default unchanged, null is the default, the file value (0.0 when the
    file carries none) gives a bit-identical model, and an explicit value
    shifts the eigenvalues by exactly -(x - x_file)."""
    has_file_reference = x_file is not None
    expected_origin = x_file if has_file_reference else 0.0
    m_default = load(None)
    assert m_default.metadata["fermi_energy_source"] == ("file" if has_file_reference else "none")
    assert m_default.metadata["fermi_energy_subtracted_eV"] == expected_origin
    m_new = load(X_NEW)
    assert m_new.metadata["fermi_energy_source"] == "run file"
    assert m_new.metadata["fermi_energy_subtracted_eV"] == X_NEW
    delta = -(X_NEW - expected_origin)
    np.testing.assert_allclose(
        _generalized_eigs(m_new, K), _generalized_eigs(m_default, K) + delta, atol=1e-10)
    m_file = load(expected_origin)
    assert np.array_equal(m_default.hoppings, m_file.hoppings)
    assert np.array_equal(m_default.overlaps, m_file.overlaps)
    assert m_file.metadata["fermi_energy_source"] == ("file" if has_file_reference else "none")


def _rng_blocks(rng, cells, n, scale, *, identity=False):
    """M(R) with M(-R) = M(R)^+, an Hermitian R = 0 block."""
    mats = {}
    for R in cells:
        neg = tuple(-x for x in R)
        if neg in mats:
            mats[R] = mats[neg].conj().T
            continue
        m = scale * (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
        if R == neg:
            m = 0.5 * (m + m.conj().T) + (np.eye(n) if identity else 0)
        mats[R] = m
    return mats


def _wannier_files(tmp_path, rng):
    from lsquant.models.orthogonal import load_wannier_metric_model
    cells = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
    mats = _rng_blocks(rng, cells, 3, 0.4)
    np.savetxt(tmp_path / "w.uc", np.diag([3.0, 3.2, 20.0]))
    (tmp_path / "w.xyz").write_text(
        "3\ncomment\nWF1 0.0 0.0 0.0\nWF2 1.0 0.2 0.0\nWF3 0.3 1.1 0.0\n")
    with open(tmp_path / "w_hr.dat", "w", encoding="utf-8") as fh:
        fh.write("written by test_connector_validation\n3\n5\n" + " ".join(["1"] * 5) + "\n")
        for R in cells:
            for i in range(3):
                for j in range(3):
                    fh.write(f"{R[0]} {R[1]} {R[2]} {i + 1} {j + 1} "
                             f"{mats[R][i, j].real:.16e} {mats[R][i, j].imag:.16e}\n")
    return lambda fermi_energy=None: load_wannier_metric_model(str(tmp_path), "w", fermi_energy=fermi_energy)


def test_fermi_policy_wannier(tmp_path):
    load = _wannier_files(tmp_path, np.random.default_rng(21))
    _assert_fermi_policy(load, x_file=None)


def _paoflow_file(tmp_path, rng):
    from lsquant.models.orthogonal import load_paoflow_model
    cells = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]
    mats = _rng_blocks(rng, cells, 4, 0.4)
    r_int = np.array(cells, int)
    spin = np.zeros((3, 4, 4))
    spin[2] = np.eye(4)
    np.savez(tmp_path / "p.npz", uc=np.diag([3.0, 3.1, 18.0]),
             tau=np.zeros((4, 3)), r_int=r_int,
             hoppings=np.array([mats[R] for R in cells]), spin=spin,
             fermi_energy=X_FILE, mesh=np.array([1, 1, 1]), edge_ratio=0.0)
    return lambda fermi_energy=None: load_paoflow_model(str(tmp_path / "p.npz"), fermi_energy=fermi_energy)


def test_fermi_policy_paoflow(tmp_path):
    _assert_fermi_policy(_paoflow_file(tmp_path, np.random.default_rng(22)), x_file=X_FILE)


def _pyscf_file(tmp_path, rng):
    from lsquant.models.pyscf_lcao import load_pyscf_model
    cells = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]
    nao = 2
    h = _rng_blocks(rng, cells, nao, 0.4)
    s = _rng_blocks(rng, cells, nao, 0.05, identity=True)
    hoppings = np.zeros((len(cells), 2 * nao, 2 * nao), complex)
    overlap_orb = np.array([s[R] for R in cells])
    position_orb = np.zeros((3, len(cells), nao, nao))
    r_int = np.array(cells, int)
    tau = np.zeros((2 * nao, 3))
    # interleaved spinors: H stored shifted by X_FILE, as write_pyscf_model does
    for i, R in enumerate(cells):
        hoppings[i] = np.kron(h[R] - X_FILE * s[R], np.eye(2))
    np.savez(tmp_path / "pyscf.npz", uc=np.diag([3.0, 3.1, 18.0]), tau=tau, r_int=r_int,
             hoppings=hoppings, overlap_orb=overlap_orb, position_orb=position_orb,
             fermi_energy=X_FILE, mesh=np.array([1, 1, 1]), edge_ratio=0.0, nelectron=2,
             basis="test")
    return lambda fermi_energy=None: load_pyscf_model(str(tmp_path / "pyscf.npz"), fermi_energy=fermi_energy)


def test_fermi_policy_pyscf(tmp_path):
    _assert_fermi_policy(_pyscf_file(tmp_path, np.random.default_rng(23)), x_file=X_FILE)


def _hamgnn_file(tmp_path, rng):
    from lsquant.models.hamgnn import FORMAT, load_hamgnn_model
    cells = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]
    n = 1                      # one atom, one s orbital, spinless
    h = _rng_blocks(rng, cells, n, 0.4)
    s = _rng_blocks(rng, cells, n, 0.05, identity=True)
    hamiltonian = np.concatenate([np.real(h[R]).ravel() for R in cells])
    overlap = np.concatenate([np.real(s[R]).ravel() for R in cells])
    np.savez(tmp_path / "h.npz", dft_code=np.array("unknown"), orbital_l=np.array([0]),
             format=np.array(FORMAT), source=np.array("prediction"), spin=np.array("none"),
             energy_unit=np.array("Hartree"), length_unit=np.array("Bohr"),
             lattice=np.diag([4.0, 4.1, 20.0]), positions=np.array([[0.0, 0.0, 10.0]]),
             atomic_numbers=np.array([14]), orbitals_per_atom=np.array([1]),
             orbital_index=np.array([0]), nao_max=np.array(1),
             cells=np.array(cells, int), atom_i=np.zeros(3, int), atom_j=np.zeros(3, int),
             hamiltonian=hamiltonian, overlap=overlap)
    return lambda fermi_energy=None: load_hamgnn_model(str(tmp_path / "h.npz"), fermi_energy=fermi_energy)


def test_fermi_policy_hamgnn(tmp_path):
    _assert_fermi_policy(_hamgnn_file(tmp_path, np.random.default_rng(24)), x_file=None)


def _deeptb_dir(tmp_path, rng):
    h5py = pytest.importorskip("h5py")
    from lsquant.models.deeptb import load_deeptb_model
    cell = np.diag([3.8, 3.8, 20.0])
    positions = np.array([[0.0, 0.0, 10.0], [1.9, 1.9, 10.3]])
    norb = np.array([2, 2])
    cells = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]

    def blocks(scale, identity=False):
        out = {}
        for R in cells:
            for i in range(2):
                for j in range(2):
                    key, partner = (R, i, j), (tuple(-x for x in R), j, i)
                    if key in out or partner in out:
                        continue
                    m = scale * rng.normal(size=(norb[i], norb[j]))
                    if key == partner:
                        m = 0.5 * (m + m.T) + (np.eye(norb[i]) if identity else 0)
                    out[key] = m
                    out[partner] = m.T
        return out

    def write_h5(name, b):
        with h5py.File(tmp_path / name, "w") as f:
            group = f.create_group("0")
            for (R, i, j), m in b.items():
                group.create_dataset(f"{i}_{j}_{R[0]}_{R[1]}_{R[2]}", data=m)

    write_h5("predicted_hamiltonians.h5", blocks(0.4))
    write_h5("predicted_overlaps.h5", blocks(0.05, identity=True))
    structure = {"basis_info": {"orbital_counts": {"Mo": 2, "S": 2}, "spinful": False},
                 "structure": {"cell": cell.tolist(), "positions": positions.tolist(),
                               "symbols": ["Mo", "S"]}}
    (tmp_path / "structure.json").write_text(json.dumps(structure))
    return lambda fermi_energy=None: load_deeptb_model(str(tmp_path), fermi_energy=fermi_energy)


def test_fermi_policy_deeptb(tmp_path):
    _assert_fermi_policy(_deeptb_dir(tmp_path, np.random.default_rng(25)), x_file=None)


def _abacus_dir(tmp_path, rng):
    from lsquant.models.abacus import load_abacus_model
    cells = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
    n = 2
    h = _rng_blocks(rng, cells, n, 0.4)
    s = _rng_blocks(rng, cells, n, 0.05, identity=True)
    (tmp_path / "STRU").write_text("""ATOMIC_SPECIES
Bi 208.98 Bi.upf

LATTICE_CONSTANT
7.5

LATTICE_VECTORS
1.0 0.0 0.0
0.5 0.8 0.0
0.0 0.0 4.0

ATOMIC_POSITIONS
Cartesian_angstrom
Bi
0.0
1
0.0 0.0 5.0 1 1 1
""")

    def write_csr(name, mats):
        out = ["STEP: 0", f"Matrix Dimension of H(R): {n}", f"Matrix number of H(R): {len(mats)}"]
        for R, m in zip(cells, mats):
            m = np.real(m)
            rows, cols = np.nonzero(m)
            out.append(f"{R[0]} {R[1]} {R[2]} {len(rows)}")
            out.append(" ".join(repr(float(v)) for v in m[rows, cols]))
            out.append(" ".join(map(str, cols)))
            out.append(" ".join(map(str, np.concatenate(
                [[0], np.cumsum(np.bincount(rows, minlength=n))]))))
        (tmp_path / name).write_text("\n".join(out) + "\n")

    write_csr("hr_tot", np.array([h[R] for R in cells]))
    write_csr("srs1_nao.csr", np.array([s[R] for R in cells]))
    return lambda fermi_energy=None: load_abacus_model(
        str(tmp_path / "hr_tot"), overlap_file=str(tmp_path / "srs1_nao.csr"),
        stru_file=str(tmp_path / "STRU"), orbitals_per_element={"Bi": n}, fermi_energy=fermi_energy)


def test_fermi_policy_abacus(tmp_path):
    _assert_fermi_policy(_abacus_dir(tmp_path, np.random.default_rng(26)), x_file=None)


def _gpaw_file(tmp_path, rng):
    from lsquant.models.gpaw import FORMAT, load_gpaw_model
    lattice = np.diag([3.0, 3.2, 15.0])
    positions = np.array([[0.1, 0.2, 7.0]])
    norb = np.array([2])
    mesh = (4, 3, 1)
    kpts = np.array([[i / mesh[0], j / mesh[1], 0.0]
                     for i in range(mesh[0]) for j in range(mesh[1])])
    cells = [(a, b, 0) for a in (-1, 0, 1) for b in (-1, 0, 1)]

    def bloch(mats):
        return np.array([sum(np.exp(2j * np.pi * np.dot(k, R)) * m for R, m in mats.items())
                         for k in kpts])

    h = _rng_blocks(rng, cells, 2, 0.3)
    s = _rng_blocks(rng, cells, 2, 0.02, identity=True)
    np.savez(tmp_path / "g.npz", format=np.array(FORMAT), kpts=kpts, mesh=np.array(mesh),
             H=np.conj(bloch(h)), S=np.conj(bloch(s)), soc=False, lattice=lattice,
             positions=positions, symbols=np.array(["Si"]), orbitals_per_atom=norb,
             fermi_energy=X_FILE)
    return lambda fermi_energy=None: load_gpaw_model(
        str(tmp_path / "g.npz"), soc=False, fermi_energy=fermi_energy)


def test_fermi_policy_gpaw(tmp_path):
    _assert_fermi_policy(_gpaw_file(tmp_path, np.random.default_rng(27)), x_file=X_FILE)


def _deeph_dir(tmp_path, rng):
    h5py = pytest.importorskip("h5py")
    from lsquant.models.deeph import load_deeph_model
    lattice = np.array([[3.0, 0.0, 0.0], [1.4, 2.6, 0.0], [0.0, 0.0, 15.0]])
    positions = np.array([[0.1, 0.2, 7.0], [1.6, 1.1, 7.4]])
    norb = np.array([2, 2])
    cells = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]
    np.savetxt(tmp_path / "lat.dat", lattice.T)
    np.savetxt(tmp_path / "site_positions.dat", positions.T)
    (tmp_path / "orbital_types.dat").write_text("0 1\n0\n")
    (tmp_path / "info.json").write_text(json.dumps({"isspinful": False, "fermi_level": X_FILE}))

    def blocks(scale, identity=False):
        out = {}
        for R in cells:
            for i in range(2):
                for j in range(2):
                    key, partner = (R, i, j), (tuple(-x for x in R), j, i)
                    if key in out or partner in out:
                        continue
                    m = scale * rng.normal(size=(norb[i], norb[j]))
                    if key == partner:
                        m = 0.5 * (m + m.T) + (np.eye(norb[i]) if identity else 0)
                    out[key] = m
                    out[partner] = m.T
        return out

    def write_h5(name, b):
        with h5py.File(tmp_path / name, "w") as f:
            for (R, i, j), m in b.items():
                f[str([R[0], R[1], R[2], i + 1, j + 1])] = m

    write_h5("hamiltonians_pred.h5", blocks(0.4))
    write_h5("overlaps.h5", blocks(0.05, identity=True))
    return lambda fermi_energy=None: load_deeph_model(str(tmp_path), fermi_energy=fermi_energy)


def test_fermi_policy_deeph(tmp_path):
    _assert_fermi_policy(_deeph_dir(tmp_path, np.random.default_rng(28)), x_file=X_FILE)


def _threebodytb_file(tmp_path, rng):
    from lsquant.models.threebodytb import load_threebodytb_model
    a = np.diag([4.6, 4.0, 30.0])      # rows, Bohr
    coords = np.array([[0.0, 0.0, 0.5], [1 / 3, 1 / 3, 0.52]])
    nwan = 3
    efermi_ry = -0.21
    cells = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
    h = _rng_blocks(rng, cells, nwan, 0.1)
    s = _rng_blocks(rng, cells, nwan, 0.02, identity=True)
    h = np.array([h[R] for R in cells])[None]          # nspin = 1
    s = np.array([s[R] for R in cells])

    def arr2str(x):
        return "\n".join(" ".join(f"{v: .16e}" for v in row) for row in np.atleast_2d(x)) + "\n"

    lines = []
    for r, R in enumerate(cells):
        for m in range(nwan):
            for n in range(nwan):
                lines.append(f"{r + 1:3d} {m + 1:3d} {n + 1:3d} "
                             f"{h[0, r, m, n].real: .10f} {h[0, r, m, n].imag: .10f} "
                             f"{s[r, m, n].real: .10f} {s[r, m, n].imag: .10f}")
    ind = "\n".join(" ".join(f"{x:7d}" for x in R) for R in cells)
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<root>
  <crystal>
    <A>{arr2str(a)}</A>
    <nat>2</nat>
    <coords>{arr2str(coords)}</coords>
    <types>{"  ".join(["C", "H"])}</types>
  </crystal>
  <nelec>5.0</nelec>
  <nspin>1</nspin>
  <dftenergy>0.0</dftenergy>
  <scf>false</scf>
  <efermi>{efermi_ry}</efermi>
  <energy>-1.0</energy>
  <tightbinding>
    <nr>{len(cells)}</nr>
    <nwan>{nwan}</nwan>
    <nspin>1</nspin>
    <nonorth>true</nonorth>
    <ind_arr>{ind}</ind_arr>
    <r_dict>Dict()</r_dict>
    <scf>false</scf>
    <h1>{arr2str(np.zeros((nwan, nwan)))}</h1>
    <scfspin>false</scfspin>
    <H>{chr(10).join(lines)}
</H>
  </tightbinding>
</root>
"""
    (tmp_path / "tb.xml").write_bytes(gzip.compress(xml.encode()))
    x_file = efermi_ry * 13.6056980659
    return (lambda fermi_energy=None: load_threebodytb_model(
        str(tmp_path / "tb.xml"), fermi_energy=fermi_energy)), x_file


def test_fermi_policy_threebodytb(tmp_path):
    load, x_file = _threebodytb_file(tmp_path, np.random.default_rng(29))
    _assert_fermi_policy(load, x_file=x_file)



def _openmx_file(tmp_path, rng):
    """A minimal OpenMX 3.9 scfout: two atoms, one orbital each (scalar
    blocks, as the layout of ``test_openmx_reader`` writes them)."""
    from lsquant.models.openmx import load_openmx_model
    cells = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]

    def herm_pair(build):
        m0 = build()
        m0 = 0.5 * (m0 + m0.conj())
        return {(0, 0, 0): m0, (1, 0, 0): build() * 0.3, (-1, 0, 0): (build() * 0.3).conj()}

    uu = herm_pair(lambda: rng.normal() + 1j * rng.normal())
    dd = herm_pair(lambda: rng.normal() + 1j * rng.normal())
    b = rng.normal() + 1j * rng.normal()
    ud = {(0, 0, 0): 0.2 * b, (1, 0, 0): 0.1 * b, (-1, 0, 0): 0.15 * b.conj()}
    s1 = rng.normal(scale=0.05)
    s = {(0, 0, 0): 1.0, (1, 0, 0): s1, (-1, 0, 0): s1}
    mu_hartree = -0.12
    natom, ncell = 2, 2
    tv = np.diag([7.0, 30.0, 30.0])
    gxyz = np.array([[0.5, 1.0, 0.0], [3.0, 2.0, 0.5]])
    nb = {a: [(a, 0)] + [(j, c) for c in range(3) for j in range(natom) if (j, c) != (a, 0)]
          for a in range(natom)}
    path = tmp_path / "toy.scfout"
    with open(path, "wb") as fh:
        w = lambda arr, dt: fh.write(np.asarray(arr, dt).tobytes())
        w([natom, 3 + 4 * 3, natom, 0, 0, ncell], "<i4")
        w([1], "<i4")                                    # order_max
        w(np.array([[0, 0, 0, 0], [0, tv[0, 0], tv[0, 1], tv[0, 2]],
                   [0, -tv[0, 0], -tv[0, 1], -tv[0, 2]]]), "<f8")
        w([[0, 0, 0, 0], [0, 1, 0, 0], [0, -1, 0, 0]], "<i4")
        w([1, 1], "<i4")                                 # Total_NumOrbs
        w([len(nb[a]) - 1 for a in range(natom)], "<i4")  # FNAN
        for a in range(natom):
            w([j + 1 for j, _ in nb[a]], "<i4")          # natn (1-based)
        for a in range(natom):
            w([c for _, c in nb[a]], "<i4")              # ncn
        w(np.c_[np.zeros(3), tv], "<f8")
        w(np.c_[np.zeros(3), 2 * np.pi * np.linalg.inv(tv).T], "<f8")
        w(np.c_[np.zeros(2), gxyz], "<f8")

        def element(m, a, j, c):
            return m[cells[c]]

        def write_all(fn):
            for a in range(natom):
                for j, c in nb[a]:
                    w([fn(a, j, c)], "<f8")
        im_split = 0.4
        write_all(lambda a, j, c: element(uu, a, j, c).real)
        write_all(lambda a, j, c: element(dd, a, j, c).real)
        write_all(lambda a, j, c: element(ud, a, j, c).real)
        write_all(lambda a, j, c: im_split * element(ud, a, j, c).imag)
        write_all(lambda a, j, c: element(uu, a, j, c).imag)
        write_all(lambda a, j, c: element(dd, a, j, c).imag)
        write_all(lambda a, j, c: (1 - im_split) * element(ud, a, j, c).imag)
        write_all(lambda a, j, c: element(s, a, j, c))
        for _ in range(3 + 4 + 2):                        # OLPpo, OLPmo, DM, iDM
            write_all(lambda a, j, c: 0.0)
        w([0], "<i4")                                    # eigenvalue method flag
        w([mu_hartree, 300.0, 0, 0, 0, 0, 0, 0, 2.0, 0.0], "<f8")
        w([0], "<i4")                                    # no input-file lines
    x_file = mu_hartree * 27.211386245988
    return (lambda fermi_energy=None: load_openmx_model(
        str(path), fermi_energy=fermi_energy)), x_file


def test_fermi_policy_openmx(tmp_path):
    load, x_file = _openmx_file(tmp_path, np.random.default_rng(30))
    _assert_fermi_policy(load, x_file=x_file)


def test_fermi_policy_siesta():
    sisl = pytest.importorskip("sisl", reason="the SIESTA reader needs sisl")
    del sisl  # the import is the skip check; the reader imports it itself
    from lsquant.models.siesta import load_siesta_model
    model_dir = os.path.join(REPO, "examples", "07_spin_hall_bi_bilayer", "siesta")
    if not os.path.isdir(model_dir):
        pytest.skip("the committed SIESTA example is absent")
    load = lambda fermi_energy=None: load_siesta_model(model_dir, fermi_energy=fermi_energy)
    m_default = load(None)
    printed = m_default.metadata["printed_fermi_level_eV"]
    assert m_default.metadata["fermi_energy_source"] == ("file" if printed is not None else "none")
    if printed is None:
        with pytest.raises(ValueError, match="no printed Fermi level"):
            load(X_NEW)
        return
    _assert_fermi_policy(load, x_file=printed)
