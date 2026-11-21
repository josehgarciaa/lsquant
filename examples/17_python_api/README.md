# Example 17: LSQuant from Python, without a run file

## The physical question

A run file covers the standard calculations. Sometimes we want to change one
step: build our own operator, choose our own random or unit vectors, loop over
a parameter, or post-process the moments in a script. The run file is a thin
layer over a Python API, and we can call that API at two levels. We compute the
DOS and the Kubo–Greenwood $\sigma_{xx}(E)$ of graphene both ways, check that
they agree, and read the DOS against the known shape of the graphene bands.

## The model

The model is the nearest-neighbour graphene $p_z$ model of
[example 10](../10_bases_compare/README.md), with $t = -2.7$ eV, stored in
Wannier90 format in `examples/10_bases_compare/graphene/`. We use a
$16 \times 16$ k mesh in the `bloch` basis, an exact trace, 128 moments, the
window $[-9, 9]$ eV and a 150 meV broadening.

## Run 1: a run file as a Python dictionary

`execute_runfile` accepts a path or a dictionary with the same keys as the JSON
file. It validates the dictionary with the same strict checks, runs, writes the
same files, and returns the list of written paths:

```python
from lsquant.runfile import execute_runfile

run = {
    "model": {"kind": "wannier", "dir": MODEL_DIR, "seed": "graphene"},
    "basis": {"name": "bloch", "supercell": [16, 16, 1]},
    "metric": {"policy": "orthogonal"},
    "word": "Tr[delta(E-H)]",
    "words": ["Tr[delta(E-H)]", "Tr[VX delta(E-H) VX delta(E-H)]"],
    "expansion": {"moments": 128, "window": [-9.0, 9.0], "trace": "exact"},
    "reconstruction": {"kernel": "jackson", "broadening_meV": 150.0,
                       "formulas": ["dos", "greenwood"]},
    "output": {"dir": out_dir, "legacy_layout": False, "write_chebmom": True},
}
result = execute_runfile(run)        # result["moment_files"], result["curve_files"], ...
```

This is the level for parameter scans. We build the dictionary in a loop and
change, for example, `basis.disorder.W` or `expansion.moments`.

## Run 2: the building blocks

We now do the same calculation step by step. We load the model, build the
basis, and take its Hamiltonian, velocity operator and exact-trace unit
vectors:

```python
from lsquant import _core, words
from lsquant.bases import make_basis
from lsquant.models.wannier import load_wannier_model

model = load_wannier_model(MODEL_DIR, seed="graphene")
basis = make_basis("bloch", model, (16, 16, 1))
window = _core.SpectralWindow(-9.0, 9.0)
states = basis.states(trace_mode="exact")          # N_orb = 2 block-local unit vectors
operators = {"VX": basis.operator("VX")}
```

We compile the word, bound to these operators, into a plan. The plan is the
list of operations that the C++ engine executes:

```python
compiled = words.compile(
    "Tr[VX delta(E-H) VX delta(E-H)]",
    basis={"name": "bloch", "operators": operators},
    params={"expansion": {"moments": 128, "window": window,
                          "states": states.number_of_states}},
)
```

We execute the plan, which returns the Chebyshev moments, and reconstruct the
curve. The two lines marked `# (*)` are the ones that `execute_runfile` adds
for an exact Bloch trace:

```python
from lsquant.moments import run_plan
from lsquant.reconstruct import kubo_greenwood

mu = run_plan(compiled.plan, basis.hamiltonian(), window, states,
              operators=compiled.operators, batch_size=128,
              symmetrise=True)                     # (*) Kubo-Greenwood moments are symmetrised
mu.number_of_states = basis.dimension              # (*) exact trace over all N_k N_orb states
energies, sigma_xx = kubo_greenwood(mu, 150.0)
```

The trace runs over $N_\mathrm{orb}$ block-local vectors, but the
normalisation must count all $N_k N_\mathrm{orb}$ states. Kubo–Greenwood
moments are also symmetrised. For the DOS, the word is `Tr[delta(E-H)]`, the
plan runs unsymmetrised, and the reconstruction is `dos_curve`.

The full script, [`python_api.py`](python_api.py), runs both levels and
compares them:

```bash
python3 examples/17_python_api/python_api.py /tmp/ex17/level1
```

## The result

The script prints:

```text
level 1: wrote 6 files to /tmp/ex17/level1
level 2: Tr[delta(E-H)]                     plan = 2 instructions, moments [128], curve 3840 points
level 2: Tr[VX delta(E-H) VX delta(E-H)]    plan = 6 instructions, moments [128, 128], curve 127 points
dos   level 1 vs level 2: max rel |d mu| = 0.00e+00, max rel |d curve| = 0.00e+00  PASS
kubo  level 1 vs level 2: max rel |d mu| = 0.00e+00, max rel |d curve| = 0.00e+00  PASS
graphene DOS per orbital at E = 0: 0.0001 1/eV; at E = |t| = 2.7 eV (van Hove): 0.1773 1/eV
CHECK PASS
```

The two levels give bit-identical moments and curves. The DOS has the
graphene shape: it vanishes at the Dirac point $E = 0$ and peaks at the van Hove
singularity $E = |t| = 2.7$ eV. The DOS curve has $30 M$ points, while the
Kubo–Greenwood reconstruction keeps $M - 1$ points.

## Runtime

On a shared 96-core node with one OpenMP thread, both levels together took
3 s.
