# Bring a model to LSQuant

LSQuant computes transport and spectral properties of a tight-binding model.
It does not run the electronic-structure calculation: you run SIESTA, OpenMX,
Quantum ESPRESSO + Wannier90 or PAOFLOW, turn on the outputs LSQuant reads,
and give LSQuant the result in a run file. The model is yours; LSQuant checks only that
its files are complete and consistent.

Every source follows the same three steps:

1. **Turn on** the outputs LSQuant reads in your code.
2. **Transfer** the model: for SIESTA, OpenMX and Wannier90 LSQuant reads the
   files of the run directly; for PAOFLOW one call in the PAOFLOW script writes a file.
3. **Run** LSQuant with a run file whose `model.kind` names the source.

| Source | 1. Turn on | 2. Transfer | 3. `model` in the run file | Guide |
|---|---|---|---|---|
| SIESTA | `Spin spin-orbit`; `TS.HS.Save true` (or `SaveHS true`) | nothing: LSQuant reads `<label>.fdf`, `<label>.TSHS` and, for the covariant velocity, the `<species>.ion.xml` files | `"kind": "siesta"`, `"dir"`, `"label"` | [Use a SIESTA calculation](use-siesta.md) |
| OpenMX | `scf.SpinPolarization NC`; `scf.SpinOrbit.Coupling on`; `HS.fileout on` | nothing: LSQuant reads `<System.Name>.scfout` | `"kind": "openmx"`, `"file"` | [Use an OpenMX calculation](use-openmx.md) |
| Quantum ESPRESSO + Wannier90 | in `<seed>.win`: `write_hr`, `write_xyz`, `write_rmn`, `write_u_matrices`; in `pw2wannier90.x`: `write_spn` | nothing: LSQuant reads the Wannier90 directory and writes `<seed>.uc` and `<seed>.xyz` on the first load | `"kind": "wannier"`, `"dir"`, `"seed"`, `params.velocity_mode`, `params.spin`, `params.fermi_energy` | [Use a Wannier90 model](use-wannier90.md) |
| Quantum ESPRESSO + PAOFLOW | `projwfc.x` with `lwrite_overlaps = .true.` | `write_paoflow_model(p, "model.npz")` at the end of the PAOFLOW script | `"kind": "paoflow"`, `"file"` | [Use a PAOFLOW model](use-paoflow.md) |

For all of them, step 3 is the same:

```bash
python -m lsquant inspect run.json   # check the run file
python -m lsquant run run.json       # compute
```

All sources go through the same executor: the same velocity modes,
spin operators, spectral projection, kernel polynomial method and exact check.
The quantities they give are listed in the
[README](../../README.md#quantities-by-model-source).
