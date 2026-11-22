# Bring a model to LSQuant

LSQuant computes transport and spectral properties of a tight-binding model.
It does not run the electronic-structure calculation: you run SIESTA, OpenMX,
Quantum ESPRESSO + Wannier90, PAOFLOW or a machine-learning Hamiltonian code,
turn on the outputs LSQuant reads, and give LSQuant the result in a run file.
The model is yours; LSQuant checks only that its files are complete and
consistent.

Every source follows the same three steps:

1. **Turn on** the outputs LSQuant reads in your code.
2. **Transfer** the model: for SIESTA, OpenMX, Wannier90 and the DeepH format
   LSQuant reads the files of the run directly; for PAOFLOW and DeePTB one
   call in a Python script writes the files; for HamGNN a script that comes with
   LSQuant writes a file in the HamGNN environment.
3. **Run** LSQuant with a run file whose `model.kind` names the source.

| Source | 1. Turn on | 2. Transfer | 3. `model` in the run file | Guide |
|---|---|---|---|---|
| SIESTA | `Spin spin-orbit`; `TS.HS.Save true` (or `SaveHS true`) | nothing: LSQuant reads `<label>.fdf`, `<label>.TSHS` and, for the covariant velocity, the `<species>.ion.xml` files | `"kind": "siesta"`, `"dir"`, `"label"` | [Use a SIESTA calculation](use-siesta.md) |
| OpenMX | `scf.SpinPolarization NC`; `scf.SpinOrbit.Coupling on`; `HS.fileout on` | nothing: LSQuant reads `<System.Name>.scfout` | `"kind": "openmx"`, `"file"` | [Use an OpenMX calculation](use-openmx.md) |
| Quantum ESPRESSO + Wannier90 | in `<seed>.win`: `write_hr`, `write_xyz`, `write_rmn`, `write_u_matrices`; in `pw2wannier90.x`: `write_spn` | nothing: LSQuant reads the Wannier90 directory and writes `<seed>.uc` and `<seed>.xyz` on the first load | `"kind": "wannier"`, `"dir"`, `"seed"`, `params.velocity_mode`, `params.spin`, `params.fermi_energy` | [Use a Wannier90 model](use-wannier90.md) |
| Quantum ESPRESSO + PAOFLOW | `projwfc.x` with `lwrite_overlaps = .true.` | `write_paoflow_model(p, "model.npz")` at the end of the PAOFLOW script | `"kind": "paoflow"`, `"file"` | [Use a PAOFLOW model](use-paoflow.md) |
| DeepH-pack, DeepH-E3, xDeepH, MACE-H | inference with the prediction written next to the overlap and structure files | nothing: LSQuant reads `hamiltonians_pred.h5` (or `hamiltonian_pred.h5`), the overlap and the structure files | `"kind": "deeph"`, `"dir"`, `params.fermi_energy` | [Use a DeepH-format model](use-deeph.md) |
| DeePTB (SK and E3, spinless) | a trained model and a structure | `TBSystem(...).to_pardiso_json("dir")` or `dptb pdso` | `"kind": "deeptb"`, `"dir"`, `params.fermi_energy` | [Use a DeePTB model](use-deeptb.md) |
| HamGNN | a prediction run (`stage: test`) on the `graph_data.npz` of the structure | `python lsquant/models/hamgnn.py graph_data.npz prediction_hamiltonian.npy model.npz --nao-max N --spin soc --dft-code openmx` in the HamGNN environment | `"kind": "hamgnn"`, `"file"`, `params.fermi_energy` | [Use a HamGNN model](use-hamgnn.md) |

For every source, step 3 is the same:

```bash
python -m lsquant inspect run.json   # check the run file
python -m lsquant run run.json       # compute
```

All sources go through the same executor: the same velocity modes, spin
operators, spectral projection, kernel polynomial method and exact check. A
machine-learning model has no position matrix, so its velocity is `bare` or
`berry_connection`, and a model without spin–orbit coupling has no spin
operator. LSQuant computes transport only; the quality of a predicted
Hamiltonian is yours to check against DFT.
The quantities they give are listed in the
[README](../../README.md#quantities-by-model-source).
