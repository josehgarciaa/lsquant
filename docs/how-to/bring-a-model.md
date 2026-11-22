# Bring a model to LSQuant

LSQuant computes transport and spectral properties of a tight-binding model
that you build with SIESTA, OpenMX, Quantum ESPRESSO + Wannier90, PAOFLOW or
a machine-learning Hamiltonian code. Turn on the outputs LSQuant reads in that
code, then give LSQuant the result in a run file. You choose and converge the
model in your own code, as usual.
LSQuant checks only that the files are complete and consistent.

Every source follows the same three steps:

1. **Turn on** the outputs LSQuant reads in your code.
2. **Transfer** the model: for SIESTA, OpenMX, Wannier90 and the DeepH format,
   LSQuant reads the output files of your calculation directly; for PAOFLOW
   and DeePTB, one call in a Python script writes them; for HamGNN, a script
   that comes with LSQuant writes one file in the HamGNN environment; ThreeBodyTB writes its
   model file with its own `write_tb_crys`; NextHAM's own scripts
   write H(R) in the ABACUS format, which LSQuant reads with the ABACUS overlap
   and `STRU`.
3. **Run** LSQuant with a run file whose `model.kind` names the source.

| Source | 1. Outputs to turn on | 2. How LSQuant gets the model | 3. `model` in the run file | Guide |
|---|---|---|---|---|
| SIESTA | `Spin spin-orbit`; `TS.HS.Save true` (or `SaveHS true`) | nothing to run: LSQuant reads `<label>.fdf`, `<label>.TSHS` and, for the covariant velocity, the `<species>.ion.xml` files | `"kind": "siesta"`, `"dir"`, `"label"` | [Use a SIESTA calculation](use-siesta.md) |
| OpenMX | `scf.SpinPolarization NC`; `scf.SpinOrbit.Coupling on`; `HS.fileout on` | nothing to run: LSQuant reads `<System.Name>.scfout` | `"kind": "openmx"`, `"file"` | [Use an OpenMX calculation](use-openmx.md) |
| Quantum ESPRESSO + Wannier90 | in `<seed>.win`: `write_hr`, `write_xyz`, `write_rmn`, `write_u_matrices`; in `pw2wannier90.x`: `write_spn` | nothing to run: LSQuant reads the Wannier90 directory and, the first time, writes `<seed>.uc` and `<seed>.xyz` there | `"kind": "wannier"`, `"dir"`, `"seed"`, `params.velocity_mode`, `params.spin`, `params.fermi_energy` | [Use a Wannier90 model](use-wannier90.md) |
| Quantum ESPRESSO + PAOFLOW | `projwfc.x` with `lwrite_overlaps = .true.` | `write_paoflow_model(p, "model.npz")` at the end of the PAOFLOW script | `"kind": "paoflow"`, `"file"` | [Use a PAOFLOW model](use-paoflow.md) |
| DeepH-pack, DeepH-E3, xDeepH, MACE-H | inference with the prediction written next to the overlap and structure files | nothing to run: LSQuant reads `hamiltonians_pred.h5` (or `hamiltonian_pred.h5`), the overlap and the structure files | `"kind": "deeph"`, `"dir"`, `params.fermi_energy` | [Use a DeepH-format model](use-deeph.md) |
| DeePTB (SK and E3, spinless) | a trained model and a structure | `TBSystem(...).to_pardiso_json("dir")` or `dptb pdso` | `"kind": "deeptb"`, `"dir"`, `params.fermi_energy` | [Use a DeePTB model](use-deeptb.md) |
| HamGNN | a prediction run (`stage: test`) on the `graph_data.npz` of the structure | `python lsquant/models/hamgnn.py graph_data.npz prediction_hamiltonian.npy model.npz --nao-max N --spin soc --dft-code openmx` in the HamGNN environment | `"kind": "hamgnn"`, `"file"`, `params.fermi_energy` | [Use a HamGNN model](use-hamgnn.md) |
| ThreeBodyTB | nothing: build the model (`scf_energy` for charge transfer) | `write_tb_crys("model.xml", tbc)` in Julia (a ThreeBodyTB function) | `"kind": "threebodytb"`, `"file"` | [Use a ThreeBodyTB model](use-threebodytb.md) |
| NextHAM (ABACUS basis) | save the predictions of the test run (`test.py`) | NextHAM's `run_bands_compare.py` writes `hr_tot`; LSQuant reads it with ABACUS's `srs1_nao.csr` and `STRU` | `"kind": "abacus"`, `"file"`, `params.overlap_file`, `params.stru_file`, `params.fermi_energy` | [Use a NextHAM model](use-nextham.md) |
| GPAW (LCAO) | `mode='lcao'`, `symmetry='off'`, `parallel={'domain': 1}`; `calc.write('calc.gpw')` | `python lsquant/models/gpaw.py calc.gpw model.npz` in the GPAW environment | `"kind": "gpaw"`, `"file"` | [Use a GPAW calculation](use-gpaw.md) |

For every source, step 3 is the same. Check the run file, then compute:

```bash
python -m lsquant inspect run.json   # check the run file
python -m lsquant run run.json       # compute
```

All sources use the same velocity modes, spin operators, spectral projection and exact check.
The kernel polynomial method (KPM) and its convergence parameters are also the
same; [Use a SIESTA calculation](use-siesta.md#number-of-moments-and-broadening)
explains them. A machine-learning model has no position matrix, so its
velocity is `bare` or `berry_connection`. A model without spin–orbit coupling
has no spin operator and gives the density of states only.

OpenMX, PAOFLOW, DeepH, DeePTB, HamGNN, ThreeBodyTB and NextHAM (ABACUS) models, and Wannier90 models with
`params.velocity_mode`, `params.spin` or `params.fermi_energy`, run as clean
crystals on a k-mesh (`basis.name` `"bloch"` or `"eigen"`). A real-space supercell or
Anderson disorder stops the run; disorder is available for SIESTA models.
<!-- source: python/lsquant/runfile.py:1576–1579 ("runs in the clean bloch or eigen basis" when basis.name = "real" or basis.disorder.W != 0); SIESTA disorder: use-siesta.md "Disordered supercells" -->

With a spectral projection, the cost of a run grows with the number of bands
of the model at each k-point, so a small Wannier90 model is much cheaper than
an atomic-orbital model of the same material. The OpenMX, PAOFLOW and
Wannier90 guides give measured costs.
<!-- source: measured, examples/20_pdse2_spin_hall_routes README Runtime / data/wannier_covariant_spn.json, AMD EPYC 9374F (12-band Wannier90 model about 5 min on 8 threads; 108-band PAOFLOW model about 2 h 30 min on 8 threads) -->

The spin Hall conductivity σ^z_xy is written in units of e²/h, with the
spin current defined by the Pauli matrix σ_z. To compare with values for
spin ħ/2, in units of (ħ/e) e²/h, divide by 2; see
[Kubo formulas](../explanation/kubo-formulas.md#spin-currents-and-the-spin-hall-conductivity).
<!-- source: docs/explanation/kubo-formulas.md (Spin currents and the spin Hall conductivity); reconstruct.py:731-743 -->

Check a predicted Hamiltonian against DFT before you use it for transport.
The [table in the README](../../README.md#quantities-by-model-source) lists
the quantities each source gives.
