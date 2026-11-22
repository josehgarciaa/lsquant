# Use a PAOFLOW model

Feed a PAOFLOW tight-binding model to LSQuant and compute its spin Hall
conductivity with the kernel polynomial method (KPM). Add one call at the end
of your PAOFLOW script, then point a run file at the file it writes. LSQuant
reads the model as it is; the quality of the PAO model (projectability, basis,
k mesh) is yours to check. You should know PAOFLOW and have read
[Running LSQuant](run-lsquant.md).

The run file sends the PAOFLOW model through the same executor as a SIESTA
model: same KPM engine, projection and exact check. Every key is listed in the
[run-file reference](../reference/run-file.md#wannier90-and-paoflow-models-in-the-metric-executors).
[Example 20](../../examples/20_pdse2_spin_hall_routes/README.md) uses it for a
PdSe₂ monolayer.

## Requirements

- LSQuant built and importable in the Python environment of PAOFLOW; see
  [Install LSQuant](install.md).
<!-- tested with PAOFLOW 3.0.0 and Quantum ESPRESSO 7.4 -->

## What to turn on

Run your usual PAOFLOW workflow with spin–orbit coupling (`pw.x` scf and nscf,
`projwfc.x` with `lwrite_overlaps = .true.`). At the end of the PAOFLOW
script, after `pao_hamiltonian()`, write the model:

```python
from lsquant.models.orthogonal import write_paoflow_model

p.read_atomic_proj_QE()
p.projectability(pthr=0.90)
p.pao_hamiltonian()
write_paoflow_model(p, "pdse2_paoflow.npz")
```

The file holds H(R), the orbital positions and PAOFLOW's on-site spin
operator, with energies measured from the DFT Fermi energy.

**Fully relativistic pseudopotentials.** PAOFLOW's spin operator assumes that
in each l shell the radial channel j = l − ½ comes before j = l + ½ in the
pseudopotential file (`PP_RELWFC` lines of the UPF). If your files store
j = l + ½ first, as PseudoDojo files do, pass the order: one entry per l
shell of every atom in basis order, with the j of the first channel.

```python
# Pd: s, p, d; Se: d, s, p; first channel j = l + 1/2
shells = {"Pd": ([0, 1, 2], [0.5, 1.5, 2.5]), "Se": ([2, 0, 1], [2.5, 0.5, 1.5])}
atoms = p.data_controller.data_dicts()[0]["atoms"]
sh_l = sum((shells[a][0] for a in atoms), [])
sh_j = sum((shells[a][1] for a in atoms), [])
write_paoflow_model(p, "pdse2_paoflow.npz", sh_l=sh_l, sh_j=sh_j)
```

## The run file

```json
{
  "model": {
    "kind": "paoflow",
    "file": "path/to/pdse2_paoflow.npz",
    "params": {"velocity_mode": "berry_connection", "oracle": true}
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "metric": {"policy": "orthogonal"},
  "projection": {"kind": "shift", "window": [-3, 3]},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 6.0,
                     "energy_grid": [-0.5, 2.0, 1251]},
  "output": {"dir": "out/paoflow_shc"}
}
```

| Key | Value |
|---|---|
| `velocity_mode` | `"berry_connection"`, the velocity of PAOFLOW's `gradient_and_momenta`, or `"bare"`. PAOFLOW has no position matrix, so `"covariant"` is refused. |
| `projection.window` | Measured from the DFT Fermi energy. Keep it below the `shift` that PAOFLOW's `projectability` prints. |

`model.params.spin` and `model.params.fermi_energy` are refused: both come
from the file.

## What LSQuant checks

- It stops when H(k) is not Hermitian or the file lacks a field.
- `write_paoflow_model` warns when H(R) has not decayed at the edge of the
  PAOFLOW k grid (`edge_ratio` above 10⁻², stored in the file). Between the
  DFT k points such a model can hold states that the DFT bands do not have.

It does not check the PAO model itself: the projectability, the projection
basis, the k mesh or the agreement with the DFT bands.
