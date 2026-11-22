# Use a PAOFLOW model

Connect a PAOFLOW tight-binding model to LSQuant and compute its spin Hall
conductivity. LSQuant reads the model as PAOFLOW built it; the PAO model
itself (projectability, basis, k mesh, agreement with the DFT bands) is yours.
You should know PAOFLOW and have read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

Every key is listed in the
[run-file reference](../reference/run-file.md#wannier90-and-paoflow-models-in-the-metric-executors).

## Requirements

- LSQuant built and importable in the Python environment of PAOFLOW; see
  [Install LSQuant](install.md).
<!-- tested with PAOFLOW 3.0.0 and Quantum ESPRESSO 7.4 -->

## 1. Turn on in Quantum ESPRESSO

Run your usual PAOFLOW workflow with spin–orbit coupling: `pw.x` scf and nscf,
then `projwfc.x` with the overlaps:

```text
&PROJWFC
  prefix = "myprefix", outdir = "./tmp"
  lwrite_overlaps = .true.
/
```

## 2. Transfer the model

At the end of the PAOFLOW script, after `pao_hamiltonian()`, write the model:

```python
from lsquant.models.orthogonal import write_paoflow_model

p.read_atomic_proj_QE()
p.projectability(pthr=0.90)
p.pao_hamiltonian()
write_paoflow_model(p, "mymodel.npz")
```

The file holds H(R), the orbital positions and PAOFLOW's on-site spin
operator, with energies measured from the DFT Fermi energy.

**Fully relativistic pseudopotentials.** PAOFLOW's spin operator assumes that
in each l shell the radial channel j = l − ½ comes before j = l + ½ in the
pseudopotential file (`PP_RELWFC` lines of the UPF). If your files store
j = l + ½ first, as PseudoDojo files do, pass the order: one entry per l
shell of every atom in basis order, with the j of the first channel.

```python
# Example: Pd with shells s, p, d and Se with shells d, s, p, first channel j = l + 1/2
shells = {"Pd": ([0, 1, 2], [0.5, 1.5, 2.5]), "Se": ([2, 0, 1], [2.5, 0.5, 1.5])}
atoms = p.data_controller.data_dicts()[0]["atoms"]
sh_l = sum((shells[a][0] for a in atoms), [])
sh_j = sum((shells[a][1] for a in atoms), [])
write_paoflow_model(p, "mymodel.npz", sh_l=sh_l, sh_j=sh_j)
```

## 3. Run LSQuant

```json
{
  "model": {
    "kind": "paoflow",
    "file": "/path/to/mymodel.npz",
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

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Value |
|---|---|
| `velocity_mode` | `"berry_connection"`, the velocity of PAOFLOW's `gradient_and_momenta`, or `"bare"`. PAOFLOW has no position matrix, so `"covariant"` is refused. |
| `projection.window` | Measured from the DFT Fermi energy. |

`model.params.spin` and `model.params.fermi_energy` are refused: both come
from the file. LSQuant stops when H(k) is not Hermitian or the file lacks a
field.

## 4. Recommendations

- Use pseudopotentials whose atomic wavefunctions include the empty orbitals
  of the conduction bands you need (for example a transition-metal s and p
  shell above the d shell). With only occupied or semicore shells the model
  reproduces the DFT bands at the DFT k points but not between them.
- `write_paoflow_model` warns when H(R) has not decayed at the edge of the
  PAOFLOW k grid (`edge_ratio` above 10⁻², stored in the file). Such a model
  can hold states between the DFT k points that the DFT bands do not have;
  use a denser nscf k mesh or a richer projection basis.
- Keep `projection.window` below the `shift` that PAOFLOW's `projectability`
  prints: above it the bands were moved by PAOFLOW.
- In a trivial gap the spin Hall conductivity must be close to zero. A clear
  signal there points to the spin operator (channel order above) or to
  spurious states.
