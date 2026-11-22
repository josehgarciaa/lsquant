# Use a PAOFLOW model

Start from a Quantum ESPRESSO calculation with spin–orbit coupling, build the
PAO tight-binding model with PAOFLOW, write it to a file LSQuant reads, and
compute the spin Hall conductivity with the kernel polynomial method (KPM).
You should know PAOFLOW and have read [Running LSQuant](run-lsquant.md).

A PAOFLOW model carries its own spin operator and runs through the same
executor as a SIESTA model. Every key is listed in the
[run-file reference](../reference/run-file.md#wannier90-and-paoflow-models-in-the-metric-executors).
[Example 20](../../examples/20_pdse2_spin_hall_routes/README.md) uses a
PAOFLOW model of a PdSe₂ monolayer.

## Requirements

- LSQuant built and importable; see [Install LSQuant](install.md).
- PAOFLOW, in the same Python environment as LSQuant.
- Quantum ESPRESSO with `projwfc.x`.
<!-- tested with PAOFLOW 3.0.0 and Quantum ESPRESSO 7.4 -->

## The DFT calculation

1. Run `pw.x` with spin–orbit coupling: an scf calculation, then an nscf
   calculation on a uniform k mesh (`K_POINTS automatic`) with enough bands
   above the energies you need.

2. Run `projwfc.x` with the overlaps:

   ```text
   &PROJWFC
     prefix = "pdse2", outdir = "./tmp"
     lwrite_overlaps = .true.
   /
   ```

PAOFLOW projects the Bloch states on the atomic wavefunctions stored in the
pseudopotential files. The PAO model is only as good as that basis: the
conduction bands need empty orbitals in it (for PdSe₂, Pd 5s and 5p, Se 4d).
Pseudopotentials whose files hold only the occupied or semicore shells give a
model that reproduces the DFT bands on the k mesh but not between the k
points; see [Check the model](#check-the-model).

## Write the model for LSQuant

At the end of the PAOFLOW script, after `pao_hamiltonian()`:

```python
from PAOFLOW import PAOFLOW
from lsquant.models.orthogonal import write_paoflow_model

p = PAOFLOW.PAOFLOW(savedir="tmp/pdse2.save", outputdir="out", npool=1)
p.read_atomic_proj_QE()
p.projectability(pthr=0.90)
p.pao_hamiltonian()
write_paoflow_model(p, "pdse2_paoflow.npz")
```

The file holds H(R), the orbital positions and PAOFLOW's on-site spin
operator, with energies measured from the DFT Fermi energy.

### Fully relativistic pseudopotentials

PAOFLOW builds its spin operator from the order of the two radial channels
j = l − ½ and j = l + ½ of each l shell, and by default assumes that j = l − ½
comes first. Look at the `PP_RELWFC` lines of each UPF file. If `jchi = l + ½`
comes first, as in PseudoDojo files, give the order explicitly: one entry per
l shell of every atom in basis order, with the j of the first channel.

```python
# Pd: s, p, d; Se: d, s, p; first channel j = l + 1/2
shells = {"Pd": ([0, 1, 2], [0.5, 1.5, 2.5]), "Se": ([2, 0, 1], [2.5, 0.5, 1.5])}
atoms = p.data_controller.data_dicts()[0]["atoms"]
sh_l = sum((shells[a][0] for a in atoms), [])
sh_j = sum((shells[a][1] for a in atoms), [])
write_paoflow_model(p, "pdse2_paoflow.npz", sh_l=sh_l, sh_j=sh_j)
```

With the wrong order the spin operator is wrong, and the spin Hall
conductivity is not zero in a trivial gap.

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

| Key | What to set |
|---|---|
| `velocity_mode` | `"berry_connection"`: the velocity of PAOFLOW's `gradient_and_momenta`. `"bare"` drops the orbital positions. PAOFLOW has no position matrix, so `"covariant"` is refused. |
| `projection.window` | Measured from the DFT Fermi energy. Keep it below the energy where PAOFLOW shifted the bands it could not project (the `shift` printed by `projectability`). |

`model.params.spin` and `model.params.fermi_energy` are refused: the spin and
the energy origin come from the file.

## Check the model

- **Locality.** `write_paoflow_model` stores `edge_ratio`, the largest |H(R)|
  on the outermost R shell of the PAOFLOW grid over the largest |H(R)| with
  R ≠ 0, and warns above 10⁻². A large value means H(R) has not decayed:
  between the DFT k points H(k) can hold states that are not in the DFT
  bands, and they show up as a spin Hall signal inside a gap. Use a denser
  nscf k mesh or a projection basis with the missing empty orbitals.
- **Bands.** Compare the bands of the model with the DFT bands, at the DFT
  k points and between them.
- **Trivial gap.** In a gap with zero spin Chern number the spin Hall
  conductivity must be close to zero.
- **PAOFLOW itself.** PAOFLOW's `spin_Hall` adds a regulator δ to the
  interband denominator, (E_n − E_m)² + δ², and uses adaptive smearing; the
  LSQuant oracle uses δ = 0 and a fixed Gaussian, and KPM broadens through its
  kernel. In a gap all three agree. Between bands compare them only with the
  same broadening: with δ = 0 those values depend strongly on the k mesh.
