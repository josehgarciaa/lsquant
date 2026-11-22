# Use a PAOFLOW model

Start from your own PAOFLOW tight-binding model and end with its spin Hall
conductivity. LSQuant reads the model as PAOFLOW built it. Check the
projectability, the basis, the k-mesh and the agreement with the DFT bands in
PAOFLOW, as usual. You should know PAOFLOW and have read
[Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

LSQuant computes the density of states and the spin Hall conductivity of the
clean crystal on a k-mesh; a real-space supercell or disorder stops the run.
Every key is listed in the
[run-file reference](../reference/run-file.md#wannier90-and-paoflow-models-in-the-metric-executors).
<!-- source: python/lsquant/runfile.py:1576–1579; words: reference/run-file.md (SIESTA executors, Words row) -->

## Requirements

You need:

- LSQuant built and importable in the Python environment of PAOFLOW; see
  [Install LSQuant](install.md).
<!-- tested with PAOFLOW 3.0.0 and Quantum ESPRESSO 7.4 -->

## 1. Turn on the overlaps in Quantum ESPRESSO

Run your usual PAOFLOW workflow with spin–orbit coupling, `pw.x` scf and
nscf. Example 20 sets `noncolin = .true.` and `lspinorb = .true.` in both
runs, with fully relativistic pseudopotentials. Then run `projwfc.x` with the
overlaps written:
<!-- source: examples/20_pdse2_spin_hall_routes/qe/scf.in:20–21, paoflow/nscf_pao.in:20–21; README (rel-pbe pseudopotentials) -->

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

`pthr` is PAOFLOW's projectability threshold: only bands with a
projectability above 0.90 enter the model, as in example 20.
The file contains H(R), the orbital positions and PAOFLOW's on-site spin
operator. Its energies are measured from the DFT Fermi energy, so give every
energy in the run file relative to it, E_F = 0.
<!-- source: examples/20_pdse2_spin_hall_routes/README.md (PAOFLOW route: "92 bands above a projectability of 0.90"); reference/run-file.md (paoflow file contents) -->

### Fully relativistic pseudopotentials

PAOFLOW builds its spin operator on the assumption that, in each l shell, the
radial channel j = l − ½ comes before j = l + ½ in the pseudopotential file
(`PP_RELWFC` lines of the UPF). If your files store j = l + ½ first, as
PseudoDojo files do, pass the channel order to `write_paoflow_model` as
`sh_l` and `sh_j`: one entry per l shell of every atom, in basis order, with
the j of the first channel.

```python
# Example: Pd with shells s, p, d and Se with shells d, s, p, first channel j = l + 1/2
shells = {"Pd": ([0, 1, 2], [0.5, 1.5, 2.5]), "Se": ([2, 0, 1], [2.5, 0.5, 1.5])}
atoms = p.data_controller.data_dicts()[0]["atoms"]
sh_l = sum((shells[a][0] for a in atoms), [])
sh_j = sum((shells[a][1] for a in atoms), [])
write_paoflow_model(p, "mymodel.npz", sh_l=sh_l, sh_j=sh_j)
```

## 3. Run LSQuant

This run file computes the spin Hall conductivity on a 24 × 24 k-mesh, with
the spectrum shifted inside [−3, 3] eV:

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

The template has the settings of the projected SIESTA run of
[example 20](../../examples/20_pdse2_spin_hall_routes/README.md). It is the
exact reference path of LSQuant, and each setting does the following:

- `word` is the observable, here the spin Hall conductivity σ^z_xy: the spin
  current with spin along z and flow along x, in response to an electric
  field along y. See [Observable words](../reference/observable-words.md).
- `basis.name = "eigen"` diagonalises H(k) at every point of the k-mesh, and
  `basis.supercell` is that mesh, Γ-centred, with 1 as the third entry for a
  2D system. With `expansion.trace = "exact"`, LSQuant sums over all states,
  so the curve has no random-vector noise.
- `metric.policy = "orthogonal"`: a PAOFLOW model has no overlap matrix.
- `projection` with `"kind": "shift"` moves every state outside [−3, 3] eV,
  relative to E_F, to the nearest edge of the window. The expansion then
  covers at most 7.2 eV instead of the full spectrum, and needs fewer moments.
- `expansion.moments` sets the energy resolution, as the smearing does for a
  DFT density of states. With `"window": "eigen"`, LSQuant pads the projected
  spectrum by 10 %, to at most [−3.6, 3.6] eV. With 1024 moments, the finest
  broadening is πa/M ≈ 11 meV, with a = 3.6 eV the half-width of the window.
  LSQuant applies this value in place of the requested 6 meV; to resolve
  6 meV, use about 1900 moments. See
  [Number of moments and broadening](use-siesta.md#number-of-moments-and-broadening).
- `reconstruction.energy_grid`, `[E_min, E_max, n]` in eV relative to E_F,
  is the grid of the exact reference: here −0.5 to 2.0 eV in 1251 points,
  2 meV apart. The KPM curve has its own energy grid.

Whether the 24 × 24 mesh is converged depends on your material: increase it
until the curve stops changing.
<!-- source: template = examples/20_pdse2_spin_hall_routes/run_siesta_projected.json (same basis, projection, word, expansion, reconstruction); eigen + exact: use-siesta.md (Basis table); orthogonal = exact Löwdin with S = 1: reference/run-file.md (Wannier90 and PAOFLOW models), runfile.py:1580–1583; shift: python/lsquant/projection/spectral.py:119–125; window pad: python/lsquant/bases.py:1410–1416; resolution: python/lsquant/reconstruct.py:35–36, pi*3.6/1024 = 0.011 eV; energy_grid: python/lsquant/apps/stage6b_siesta.py:132–143, projection/runfile.py:214 -->

Check the file, then run it:

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Values and meaning |
|---|---|
| `model.params.velocity_mode` | `"berry_connection"` (the default), the velocity of PAOFLOW's `gradient_and_momenta`, or `"bare"`. PAOFLOW has no position matrix, so `"covariant"` stops the run. |
| `model.params.oracle` | `true` adds an exact reference in a run without a projection (`exact_shc.json`). A projection run writes its exact reference in any case. |
| `projection.window` | In eV, measured from the DFT Fermi energy. |
<!-- source: default velocity runfile.py:1515; covariant refused models/orthogonal.py:418, runfile.py:1515–1516; oracle: apps/stage6b_siesta.py:322, projection/runfile.py:214–228 -->

Leave out `model.params.spin` and `model.params.fermi_energy`: both come from
the file, and either key stops the run with
`a PAOFLOW model carries its own on-site spin and energies from E_F; remove model.params.spin and model.params.fermi_energy`.
LSQuant also stops the run when H(k) is not Hermitian or the file lacks a
field.
<!-- source: python/lsquant/runfile.py:1512–1514 -->

With a projection, LSQuant writes the KPM curve to `word_0.dat` and the exact
result to `word_0_oracle.json`, which contains the exact curve both before and
after the projection; `projection.json` gives the difference between the KPM
curve and the exact one. Without a projection, it writes `shc.dat` and
`exact_shc.json`. The curves are σ^z_xy in units of e²/h, with the Pauli
matrix σ_z in the spin current; divide by 2 to compare with values for spin
ħ/2. See [Output files](../reference/output-files.md#projection-runs) and
[Results and units](use-siesta.md#results-and-units).
<!-- source: reference/output-files.md (Projection runs; Clean SIESTA runs); docs/explanation/kubo-formulas.md -->

PAOFLOW's own spin Hall module gives σ^z_xy in units of e²/ħ with the spin
operator σ_z/2. Multiply its values by 4π to compare them with LSQuant: 2π from
e²/ħ to e²/h, and 2 from σ_z/2 to σ_z. On the PdSe₂ monolayer of example 20,
LSQuant on a Wannier90 model and PAOFLOW on the same Quantum ESPRESSO bands
give spin Hall plateaus 0.1 % apart.
<!-- source: examples/20_pdse2_spin_hall_routes/README.md (PAOFLOW route, factor 4 pi; The result: 1.944 vs 1.942 e^2/h) -->

The cost grows with the number of bands of the model at each k-point, because
the shift projection keeps every band. The template run on a PAOFLOW model of
PdSe₂ with 108 bands per k-point took about 2 h 30 min on 8 threads, against
about 5 min for a 12-band Wannier90 model of the same material.
<!-- source: measured, examples/20_pdse2_spin_hall_routes README Runtime / data/wannier_covariant_spn.json, AMD EPYC 9374F -->

## 4. Recommendations

- Use pseudopotentials whose atomic wavefunctions include the empty orbitals
  of the conduction bands you need (for example a transition-metal s and p
  shell above the d shell). With only occupied or semicore shells, the model
  reproduces the DFT bands at the DFT k-points but not between them.
- `write_paoflow_model` warns when H(R) has not decayed at the edge of the
  PAOFLOW k-mesh (`edge_ratio` above 10⁻², stored in the file). Such a model
  can contain states between the DFT k-points that the DFT bands do not have;
  use a denser nscf k-mesh or a richer projection basis.
- Keep `projection.window` below the `shift` that PAOFLOW's `projectability`
  prints: above it, the bands were moved by PAOFLOW.
- In a trivial gap, the spin Hall conductivity is expected to be close to
  zero. A clearly non-zero value there comes from the spin operator (see the
  channel order above) or from spurious states.
