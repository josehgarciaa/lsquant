# Example 20: set up and run the PdSe₂ spin Hall routes

This page tells you how to reproduce the spin Hall conductivity of pentagonal
PdSe₂ from SIESTA, OpenMX, Quantum ESPRESSO + Wannier90 and Quantum ESPRESSO +
PAOFLOW. The result and what it shows are on the
[example page](PAGE.md).

The figure is made from the curves in [`data/`](data/), so you can plot it
without running anything (step 6). Steps 1 to 5 regenerate those curves.

## Requirements

- LSQuant built and importable; see [Install LSQuant](../../docs/how-to/install.md).
- For the Quantum ESPRESSO routes: Quantum ESPRESSO with `pw2wannier90.x`,
  Wannier90, and PAOFLOW in the Python environment of LSQuant.
- For the OpenMX route: OpenMX 3.9 with the DFT_DATA19 database.
- About 1 h on 32 cores for the DFT steps; see [Runtime](#runtime).
<!-- tested with Quantum ESPRESSO 7.4, Wannier90 3.1.0, PAOFLOW 3.0.0, OpenMX 3.9.9 -->

The four connectors are described in
[Bring a model to LSQuant](../../docs/how-to/bring-a-model.md).

## The inputs

| | SIESTA | OpenMX | Quantum ESPRESSO |
|---|---|---|---|
| Files | `siesta/` of [example 8](../08_spin_hall_pdse2_monolayer/README.md) | [`openmx/pdse2.dat`](openmx/pdse2.dat) | [`qe/`](qe/), [`paoflow/`](paoflow/) |
| Functional | PBE + SOC | PBE + SOC | PBE + SOC |
| Pseudopotentials | PseudoDojo NC FR (psml) | DFT_DATA19 `*_PBE19` | pslibrary 1.0.0 PAW: `Pd.rel-pbe-n-kjpaw_psl.1.0.0.UPF`, `Se.rel-pbe-dn-kjpaw_psl.1.0.0.UPF` |
| Basis | DZP | Pd7.0-s3p2d2, Se7.0-s3p2d2 | plane waves, 80 / 800 Ry |
| Self-consistent k mesh | 12 × 12 | 12 × 12 | 4 × 4 |

All three use the same cell ($a = 5.750$ Å, $b = 5.920$ Å, 16.45 Å along $z$)
and atomic positions.
<!-- Geometry of wannier2sparse example 05. In the QE inputs the slab sits at
z = 0: at z = c/2 the Wannier functions cross the cell boundary along z and
Wannier90 returns spreads of about 100 A^2. Wannier model: 12 spinor functions
(Pd dx2-y2, Se pz), nscf 12x12 with 96 bands, outer window -3.5..3.0 eV, frozen
window -2.8..1.0 eV (absolute; E_F = -1.316 eV), 1000 + 400 iterations, spreads
12.6 A^2 (Pd) and 10.9 A^2 (Se, centred near the Se-Se dimers), Kramers pairs
degenerate to 1.4 meV. The 4x4 nscf mesh and 10 + 10 iterations of wannier2sparse
example 05 give no global gap near 1.2 eV and Kramers splittings up to 24 meV. -->

## Step 1: Quantum ESPRESSO and Wannier90

Put the two pseudopotentials in `pseudo/`, next to `qe/`, as `Pd.UPF` and
`Se.UPF`. In `qe/`:

```bash
mpirun -np 32 pw.x -nk 8 -in scf.in  > scf.out
mpirun -np 32 pw.x -nk 8 -in nscf.in > nscf.out
wannier90.x -pp pdse2_proj
mpirun -np 16 pw2wannier90.x -in pw2wan_spn.in > pw2wan_spn.out
wannier90.x pdse2_proj
```

Then copy the Wannier90 files to `w90/`, the directory named in
`run_wannier_kpm.json`:

```bash
mkdir -p ../w90 && cp pdse2_proj.win pdse2_proj_* pdse2_proj.spn pdse2_proj.eig ../w90/
```

## Step 2: PAOFLOW

Copy the directory `tmp/` of the self-consistent run to `paoflow/tmp/`. In
`paoflow/`:

```bash
mpirun -np 32 pw.x -nk 8 -in nscf_pao.in > nscf_pao.out
projwfc.x -in projwfc.in > projwfc.out
python run_pao.py
```

`run_pao.py` writes PAOFLOW's own spin Hall curve and `pdse2_paoflow.npz`, the
model that LSQuant reads.

## Step 3: OpenMX

Set `DATA.PATH` in `openmx/pdse2.dat` to your DFT_DATA19 directory. In
`openmx/`:

```bash
mpirun -np 32 openmx pdse2.dat > pdse2.std
```

`HS.fileout on` makes OpenMX write `pdse2.scfout`.

## Step 4: LSQuant

From the repository root:

```bash
python -m lsquant run examples/20_pdse2_spin_hall_routes/run_siesta_direct.json
python -m lsquant run examples/20_pdse2_spin_hall_routes/run_siesta_projected.json
python -m lsquant run examples/20_pdse2_spin_hall_routes/run_wannier_kpm.json
python -m lsquant run examples/20_pdse2_spin_hall_routes/run_paoflow_kpm.json
python -m lsquant run examples/20_pdse2_spin_hall_routes/run_openmx_kpm.json
```

Each run writes the KPM curve (`word_0.dat`) and the exact sum on the same
k mesh (`word_0_oracle.json`) to its `output.dir`.

| Run file | Model | Velocity | k mesh | Moments |
|---|---|---|---|---|
| `run_siesta_direct.json` | SIESTA, all states | covariant | 24 × 24 | exact sum only |
| `run_siesta_projected.json` | SIESTA, states in $[-3, 3]$ eV | covariant | 24 × 24 | 1024 |
| `run_wannier_kpm.json` | Wannier90, spin projected from DFT | covariant | 24 × 24 | 1024 |
| `run_paoflow_kpm.json` | PAOFLOW model | berry_connection | 24 × 24 | 1024 |
| `run_openmx_kpm.json` | OpenMX | covariant | 24 × 24 | 1024 |

## Step 5: the exact Wannier90 and OpenMX curves

The exact curves in the figure use denser meshes than the run files:

```bash
python examples/20_pdse2_spin_hall_routes/wannier_route.py \
    examples/20_pdse2_spin_hall_routes/w90 pdse2_proj -1.3162 wannier_covariant_spn
python examples/20_pdse2_spin_hall_routes/openmx_exact.py \
    examples/20_pdse2_spin_hall_routes/openmx/pdse2.scfout openmx_exact
```

`wannier_route.py` uses a 300 × 300 mesh; `openmx_exact.py` a 60 × 60 mesh.

## Step 6: the figure

```bash
python examples/20_pdse2_spin_hall_routes/plot.py
```

It writes `fig_routes.png` from the curves in `data/`.

## Runtime

On an AMD EPYC 9374F:

| Step | Wall time | Peak memory |
|---|---|---|
| Quantum ESPRESSO, scf and 12 × 12 nscf (32 MPI processes) | 3 min + 12 min | — |
| OpenMX scf (32 MPI processes) | 66 s | — |
| PAOFLOW | 8 min (CPU) | — |
| `wannier_route.py`, 300 × 300 | 73 s (16 threads) | 2.8 GB |
| `run_siesta_projected.json` | 3 h 10 min (16 threads) | 23.3 GB |
| `run_wannier_kpm.json` | 5 min (8 threads) | — |
| `run_paoflow_kpm.json` | 2 h 30 min (8 threads) | — |
| `run_openmx_kpm.json` | more than 4 h (16 threads) | — |

The KPM runs of PAOFLOW and OpenMX are slow because the `shift` projection keeps every band of the model at each k point (108 and 228); `openmx_exact.py` gives the OpenMX curve of the figure in a fraction of that time.

## Files

| File | Content |
|---|---|
| `qe/`, `paoflow/`, `openmx/` | DFT and PAOFLOW inputs |
| `run_*.json` | LSQuant run files |
| `wannier_route.py`, `openmx_exact.py` | exact curves on dense meshes |
| `data/*.dat` | $\sigma^z_{xy}(E - E_F)$ of each route, in $e^2/h$ (Pauli $\sigma_z$) |
| `plot.py` | the figure |

## References

- A. A. Mostofi *et al.*, Comput. Phys. Commun. **185**, 2309 (2014).
- M. Buongiorno Nardelli *et al.*, Comput. Mater. Sci. **143**, 462 (2018).
- T. Ozaki, Phys. Rev. B **67**, 155108 (2003).
- J. H. García, L. Covaci and T. G. Rappoport, Phys. Rev. Lett. **114**, 116602 (2015).
