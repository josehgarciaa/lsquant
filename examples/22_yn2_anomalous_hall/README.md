# Example 22: set up and run the YN₂ anomalous Hall calculation

Reproduce the anomalous Hall conductivity of the 1T-YN₂ monolayer from
SIESTA and from Quantum ESPRESSO + Wannier90, and its split into Fermi-sea
and Fermi-surface parts. The [example page](PAGE.md) shows the result.

The figure is made from the files in [`data/`](data/), so you can plot it
without running anything (step 6). Steps 1 to 5 regenerate those files.

## Requirements

- LSQuant built and importable; see [Install LSQuant](../../docs/how-to/install.md).
- SIESTA with spin–orbit coupling.
- Quantum ESPRESSO with `pw2wannier90.x`, and Wannier90.
- For the KPM run (step 4): about 150 GB of memory and 2 h on 48 cores; see
  [Runtime](#runtime).
<!-- tested with SIESTA 5.4.2, Quantum ESPRESSO 7.4, Wannier90 3.1.0 -->

## The inputs

| | SIESTA | Quantum ESPRESSO |
|---|---|---|
| Files | [`siesta/yn2.fdf`](siesta/yn2.fdf) | [`qe/`](qe/) |
| Functional | PBE + SOC | PBE + SOC |
| Pseudopotentials | PseudoDojo NC FR 0.4, standard (psml) | PseudoDojo NC FR 0.4, standard (upf) |
| Basis | DZP; Y with 4s and 4p semicore | plane waves, 80 Ry |
| Self-consistent k mesh | 36 × 36 | 24 × 24 |

Both use the same cell ($a = 3.776$ Å, 24 Å along $z$) and atomic positions:
Y at the centre, one N 0.877 Å above and one 0.877 Å below, Y–N bond
2.350 Å. The initial moments on N point along $z$; the self-consistent
moment is 3 μB per cell.

## Step 1: SIESTA

Put `Y.psml` and `N.psml` in `siesta/` and run, in `siesta/`:

```bash
mpirun -np 16 siesta < yn2.fdf > yn2.out
```

`SaveHS true` writes `yn2.TSHS`; SIESTA also writes the `*.ion.xml` files of
the basis, which LSQuant reads for the covariant velocity.

## Step 2: Quantum ESPRESSO and Wannier90

Put `Y.upf` and `N.upf` in `pseudo/`, next to `qe/`. In `qe/`:

```bash
mpirun -np 32 pw.x -nk 8 -in scf.in  > scf.out
mpirun -np 32 pw.x -nk 8 -in nscf.in > nscf.out
wannier90.x -pp y
mpirun -np 4 pw2wannier90.x -in pw2wan.in > pw2wan.out
wannier90.x y
```

The model has 12 spinor Wannier functions from the N p orbitals; the 12
lowest bands (Y 4s, 4p and N 2s) are excluded. The energy windows in
[`qe/y.win`](qe/y.win) are set from the scf Fermi energy, $E_F = -3.504$ eV
with these pseudopotentials: the disentanglement window is
$E_F - 8$ eV to $E_F + 5$ eV and the frozen window ends $0.4$ eV above
$E_F$. If your scf run gives a different $E_F$, shift the windows by the
difference. `write_rmn = true` writes `y_r.dat`, the position matrix that
LSQuant needs for the covariant velocity.

Copy the Wannier90 files to `w90/`, the directory named in
[`run_wannier_ahc.json`](run_wannier_ahc.json):

```bash
mkdir -p ../w90 && cp y.win y_hr.dat y_r.dat y_centres.xyz ../w90/
```

## Step 3: exact anomalous Hall conductivity

[`exact_hall.py`](exact_hall.py) sums the Berry curvature of the occupied
states on a 300 × 300 mesh, with the covariant velocity of LSQuant for both
models. From the repository root:

```bash
python examples/22_yn2_anomalous_hall/exact_hall.py siesta \
    examples/22_yn2_anomalous_hall/siesta siesta_exact --workers 16
python examples/22_yn2_anomalous_hall/exact_hall.py wannier \
    examples/22_yn2_anomalous_hall/w90 wannier_exact --fermi-energy -3.504 --workers 16
```

Each writes `<name>.dat` (energy from $E_F$, σ_xy in e²/h) and prints
σ_xy(E_F). `--fermi-energy` is the scf Fermi energy of step 2.

## Step 4: KPM and the Fermi-sea and Fermi-surface parts

[`run_wannier_ahc.json`](run_wannier_ahc.json) asks for the anomalous Hall
word and three reconstructions:

```json
"word": "Tr[VX delta(E-H) VY dG-(E)/dE]",
"reconstruction": {"kernel": "jackson", "formulas": ["ahc", "sea", "surf"], "broadening_meV": 4.0}
```

`ahc` is the full Kubo–Bastin curve in e²/h; `sea` and `surf` are its two
parts with the same normalisation. From the repository root:

```bash
lsquant run examples/22_yn2_anomalous_hall/run_wannier_ahc.json
```

The run writes `ahc.dat`, `ahc_sea.dat` and `ahc_surf.dat` in `out/yn2_ahc/`,
with energies from $E_F$ (the `fermi_energy` of the run file).

## Step 5: bands

```bash
python examples/22_yn2_anomalous_hall/bands.py examples/22_yn2_anomalous_hall/siesta \
    examples/22_yn2_anomalous_hall/w90 bands.npz --fermi-energy -3.504
```

## Step 6: plot

Copy the files of steps 3 to 5 to `data/` with these names, or keep the
committed ones:

| File in `data/` | From |
|---|---|
| `siesta_exact.dat`, `wannier_exact.dat` | step 3 |
| `wannier_kpm_ahc.dat`, `wannier_kpm_sea.dat`, `wannier_kpm_surf.dat` | step 4: `ahc.dat`, `ahc_sea.dat`, `ahc_surf.dat` |
| `bands.npz` | step 5 |

```bash
python examples/22_yn2_anomalous_hall/plot.py
```

writes `fig_yn2.png`.

## Runtime

The KPM run of step 4 took 2 h on 48 cores and 142 GB of memory: it keeps
300 × 300 × 12 = 1.08 million states. A 120 × 120 mesh takes 11 min and
23 GB, but its curve overshoots the plateau (−3.7 e²/h at $E_F$): the Berry
curvature is concentrated near the M–K line and needs the dense mesh. The
exact sums of step 3 take minutes for the Wannier90 model and under an hour
for SIESTA on 12 processes.

## Sign convention

LSQuant gives σ_xy = −3 e²/h here. Wannier90's `postw90.x` (`berry_task =
ahc`) gives +3 e²/h for the same Wannier90 model, with the same magnitude at
every energy: the two codes use opposite sign conventions for σ_xy.
