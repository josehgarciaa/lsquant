# Use a GPAW calculation

Connect a GPAW calculation in the LCAO mode to LSQuant and compute its
density of states or its spin Hall conductivity. GPAW adds spin–orbit
coupling after the self-consistent cycle; LSQuant reads the scalar-relativistic
H(R) and S(R) of the atomic orbitals and the same spin–orbit term, written as
a matrix in those orbitals. You should know GPAW and have read
[Running LSQuant](run-lsquant.md) and [Bring a model to LSQuant](bring-a-model.md).

Every key is listed in the
[run-file reference](../reference/run-file.md#machine-learning-hamiltonians-in-the-metric-executors).

## Requirements

- GPAW with its Python (the old calculator, `gpaw.calculator.GPAW`, is used
  for the export). LSQuant does not need to be installed there.
- LSQuant built and importable where you run it; see [Install LSQuant](install.md).
<!-- tested with GPAW 26.7.0 and the gpaw-setups-24.11.0 PAW setups -->

## 1. Turn on in GPAW

| Setting | Value | Why |
|---|---|---|
| `mode` | `'lcao'` | LSQuant needs the atomic-orbital matrices. |
| `symmetry` | `'off'` | The export Fourier transforms the k matrices; the mesh must cover the whole Brillouin zone. |
| `kpts` | a mesh larger than twice the range of the orbitals, e.g. `{'size': (12, 12, 1), 'gamma': True}` for a 6 Å cell | The transform is exact when the mesh is larger than that range. |
| `parallel` | `{'domain': 1}` | The export needs every atom on every rank; k points may be distributed. |
| spin | spin paired (no magnetic moments) | Spin–orbit coupling is added by the export. |

Write the restart file after the self-consistent cycle:

```python
atoms.get_potential_energy()
atoms.calc.write("pdse2.gpw")
```

## 2. Transfer the model

In the GPAW environment, run the export script that comes with LSQuant:

```bash
python /path/to/lsquant/python/lsquant/models/gpaw.py pdse2.gpw pdse2_lsquant.npz
```

It also runs under `mpirun -np N gpaw python ...`. The script writes H(k),
S(k) and the spin–orbit matrix on the k mesh of the calculation, the
structure and the Fermi level of GPAW. `--no-soc` leaves out the spin–orbit
term. In a running GPAW script the same step is
`export_gpaw_lcao(calc, "pdse2_lsquant.npz")` from `lsquant.models.gpaw`.

## 3. Run LSQuant

```json
{
  "model": {
    "kind": "gpaw",
    "file": "/path/to/pdse2_lsquant.npz",
    "params": {"velocity_mode": "covariant", "dipole_cache": "/path/to/pdse2_dipoles.npz"}
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "projection": {"kind": "shift", "window": [-3, 3]},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 10.0,
                     "energy_grid": [-1.0, 3.0, 801]},
  "output": {"dir": "out/gpaw_shc"}
}
```

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Value |
|---|---|
| `velocity_mode` | `"covariant"` (recommended), `"berry_connection"` or `"bare"`. For `"covariant"` LSQuant builds the position matrix of the orbitals from the radial functions in the export (a few minutes; see below). |
| `dipole_cache` | Optional `.npz` path where LSQuant stores that position matrix and reuses it for the same export. |
| `quadrature_order` | Gauss–Legendre order of that integration, default 48. LSQuant checks the overlap it integrates against GPAW's S(R) and stops if they differ by more than 10⁻⁴. |
| `soc` | `true` (default): spinful model, H ⊗ 1 + spin–orbit. `false`: spinless model. |
| `fermi_energy` | In eV, removed as H − E_F S. Default: the Fermi level of GPAW. |

## What LSQuant builds and what you should know

| Quantity | From GPAW | Built by LSQuant |
|---|---|---|
| H(R), S(R) | k matrices of the LCAO basis | inverse Bloch transform, each element at its nearest lattice image |
| Spin–orbit coupling | the on-site PAW term of `soc_eigenstates` | its matrix in the orbitals, `P† (dV/dr / r) L·σ P` |
| Spin | — | Pauli matrices on interleaved spinors |
| Position matrix | radial functions of the basis and of the PAW projectors, PAW corrections | `<M,0|r|N,R>` with the PAW terms, by numerical integration; the covariant velocity |
| Velocity | — | covariant (or the derivative of the Bloch sums, `berry_connection`) |

- **Spin–orbit coupling is not self-consistent.** This is how GPAW treats it.
  With every LCAO band, the eigenvalues of the LSQuant model equal those of
  `soc_eigenstates`.
- **Basis.** Check the bands of the LCAO basis (`dzp` by default) against a
  plane-wave or a larger-basis calculation before you trust the transport.
- **Velocity.** Use `covariant`. `berry_connection` leaves out the dipoles
  inside each atom; for the spin Hall plateau of PdSe₂ that changes the
  result by about 1.5 %.
