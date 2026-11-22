# Use an OpenMX calculation

Connect an OpenMX calculation with spin–orbit coupling to LSQuant and compute
its spin Hall conductivity. LSQuant reads the Hamiltonian, the overlap and the
position matrix that OpenMX writes; the DFT calculation itself (basis,
pseudopotentials, k mesh, convergence) is yours. You should know OpenMX and
have read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

An OpenMX model is non-orthogonal, like a SIESTA model, and runs through the
same executors with the same keys. Every key is listed in the
[run-file reference](../reference/run-file.md#openmx-models).

## Requirements

- LSQuant built and importable; see [Install LSQuant](install.md).
- OpenMX 3.9 or later; LSQuant reads the scfout format of OpenMX 3.9.
<!-- tested with OpenMX 3.9.9 and the DFT_DATA19 basis and pseudopotentials -->

## 1. Turn on in OpenMX

In the input file:

```text
scf.SpinPolarization     NC     # non-collinear spinors
scf.SpinOrbit.Coupling   on     # with j-dependent pseudopotentials (*_PBE19, ...)
HS.fileout               on     # writes <System.Name>.scfout
```

## 2. Transfer the model

Nothing to run. LSQuant reads `<System.Name>.scfout`: H(R), S(R), the
chemical potential and the position matrix of the basis orbitals, which the
covariant velocity uses. Energies are measured from the chemical potential,
so E_F = 0 in the run file.

## 3. Run LSQuant

```json
{
  "model": {
    "kind": "openmx",
    "file": "/path/to/mysystem.scfout",
    "params": {"velocity_mode": "covariant", "oracle": true}
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "metric": {"policy": "nonorthogonal", "route": "lowdin", "applier": "exact"},
  "projection": {"kind": "shift", "window": [-3, 3]},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 6.0,
                     "energy_grid": [-0.5, 2.0, 1251]},
  "output": {"dir": "out/openmx_shc"}
}
```

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Value |
|---|---|
| `velocity_mode` | `"covariant"` (uses the position matrix of the scfout file), `"berry_connection"` or `"bare"`, as for SIESTA. |
| `metric.policy` | `"nonorthogonal"`, as for SIESTA. |
| `oracle` | `true` also writes the exact diagonalisation on the same k mesh. |

With a projection the KPM curve is `word_0.dat` and the exact one
`word_0_oracle.json`; without it, `shc.dat` and `exact_shc.json`.

LSQuant stops on a scfout file of another format version
(`only version 3 (OpenMX 3.9) is supported`) and on a calculation without
non-collinear spin–orbit coupling. OpenMX integrates part of H, S and the
position matrix on its real-space grid, so their Hermitian symmetry holds
only to round-off; LSQuant folds them Hermitian and records the removed
part in the model metadata (`hermitian_fold_relative`, copied to `basis.json`
by the clean executor).

## 4. Recommendations

- Use a basis with enough radial functions for the conduction bands you need
  (for example the precise `s3p2d2` sets of DFT_DATA19 for 4d metals and
  chalcogens).
- For the spin Hall conductivity use `velocity_mode = "covariant"`.
- For an insulator, the chemical potential lies somewhere in the gap, at a
  position set by the electronic temperature. To compare with another code,
  measure energies from the valence-band maximum.
