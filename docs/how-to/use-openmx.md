# Use an OpenMX calculation

Start from your own OpenMX calculation with spin–orbit coupling and end with
its spin Hall conductivity. LSQuant reads the Hamiltonian, the overlap and the
position matrix that OpenMX writes. Choose and converge the basis, the
pseudopotentials and the k-mesh in OpenMX, as usual. You should know OpenMX
and have read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

> **Known limitation: non-collinear spin–orbit calculations only.** LSQuant
> reads only an OpenMX calculation with non-collinear spin and spin–orbit
> coupling. For a collinear or spin-unpolarised calculation, including one
> meant for the density of states alone, the run stops with
> `load_openmx_model requires a non-collinear spin-orbit calculation (scf.SpinPolarization NC, scf.SpinOrbit.Coupling on)`.
<!-- source: python/lsquant/models/openmx.py:109–111 -->

An OpenMX model is non-orthogonal, like a SIESTA model, and LSQuant runs it
with the same settings and the same keys. It computes the density of states
and the spin Hall conductivity of the clean crystal on a k-mesh; a real-space
supercell or disorder stops the run. Every key is listed in the
[run-file reference](../reference/run-file.md#openmx-models).
<!-- source: clean only: python/lsquant/runfile.py:1576–1579; words: reference/run-file.md (SIESTA executors, Words row) -->

## Requirements

You need:

- LSQuant built and importable; see [Install LSQuant](install.md).
- OpenMX 3.9. LSQuant reads version 3 of the scfout format, which OpenMX 3.9
  writes, and a file of another version stops the run with
  `only version 3 (OpenMX 3.9) is supported`.
<!-- tested with OpenMX 3.9.9 and the DFT_DATA19 basis and pseudopotentials -->
<!-- source: python/lsquant/models/openmx.py:59–61 (scfout version must be 3). The earlier "3.9 or later" was not verified: whether later OpenMX versions write scfout version 3 is open. -->

## 1. Turn on the outputs in OpenMX

Set these keys in the OpenMX input file:

```text
scf.SpinPolarization     NC     # non-collinear spinors
scf.SpinOrbit.Coupling   on     # with j-dependent pseudopotentials (*_PBE19, ...)
HS.fileout               on     # writes <System.Name>.scfout
```

## 2. Transfer the model

There is nothing to run. LSQuant reads `<System.Name>.scfout` directly: H(R),
S(R), the chemical potential, and the position matrix of the basis orbitals
needed by `model.params.velocity_mode = "covariant"`.

Give every energy in the run file relative to the chemical potential of
OpenMX, E_F = 0; LSQuant measures the stored energies from it. The model takes
no `model.params.fermi_energy` and no `model.params.spin`: either key stops
the run with
`an OpenMX model carries its spinors and stores energies from the chemical potential; remove model.params.spin and model.params.fermi_energy`.
<!-- source: python/lsquant/models/openmx.py:9–11, 151 (H − mu S); runfile.py:1509–1511 -->

OpenMX integrates part of H, S and the position matrix on its real-space grid,
so these matrices are Hermitian only to round-off. LSQuant replaces each by
its Hermitian part and records the relative size of the removed part, about
10⁻¹⁰ for H and S and about 10⁻⁴ for the position matrix. A run without a
projection writes these values to `basis.json` as `hermitian_fold_relative`.
<!-- source: python/lsquant/models/openmx.py:152–171 (fold and comments on typical sizes); reference/run-file.md (OpenMX models: model metadata, copied to basis.json by the clean executor) -->

## 3. Run LSQuant

This run file computes the spin Hall conductivity of the clean crystal on a
24 × 24 k-mesh, with the spectrum shifted inside [−3, 3] eV:

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
- `metric` orthogonalises the atomic-orbital basis with S^(−1/2) (Löwdin),
  computed by a dense factorisation.
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
until the curve stops changing, as for the k-mesh in OpenMX.
<!-- source: template = examples/20_pdse2_spin_hall_routes/run_siesta_projected.json (same basis, metric, projection, word, expansion, reconstruction); eigen + exact: use-siesta.md (Basis table), python/lsquant/projection/runfile.py:165; Γ-centred mesh: use-siesta.md (Adapt the template, 2); metric: reference/run-file.md (metric table); shift: python/lsquant/projection/spectral.py:119–125 (clip E to [lo, hi]); window pad: python/lsquant/bases.py:1410–1416 (pad 10 % of width each side); resolution: python/lsquant/reconstruct.py:35–36 (M_eff = min(ceil(pi/eta), M)), pi*3.6/1024 = 0.011 eV, pi*3.6/0.006 = 1885; energy_grid: python/lsquant/apps/stage6b_siesta.py:132–143, projection/runfile.py:214 -->

Check the file, then run it:

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Values and meaning |
|---|---|
| `model.params.velocity_mode` | `"covariant"` (uses the position matrix of the scfout file), `"berry_connection"` (the default) or `"bare"`, as for SIESTA. |
| `metric.policy` | `"nonorthogonal"`, as for SIESTA. |
| `model.params.oracle` | `true` adds an exact reference in a run without a projection: LSQuant also diagonalises H(k) on the same k-mesh and computes the spin Hall conductivity exactly, for comparison with the KPM curve. A projection run writes its exact reference in any case. |
<!-- source: default velocity: python/lsquant/runfile.py:1515, reference/run-file.md (OpenMX models); oracle: apps/stage6b_siesta.py:322 (clean run reads params.oracle), projection/runfile.py:214–228 (projection run always computes the oracle, params.oracle not read); reference/run-file.md (SIESTA executors: oracle "not used" in projection) -->

With a projection, LSQuant writes the KPM curve to `word_0.dat` and the exact
result to `word_0_oracle.json`, which contains the exact curve both before and
after the projection; `projection.json` gives the difference between the KPM
curve and the exact one. Without a projection, it writes the curve to
`shc.dat` and the exact result to `exact_shc.json`. The curves are σ^z_xy in
units of e²/h, with the Pauli matrix σ_z in the spin current; divide by 2 to
compare with values for spin ħ/2. See
[Output files](../reference/output-files.md#projection-runs) and
[Results and units](use-siesta.md#results-and-units).
<!-- source: reference/output-files.md (Projection runs; Clean SIESTA runs); docs/explanation/kubo-formulas.md (factor 1/2) -->

The KPM curve uses the Jackson kernel and the exact curve a Gaussian, so they
differ slightly even when both are right. In example 20, the KPM plateau of
the SIESTA model lies 2 % above the exact result for the same projected
states.
<!-- source: projection/runfile.py:257 ("KPM/Jackson versus Gaussian oracle; kernels differ"); examples/20_pdse2_spin_hall_routes/README.md (The result) -->

The cost grows with the number of bands of the model at each k-point, and an
OpenMX model has many: the test model has 228 bands per k-point. The exact
reference on the 24 × 24 mesh takes minutes; the KPM part takes much longer.
The same settings on the SIESTA model of PdSe₂, with 220 spinor
orbitals per cell, took 3 h 10 min and 23.3 GB on 16 threads.
<!-- source: measured, examples/20_pdse2_spin_hall_routes README Runtime / data/wannier_covariant_spn.json, AMD EPYC 9374F (OpenMX, 228 bands per k: not finished after 4 h on 16 threads; oracle on 24x24 takes minutes); SIESTA figure: examples/20_pdse2_spin_hall_routes/README.md (Runtime) -->

## 4. Recommendations

- Use a basis with enough radial functions for the conduction bands you need
  (for example the precise `s3p2d2` sets of DFT_DATA19 for 4d metals and
  chalcogens).
- For the spin Hall conductivity, use
  `model.params.velocity_mode = "covariant"`.
- For an insulator, the chemical potential lies somewhere in the gap, at a
  position set by the electronic temperature. To compare with another code,
  measure energies from the valence-band maximum: subtract it from the energy
  column of the output.
