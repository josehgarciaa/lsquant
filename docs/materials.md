\page materials Materials: SIESTA Hamiltonians, spin Hall conductivity, disorder and projection

This part brings DFT Hamiltonians into the engine. It combines the development stages
6b (SIESTA reader and spin Hall conductivity), 6c (disorder in a non-orthogonal basis),
6d (spectral projection) and 8 (orthogonal downfolding). It builds on the overlap metric of
\ref nonorthogonal, the bases of \ref bases and the run files of \ref run_files.

## SIESTA models

```json
"model": {"kind": "siesta", "dir": "examples/07_spin_hall_bi_bilayer/siesta", "label": "bi",
          "params": {"velocity_mode": "covariant"}}
```

`lsquant.models.siesta.load_siesta_model(dir, label)` reads `<label>.fdf` with the
`.TSHS` or `.HSX` file through sisl. It needs spin-orbit coupling, and every orbital
becomes a spinor pair, so W = 2 N_orb. The real-space blocks H(R) and S(R) are
reconstructed from sisl's storage and checked against sisl's own H(k) and S(k) in both
gauges at a generic k-point; an unknown storage layout is refused. TSHS files store rounded
values: the reader measures the non-Hermitian part (up to 1e-7 relative), records it, and
folds it out in the Bloch sums.

| `params` key | Values | Effect |
|---|---|---|
| `velocity_mode` | `bare`, `berry_connection` (default), `covariant` | the velocity matrix (below) |
| `oracle` | `true` / `false` | also evaluate the exact-diagonalisation spin Hall conductivity on the same k mesh |
| `oracle_degeneracy_tolerance` | eV | grouping of degenerate states in the oracle |
| `bloch_memory_mb` | MiB | memory budget of the Bloch blocks |
| `nonorthogonal_disorder` | `true` / `false` | Anderson disorder in the non-orthogonal basis (below) |

**Velocities.** The Hermitian velocity in a non-orthogonal Bloch basis is

V(k) = ∂H − ½(∂S S⁻¹H + H S⁻¹∂S) [+ i(H S⁻¹A − A S⁻¹H)].

- `bare` uses cell derivatives and no connection A;
- `berry_connection` uses atom-resolved derivatives, i(R + τ_j − τ_i), and no A;
- `covariant` adds the midpoint dipole connection A, integrated from the radial functions
  of the `.ion.xml` files (`dipole_quadrature`). The result is cached in a
  `.lsquant-cache` directory outside the input directory and is recomputed when the inputs
  change.

`lsquant.metric_forms` derives these forms symbolically (SymPy). It proves that the
Löwdin and generalised routes give the same trace for every word shape, and it states the
Hermiticity of each operator form.

## Spin Hall conductivity

The spin Hall conductivity σ^z_xy is the Bastin word with a spin current,

```json
"word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
"reconstruction": {"formulas": ["shc"], "broadening_meV": 10.0, "fermi_energies": [0.0]}
```

where J^z_x = {V_x, S_z}/2 in Pauli units. The formula `shc` returns
σ = 2π·KuboBastin/(N_k·A_cell) in units of e²/h (Pauli σ_z). For SIESTA models the run takes
N_k and A_cell from the supercell and the lattice. Three metric routes are available:
`lowdin`, `generalised` and `anticommutator`. The last uses the symmetrised operators
(S⁻¹O + OS⁻¹)/2.

The oracle (`lsquant.oracle.spin_hall`) diagonalises H(k) c = E S(k) c on the same mesh. It
returns the spin Berry curvature, σ^z_xy(E_F) and the spin Chern number, so that every
KPM result has an exact reference.

**Known limitation.** With `model.kind` = `csr` or `wannier`, `python -m lsquant run` writes
the moments of an `shc` word and then exits with status 1: the three run paths call the
formula with `num_cells` and `cell_area` unset. To get the curve, give `num_cells` and
`cell_area` in `reconstruction`, run without `shc`, and call `python -m lsquant reconstruct
<file.chebmom2D> --formula shc --broadening-meV <b>`. SIESTA models are not affected.

## Disorder in a non-orthogonal basis

With `nonorthogonal_disorder = true` and `basis.disorder = {"W": ..., "seed": ...}`, the
on-site Anderson term becomes, in the non-orthogonal basis,

V_ij = (w_i + w_j) S_ij / 2

(`lsquant.models.disorder.anderson_potential`). For a uniform w this equals w·S, the
identity operator in this basis. The values w_i are drawn uniformly in [−W/2, W/2], one per
spinor coefficient, in the layout of the real-space basis. Up and down components of an
orbital therefore get independent values: the disorder carries a random spin-dependent
part, a Zeeman-like field along z. `lsquant.disorder_bases` provides the disordered Bloch
basis (FFT sandwich) and the real-space basis with the overlap, and both give the same
moments for the same seed.

## Projection and downfolding

The optional `projection` section (`kind` = `none`, `hard`, `shift` or `orthogonal`, with
`window` in eV and `radius` in Å) changes a SIESTA model so that a word sees only the
states of an energy window. The coordinates and the overlap stay the same: the Bloch blocks
are diagonalised once per k, and the eigenvalues outside the window are moved:
- `hard` sets them to 0;
- `shift` moves each to the nearest window edge;
- `orthogonal` reduces the model to an orthogonal tight-binding model of the window
  manifold with the downfolding tools of `tools/downfold_origin/`, driven by
  `lsquant.projection.orthogonal`. It needs a window of constant rank around an isolated
  manifold; `radius` truncates the reduced hoppings.

`lsquant.projection.diagnostics` compares the projected traces with the full ones, and
`runs/projection/` has two run files for the Bi and PdSe₂ examples.

## Examples

`examples/07_spin_hall_bi_bilayer/` and `examples/08_spin_hall_pdse2_monolayer/` hold the
SIESTA inputs and checksums (`siesta/`) and quick and production run files for the DOS and
the spin Hall conductivity. Each README gives the commands, the run times and the results.
For the Bi bilayer, the quick run (24 × 24 mesh, 20 meV) gives an exact-diagonalisation
σ^z_xy(E_F = 0) = 1.144 e²/h and the spin Chern number 2 (Z₂ = 1) across a gap of 0.304 eV.
The KPM curve of the quick run (M = 256) resolves about 0.6 eV on the SIESTA window, more
than the gap, so it stays far from the plateau; the converged KPM value of the production
setting is still pending.

## Tests

| Test | What it checks | Run |
|---|---|---|
| `tests/unit/test_python_stage6b.py` | the symbolic route equivalence for every word shape; Hermiticity of the velocity forms (S-Hermiticity of S⁻¹V); the Bi SIESTA reader: spinor blocks, S(0) Hermitian to the TSHS rounding and positive definite, H(k) and S(k) Hermitian, real generalised spectrum; `anderson_potential` (formula, uniform case w·S, refusals) | `PYTHONPATH=build/lsquant_core:python python3 -m pytest tests/unit` |
| `tests/dev/stage_6b_check.py` | the SIESTA reader and velocity ladder; the spin Hall conductivity of a Kane-Mele model against the oracle, by every route and basis, with the noise of random-phase traces | `PYTHONPATH=python:build:build/lsquant_core python3 tests/dev/stage_6b_check.py` |
| `tests/dev/stage_6c_check.py` | disorder in the non-orthogonal Bloch and real bases: identity of the moments, the uniform-w limit | same |
| `tests/dev/stage_6d_check.py` | spectral projection and downfolding: projected traces against the full ones | same |

The stage-7 gate covers the regression of every golden (G1–G10) after these changes.
