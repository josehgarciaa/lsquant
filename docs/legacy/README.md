# Legacy pages

The pages in this directory describe the legacy drivers of LSQuant, built from
`src/`, `include/` and `test/`. They stay here, unchanged, while the new code is
tested. Statements such as the `RECON_CUTOFF` constant, the Lorentz kernel and
the division by Z hold only for the legacy drivers.

For the current code, start at the [documentation index](../index.md). Each
legacy page has a current counterpart:

| Legacy page | Current page |
|---|---|
| [`spectral_scales.md`](spectral_scales.md) | [The kernel polynomial method](../explanation/kpm.md) and [Kubo formulas](../explanation/kubo-formulas.md) |
| [`reconstruction.md`](reconstruction.md) | [Architecture](../explanation/architecture.md) |
| [`time_evolution.md`](time_evolution.md) | [Time evolution](../explanation/time-evolution.md) |
| [`spinrelaxation.md`](spinrelaxation.md) | [Time evolution](../explanation/time-evolution.md) and [example 06](../../examples/06_spin_relaxation_square/README.md) |
