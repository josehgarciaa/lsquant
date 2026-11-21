# Non-orthogonal basis: the Haldane model rewritten with an overlap

`haldane_rotated.py` takes the golden Haldane model (512 sites, Chern number C = +1),
rewritten in a local non-orthogonal basis, from `goldens/stage0/nonorth/haldane/complex`:
- H' = R†HR and V' = R†VR;
- the overlap S' = R†R, with cond(S') = 5.2.

The physics is the same as in the orthogonal model, so every correct method must give its
Chern number.

The script computes the VX–VY moments (M = 128, exact trace) and the Kubo-Bastin σ_xy(E),
and takes C = 2π⟨σ_xy⟩/A on the gap plateau |E| < 0.5 eV. It does this in three ways:

| Route | Recursion | Velocities | Bra side |
|---|---|---|---|
| Löwdin | S^{−1/2} H' S^{−1/2} | S^{−1/2} V' S^{−1/2} | the orthogonal plan |
| generalised | S⁻¹(aH' + bS) | S⁻¹V' | S T_m(S⁻¹H̃) S⁻¹V_x S⁻¹ r (covariant) |
| naive | S⁻¹(aH' + bS) | S⁻¹V' on both sides | no S (the legacy estimator) |

The generalised route also runs with the polynomial (Chebyshev, order 60) and the
iterative (Ischo, tolerance 1e-13) appliers of S⁻¹.

## Run

```bash
PYTHONPATH=build/lsquant_core:python python3 examples/nonorthogonal/haldane_rotated.py
```

Output on a 64-core workstation (Ubuntu 22.04, g++ 11.4), 8 threads, 3 min 54 s wall time:

```
rotated Haldane model: 512 sites, M = 128, exact trace, window (-4.5, 4.5)
orthogonal golden:          C = +1.0170
lowdin, exact S^-1/2                   C = +1.0170   max|mu - mu_golden|/max|mu_golden| = 4.0e-14    28.1 s
generalised, exact S^-1                C = +1.0170   max|mu - mu_golden|/max|mu_golden| = 1.2e-14    18.0 s
generalised, Chebyshev S^-1 (P = 60)   C = +1.0170   max|mu - mu_golden|/max|mu_golden| = 5.3e-13    69.9 s
generalised, Ischo S^-1 (tol 1e-13)    C = +1.0170   max|mu - mu_golden|/max|mu_golden| = 4.4e-14   100.1 s
naive (legacy estimator)               C = +1.1935   max|mu - mu_golden|/max|mu_golden| = 1.7e-01    16.4 s
checks: PASS
```

The script exits with status 1 if a correct route misses the orthogonal moments by more
than 1e-8, or if the naive route comes within 1 %.

## Reading the result

- **Covariance.** The Löwdin and generalised routes reproduce the orthogonal moments to
  rounding, about 1e-14, and give the orthogonal C = +1.0170. The basis change leaves the
  physics unchanged, as it must.
- **Appliers.** The Chebyshev applier of order 60 reaches 5e-13 here: S' is well
  conditioned, so 60 terms of the series of 1/s are enough. The Ischo applier reaches the
  exact result. Both cost more than the exact dense S⁻¹ at this small size; at large sizes
  the dense applier is not available (`max_dimension`, 4096 by default).
- **The naive estimator.** Omitting S on the bra computes Tr[V T_m S⁻² V T_n] instead of
  Tr[V T_m S⁻¹ V T_n S⁻¹]. The moments are off by 17 %, and the Hall plateau gives
  C = 1.19, wrong by 17 %. The legacy non-orthogonal drivers used this estimator.
  `docs/nonorthogonal.md` explains the derivation.
