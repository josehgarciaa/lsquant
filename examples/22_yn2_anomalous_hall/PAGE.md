# Quantum anomalous Hall plateau of a YN₂ monolayer, from two DFT codes

A Chern insulator has a Hall conductivity σ_xy = C e²/h with E_F in its gap,
and C is fixed by the occupied bands alone. A calculation that gets C right
has to get the Berry curvature of the bands, and so their velocity matrix
elements, right everywhere in the Brillouin zone. The ferromagnetic 1T-YN₂
monolayer is predicted to have C = 3 [1]. We compute its anomalous Hall
conductivity with LSQuant from SIESTA and from Quantum ESPRESSO through
Wannier90, and we split the Kubo–Bastin curve into its Fermi-sea and
Fermi-surface parts.

![Bands and anomalous Hall conductivity of YN2 from SIESTA and Quantum ESPRESSO](fig_yn2.png)

**FIG. 1.** 1T-YN₂, PBE with spin–orbit coupling, magnetisation along z.
(a) Bands from SIESTA (black solid line, numerical atomic orbitals) and from
Quantum ESPRESSO with a 12-orbital Wannier90 model (green dashed line, plane
waves). (b) σ_xy in e²/h. Black solid line: SIESTA, exact sum over occupied
states. Green dashed line: Quantum ESPRESSO + Wannier90, exact sum. Both on
a 300 × 300 k mesh, with a 1 meV interband regulator and a 2 meV Gaussian.
Blue dotted line: the same Wannier90 model, KPM: 300 × 300 k mesh,
M = 2048 moments, Jackson kernel, 4 meV broadening. Grey line: −3 e²/h.
(c) The KPM curve of (b) (blue dotted line) and its two parts: Fermi sea
(pink solid line) and Fermi surface (orange dash-dotted line).

| Route | σ_xy(E_F), exact (e²/h) | σ_xy(E_F), KPM (e²/h) |
|---|---|---|
| SIESTA | −3.01 | — |
| QE + Wannier90 | −3.00 | −3.13 |

<!-- Exact SIESTA 600 x 600 gives -2.999. Global gap of the SIESTA bands: 36 meV, on the M-K line. postw90 (berry_task = ahc) on the same Wannier90 model gives +3.00: the same magnitude, opposite sign convention. KPM with 120 x 120 k points gives -3.74 at E_F: the Berry curvature is concentrated near M-K and needs the dense mesh. -->

The two DFT codes give the same bands near E_F and the same plateau,
C = −3 in the sign convention of LSQuant. The KPM curve reproduces the exact
one within 4 % at E_F; the Jackson kernel of 2048 moments resolves the
plateau, which is about 25 meV wide, only just.

In the gap the Fermi-surface part vanishes and the Fermi-sea part carries the
whole plateau, as for a Chern insulator it must: the quantised Hall current
comes from all occupied states, not from states at E_F. Outside the gap, within
0.3 eV of E_F, the Fermi-surface part stays below 0.1 e²/h. [Example 13](../13_bastin_sea_surface/README.md)
defines the two parts.

[1] X. Kong, L. Li, O. Leenaerts, W. Wang, X. Liu and F. M. Peeters,
Quantum anomalous Hall effect in a stable 1T-YN₂ monolayer with a large
nontrivial bandgap and a high Chern number, *Nanoscale* **10**, 8153 (2018),
[doi:10.1039/C8NR00571K](https://doi.org/10.1039/C8NR00571K).

**Set it up and run it:** [Example 22 README](README.md).

## Check

Check the committed Hall curves and path band gaps, from the repository root:

```bash
python3 examples/22_yn2_anomalous_hall/check.py --data-only
```

Compare completed anomalous Hall outputs with the committed curves:

```bash
python3 examples/22_yn2_anomalous_hall/check.py out/yn2_ahc
```

The signed Chern plateau is negative in this convention. The gap inferred from
the committed band path is checked separately from the conductivity plateau;
it does not establish a global gap over the full Brillouin zone.
