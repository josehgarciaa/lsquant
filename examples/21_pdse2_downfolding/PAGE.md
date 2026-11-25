# Downfolding: the same spin Hall plateau from 220 orbitals or 12

A DFT basis carries many more orbitals than transport needs. LSQuant can
reduce a SIESTA Hamiltonian to the bands near the Fermi level, and it reads
the reduced models that Wannier90 and PAOFLOW build from Quantum ESPRESSO.
For a pentagonal PdSe₂ monolayer, the 12-orbital models give the quantum spin
Hall plateau of the full 220-orbital calculation within 2 %, with 12 bands per
k point instead of 220 in the expansion.

![Spin Hall conductivity of PdSe2 from five models of different size](fig_downfolding.png)

**FIG. 1.** σᶻ_xy of pentagonal PdSe₂ against the energy from the
valence-band maximum, in e²/h with the Pauli σ_z. Black solid line: SIESTA
with all 220 spinor orbitals, exact diagonalisation. Black dashed line: the
same SIESTA Hamiltonian with the states outside ±3 eV of E_F moved to the
window edges. Thin blue line with circles: SIESTA reduced to its 12 isolated
bands around E_F. Green dash-dotted line: Quantum ESPRESSO with a 12-orbital
Wannier90 model. Pink dotted line: Quantum ESPRESSO with a 108-orbital
PAOFLOW model. All but the first are KPM curves: 24 × 24 k mesh,
M = 1024 moments, Jackson kernel, 6 meV broadening. (a) The full energy
range. (b) The quantum spin Hall gap; dark grey: the gap of Quantum ESPRESSO,
2.01–2.06 eV; light grey: the gap of SIESTA, 1.99–2.06 eV.

| Model | Bands per k point | Plateau (e²/h) | σᶻ_xy in the trivial gap (e²/h) |
|---|---|---|---|
| SIESTA, all states (exact) | 220 | 1.98 | +0.01 |
| SIESTA, window ±3 eV | 220 | 2.01 | −0.01 |
| SIESTA, reduced to 12 bands | 12 | 1.98 | −0.06 |
| QE + Wannier90 | 12 | 2.02 | −0.07 |
| QE + PAOFLOW | 108 | 2.00 | +0.02 |

<!-- Plateau: mean over the central half of each topological gap; trivial gap: E - E_VBM = 0.5 eV. Reduced SIESTA, exact on the same 24x24 mesh: plateau 1.964, trivial gap -0.061, the same as its KPM value, so the trivial-gap shift comes from the model and not from KPM. -->

In the quantum spin Hall gap all models agree within 2 %. In the trivial gap
the two 12-band models give −0.06 e²/h instead of +0.01: they leave out the
bands outside their window, which add about +0.07 e²/h to the Fermi-sea sum
there.

**Set it up and run it:** [Example 21 README](README.md).

## Check

Verify the published plateaus from the committed curves, from the repository root:

```bash
python3 examples/21_pdse2_downfolding/check.py --data-only
```

Compare completed SIESTA outputs with their matching committed curves:

```bash
python3 examples/21_pdse2_downfolding/check.py --siesta-direct out/pdse2_routes --siesta-projected out/pdse2_routes_projected --siesta-reduced12 out/pdse2_reduced12
```

Supply only the output directories available. Omitted routes are reported as
not checked. Every supplied route requires its topology JSON; the current
projection execution path does not produce that file, which is reported as a
failed input requirement. The direct reconstruction has no committed matching
reference, so its curve comparison is explicitly reported as not checked.
