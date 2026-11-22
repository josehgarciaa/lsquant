# Downfolding: the same spin Hall plateau from 220 orbitals or 12

A DFT basis carries many more orbitals than transport needs. LSQuant can
reduce a SIESTA Hamiltonian to the bands near the Fermi level, and it reads
the reduced models that Wannier90 and PAOFLOW build from Quantum ESPRESSO.
For a pentagonal PdSe₂ monolayer, the 12-orbital models give the quantum spin
Hall plateau of the full 220-orbital calculation within 2 %, and the reduced
SIESTA model computes it about 14 times faster than an energy window on all
220 orbitals.

![Spin Hall conductivity of PdSe2 from five models of different size](fig_downfolding.png)

**FIG. 1.** σᶻ_xy of pentagonal PdSe₂ against the energy from the
valence-band maximum, in e²/h with the Pauli σ_z. Black solid line: SIESTA
with all 220 spinor orbitals, exact diagonalisation. Black dashed line: the
same SIESTA Hamiltonian with the states outside ±3 eV of E_F moved to the
window edges. Blue solid line: SIESTA reduced to its 12 isolated bands around
E_F. Green dash-dotted line: Quantum ESPRESSO with a 12-orbital Wannier90
model. Pink dotted line: Quantum ESPRESSO with a 108-orbital PAOFLOW model.
All but the first are KPM curves: 24 × 24 k mesh, M = 1024 moments, Jackson
kernel, 6 meV broadening. (a) The full energy range. (b) The quantum spin Hall
gap; the grey bands mark the gap of SIESTA (1.99–2.06 eV) and of Quantum
ESPRESSO (2.01–2.06 eV).

| Model | Orbitals per cell | Plateau (e²/h) | σᶻ_xy in the trivial gap (e²/h) | Moments |
|---|---|---|---|---|
| SIESTA, all states (exact) | 220 | 1.98 | +0.01 | — |
| SIESTA, window ±3 eV | 220 | 2.01 | −0.01 | 3 h 10 min |
| SIESTA, reduced to 12 bands | 12 | 1.98 | −0.06 | 14 min |
| QE + Wannier90 | 12 | 2.02 | −0.07 | 2.5 min |
| QE + PAOFLOW | 108 | 2.00 | +0.02 | 2 h 30 min |

<!-- Plateau: mean over the central half of each topological gap; trivial gap: E - E_VBM = 0.5 eV. Times: AMD EPYC 9374F, 16 threads (SIESTA) and 8 threads (Wannier90, PAOFLOW). The reduced SIESTA run also spends 58 min once on the orbital dipoles of the covariant velocity. Reduced SIESTA, exact on the same mesh: plateau 1.964, trivial gap -0.061. -->

In the quantum spin Hall gap all models agree. In the trivial gap the two
12-band models give −0.06 e²/h instead of +0.01: the bands they leave out add
a small Fermi-sea contribution there.

**Set it up and run it:** [Example 21 README](README.md).
