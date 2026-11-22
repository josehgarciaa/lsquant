# One quantum spin Hall plateau, four codes

A pentagonal PdSe₂ monolayer has a quantum spin Hall gap about 2 eV above its
valence band. LSQuant computes the spin Hall conductivity in that gap from
four different electronic structures: SIESTA, OpenMX, and Quantum ESPRESSO
through Wannier90 and through PAOFLOW. All four give the same plateau.

![Spin Hall conductivity of PdSe2 from SIESTA, OpenMX, Wannier90 and PAOFLOW](fig_routes.png)

σᶻ_xy against energy from the valence-band maximum, in e²/h with the Pauli
σ_z. (a) Kernel polynomial method (KPM) for SIESTA, Wannier90 and PAOFLOW;
exact diagonalisation for OpenMX. (b) The quantum spin Hall gap, with the
exact curves (dotted). Grey bands: the gap of each electronic structure.

| Electronic structure | Plateau, KPM (e²/h) | Plateau, exact (e²/h) |
|---|---|---|
| SIESTA | 2.01 | 1.98 |
| QE + Wannier90 | 2.02 | 1.94 |
| QE + PAOFLOW | 2.00 | 1.94 (PAOFLOW's own) |
| OpenMX | — | 1.97 |

<!-- Plateau: mean over the central half of each topological gap, from the curves in data/ (plot.py constants VBM and GAP). KPM: 24x24, 1024 moments, Jackson kernel, 6 meV. Exact: SIESTA 24x24, Wannier90 300x300, OpenMX 60x60 (24x24 gives 1.850, not converged), PAOFLOW its own spin_Hall (120x120, adaptive smearing). -->

OpenMX places the gap 0.07 eV lower than the other codes, because its band
structure differs; its plateau has the same value.

**Set it up and run it:** [Example 20 README](README.md).
