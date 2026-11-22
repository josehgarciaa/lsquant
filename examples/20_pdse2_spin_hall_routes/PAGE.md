# One quantum spin Hall plateau, four routes from DFT

A pentagonal PdSe₂ monolayer is a semiconductor with a quantum spin Hall gap
about 2 eV above its valence-band maximum, in the conduction bands; electron
doping or gating would bring the Fermi level there. We compute the spin Hall
conductivity in that gap from three DFT codes and four routes: SIESTA, OpenMX,
and Quantum ESPRESSO through Wannier90 and through PAOFLOW. The plateaus agree
within 4 %.

![Spin Hall conductivity of PdSe2 from SIESTA, OpenMX, Wannier90 and PAOFLOW](fig_routes.png)

**FIG. 1.** σᶻ_xy of pentagonal PdSe₂ against the energy from the
valence-band maximum, in e²/h with the Pauli σ_z. (a) Black solid line:
SIESTA, states within ±3 eV of E_F. Green dashed line: Quantum ESPRESSO with a
Wannier90 model. Pink dotted line: Quantum ESPRESSO with a PAOFLOW model.
These three are KPM curves: 24 × 24 k mesh, M = 1024 moments, Jackson kernel,
6 meV broadening. Blue dash-dotted line: OpenMX, exact diagonalisation on a
60 × 60 mesh with a 6 meV interband regulator and a 6 meV Gaussian. (b) The
quantum spin Hall gap, with the same curves and the exact results (thin lines
with markers): SIESTA (circles, 24 × 24), Wannier90 (squares, 300 × 300) and
PAOFLOW's own calculation (triangles, 120 × 120). The grey bands mark the gap
of each electronic structure: OpenMX 1.92–1.98 eV, SIESTA 1.99–2.06 eV,
Quantum ESPRESSO 2.01–2.06 eV.

| Route | Plateau, KPM (e²/h) | Plateau, exact (e²/h) |
|---|---|---|
| SIESTA | 2.01 | 1.98 |
| QE + Wannier90 | 2.02 | 1.94 |
| QE + PAOFLOW | 2.00 | 1.94 (PAOFLOW's own) |
| OpenMX | — | 1.97 |

Only for SIESTA do the KPM and exact columns use the same k mesh; for the
other routes the exact result uses a denser mesh, so the difference between
the columns mixes the mesh and the KPM. The KPM run of OpenMX is not shown: it
keeps every band of the model and was stopped after 4 h, so its curve comes
from the exact calculation.

<!-- Plateau: mean over the central half of each topological gap, from the curves in data/ (plot.py constants VBM and GAP). OpenMX on a 24x24 mesh gives 1.850, not converged; 60x60 gives 1.974 without the regulator and 1.971 with it. -->

OpenMX places the gap 0.07 eV lower than SIESTA and Quantum ESPRESSO; its
plateau, 1.97 e²/h, agrees with theirs.

**Set it up and run it:** [Example 20 README](README.md).
