# One quantum spin Hall plateau, four routes from DFT

The quantum spin Hall effect tests a spin transport calculation as a whole.
Its plateau comes from the Fermi sea, so the occupied bands must carry the
right topology, as for the anomalous Hall plateau of
[example 22](../22_yn2_anomalous_hall/PAGE.md). It is also the response of a
spin current, so the velocity, the spin operator and the spin current built
from them, $J = \{v_x, \sigma_z\}/2$, must be right in the basis of each code. A
connector that passes this test reads the Hamiltonian, the overlap, the
positions and the spin of its code correctly.

A pentagonal PdSe₂ monolayer is a semiconductor whose lowest conduction bands
are topologically nontrivial [1]: a quantum spin Hall gap opens about 2 eV
above the valence-band maximum (VBM); electron doping or gating would bring
the Fermi level there. We compute the spin Hall conductivity in that gap from
three DFT codes and four routes: SIESTA, OpenMX, and Quantum ESPRESSO through
Wannier90 and through PAOFLOW. All four find the plateau, and the plateaus
agree within 4 %.

![Bands of PdSe2 from Quantum ESPRESSO and OpenMX](fig_bands.png)

**FIG. 1.** Bands of pentagonal PdSe₂ from the VBM of each code, PBE with
spin–orbit coupling. Green circles: Quantum ESPRESSO (plane waves) at the
points of its 12 × 12 mesh that lie on the path; green dashed line: its
12-orbital Wannier90 model. Blue solid line: OpenMX (numerical atomic
orbitals). Grey bands: the quantum spin Hall gap of each code,
Quantum ESPRESSO 2.01–2.06 eV and OpenMX 1.92–1.98 eV.

![Spin Hall conductivity of PdSe2 from SIESTA, OpenMX, Wannier90 and PAOFLOW](fig_routes.png)

**FIG. 2.** σᶻ_xy of pentagonal PdSe₂ against the energy from the
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

[1] S. Bravo, M. Pacheco, J. D. Correa and L. Chico, Topological bands in
the PdSe₂ pentagonal monolayer, *Phys. Chem. Chem. Phys.* (2022),
[doi:10.1039/D2CP01822E](https://doi.org/10.1039/D2CP01822E).

**Set it up and run it:** [Example 20 README](README.md).

## Check

Verify the published plateaus from the committed curves, from the repository root:

```bash
python3 examples/20_pdse2_spin_hall_routes/check.py --data-only
```

Compare completed SIESTA outputs with their matching committed curves:

```bash
python3 examples/20_pdse2_spin_hall_routes/check.py --siesta-direct out/pdse2_routes --siesta-projected out/pdse2_routes_projected --siesta-reduced12 out/pdse2_reduced12
```

Supply only the output directories available. Omitted routes are reported as
not checked. Every supplied route requires its topology JSON; the current
projection execution path does not produce that file, which is reported as a
failed input requirement. The direct reconstruction has no committed matching
reference, so its curve comparison is explicitly reported as not checked.
