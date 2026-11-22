# Use a ThreeBodyTB model

Connect a ThreeBodyTB (TB3) tight-binding model to LSQuant and compute its
bands and its density of states. ThreeBodyTB has no spin–orbit coupling, so
use it for bands and DOS, not for spin Hall transport. ThreeBodyTB builds a non-orthogonal tight-binding model of
a crystal from two- and three-body terms fitted to DFT, for most of the
periodic table, without a DFT calculation of the structure. LSQuant reads the
model file that ThreeBodyTB writes; the quality of the model (fit, range of
validity, self-consistency) is yours. You should know ThreeBodyTB and have
read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

Every key is listed in the
[run-file reference](../reference/run-file.md#machine-learning-hamiltonians-in-the-metric-executors).

## Requirements

- ThreeBodyTB.jl in Julia, to build the model.
- LSQuant built and importable; see [Install LSQuant](install.md). Reading the
  file needs only the Python standard library and NumPy.
<!-- format read from the sources of ThreeBodyTB.jl 0752192 (src/TB.jl write_tb_crys, read_tb_crys, Hk; src/Crystal.jl orbital_index; src/Atomdata.jl; src/GlobalUnits.jl); unit tests use synthetic XML files, no real model yet -->

## 1. Turn on in ThreeBodyTB

Nothing to change. Build the model of your structure as usual. For a system
with charge transfer, run the self-consistent step first, so that the file
holds the charge term:

```julia
using ThreeBodyTB
crys = makecrys([4.0 0 0; 0 4.0 0; 0 0 4.0], [0 0 0; 0.5 0.5 0.5], ["Na", "Cl"])
energy, tbc, flag = scf_energy(crys)
```

## 2. Transfer the model

Write the model with ThreeBodyTB's function `write_tb_crys`:

```julia
write_tb_crys("mymodel.xml", tbc)
```

The XML file holds H(R) and S(R), the lattice, the atomic positions, the
self-consistent charge term and the Fermi energy. LSQuant also reads it
gzipped (`mymodel.xml.gz`). ThreeBodyTB's `write_hr_dat` writes a
Löwdin-orthogonalized Wannier90 file on a k grid instead; use the XML file,
which keeps the model exactly.

## 3. Bands and DOS with LSQuant

**Bands.** Read the model in Python and solve H(k) c = E S(k) c along a path
of fractional k points:

```python
import numpy as np
from scipy.linalg import eigh
from lsquant.models.threebodytb import load_threebodytb_model

model = load_threebodytb_model("mymodel.xml")          # energies in eV from E_F
nodes = np.array([[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0], [0, 0, 0]])   # Γ X M Γ
path = np.concatenate([np.linspace(a, b, 60, endpoint=False) for a, b in zip(nodes[:-1], nodes[1:])]
                      + [nodes[-1:]])
k_cart = path @ (2 * np.pi * np.linalg.inv(model.uc).T)
H, S = model.Hk(k_cart), model.Sk(k_cart)
bands = np.array([eigh(h, s, eigvals_only=True) for h, s in zip(H, S)])   # (n_k, n_bands)
```

The bands are those of ThreeBodyTB's `Hk` at the same k points: the charge
term is included and the Fermi energy of the file is the zero.

**DOS.** Run LSQuant with a run file:


```json
{
  "model": {
    "kind": "threebodytb",
    "file": "/path/to/mymodel.xml",
    "params": {"velocity_mode": "berry_connection"}
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 24]},
  "word": "Tr[delta(E-H)]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["dos"], "broadening_meV": 20.0,
                     "energy_grid": [-6.0, 6.0, 601]},
  "output": {"dir": "out/tb3_dos"}
}
```

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Value |
|---|---|
| `velocity_mode` | `"berry_connection"` or `"bare"`. There is no position matrix, so `"covariant"` is refused. |
| `fermi_energy` | In eV, removed as H − E_F S. Default: the Fermi energy stored in the file (`efermi`, converted from Rydberg). |
| `orbitals_per_element` | Number of orbitals of an element, for example `{"Xx": 9}`. Needed only for an element missing from ThreeBodyTB's table, or for a modified database. |

LSQuant converts Rydberg to eV and Bohr to Å with ThreeBodyTB's own
constants, puts every orbital on its atom, and adds the self-consistent charge
term as ThreeBodyTB does when it builds H(k). ThreeBodyTB sums
conj(H(R)) e^{−2πik·R}; LSQuant stores conj(H(−R)), which gives the same H(k)
with its convention e^{+2πik·R}.

## Basis, positions and training data

The velocity that LSQuant builds is only as good as the basis and the orbital
positions:

- **Basis.** ThreeBodyTB uses a minimal set of atomic orbitals per element
  (s and p for main-group elements, s, d and p for transition metals; the table
  of `Atomdata.jl`). The labels of the orbitals (s, pz, px, py, dz2, ...) are in
  `basis.json` as `orbital_labels`.
- **Positions.** Every orbital sits on its atom; the position operator is
  diagonal and the intra-atomic dipoles are left out (`berry_connection`).
- **Training data.** ThreeBodyTB is not a neural network: its two- and
  three-body parameters are fitted to a large set of Quantum ESPRESSO
  calculations (PBEsol, GBRV pseudopotentials) projected on atomic orbitals.
  Its energies, band gaps and Fermi energy are those of that fit. For
  transport, prefer an equivariant network trained on your material when you
  have one; ThreeBodyTB is the choice when no trained network covers your
  elements.

## What the file gives and what LSQuant builds

| Quantity | From the ThreeBodyTB file | Built by LSQuant |
|---|---|---|
| H(R) | yes, with the charge term | — |
| S(R) | yes (`nonorth = true`, the default of ThreeBodyTB); without it the model is orthogonal | — |
| Velocity | no | from ∂H/∂k and ∂S/∂k with every orbital on its atom (`bare` or `berry_connection`); `covariant` needs a position matrix and is refused |
| Spin | no; ThreeBodyTB has no spin–orbit coupling | σ on the spinor labels, only for a file with two spin channels; otherwise the spin words stop with an error |
| Spin current | no | from the velocity and the spin, so only when the spin exists |

## 4. Recommendations

- ThreeBodyTB has no spin–orbit coupling: compare bands and DOS, not spin
  transport. The spin Hall word stops with an error. A file with two spin
  channels (a magnetic model) gives a model with spin up and spin down and no
  mixing.
- A file written with `scfspin = true` is refused: `write_tb_crys` does not
  store the magnetic charge term, so the model in the file is incomplete.
- Compare the bands with ThreeBodyTB's own `plot_bandstr` on the same path,
  and the LSQuant DOS with ThreeBodyTB's `dos` on the same k mesh.
- ThreeBodyTB warns when a structure lies outside its fitting range
  (`within_fit = false`). Trust the transport only for structures within the
  fit.

## References

- K. F. Garrity and K. Choudhary, "Fast and accurate prediction of material
  properties with three-body tight-binding model for the periodic table",
  Phys. Rev. Materials 7, 044603 (2023).
- ThreeBodyTB.jl source and documentation:
  <https://github.com/usnistgov/ThreeBodyTB.jl>.
