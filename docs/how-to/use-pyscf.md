# Use a PySCF periodic calculation

Start from your own periodic PySCF calculation with spin–orbit coupling and
end with its spin Hall conductivity. LSQuant reads the Kohn–Sham matrix, the
overlap and the position matrix of the Gaussian basis. Choose and converge the
basis, the k-mesh and the integration grids in PySCF, as usual. You should
know PySCF and have read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

A PySCF model is non-orthogonal, like a SIESTA or OpenMX model, and LSQuant
runs it with the same settings and the same keys. It computes the density of
states and the spin Hall conductivity of the clean crystal on a k-mesh; a
real-space supercell or disorder stops the run. Every key is listed in the
[run-file reference](../reference/run-file.md#pyscf-models).

## Requirements

You need:

- LSQuant importable in the Python environment of PySCF; see
  [Install LSQuant](install.md).
- PySCF 2.x with periodic two-component DFT (`pyscf.pbc.dft.KGKS`).
<!-- tested with PySCF 2.15.1 -->

## 1. Turn on spin–orbit coupling in PySCF

Run a two-component calculation on a full Γ-centred Monkhorst–Pack mesh, for
example all-electron with the X2C Hamiltonian:

```python
from pyscf.pbc import gto, dft, df
from pyscf.pbc.dft import gen_grid

cell = gto.Cell(a=..., atom=..., basis="x2c-svpall", unit="A")
cell.build()
kpts = cell.make_kpts([6, 6, 1])           # 2D: one k point along the vacuum
mf = dft.KGKS(cell, kpts).x2c1e()          # two-component, spin-orbit X2C
mf.xc = "pbe"
mf.grids = gen_grid.BeckeGrids(cell)       # atom-centred grids for all-electron
mf.with_df = df.GDF(cell, kpts)
mf.kernel()
```

A slab is computed as a 3D cell with vacuum and one k point along it.

## 2. Transfer the model

After the SCF, write the model:

```python
from lsquant.models.pyscf import write_pyscf_model

write_pyscf_model(mf, "mymodel.npz")
```

The file contains H(R), S(R) and the position matrix Q(R) = ⟨i,0|r|j,R⟩.
S and Q are analytic Gaussian integrals, so the covariant velocity needs no
quadrature. H(R) is the Fourier transform of the Kohn–Sham matrix on the SCF
mesh: H(k) is exact on that mesh and interpolated between its points. The
file stores `edge_ratio`, the size of H(R) at the edge of the supercell of the
mesh relative to its largest block; above 10⁻³ a warning asks for a denser
mesh. `mesh=(n1, n2, 1)` builds H on a denser mesh from the converged density,
at the cost of the core Hamiltonian at each new k point; this needs a
non-magnetic density and a pure functional.

Energies are measured from `fermi_energy`, by default the middle of the gap;
for a metal, pass it. Give every energy in the run file relative to it. The
model takes no `model.params.fermi_energy` and no `model.params.spin`.

## 3. Run LSQuant

The run file is that of an OpenMX model with `"kind": "pyscf"` and `"file"`
set to the model; see
[Use an OpenMX calculation](use-openmx.md#3-run-lsquant):

```json
{
  "model": {"kind": "pyscf", "file": "mymodel.npz",
            "params": {"velocity_mode": "covariant", "oracle": true, "bloch_memory_mb": 4096}},
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "metric": {"policy": "nonorthogonal", "route": "lowdin", "applier": "exact"},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 6.0,
                     "energy_grid": [-0.5, 2.0, 1251]},
  "output": {"dir": "out/pyscf_shc"}
}
```

`bloch_memory_mb` raises the workspace of the Bloch sums: an all-electron
basis has large blocks (468 spinor functions per cell for PdSe₂ with
x2c-SVPall), and the default of 256 MiB cannot hold them.

## 4. Recommendations

- Check the bands of the model against PySCF's own (`mf.get_bands`) at k points
  between those of the mesh before computing transport.
- All-electron Gaussian bases built for molecules carry diffuse functions;
  check the smallest eigenvalue of S(k) and, if it is small, remove them or use
  a basis made for solids.
- The X2C core Hamiltonian is the expensive step of an all-electron
  spin–orbit calculation, minutes per k point for 4d elements; keep the mesh of
  the model equal to the SCF mesh unless the bands between its points are off.
