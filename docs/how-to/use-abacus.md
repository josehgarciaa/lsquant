# Use an ABACUS calculation

Start from your own ABACUS calculation in numerical atomic orbitals with
spin–orbit coupling and end with its spin Hall conductivity. LSQuant reads the
Hamiltonian, the overlap and the position matrix that ABACUS writes, with the
structure file `STRU`. Choose and converge the pseudopotentials, the orbitals
and the k-mesh in ABACUS, as usual. You should know ABACUS and have read
[Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

The model kind is `abacus`, the same one that reads NextHAM predictions in the
ABACUS format ([Use a NextHAM model](use-nextham.md)). An ABACUS model is
non-orthogonal, like a SIESTA or OpenMX model, and LSQuant runs it with the
same settings. It computes the density of states and the spin Hall
conductivity of the clean crystal on a k-mesh. Every key is listed in the
[run-file reference](../reference/run-file.md#machine-learning-hamiltonians-in-the-metric-executors).
<!-- source: python/lsquant/runfile.py ML_MODEL_KINDS["abacus"]; models/abacus.py (load_abacus_model, position_file) -->

## Requirements

- LSQuant built and importable; see [Install LSQuant](install.md).
- ABACUS with the LCAO basis (`basis_type lcao`).
<!-- tested with ABACUS 3.10.1, PseudoDojo NC FR 0.4 PBE standard and the ABACUS-orbitals Dojo-NC-FR DZP orbitals (Orbitals_v2.0) -->

## 1. Turn on the outputs in ABACUS

Set these keys in `INPUT`:

```text
basis_type      lcao
nspin           4       # spinors
lspinorb        1       # spin-orbit coupling, with fully relativistic pseudopotentials
noncolin        0       # non-magnetic: spin-orbit coupling without magnetisation
out_mat_hs2     1       # H(R) and S(R)
out_mat_r       1       # position matrix r(R), for the covariant velocity
```

For a magnetic material, set `noncolin 1` and the initial moments in `STRU`.
For a non-magnetic one keep `noncolin 0`: with `noncolin 1`, ABACUS can
converge to a solution with a spurious magnetisation.
<!-- observed on PdSe2: noncolin 1 with mag 0 0 0 in STRU gave TMAGZ 0.56 at the first SCF step; noncolin 0 keeps |m| ~ 1e-14 -->

ABACUS writes three files to `OUT.<suffix>/`: `data-HR-sparse_SPIN0.csr`,
`data-SR-sparse_SPIN0.csr` and `data-rR-sparse.csr`. Newer versions may
name the first two `hrs1_nao.csr` and `srs1_nao.csr`; LSQuant reads the
format, not the name.

## 2. Transfer the model

There is nothing to run. LSQuant reads the three matrix files, `STRU`, and
the numerical orbital files named in `STRU`, which give the number of
orbitals of each atom.

Give the Fermi energy of the ABACUS run as `model.params.fermi_energy`, in eV:
the second number of the last `E_Fermi` line of `OUT.<suffix>/running_scf.log`
(or `EFERMI = ... eV` at its end). LSQuant subtracts it, so every energy in
the run file is measured from it.
<!-- source: models/abacus.py docstring (energy reference: the absolute scale of ABACUS) -->

## 3. Run LSQuant

This run file computes the spin Hall conductivity of the clean crystal on a
24 × 24 k-mesh, with the spectrum shifted inside [−3, 3] eV:

```json
{
  "model": {
    "kind": "abacus",
    "file": "OUT.mysystem/data-HR-sparse_SPIN0.csr",
    "params": {
      "overlap_file": "OUT.mysystem/data-SR-sparse_SPIN0.csr",
      "position_file": "OUT.mysystem/data-rR-sparse.csr",
      "stru_file": "STRU",
      "orbital_dir": "/path/to/orbitals",
      "fermi_energy": -0.967,
      "velocity_mode": "covariant"
    }
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "metric": {"policy": "nonorthogonal", "route": "lowdin", "applier": "exact"},
  "projection": {"kind": "shift", "window": [-3, 3]},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 6.0,
                     "energy_grid": [-0.5, 2.0, 1251]},
  "output": {"dir": "out/abacus_shc"}
}
```

The settings are those of [Use an OpenMX calculation](use-openmx.md#3-run-lsquant),
which explains each of them. Paths are relative to the run file.

| Key | Values and meaning |
|---|---|
| `model.params.velocity_mode` | `"covariant"` needs `position_file`; without it, `"berry_connection"` (the default) or `"bare"`. |
| `model.params.position_file` | ABACUS's `data-rR-sparse.csr`. LSQuant builds the connection of the covariant velocity from it. |
| `model.params.orbital_dir` | Where the orbital files named in `STRU` are; by default the directory of `STRU`. |

Check the file, then run it:

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

The output files and units are those of the OpenMX and SIESTA runs, named
after the stem of `model.file` (ABACUS `hr_tot` writes `hr_totKPM_...`): σ^z_xy in
e²/h with the Pauli matrix σ_z; divide by 2 for spin ħ/2. See
[Results and units](use-siesta.md#results-and-units).

## 4. Recommendations

- Use orbitals generated for your pseudopotentials. The ABACUS orbital
  library has sets for the PseudoDojo fully relativistic pseudopotentials
  (`Dojo-NC-FR`); use at least DZP, with the cut-off radius the library
  recommends for each element, and `ecutwfc` equal to the energy cut-off of
  the orbitals (100 Ry for the v2.0 sets).
- For the spin Hall conductivity, use `velocity_mode = "covariant"`. ABACUS
  computes the position matrix with the same two-centre integrals as the
  overlap, so it is exact in the basis: LSQuant replaces it by its Hermitian
  part, and the removed part, written to `basis.json` as
  `connection_hermiticity_residual`, is about 10⁻⁸.
- For an insulator, the Fermi energy lies somewhere in the gap. To compare
  with another code, measure energies from the valence-band maximum.
- The cost grows with the number of bands per k-point: a DZP basis with
  spin–orbit coupling has 50 spinor orbitals per transition-metal or
  chalcogen atom (300 for the six atoms of PdSe₂). For that model, the exact
  spin Hall conductivity took about 0.7 s per k-point on 8 threads, and its
  plateau needed a 48 × 48 mesh: 24 × 24 gives a value 4 % too low.
<!-- measured: PdSe2 DZP (300 bands), BlochMetricBasis covariant + JXSZ/VY per k, 48x48 in 1667 s on 8 threads (AMD EPYC 9374F, shared node); plateau 1.870 (24x24), 1.948 (48x48), 1.948 (60x60), 1.948 (72x72) -->
