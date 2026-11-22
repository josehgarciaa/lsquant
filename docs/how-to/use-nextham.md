# Use a NextHAM model

Connect a Hamiltonian predicted by NextHAM to LSQuant and compute its density
of states or its spin Hall conductivity. NextHAM is an equivariant network
that predicts the Hamiltonian of ABACUS in its numerical atomic orbitals, with
spin–orbit coupling; the overlap S(R) comes from ABACUS. NextHAM writes its
prediction in the sparse-matrix format of ABACUS, and LSQuant reads that
format (`model.kind = "abacus"`), so the same run file also reads the DFT
Hamiltonian of ABACUS. The quality of the model (training data,
transferability, agreement with DFT) is yours. You should know NextHAM and
ABACUS and have read [Running LSQuant](run-lsquant.md) and
[Bring a model to LSQuant](bring-a-model.md).

Every key is listed in the
[run-file reference](../reference/run-file.md#machine-learning-hamiltonians-in-the-metric-executors).

## Requirements

- NextHAM and its ABACUS pipeline (`ABACUS2NextHAM`, `NextHAM2Bands`).
- LSQuant built and importable; see [Install LSQuant](install.md). Reading the
  files needs only NumPy.
<!-- format read from the sources of NextHAM 4397dc4 (NextHAM2Bands/output_hs_data.py, add_element.py, run_bands_compare.py) and pyatb d688b17 (io/abacus_read_xr.py, io/abacus_read_stru.py, core/base_data.cpp); unit tests use synthetic files, no real model yet -->

## 1. Turn on in NextHAM

Prepare the structure with NextHAM's ABACUS pipeline (`ABACUS2NextHAM`): it
runs ABACUS and writes, for each case, `abacus_scf/STRU` and, in
`abacus_geths/OUT.ABACUS/`, the overlap `srs1_nao.csr` and the ABACUS
Hamiltonian `hrs1_nao.csr`. Then save the predictions of the test run:
uncomment the line marked in NextHAM's README in `test.py` (line 440 in the
tested version). The test run writes one `.pth` file per case.

## 2. Transfer the model

NextHAM's `NextHAM2Bands` scripts turn the `.pth` file into H(R) in the
ABACUS format. `run_bands.sh` (or `run_bands_compare.py`) does it before it
plots the bands:

```bash
cd NextHAM2Bands
python run_bands_compare.py --pth test_res/case-0_out.pth \
    --case-dir ../ABACUS2NextHAM/simulate/case-0 \
    --efermi-pth ../ABACUS2NextHAM/processed_data/efermi_case-0.pth --out-dir band_res
```

It writes `band_res/pred_band/hr_tot` (the prediction) and
`band_res/gt_band/hr_tot` (the DFT Hamiltonian of the same case, through the
same steps). The Fermi energy, in eV, is the number stored in
`efermi_case-0.pth` (`python -c "import torch; print(float(torch.load('efermi_case-0.pth')))"`).

LSQuant reads three files: `hr_tot`, `srs1_nao.csr` and `STRU`. It also reads
the numerical orbital files named in `STRU` (block `NUMERICAL_ORBITAL`), which
give the number of orbitals and the angular momentum of each one.

## 3. Run LSQuant

```json
{
  "model": {
    "kind": "abacus",
    "file": "band_res/pred_band/hr_tot",
    "params": {
      "overlap_file": "../ABACUS2NextHAM/simulate/case-0/abacus_geths/OUT.ABACUS/srs1_nao.csr",
      "stru_file": "../ABACUS2NextHAM/simulate/case-0/abacus_scf/STRU",
      "orbital_dir": "/path/to/orbitals",
      "fermi_energy": 5.71,
      "velocity_mode": "berry_connection"
    }
  },
  "basis": {"name": "eigen", "supercell": [24, 24, 1]},
  "projection": {"kind": "shift", "window": [-3, 3]},
  "word": "Tr[JXSZ delta(E-H) VY dG-(E)/dE]",
  "expansion": {"moments": 1024, "window": "eigen", "trace": "exact"},
  "reconstruction": {"kernel": "jackson", "formulas": ["shc"], "broadening_meV": 10.0,
                     "energy_grid": [-1.0, 1.0, 401]},
  "output": {"dir": "out/nextham_shc"}
}
```

```bash
python -m lsquant inspect run.json
python -m lsquant run run.json
```

| Key | Value |
|---|---|
| `file` | H(R) in the ABACUS sparse format, text or binary: NextHAM's `hr_tot`, or ABACUS's own `hrs1_nao.csr` / `data-HR-sparse_SPIN0.csr`. |
| `overlap_file` | S(R) of the same structure (`srs1_nao.csr` / `data-SR-sparse_SPIN0.csr`). Required. |
| `stru_file` | The ABACUS `STRU` of the structure. Required. |
| `orbital_dir` | Directory of the numerical orbital files named in `STRU`. Default: the directory of `STRU`. |
| `orbitals_per_element` | Orbitals per atom, for example `{"Bi": 25, "Se": 13}`, if the orbital files are not available. Then `basis.json` holds no angular momenta. |
| `fermi_energy` | In eV, removed as H − E_F S. ABACUS energies are absolute; take the Fermi energy of the ABACUS run of the structure. Default 0. |
| `velocity_mode` | `"berry_connection"` or `"bare"`. There is no position matrix, so `"covariant"` is refused. |

Relative paths in `file`, `overlap_file`, `stru_file` and `orbital_dir` are
resolved like `model.file`. LSQuant converts Rydberg to eV and Bohr to Å with
the constants of pyatb, NextHAM's band code. ABACUS and NextHAM store
⟨i,0|H|j,R⟩ and sum e^{+ik·R}, as LSQuant does, and with spin–orbit coupling
the spin index runs fastest, which is LSQuant's interleaved layout; no sign or
order changes.

## Basis, positions and training data

- **Basis.** The ABACUS numerical atomic orbitals of the training data (the
  orbital files in `STRU`). A prediction is meaningful only for the elements
  and the orbital files of the training set; NextHAM's own table in
  `add_element.py` lists the number of orbitals it expects per element.
- **Positions.** The orbitals are centred on the atoms of `STRU`, so LSQuant
  puts every orbital on its atom. The intra-atomic dipoles are left out
  (`berry_connection`).
- **Training data.** Note the dataset of your model (the published
  Materials-SOC-HAM set or your own ABACUS data: functional, pseudopotentials,
  orbitals). The energy scale, the Fermi energy and the spin–orbit strength are
  those of that dataset.
- **Cut-off.** NextHAM predicts the blocks within 7 Å; `hr_tot` takes the
  blocks beyond it from ABACUS's `hrs1_nao.csr`. The model is partly predicted
  and partly DFT.

## What the files give and what LSQuant builds

| Quantity | From NextHAM and ABACUS | Built by LSQuant |
|---|---|---|
| H(R) | yes (`hr_tot`) | — |
| S(R) | yes, from ABACUS (`srs1_nao.csr`) | — |
| Velocity | no | from ∂H/∂k and ∂S/∂k with every orbital on its atom (`bare` or `berry_connection`); `covariant` needs the position matrix and is refused |
| Spin | no operator; the spinor layout of an `nspin = 4` file | σ on the spinor labels; a spinless file has no spin and the spin words stop with an error |
| Spin current | no | from the velocity and the spin |

## 4. Recommendations

- Run the same run file on `pred_band/hr_tot` and `gt_band/hr_tot`. The
  difference is the error of the prediction in the quantity you compute.
- Use `"oracle": true` on a small supercell to compare the kernel polynomial
  result with exact diagonalization.
- If LSQuant warns that H(R) differs from H(−R)† by more than 10⁻³, the
  prediction is not Hermitian: check the model before transport.

## References

- S. Yin, Z. Dai, X. Pan and L. He, "Advancing Universal Deep Learning for
  Electronic-Structure Hamiltonian Prediction of Materials", ICLR 2026.
- NextHAM source: <https://github.com/DavidYin94/NextHAM>.
- ABACUS sparse matrix output (`out_mat_hs2`): <https://abacus.deepmodeling.com>.
