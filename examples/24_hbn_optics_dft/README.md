# Example 24: set up and run the optical conductivity of monolayer hBN from DFT

This example computes σ(ω) of monolayer hexagonal boron nitride from seven
DFT Hamiltonians: SIESTA (PBE), PySCF (PBE and HSE06), and Quantum ESPRESSO
(PBE and HSE06) projected on atomic orbitals with Wannier90 and with
PAOFLOW. Each result from the Chebyshev moments is compared with the exact
sum over states of the same model. The [page](PAGE.md) shows the results;
this file says how to reproduce them.

## Requirements

- LSQuant built with its Python package.
- Quantum ESPRESSO 7.4 (`pw.x`, `open_grid.x`, `projwfc.x`, `pw2wannier90.x`).
- Wannier90 4.0 or later: the projectability disentanglement
  (`dis_froz_proj`) and `auto_projections` are new in version 4.
- PAOFLOW 3.0.
- SIESTA 5.4 with `sisl` in the LSQuant environment.
- PySCF 2.15 (a separate environment works; the writer needs `lsquant` importable).
- Pseudopotentials: PseudoDojo nc-sr-04 PBE standard for B and N, as
  `pseudo/B.upf`, `pseudo/N.upf` (QE) and `pseudo/B.psml`, `pseudo/N.psml`
  (SIESTA). They are not distributed here.

## The structure

All codes use the same cell: a = 2.504 Å, 20 Å along z, B at (1/3, 2/3, 1/2)
and N at (2/3, 1/3, 1/2) in fractional coordinates. No spin-orbit coupling.

## Step 1: Quantum ESPRESSO

In `qe/pbe` and `qe/hse` (the inputs differ only in the functional: HSE06
adds `input_dft = 'hse'`, a 6 × 6 q mesh for the exchange and `ecutfock`):

```bash
mpirun -np 16 pw.x -nk 4 -in scf.in > scf.out          # 12 x 12 k, 40 bands, 90 Ry
mpirun -np 4 open_grid.x -in open_grid.in > open_grid.out
mpirun -np 4 projwfc.x -in projwfc.in > projwfc.out    # PAOFLOW reads the overlaps
```

`open_grid.x` unfolds the symmetry-reduced scf states to the full 12 × 12
mesh, which Wannier90 and PAOFLOW need. The run prints the highest occupied
level (PBE −4.0569 eV, HSE06 −4.8265 eV); `make_runs.py` uses it as the
energy origin.

The lowest conduction states at Γ are not π* states. They are
nearly-free-electron states of the vacuum, 0.65 eV above the valence-band
maximum in PBE, and they lie between the π* and σ* bands. No energy window
separates them from the atomic-orbital bands, which is why both projections
select by projectability.

## Step 2: projected Wannier functions

In `qe/pbe` (FROZEN = −3.0) and `qe/hse` (FROZEN = −3.8, both in the gap):

```bash
python ../make_win.py open_grid.out FROZEN
wannier90.x -pp hbn
mpirun -np 8 pw2wannier90.x -in pw2wan.in > pw2wan.out
wannier90.x hbn
```

Eight functions (s and p on B and N) from the projections on the
pseudo-atomic orbitals (`auto_projections`, `atom_proj`), disentangled by
projectability, without maximal localisation (`num_iter = 0`):
[`make_win.py`](qe/make_win.py) explains the choices. Check the result with

```bash
python ../check_wannier.py hbn
```

which compares the bands with the DFT ones on the mesh and prints the decay
of H(R). In PBE the valence bands agree to 2 × 10⁻⁵ eV, the π* band to
0.01 eV at K and M and within 2 meV on most of the mesh, and H(R) is
below 5 meV beyond 11 Å. Around Γ, where the π* band lies 10 eV above the
valence-band maximum among the vacuum states, it deviates by up to 0.8 eV
(38 % of the mesh points by more than 50 meV).

## Step 3: projected PAOFLOW

In `qe/pbe` and `qe/hse`:

```bash
python ../run_paoflow.py hbn_paoflow_pbe.npz --pthr 0.6     # hbn_paoflow_hse.npz in qe/hse
```

PAOFLOW keeps the bands whose projectability, averaged over k, exceeds
`pthr`. The vacuum state at Γ lowers the average of the π* band to 0.70, so
the usual 0.95 keeps only the valence bands; 0.6 keeps the π* band as well.
The script prints the PAO and DFT bands at Γ, M and K: valence and π* agree
at M and K; the other three bands sit at PAOFLOW's shift (5.6 eV in PBE,
6.5 eV in HSE06 above the valence-band maximum).

## Step 4: SIESTA and PySCF

```bash
cd siesta && siesta < hbn.fdf > hbn.out                  # PBE, DZP, 24 x 24 k
cd pyscf && python hbn_pyscf.py pbe hbn_pyscf_pbe.npz    # gth-dzvp, 12 x 12 k
cd pyscf && python hbn_pyscf.py hse hbn_pyscf_hse.npz
```

SIESTA writes `hbn.HSX`; LSQuant reads it through `sisl`. The PySCF script
writes the LSQuant model directly.

## Step 5: LSQuant runs

From the repository root:

```bash
python examples/24_hbn_optics_dft/make_runs.py
for r in siesta_pbe pyscf_pbe pyscf_hse w90_pbe w90_hse pao_pbe pao_hse; do
  python -m lsquant run examples/24_hbn_optics_dft/run_$r.json
done
```

Each run computes the moments of `Tr[VX delta(E-H) VX dG-(E)/dE]` and
`Tr[VY delta(E-H) VY dG-(E)/dE]` on a 120 × 120 k mesh with an exact trace,
with a resolution of 60 meV. [`make_runs.py`](make_runs.py) explains the
windows and the velocity of each model. The two local-basis codes need a
spectral projection: their basis reaches 100–190 eV, and only one gap
separates their bands below 60 eV (the band gap itself), so an isolated
window that keeps conduction bands does not exist. The `shift` projection
on [−21, 22] eV projects the velocity on the states inside and clips the
others to the window edges; transitions to states above 22 eV are left
out.

## Step 6: exact references

From the repository root:

```bash
E=examples/24_hbn_optics_dft; mkdir -p $E/sos
for r in siesta_pbe pyscf_pbe pyscf_hse; do python $E/sos.py $E/run_$r.json 1 $E/sos/$r.npz; done
for r in w90_pbe w90_hse pao_pbe pao_hse; do python $E/sos.py $E/run_$r.json 2 $E/sos/$r.npz; done
```

[`sos.py`](sos.py) diagonalises each model on the same mesh, with the same
velocity and projection window, and evaluates the interband Kubo formula
with a Gaussian of the width of one Chebyshev transition. The second
argument is the spin degeneracy: 1 for the spinor-lifted SIESTA and PySCF
models, 2 for the spinless QE models.

## Step 7: analysis and figure

```bash
mkdir -p data && cd data
python ../analyze.py ../../../out ../sos
cd .. && python plot.py
```

[`analyze.py`](analyze.py) reconstructs σ_xx and σ_yy at a chemical
potential in the middle of the gap and runs the checks of the page. It
writes `analysis.json` and `curves.npz`; [`plot.py`](plot.py) writes
`fig_hbn_dft.png`.

The lower limit of the integration range passed to `optical_conductivity`
(−22 eV here) must lie below every state of the moments, including the
states a `shift` projection piles at its lower edge. In HSE06 the N 2s
band of PySCF starts below −21 eV; with −21 eV as the limit the Fermi sea
is cut through that pile, and α_2D(ω → 0) diverges.

## Runtime

On one shared 64-core node: QE PBE 2 min and QE HSE06 2 h (scf with 40
bands); SIESTA 8 min; PySCF 10 min (PBE) and 2.4 h (HSE06). Each QE
model run of step 5 takes about 25 min per word on 4 cores. The SIESTA and
PySCF runs keep all 52 spinor states per k and take several hours per word
on 16 cores; their metric (S⁻¹ᐟ²) and projection take 15 and 7 min.

## What is not covered

- Excitons. hBN absorbs mainly through a bound exciton below the
  quasiparticle gap; independent-particle σ(ω) has no electron-hole
  interaction.
- Local fields, and the out-of-plane polarisation.
- Transitions to states outside each model: the vacuum states for every
  model, states above 22 eV for SIESTA and PySCF, everything above π* except
  σ* for Wannier90, and everything above π* for PAOFLOW.
