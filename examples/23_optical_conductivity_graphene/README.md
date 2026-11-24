# Example 23: set up and run the optical tests of graphene and hBN

Compute the optical conductivity σ(μ, ω) of nearest-neighbour graphene and
hBN from Chebyshev moments, and check it against exact results. The
[example page](PAGE.md) shows the result and states the formula.

The figure is made from the files in [`data/`](data/), so you can plot it
without running anything (step 5). Steps 1 to 4 regenerate those files. Run
the commands from the repository root unless a step says otherwise.

## Requirements

- LSQuant built and importable; see [Install LSQuant](../../docs/how-to/install.md).
- SymPy, for step 1.
- For step 3: about 21 GB of memory per run, and about 4.5 min per velocity
  operator on 30 cores; see [Runtime](#runtime).

## Step 1: the derivation

```bash
python examples/23_optical_conductivity_graphene/derivation.py
```

The script starts from the finite-frequency Kubo–Bastin formula and derives,
with SymPy, each element of the formula on the page:

1. The Chebyshev expansion of δ(E − H).
2. The expansion of G^±(E), and which sign of the exponent belongs to G^+.
3. That the Bastin form equals the Kubo–Greenwood conductivity with time
   dependence e^{−iωt}, and its limits: Re σ ≥ 0, Im σ < 0 for an insulator,
   and the static value as ω → 0.
4. The Γ matrix, the complex conjugate and the prefactor 4 g_s/(π A D²).
5. The moment normalisation and the order of the velocities.
6. The units and the dielectric function.

It prints `PASS` in a few seconds.

## Step 2: the exact-trace check

```bash
python examples/23_optical_conductivity_graphene/exactness.py out/optical_exactness
```

[`exactness.py`](exactness.py) builds graphene, hBN and the Haldane model on
6 × 6 cells (72 orbitals). It computes their moments with an exact trace
(128 moments) and compares σ_xx and σ_xy from `optical_conductivity` with the
Kubo–Bastin formula evaluated on the eigenpairs of the same supercell. It also
checks the following:

- **The moment file:** μ_mn = 4 h_m h_n Tr[v_a T_m(H) v_b T_n(H)]/N, which is
  where the N/4 of the formula comes from.
- **The two polarisations:** the v_x v_x and v_y v_y moments are equal.
- **The Hall limit:** σ_xy(ω → 0) of the Haldane model equals the static
  Kubo–Bastin value.

Every comparison agrees to 10⁻¹¹ or better, and the script runs in seconds.

## Step 3: the moments of the large supercells

[`make_honeycomb.py`](make_honeycomb.py) writes H, v_x and v_y of the
honeycomb model as CSR files. The on-site energy is +Δ on one sublattice (B)
and −Δ on the other (N), the hopping is t, and the supercell has n × n cells
with periodic boundaries.

```bash
cd examples/23_optical_conductivity_graphene
python make_honeycomb.py graphene 400 operators
python make_honeycomb.py hbn 400 operators --delta 3.625 --t -2.33 --a 2.50
cd ../..
python -m lsquant run examples/23_optical_conductivity_graphene/run_graphene.json
python -m lsquant run examples/23_optical_conductivity_graphene/run_hbn.json
```

[`run_graphene.json`](run_graphene.json) and [`run_hbn.json`](run_hbn.json)
compute the moments of `Tr[VX delta(E-H) VX dG-(E)/dE]` and
`Tr[VY delta(E-H) VY dG-(E)/dE]`: 1400 × 1400 moments, 8 random-phase
states, and the window [−26, 26] eV (resolution 58 meV).

The window needs care. The Green functions are expanded at E ± ħω, so the
window must contain the spectrum widened by the largest photon energy on each
side. The spectra are |E| ≤ 8.10 eV (graphene) and 7.87 eV (hBN). The
reconstruction integrates over [−8.6, 8.6] eV, so photon energies up to
17 eV fit. The Kramers–Kronig test needs Re σ up to the end of the
absorption, above 2√(Δ² + t²) ≈ 8.6 eV, and the hBN polarisability is
plotted up to the gap; both need the wide window. `optical_conductivity` refuses a window that is too narrow and
names the range it needs.

## Step 4: the tests

```bash
cd examples/23_optical_conductivity_graphene
python benchmark.py ../../out/optical_graphene ../../out/optical_hbn
```

[`benchmark.py`](benchmark.py) reconstructs σ(μ, ω) for both velocities and
both materials and runs the tests the [page](PAGE.md) lists. The exact
references come from [`exact_honeycomb.py`](exact_honeycomb.py), which uses
the same model on a 1200 × 1200 k mesh. The script prints one line per test
and the count of passing tests (18 of 18), writes every number to
`benchmark.json`, and writes the curves to `bench_*.npz` in the current
directory.

The tolerances are as follows:

- 20 meV for the van Hove peaks.
- One resolution (58 meV) for the Pauli edge at 2μ.
- A quarter of a resolution for the hBN step, against the broadened exact
  step.
- 2 % for the graphene plateau and the polarisability.
- 1 % for the static limit.
- 3 % for the difference between the two polarisations.
- 3 % of max |σ| for the difference between Im σ and the Kramers–Kronig
  transform of Re σ, from 0.5 to 16 eV.

## Step 5: plot

Copy `bench_graphene.npz`, `bench_hbn.npz`, `bench_map.npz` and
`benchmark.json` to `data/`, or keep the committed ones, then run:

```bash
python examples/23_optical_conductivity_graphene/plot_benchmark.py
```

This writes `fig_benchmark.png`.

## The reconstruction in your own script

```python
from lsquant import _core
from lsquant.reconstruct import optical_conductivity, dielectric_function

moments = _core.ReadChebmom("out/optical_hbn/NonEqOpVX-VXhbnKPM_M1400x1400_statedefault.chebmom")
energies, sigma = optical_conductivity(
    moments, photon_energies, energy_window=(-8.6, 8.6),
    num_cells=400 * 400, cell_area=5.412659, spin_degeneracy=2)
# sigma[w, i]: complex, in e^2/h, at photon energy photon_energies[w] and chemical potential energies[i]
eps = dielectric_function(sigma, photon_energies, thickness=3.33)   # thickness in Angstrom
```

There are three conventions to keep in mind:

- **Time dependence:** e^{−iωt}. Re σ ≥ 0 is absorption, and an insulator
  has σ = −iωα with α > 0 below its gap.
- **Sheet conductivity:** σ is a 2D quantity. A dielectric function needs a
  thickness, which is a convention for a monolayer and not the height of the
  cell. The 2D polarisability α_2D = −e² Im σ/(ħω), shown in the figure,
  needs none.
- **Two velocities:** with moments of (v_a, v_b), the function returns σ_ab.
  The reference script returns σ_ba for these moments (see
  [Agreement with the reference script](#agreement-with-the-reference-script)).

## Agreement with the reference script

This reconstruction was ported from the post-processing script of T. Galvani
(private communication, 2026). The script takes `conj` of the final sum. For
one velocity the moments are real, and LSQuant and the script agree term by
term: on his 100 × 100-cell polycrystalline-graphene moments (256 × 256,
1024 energies, 500 photon energies) the largest difference is 2.4 × 10⁻¹⁰
e²/h, on values up to 21 e²/h. For two velocities, conj(μ^ab) = μ^ba, so the
script's `conj` returns σ_ba. LSQuant conjugates the moments first and
returns σ_ab, consistent with the static Hall conductivity
(`derivation.py`, step 5; `exactness.py`, Haldane model).

The script evaluates the double sum for every energy and photon energy. The
LSQuant kernel forms two half-sums once per energy:
O(N_E M² + N_E N_ω M) operations instead of O(N_E N_ω M²), which took 0.07 s
against 92 s on the same moments. LSQuant also accepts a window that is not
centred on zero; the script assumes a symmetric one.

## Runtime

Each run of step 3 keeps 320 000 orbitals and 1400 moments per axis. The
graphene and hBN runs ran side by side on 30 cores each of a shared node, in
about 16 min per velocity operator. Steps 1, 2, 4 and 5 take seconds
to a minute; most of step 4 is the 1200 × 1200 k mesh of the exact
references.

## What is not covered

- Finite temperature.
- A run-file formula for this reconstruction: it is a Python function.
- Excitons: hBN has strong excitonic effects, so the single-particle
  absorption tested here is not its measured spectrum.
