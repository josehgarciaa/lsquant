# Example 16: spectral projection of a SIESTA Hamiltonian onto an energy window

## The physical question

A SIESTA Hamiltonian in a DZP basis contains many bands far from the energies
that a transport calculation probes. For the Bi(111) bilayer of example 7, the
72 bands per $\mathbf k$ span $[-24.7, 56.4]$ eV, while transport happens within
a few eV of $E_F$. The Chebyshev expansion must nevertheless cover the whole
range. The energy resolution of $M$ moments is about $\pi a/M$, with $a$ the
half-width of the spectral window. An 80 eV window therefore needs about 17
times more moments than a 4.8 eV one for the same broadening.

A spectral projection removes the far bands from the matrix. It keeps the
basis, the overlap and every state inside the window unchanged. Take the
generalised eigenpairs $H(\mathbf k)\,c_n = \varepsilon_n S(\mathbf k)\,c_n$
and a window $[E_\mathrm{lo}, E_\mathrm{hi}]$. The `shift` projection replaces

$$ H \;\to\; H_s = S\,C\,\mathrm{diag}\big(\mathrm{clip}(\varepsilon_n,
E_\mathrm{lo}, E_\mathrm{hi})\big)\,C^\dagger S . $$

Every band outside the window moves to the nearest window edge, and the bands
inside keep their exact energies and eigenvectors. Operators are projected
consistently, $O \to P^\dagger O P$, where $P = C_W C_W^\dagger S$ is the
$S$-orthogonal projector on the kept bands. The `hard` kind moves the outside
bands to zero instead. The `orthogonal` kind builds a smaller orthonormal,
Wannier-like basis for an isolated group of bands.

We project the Bi bilayer onto $[-2, 2]$ eV and check that the spectrum inside
the window is unchanged.

## The model

We use the committed SIESTA files of example 7 on an $8\times8$ k mesh in the
`bloch` basis, with the Löwdin route. The projection is one section of
[`run_projection_shift.json`](run_projection_shift.json):

```json
"projection": {"kind": "shift", "window": [-2.0, 2.0]}
```

The window is given in eV, measured from the SIESTA Fermi level at $E = 0$.

## The run

We compute two words. `Tr[delta(E-H)]` is the full DOS. It still counts every
state, and the projected-out ones pile up at the window edges.
`Tr[WINDOW delta(E-H)]` is the DOS of the kept bands only. `WINDOW` is the
projected identity, and it must be declared in `model.operators`:

```json
"model": {"kind": "siesta", "dir": "examples/07_spin_hall_bi_bilayer/siesta", "label": "bi",
          "params": {"velocity_mode": "berry_connection"}, "operators": ["WINDOW"]},
"words": ["Tr[delta(E-H)]", "Tr[WINDOW delta(E-H)]"]
```

Every projected run also evaluates each word exactly, by diagonalisation with
a Gaussian broadening, before and after projection. This happens on an energy
grid that we choose inside the window:

```json
"reconstruction": {"kernel": "jackson", "broadening_meV": 50.0, "formulas": ["dos"],
                   "energy_grid": [-1.5, 1.5, 61]}
```

The run is started with:

```bash
python3 -m lsquant run examples/16_projection_downfold/run_projection_shift.json --out /tmp/ex16/projection_shift
```

For word $i$ the run writes the moments `word_<i>.chebmom1D` and `.chebmom`,
and the KPM DOS in states/eV/cell as `word_<i>.dat`. The file
`word_<i>_oracle.json` holds the exact curves by diagonalisation, before and
after projection. `projection.json` records the window, the number of kept
bands at each $\mathbf k$, and the spectral bounds before and after. It also
records the residuals of the projector: $P^2 = P$, $P^\dagger S = SP$, the
eigenpairs and $C^\dagger S C = 1$, all at $10^{-13}$ or below.

## The result

After projection the spectrum spans exactly $[-2, 2]$ eV. With the default
10 % padding, the KPM window is $[-2.4, 2.4]$ eV instead of about
$[-33, 64]$ eV. Between 2 and 8 bands per $\mathbf k$ fall inside the window,
4.53 on average. The number varies with $\mathbf k$, so this window cuts
through bands instead of isolating a fixed group. The `orthogonal` reduction
needs a constant number and would refuse this window. The `shift` and `hard`
kinds accept a number that varies.

The DOS inside the window shows the gap of the Bi bilayer at $E_F$. At
$E = 0$ the exact DOS is $8\times10^{-4}$ states/eV/cell at 50 meV
broadening. The KPM curve reads $3\times10^{-3}$ there, because the Jackson
kernel has longer tails than the Gaussian of the exact curve. The two curves
use different broadenings, and the run therefore labels its comparison of KPM
and exact curves "kernels differ" rather than giving a verdict.

## Checks

[`check_projection.py`](check_projection.py) reads the run directory. We run it
with:

```bash
python3 examples/16_projection_downfold/check_projection.py /tmp/ex16/projection_shift
```

Output on the reference run:

```text
[1] Tr[delta(E-H)]           exact DOS projected vs unprojected on |E| <= 1.5 eV: max rel diff = 2.79e-14  PASS
[1] Tr[WINDOW delta(E-H)]    exact DOS projected vs unprojected on |E| <= 1.5 eV: max rel diff = 2.78e-14  PASS
[2] spectrum before [-24.69, 56.37] eV, after [-2.000000, 2.000000] eV  PASS
[3] int Tr[delta] dE = 71.998 (72 spin-orbitals); int Tr[WINDOW delta] dE = 4.5312 (mean kept bands 4.53125)  PASS
[4] |E| < 1.5 eV: max |DOS - WINDOW DOS| = 0.0064 (max DOS 5.805) states/eV/cell  PASS
CHECK PASS
```

1. Inside the window the projection changes nothing: the exact DOS before and
   after agree to $3\times10^{-14}$.
2. Line [2] prints the spectral bounds, $[-24.69, 56.37]$ eV before
   projection and exactly $[-2, 2]$ eV after. The projection has moved every
   band outside the window onto the window edges, as the `shift` formula of
   [the physical question](#the-physical-question) prescribes.
   [The result](#the-result) gives the KPM window that follows from these
   bounds.
3. The full DOS still integrates to the 72 spin-orbitals per cell, and the
   projected identity integrates to the average number of kept bands.
4. Away from the window edges, where the full DOS holds the piled-up outside
   states, the two words give the same curve.

## Runtime

On an oversubscribed 96-core node, where the process received about 10 % of a
core, the run took 856 s with one OpenMP thread and the check under 1 s. Of
the 856 s, 150 s went to the moments of the first word and 659 s to those of
the `WINDOW` word.
