# Example 19: Five observables of a quantum spin Hall insulator from one run file

## The physical question

The run file names the quantity to compute in one line, the `word` key. We take
one model and change only that line to obtain five observables: the density of
states, the spin-resolved density of states, the longitudinal conductivity
$\sigma_{xx}$, the Hall conductivity $\sigma_{xy}$ and the spin Hall conductivity
$\sigma^z_{xy}$. The model is a quantum spin Hall insulator, so the five curves
answer one physical question: what remains of the quantised spin Hall response
when the spin is no longer conserved and time-reversal symmetry is broken?
Every curve is compared with an exact calculation on the same k-point mesh.

## The model

We use the Kane–Mele model of a honeycomb crystal, with the bond length
$a = 1$ Å and two $p_z$-like orbitals per cell, each with two spin states
(four states per cell):

$$
H = -t\sum_{\langle ij\rangle} c_i^\dagger c_j
  + i\lambda_{\mathrm{SO}}\sum_{\langle\langle ij\rangle\rangle}\nu_{ij}\,c_i^\dagger\sigma_z c_j
  + i\lambda_{\mathrm R}\sum_{\langle ij\rangle} c_i^\dagger(\boldsymbol\sigma\times\hat{\mathbf d}_{ij})_z c_j
  + M\sum_i c_i^\dagger\sigma_z c_i ,
$$

with $t = 1$ eV, intrinsic spin–orbit coupling $\lambda_{\mathrm{SO}} = 0.1$ eV,
Rashba coupling $\lambda_{\mathrm R} = 0.05$ eV and an exchange field
$M = 0.05$ eV along $z$. The Rashba term mixes the two spins, so $\sigma^z_{xy}$
is no longer quantised. The exchange field breaks time-reversal symmetry, so the
spin density of states and the Hall conductivity can be non-zero.

The script `make_kane_mele.py` writes the model as sparse matrices for a
$16 \times 16$ supercell, periodic in both directions, with 1024 states. It
writes five matrices: the Hamiltonian, the velocities $V_x$ and $V_y$, the spin
operator $\sigma_z$ and the spin current $J^z_x = \{V_x, \sigma_z\}/2$. A periodic
$16 \times 16$ supercell has exactly the k points of a $16 \times 16$ mesh. The
script `oracle.py` diagonalises the Bloch Hamiltonian on that mesh and gives the
exact reference for every curve.
<!-- source: examples/19_observable_words/make_kane_mele.py (parameters at the top), oracle.py -->

## The run

The run file `run_words.json` loads the five matrices as a `csr` model and lists
five values of the `word` key:

| `word` | Observable | Curve file prefix |
|---|---|---|
| `Tr[delta(E-H)]` | density of states $\rho(E)$ | `meanSpectralOp1` |
| `Tr[SZ delta(E-H)]` | spin density of states $\rho_z(E) = \rho_\uparrow - \rho_\downarrow$ | `meanSpectralOpSZ` |
| `Tr[VX delta(E-H) VX delta(E-H)]` | Kubo–Greenwood $\sigma_{xx}(E_F)$ | `KuboGreenWood_` |
| `Tr[VX delta(E-H) VY dG-(E)/dE]` | Kubo–Bastin $\sigma_{xy}(E_F)$ | `KuboBastin_NonEqOpVX-VY` |
| `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` | spin Hall $\sigma^z_{xy}(E_F)$ | `SpinHall_` |

A single calculation uses the `word` key alone. Here the `words` list runs all
five words in one call, and every other key of the run file stays the same:

```json
"word": "Tr[delta(E-H)]",
"words": ["Tr[delta(E-H)]", "Tr[SZ delta(E-H)]",
          "Tr[VX delta(E-H) VX delta(E-H)]",
          "Tr[VX delta(E-H) VY dG-(E)/dE]",
          "Tr[JXSZ delta(E-H) VY dG-(E)/dE]"],
"expansion": {"moments": 512, "window": [-4.0, 4.0], "trace": "exact"},
"reconstruction": {"kernel": "jackson", "broadening_meV": 20.0,
                   "formulas": ["dos", "greenwood", "bastin"],
                   "num_cells": 256, "cell_area": 2.598076211353316}
```

The number of moments sets the energy resolution, like a smearing: with 512
moments on the window of half-width 4 eV the resolution is about
$\pi \cdot 4\ \mathrm{eV}/512 \approx 25$ meV. We ask for a 20 meV broadening.
The reconstruction keeps $\lceil \pi W/(2\eta) \rceil$ moments for a window of
width $W = 8$ eV and a broadening $\eta$, which would be 628 moments here. The run
has 512, so all 512 are used and the curves keep the resolution of 512 moments,
about 25 meV, rather than 20 meV.
<!-- source: python/lsquant/reconstruct.py dos_curve (eta = 2 b / W, M_eff = min(ceil(pi/eta), M)) and _truncate_and_damp_2d (same rule per axis) -->
The trace is exact, so the curves carry no random-vector
noise. The supercell size sets the finite-size error, like a k-mesh.

From the repository root we generate the model, run the five words, build the
spin Hall curve, compute the reference and plot:

```bash
python3 examples/19_observable_words/make_kane_mele.py 16
python3 -m lsquant run examples/19_observable_words/run_words.json
python3 -m lsquant reconstruct \
    examples/19_observable_words/out/NonEqOpJXSZ-VYkm_L16KPM_M512x512_stateexact.chebmom \
    --formula shc --broadening-meV 20
cd examples/19_observable_words && python3 oracle.py 16 20 && cd ../..
python3 examples/19_observable_words/plot.py
```

`reconstruct` reads `num_cells` and `cell_area` from `run.resolved.json` in the
output directory and writes $\sigma^z_{xy}$ in units of $e^2/h$, with the Pauli
matrix $\sigma_z$. Codes that define the spin current with
$S_z = (\hbar/2)\sigma_z$ report half this value.

> **Known limitation: the spin Hall curve of a `csr` model needs a second step.**
> With the formula `shc` in `reconstruction.formulas`, `python -m lsquant run`
> writes the moments of a `csr` model and then stops with exit status 1. We
> therefore leave `shc` out of the run file, give `num_cells` and `cell_area`,
> and call `python -m lsquant reconstruct --formula shc` on the moment file, as
> above. SIESTA models do not need the second step.
<!-- source: docs/explanation/kubo-formulas.md, known limitation "shc in csr and wannier runs"; python/lsquant/runfile.py:1714, 2239, 2318 (the three run paths call the formula without num_cells and cell_area); runfile.py:2431-2440 (reconstruct reads num_cells and cell_area from run.resolved.json) -->

## The result

![Five observables of the Kane–Mele model with Rashba coupling and exchange, KPM against the exact k-mesh result](fig_words.png)

**FIG. 1.** Five observables of the Kane–Mele model with Rashba coupling
($\lambda_{\mathrm R} = 0.05$ eV) and exchange ($M = 0.05$ eV), on a
$16 \times 16$ supercell. (a) Density of states per cell. (b) Spin density of
states per cell. (c) Longitudinal conductivity $\sigma_{xx}$. (d) Hall
conductivity $\sigma_{xy}$. (e) Spin Hall conductivity $\sigma^z_{xy}$ (Pauli
$\sigma_z$). Solid lines: KPM with 512 moments, exact trace and 20 meV broadening.
Dashed lines: exact diagonalisation on the same $16 \times 16$ k mesh, with a
20 meV Gaussian. The grey band marks the gap, from −0.510 to +0.412 eV.

The exact calculation gives a global gap of 0.922 eV, from −0.510 to
+0.412 eV, and a spin Chern number of −2 ($Z_2 = 1$). The run prints the
following comparison at $E_F = 0$ and over the range −3.5 to 3.5 eV:

| Observable | KPM at $E_F = 0$ | Exact at $E_F = 0$ | Largest difference over the range |
|---|---|---|---|
| $\rho(E)$, eV$^{-1}$ per cell | 0.0000 | 0.0000 | 0.37 (peak height 4.0) |
| $\rho_z(E)$, eV$^{-1}$ per cell | 0.0000 | 0.0000 | 0.33 (peak height 1.8) |
| $\sigma_{xx}$, $e^2/h$ | 0.0000 | 0.0000 | 92 (peak height 320) |
| $\sigma_{xy}$, $e^2/h$ | +0.0001 | +0.0001 | 0.16 |
| $\sigma^z_{xy}$, $e^2/h$ | −2.021 | −2.018 | 0.047 |

The two Hall curves agree with the exact result over the whole spectrum. In the
gap, $\sigma_{xy}$ vanishes and $\sigma^z_{xy}$ forms a plateau at
−2.02 $e^2/h$. Without the Rashba term the plateau would be exactly
$-2\,e^2/h$; the spin mixing shifts it by 1%. In the bands the exchange field
gives an anomalous Hall conductivity that reaches −0.8 $e^2/h$ near +1 eV, and it
splits the two spins in $\rho_z(E)$.

The density of states integrates to 4.000 states per cell in both calculations,
one per state. The spin density of states integrates to zero in both. The peak
heights of $\rho$, $\rho_z$ and $\sigma_{xx}$ differ for two reasons: the KPM
curves have the resolution of 512 moments, about 25 meV, against the 20 meV of
the exact calculation, and the Jackson kernel and the Gaussian have different
shapes. A
$16 \times 16$ mesh has sharp peaks, and the shape difference shows most there.
The clean $\sigma_{xx}$ depends on the square of the broadened delta function and
grows as the broadening decreases. Its integral over the range is 297 $e^2/h\cdot$eV
in KPM. The exact calculation gives 313 $e^2/h\cdot$eV with a 20 meV Gaussian
and 251 $e^2/h\cdot$eV with a 25 meV Gaussian, so the KPM value lies between the
two.

> **Known limitation: Hall curves of a clean metal are noisy with few random vectors.**
> On a $64 \times 64$ supercell with 4 random vectors instead of an exact trace,
> the plateau stays right: −2.014 and −2.031 $e^2/h$ for two seeds, against
> −2.018 $e^2/h$ exact. In the bands the two seeds differ by several $e^2/h$, for
> example +3.09 and +6.54 $e^2/h$ at −1.5 eV, against −0.40 $e^2/h$ exact. In this
> run the seeds differ most in the bands, where $\sigma_{xx}$ is large, and agree
> in the gap, where $\sigma_{xx}$ vanishes. For small clean cells we therefore use
> an exact trace; otherwise many random vectors, or disorder that lowers
> $\sigma_{xx}$.

## Runtime

On 8 threads of an AMD Ryzen Threadripper PRO 3995WX, the run of the five words
took 2 min 52 s, with a peak memory of 0.41 GB. The spin Hall reconstruction
took 4 s, the exact calculation 2 s and the model 1 s. The exact trace costs one
pass per orbital, so its time grows as the square of the number of orbitals. The
$64 \times 64$ run with 4 random vectors took 28 s.

## Output files

| File | Content |
|---|---|
| `out/*.chebmom` | Chebyshev moments of each word, for later reconstruction |
| `out/meanSpectralOp1km_L16*JACKSON.dat` | $\rho(E)$, states per eV for the whole supercell |
| `out/meanSpectralOpSZkm_L16*JACKSON.dat` | $\rho_z(E)$, same units |
| `out/KuboGreenWood_NonEqOpVX-VX*JACKSON.dat` | $\sigma_{xx}$ before the factor $2\pi/(N_{\mathrm{cells}}A_{\mathrm{cell}})$ |
| `out/KuboBastin_NonEqOpVX-VY*JACKSON.dat` | $\sigma_{xy}$ before the same factor |
| `out/SpinHall_NonEqOpJXSZ-VY*JACKSON.dat` | $\sigma^z_{xy}$ in $e^2/h$ |
| `out/run.resolved.json` | the resolved run file, with `num_cells` and `cell_area` |
| `operators/km_L16.oracle.json` | the exact reference |
| `fig_words.png` | FIG. 1 |

## References

- C. L. Kane and E. J. Mele, *Quantum spin Hall effect in graphene*,
  Phys. Rev. Lett. **95**, 226801 (2005).
- E. Prodan, *Robustness of the spin-Chern number*, Phys. Rev. B **80**, 125327 (2009).
- J. H. García, L. Covaci and T. G. Rappoport, *Real-space calculation of the
  conductivity tensor for disordered topological matter*, Phys. Rev. Lett. **114**,
  116602 (2015).

## Check

After the run, spin Hall reconstruction and oracle calculation above, check all
five saved curves without generating figures. The script states each tolerance,
prints CHECK PASS or CHECK FAIL and exits with status 0 or 1 respectively:

```bash
python3 examples/19_observable_words/check.py examples/19_observable_words/out --oracle examples/19_observable_words/operators/km_L16.oracle.json
```

An additional 25 meV oracle can be supplied with `--oracle-25` to check the
conductivity-integral bracket as well.
