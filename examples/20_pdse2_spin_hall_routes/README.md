# Example 20: Spin Hall plateau of pentagonal $\mathrm{PdSe_2}$ from SIESTA, Wannier90 and PAOFLOW

## The physical question

A pentagonal $\mathrm{PdSe_2}$ monolayer with spin–orbit coupling is a
semiconductor. Its gap at the Fermi level $E_F$ is trivial ($Z_2 = 0$), and
$\sigma^z_{xy}$ there is close to zero. About 1 eV above $E_F$ the conduction
bands open a second, narrower gap with $Z_2 = 1$, a quantum spin Hall gap. With
the Fermi level moved into that gap by doping or gating, $\sigma^z_{xy}$ forms a
plateau. We ask whether the value of that plateau depends on how the
Hamiltonian is obtained. We compute it from two density-functional
calculations and four routes from the electronic structure to $\sigma^z_{xy}$:

1. SIESTA Hamiltonian and overlap, all states, exact trace;
2. the same, with the states restricted to $[-3, 3]$ eV, exact trace and KPM;
3. Quantum ESPRESSO bands, a 12-orbital Wannier90 model, LSQuant;
4. the same Quantum ESPRESSO bands, the spin Hall module of PAOFLOW 3.0.

Routes 1 and 2 use the SIESTA files of
[example 8](../08_spin_hall_pdse2_monolayer/README.md). Route 4 is an
independent code and serves as the external check of route 3.

All conductivities are in $e^2/h$ with the Pauli matrix $\sigma_z$, as in
[Kubo formulas](../../docs/explanation/kubo-formulas.md); divide by 2 for
units of $\hbar/2$.

## The electronic structures

| | SIESTA | Quantum ESPRESSO 7.4 |
|---|---|---|
| Functional | PBE, fully relativistic | PBE, fully relativistic |
| Pseudopotentials | PSML (`siesta/` of example 8) | PAW, pslibrary 1.0.0: `Pd.rel-pbe-n-kjpaw_psl.1.0.0.UPF`, `Se.rel-pbe-dn-kjpaw_psl.1.0.0.UPF` |
| Basis | DZP numerical orbitals, 220 spinor orbitals per cell | plane waves, 80 Ry (800 Ry for the charge density) |
| Self-consistent k mesh | $12 \times 12 \times 1$ | $4 \times 4 \times 1$ |
| Cell | $a = 5.750$ Å, $b = 5.920$ Å, 16.45 Å along $z$ | same |

The geometry is that of
[wannier2sparse](https://github.com/adamecius/wannier2sparse) example 05. In
the Quantum ESPRESSO inputs the slab is centred at $z = 0$: a slab at $z = c/2$ puts the Wannier
functions across the boundary of the cell along $z$, and Wannier90 then returns
spreads of about 100 Å$^2$.
<!-- source: qe/scf.in; ../08_spin_hall_pdse2_monolayer/siesta/pdse2.fdf -->

The Wannier90 model has 12 spinor Wannier functions, Pd $d_{x^2-y^2}$ and Se
$p_z$. It is built from a non-self-consistent run on a $12 \times 12 \times 1$
mesh with 96 bands, an outer window from −3.5 to 3.0 eV and a frozen window from
−2.8 to 1.0 eV (absolute energies; $E_F = -1.316$ eV), with 1000
disentanglement and 400 localisation iterations. The final spreads are
12.6 Å$^2$ for the Pd functions, centred on the Pd atoms, and 10.9 Å$^2$ for
the Se functions, centred near the middle of the Se–Se dimers. Kramers
pairs of the model stay degenerate to 1.4 meV.
<!-- source: qe/pdse2_proj.win; Wannier90 .wout of the run -->

## The routes

**SIESTA, all states.** [`run_siesta_direct.json`](run_siesta_direct.json)
diagonalises $H(\mathbf k)$ and $S(\mathbf k)$ on a $24 \times 24$ mesh and
evaluates the word `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` with an exact trace,
the `covariant` velocity and the `lowdin` metric. With `oracle = true` it also
writes the exact occupied-band sum, `exact_shc.json`, with a 6 meV Gaussian.
The rows marked "exact" in the table and in FIG. 1 are these occupied-band sums. The run also writes a KPM
curve with 64 moments, which we do not use; the moment count is kept low
because this run serves only the exact sum.

**SIESTA, $[-3, 3]$ eV.** [`run_siesta_projected.json`](run_siesta_projected.json)
adds `"projection": {"kind": "shift", "window": [-3, 3]}`: the states outside the
window are moved to its nearest edge, so the expansion covers 6 eV instead of
the full spectrum. It writes the KPM curve with 1024 moments (Jackson kernel)
and, as `word_0_oracle.json`, the exact occupied-band sum of the same projected
states with a 6 meV Gaussian.

**Wannier90 and LSQuant.** [`wannier_route.py`](wannier_route.py) reads the
Wannier90 output and computes $\sigma^z_{xy}(E)$ on a $300 \times 300$ mesh with a
6 meV Gaussian. It uses two operators that the run file of a `wannier` model
does not provide:

- the covariant velocity, `WannierModel.velocity_blocks(k, d, mode="covariant")`,
  which needs `<seed>_r.dat` (see the
  [Python API reference](../../docs/reference/python-api.md));
- the spin operator projected from the DFT bands,
  $S_W(\mathbf k) = V^\dagger S V$ with $V = U_{\mathrm{dis}}U$, from
  `<seed>.spn`, `<seed>_u.mat`, `<seed>_u_dis.mat` and `<seed>.eig`.

With spin–orbit coupling the Wannier functions are not spin eigenstates: the
diagonal elements of $S_W$ average 0.58 and 0.75 for the two kinds of function,
not ±1. A diagonal $\sigma_z$ on the spinor labels is therefore a different
operator, and it gives a different plateau (see the result).

**PAOFLOW.** [`paoflow/run_pao.py`](paoflow/run_pao.py) projects the Quantum
ESPRESSO bands on the pseudo-atomic orbitals (92 bands above a projectability of
0.90), interpolates on a $120 \times 120$ mesh and calls
`spin_Hall(twoD=True, delta=0.006)` after `adaptive_smearing()`. PAOFLOW gives
$\sigma^z_{xy}$ in units of $e^2/\hbar$ with the spin operator $\sigma_z/2$; the
factor $4\pi$, $2\pi$ from $e^2/\hbar$ to $e^2/h$ and 2 from $\sigma_z/2$ to
$\sigma_z$, converts it to $e^2/h$ with the Pauli $\sigma_z$.
<!-- source: PAOFLOW 3.0.0, PAOFLOW.py:1476 (spin operator sP = 0.5 * Pauli);
response/do_Hall.py:77 (spin current 0.5*(Sj dHk + dHk Sj)), :447 (Omega_n =
-2 Im sum j_nm v_mn / dE^2, energies in eV, so Omega in bohr^2), :486 (sum over
k divided by N_k), :96 (twoD: divided by the cell area in bohr^2, no other
factor). The result is (1/(N_k A)) sum f Omega, i.e. sigma in units of e^2/hbar,
since int d^2k/(2pi)^2 = 1/(N_k A). LSQuant: sigma = 2pi/(N_k A) sum f Omega with
Pauli sigma_z, in e^2/h (docs/explanation/kubo-formulas.md). Ratio: 2pi * 2 = 4pi.
The 3D branch, :98, applies e^2/hbar in SI and 1/omega, the same convention. -->

## The run

The density-functional steps are the slow part. The figure and the table below
are made from the curves in [`data/`](data/), so Step 4 runs without them.

**Step 1: Quantum ESPRESSO and Wannier90.** We put the two pseudopotentials in
`pseudo/`, next to `qe/`, as `Pd.UPF` and `Se.UPF`, and run in `qe/`:

```bash
mpirun -np 32 pw.x -nk 8 -in scf.in  > scf.out
mpirun -np 32 pw.x -nk 8 -in nscf.in > nscf.out
wannier90.x -pp pdse2_proj
mpirun -np 16 pw2wannier90.x -in pw2wan_spn.in > pw2wan_spn.out
wannier90.x pdse2_proj
```

**Step 2: the LSQuant routes.** We run from the repository root:

```bash
python examples/20_pdse2_spin_hall_routes/wannier_route.py \
    examples/20_pdse2_spin_hall_routes/qe pdse2_proj -1.3162 wannier_covariant_spn \
    --compare-operators
python -m lsquant run examples/20_pdse2_spin_hall_routes/run_siesta_direct.json
python -m lsquant run examples/20_pdse2_spin_hall_routes/run_siesta_projected.json
```

`wannier_route.py` writes `<seed>.uc` and `<seed>.xyz` (cell and Wannier
centres) into the Wannier90 directory, then `wannier_covariant_spn.dat` and a
`.json` with the gaps, $Z_2$ and the operator comparison.

**Step 3: PAOFLOW.** We copy the directory `tmp/` written by the
self-consistent run of Step 1 to `paoflow/tmp/` and run in `paoflow/`:

```bash
mpirun -np 32 pw.x -nk 8 -in nscf_pao.in > nscf_pao.out
projwfc.x -in projwfc.in > projwfc.out
python run_pao.py
```

**Step 4: the figure.** We plot the curves of `data/`:

```bash
python examples/20_pdse2_spin_hall_routes/plot.py
```

## The result

![Spin Hall conductivity of pentagonal PdSe2 from four routes](fig_routes.png)

**FIG. 1.** $\sigma^z_{xy}(E)$ of pentagonal $\mathrm{PdSe_2}$, in $e^2/h$ with
Pauli $\sigma_z$. (a) From $E_F - 0.5$ eV to $E_F + 2$ eV. (b) The two
quantum spin Hall gaps. Black solid line: SIESTA, all states, exact
occupied-band sum on a $24 \times 24$ mesh with a 6 meV Gaussian. Blue dashed
line: SIESTA, states in $[-3, 3]$ eV, exact sum, same mesh and Gaussian.
Orange dotted line: the same projected states, KPM with $M = 1024$ moments, the
Jackson kernel and an exact trace. Green dash-dotted line: Quantum ESPRESSO
bands, Wannier90 model, covariant velocity and spin projected from DFT,
$300 \times 300$ mesh, 6 meV Gaussian. Pink dashed line: the same Quantum
ESPRESSO bands, PAOFLOW 3.0, $120 \times 120$ mesh, adaptive smearing. Dark
grey band: gap of the SIESTA electronic structure; light grey band: gap of the
Quantum ESPRESSO one, as reproduced by the Wannier90 model.

| Route | Gap with $Z_2 = 1$, $E - E_F$ (eV) | Plateau ($e^2/h$) | At $E_F$ ($e^2/h$) |
|---|---|---|---|
| SIESTA, all states, exact | 1.014 – 1.090 | 1.979 | +0.011 |
| SIESTA, $[-3, 3]$ eV, exact | 1.014 – 1.090 | 1.968 | −0.010 |
| SIESTA, $[-3, 3]$ eV, KPM, 1024 moments | 1.014 – 1.090 | 2.012 | −0.010 |
| Quantum ESPRESSO, Wannier90 + LSQuant | 1.188 – 1.233 | 1.944 | −0.066 |
| Quantum ESPRESSO, PAOFLOW 3.0 | flat from 1.196 to 1.246 | 1.942 | +0.003 |

The plateau is the mean over the central half of the gap. The gaps of the two
SIESTA routes and of the Wannier90 model are global gaps on their k meshes; for
PAOFLOW the table gives the range where its curve stays within 0.02 $e^2/h$ of
the plateau. $Z_2$ comes from the spin Chern number of the occupied bands, with
the spin operator of each model: $Z_2 = 1$ in the upper gap and $Z_2 = 0$ at
$E_F$, for SIESTA on a $24 \times 24$ mesh and for the Wannier90 model on a
$60 \times 60$ mesh.
<!-- source: plot.py output; data/wannier_covariant_spn.json -->

The five plateaus lie between 1.94 and 2.01 $e^2/h$. With one electronic
structure, the routes agree closely: the two exact SIESTA routes differ by 0.6%,
and LSQuant and PAOFLOW on the Quantum ESPRESSO bands differ by 0.1%. The KPM
plateau lies 2% above the exact result for the same projected states. The two
electronic structures place the gap 0.16 eV apart and give plateaus 1–2%
apart; they come from different codes, pseudopotentials and basis sets, and
this example does not separate those effects.

Outside the gaps the curves differ. The band structures of the two codes are
not the same, and the four routes use different k meshes and broadenings
($24 \times 24$ with a 6 meV Gaussian, $300 \times 300$ with a 6 meV Gaussian,
$120 \times 120$ with adaptive smearing in PAOFLOW). The Wannier90 model
contains only the 12 bands of its windows and reproduces the DFT bands inside
its frozen window, −1.5 to +2.3 eV around $E_F$; the Fermi-sea sum of the other
routes also runs over all lower bands. We have not converged the curves in the
bands, and the example compares only the plateaus.

**The spin operator and the Wannier90 plateau.** For the Wannier90 model,
`--compare-operators` evaluates the centre of the upper gap with three
velocities and two spin operators:

| Velocity | Diagonal $\sigma_z$ | Spin projected from DFT |
|---|---|---|
| `bare` | 1.583 | 1.824 |
| `berry_connection` | 1.507 | 1.947 |
| `covariant` | 1.577 | 1.946 |

With the projected spin, the covariant and Berry-connection velocities agree to
0.1%, and both agree with PAOFLOW; the bare velocity is 6% lower. With the
diagonal $\sigma_z$, the three velocities give 1.51 to 1.58 $e^2/h$; for the
covariant velocity, 1.577 against 1.946 $e^2/h$ is 19% lower.
<!-- source: data/wannier_covariant_spn.json, key operators_at_upper_gap -->

The same inputs with the $4 \times 4$ non-self-consistent mesh and 10 + 10
iterations of [wannier2sparse](https://github.com/adamecius/wannier2sparse)
example 05 give a Wannier90 model without a global gap near 1.2 eV (bands overlap by 13 meV) and with Kramers splittings up to
24 meV. The denser mesh and the longer minimisation above remove both.

## Runtime

On 16 threads of an AMD EPYC 9374F:

| Step | Wall time | Peak memory |
|---|---|---|
| Quantum ESPRESSO, self-consistent and $12 \times 12$ non-self-consistent (32 MPI processes) | 3 min + 12 min | — |
| `wannier_route.py`, $300 \times 300$ mesh | 73 s | 2.8 GB |
| SIESTA, all states, exact trace, 64 moments | 46 min | 10.8 GB |
| SIESTA, $[-3, 3]$ eV, exact trace, 1024 moments | 3 h 10 min | 23.3 GB |
| PAOFLOW | 8 min (CPU time) | — |

The exact trace costs one pass per state of the $24 \times 24$ mesh, 126 720
states, and the projected run does it for 1024 moments.

## Files

| File | Content |
|---|---|
| `qe/scf.in`, `qe/nscf.in`, `qe/pw2wan_spn.in`, `qe/pdse2_proj.win` | Quantum ESPRESSO and Wannier90 inputs |
| `paoflow/nscf_pao.in`, `paoflow/projwfc.in`, `paoflow/run_pao.py` | PAOFLOW inputs |
| `run_siesta_direct.json`, `run_siesta_projected.json` | LSQuant run files on the SIESTA files of example 8 |
| `wannier_route.py` | Wannier90 route: covariant velocity and projected spin |
| `data/*.dat` | $\sigma^z_{xy}(E - E_F)$ of each route, in $e^2/h$ (Pauli $\sigma_z$) |
| `data/wannier_covariant_spn.json` | gaps, $Z_2$ and operator comparison of the Wannier90 model |
| `plot.py` | FIG. 1 and the table |
| `fig_routes.png` | FIG. 1 |

## References

- A. A. Mostofi *et al.*, *An updated version of wannier90*, Comput. Phys.
  Commun. **185**, 2309 (2014).
- G. Pizzi *et al.*, *Wannier90 as a community code*, J. Phys.: Condens. Matter
  **32**, 165902 (2020).
- M. Buongiorno Nardelli *et al.*, *PAOFLOW: A utility to construct and operate
  on ab initio Hamiltonians from the projections of electronic wavefunctions on
  atomic orbital bases*, Comput. Mater. Sci. **143**, 462 (2018).
- E. Prodan, *Robustness of the spin-Chern number*, Phys. Rev. B **80**, 125327 (2009).
- J. H. García, L. Covaci and T. G. Rappoport, *Real-space calculation of the
  conductivity tensor for disordered topological matter*, Phys. Rev. Lett. **114**,
  116602 (2015).
