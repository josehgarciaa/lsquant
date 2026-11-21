# Example 14: the SSH chain from a hand-written Wannier90 `_hr.dat`

## The physical question

The Su–Schrieffer–Heeger (SSH) chain is the simplest crystal with a band gap
that comes from the lattice alone. Each cell holds two orbitals, A and B. An
electron hops with $t_1$ inside a cell (A–B) and with $t_2$ between cells (B to
the next A). In the cell-periodic gauge the Bloch Hamiltonian is

$$ H(k) = \begin{pmatrix} 0 & t_1 + t_2 e^{-ika} \\ t_1 + t_2 e^{ika} & 0 \end{pmatrix} , $$

with the two bands

$$ E_\pm(k) = \pm\,|t_1 + t_2 e^{ika}| . $$

The bands reach $\pm|t_1 + t_2|$ at $k = 0$ and $\pm|t_1 - t_2|$ at
$k = \pi/a$. When $t_1 \neq t_2$ a gap of width $2|t_1 - t_2|$ opens around
$E = 0$, and the DOS has a 1D van Hove divergence at each of the four band
edges. Inside a band the DOS per orbital follows from
$\rho(E) = (1/2\pi)\,|dk/dE|$:

$$ \rho(E) = \frac{|E|}{2\pi\,|t_1 t_2|\,|\sin ka|}, \qquad
\cos ka = \frac{E^2 - t_1^2 - t_2^2}{2\,t_1 t_2} . $$

We write this model by hand in the Wannier90 format, the format that most
DFT-derived tight-binding models use, and let LSQuant read it with
`model.kind = "wannier"`. We then check that LSQuant reproduces the analytic
bands and that the DOS vanishes in the gap.

LSQuant reads the Wannier90 files itself. It assembles $H(\mathbf k)$,
$\mathbf v(\mathbf k)$ or the real-space supercell CSR in Python, in the
`real`, `bloch` or `eigen` basis, so the `wannier2sparse` converter is
optional. When `model.csr_dir` points at its output, the `real` basis reads
those CSR files instead and checks that they match the model.

## The model

We take $t_1 = -1.0$ eV, $t_2 = -0.6$ eV and $a = 3$ Å. The bands then span
$0.4 \le |E| \le 1.6$ eV, with a gap $2|t_1 - t_2| = 0.8$ eV.

A Wannier model is a directory with three files that share a seed name (here
`ssh`). The run file gives them as `model.dir` and `model.seed`.

[`ssh/ssh_hr.dat`](ssh/ssh_hr.dat) holds the hoppings
$H_{ij}(\mathbf R) = \langle i,\mathbf 0|H|j,\mathbf R\rangle$ in eV, in the
standard Wannier90 layout. The file starts with a comment line, the number of
orbitals $N_\mathrm{orb}$ and the number of lattice vectors $n_R$. The $n_R$
degeneracy weights follow, 15 per line in files written by Wannier90. Each
remaining line, `R1 R2 R3 i j Re Im`, gives one matrix element, with orbital
indices counted from 1. For the SSH chain, $t_1$ appears in the
$\mathbf R = 0$ block and $t_2$ in the $\mathbf R = \pm 1$ blocks:

```text
SSH chain, t1 = -1.0 eV (intra-cell A-B), t2 = -0.6 eV (inter-cell B-A1), hand-written
2
3
1 1 1
   -1    0    0    1    1    0.000000    0.000000
   -1    0    0    2    1    0.000000    0.000000
   -1    0    0    1    2   -0.600000    0.000000
   ...
    0    0    0    2    1   -1.000000    0.000000
    0    0    0    1    2   -1.000000    0.000000
   ...
    1    0    0    2    1   -0.600000    0.000000
   ...
```

LSQuant divides each $H(\mathbf R)$ by its degeneracy weight, accepts files that
list only the non-zero elements, and reads `<seed>_wsvec.dat` if it is present.

[`ssh/ssh.uc`](ssh/ssh.uc) holds the three lattice vectors (rows, Å). The chain
runs along $x$ with $a = 3$ Å. The other two vectors only need to be long.

[`ssh/ssh.xyz`](ssh/ssh.xyz) holds the orbital count and one `label x y z` line
per orbital (Å), in the order of the orbital indices of `_hr.dat`. The centres
enter the Bloch phases (atom gauge) and the velocity

$$ v_x(\mathbf k) \propto \sum_{\mathbf R}(\mathbf R + \tau_j - \tau_i)_x\,
e^{i\mathbf k\cdot(\mathbf R+\tau_j-\tau_i)}\,H_{ij}(\mathbf R) , $$

so they must be the Wannier centres reported by Wannier90
(`<seed>_centres.xyz`), not the atomic positions, if the two differ.

For a real Wannier90 calculation we copy `<seed>_hr.dat`, write `<seed>.uc`
from the `unit_cell_cart` block, and write `<seed>.xyz` from the Wannier
centres.

## The run

[`run_ssh_dos.json`](run_ssh_dos.json) computes the DOS and $\sigma_{xx}$ in the
`bloch` basis on a 256-point k mesh, set by `basis.supercell = [256, 1, 1]`. It
uses $M = 256$ moments, an exact trace, a 20 meV Jackson broadening and the
spectral window from the band extrema (`"eigen"`). The run is started with:

```bash
python3 -m lsquant run examples/14_wannier/run_ssh_dos.json --out /tmp/ex14/ssh_dos
```

The `"eigen"` window is the band extrema $\pm 1.6$ eV padded by 10 % of the
bandwidth on each side, which gives $[-1.92, 1.92]$ eV. With this half-width,
a 20 meV broadening would need $\lceil \pi \cdot 1.92 / 0.020 \rceil = 302$
moments, more than the 256 we computed. The reconstruction then uses all 256,
and the actual resolution is $\pi \cdot 1.92/256 \approx 24$ meV.

The run directory holds the DOS moments and curve as
`SpectralOp1ssh_bloch_*`, the $\sigma_{xx}$ moments and curve as
`NonEqOpVX-VXssh_bloch_*`, and one `*.basis.json` per word. Each
`*.basis.json` records the basis, the mesh, the window and the origin of each
operator. Here the origin reads `"built by lsquant"`. For operators read from
`wannier2sparse` output, it is the file name.

Setting `basis.name` to `real` runs the same model as a 512-orbital
real-space chain. This is the starting point for adding disorder with
`basis.disorder = {"W": ..., "seed": ...}`.

## The result

The analytic DOS is zero in the gap $|E| < 0.4$ eV and has two bands with van
Hove peaks at $\pm 0.4$ and $\pm 1.6$ eV. Inside the bands, away from the
edges, its average over $0.6 < |E| < 1.4$ eV is

$$ \langle\rho\rangle = \frac{1}{0.8\ \mathrm{eV}}\int_{0.6}^{1.4}\rho(E)\,dE
= \frac{\arccos(-5/6) - \pi/3}{1.6\,\pi}\ \mathrm{eV}^{-1} = 0.300\ \mathrm{eV}^{-1} , $$

and the reconstructed DOS gives 0.300 1/eV on the same interval. In the gap,
the reconstructed DOS is $5\times10^{-4}$ of its value in the bands. The residue
is the tail of the broadened band edges at $\pm 0.4$ eV.

The bands depend only on $|t_1 + t_2 e^{ika}|$, which is unchanged when we swap
$t_1$ and $t_2$. The DOS of the periodic chain therefore cannot tell the
trivial phase, $|t_1| > |t_2|$ as here, from the topological one,
$|t_1| < |t_2|$. The difference shows up only at the ends of an open chain.

## Checks

[`check_ssh.py`](check_ssh.py) computes the exact DOS moments from the analytic
bands on the same 256-point mesh,

$$ \mu_m = \frac{2 - \delta_{m0}}{2N_k}\sum_{k,\pm} T_m\!\Big(\frac{E_\pm(k) - b}{a_w}\Big) , $$

with $a_w$ and $b$ the half-width and centre of the window stored in the moment
file. It compares them with the moments that LSQuant computed from
`ssh_hr.dat`, and it checks that the reconstructed DOS vanishes in the gap. We
run it with:

```bash
python3 examples/14_wannier/check_ssh.py /tmp/ex14/ssh_dos
```

Output on the reference run:

```text
moment file: SpectralOp1ssh_bloch_KPM_M256_stateexact.chebmom
window [-1.920000, 1.920000] eV; analytic band edges +/-0.400, +/-1.600 eV
max |mu_m - mu_m(analytic)| over 256 moments = 6.758e-15
DOS per orbital: max inside |E|<0.25 eV = 1.43e-04; mean inside 0.6<|E|<1.4 eV = 0.300 1/eV
CHECK PASS
```

1. The second line prints the window, $[-1.92, 1.92]$ eV: the band extrema
   $\pm 1.6$ eV padded by 10 % of the bandwidth, as computed in
   [the run](#the-run). It also prints the analytic band edges on the mesh,
   $\pm 0.4$ and $\pm 1.6$ eV. They equal $|t_1 \mp t_2|$, as predicted,
   because $k = \pi/a$ lies on the 256-point mesh.
2. The 256 moments computed from the hand-written `_hr.dat` agree with the
   analytic ones to $6.8\times10^{-15}$, which is round-off. LSQuant reads the
   Wannier90 file as the intended model.
3. The DOS per orbital is at most $1.4\times10^{-4}$ 1/eV in the gap, against
   0.30 1/eV in the bands.

## Runtime

On a shared 96-core node with a load of about 160 and four OpenMP threads,
the run took 25.5 s and the check 1.5 s.
