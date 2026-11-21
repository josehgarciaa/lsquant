# Example 15: layer projectors and the mean-square displacement

## The physical question

In a layered or nanostructured system we often need to know where the states
at a given energy live, and how fast electrons at that energy spread. LSQuant
answers both with operators built from the orbital positions:

- a **projector** $P_\ell$ on a region, such as a layer, a ribbon edge or an
  impurity site, turns the DOS into a local DOS,
  $\rho_\ell(E) = \frac1N\mathrm{Tr}[P_\ell\,\delta(E-H)]$;
- the **position operator** $X$ enters through the commutator $[X, U(t)]$ with
  the time-evolution operator $U(t) = e^{-iHt/\hbar}$, which gives the
  energy-resolved mean-square displacement

$$ \rho(E)\,\Delta X^2(E,t) = \frac1N\,\mathrm{Tr}\big[[X,U(t)]^\dagger\,
\delta(E-H)\,[X,U(t)]\big] . $$

$X$ itself is unbounded on a periodic system, and LSQuant never builds it.
LSQuant evaluates the commutator from the velocity
$V_X = \frac{i}{\hbar}[H, X]$, and we build this velocity from the positions.

Our system is a biased bilayer chain, a one-dimensional model of a gated
bilayer. We ask how the bias moves weight between the layers, and we check
that electrons move ballistically with the band velocity.

## The model

[`make_bilayer_chain.py`](make_bilayer_chain.py) writes two chains of $L$
sites, the top and the bottom layer, with lattice constant $a = 1$ Å. The
intra-chain hopping is $-t$, the inter-layer hopping $-t_\perp$, and the layer
bias $\pm\Delta/2$:

$$ H = -t\sum_{\ell,i}\big(|\ell,i\rangle\langle \ell,i+1| + \mathrm{h.c.}\big)
- t_\perp\sum_i\big(|\mathrm{T},i\rangle\langle \mathrm{B},i| + \mathrm{h.c.}\big)
+ \frac\Delta2\sum_i\big(|\mathrm{T},i\rangle\langle \mathrm{T},i| -
|\mathrm{B},i\rangle\langle \mathrm{B},i|\big), $$

with $t = 1$ eV, $t_\perp = 0.5$ eV, $\Delta = 0.6$ eV, periodic along $x$. The
bands are $E_\pm(k) = -2t\cos ka \pm s$ with $s = \sqrt{t_\perp^2 +
\Delta^2/4}$, and both have the velocity $v(k) = 2ta\sin ka$.

Next to the Hamiltonian, the script writes the layer projectors as diagonal
CSR matrices, with 1 on the orbitals of one layer and 0 elsewhere. It also
writes the velocity $V_X = i[H, X]$ from the orbital positions $x_i$, taking
the minimum image across the periodic boundary:
$(V_X)_{ij} = i\,H_{ij}\,(x_j - x_i)$. We generate the small system for the
layer-resolved DOS and the long one for the mean-square displacement:

```bash
python3 examples/15_spatial_operators/make_bilayer_chain.py 256     # layer-resolved DOS
python3 examples/15_spatial_operators/make_bilayer_chain.py 4096    # mean-square displacement
```

Each call writes `operators/bilayer_L<L>.{HAM,VX,A,B}.CSR` and the positions
`bilayer_L<L>.xyz`. `A` is the projector on the top layer and `B` the projector
on the bottom layer.

## Run 1: the layer-resolved density of states

[`run_layer_dos.json`](run_layer_dos.json) computes three words on the
$L = 256$ bilayer, which has 512 orbitals, with an exact trace and 512
moments:

```json
"words": ["Tr[delta(E-H)]", "Tr[A delta(E-H)]", "Tr[B delta(E-H)]"]
```

> **Known limitation: operator names.** A projector name outside the word
> catalogue, such as `PTOP`, passes validation and `lsquant inspect`, but
> `lsquant run` fails with `Unknown operator`. The cause is that the per-word
> compile does not receive the model's operators. Catalogue names such as `A`,
> `B`, `O`, `P`, `PX`, `SX` and `VX` work.
> <!-- source: python/lsquant/runfile.py:1546; words/compiler.py:28, 393-394, 405 -->

The run is started with:

```bash
python3 -m lsquant run examples/15_spatial_operators/run_layer_dos.json --out /tmp/ex15/layer_dos
```

It writes one curve per word: the total DOS in
`meanSpectralOp1bilayer_L256…JACKSON.dat`, the top layer in
`meanSpectralOpAbilayer_L256…JACKSON.dat` and the bottom layer in
`meanSpectralOpBbilayer_L256…JACKSON.dat`.

The bias pushes the top layer up in energy and the bottom layer down. The
bonding band $E_-$ lives mostly on the bottom layer, and the antibonding band
$E_+$ mostly on the top layer. The first moment makes this exact. The mean
energy of the states in layer $\ell$ is $\mathrm{Tr}[P_\ell H]/L = \pm\Delta/2$,
and the layer DOS give $\langle E\rangle_\mathrm{top} = +0.3000$ eV and
$\langle E\rangle_\mathrm{bottom} = -0.3000$ eV.

## Run 2: the mean-square displacement

[`run_msd.json`](run_msd.json) uses the $L = 4096$ bilayer with 8192 orbitals.
The chain is long enough that the electrons stay clear of the boundary for
40 fs. The run uses a random-phase trace with 2 vectors, 256 moments, and 41
time steps up to 40 fs:

```json
"words": ["Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]", "Tr[delta(E-H)]"],
"time": {"steps": 41, "tmax": 40.0}
```

LSQuant uses $\hbar = 0.6582119624$ eV fs, so with $H$ in eV the times
`time.dt` and `time.tmax` are in fs.
<!-- source: python/lsquant/moments.py:22 -->
The reconstruction evaluates the mean-square displacement at two Fermi
energies:

```json
"reconstruction": {"kernel": "jackson", "broadening_meV": 60.0,
                   "fermi_energies": [0.0, 1.0], "formulas": ["time_correlation", "dos"]}
```

The run is started with:

```bash
python3 -m lsquant run examples/15_spatial_operators/run_msd.json --out /tmp/ex15/msd
```

It writes `meanCorrelationVXbilayer_L4096KPM_M256_state40defaultEF<E_F>JACKSON.dat`
and the DOS from the same random vectors. The first file has two columns, $t$
in fs and $\rho(E_F)\,\Delta X^2(E_F,t)$ per orbital in Å²/eV. Dividing it by
the DOS per orbital gives $\Delta X^2(E_F,t)$ in Å².

The clean chain has no scattering, so the spreading is ballistic,
$\Delta X^2 = \langle v^2\rangle_E\,t^2/\hbar^2$. At energy $E$ both bands
contribute, with $|v_\pm| = \sqrt{4t^2 - (E\mp s)^2}$. Each band contributes
with its DOS $\propto 1/|v_\pm|$, so that
$\langle v^2\rangle_E = \sum_\pm |v_\pm| / \sum_\pm |v_\pm|^{-1}$. At
$E_F = 0$ electrons spread to $\Delta X^2 = 13\,500$ Å² in 40 fs, about 116 Å.

## Checks

[`check_spatial.py`](check_spatial.py) reads both runs. We run it with:

```bash
python3 examples/15_spatial_operators/check_spatial.py /tmp/ex15
```

Output on the reference run:

```text
[1] max |rho_A + rho_B - rho| / max rho = 2.92e-06  PASS
[2,3] layer A (top): weight = 256.00 (L = 256), <E> = +0.3000 eV (expected +0.300)  PASS
[2,3] layer B (bottom): weight = 256.00 (L = 256), <E> = -0.3000 eV (expected -0.300)  PASS
[4,5] E_F = 0.0 eV: slope = 2.00000, DeltaX^2/t^2 = 8.435 A^2/fs^2, <v^2>/hbar^2 = 8.448; DeltaX^2(t=40 fs) = 13495 A^2  PASS
[4,5] E_F = 1.0 eV: slope = 2.00000, DeltaX^2/t^2 = 5.572 A^2/fs^2, <v^2>/hbar^2 = 5.518; DeltaX^2(t=40 fs) = 8916 A^2  PASS
CHECK PASS
```

1. The two projectors resolve the identity, so the layer DOS add up to the
   total DOS to the 6 significant digits of the curve files.
2. Each layer holds $L = 256$ states.
3. Lines [2,3] also print the mean energy of each layer, $+0.3000$ eV for the
   top layer and $-0.3000$ eV for the bottom one. These equal the exact first
   moments $\mathrm{Tr}[P_\ell H]/L = \pm\Delta/2$ of
   [run 1](#run-1-the-layer-resolved-density-of-states): the layer DOS place
   each layer at its bias.
4. $\Delta X^2(t)$ grows exactly as $t^2$, with a log–log slope of 2.00000.
5. $\Delta X^2/t^2$ matches $\langle v^2\rangle_E/\hbar^2$ to 0.2 % at
   $E_F = 0$ and to 1 % at $E_F = 1$ eV. The difference comes from the 60 meV
   broadening and the 2-vector stochastic trace.

See also [running on a cluster](../../docs/how-to/run-on-a-cluster.md) for the
atom-block spatial layout of the native build.

## Runtime

On a shared 96-core node with one OpenMP thread, the layer DOS took about
19 s, the mean-square displacement run 9 s and the check under 1 s.
