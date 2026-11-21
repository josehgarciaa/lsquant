# Example 9: one set of moments, many broadenings

## The physical question

A KPM calculation has two steps. The first evaluates the Chebyshev moments of
the rescaled Hamiltonian $\tilde H = (H - b)/a$, where $a$ is the half-width
and $b$ the centre of the spectral window. LSQuant stores the moments with the
Chebyshev weight already included:

$$ \mu_m = \frac{2 - \delta_{m0}}{N}\,\mathrm{Tr}\,T_m(\tilde H) . $$

The moments depend only on the Hamiltonian and on the window. This is the
expensive step. The second step turns the stored moments into a
curve:

$$ \rho(E) = \sum_{m=0}^{M_\mathrm{eff}-1} g_m\,\mu_m\,
\frac{T_m(\tilde E)}{\pi\sqrt{a^2 - (E-b)^2}} . $$

Here $g_m$ are the Jackson kernel coefficients. They turn each delta peak into
a peak of width about $\eta$. The reconstruction keeps only

$$ M_\mathrm{eff} = \min\!\big(\lceil \pi a/\eta \rceil,\; M\big) $$

moments, so the finest broadening that $M$ moments can give is $\pi a/M$.
<!-- source: python/lsquant/reconstruct.py:22-120 (dos_curve: no (2 - delta) factor at reconstruction); moments.py:25-32 (weights 1/N, 2/N folded in); docs/explanation/kpm.md "Chebyshev moments" -->
The broadening $\eta$ is therefore a resolution that we choose after the
moments exist, not a property of the system. We can change it as often as we
like at almost no cost.

We show this on the periodic 1D chain of
[tutorial 1](../01_dos_1d_chain/README.md), with hopping $t = -1$ eV. Its bulk
density of states is

$$ \rho(E) = \frac{1}{\pi\sqrt{4t^2 - E^2}} , $$

with van Hove divergences at the band edges $E = \pm 2|t| = \pm 2$ eV. A finite
broadening rounds each divergence into a peak. As $\eta$ shrinks, the peaks grow
and move towards $\pm 2$ eV, and the weight under the curve stays the same. We
watch this with two commands. `lsquant inspect` shows what a finished run
stored, and `lsquant reconstruct` rebuilds curves from the stored moments alone.

## The moments

We reuse the run of tutorial 1: chains of $N = 256, 1024, 4096$ sites, the
window $[-2, 2]$ eV and $M = 512$ exact-trace moments. The run is started with:

```bash
python3 -m lsquant run examples/01_dos_1d_chain/run_dos.json --out /tmp/ex09/run
```

This is the slow step, and the exact trace of the $N = 4096$ chain takes most
of its time. Every command below takes about a second.

## What the run stored

The header of a moment file holds the numbers that enter the formulas above:
the rank, the number of moments $M$, the system size $N$ and the window that
fixes $a$ and $b$. We read it with `inspect`, which always prints JSON:

```bash
python3 -m lsquant inspect /tmp/ex09/run/SpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256.chebmom1D
```

```json
{
  "kind": "moments",
  "rank": 1,
  "dims": [512],
  "system_size": 256,
  "number_of_states": 0,
  "state_label": "exact_chain1d_N256",
  "word": "",
  "window": [-2.0, 2.0],
  "band_width": 4.0,
  "half_width": 2.0,
  "band_center": 0.0,
  "kpm_weighted": true,
  "symmetrised": false,
  "number_of_time_steps": 0,
  "time_step": 0.0
}
```

The flag `kpm_weighted: true` says that the factor $(2 - \delta_{m0})/N$ is
already in the stored moments. The text moment files (`.chebmom1D`,
`.chebmom2D`, `.chebmomTD`) keep the format of the earlier codes. They do not
store the word or the number of random vectors, so these fields read `""` and
`0`. The binary `.chebmom` format, written when `output.write_chebmom` is true,
keeps them.

The same command also takes a run directory, where it prints the resolved run
file and the header of every moment file, and a run file, which it only
validates. The
[command-line reference](../../docs/reference/cli.md) describes these forms.

## Choosing the resolution

`reconstruct` takes a moment file, a formula and a broadening, as named options
or as positional arguments, `reconstruct <target> [formula] [broadening_meV]`.
We rebuild the $N = 256$ DOS at 60 meV with named options, and at 120 meV with
positional ones:

```bash
M=/tmp/ex09/run/SpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256.chebmom1D
python3 -m lsquant reconstruct $M --formula dos --broadening-meV 60 --out /tmp/ex09/eta60 --json
python3 -m lsquant reconstruct $M dos 120 --out /tmp/ex09/eta120
```

The second command prints only the path of the curve it wrote. The first, with
`--json`, prints a summary:

```json
{
  "moment_file": ".../SpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256.chebmom1D",
  "output_dir": "/tmp/ex09/eta60",
  "formula": "dos",
  "broadening_meV": 60.0,
  "fermi_energies": [0.0],
  "curve_files": ["/tmp/ex09/eta60/meanSpectralOp1chain1d_N256KPM_M512_stateexact_chain1d_N256JACKSON.dat"]
}
```

The curve file has two columns: $E$ in eV and the extensive DOS $N\rho(E)$ in
1/eV. Its $30 M_\mathrm{eff}$ points cover 95 % of the spectral window, that is
$|E| \le 1.9$ eV. With $a = 2$ eV, the rule for $M_\mathrm{eff}$ gives the
number of moments that each broadening uses:

| $\eta$ | $M_\mathrm{eff}$ | Where the curve comes from |
|---|---|---|
| 120 meV | 53 | `reconstruct ... dos 120`, and the curve that `run` wrote |
| 60 meV | 105 | `reconstruct ... dos 60` |
| 20 meV | 315 | `eta20/` subfolder written by `run` |

All three use only part of the 512 stored moments. A smaller $\eta$ reads
further into the same list.

Given a directory, `reconstruct` processes every moment file in it. It takes
the broadening from the command line, or from `run.resolved.json` when we omit
it. It writes the curves to `--out`, or next to the moments by default. We
rebuild all three chain lengths at 60 meV with:

```bash
python3 -m lsquant reconstruct /tmp/ex09/run -b 60 --out /tmp/ex09/dir_eta60
```

## Another formula from the same moments

For transport the formula is also a choice made after the moments exist, as
long as it matches the word that produced them. A rank-1 file such as
`Tr[delta(E-H)]` or `Tr[A delta(E-H)]` takes `dos`. A rank-2 Kubo file takes
`greenwood`, `greenwood_fft`, `bastin`, `sea`, `surf`, `bastin_fft` or `shc`. A
time-dependent file takes `time_correlation` with `--fermi-energies`. The
default `--formula` is `dos` for rank 1, `bastin` for rank 2 and
`time_correlation` for time-dependent moments. For a Kubo–Greenwood file we
therefore give the formula explicitly.

From the Hall moments of the Haldane model in
[example 13](../13_bastin_sea_surface/README.md), for instance, we rebuild only
the Fermi-sea term at 100 meV:

```bash
python3 -m lsquant reconstruct <out>/NonEqOpVX-VYhaldaneKPM_M128x128_stateexact.chebmom2D sea 100
```

It writes `KuboBastinSea_…JACKSON.dat` in about 3 s.

## The result

The van Hove peaks behave as the kernel predicts. On the $N = 256$ chain, the
largest value of $\rho(E)$ on the grid is 0.510 1/eV at 60 meV and 0.530 1/eV
at 20 meV. The bulk DOS grows towards the band edges, so it reaches its largest
value on the grid at the grid edge, $\rho(1.9\ \mathrm{eV}) = 1/(\pi\sqrt{4 -
1.9^2}) = 0.510$ 1/eV. The 60 meV curve agrees with it. At 20 meV the curve
uses $M_\mathrm{eff} = 315 > N$ moments. It starts to resolve the discrete
levels of the 256-site ring of tutorial 1, and its maximum rises above the bulk
value.

The weight is unchanged. On $|E| \le 1.9$ eV the $N = 256$ curve holds 0.7979
states per site at 20 meV and 0.7978 at 60 meV. The kernel redistributes
spectral weight but conserves the number of states.

## Checks

[`check_reconstruct.py`](check_reconstruct.py) reads the run directory and the
reconstructions. We run it with:

```bash
python3 examples/09_inspect_reconstruct/check_reconstruct.py /tmp/ex09/run /tmp/ex09
```

Output on the reference run:

```text
[1] reconstruct at 120 meV == curve written by run: True  PASS
[2] N=4096, 60 meV: weight on |E|<1.9 eV = 0.7978 per site (bulk 0.7978)  PASS
[3] N=256: max rho 20 meV = 0.530, 60 meV = 0.510 1/eV; weights 0.7979, 0.7978  PASS
CHECK PASS
```

1. Reconstructing the stored moments at 120 meV reproduces, byte for byte, the
   120 meV curve that `run` wrote. The second step is deterministic: the same
   moments and the same $\eta$ give the same curve.
2. The bulk chain DOS carries the weight

   $$ \int_{-1.9}^{1.9} \frac{dE}{\pi\sqrt{4 - E^2}} = \frac{2}{\pi}\arcsin(0.95) = 0.7978 $$

   on the reconstruction interval. The $N = 4096$ curve at 60 meV gives 0.7978:
   the long chain reaches the bulk limit.
3. Line [3] prints the largest DOS of the $N = 256$ curve, 0.530 1/eV at
   20 meV and 0.510 1/eV at 60 meV, and its weights, 0.7979 and 0.7978 states
   per site. A larger broadening lowers the maximum and keeps the weight.
   [The result](#the-result) discusses these numbers.

## Runtime

On a shared 96-core node with all cores busy, the `run` took 149 s with one
OpenMP thread and 358 s with four, and each `inspect` or `reconstruct` call
took 1–3 s.
