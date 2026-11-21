# Observable words

An **observable word** is the trace expression that names the quantity LSQuant
computes, for example `Tr[VX delta(E-H) VY dG-(E)/dE]`. You give it in the
`word` (or `words`) key of a [run file](run-file.md#word-and-words). Every word
must match one of seven **word shapes**. The shape fixes the moment tensor that
LSQuant computes and the reconstruction formulas that can turn it into a
curve.

The source of truth is `python/lsquant/words/`: the parser (`grammar.py`), the
shape catalogue (`catalogue.py`) and the compiler (`compiler.py`). To see the
catalogue and how one word is compiled:

```bash
python -m lsquant words list
python -m lsquant words describe "Tr[VX delta(E-H) VY dG-(E)/dE]"
```

See the [command-line reference](cli.md#words).

## Syntax

A word is `Tr[`, a product of factors, and `]`. Separate the factors with
spaces; the parser accepts and ignores a `*` between two factors.

| Factor | Written as | Meaning |
|---|---|---|
| Spectral function | `delta(E-H)` | δ(E − H), expanded in Chebyshev polynomials. |
| | `dG-(E)/dE` | Energy derivative of the advanced Green's function, expanded. |
| | `G+(E)`, `G-(E)`, `dG+(E)/dE` | Parsed, but no shape uses them: a word that contains them is refused. |
| Time evolution | `U(t)`, `U(t)^+` | e^(−iHt/ħ) and its adjoint, applied to vectors, not expanded. |
| Commutator | `[X,U(t)]`, `[X,U(t)]^+` | Position commutator of the mean-square displacement. Axis `X`, `Y`, `Z`, `VX`, `VY` or `VZ`. |
| Operator | `NAME` | A static operator: a letter, then letters, digits or `_`. See [Operator names](#operator-names). |
| | `1` | The identity. |
| | `P:NAME` | A projector; LSQuant looks up the operator `NAME`. |
| | `NAME^+` | Adjoint. Every shape refuses it: the operators of a word must be Hermitian. |
| | `NAME@real`, `NAME@bloch`, `NAME@eigen` | Basis tag of the operator. |
| | `B^+ NAME B`, `U^+ NAME U` | Explicit basis transform (Bloch to real space, eigenbasis to Bloch). |
| Fermi-sea term | `Tr[sea[ … ]]` or `sea[Tr[ … ]]` | Marks the Fermi-sea term of a Kubo–Bastin word. |

The **rank** of a word is its number of expanded spectral factors
(`delta(E-H)`, `dG-(E)/dE`). A word of rank 3 or more is refused with
`rank >= 3 moment tensors are out of scope`.

A word with `U(t)` or `[X,U(t)]` needs a [`time`](run-file.md#time) section in
the run file.

## Word shapes

From `CATALOGUE` and `match_catalogue` in `words/catalogue.py`:

| Shape | Pattern | Example | Rank | Time | Moment file prefix, extension |
|---|---|---|---|---|---|
| `dos_identity` | `Tr[delta(E-H)]` | `Tr[delta(E-H)]` | 1 | no | `SpectralOp1`, `.chebmom1D` |
| `dos_operator` | `Tr[A delta(E-H)]` | `Tr[SZ delta(E-H)]` | 1 | no | `SpectralOp<A>`, `.chebmom1D` |
| `kubo_greenwood` | `Tr[A delta(E-H) B delta(E-H)]` | `Tr[VX delta(E-H) VX delta(E-H)]` | 2 | no | `NonEqOp<A>-<B>`, `.chebmom2D` |
| `kubo_bastin` | `Tr[A delta(E-H) B dG-(E)/dE]`, optionally in `sea[…]` | `Tr[VX delta(E-H) VY dG-(E)/dE]` | 2 | no | `NonEqOp<A>-<B>`, `.chebmom2D` |
| `msd` | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | 1 | yes | `Correlation<V>`, `.chebmomTD` |
| `projected_evolution` | `Tr[P U(t)^+ A U(t) delta(E-H) P]` | `Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]` | 1 | yes | `TimeProj<A>-<P>`, `.chebmomTD` |
| `time_correlation` | `Tr[A delta(E-H) U(t)^+ B U(t)]` | `Tr[A delta(E-H) U(t)^+ B U(t)]` | 1 | yes | `TimeCorr<A>-<B>`, `.chebmomTD` |

The moment-file names use the operator names in place of `<A>`, `<B>`, `<P>`
and `<V>`. With `output.legacy_layout = false`, the unified `.chebmom` file is
written instead; see [Output files](output-files.md#moment-files).

Rules inside the shapes:

| Shape | Rule |
|---|---|
| `dos_operator` | `Tr[1 delta(E-H)]` gives the plain DOS. |
| `kubo_greenwood` | A `sea[…]` word with two `delta(E-H)` matches no shape. |
| `kubo_bastin` | The second spectral factor is `dG-(E)/dE`. |
| `msd` | Both commutators have the same axis. The axis selects the velocity operator: `X` or `VX` gives `VX`, `Y` or `VY` gives `VY`, `Z` or `VZ` gives `VZ`. |
| `projected_evolution` | The first and the last operator are the same. |

A word that matches no shape is refused with the nearest shape, for example:
`Unrecognised word shape 'Tr[VX G+(E) VY delta(E-H)]'; nearest recognised shape is 'Tr[VX delta(E-H) VY dG-(E)/dE]' (kubo_bastin)`.

## Formulas by shape

`reconstruction.formulas` in the run file, and `--formula` of
[`reconstruct`](cli.md#reconstruct), must be allowed by the shape of the word
(`SHAPE_ALLOWED_FORMULAS` in `words/catalogue.py`). Matching is
case-insensitive. Aliases are in parentheses.

| Shape | Allowed formulas | Default |
|---|---|---|
| `dos_identity` | `dos` (`dos_curve`) | `dos` |
| `dos_operator` | `dos` (`dos_curve`) | `dos` |
| `kubo_greenwood` | `greenwood` (`kubo_greenwood`), `greenwood_fft` (`kubo_greenwood_fft`) | `greenwood` |
| `kubo_bastin` | `bastin` (`kubo_bastin`), `sea` (`kubo_bastin_sea`), `surf` (`kubo_bastin_surf`), `bastin_fft` (`kubo_bastin_fft`), `shc` | `bastin`; `sea` for a `sea[…]` word |
| `msd` | `time_correlation` (`time_correlation_curve`) | `time_correlation` |
| `projected_evolution` | `row_m0` (`row0`), `time_correlation` (`time_correlation_curve`) | `row_m0` |
| `time_correlation` | `time_correlation` (`time_correlation_curve`) | `time_correlation` |

| Formula | Curve |
|---|---|
| `dos` | DOS, or the operator-weighted DOS for `dos_operator`. |
| `greenwood` | Kubo–Greenwood conductivity σ(E). |
| `greenwood_fft` | Kubo–Greenwood σ(E) on the Chebyshev-node grid. |
| `bastin` | Kubo–Bastin σ(E), Fermi-sea and Fermi-surface terms together. |
| `sea`, `surf` | The Fermi-sea and the Fermi-surface term of Kubo–Bastin. |
| `bastin_fft` | Kubo–Bastin σ(E) on the Chebyshev-node grid. |
| `shc` | Kubo–Bastin with the spin Hall normalisation; needs `num_cells` and `cell_area`. |
| `time_correlation` | One time curve per entry of `reconstruction.fermi_energies`. |
| `row_m0` | No curve file; the time series stays in the moment file. |

[Output files](output-files.md#curve-files) gives the columns and units of
each curve. The SIESTA executors accept only `dos`, `dos_curve` and `shc`; see
[SIESTA executors](run-file.md#siesta-executors).

## Operator names

### Names the compiler accepts

When the run file leaves out `model.operators`, the compiler accepts these
names (`DEFAULT_KNOWN_OPERATORS`, `words/compiler.py:28`):

| Group | Names |
|---|---|
| Identity | `1` |
| Velocity | `VX`, `VY`, `VZ` |
| Spin | `SX`, `SY`, `SZ` |
| Spin current | `JXSX`, `JXSY`, `JXSZ`, `JYSX`, `JYSY`, `JYSZ`, `JZSX`, `JZSY`, `JZSZ` |
| Velocity times spin | `VXSX`, `VXSY`, `VXSZ`, `VYSX`, `VYSY`, `VYSZ`, `VZSX`, `VZSY`, `VZSZ` |
| Projectors | `P`, `PX`, `PY`, `PZ` |
| Generic | `A`, `B`, `O` |

With `model.operators` in the run file, the validator accepts only `1` and
the listed names: the list replaces the default list
(`words/compiler.py:294–300`).
If `P` is not accepted but `PX` is, `P` stands for `PX`
(`words/compiler.py:397–398`).

### Names each model provides

An accepted name must also exist in the model:

| `model.kind` | Operators |
|---|---|
| `csr` | One per file `<label>.<NAME>.CSR` in `model.dir`. Without `model.operators`, every such file except `HAM` and `S` is loaded. |
| `wannier` | `VX`, `VY`, `VZ`, built from the hoppings and the positions. |
| `siesta` | `VX`, `VY`, `VZ`, `SX`, `SY`, `SZ`, the spin currents `J<a>S<b>`, and `WINDOW` in projection runs. The executors restrict the words further; see [Use a SIESTA calculation](../how-to/use-siesta.md). |
| `generator` | Built by the tutorial driver. |

For `csr` models, give your operator files names from the default list, for
example `<label>.A.CSR` for a layer projector, as
`examples/15_spatial_operators` does.

> **Known limitation: custom operator names fail in `run`.** A name outside
> `DEFAULT_KNOWN_OPERATORS`, for example `PTOP`, passes validation and
> `inspect` when it is listed in `model.operators`, but
> `python -m lsquant run` stops with `Unknown operator 'PTOP'` for `csr` and
> `wannier` models. Before it runs each word, `execute_runfile` compiles it
> again without the model's operators (`python/lsquant/runfile.py:1546–1554`),
> so the compiler falls back to the default list (`words/compiler.py:28`,
> `393–394`) and raises at `words/compiler.py:405`. `words describe` refuses
> such a name too. Use names from the default list, as above.

> **Known limitation: `WINDOW` needs `model.operators`.** `WINDOW` is not in
> the default list, so the validator refuses `Tr[WINDOW delta(E-H)]` unless the
> run file has `"operators": ["WINDOW"]` in `model`. The word runs only in
> SIESTA projection runs, which skip the compile step of the previous box; see
> `examples/16_projection_downfold`. Because of the replacement rule in
> [Names the compiler accepts](#names-the-compiler-accepts), list every
> operator of every word.

## Examples

| Quantity | Word | Formula |
|---|---|---|
| Density of states | `Tr[delta(E-H)]` | `dos` |
| DOS projected on a region or orbital set | `Tr[A delta(E-H)]` | `dos` |
| Spin-resolved DOS | `Tr[SZ delta(E-H)]` | `dos` |
| Longitudinal conductivity σ_xx | `Tr[VX delta(E-H) VX delta(E-H)]` | `greenwood` |
| Hall conductivity σ_xy | `Tr[VX delta(E-H) VY dG-(E)/dE]` | `bastin` |
| Fermi-sea term of σ_xy | `Tr[sea[VX delta(E-H) VY dG-(E)/dE]]` | `sea` |
| Spin Hall conductivity σ^z_xy | `Tr[JXSZ delta(E-H) VY dG-(E)/dE]` | `shc` |
| Mean-square displacement along x | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | `time_correlation` |
| Spin precession of a projected state | `Tr[P:PX U(t)^+ SX U(t) delta(E-H) P:PX]` | `row_m0`, `time_correlation` |
| Time correlation of two operators | `Tr[A delta(E-H) U(t)^+ B U(t)]` | `time_correlation` |

[How we run LSQuant](../how-to/run-lsquant.md#the-quantity)
points to a use case for each of these.
