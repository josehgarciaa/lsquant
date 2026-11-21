# Command-line interface

`python -m lsquant` is the command-line front end of LSQuant.
`python/pyproject.toml` also declares an `lsquant` console script, but
`pip install -e python/` currently fails (see
[Install LSQuant](../how-to/install.md)), so run the commands with
`python -m lsquant`. The compiled module must be importable:

```bash
export PYTHONPATH="$PWD/build/golden/lsquant_core:$PWD/python"
```

Every command prints `lsquant error: <message>` to standard error and exits
with status 1 if it fails. On success it exits with status 0.

## Synopsis

```text
python -m lsquant run <runfile> [--out DIR] [--json]
python -m lsquant reconstruct <target> [FORMULA] [BROADENING_MEV]
                              [--formula NAME] [--broadening-meV MEV]
                              [--fermi-energies E [E ...]] [--out DIR] [--json]
python -m lsquant inspect <path>
python -m lsquant words [list]
python -m lsquant words describe "<word>"
python -m lsquant words "<word>"
```

## `run`

Validates a JSON run file, computes the Chebyshev moments, and reconstructs the
curves the run file asks for.

| Argument | Meaning |
|---|---|
| `runfile` | Path to the JSON run file. See [Run-file reference](run-file.md). |
| `--out DIR`, `-o DIR` | Output directory. Overrides `output.dir` in the run file. |
| `--json` | Print the execution summary as JSON (`output_dir`, `resolved_runfile`, `moment_files`, `curve_files`, `written_files`) instead of one written path per line. |

The run file is validated before any input is read. An unknown key, at the top
level or inside any section, stops the run with an error.

## `reconstruct`

Rebuilds curves from the moments that a previous `run` wrote, without
calculating the moments again. Use it to try another broadening or another
formula on the same moments.

| Argument | Meaning |
|---|---|
| `target` | A moment file (`*.chebmom1D`, `*.chebmom2D`, `*.chebmomTD`, `*.chebmom`) or a run output directory. For a directory, every moment file in it is reconstructed. |
| `FORMULA` (positional, optional) | Same as `--formula`. |
| `BROADENING_MEV` (positional, optional) | Same as `--broadening-meV`. |
| `--formula NAME`, `-f NAME` | Reconstruction formula; see the table below. |
| `--broadening-meV MEV`, `-b MEV` | Target energy resolution of the Jackson kernel, in meV. Must be positive. |
| `--fermi-energies E ...`, `--fermi-energy`, `-e` | Fermi energies in eV, used by `time_correlation`. |
| `--out DIR`, `-o DIR` | Directory for the curve files. Default: the directory of the moments. |
| `--json` | Print the summary as JSON instead of one curve path per line. |

Defaults when an option is omitted:

- **Formula.** `time_correlation` if the moments have time steps, `bastin`
  if they are two-dimensional, and `dos` otherwise.
- **Broadening.** For a directory, the value of
  `reconstruction.broadening_meV` in its `run.resolved.json`. Otherwise
  10 meV.
- **Fermi energies.** For a directory, `reconstruction.fermi_energies` in
  `run.resolved.json`. Otherwise `[0.0]`.

Formulas (aliases in parentheses):

| Formula | Moments | Result |
|---|---|---|
| `dos` (`dos_curve`) | 1D | Density of states. |
| `greenwood` (`kubo_greenwood`) | 2D | Kubo–Greenwood conductivity, for σ_xx. |
| `bastin` (`kubo_bastin`) | 2D | Kubo–Bastin conductivity, for σ_xy. |
| `sea` (`kubo_bastin_sea`) | 2D | Fermi-sea part of the Kubo–Bastin formula. |
| `surf` (`kubo_bastin_surf`) | 2D | Fermi-surface part of the Kubo–Bastin formula. |
| `greenwood_fft` (`kubo_greenwood_fft`) | 2D | Kubo–Greenwood on the Chebyshev-node grid, computed with an FFT. |
| `bastin_fft` (`kubo_bastin_fft`) | 2D | Kubo–Bastin on the Chebyshev-node grid, computed with an FFT. |
| `shc` | 2D | Spin Hall conductivity σ^z_xy in units of e²/h. It needs `num_cells` and `cell_area` from the `run.resolved.json` next to the moments. |
| `time_correlation` (`time_correlation_curve`) | time-dependent | Correlation against time, one curve per Fermi energy. |
| `row_m0` (`row0`) | 2D | First row of the moment matrix (diagnostic). |

## `inspect`

Prints JSON metadata for a path. Its output is always JSON, so the command has
no `--json` option.

| Path | Output |
|---|---|
| Run output directory | `kind: "directory"`, the resolved run file, the basis defaults, and the metadata of every moment file in it. |
| Run file (`*.json`) | `kind: "runfile"`, the run file after validation with every default filled in, and the basis defaults. The command stops after validation, so use it to check a run file before you submit it. |
| Moment file | `kind: "moments"`: rank, dimensions, system size, number of states, word, spectral window, band centre and width, kernel and symmetrisation flags, and the number of time steps and the time step. |

## `words`

Lists the observable word shapes, or parses and compiles one word.

| Form | Output |
|---|---|
| `words` or `words list` | The catalogue of supported word shapes, as JSON. |
| `words describe "<word>"` | The compiled word, as JSON: `word` (the parse tree), `catalogue` (the matched shape, with `shape_id`, `rank` and `time_dependent`), `plan` (the engine instructions), `operators`, `reconstruction` (default formula, kernel, broadening) and `dt`. |
| `words "<word>"` | Same as `describe`. |

Quote the word in the shell, because it contains brackets and spaces:

```bash
python -m lsquant words describe "Tr[VX delta(E-H) VY dG-(E)/dE]"
```

## See also

- [Run-file reference](run-file.md)
- [Observable words](observable-words.md)
- [MCP tools](mcp-tools.md), which expose the same operations to an AI assistant.
