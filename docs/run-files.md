\page run_files Run files, observable words, the CLI and the MCP server

Stage 7 puts one interface on top of the stages before it. A calculation is a JSON **run
file**; the quantity is an **observable word** such as `Tr[VX delta(E-H) VY dG-(E)/dE]`; the
command line tool `python3 -m lsquant` runs, reconstructs and inspects; an MCP server
offers the same operations to AI assistants. Every golden of stages 0–6 is reproduced
byte for byte from a run file under `runs/goldens/`.

## Observable words

A word is a trace of operators and spectral functions of H:

| Shape | Example | Reconstruction formulas |
|---|---|---|
| `Tr[delta]` | `Tr[delta(E-H)]` | `dos` |
| `Tr[A delta]` | `Tr[SZ delta(E-H)]` | `dos` |
| `Tr[A delta B delta]` | `Tr[VX delta(E-H) VX delta(E-H)]` | `greenwood`, `greenwood_fft` |
| `Tr[A delta B dG/dE]` | `Tr[VX delta(E-H) VY dG-(E)/dE]` | `bastin`, `sea`, `surf`, `bastin_fft` |
| `Tr[A delta U(t)^+ B U(t)]` | `Tr[VX delta(E-H) U(t)^+ VX U(t)]` | `time_correlation` |
| `Tr[P U(t)^+ A U(t) delta P]` | `Tr[PX U(t)^+ SZ U(t) delta(E-H) PX]` | `time_correlation`, `row0` |
| `Tr[[X,U(t)]^+ delta [X,U(t)]]` | `Tr[[X,U(t)]^+ delta(E-H) [X,U(t)]]` | `time_correlation` (MSD) |

`python3 -m lsquant words list` prints the catalogue, and `python3 -m lsquant words
describe "<word>"` prints what a word compiles to: the plan instructions, the operators
and the allowed formulas. A word outside the catalogue is refused with the nearest
catalogue entry; the compiler never approximates a word by another one. Operators are
named after the CSR files of the model (`VX`, `SZ`, `PX`, …) and must be Hermitian.

## Run files

A run file has up to nine sections. Unknown keys, at the top level or inside any section,
are refused before anything is loaded.

| Section | Keys | Content |
|---|---|---|
| `model` | `kind` = `csr` (`dir`, `label`, `operators`), `wannier` (`dir`, `seed`, `csr_dir`) or `generator` (`name`, `params`) | where H and the operators come from; `generator` runs one of the model scripts of `examples/` |
| `basis` | `name` (`real`, `bloch`, `eigen`), `supercell`, `disorder` (`W`, `seed`) | the basis of \ref bases |
| `metric` | `policy` (`orthogonal`, `nonorthogonal`), `route`, `applier`, `order`, `coeff_tol`, `tol`, `max_iter`, `max_dim` | the overlap of \ref nonorthogonal |
| `word`, `words` | one word, or a list | what to compute |
| `expansion` | `moments`, `window` (`[lo, hi]`, `gershgorin`, `eigen`), `trace` (`exact`, `random`), `states`, `seed`, `batch`, `state_label` | the Chebyshev expansion |
| `time` | `steps`, `tmax`, `dt`, `bessel` | the time axis of \ref time_evolution |
| `reconstruction` | `kernel`, `broadening_meV`, `fermi_energies`, `formulas`, `num_points`, `fft_grid`, `legacy_text_round_trip`, `per_site_dos`, `targets` | the curves |
| `output` | `dir`, `legacy_names`, `legacy_layout`, `write_chebmom`, … | where and how files are written |

`model`, `word` and `expansion` are required. `window = "eigen"` needs a Bloch or eigen
basis and is refused on a `csr` model. Each run writes `run.resolved.json`, the run file with
every default filled in.

## The command line

```bash
export PYTHONPATH=$PWD/build/lsquant_core:$PWD/python
python3 -m lsquant run <runfile> [--out <dir>] [--json]
python3 -m lsquant reconstruct <moment file | output dir> [--formula f] [--broadening-meV b] [--fermi-energies e ...] [--out <dir>]
python3 -m lsquant inspect <moment file | run file | output dir>        # always prints JSON
python3 -m lsquant words list | describe "<word>"
```

`reconstruct` reads stored moments, so a new broadening or formula costs no new moments. A
refused run file exits with a non-zero status, prints the reason on stderr and creates no
output directory.

## The MCP server

```bash
python3 -m lsquant.mcp_server --root <project dir>
```

starts a stdio server (it needs the `mcp` package; `--root` defaults to `LSQUANT_MCP_ROOT`, then
the repository root) with five tools: `words_describe`,
`run` (validates the run file and starts the run in the background), `status`,
`reconstruct` and `inspect`. Every path in a request must lie inside `--root`; a model or
output directory outside it is refused. The client configuration is described in
`python/lsquant/README.md`.

## Example: the Haldane model from a run file

`examples/run_files/` holds a run file for the Hall and longitudinal conductivities of the
golden Haldane model and `run.sh`, which goes through the whole workflow: describe the word,
inspect the run file, run it, reconstruct again at 40 meV from the stored moments, and
extract the Chern number from both curves (`chern.py`):

```bash
examples/run_files/run.sh
```

The output and the timings are in `examples/run_files/README.md`.

## Tests

| Test | What it checks | Run |
|---|---|---|
| `tests/unit/test_python_stage7.py` | words: parse, catalogue, compile (one `Expand` per spectral slot), refusals, `describe`; run files g1 and g5 reproduce G1 and the four G5 curves byte for byte; validation refuses unknown keys, a missing word and `window = "eigen"` on a csr model; CLI `inspect` and a refused run; MCP `words_describe`, `run`/`status` inside a temporary root, and a model directory outside the root refused | `PYTHONPATH=build/lsquant_core:python python3 -m pytest tests/unit` |
| `tests/dev/stage_7_check.py` | the run files of G1–G7 and the tutorials, the CLI and the MCP tools, as a standalone script | `PYTHONPATH=... python3 tests/dev/stage_7_check.py` |
| stage-7 gate | every run file of `runs/goldens/` at 16 and 1 threads (G1–G7, G4c, tutorials 01–06 byte for byte; G9 and G10 against the stage-5 and stage-6 programs and the tolerance tables), then the MCP server, the CLI refusals and the greps of stages 1–7. `S7_SET=quick` runs the closure subset | `W2S_BIN=... W2S_TOOLS=... scripts/check_goldens.sh 7 build` |

The tutorial-06 goldens were rewritten from the new code in this stage because the legacy
pipeline never built the identity operator of its denominator (legacy defect D5); the old
files are kept in `goldens/stage0/legacy_defects/tutorial06_identity_operator/`, and
`goldens/stage0/PROVENANCE.md` (G8) explains the change.

The complete run-file reference, the observable-word reference and the user guide come
with the documentation pull requests that follow the stage series.
