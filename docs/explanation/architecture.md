# Architecture

A KPM code spends almost all its time in one loop: sparse matrix–vector
products and inner products on long vectors. The choice of observable, the
basis, the overlap matrix, the Kubo formula and its units are what make a
calculation physical. They cost little by comparison, but they change far
more often. LSQuant separates the two. A small C++ engine runs the loop. A
Python front end holds the physics and drives the engine. Besides the loop,
the C++ side keeps only the contraction helpers of the reconstruction and the
legacy table evaluators in `lsquant/compat/`.

This page describes the parts and why they are arranged this way. How
observables become traces is explained in [The kernel polynomial method](kpm.md).

## The C++ moment engine

The engine is in `lsquant/` and has three layers with one-way
dependencies:

- `lsquant/core/algebra/` provides linear operators ($y \leftarrow a\,O x + b\,y$),
  the sparse CSR matrix, composite operators (sums, products, the Löwdin and
  metric-inverse products), the overlap appliers, the Bloch transform, the
  spectral bounds and the generators of trace vectors (unit, random-phase and block
  vectors);
- `lsquant/core/chebyshev/` provides the spectral window, the rescaled
  Hamiltonian, the instruction set, the `MomentEngine` that executes it, the
  `ChebyshevMoments` tensor and the contraction kernels used by the
  reconstruction;
- `lsquant/io/` reads CSR text files and writes moment and curve files.

The engine executes a `Plan`, the list of operations it carries out on
vectors. There are seven instructions. `ApplyOperator`
multiplies the live vectors by a registered operator. `Expand` opens a
Chebyshev index, replacing a vector by $T_0(\tilde H)v, \dots, T_{M-1}(\tilde H)v$.
`ApplySeries` and `ApplyDerivedSeries` apply a Chebyshev series without
opening an index; time evolution uses them. `MirrorToBra` starts the bra side
from a fresh copy of the starting vector, `TimeAxis` opens a time index, and
`Trace` contracts bra and ket into the moment tensor.
<!-- lsquant/core/chebyshev/plan.hpp -->

The rank of the result is the number of `Expand` instructions. None of the
instructions names a velocity, a spin or a delta function. The test script
described below checks this by searching
`lsquant/core` and `lsquant/io` for physical words such as `delta`, `Kubo`,
`velocity`, `spin` and `basis`.
<!-- scripts/check_goldens.sh, stages 1 to 4 -->

For example, the Kubo–Bastin word `Tr[VX delta(E-H) VY dG-(E)/dE]` becomes
the following plan:
<!-- python/lsquant/moments.py, kubo_plan -->

```python
[Expand(M), ApplyOperator(VY),      # ket: VY T_n(H~) |r>
 MirrorToBra(),
 ApplyOperator(VX), Expand(M),      # bra: T_m(H~) VX |r>
 Trace([1])]                        # mu[m, n] += <bra_m | ket_n>
```

The same plan serves the Kubo–Greenwood word with `VX` on both sides. The
difference between the two formulas appears only later, in the
reconstruction.

`lsquant/compat/` holds evaluators of the Green's-function and Bessel tables
in the form used by the earlier codes. They are kept only so that the
reference curves match byte for byte, and the core never includes them.
<!-- lsquant/compat/legacy_kernels.hpp:1-8 -->

## The Python front end

The pybind11 module `lsquant._core` exposes the engine to Python.
<!-- lsquant/api/bindings.cpp --> Everything above it is Python in `python/lsquant/`:

- `words/` parses an observable word, matches it to one of the shapes in the
  catalogue, and compiles it into a `Plan` and a list of operators;
- `bases.py`, `disorder_bases.py` and `models/` build the operators and trace
  states for the `real`, `bloch` and `eigen` bases, for Wannier and SIESTA
  models and for non-orthogonal metrics. A basis is only an adapter that hands
  the engine linear operators;
- `runfile.py` validates a run file, fills in defaults, writes
  `run.resolved.json`, and passes the run to the code for its model kind;
- `reconstruct.py`, `kernels.py` and `spectral_functions.py` turn moments into
  curves. The coefficient tables are computed in Python and contracted with
  helper routines of `_core`;
- `metric_forms.py` states the non-orthogonal formulas symbolically with SymPy,
  as the source that the numerical code implements.

The engine stops at the moment tensor. Python drives the reconstruction and
calls only small contraction helpers of `_core`. The run writes the tensor to
disk next to the curves. This split is what makes
`python -m lsquant reconstruct` possible: a new broadening or a new formula
needs only the moment file.

## The command line and the MCP server

Three entry points drive the same Python code. The command line,
`python -m lsquant`, has four commands ([CLI reference](../reference/cli.md)):

- `words` lists and describes the accepted word shapes;
- `run` executes a run file;
- `reconstruct` turns a moment file into a curve;
- `inspect` prints the metadata of a moment file, run file or output directory.

Python scripts can call `execute_runfile` or build plans themselves.

The MCP server, `python -m lsquant.mcp_server`, exposes the same operations as
five tools for an AI assistant: `words_describe`, `run`, `status`,
`reconstruct` and `inspect`. The `run` tool validates the run file and starts
`python -m lsquant run` in a background process, and `status` reports on it.
Every path must lie inside one root directory, set by `LSQUANT_MCP_ROOT` or
the command line ([MCP reference](../reference/mcp-tools.md)).

## Build profiles

LSQuant builds in two arithmetic profiles. The `golden` profile uses portable
IEEE arithmetic with no machine-specific flags, so its results are
reproducible across machines. The `native` profile turns on `-march=native`,
SIMD and cache-blocked kernels, and its results differ from the golden ones in
the last digits. The [build-presets reference](../reference/build-presets.md)
lists the presets built on each profile.

## Golden files instead of unit tests

LSQuant unifies three earlier codes. The question that matters for each change
is whether the new code still gives the numbers the old code gave, for the
whole pipeline from operator files to curves. A KPM curve can be wrong in
many small ways that escape a unit test of one function: a factor
$2 - \delta_{m0}$ applied twice, a window rescaled in the wrong order, a
missing normalisation by the number of random vectors.

The project therefore relies on regression testing of the whole pipeline
against golden masters. A golden master, or golden file, is a reference
output produced once by a trusted code and kept unchanged. LSQuant's golden
files were generated from the earlier codes and committed under
`goldens/stage0/`, with their provenance in `goldens/stage0/PROVENANCE.md`.
They are never regenerated from the new code.

The new code was written in stages. Each stage added one part: the moment
core, the Python binding, the rank-2 reconstruction, time evolution, the
Bloch and eigen bases, the non-orthogonal metric, and finally the words, run
files, command line and MCP server. The gate is the script
`scripts/check_goldens.sh <stage>`. For a given stage it runs the checks of
that stage and of the stages below it on the committed inputs, and it
compares every output file with `cmp`, byte for byte. A difference that is a
pure floating-point reordering, below $10^{-12}$ relative, is reported
separately and has to be accepted explicitly.

Each golden also has a physical check where one exists: the chain DOS moments
$\mu_m = \delta_{m0}$, the quantised Hall plateau $C = 1$ of the Haldane model,
the ballistic $\Delta X^2 \propto t^2$ of a clean chain, exact
diagonalisation for the Bloch and eigen bases, and agreement of the
non-orthogonal routes with the orthogonal goldens. Byte identity shows
that the output is unchanged, and the physical check shows that this
unchanged output is right.

The choice has a visible cost in the code. Reproducing the output of the
earlier codes bit for bit fixes the order of floating-point operations.
Several routines therefore use scalar `math` functions where NumPy would be
faster, and their docstrings warn against vectorising them. A `native` build
differs from the golden files in the last digits (see
[Build profiles](#build-profiles)), so it is checked with the tolerance mode
of the same script (`GOLDEN_TOL=1`).
<!-- python/lsquant/kernels.py:49-51 -->

A stage may carry one developer self-check in `tests/dev/`, used while
implementing that stage. These checks stay outside the gate.
