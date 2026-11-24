# Core module `lsquant._core`

`lsquant._core` is the compiled C++ engine, bound to Python with pybind11. The
Python package builds every calculation from it: operators, trace states, the
bytecode plan of a word, the moment engine, and the moment and curve files.
Most users need the [Python API](python-api.md) or a [run file](run-file.md)
instead. Use this page when you write your own operator, state set or plan.

```python
from lsquant import _core
```

Import it through `lsquant`, as above, rather than as a top-level `_core`.
Importing `lsquant` first sets `OMP_NUM_THREADS`, and, only when the process is
the `lsquant` command with its own cores (or `LSQUANT_BIND_THREADS=1`),
`OMP_PROC_BIND` and `OMP_PLACES`; each only if it is
unset (`python/lsquant/runtime.py`, `configure_threads`). In every build it then calls
`ConfigureEngineDiagnostics` with the settings of `LSQ_ENGINE_MEM`,
`LSQ_ENGINE_TIMING` and `LSQ_CHECKPOINT_ABORT_AFTER` (`runtime.py:44–49`). In
native builds (`build_profile == "native"`), it also calls
`ConfigurePerformance` with the values of `LSQUANT_SPMV_MIN_NNZ` and
`LSQUANT_KERNEL_BATCH` (defaults `16384` and `0`; `runtime.py:50–54`). The
engine receives these settings only through the two calls.

The source of truth is `lsquant/api/bindings.cpp`. The type stub
`python/lsquant/_core.pyi` is incomplete; differences are listed in
[Differences from the type stub](#differences-from-the-type-stub).

Conventions used on this page:

- `complex` arrays are NumPy `complex128`; `int64` arrays are NumPy `int64`.
- "1D contiguous" means a one-dimensional C-contiguous array of exactly that
  dtype. Functions marked *no conversion* raise `TypeError` for any other
  dtype or layout instead of converting the array.
- Every operator applies `y = a·O·x + b·y` in place internally. Errors raised
  in C++ as `std::invalid_argument` arrive in Python as `ValueError`, and
  `std::runtime_error` as `RuntimeError`.
- *Keeps alive* means the object holds a reference to its arguments, so they
  live at least as long as the object.

## Module attributes

| Name | Type | Meaning |
|---|---|---|
| `build_profile` | `str` | `"golden"` or `"native"`: the `LSQUANT_BUILD_PROFILE` of the build ([build presets](build-presets.md)). |
| `ginkgo_available` | `bool` | `True` if built with `LSQUANT_WITH_GINKGO=ON`. |
| `cuda_available` | `bool` | `True` if built with `LSQUANT_WITH_CUDA=ON`. |
| `mpi_available` | `bool` | `True` if built with `LSQUANT_WITH_MPI=ON`. No MPI class or function is bound; the flag only reports the build. |

## Performance configuration

Call these once, before any engine runs.

| Function | Meaning |
|---|---|
| `ConfigurePerformance(spmv_min_nnz=16384, kernel_batch=0)` | Sets two process-wide kernel settings. `spmv_min_nnz`: a sparse product runs in an OpenMP parallel region only when the operator has at least this many non-zeros (atom pairs for `AtomBlockOperator`). `kernel_batch`: maximum number of vectors per cache micro-batch in block products; `0` lets the kernel choose from a 512 KiB L2 budget. Both must be ≥ 0, otherwise `ValueError`. Returns `None`. |
| `ConfigureEngineDiagnostics(memory_report=False, phase_timing=False, checkpoint_abort_after=-1)` | Sets three process-wide diagnostic settings. `memory_report` and `phase_timing` print the allocation and the phase timings of the native rank-2 path to stderr. `checkpoint_abort_after = k` stops a checkpointed run right after the k-th checkpoint write; `-1` never stops. A value below `-1` raises `ValueError`. Returns `None`. |
| `ConfigureBenchmarkThreads(count)` | Calls `omp_set_num_threads(count)` and returns `omp_get_max_threads()`. `count` ≤ 0 raises `ValueError`. |

`import lsquant` makes this call in native builds, as described at the top of
this page. The engine receives these values only through this call. Only the
native code paths of the kernels read the settings (for example
`lsquant/core/algebra/sparse_matrix.cpp:136`).

## Sparse matrices

### `SparseMatrix`

A square complex matrix in CSR format.

| Member | Signature | Meaning |
|---|---|---|
| constructor | `SparseMatrix()` | Empty 0 × 0 matrix. |
| static | `SparseMatrix.from_csr(dim, row_pointers, columns, values)` | Builds a matrix from CSR arrays. *No conversion*: `row_pointers` and `columns` 1D contiguous `int64`, `values` 1D contiguous `complex128`. |
| property | `dimension` → `int` | N. |
| property | `non_zeros` → `int` | Stored entries, explicit zeros included. |
| method | `rescaled(a, b)` → `SparseMatrix` | New matrix `a·M + b·I`. |
| method | `hermitian_fold()` → `SparseMatrix` | New matrix `(M + M†)/2`. |
| method | `is_hermitian(tol=1e-12)` → `bool` | `True` if the relative Frobenius norm of `M − M†` is within `tol`. |
| method | `csr_arrays()` → `(row_pointers, columns, values)` | Copies of the CSR arrays (`int64`, `int64`, `complex128`). |

`from_csr` parameters:

| Parameter | Type | Meaning |
|---|---|---|
| `dim` | `int` | Matrix dimension N ≥ 0. |
| `row_pointers` | `int64[N+1]` | Row offsets, starting at 0 and ending at nnz. |
| `columns` | `int64[nnz]` | Column indices in `[0, N)`. |
| `values` | `complex128[nnz]` | Non-zero values. |

### `ReadCsrText(path, legacy_fold=False)` → `SparseMatrix`

Reads an LSQuant `.CSR` text file: `dim`, `nnz`, then `nnz` value pairs
`(re, im)`, `nnz` column indices, and `dim + 1` row offsets, in that order. It
computes ‖(M − M†)/2‖_F / ‖M‖_F, refuses the file if this exceeds `1e-10`
unless `legacy_fold=True`, and returns the folded matrix `(M + M†)/2`.

## Linear operators

### `LinearOperator` (base class)

Not constructible. Every operator below derives from it.

| Member | Signature | Meaning |
|---|---|---|
| property | `dimension` → `int` | N. |
| property | `is_identity` → `bool` | `True` only for a structural identity. |
| property | `is_hermitian` → `bool` | As reported by the operator. |
| method | `multiply(x)` → `complex128[N]` | Returns `O·x`. `x` is converted to 1D `complex128`. |
| method | `multiply(a, x, b, y)` → `None` | `y ← a·O·x + b·y`. `y` must be a 1D contiguous `complex128` array of length N; it is updated in place. |
| method | `multiply_block(a, x, b, y)` → `None` | `Y ← a·O·X + b·Y` on N × S blocks. `x`, `y` must be 2D Fortran-ordered `complex128`. |

### Operator classes

| Class | Constructor | Meaning | Extra members |
|---|---|---|---|
| `SparseMatrixOperator` | `(matrix)` | Operator view of a `SparseMatrix`. Keeps `matrix` alive. | — |
| `IdentityOperator` | `(dim)` | I of dimension `dim`. | — |
| `DiagonalOperator` | `(values)` | diag(`values`); `values` is converted to 1D `complex128`. | `adjoint()` → diag(conj(values)) |
| `ScaledOperator` | `(alpha, op)` | `alpha·op`. Keeps `op` alive. | — |
| `OperatorSum` | `(first, second)` | `first + second`. Same dimension required. Keeps both alive. | — |
| `OperatorProduct` | `(first, second)` | `first·second` (`second` applied first). Keeps both alive. | — |
| `AdjointOperator` | `(op)` or `(matrix)` | `op†`, or `matrix†` stored as a copy. From an operator: raises `ValueError` if `op` is not Hermitian and not of a type whose adjoint is known. | — |
| `SandwichOperator` | `(B, V)` | `B·V·B†`. Keeps both alive. | — |
| `LowdinOperator` | `(inverse_sqrt, operator)` | `S^(−1/2)·O·S^(−1/2)`. Keeps both alive. | — |
| `MetricInverseProduct` | `(inverse, operator)` | `S^(−1)·O`. Keeps both alive. | — |
| `HermitianTripleProduct` | `(first, middle, last, coefficient=0.5)` | c·A·M·B + conj(c)·B·M·A with A = `first`, M = `middle`, B = `last`, c = `coefficient`, for Hermitian factors. c = ½ gives the symmetric product, c = i a Hermitian commutator. Keeps all three alive. | — |
| `BlockDiagonalOperator` | `(blocks)` or `(num_blocks, block_size, blocks)` | N_k dense W × W blocks. First form: array of shape `(N_k, W, W)`. Second form: flat row-major array of length N_k·W·W. | `num_blocks`, `block_size`, `adjoint()` |
| `BlochTransform` | `(n1, n2, n3, num_orbitals, atom_phases=None, direction=BlochDirection.RealToBloch)` | Unitary lattice Fourier transform B (or B† with `BlochToReal`) by FFTW, on an `n1 × n2 × n3` supercell with `num_orbitals` orbitals per cell. `atom_phases`: unscaled table `exp(−i k·τ_α)` of length N_k·W; the constructor multiplies it by 1/√N_k. `None` means τ_α = 0. | `n1`, `n2`, `n3`, `num_cells`, `num_orbitals`, `direction`, `adjoint()` |
| `AtomBlockOperator` | see below | Sparse operator stored as atom-pair blocks with a dictionary of distinct dense blocks. | `num_atoms`, `num_pairs`, `dictionary_size` |

`BlochDirection` is an enumeration with values `RealToBloch` and
`BlochToReal`.

### `AtomBlockOperator`

```text
AtomBlockOperator(atom_offset, pair_offset, nbr, block_id, dict_offset,
                  dict_rows, dict_cols, dict_values, diagonal, hermitian=True)
```

| Parameter | Type | Meaning |
|---|---|---|
| `atom_offset` | `int64[A+1]` | First orbital of each atom; `atom_offset[A]` = N. |
| `pair_offset` | `int64[A+1]` | First pair of each atom; `pair_offset[A]` = P. |
| `nbr` | integer `[P]` | Neighbour atom of each pair. 32- or 64-bit; stored as `int32`. |
| `block_id` | integer `[P]` | Dictionary block of each pair. |
| `dict_offset` | `int64[D+1]` | First value of each dictionary block in `dict_values`. |
| `dict_rows`, `dict_cols` | integer `[D]` | Shape of each block. |
| `dict_values` | `complex128` | Row-major values of all blocks. |
| `diagonal` | `complex128[N]` or empty | On-site term (for example Anderson disorder), applied after the pair sum. |
| `hermitian` | `bool` | Value reported by `is_hermitian`; not checked. |

Results agree with the CSR product of the same matrix to rounding (≈ 10⁻¹⁶
per row), not bit for bit, because the summation order differs.

### `AtomBlocksFromCsr(rowptr, columns, values, atom_map, permutation)` → `dict`

Builds the arrays of `AtomBlockOperator` from a CSR matrix. Returns a dict with
keys `atom_offset`, `pair_offset`, `nbr`, `block_id`, `dict_offset`,
`dict_rows`, `dict_cols`, `dict_values`, `diagonal`, ready to pass as keyword
arguments. *No conversion*: all inputs 1D contiguous, `values` `complex128`,
the others `int64`.

| Parameter | Meaning |
|---|---|
| `rowptr`, `columns`, `values` | CSR arrays of the matrix. |
| `atom_map` | Atom of each orbital (length N). |
| `permutation` | New-to-old orbital order; must keep each atom contiguous. |

The function rejects duplicate CSR entries. It applies no Hermitian folding.

### Optional operators

| Class | Built when | Constructor | Members |
|---|---|---|---|
| `GinkgoOperator` | `LSQUANT_WITH_GINKGO=ON` | `(matrix, executor="omp", format="csr")`. `executor`: `"omp"` or `"cuda"`. `format`: `"csr"`, `"sellp"`, `"ell"`, `"hybrid"` or `"coo"`. | `executor`, `format`, `multiply_block_row_major(a, x, b, y)` (row-major blocks, *no conversion*) |
| `CudaAtomBlockOperator` | `LSQUANT_WITH_CUDA=ON` | Same arguments as `AtomBlockOperator`. | `num_atoms`, `num_pairs`, `dictionary_size`, `multiply_block_device(a, x, b, y)` on `CudaVectorBlock` |
| `CudaVectorBlock` | `LSQUANT_WITH_CUDA=ON` | `(rows, cols)`: zero-initialised N × S block on the GPU. | `rows`, `cols`, `from_numpy(array)` (Fortran-ordered `(N, S)` `complex128`), `to_numpy()` |

With `LSQUANT_WITH_CUDA=ON` the module also has `synchronize_cuda()`, which
waits for the GPU and releases the GIL meanwhile.

## Overlap metric

These classes treat a non-orthogonal basis with overlap matrix S. See
the [`metric` keys of the run file](run-file.md#metric) for how the run file selects them.

### Metric policies

| Class | Constructor | Meaning |
|---|---|---|
| `MetricPolicy` | — (base class) | Read-only properties `is_orthogonal`, `overlap` (S), `inverse_overlap` (S⁻¹), `inverse_sqrt_overlap` (S^(−1/2)); each operator property returns a `LinearOperator`. |
| `OrthogonalMetric` | `(dim=0)` | S = S⁻¹ = S^(−1/2) = I. `is_orthogonal` is `True`. |
| `NonOrthogonalMetric` | `(s, inverse, inverse_sqrt)` | S from a `SparseMatrix`, with the two given appliers. `is_orthogonal` is `False`, even if S = I. Keeps all three alive. |

### Inverse appliers

| Class | Constructor | Method | Members |
|---|---|---|---|
| `ExactInverseOverlap` | `(s, max_dimension=4096)` | Dense Cholesky inverse of S. Raises `ValueError` if N > `max_dimension` or S is not Hermitian. | — |
| `ExactInverseSqrtOverlap` | `(s, max_dimension=4096)` | Dense S^(−1/2) from the eigendecomposition of S. Same limits. | — |
| `IschoInverseOverlap` | `(s, tolerance=1e-13, max_iterations=200)` | Solves S z = x by conjugate gradients with an incomplete Cholesky preconditioner, at every application. Keeps `s` alive. Raises `RuntimeError` if it does not converge within `max_iterations`. | `last_iterations`, `last_error` (relative residual of the last application) |
| `ChebyshevInverseOverlap` | `(s, bounds_of_s, order, coeff_tol=1e-12)` | Chebyshev polynomial of degree ≤ `order` for S⁻¹ on `bounds_of_s` (a `SpectralBounds` with 0 < lower < upper that encloses spec(S)). Terms with \|c_m\| < `coeff_tol`·\|c_0\| are dropped. Keeps `s` alive. | `effective_order` |
| `ChebyshevInverseSqrtOverlap` | `(s, bounds_of_s, order, coeff_tol=1e-12)` | The same for S^(−1/2). | `effective_order` |

There is no S^(−1/2) applier based on conjugate gradients, so the Löwdin route
cannot use `ischo`.

## Spectral bounds and window

| Name | Signature | Meaning |
|---|---|---|
| `SpectralBounds` | `SpectralBounds(lower=0.0, upper=0.0)` | Interval `[lower, upper]`. Attributes `lower`, `upper` are writable. |
| `GershgorinBounds` | `GershgorinBounds(matrix, pad=0.10)` → `SpectralBounds` | Symmetric enclosure `[−ρ(1+pad), ρ(1+pad)]` with ρ = maxᵢ Σⱼ \|Hᵢⱼ\| (ρ = 1 for a zero matrix). Always symmetric about 0, so an off-centre spectrum gets a wider window than needed. |
| `LanczosBounds` | `LanczosBounds(operator, iterations=40, seed=12345)` → `SpectralBounds` | Extreme Ritz values after `iterations` Lanczos steps from a random-phase start. They lie inside the true spectrum; add a margin before you use them as a window. |
| `SpectralWindow` | `SpectralWindow()` or `SpectralWindow(lower, upper)` | Window onto which H is rescaled. Writable `lower`, `upper` (eV); read-only `half_width` = (upper − lower)/2, `center` = (upper + lower)/2, `scale` = 1/half_width, `shift` = −center/half_width. |

### Rescaled Hamiltonians

| Class | Constructor | Meaning |
|---|---|---|
| `RescaledOperator` | `(matrix, window)` or `(operator, window)` | H̃ = scale·H + shift·I. Raises `ValueError` unless `window.upper > window.lower`. Property `window`. |
| `GeneralisedRescaledOperator` | `(h, s, inverse, window)` | H̃ = S⁻¹(scale·H + shift·S), the generalised route. `is_hermitian` is always `False`. Keeps `h`, `s`, `inverse` alive. Property `window`. |

## Trace states

### `StateGenerator` (base class)

| Member | Signature | Meaning |
|---|---|---|
| property | `dimension` → `int` | Length of each state. |
| property | `number_of_states` → `int` | Number of states in the sequence. |
| property | `label` → `str` | State-set label written to file names and headers. |
| method | `reset()` | Rewinds the sequence; random generators are re-seeded. |
| method | `next()` → `complex128[N]` or `None` | Copy of the next state, or `None` at the end. |

### Generator classes

| Class | Constructor | States | `label` |
|---|---|---|---|
| `RandomPhaseStates` | `(dim, num_states, seed=12345)` | `num_states` normalised random-phase vectors from `std::mt19937`. | `"default"` |
| `RangedRandomPhaseStates` | `(dim, num_states, seed, row_begin, row_end)` | Rows `[row_begin, row_end)` of the same random-phase vectors. | `"default"` |
| `LocalStates` | `(dim, sites)` | Unit vectors \|e_s⟩ for each index in `sites` (`int64` array or list). | `"exact"` if `sites` is 0…N−1, else `"local"` |
| `BlockLocalStates` | `(num_blocks, block_size)` or `(n1, n2, n3, block_size)` | W states; state r has 1/√N_k on orbital r of every block. An exact trace only when every operator of the word is block diagonal. | `"exact"` |
| `TransformedStates` | `(transform, inner)` | `transform·|r⟩` for each state of `inner`. Keeps both alive. | label of `inner` |
| `PermutedStates` | `(permutation, inner)` | States of `inner` with components reordered: new index i takes old index `permutation[i]`. Keeps `inner` alive. | label of `inner` |
| `UserStates` | `(states)` | The given list of 1D `complex128` arrays. | — |

`BlockLocalStates.block_index(m1, m2, m3, n1, n2)` (static) returns the flat
cell index m1 + n1·m2 + n1·n2·m3.

## Plans

A `Plan` is the bytecode that `words.compile` produces from an observable word.
The instructions act on the ket side until `MirrorToBra`, then on the bra side,
and `Trace` contracts the two.

| Class | Constructor | Attributes | Meaning |
|---|---|---|---|
| `ApplyOperator` | `(operator_index)` | `operator_index` | Replace each live vector v by O_k v, with O_k the `operator_index`-th entry of the engine's `operators` list. |
| `Expand` | `(num_moments)` | `num_moments` | Replace v by the sequence T₀(H̃)v … T_{M−1}(H̃)v. Each `Expand` adds one axis to the moment tensor; at most one per side. |
| `ApplySeries` | `(coeffs, scale=1+0j)` | `coeffs`, `scale` | Replace v by `scale`·Σₙ cₙ Tₙ(H̃) v (for example U(t)). Empty `coeffs` gives v = 0. Must come before `Expand` on its side. |
| `ApplyDerivedSeries` | `(coeffs, operator_index, scale=1+0j)` | `coeffs`, `operator_index`, `scale` | Replace v by `scale`·Σₙ cₙ Dₙ v, with D₀ = 0, D₁ = O_k T₀, D_{n+1} = 2H̃Dₙ − D_{n−1} + 2O_k Tₙ (the commutator [O_k, U(t)]). |
| `MirrorToBra` | `()` | — | Start the bra side from a fresh copy of the trace state. At most once. |
| `TimeAxis` | `(num_steps, evolve=[])` | `num_steps`, `evolve` | Add a time axis of `num_steps` slices; the steps in `evolve` run between slices. |
| `Trace` | `(weights)` | `weights` | Final instruction: μ[i_bra…, i_ket…] += w·⟨bra\|ket⟩ for every trace state. `weights` has length 1, the moment-axis length, or the number of tensor elements. |
| `Plan` | `()` or `(code=[], word="")` | `code`, `word` | Instruction list (must end with one `Trace`) and the word string stored in the moments. |

Steps of a `TimeAxis` (`Register` is an enumeration with `Ket`, `Bra`, `Aux`):

| Class | Constructor | Meaning |
|---|---|---|
| `EvolveRegister` | `(target, coeffs, scale=1+0j)` | seed[target] ← `scale`·Σₙ cₙ Tₙ(H̃) seed[target]. |
| `EvolveDerived` | `(target, source, coeffs, operator_index, scale=1+0j)` | seed[target] += `scale`·Σₙ cₙ Dₙ seed[source]. |
| `CopyRegister` | `(from_reg, to_reg)` | seed[to_reg] ← seed[from_reg]. |

All `coeffs` and `weights` accept a 1D `complex128` array or a list of
complex numbers.

## Moment engine

### `MomentEngine`

```text
MomentEngine(rescaled_h, operators, states, batch_size,
             state_block=1, left_window=32, right_window=0)
```

| Parameter | Type | Meaning |
|---|---|---|
| `rescaled_h` | `LinearOperator` | H̃, with spectrum in [−1, 1] (a `RescaledOperator` or `GeneralisedRescaledOperator`). |
| `operators` | list of `LinearOperator` | Operators referenced by index in the plan. |
| `states` | `StateGenerator` | Trace states. |
| `batch_size` | `int` > 0 | Chebyshev vectors kept per batch during `Expand`. If M is not a multiple of `batch_size`, the recursion runs to the next multiple and the extra moments are dropped. |
| `state_block` | `int` ≥ 1 | Trace states advanced together through one block product. The result is bit-identical to `state_block=1`. |
| `left_window` | `int` > 0 | Native builds, rank-2 words without time axis: bra vectors streamed per window. |
| `right_window` | `int` ≥ 0 | Same path: ket vectors stored per block; `0` stores all M (fastest, memory M·N·16 bytes). A value B < M stores B vectors and restarts the ket recursion per block (M + M²/B products per state instead of 2M). |

All operators and `states` must have the dimension of `rescaled_h`. The engine
keeps `rescaled_h`, `operators` and `states` alive.

| Method | Signature | Meaning |
|---|---|---|
| `run` | `run(plan, checkpoint_path="")` → `ChebyshevMoments` | Rewinds `states`, runs `plan` over every trace state and returns the tensor, whose rank is the number of `Expand` instructions. Fills `system_size`, `number_of_states`, `state_label`, `word`; you set `window`. Releases the GIL. With a non-empty `checkpoint_path`, it rewrites the partial tensor after each finished state (or state block) and resumes from a matching file, with a bit-identical result. |

> **Known limitation.** In a `golden` build `left_window` and `right_window`
> are validated but never used: the windowed rank-2 path is compiled only with
> `LSQUANT_NATIVE_BUILD` (`lsquant/core/chebyshev/moment_engine.cpp:1486–1492`).
> In a `native` build that path ignores `batch_size` for static rank-2 words.
> The defaults also differ between layers: `MomentEngine` uses
> `left_window=32, right_window=0`, while `lsquant.moments.run_plan` uses 16
> and 16 (`python/lsquant/moments.py:504–505`).

> **Known limitation.** A rank-2 plan with `Trace([1])` returns raw sums over
> the states. Call `ApplyKpmWeights` afterwards, as `run_plan` does
> (`python/lsquant/moments.py:597–598`).

### `ApplyKpmWeights(moments, num_states, symmetrise)` → `None`

Applies the rank-2 KPM factors in place: 4/`num_states`, times ½ for each
index equal to 0. With `symmetrise=True` it also sets μ(m, n) = conj μ(n, m)
from their average. Raises `ValueError` if the rank is not 2, `num_states` ≤ 0,
the tensor is already weighted, or `symmetrise` is requested on a non-square
tensor. Sets `kpm_weighted` and `symmetrised`.

### `ChebyshevMoments`

| Member | Type | Access | Meaning |
|---|---|---|---|
| constructor | `ChebyshevMoments(dims=[])` | — | Zero tensor of shape `dims`; every dimension > 0. |
| `rank` | `int` | read | Number of axes. |
| `dims` | `list[int]` | read | Shape. |
| `data` | `complex128` array | read | Row-major view of the tensor, without copy. Writing to it changes the moments. |
| `window` | `SpectralWindow` | read/write | Window used for H̃. |
| `system_size` | `int` | read/write | N. |
| `number_of_states` | `int` | read/write | Trace states used. |
| `state_label` | `str` | read/write | `"exact"`, `"default"`, `"local"`, … |
| `word` | `str` | read/write | Observable word. |
| `build_profile` | `str` | read | Profile of the build that created the object. |
| `kpm_weighted` | `bool` | read/write | KPM factors applied. |
| `symmetrised` | `bool` | read/write | Hermitian symmetrisation applied. |
| `time_step` | `float` | read/write | Time step of a time axis. |
| `number_of_time_steps` | `int` | read/write | Length of the time axis; 0 if none. |

## Moment files

| Function | Meaning |
|---|---|
| `ReadChebmom(path)` → `ChebyshevMoments` | Reads a unified `.chebmom` file (first token `chebmom`) or a legacy `.chebmom1D` / `.chebmom2D` file. For a legacy file it sets `kpm_weighted=True`, `symmetrised = (rank == 2)`, and takes `state_label` from the token after `_state` in the file name. |
| `WriteChebmom(path, moments)` | Writes the unified `.chebmom` version-1 text format with its metadata. |
| `WriteLegacyChebmom1D(path, moments)` | Rank-1 tensor in the legacy `.chebmom1D` layout. |
| `WriteLegacyChebmom2D(path, moments)` | Rank-2 tensor in the legacy `.chebmom2D` layout. |
| `WriteLegacyChebmomTD(path, moments)` | Rank-2 `(moments, times)` tensor in the legacy `.chebmomTD` layout. |

The [output-file reference](output-files.md) describes the formats.

## Contractions and curves

These functions evaluate Chebyshev series on a grid. They work on dimensionless
quantities in [−1, 1]; the caller applies kernels, tables and physical
prefactors (see `python/lsquant/reconstruct.py`). The sum over moments is
sequential for each grid point, so the result is the same for any number of
threads.

| Function | Result for each grid point i |
|---|---|
| `ContractRank1(moments, table)` | Σ_m table[m, i]·μ[m]. `moments`: 1D `float64` array, or a rank-1 `ChebyshevMoments` (its real part is used). `table`: `float64` of shape `(M, N)`. |
| `ContractRank2Imag(moments, table_a, table_b)` | Σ_{m,n} A[m, i]·Im(B[n, i]·μ[m, n]). `table_a` `float64`, `table_b` `complex128`. |
| `ContractRank2Real(moments, table_a, table_b)` | Σ_{m,n} Re(A[m, i]·B[n, i]·μ[m, n]). Both tables `complex128`. |
| `ContractRank2RealReal(moments, table_a, table_b)` | Σ_{m,n} A[m, i]·B[n, i]·Re μ[m, n]. Both tables `float64`. |
| `FftMomentContraction(moments, points)` | Returns a `(2, points)` `complex128` array: row 0 = Σ μ[m, n] cos(mθ) e^(−inθ), row 1 = Σ n·μ[m, n] cos(mθ) e^(−inθ), on θ_j = (2j + 0.5)π/`points`. `moments`: 2D `complex128`. |

For the rank-2 functions, `moments` is a 2D `complex128` array or a rank-2
`ChebyshevMoments`; `table_a` has shape `(M0, N)` and `table_b` `(M1, N)`.
Each returns a 1D `float64` array of length N.

| Function | Meaning |
|---|---|
| `WriteCurveText(path, x, y)` | Writes two 1D `float64` arrays of equal length as two columns, `x y` per line, in C++ default stream format (6 significant digits). |

## `compat` submodule

`_core.compat` holds kernels that reproduce the rounding of the legacy code.
The golden gate uses them; keep them out of new code.

| Function | Meaning |
|---|---|
| `compat.legacy_green_tables(x_grid, num_moments_0, num_moments_1=None)` → `(green, dgreen)` | Legacy Chebyshev tables of the retarded Green function and its derivative on `x_grid`, of shapes `(num_moments_0, len(x_grid))` and `(num_moments_1, len(x_grid))`. `num_moments_1` defaults to `num_moments_0`. |
| `compat.bessel_series_coefficients(x, tol=1e-15)` → `complex128` array | Coefficients 2(−i)ⁿ Jₙ(x) of the time-evolution series (first term J₀(x)), computed with GSL until \|Jₙ(x)\| ≤ `tol`. |

## Differences from the type stub

`python/lsquant/_core.pyi` differs from the bindings in these points:

| Binding (`lsquant/api/bindings.cpp`) | Stub (`python/lsquant/_core.pyi`) |
|---|---|
| `RangedRandomPhaseStates` (line 661) | missing |
| `BlockLocalStates(n1, n2, n3, block_size)` and `BlockLocalStates.block_index` (lines 678–681) | only `(num_blocks, block_size)` (line 268) |
| `synchronize_cuda()` in CUDA builds (line 465) | missing |
| `ChebyshevMoments.build_profile` is read-only (line 901) | declared as a writable attribute (line 427) |
| `GinkgoOperator`, `CudaAtomBlockOperator`, `CudaVectorBlock` exist only in builds with Ginkgo or CUDA (lines 418, 468) | declared unconditionally |
