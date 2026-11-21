# Orthogonal downfolding tools

Two standalone scripts that reduce a SIESTA Hamiltonian to an orthogonal tight-binding
model on a smaller basis. `python/lsquant/projection/origin.py` loads them; the LSQUANT
wrapper, the FFT and Wigner-Seitz reuse and the run-file interface live in
`python/lsquant/projection/`, and changes belong there.

| Script | Method |
|---|---|
| `downfold.py` | shell/polar and frozen-window reduction |
| `downfold_scdm.py` | SCDM (selected columns of the density matrix); imports `downfold.py` |

Both validate the reduced bands on deterministic off-grid k-points and on a symmetry path.
The literature for each method is cited in the source docstrings. `SHA256SUMS` records the
bytes of the two scripts.

Standalone usage (NumPy 2.2, SciPy 1.15, sisl 0.16):

```sh
python tools/downfold_origin/downfold.py RUN_DIRECTORY LABEL --out OUTPUT
python tools/downfold_origin/downfold_scdm.py RUN_DIRECTORY LABEL --out OUTPUT
```
