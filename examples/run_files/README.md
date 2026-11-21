# Run files: the Haldane model from the command line

This example drives a whole calculation with a JSON run file and the `lsquant` command
line, with no Python code of your own. The model is the golden Haldane model of
`goldens/stage0/haldane/operators`: 512 sites, a Chern insulator with C = +1.

`haldane_hall.json` asks for two observable words on the golden operators (M = 128, exact
trace, window [−4.5, 4.5] eV):

| Word | What it is | Reconstruction |
|---|---|---|
| `Tr[VX delta(E-H) VY dG-(E)/dE]` | the Kubo-Bastin Hall conductivity σ_xy | `bastin` |
| `Tr[VX delta(E-H) VX delta(E-H)]` | the Kubo-Greenwood longitudinal conductivity σ_xx | `greenwood` |

`run.sh` goes through five steps:
1. `words describe`: what the Bastin word compiles to;
2. `inspect`: the run file after validation, with every default filled in;
3. `run`: moments and curves at 10 meV;
4. `reconstruct`: σ_xy again at 150 meV from the stored VX–VY moments, without computing new
   moments;
5. `chern.py`: C = 2π⟨σ_xy⟩/A on the plateau |E| < 0.5 eV, for both broadenings.

## Run

```bash
export PYTHONPATH=$PWD/build/lsquant_core:$PWD/python
examples/run_files/run.sh            # outputs under out/run_files (or give a directory)
```

On a 64-core workstation (Ubuntu 22.04, g++ 11.4) the five steps take 12 s in total. The last
lines are:

```
== 5. Chern number from the plateau of sigma_xy at 10 and 150 meV
out/run_files/haldane/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: C = +1.0170
out/run_files/haldane_150meV/KuboBastin_NonEqOpVX-VYhaldaneKPM_M128x128_stateexactJACKSON.dat: C = +1.0306
checks: PASS
```

## Reading the result

- **Step 2** prints the resolved run file. The keys the run file leaves out appear with
  their defaults, for example `basis.name = "real"` and `basis.supercell = [1, 1, 1]`. The
  run writes the same content to `run.resolved.json` in its output directory.
- **Step 3** writes, for each word, the moments (`.chebmom2D`) and the curve (`.dat`) under
  their legacy names, because the run file sets `output.legacy_names = true`.
- **Step 4** reads the stored moments only. At 150 meV the reconstruction uses
  M_eff = ⌈π/η⌉ = 95 of the 128 moments, with η = 2 × 0.150/9 in rescaled units. At
  10 meV, ⌈π/η⌉ = 1414 exceeds M, so all 128 moments are used: every broadening up to
  110 meV gives the same curve, the one set by the resolution of 128 moments (πa/M = 0.11 eV).
- **Step 5.** C = +1.017 at 10 meV, the value of the golden model. At 150 meV the
  plateau is still quantised, C = +1.031, but the broadening rounds the edges of the gap at
  ±0.78 eV and lifts the average over |E| < 0.5 eV slightly. `chern.py` accepts |C − 1| < 0.05.

The words, the run-file sections and the commands are described in the manual page *Run
files, observable words, the CLI and the MCP server* (`docs/run-files.md`).
