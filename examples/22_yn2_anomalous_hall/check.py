#!/usr/bin/env python3
"""Check YN2 Hall plateaus, band gaps, and optional completed AHC outputs.

Conductivities are in e^2/h and energies in eV from each model's E_F.
LSQuant's signed C is -3; the invariant magnitude is 3. The gap from the
committed path bands is a path gap, not a proof of a global insulating gap.
Warning: no electronic-structure calculation is attempted here.
Pitfall: cross-code curves agree in the plateau, not at every band energy.
Work is linear in the small tables and bands; files are limited to 8 MB.
To-do: establish a global gap from a full-mesh band table when available.
"""
import argparse
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from heavy_checks import Checks, at, compare, guarded, small, table


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    p.add_argument('output', nargs='?', type=Path)
    p.add_argument('--data-only', action='store_true')
    args = p.parse_args()
    if args.data_only and args.output:
        raise ValueError('--data-only cannot be combined with output')
    data = Path(__file__).resolve().parent / 'data'
    c = Checks()
    curves = {key: table(data / (key + '.dat')) for key in
              ('siesta_exact', 'wannier_exact', 'wannier_kpm_ahc',
               'wannier_kpm_sea', 'wannier_kpm_surf')}
    for key, expected, tol in [('siesta_exact', -3.01, 0.02),
                               ('wannier_exact', -3.00, 0.02),
                               ('wannier_kpm_ahc', -3.13, 0.02)]:
        c.near(key + ' sigma(E_F)', float(at(curves[key], 0)), expected, tol,
               'PAGE rounding and interpolation to E_F')
        c.near(key + ' signed C plateau', float(at(curves[key], 0)), -3, 0.16,
               'exact mesh error; KPM finite resolution allowance about 5%')
    points = np.linspace(-0.005, 0.005, 11)
    c.near('exact inter-code plateau max difference',
           float(np.max(np.abs(at(curves['siesta_exact'], points) -
                               at(curves['wannier_exact'], points)))), 0, 0.08,
           'different electronic structures and finite 300x300 meshes')
    total = curves['wannier_kpm_ahc']
    c.near('max |sea+surface-total|',
           float(np.max(np.abs(at(curves['wannier_kpm_sea'], total[0]) +
                               at(curves['wannier_kpm_surf'], total[0]) - total[1]))),
           0, 0.01, 'the total is integrated with a rectangle rule, the two parts are '
           'evaluated separately: they agree to the quadrature error (measured 0.0047)')
    c.near('surface part at E_F', float(at(curves['wannier_kpm_surf'], 0)),
           0, 0.02, 'finite Jackson resolution in the gap')
    with np.load(small(data / 'bands.npz'), allow_pickle=False) as bands:
        for code in ('siesta', 'wannier'):
            energy = np.asarray(bands['E_' + code])
            if energy.ndim != 2 or not np.all(np.isfinite(energy)):
                raise ValueError('invalid path band energies')
            # Bands are stored as (path points, band index), as plot.py uses them.
            if energy.shape[0] != len(bands['x_' + code]):
                raise ValueError('path band shape does not match path coordinates')
            occupied = np.sum(energy < 0, axis=1)
            if not np.all(occupied == occupied[0]) or occupied[0] == 0 or occupied[0] == energy.shape[1]:
                raise ValueError('path does not have a common occupied-band count')
            ordered = np.sort(energy, axis=1)
            n = occupied[0]
            gap = float(np.min(ordered[:, n]) - np.max(ordered[:, n - 1]))
            c.near(code + ' path gap (eV)', gap, 0.036 if code == 'siesta' else 0.025,
                   0.015, 'path sampling; 36 meV global SIESTA and about 25 meV plateau statements')
    if args.output:
        for produced, reference in [('ahc.dat', 'wannier_kpm_ahc'),
                                     ('ahc_sea.dat', 'wannier_kpm_sea'),
                                     ('ahc_surf.dat', 'wannier_kpm_surf')]:
            compare(c, produced, table(args.output / produced), curves[reference])
    else:
        print('new AHC output: not checked (no output directory supplied)')
    return c.finish()


if __name__ == '__main__':
    sys.exit(guarded(main))
