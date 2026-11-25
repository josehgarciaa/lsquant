"""Read small Hall-response outputs without importing the calculation package.

The checks use eV and e^2/h (Pauli spin convention for spin Hall curves).
They compare tabulated observables, never moments or Hamiltonians. Work and
storage are linear in the curve length; each input is limited to 8 MB.
Warning: interpolation is permitted only inside the sampled energy interval.
Pitfall: topology is independent of a nonquantised spin Hall conductivity;
rounding that conductivity cannot establish a spin Chern number.
To-do: validate the stated numerical tolerances against additional builds.
"""
import json
from pathlib import Path

import numpy as np


class Checks:
    """Accumulate failures while displaying the numerical acceptance criteria."""

    def __init__(self):
        self.ok = True

    def near(self, label, value, expected, tolerance, reason):
        good = bool(np.isfinite(value) and abs(value - expected) <= tolerance)
        self.ok &= good
        print(f"{label}: {value:.9g}; expected {expected:.9g} +/- {tolerance:g}; "
              f"{reason}: {'PASS' if good else 'FAIL'}")

    def require(self, condition, reason):
        self.ok &= bool(condition)
        print(f"{reason}: {'PASS' if condition else 'FAIL'}")

    def finish(self):
        print('CHECK PASS' if self.ok else 'CHECK FAIL')
        return 0 if self.ok else 1


def small(path):
    path = Path(path)
    if path.stat().st_size > 8_000_000:
        raise ValueError(f"input exceeds 8 MB: {path}")
    return path


def read_json(path):
    with small(path).open() as stream:
        return json.load(stream)


def curve_arrays(energy, values):
    e, v = np.asarray(energy, dtype=float), np.asarray(values, dtype=float)
    if (e.ndim != 1 or v.shape != e.shape or len(e) < 2 or
            not np.all(np.isfinite(e)) or not np.all(np.isfinite(v))):
        raise ValueError('curve must contain at least two finite energy/value pairs')
    if np.all(np.diff(e) < 0):
        e, v = e[::-1], v[::-1]
    if not np.all(np.diff(e) > 0):
        raise ValueError('energy samples must be unique and strictly ordered')
    return e, v


def table(path):
    a = np.loadtxt(small(path), ndmin=2)
    if a.shape[1] != 2:
        raise ValueError(f"expected exactly two columns: {path}")
    return curve_arrays(a[:, 0], a[:, 1])


def oracle(path, field='values'):
    a = read_json(path)
    return curve_arrays(a['energies'], a[field])


def at(curve, energy):
    e, v = curve
    if np.any(np.asarray(energy) < e[0]) or np.any(np.asarray(energy) > e[-1]):
        raise ValueError(f"requested energy outside [{e[0]}, {e[-1]}] eV")
    return np.interp(energy, e, v)


def plateau(curve, low, high):
    """Average samples in the central half of a documented gap, as in PAGE."""
    width = (high - low) / 4
    e, v = curve
    at(curve, [low + width, high - width])
    mask = (e >= low + width) & (e <= high - width)
    if mask.sum() < 3:
        raise ValueError('fewer than three plateau samples')
    return float(np.mean(v[mask]))


def compare(c, label, produced, reference, tolerance=5e-5):
    """Compare every reference sample; reject incomplete energy coverage.

    Text curves round energy to six decimals and conductivity to about seven
    significant digits. The absolute tolerance also allows interpolation of
    those rounded grids; it is not a tolerance for changing a physical model.
    """
    error = float(np.max(np.abs(at(produced, reference[0]) - reference[1])))
    c.near(label + ' max curve error', error, 0, tolerance,
           'text rounding and interpolation of rounded energy grids')


def topology(c, path, *, nontrivial=True, occupied=104, gap=0.0766):
    a = read_json(path)
    target = 1 if nontrivial else 0
    for key, value in [('chern_plus', target), ('chern_minus', -target),
                       ('spin_chern', 2 * target), ('z2', target),
                       ('occupied_bands', occupied)]:
        c.near(key, a[key], value, 1e-6 if 'chern' in key else 0,
               'integer topology; numerical FHS roundoff only')
    c.near('energy gap (eV)', a['energy_gap'], gap, 0.003,
           'rounded calibration and finite mesh band edges')
    c.require(np.isfinite(a['projected_spin_gap']) and a['projected_spin_gap'] > 1e-7,
              'projected spin gap open (> 1e-7)')
    return a


def guarded(main):
    """Make absent, malformed, and ambiguous inputs explicit failed checks."""
    try:
        return main()
    except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
        print(f"input error: {error}\nCHECK FAIL")
        return 1


def spin_example(material):
    """Apply mesh-specific calibration and a production plateau gate.

    The quick runs use different meshes: 24x24 for Bi, 12x12 for PdSe2 (their
    run files). The PdSe2 quick value at E_F is small because E_F = 0 lies in
    the trivial gap; its calibration (0.0105 e^2/h) was measured on the
    reference build.
    """
    import argparse
    p = argparse.ArgumentParser(description='Check completed Hall outputs')
    p.add_argument('output', type=Path)
    mode = p.add_mutually_exclusive_group()
    mode.add_argument('--quick', action='store_true')
    mode.add_argument('--production', action='store_true')
    args = p.parse_args()
    d = args.output
    r = read_json(d / 'run.resolved.json')
    read_json(d / 'basis.json')
    m = r['expansion']['moments']
    if m not in (256, 1024):
        raise ValueError('only M=256 quick and M=1024 production are calibrated')
    quick = m == 256
    if (args.quick and not quick) or (args.production and quick):
        raise ValueError('requested mode disagrees with resolved moments')
    c = Checks()
    mesh = r['basis']['supercell']
    c.require(r['model']['label'] == material, 'resolved model matches example')
    quick_mesh = [24, 24, 1] if material == 'bi' else [12, 12, 1]
    c.require(mesh == (quick_mesh if quick else [48, 48, 1]),
              f'mesh matches the calibration ({quick_mesh[0]}x{quick_mesh[1]} quick, 48x48 production)')
    c.require(r['model']['params']['velocity_mode'] == 'covariant' and
              r['metric']['route'] == 'lowdin', 'covariant velocity and lowdin metric')
    bi = material == 'bi'
    a = topology(c, d / 'spin_chern.json', nontrivial=bi,
                 occupied=30 if bi else 100,
                 gap=(0.3044 if quick else 0.30325) if bi else 1.1746)
    ef = a['fermi_energy_eV']
    c.near('Fermi energy (eV)', ef, 0, 0, 'run files specify E_F=0')
    exact = oracle(d / 'exact_shc.json')
    c.near('oracle sigma(E_F) vs topology record', float(at(exact, ef)),
           a['sigma_xy'], 0.005, 'Gaussian curve interpolation at E_F')
    if bi:
        c.near('exact sigma(E_F) (e^2/h)', a['sigma_xy'],
               1.144 if quick else 1.1747, 0.005, 'rounded reference calibration (covariant velocity)')
        if not quick:
            c.near('projected spin gap', a['projected_spin_gap'], 0.5131, 0.005,
                   'rounded production calibration')
    elif quick:
        c.near('exact sigma(E_F) (e^2/h)', a['sigma_xy'], 0.0105, 0.005,
               'trivial gap at E_F = 0; value measured on the 12x12 quick mesh')
    else:
        print('Independent PdSe2 production sigma(E_F=0) calibration: not checked; '
              '1.987 describes a different energy and projection')
    if quick:
        print('Finite-M reconstruction: not checked; M=256 is not expected to reach the plateau')
    else:
        reconstructed = table(d / 'shc.dat')
        resolutions = r['receipt']['curve_resolution']
        if not resolutions:
            raise ValueError('missing production resolution diagnostics')
        for entry in resolutions:
            widths = entry['effective_broadening_meV']
            c.require(bool(widths) and all(np.isfinite(x) and 0 < x < 1000 * a['energy_gap']
                                          for x in widths),
                      'production effective resolution below the insulating gap')
        # Reported, not asserted. The page of example 07 states that the converged
        # reconstructed value of the production setting is still under study; measured on
        # the reference build (48x48, M = 1024, covariant, lowdin): reconstructed 1.367 at
        # E_F against the exact 1.175. A gate here would assert something the page does not
        # claim. To-do: turn this into an assertion once the difference is understood.
        points = np.linspace(-0.02, 0.02, 5)
        error = float(np.max(np.abs(at(reconstructed, points) - at(exact, points))))
        print(f'production reconstructed sigma(E_F): {float(at(reconstructed, ef)):.6g} e^2/h; '
              f'exact: {float(at(exact, ef)):.6g}; largest residual in the central 40 meV: '
              f'{error:.3g}: reported, not checked (the page states it is under study)')
    return c.finish()


def routes_example(number):
    """Check page plateaus and explicitly supplied SIESTA results.

    The page mean uses the central half-gap on each code's own energy axis.
    Warning: a missing route is reported as not checked, never as passed.
    Projected runs do not currently emit topology: supplying such a run still
    requires spin_chern.json and fails clearly if it is absent.
    """
    import argparse
    root = Path(__file__).resolve().parent
    ex20 = root / '20_pdse2_spin_hall_routes' / 'data'
    ex21 = root / '21_pdse2_downfolding' / 'data'
    p = argparse.ArgumentParser(description='Check PdSe2 curves and topology')
    p.add_argument('--data-only', action='store_true')
    for route in ('siesta-direct', 'siesta-projected', 'siesta-reduced12'):
        p.add_argument('--' + route, type=Path)
    args = p.parse_args()
    outputs = {name: getattr(args, name.replace('-', '_')) for name in
               ('siesta-direct', 'siesta-projected', 'siesta-reduced12')}
    if args.data_only and any(outputs.values()):
        raise ValueError('--data-only cannot be combined with output directories')
    if not args.data_only and not any(outputs.values()):
        raise ValueError('give an output directory or explicitly select --data-only')
    c = Checks()
    # Constants and filenames follow plot.py, without importing plotting code.
    specifications = [
        (ex20, 'siesta_direct_exact.dat', -0.972, (1.986, 2.061), 1.98, 0.01),
        (ex20, 'siesta_projected_kpm.dat', -0.972, (1.986, 2.061), 2.01, -0.01),
        (ex20, 'wannier_kpm.dat', -0.822, (2.009, 2.055), 2.02, -0.07),
        (ex20, 'paoflow_kpm.dat', -0.827, (2.009, 2.055), 2.00, 0.02),
    ]
    if number == 20:
        specifications += [
            (ex20, 'wannier_covariant_spn.dat', -0.822, (2.009, 2.055), 1.94, None),
            (ex20, 'paoflow.dat', -0.827, (2.009, 2.055), 1.94, None),
            (ex20, 'openmx_exact.dat', -0.606, (1.918, 1.982), 1.97, None),
        ]
    else:
        specifications += [(ex21, 'siesta_reduced12_kpm.dat', -0.972,
                            (1.986, 2.061), 1.98, -0.06)]
    for directory, filename, vbm, gap, expected, trivial in specifications:
        curve = table(directory / filename)
        c.near(filename + ' plateau (e^2/h)',
               plateau(curve, gap[0] + vbm, gap[1] + vbm), expected, 0.015,
               'PAGE rounded to two decimals; discrete central-half-gap mean')
        if number == 21 and trivial is not None:
            c.near(filename + ' trivial gap (e^2/h)', float(at(curve, vbm + 0.5)),
                   trivial, 0.015, 'PAGE rounded to two decimals at E-E_VBM=0.5 eV')
    errors = []
    for route, output in outputs.items():
        if output is None:
            print(route + ': not checked (no output directory supplied)')
            continue
        try:
            read_json(output / 'run.resolved.json')
            if route == 'siesta-direct':
                compare(c, route + ' oracle', oracle(output / 'exact_shc.json'),
                        table(ex20 / 'siesta_direct_exact.dat'))
                # No committed direct M=64 KPM curve exists. Never substitute
                # the projected M=1024 curve, which describes another run.
                direct = table(output / 'shc.dat')
                print(f'direct reconstructed curve: {len(direct[0])} finite samples; '
                      'reference comparison not checked (no committed same-route KPM curve)')
            else:
                stem = route.replace('-', '_')
                directory = ex21 if route == 'siesta-reduced12' else ex20
                compare(c, route + ' projected oracle',
                        oracle(output / 'word_0_oracle.json', 'projected'),
                        table(directory / (stem + '_exact.dat')))
                compare(c, route + ' original oracle',
                        oracle(output / 'word_0_oracle.json', 'unprojected'),
                        table(ex20 / 'siesta_direct_exact.dat'))
                compare(c, route + ' reconstruction', table(output / 'word_0.dat'),
                        table(directory / (stem + '_kpm.dat')))
            a = topology(c, output / 'spin_chern.json')
            c.near(route + ' E_F (eV)', a['fermi_energy_eV'], 1.052, 1e-9,
                   'upper-gap calibration')
            if route == 'siesta-direct':
                c.near('direct exact sigma(E_F)', a['sigma_xy'], 1.9791, 0.005,
                       'rounded reference calibration')
        except (OSError, ValueError, KeyError, TypeError, IndexError) as error:
            errors.append(f'{route}: {error}')
    for error in errors:
        c.require(False, error)
    for route in ('Wannier90', 'PAOFLOW', 'OpenMX'):
        print(route + ' new calculation: not checked (committed data checked separately)')
    return c.finish()
