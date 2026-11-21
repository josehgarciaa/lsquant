#!/usr/bin/env python3
"""Checks for the layer projectors and the [X,U(t)] mean-square displacement.

Layer DOS (run_layer_dos.json, L = 256, exact trace):
  1. Tr[A delta] + Tr[B delta] = Tr[delta]: the two layer projectors resolve the
     identity, so the layer DOS add up to the total DOS (to the 6 significant
     digits of the curve files).
  2. Each layer holds L states:  int rho_A dE = int rho_B dE = L.
  3. The layer bias +/-Delta/2 sets the mean energy of each layer exactly:
     <E>_A = Tr[P_A H]/L = +Delta/2,  <E>_B = -Delta/2.

Mean-square displacement (run_msd.json, L = 4096, random-phase trace):
  The MSD curve LSQuant writes is rho(E) DeltaX^2(E,t) per orbital, so
  DeltaX^2(E,t) = curve / (DOS/N). In the clean bilayer chain the motion is
  ballistic, DeltaX^2 = <v^2>_E t^2 / hbar^2, with both bands having
  v(k) = 2 t a sin(ka) and, at energy E, |v_b| = sqrt(4t^2 - (E -/+ s)^2),
  s = sqrt(t_p^2 + Delta^2/4); the DOS-weighted average is
  <v^2>_E = sum_b |v_b| / sum_b |v_b|^-1.
  4. The log-log slope of DeltaX^2(t) is 2 (ballistic).
  5. DeltaX^2/t^2 matches <v^2>_E/hbar^2 (t in fs, hbar = 0.6582 eV fs) within 3 %.

Usage:
    python check_spatial.py <out_root>    # holds layer_dos/ and msd/
"""
import glob
import os
import sys

import numpy as np

T, T_PERP, DELTA = 1.0, 0.5, 0.6
L_LAYER, L_MSD = 256, 4096
HBAR = 0.6582119569  # eV fs


def curve(pattern):
    return np.loadtxt(sorted(glob.glob(pattern))[0])


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "out"
    ok = True
    lay = os.path.join(root, "layer_dos")
    tot = curve(os.path.join(lay, "meanSpectralOp1bilayer_L*JACKSON.dat"))
    pa = curve(os.path.join(lay, "meanSpectralOpAbilayer_L*JACKSON.dat"))
    pb = curve(os.path.join(lay, "meanSpectralOpBbilayer_L*JACKSON.dat"))
    dsum = float(np.max(np.abs(pa[:, 1] + pb[:, 1] - tot[:, 1])) / tot[:, 1].max())
    good = dsum < 1e-5         # curve files carry 6 significant digits
    ok &= good
    print(f"[1] max |rho_A + rho_B - rho| / max rho = {dsum:.2e}  {'PASS' if good else 'FAIL'}")
    for name, c, sign in (("A (top)", pa, +1), ("B (bottom)", pb, -1)):
        w = float(np.trapezoid(c[:, 1], c[:, 0]))
        e_mean = float(np.trapezoid(c[:, 0] * c[:, 1], c[:, 0]) / w)
        good = abs(w - L_LAYER) < 1e-3 * L_LAYER and abs(e_mean - sign * DELTA / 2) < 5e-3
        ok &= good
        print(f"[2,3] layer {name}: weight = {w:.2f} (L = {L_LAYER}), <E> = {e_mean:+.4f} eV "
              f"(expected {sign * DELTA / 2:+.3f})  {'PASS' if good else 'FAIL'}")
    msd_dir = os.path.join(root, "msd")
    dos = curve(os.path.join(msd_dir, "meanSpectralOp1bilayer_L*JACKSON.dat"))
    s = np.hypot(T_PERP, DELTA / 2)
    for ef in (0.0, 1.0):
        m = curve(os.path.join(msd_dir, f"meanCorrelationVX*EF{ef:.6f}JACKSON.dat"))
        rho = np.interp(ef, dos[:, 0], dos[:, 1]) / (2 * L_MSD)
        dx2 = m[:, 1] / rho
        t = m[:, 0]
        slope = float(np.polyfit(np.log(t[1:]), np.log(dx2[1:]), 1)[0])
        v = np.sqrt([4 * T**2 - (ef - s) ** 2, 4 * T**2 - (ef + s) ** 2])
        v2 = float(v.sum() / (1 / v).sum())
        ratio = float(np.mean(dx2[1:] / t[1:] ** 2))
        good = abs(slope - 2) < 1e-3 and abs(ratio / (v2 / HBAR**2) - 1) < 0.03
        ok &= good
        print(f"[4,5] E_F = {ef:.1f} eV: slope = {slope:.5f}, DeltaX^2/t^2 = {ratio:.3f} A^2/fs^2, "
              f"<v^2>/hbar^2 = {v2 / HBAR**2:.3f}; DeltaX^2(t=40 fs) = {dx2[-1]:.0f} A^2  "
              f"{'PASS' if good else 'FAIL'}")
    print("CHECK", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
