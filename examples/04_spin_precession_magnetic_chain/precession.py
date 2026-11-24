"""Analytic Larmor gate shared by plotting and the standalone check.

The exact trace uses dimensionless spin normalized to SX at t=0.
Only the m=0 time moments are needed, so the cost is O(N_time).
Pitfall: the published 96-step, 95 fs grid and J=0.1 eV are assumed.
The 1e-10 budget allows round-off; the reference error is 4e-15.
"""
import glob
import os
import numpy as np

NTIME = 96
TMAX = 95.0
HBAR = 0.6582119624
J_EX = 0.1


def _m0(path: str) -> tuple[np.ndarray, int]:
    """Read the m=0 row (energy-integrated spin expectation) from a .chebmomTD file."""
    with open(path, encoding="utf-8") as f:
        f.readline()
        nt = int(f.readline().split()[1])
        return np.array([float(f.readline().split()[0]) for _ in range(nt)]), nt


def check_precession(run_dir: str, tol: float = 1e-10) -> bool:
    """Compare the KPM spin with the analytic Larmor precession and print the result.

    The precession frequency is fitted from the unwrapped phase of
    <S_x> - i<S_y> and compared with 2 J_ex / hbar. The normalised spin
    components are compared point by point with cos(wt), -sin(wt) and 0.
    """
    series = {}
    for op in ("SX", "SY", "SZ"):
        mom = sorted(
            glob.glob(os.path.join(run_dir, f"TimeProj{op}-PX*stateexact*chebmomTD"))
        )[-1]
        series[op], nt = _m0(mom)
    t = np.arange(nt) * (TMAX / (NTIME - 1))
    norm = series["SX"][0]
    sx, sy, sz = series["SX"] / norm, series["SY"] / norm, series["SZ"] / norm

    omega = 2.0 * J_EX / HBAR
    phase = np.unwrap(np.angle(sx - 1j * sy))
    omega_fit = np.polyfit(t, phase, 1)[0]
    rel_omega = abs(omega_fit - omega) / omega
    dev_x = np.max(np.abs(sx - np.cos(omega * t)))
    dev_y = np.max(np.abs(sy + np.sin(omega * t)))
    dev_z = np.max(np.abs(sz))
    amp = np.max(np.abs(np.hypot(sx, sy) - 1.0))

    print(f"[1] fitted omega = {omega_fit:.12f} rad/fs, 2 J_ex/hbar = {omega:.12f} rad/fs, "
          f"relative difference {rel_omega:.2e}")
    print(f"[2] max |<S_x> - cos wt| = {dev_x:.2e}, max |<S_y> + sin wt| = {dev_y:.2e}, "
          f"max |<S_z>| = {dev_z:.2e}")
    print(f"[3] max |sqrt(<S_x>^2 + <S_y>^2) - 1| = {amp:.2e} over {nt} time steps")
    ok = max(rel_omega, dev_x, dev_y, dev_z, amp) < tol
    print("CHECK PASS" if ok else "CHECK FAIL")
    return ok
