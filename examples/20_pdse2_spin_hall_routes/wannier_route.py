"""Spin Hall conductivity of a Wannier90 model with the exact Wannier spin.

Reads a finished Wannier90 run (``<seed>_hr.dat``, ``<seed>_r.dat``,
``<seed>.spn``, ``<seed>_u.mat``, ``<seed>_u_dis.mat``, ``<seed>.eig`` and
``<seed>.win``) and writes sigma^z_xy(E) in e^2/h with Pauli sigma_z.

The spin operator is projected from the DFT bands, S_W(k) = V^+ S V with
V = U_dis U restricted to the bands inside the outer window at each k. A
diagonal sigma_z on the Wannier spinor labels is not used: with spin-orbit
coupling the Wannier functions are spin-mixed, and that diagonal operator
misses part of the spin Hall response.

Usage:
    python wannier_route.py W90_DIR SEED E_FERMI OUT_PREFIX [--mesh 300]
"""
import argparse
import json
import os
import sys

import numpy as np
from scipy.io import FortranFile

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "python"))
from lsquant.models.wannier import load_wannier_model  # noqa: E402
from lsquant.oracle import spin_hall  # noqa: E402


def read_win(path):
    """Return the cell (rows in Angstrom), num_wann and the outer window."""
    text = open(path, encoding="utf-8").read().splitlines()
    keys, cell = {}, None
    for i, line in enumerate(text):
        low = line.split("!")[0].strip().lower()
        if low.startswith("begin unit_cell_cart"):
            rows = [ln.split() for ln in text[i + 1:i + 5] if ln.strip()]
            rows = [r for r in rows if len(r) == 3]
            cell = np.array(rows[:3], float)
        elif "=" in low:
            k, v = (s.strip() for s in low.split("=", 1))
            keys[k] = v
    return cell, int(keys["num_wann"]), float(keys["dis_win_min"]), float(keys["dis_win_max"])


def write_lsquant_geometry(w90_dir, seed, cell):
    """Write ``<seed>.uc`` and ``<seed>.xyz`` (Wannier centres) for LSQuant."""
    centres = []
    for line in open(os.path.join(w90_dir, seed + "_centres.xyz"), encoding="utf-8"):
        parts = line.split()
        if len(parts) == 4 and parts[0] == "X":
            centres.append([float(x) for x in parts[1:]])
    np.savetxt(os.path.join(w90_dir, seed + ".uc"), cell, fmt="%.14f")
    with open(os.path.join(w90_dir, seed + ".xyz"), "w", encoding="utf-8") as fh:
        fh.write(f"{len(centres)}\n")
        for i, c in enumerate(centres):
            fh.write(f"WF{i + 1} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}\n")


def read_spn_z(path):
    """Read the z component of ``<psi_mk|sigma|psi_nk>`` from an unformatted .spn."""
    f = FortranFile(path, "r")
    f.read_record("S60")
    nb, nk = f.read_ints(np.int32)
    S = np.zeros((nk, nb, nb), complex)
    for k in range(nk):
        # Upper triangle, m <= n, n fastest outside; three components per element.
        rec = f.read_reals(np.float64).view(np.complex128).reshape(-1, 3)
        c = 0
        for n in range(nb):
            for m in range(n + 1):
                S[k, m, n] = rec[c, 2]
                S[k, n, m] = np.conj(rec[c, 2])
                c += 1
    return S


def read_u(path, nk, n1, n2):
    """Read ``_u.mat`` or ``_u_dis.mat``: k points and (nk, n1, n2) matrices."""
    lines = open(path, encoding="utf-8").read().splitlines()
    i, kpts, mats = 2, [], []
    for _ in range(nk):
        while not lines[i].strip():
            i += 1
        kpts.append([float(x) for x in lines[i].split()])
        i += 1
        v = np.array([[float(x) for x in lines[i + j].split()] for j in range(n1 * n2)])
        i += n1 * n2
        mats.append((v[:, 0] + 1j * v[:, 1]).reshape(n1, n2, order="F"))
    return np.array(kpts), np.array(mats)


def wannier_spin(w90_dir, seed, num_wann, win_min, win_max):
    """Return the k points (fractional) and S_W(k) of the Wannier90 mesh."""
    S = read_spn_z(os.path.join(w90_dir, seed + ".spn"))
    nk, nb, _ = S.shape
    eig = np.loadtxt(os.path.join(w90_dir, seed + ".eig")).reshape(nk, nb, 3)[:, :, 2]
    kpts, U = read_u(os.path.join(w90_dir, seed + "_u.mat"), nk, num_wann, num_wann)
    _, Ud = read_u(os.path.join(w90_dir, seed + "_u_dis.mat"), nk, nb, num_wann)
    V = np.zeros((nk, nb, num_wann), complex)
    for k in range(nk):
        # _u_dis rows list only the bands inside the outer window, packed at the top.
        inside = np.where((eig[k] >= win_min) & (eig[k] <= win_max))[0]
        V[k, inside, :] = Ud[k, :len(inside), :] @ U[k]
    err = np.abs(np.einsum("kba,kbc->kac", V.conj(), V) - np.eye(num_wann)).max()
    if err > 1e-6:
        raise RuntimeError(f"V^+V is not the identity (error {err:.1e}); check the window")
    return kpts, np.einsum("kba,kbc,kcd->kad", V.conj(), S, V)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("w90_dir")
    ap.add_argument("seed")
    ap.add_argument("e_fermi", type=float, help="DFT Fermi energy in eV")
    ap.add_argument("out_prefix")
    ap.add_argument("--mesh", type=int, default=300, help="N for an N x N x 1 k mesh")
    ap.add_argument("--broadening", type=float, default=0.006, help="Gaussian width, eV")
    ap.add_argument("--compare-operators", action="store_true",
                    help="also print the plateau for three velocities and two spin operators")
    ap.add_argument("--eta", type=float, default=0.0,
                    help="interband regulator in eV: denominators (E_n - E_m)^2 + eta^2, as in "
                         "openmx_exact.py; 0 is the strict sum")
    ap.add_argument("--filling", type=int, default=4,
                    help="occupied Wannier bands at E_F (4 for the 12-function model); "
                         "the quantum spin Hall gap is 4 bands higher")
    ap.add_argument("--grid", type=float, nargs=3, default=[-0.5, 2.0, 1251],
                    help="E - E_F grid: start stop points")
    args = ap.parse_args()

    cell, nw, wmin, wmax = read_win(os.path.join(args.w90_dir, args.seed + ".win"))
    write_lsquant_geometry(args.w90_dir, args.seed, cell)
    model = load_wannier_model(args.w90_dir, args.seed)
    recip = 2 * np.pi * np.linalg.inv(cell).T
    area = abs(np.linalg.det(cell[:2, :2]))

    kpts, SWk = wannier_spin(args.w90_dir, args.seed, nw, wmin, wmax)
    # Fourier transform to the cell gauge on the R range of the n x n Wannier90 mesh.
    n = int(round(np.sqrt(len(kpts))))
    Rs = np.array([[a, b, 0] for a in range(-n // 2, n - n // 2) for b in range(-n // 2, n - n // 2)])
    SR = np.einsum("kr,kab->rab", np.exp(-2j * np.pi * (kpts @ Rs.T)), SWk) / len(kpts)

    def spin_atom_gauge(kc):
        kf = kc @ cell.T / (2 * np.pi)
        Sc = np.einsum("kr,rab->kab", np.exp(2j * np.pi * (kf @ Rs.T)), SR)
        p = np.exp(1j * (kc @ model.tau.T))
        Sa = p.conj()[:, :, None] * Sc * p[:, None, :]
        return 0.5 * (Sa + np.conj(np.swapaxes(Sa, 1, 2)))

    N = args.mesh
    g = (np.arange(N) + 0.5) / N
    kc = np.array([[i, j, 0.0] for j in g for i in g]) @ recip
    E, C = np.linalg.eigh(model.Hk(kc))
    Sz = spin_atom_gauge(kc)
    VX = model.velocity_blocks(kc, "x", mode="covariant")
    VY = model.velocity_blocks(kc, "y", mode="covariant")
    if args.eta > 0:
        omega = np.empty_like(E)
        for c in np.array_split(np.arange(len(E)), max(1, len(E) // 2000)):
            Cc, Cd = C[c], np.conj(np.swapaxes(C[c], 1, 2))
            J = Cd @ (0.5 * (VX[c] @ Sz[c] + Sz[c] @ VX[c])) @ Cc
            V = Cd @ VY[c] @ Cc
            dE = E[c][:, :, None] - E[c][:, None, :]
            omega[c] = -2 * np.sum(np.where(np.abs(dE) > 1e-6, np.imag(J * np.swapaxes(V, 1, 2))
                                            / (dE**2 + args.eta**2), 0), axis=2)
    else:
        omega = spin_hall.spin_berry_curvature(E, C, 0.5 * (VX @ Sz + Sz @ VX), VY)

    energies = np.linspace(args.grid[0], args.grid[1], int(args.grid[2]))
    sigma = np.concatenate([
        spin_hall.sigma_xy(E, omega, args.e_fermi + chunk, area, broadening=args.broadening)
        for chunk in np.array_split(energies, max(1, len(energies) // 25))])
    np.savetxt(args.out_prefix + ".dat", np.column_stack([energies, sigma]),
               header="E-E_F (eV)  sigma^z_xy (e^2/h, Pauli sigma_z); Wannier90 + LSQuant, "
                      f"covariant velocity, exact Wannier spin, {N}x{N}, interband regulator {args.eta * 1e3:g} meV, "
                      f"Gaussian {args.broadening * 1e3:g} meV")

    # Topology of the two gaps: the gap at E_F and the one 4 bands higher (about 2 eV above the VBM).
    Nc = 60
    gc = np.arange(Nc) / Nc
    kc2 = np.array([[i, j, 0.0] for j in gc for i in gc]) @ recip
    E2, C2 = np.linalg.eigh(model.Hk(kc2))
    S2 = spin_atom_gauge(kc2)
    report = {}
    for name, nocc in [("gap at E_F", args.filling), ("upper gap", args.filling + 4)]:
        lo, hi = E2[:, nocc - 1].max(), E2[:, nocc].min()
        r = spin_hall.spin_chern(E2, C2, np.broadcast_to(np.eye(nw), C2.shape), S2,
                                 0.5 * (lo + hi), (Nc, Nc, 1), k_points=kc2, centres=model.tau)
        mid = 0.5 * (lo + hi) - args.e_fermi
        report[name] = dict(gap_eV=[lo - args.e_fermi, hi - args.e_fermi],
                            spin_chern=float(r["spin_chern"]), z2=int(r["z2"]),
                            sigma_mid=float(np.interp(mid, energies, sigma)))
        print(f"{name}: E-E_F {lo - args.e_fermi:.3f} to {hi - args.e_fermi:.3f} eV, "
              f"spin Chern {r['spin_chern']:.3f}, Z2 {r['z2']}, "
              f"sigma {report[name]['sigma_mid']:.4f} e^2/h")
    if args.compare_operators:
        # Diagonal sigma_z on the spinor labels: Wannier90 orders each spinor
        # projection as (up, down).
        diag = np.broadcast_to(np.kron(np.eye(nw // 2), np.diag([1.0, -1.0])).astype(complex), Sz.shape)
        mid = [args.e_fermi + 0.5 * sum(report["upper gap"]["gap_eV"])]
        report["operators_at_upper_gap"] = {}
        for mode in ("bare", "berry_connection", "covariant"):
            vx = model.velocity_blocks(kc, "x", mode=mode)
            vy = model.velocity_blocks(kc, "y", mode=mode)
            for name, S in (("diagonal sigma_z", diag), ("exact Wannier spin", Sz)):
                om = spin_hall.spin_berry_curvature(E, C, 0.5 * (vx @ S + S @ vx), vy)
                val = float(spin_hall.sigma_xy(E, om, mid, area)[0])
                report["operators_at_upper_gap"][f"{mode}, {name}"] = val
                print(f"{mode:17s} {name:19s} {val:.4f} e^2/h")
    json.dump(report, open(args.out_prefix + ".json", "w"), indent=1)


if __name__ == "__main__":
    main()
