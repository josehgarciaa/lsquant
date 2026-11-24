"""Exact anomalous Hall conductivity of a SIESTA or Wannier90 model on a dense k mesh.

Usage:
    python exact_hall.py siesta SIESTA_DIR OUT [--label yn2] [--mesh 300] ...
    python exact_hall.py wannier W90_DIR OUT --fermi-energy EF [--seed y] [--mesh 300] ...

Options: --eta (interband regulator, eV), --broadening (Gaussian standard
deviation of the occupations, eV), --workers, --chunk, --emin/--emax (energy
grid relative to E_F, eV).

    sigma_xy(E) = (e^2/h) 2 pi/(N_k A) sum_k sum_n f(E_nk - E) Omega_n(k),
    Omega_n = -2 sum_{m != n} Im(<n|VX|m><m|VY|n>) / ((E_n - E_m)^2 + eta^2),

with the covariant velocity of LSQuant in both cases: for SIESTA, the dipoles
of the radial orbitals (by quadrature); for Wannier90, the position matrix
of seed_r.dat (write_rmn = true). This is the sign convention of the LSQuant
Kubo-Bastin curves. Energies are relative to E_F: the E_F of the .TSHS for
SIESTA, --fermi-energy (the scf Fermi energy, eV) for Wannier90.

Writes OUT.dat (E - E_F, sigma_xy in e^2/h). The k mesh is split in chunks
over --workers processes to bound the memory.
"""
import argparse
import multiprocessing as mp
import time

import numpy as np
from scipy.special import erf

from lsquant.bases import BlochMetricBasis

STATE = {}


def _siesta_chunk(k):
    model, mesh = STATE["model"], STATE["args"].mesh
    model.k_mesh = lambda supercell: k
    basis = BlochMetricBasis(model, [mesh, mesh, 1], velocity_mode="covariant")
    E, C = basis.physical_eigenvalues, basis.physical_eigenvectors
    Cd = np.conj(np.swapaxes(C, 1, 2))
    return E, Cd @ basis.operator_blocks("VX") @ C, Cd @ basis.operator_blocks("VY") @ C


def _wannier_chunk(k):
    model = STATE["model"]
    E, C = np.linalg.eigh(model.Hk(k))
    Cd = np.conj(np.swapaxes(C, 1, 2))
    X = Cd @ model.velocity_blocks(k, "x", mode="covariant") @ C
    Y = Cd @ model.velocity_blocks(k, "y", mode="covariant") @ C
    return E - STATE["args"].fermi_energy, X, Y


def _part(start):
    a, grid = STATE["args"], STATE["grid"]
    k = STATE["kpoints"][start:start + a.chunk]
    E, X, Y = (_siesta_chunk if a.kind == "siesta" else _wannier_chunk)(k)
    dE = E[:, :, None] - E[:, None, :]
    # Degenerate pairs carry no interband term; eta regularises near-degeneracies.
    curvature = -2 * np.sum(np.where(np.abs(dE) > 1e-6,
                                     np.imag(X * np.swapaxes(Y, 1, 2)) / (dE**2 + a.eta**2), 0), axis=2)
    occupation = 0.5 * (1 - erf((E[None] - grid[:, None, None]) / (np.sqrt(2) * a.broadening)))
    return np.einsum("gkn,kn->g", occupation, curvature)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("kind", choices=["siesta", "wannier"])
    ap.add_argument("dir")
    ap.add_argument("out")
    ap.add_argument("--label", default="yn2")
    ap.add_argument("--seed", default="y")
    ap.add_argument("--fermi-energy", type=float, default=None)
    ap.add_argument("--mesh", type=int, default=300)
    ap.add_argument("--eta", type=float, default=0.001)
    ap.add_argument("--broadening", type=float, default=0.002)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--chunk", type=int, default=500)
    ap.add_argument("--emin", type=float, default=-0.3)
    ap.add_argument("--emax", type=float, default=0.3)
    a = ap.parse_args()
    if a.kind == "siesta":
        from lsquant.models.siesta import load_siesta_model
        model = load_siesta_model(a.dir, a.label, dipoles=True)
        uc = model.uc
        kpoints = model.k_mesh([a.mesh, a.mesh, 1])
    else:
        if a.fermi_energy is None:
            ap.error("wannier needs --fermi-energy")
        from lsquant.models.orthogonal import write_wannier_geometry
        from lsquant.models.wannier import load_wannier_model
        write_wannier_geometry(a.dir, a.seed)  # seed.uc and seed.xyz from seed.win and the centres
        model = load_wannier_model(a.dir, a.seed)
        uc = np.asarray(model.uc, float)
        recip = 2 * np.pi * np.linalg.inv(uc).T
        g = np.arange(a.mesh) / a.mesh
        kpoints = np.array([[i, j, 0.0] for j in g for i in g]) @ recip
    area = abs(np.linalg.det(uc[:2, :2]))
    grid = np.arange(a.emin, a.emax + 1e-9, 0.001)
    STATE.update(model=model, args=a, grid=grid, kpoints=kpoints)
    total, t0 = np.zeros(len(grid)), time.time()
    starts = list(range(0, len(kpoints), a.chunk))
    with mp.get_context("fork").Pool(a.workers) as pool:
        for n, part in enumerate(pool.imap_unordered(_part, starts), 1):
            total += part
            if n % max(1, len(starts) // 10) == 0:
                print(f"{n}/{len(starts)} chunks, {time.time() - t0:.0f} s", flush=True)
    sigma = 2 * np.pi * total / (len(kpoints) * area)
    np.savetxt(a.out + ".dat", np.column_stack([grid, sigma]), fmt="%.4f %.6e",
               header=f"E-E_F (eV)  sigma_xy (e^2/h); {a.kind}; covariant velocity; "
                      f"{a.mesh}x{a.mesh} mesh; eta {a.eta} eV; Gaussian {a.broadening} eV")
    print("sigma_xy(E_F) = %.3f e^2/h" % np.interp(0.0, grid, sigma))


if __name__ == "__main__":
    main()
