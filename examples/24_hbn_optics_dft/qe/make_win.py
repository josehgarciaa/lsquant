#!/usr/bin/env python3
"""Write the Wannier90 input hbn.win from the k list that open_grid.x prints.

Usage (in qe/pbe or qe/hse, after open_grid.x):
    python ../make_win.py open_grid.out FROZEN_MAX_eV

Eight projected Wannier functions, s and p on B and on N, disentangled from
the 40 bands of the scf run. The lowest conduction states at Gamma are
nearly-free-electron states of the vacuum, which s and p orbitals on the
atoms cannot represent, and they lie between the pi* and sigma* bands, so no
energy window separates them. The disentanglement therefore selects by
projectability (Wannier90 4, dis_froz_proj; Qiao, Pizzi and Marzari, npj
Comput. Mater. 9, 208 (2023)): states with projectability above 0.95 on the
pseudo-atomic orbitals are frozen, below 0.01 excluded. The energy frozen
window keeps the valence bands (FROZEN_MAX, absolute QE energy in eV, in
the gap). pw2wannier90 must write the projections on the pseudo-atomic
orbitals (atom_proj = .true. in pw2wan.in), which auto_projections requests
in hbn.nnkp; with two orbitals (s, p) per species in the pseudopotentials they
are the eight Wannier functions. No maximal localisation follows
(num_iter = 0): the functions are the projected ones of the disentangled
subspace. After 500 minimisation steps H(R) was still 0.1-0.2 eV at
11-17 Angstrom; the projected functions give at most 6 meV beyond 11 Angstrom.
The spreads Wannier90 prints stay near 95 Angstrom^2 per function in both
cases, presumably from the small vacuum weight of the subspace; only H(R)
and the in-plane position matrix enter LSQuant. write_rmn writes hbn_r.dat
(position matrix) for the covariant velocity of LSQuant.
"""
import sys

A = 2.504
C = 20.0


def main():
    out, froz = sys.argv[1], float(sys.argv[2])
    lines = open(out).read().splitlines()
    start = next(i for i, l in enumerate(lines) if "List to be put in the .win file" in l) + 1
    kpts = []
    for l in lines[start:]:
        parts = l.split()
        if len(parts) >= 3:
            try:
                kpts.append([float(x) for x in parts[:3]])
                continue
            except ValueError:
                pass
        if kpts:
            break
    assert len(kpts) == 144, len(kpts)
    win = f"""num_bands       = 40
num_wann        = 8
dis_win_min     = -40.0
dis_win_max     =  40.0
dis_froz_min    = -40.0
dis_froz_max    = {froz:.4f}
dis_froz_proj   = true
dis_proj_min    = 0.01
dis_proj_max    = 0.95
dis_num_iter    = 2000
num_iter        = 0
write_hr        = true
write_rmn       = true
write_xyz       = true
translation_centre_frac = 0.0 0.0 0.5
use_ws_distance = true

begin unit_cell_cart
ang
  {A:.6f}  0.000000  0.000000
  {-A/2:.6f}  {A*3**0.5/2:.6f}  0.000000
  0.000000  0.000000  {C:.6f}
end unit_cell_cart

begin atoms_frac
B  0.333333333333  0.666666666667  0.5
N  0.666666666667  0.333333333333  0.5
end atoms_frac

auto_projections = true

mp_grid = 12 12 1
begin kpoints
""" + "\n".join(f"  {k[0]:.10f} {k[1]:.10f} {k[2]:.10f}" for k in kpts) + "\nend kpoints\n"
    open("hbn.win", "w").write(win)
    print(f"hbn.win: {len(kpts)} k points, frozen window up to {froz} eV")


if __name__ == "__main__":
    main()
