"""Unit conversion constants shared by the model readers.

Every reader converts its electronic-structure code's units to the model
units (eV, Angstrom). Bohr and Hartree are the CODATA 2018 values; the
Rydberg values are NOT: ABACUS (through pyatb) and ThreeBodyTB (Atomdata.jl)
converted with their own constants, and the readers must undo each
conversion with the same number the code used. Both values therefore stay,
under distinct names, even though they differ from each other and from
``HARTREE_EV / 2``. Likewise ``ANGSTROM_BOHR`` is pyatb's Angstrom-per-Bohr,
kept as that code wrote it rather than as ``1 / BOHR_ANGSTROM``.

Pitfall: do not "fix" these numbers to CODATA or unify them: the readers
multiply by exactly the constant the producing code divided by, and a
1-digit change shifts every energy of that model.
"""

#: Bohr radius in Angstrom (CODATA 2018).
BOHR_ANGSTROM = 0.529177210903

#: Hartree in eV (CODATA 2018).
HARTREE_EV = 27.211386245988

#: Rydberg in eV as ABACUS/pyatb's ``constants.py`` defines it: ABACUS
#: divides its eV by this number, and the reader multiplies back by it.
RYDBERG_EV_ABACUS = 13.605698066

#: Rydberg in eV as ThreeBodyTB's Atomdata.jl defines it: ThreeBodyTB
#: converted with its own constant, kept distinct from the ABACUS value.
RYDBERG_EV_THREEBODYTB = 13.6056980659

#: Angstrom per Bohr as pyatb's ``constants.py`` defines it (ABACUS
#: ``LATTICE_CONSTANT`` conversion); kept as that code's number, not as
#: ``1 / BOHR_ANGSTROM``.
ANGSTROM_BOHR = 1.8897259886
