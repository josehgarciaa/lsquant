# LSQUANT: A tool for linear quantum transport

LSQUANT is a tool for quantum transport calculations in tight-binding models. It provides interfaces with other tight-binding modelling packages such as **Wannier90**, **KWANT**, and **TBModels**. Its core is written in C++ and uses OpenMP for optimized calculations.

LSQUANT allows for the following kind of calculations:

1. Nonequilibrium electrical response of different operators within the linear regime. Some particular examples are:
  1. The conductivity.
  2. The spin Hall conductivity.
  3. The spin susceptibility.
2. Time-correlation which are the basis for computing:
  1. The Mean-square spreading.
  2. \subpage spinrelaxation.

## Getting started

Build LSQUANT with CMake (`cmake -B build && cmake --build build -j`); the
installer helper is `install.py` at the repository root, and the full
prerequisites and build options are described in the top-level `README.md`. This
page does not repeat those steps so they cannot drift out of sync.

Worked, runnable calculations live in the `examples/` directory, one
self-contained subdirectory per topic (1D density of states, localization,
Haldane Hall, spin precession, charge transport, and spin relaxation). Each
example carries its own `README.md` explaining the physics it demonstrates and
how to run it.


## Core engine

The new C++ core computes every KPM quantity by running a small instruction list (a
*plan*) over trace vectors; see \subpage core_engine.

## Python layer

The `lsquant` Python package drives the core from Python:

- \subpage python_layer: the pybind11 module, the plans and the density of states;
- \subpage kubo: Kubo-Greenwood and Kubo-Bastin conductivities from rank-2 moments;
- \subpage time_evolution: mean-square displacement, projected operators, time correlations;
- \subpage bases: real-space, Bloch and eigen bases built from one Wannier model;
- \subpage nonorthogonal: an overlap matrix S, the Löwdin and generalised routes, the
  overlap appliers.

## Run files and the command line

Calculations are JSON run files written with observable words; the command line tool and
the MCP server are covered in \subpage run_files.

## Testing

Every numerical result is checked against committed golden masters; see
\subpage testing for the gate, how to read its output and the unit tests.
