#!/usr/bin/env python3
"""End-to-end test of the H2 requirement of molecules.input (B-62).

usage: molecules_input.py --rhf1d PATH [--src DIR] [--work DIR] [-v]

RH uses atmos.H2, the molecule with ID H2, in the background opacities
(H2- free-free, Rayleigh scattering by H2), the H2 collisions and
hydrostatic equilibrium, and takes molecules[0] for H2 where Background()
frees the passive molecules and (until B-62) where Hydrostatic() decides
whether H2 is present at a depth point. readMolecules required the first
molecule to be H2 with a substring test, so H2+ or H2O first passed, and an
empty list was not checked:

 1. H2+ first, no H2 in the list: atmos.H2 was NULL, segmentation fault in
    the background opacities (H2minus_ff).
 2. H2O first, H2 later: Background() freed the passive H2 that atmos.H2
    points to (read afterwards by the collisions of an active CO or H2);
    with HYDROSTATIC = TRUE, H2 was switched on and off by H2O's
    temperature range (converged nHtot of FAL-C off by 2.7e-4).
 3. No molecules: as 1.

Each now stops at input with a message and status 70. The shipped lists
(H2 first) run in every other test.

Python 3 standard library only. Exit status 0 = pass.
"""
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.run import EXIT_FATAL, execute, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d",), verbose=True)
work = WorkDir(a.src, a.work)
ck = Checker(a.verbose)

KWS = [("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"), ("N_THREADS", "1"),
       ("PRINT_CPU", "FALSE"), ("KURUCZ_DATA", "none"), ("N_MAX_ITER", "1")]
ATOMS = [("../../Atoms/H_6.atom", "PASSIVE", "LTE_POPULATIONS")]


def shipped_list():
    """The molecule lines of the shipped rhf1d/run/molecules.input."""
    with open(os.path.join(a.src, "rhf1d", "run", "molecules.input")) as f:
        return [l for l in f.read().splitlines() if ".molecule" in l and not l.lstrip().startswith("#")]


def case(name, lines, text):
    d = prepare(work, name, "rhf1d", KWS, ATOMS, copy=("ray.input",))
    with open(os.path.join(d, "molecules.input"), "w") as f:
        f.write("  %d\n\n%s\n" % (len(lines), "\n".join(lines)))
    r = execute(a.rhf1d, d, timeout=300)
    flat = " ".join(r.stdout.split())
    ck.check(r.returncode == EXIT_FATAL and text in flat, f"{name}: stops ({text!r})",
             f"(status {r.returncode})" if r.returncode != EXIT_FATAL else "")


shipped = shipped_list()
ck.check(shipped and "/H2.molecule" in shipped[0] and any("/H2+.molecule" in l for l in shipped)
         and any("/H2O.molecule" in l for l in shipped),
         "shipped molecules.input lists H2 first, and H2+ and H2O")

no_h2 = [l for l in shipped if "/H2.molecule" not in l]
case("h2plus_first", no_h2, "The first molecule in molecules.input must be H2, not H2+")

h2o = [l for l in shipped if "/H2O.molecule" in l]
h2o_first = h2o + [l for l in shipped if "/H2O.molecule" not in l]
case("h2o_first", h2o_first, "The first molecule in molecules.input must be H2, not H2O")

case("none", [], "molecules.input lists 0 molecules: it must list at least H2")

ck.finish("molecules_input", work)
