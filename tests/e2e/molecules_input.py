#!/usr/bin/env python3
"""End-to-end test of the H2 requirement of molecules.input (B-62).

usage: molecules_input.py --rhf1d PATH [--src DIR] [--work DIR] [-v]

RH uses atmos.H2, the molecule with ID H2, in the background opacities
(H2- free-free, Rayleigh scattering by H2), the H2 collisions and
hydrostatic equilibrium, and takes molecules[0] for H2 where Background()
frees the passive molecules and (until B-62) where Hydrostatic() decides
whether H2 is present at a depth point. readMolecularModels (readmolecule.c) required the first
molecule to be H2 with a substring test, so H2+ or H2O first passed, and an
empty list was not checked:

 1. H2+ first, no H2 in the list: atmos.H2 was NULL, segmentation fault in
    the background opacities (H2minus_ff).
 2. H2O first, H2 later: Background() freed the passive H2 that atmos.H2
    points to (read afterwards by the collisions of an active CO or H2);
    with HYDROSTATIC = TRUE, H2 was switched on and off by H2O's
    temperature range (converged nHtot of FAL-C off by 2.7e-4).
 3. No molecules: as 1.
 4. H2 first with its ID written "h2" in the file: the first-molecule
    check (on getMoleculeID's upper case) passed, atmos.H2 stayed NULL.

Each now stops at input with a message and status 70. The shipped lists
(H2 first) run in every other test.

Python 3 standard library only. Exit status 0 = pass.
"""
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.run import EXIT_FATAL, execute, prepare, template
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d",), verbose=True)
work = WorkDir(a.src, a.work)
ck = Checker(a.verbose)

KWS = [("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"), ("N_THREADS", "1"),
       ("PRINT_CPU", "FALSE"), ("KURUCZ_DATA", "none"), ("N_MAX_ITER", "1")]
ATOMS = [("../../Atoms/H_6.atom", "PASSIVE", "LTE_POPULATIONS")]


def shipped_list():
    """The molecule lines of the shipped examples/rhf1d/molecules.input."""
    with open(os.path.join(template(a.src, "rhf1d"), "molecules.input")) as f:
        return [l for l in f.read().splitlines() if ".molecule" in l and not l.lstrip().startswith("#")]


def case(name, lines, text, files=None):
    d = prepare(work, name, "rhf1d", KWS, ATOMS, copy=("ray.input",))
    for fname, content in (files or {}).items():
        with open(os.path.join(d, fname), "w") as f:
            f.write(content)
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

h2plus = [l for l in shipped if "/H2+.molecule" in l]
no_h2 = h2plus + [l for l in shipped if "/H2.molecule" not in l and "/H2+.molecule" not in l]
case("h2plus_first", no_h2, "The first molecule in molecules.input must be H2, not H2+")

h2o = [l for l in shipped if "/H2O.molecule" in l]
h2o_first = h2o + [l for l in shipped if "/H2O.molecule" not in l]
case("h2o_first", h2o_first, "The first molecule in molecules.input must be H2, not H2O")

# H2 first, but its file's ID in lower case: the check on the ID line
# (upper-cased by getMoleculeID) passes, atmos.H2 (ID compared as written)
# is not set: caught by the check after the list is read.
with open(os.path.join(a.src, "Molecules", "H2.molecule")) as f:
    h2_lower = f.read().replace("\n  H2\n", "\n  h2\n", 1)
lower = ["  h2.molecule   PASSIVE    LTE_POPULATIONS"] + shipped[1:]
case("h2_lowercase", lower, "the first molecule must be the only model with ID H2 (found none)",
     {"h2.molecule": h2_lower})

case("none", [], "Number of molecules in molecules.input is 0: the list must contain at least H2")

ck.finish("molecules_input", work)
