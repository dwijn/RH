#!/usr/bin/env python3
"""End-to-end test of which molecules can be active (B-69).

usage: active_molecules.py --rhf1d PATH [--src DIR] [--work DIR] [-v]

Collisional rates exist only for CO (cocollisions.c, Ayres & Wiedemann
1989). Until B-69 the collisions were chosen by substring: an ID containing
"CO" got CO's rates, one containing "H2" (H2, H2O, H2+) got a copy of them
with only the names changed (h2collisions.c: CO's vibrational constant and
CO's Landau-Teller constants), anything else stopped at the first
iteration's setup. Now readMolecule stops at input for an active molecule
other than CO:

 1. H2 active (its own line list): stops, naming H2. Before, it ran with
    CO's rates.
 2. H2O active (no line list): stops with the same message, not the
    "No lines" one.
 3. CO active (the v <= 3, J <= 49 line list, one iteration): runs to its
    end. No other test makes a molecule active.

Python 3 standard library only. Exit status 0 = pass.
"""
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.run import EXIT_FATAL, RAN, execute, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d",), verbose=True)
work = WorkDir(a.src, a.work)
ck = Checker(a.verbose)

KWS = [("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"), ("N_THREADS", "1"),
       ("PRINT_CPU", "FALSE"), ("KURUCZ_DATA", "none"), ("N_MAX_ITER", "1")]
ATOMS = [("../../Atoms/H_6.atom", "PASSIVE", "LTE_POPULATIONS")]
CO_LINES = "../../Molecules/CO/vmax=3_Jmax=49_dv=1_26"


def shipped_list():
    """The molecule lines of the shipped rhf1d/run/molecules.input."""
    with open(os.path.join(a.src, "rhf1d", "run", "molecules.input")) as f:
        return [l for l in f.read().splitlines() if ".molecule" in l and not l.lstrip().startswith("#")]


def run(name, lines, files=None):
    d = prepare(work, name, "rhf1d", KWS, ATOMS, copy=("ray.input",))
    for fname, content in (files or {}).items():
        with open(os.path.join(d, fname), "w") as f:
            f.write(content)
    with open(os.path.join(d, "molecules.input"), "w") as f:
        f.write("  %d\n\n%s\n" % (len(lines), "\n".join(lines)))
    return execute(a.rhf1d, d, timeout=600)


def activate(lines, name, molfile=None):
    """The list with molecule file <name>.molecule ACTIVE (or replaced by
    the local file molfile, ACTIVE)."""
    out = []
    for l in lines:
        if "/%s.molecule" % name in l:
            l = "  %s   ACTIVE   LTE_POPULATIONS" % (molfile or l.split()[0])
        out.append(l)
    assert out != lines, name
    return out


def stops(name, lines, mol, files=None):
    r = run(name, lines, files)
    flat = " ".join(r.stdout.split())
    text = "Molecule %s cannot be active: collisional rates exist only for CO" % mol
    ck.check(r.returncode == EXIT_FATAL and text in flat, f"{name}: stops at input ({text!r})",
             f"(status {r.returncode})" if r.returncode != EXIT_FATAL else "")


shipped = shipped_list()
ck.check(shipped and "/H2.molecule" in shipped[0] and any("/H2O.molecule" in l for l in shipped)
         and any("/CO.molecule" in l for l in shipped),
         "shipped molecules.input lists H2 first, and H2O and CO")

# H2 with a line list (CO's, to get past the "No lines" check of the old code)
with open(os.path.join(a.src, "Molecules", "H2.molecule")) as f:
    h2 = f.read().rstrip("\n") + "\n  %s\n" % CO_LINES
stops("h2_active", activate(shipped, "H2", "h2_lines.molecule"), "H2", {"h2_lines.molecule": h2})

stops("h2o_active", activate(shipped, "H2O"), "H2O")

with open(os.path.join(a.src, "Molecules", "CO_NLTE.molecule")) as f:
    co = [l for l in f.read().splitlines() if "Molecules/CO/" not in l]
co.append("  " + CO_LINES)
r = run("co_active", activate(shipped, "CO", "co_small.molecule"), {"co_small.molecule": "\n".join(co) + "\n"})
ck.check(r.returncode in RAN and os.path.exists(os.path.join(r.dir, "pops_mol.CO.out")),
         "co_active: an active CO runs to its end and writes its populations",
         f"(status {r.returncode})\n{r.tail(15)}")

ck.finish("active_molecules", work)
