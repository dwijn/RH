#!/usr/bin/env python3
"""End-to-end test of ANGLE_SET = SET_EDDINGTON, which is 1-D only (Q-41).

usage: eddington.py --rhf1d P --rhsc2d P [--rhsc3d P] [--solveray2d P]
                    [--src DIR] [--work DIR] [-v]

SET_EDDINGTON is one ray at mu = 1/sqrt(3) with weight 1. The 2-D solver
of B-57, which solved each column on its own without horizontal transfer,
was retired in Q-41: a column of it is a 1-D run.

 1. rhf1d with SET_EDDINGTON converges (Ca II with PRD on FAL-C); NRAYS
    other than 1 is set to 1 with a warning.
 2. rhsc2d, rhsc3d and the 2-D solveray stop at readInput with status 70
    and a message that says SET_EDDINGTON is 1-D only.

Python 3 standard library only. Exit status 0 = pass.
"""
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.keywords import H6_CAII
from rhtest.run import EXIT_FATAL, execute, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d", "rhsc2d"), optional=("rhsc3d", "solveray2d"),
          verbose=True, src_help="RH source tree (Atoms/, Atmos/, Molecules/, examples/)")
work = WorkDir(a.src, a.work)
ck = Checker(a.verbose)

BASE = [("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"), ("N_THREADS", "1"),
        ("PRINT_CPU", "FALSE"), ("STARTING_J", "NEW_J"), ("KURUCZ_DATA", "none"),
        ("SOLVE_NE", "NONE"), ("PRD_N_MAX_ITER", "2"), ("PRD_ANGLE_DEP", "PRD_ANGLE_AVER"),
        ("LAMBDA_REF", "500.0")]


def run(name, program, geom, kws):
    d = prepare(work, name, geom, BASE + kws, H6_CAII, copy=(),
                required=("molecules.input", "ray.input"))
    return execute(program, d, timeout=300)


# 1: 1-D runs, NRAYS = 3 is set to 1 with a warning
r = run("e1d", a.rhf1d, "rhf1d", [("ANGLE_SET", "SET_EDDINGTON"), ("NRAYS", "3"),
                                  ("ITER_LIMIT", "1.0E-4"), ("N_MAX_ITER", "300"),
                                  ("ATMOS_FILE", "../../Atmos/FALC_82.atmos")])
ck.check(r.returncode == 0 and "Converged in" in r.stdout, "1-D SET_EDDINGTON run converges",
         f"(status {r.returncode})")
ck.check("Ignoring value of NRAYS in 1-D plane geometry when SET_EDDINGTON"
         in " ".join(r.stdout.split()), "1-D SET_EDDINGTON: NRAYS = 3 set to 1 with a warning")

# 2: refused in 2-D and 3-D, at readInput
MSG = "ANGLE_SET = SET_EDDINGTON is available in 1-D plane geometry only"
refusals = [("rhsc2d", a.rhsc2d, "rhsc2d")]
if a.rhsc3d:
    refusals.append(("rhsc3d", a.rhsc3d, "rhsc3d"))
if a.solveray2d:
    refusals.append(("2-D solveray", a.solveray2d, "rhsc2d"))
for what, prog, geom in refusals:
    r = run("refuse_" + what.replace(" ", "_"), prog, geom,
            [("ANGLE_SET", "SET_EDDINGTON"), ("N_MAX_ITER", "2")])
    flat = " ".join(r.stdout.split())
    ck.check(r.returncode == EXIT_FATAL and MSG in flat and "readInput" in flat,
             f"SET_EDDINGTON refused at input: {what}", f"(status {r.returncode})")

ck.finish("eddington", work)
