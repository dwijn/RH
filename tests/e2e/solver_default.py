#!/usr/bin/env python3
"""End-to-end test of the default solver (S-10, docs/SOLVER.md, Anderson as
the default solver).

usage: solver_default.py --rhf1d PATH --rhsc2d PATH --rhsc3d PATH
                         [--src DIR] [--work DIR]

Each program runs its shipped example (run/) for one or two iterations:
 1. rhf1d, rhsc2d and rhsc3d iterate with ANDERSON when the keyword
    file does not set SOLVER, and say so (" -- SOLVER = ANDERSON (the
    default)"). rhf1d extrapolates the whole state
    (ANDERSON_EXTRAPOLATE = ALL), rhsc2d and rhsc3d only the populations
    and rho (POPULATIONS_RHO, S-11); the banner marks it "(the default)".
    Without PRD lines (PRD_N_MAX_ITER = 0) rhsc2d's default runs as
    POPULATIONS without the warning that a keyword file setting
    POPULATIONS_RHO gets.
 2. A 1-D keyword file written for the old default (Ng keywords, no SOLVER)
    stops at input (status 70) with messages that name the default and
    SOLVER = NG as a fix; with SOLVER = NG added it runs, and the banner
    does not call NG the default.
Exit status 0 = pass.
"""
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.keywords import edit_file
from rhtest.run import EXIT_FATAL, copy_dir, execute, template
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d", "rhsc2d", "rhsc3d"))
work = WorkDir(a.src, a.work)
ck = Checker()


def run(name, program, rundir, keywords):
    """The shipped example of rundir, with keywords."""
    d = copy_dir(work, name, template(a.src, rundir))
    edit_file(os.path.join(d, "keyword.input"),
              [("N_THREADS", "1"), ("PRINT_CPU", "FALSE")] + keywords)
    return execute(program, d)


def report_run(label, ok, r):
    ck.report(ok, f"{label}: {'ok' if ok else 'FAILED'} (status {r.returncode})", label)


# --- 1. The default of each geometry
for name, program, rundir, want in (
        ("rhf1d", a.rhf1d, "rhf1d", "ALL"),
        ("rhsc2d", a.rhsc2d, "rhsc2d", "POPULATIONS_RHO"),
        ("rhsc3d", a.rhsc3d, "rhsc3d", "POPULATIONS_RHO")):
    r = run(f"default_{name}", program, rundir, [("N_MAX_ITER", "2")])
    ok = r.ran and " -- SOLVER = ANDERSON (the default)\n" in r.stdout \
        and "-- Anderson: extrapolated" in r.stdout \
        and f"ANDERSON_EXTRAPOLATE = {want} (the default):" in r.stdout
    report_run(f"{name} without SOLVER iterates with ANDERSON, {want} (the default)", ok, r)

# --- Without PRD lines: the 2-D default is silent, an explicit setting warns
WARN = "ANDERSON_EXTRAPOLATE = POPULATIONS_RHO, but there are no PRD"
for label, extra, warn in (("default", [], False),
                           ("explicit", [("ANDERSON_EXTRAPOLATE", "POPULATIONS_RHO")], True)):
    r = run(f"noprd_{label}", a.rhsc2d, "rhsc2d", [("N_MAX_ITER", "2"), ("PRD_N_MAX_ITER", "0")] + extra)
    ok = r.ran and (WARN in r.stdout) == warn \
        and "ANDERSON_EXTRAPOLATE = POPULATIONS_RHO" in r.stdout
    report_run(f"rhsc2d without PRD lines, POPULATIONS_RHO {label}: "
               f"{'warns' if warn else 'no warning'}", ok, r)

# --- 2. A keyword file written for the old default
OLD = [("NG_ORDER", "2"), ("NG_DELAY", "10"), ("NG_PERIOD", "3"),
       ("PRD_ITER_LIMIT", "1.0E-2")]
r = run("old_file", a.rhf1d, "rhf1d", OLD)
ok = r.returncode == EXIT_FATAL and "SOLVER = ANDERSON (the default)" in r.stdout and \
    "Set SOLVER = NG" in r.stdout and "4 keyword settings conflict" in r.stdout
report_run("rhf1d, old keyword file without SOLVER stops at input", ok, r)
r = run("old_file_ng", a.rhf1d, "rhf1d", OLD + [("SOLVER", "NG"), ("N_MAX_ITER", "2")])
ok = r.ran and " -- SOLVER = NG\n" in r.stdout
report_run("rhf1d, old keyword file with SOLVER = NG runs with NG", ok, r)

ck.finish("solver_default", work)
