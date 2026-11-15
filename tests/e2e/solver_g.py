#!/usr/bin/env python3
"""End-to-end test of the map G and SOLVER = FIXED_POINT (S-2, docs/SOLVER.md).

usage: solver_g.py --rhf1d PATH --gcheck1d PATH [--gcheck2d PATH]
                   [--src DIR] [--work DIR]

Ca II with angle-averaged PRD (two PRD lines) on FAL-C, H_6 passive:

 1. history: SOLVER = FIXED_POINT with NG_ORDER = 0 reproduces SOLVER = NG
    with NG_ORDER = PRD_NG_ORDER = PRD_ITER_LIMIT = 0 and PRD_DAMPING = 0.7
    bit for bit (PRD_DAMPING is unset in the FIXED_POINT runs: its default
    there is 0.7, under NG 1.0, N-38)
    (populations, J, spectrum) after 5 iterations, and so does
    FIXED_POINT without PRD_ITER_LIMIT in the file (its default is 0
    under a G solver). With PRD_N_MAX_ITER = 0 (every line in CRD, B-46)
    NG and FIXED_POINT agree as well, and differ from the runs with PRD.
 2. gcheck (rhf1d, and rhsc2d at 2 threads with --gcheck2d): G(x_a) twice
    gives bit-identical results, and G(x_b) between two calls of G(x_a)
    changes nothing.
 3. validation: SOLVER = FIXED_POINT with PRD_ITER_LIMIT = 1e-2 (NG's
    default, set in the file) stops at input with a message that names PRD_ITER_LIMIT; so do
    PRD_NG_DELAY and PRD_NG_PERIOD (set at all), COLL_RELAX_FACTOR other
    than 1, and SOLVE_NE = ITERATION with BACKGROUND_DOUBLE = FALSE. With
    SOLVER = NG, a negative NG_ORDER, NG_DELAY, PRD_NG_ORDER or
    PRD_NG_DELAY stops in readInput (B-45; negative orders crashed the run).

Exit status 0 = pass.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker, same_file
from rhtest.cli import parse
from rhtest.keywords import H6_CAII
from rhtest.run import EXIT_FATAL, execute, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d", "gcheck1d"), optional=("gcheck2d",))
work = WorkDir(a.src, a.work)
ck = Checker()

G_KEYWORDS = [("SOLVER", "FIXED_POINT"), ("NG_ORDER", "0"),
              ("PRD_NG_ORDER", "0"), ("PRD_ITER_LIMIT", "0")]


def prepare_g(name, geom, keywords):
    rundir = {"1d": "rhf1d", "2d": "rhsc2d"}[geom]
    atmos = {"1d": "../../Atmos/FALC_82.atmos", "2d": "../../Atmos/FALC_03x82.atmos"}[geom]
    return prepare(work, name, rundir,
                   [("ATMOS_FILE", atmos), ("STOKES_MODE", "NO_STOKES"),
                    ("STOKES_INPUT", "none"), ("PRD_N_MAX_ITER", "2"),
                    ("PRD_ANGLE_DEP", "PRD_ANGLE_AVER"), ("N_THREADS", "1"),
                    ("PRINT_CPU", "FALSE"), ("STARTING_J", "NEW_J"),
                    ("KURUCZ_DATA", "none"), ("SOLVE_NE", "NONE"),
                    ("COLL_RELAX_NSTEP", "0")] + keywords, H6_CAII)


# --- 1. Iteration history of FIXED_POINT = SOLVER = NG without Ng
runs = {}
UNSET = [k for k in G_KEYWORDS if k[0] != "PRD_ITER_LIMIT"] + [("PRD_ITER_LIMIT", None)]
CRD = [("PRD_N_MAX_ITER", "0")]
NG_DAMPED = [("SOLVER", "NG")] + G_KEYWORDS[1:] + [("PRD_DAMPING", "0.7")]
for name, extra in (("hist_ng", NG_DAMPED),
                    ("hist_fp", G_KEYWORDS), ("hist_fp_unset", UNSET),
                    ("hist_ng_crd", NG_DAMPED + CRD),
                    ("hist_fp_crd", G_KEYWORDS + CRD)):
    d = prepare_g(name, "1d", extra + [("N_MAX_ITER", "5"), ("ITER_LIMIT", "1e-14")])
    execute(a.rhf1d, d).require()             # 1: not converged, expected
    runs[name] = d
for ref, other in (("hist_ng", "hist_fp"), ("hist_ng", "hist_fp_unset"),
                   ("hist_ng_crd", "hist_fp_crd")):
    for f in ("pops.CA.out", "J.dat", "spectrum.out"):
        same = same_file(runs[ref], runs[other], f)
        ck.report(same, f"history  {other:13s} {f:13s} {'bit-identical' if same else 'DIFFERS'}",
                  f"history {other} {f}")
# The CRD runs must really be CRD: their populations differ from the PRD runs'
crd = not same_file(runs["hist_fp_crd"], runs["hist_fp"], "pops.CA.out")
ck.report(crd, f"history  PRD_N_MAX_ITER = 0 {'differs from PRD (CRD)' if crd else 'SAME AS PRD'}",
          "history: PRD_N_MAX_ITER = 0 did not change the run")

# --- 2. gcheck: twice and history independence
checks = [("gcheck_1d", "1d", a.gcheck1d, "1")]
if a.gcheck2d:
    checks.append(("gcheck_2d", "2d", a.gcheck2d, "2"))
for name, geom, program, threads in checks:
    d = prepare_g(name, geom, G_KEYWORDS + [("N_THREADS", threads)])
    r = execute(program, d, ["--steps", "3", "--lag", "1"])
    lines = [l for l in r.stdout.splitlines() if l.startswith("gcheck: twice")
             or l.startswith("gcheck: history")]
    r.require(ok=(0,), tail=30)
    if len(lines) != 2:
        sys.exit(f"{os.path.basename(program)} failed in {d}: expected 2 gcheck lines, "
                 f"end of out.log:\n{r.tail(30)}")
    for l in lines:
        ck.report("bit-identical" in l,
                  f"{name} ({threads} thread{'s' if threads != '1' else ''}): {l[8:]}",
                  f"{name}: {l}")

# --- 3. Conflicting keywords stop the run at input
d = prepare_g("validate", "1d", [("SOLVER", "FIXED_POINT"), ("PRD_ITER_LIMIT", "1e-2")])
r = execute(a.rhf1d, d)
ok = r.returncode == EXIT_FATAL and "PRD_ITER_LIMIT" in r.stdout and \
    "reading input file" not in r.stdout.split("TERMINATING")[0].split("readInput")[-1]
ck.report(ok, f"validate SOLVER = FIXED_POINT with PRD_ITER_LIMIT = 1e-2: "
              f"{'stops at input' if ok else 'NOT STOPPED'} (status {r.returncode})",
          "validation")

# Keywords without an "off" value stop the run when set at all
for kw, val, want in (("PRD_NG_DELAY", "2", "PRD_NG_DELAY is set"),
                      ("PRD_NG_PERIOD", "3", "PRD_NG_PERIOD is set"),
                      ("COLL_RELAX_FACTOR", "0.01", "COLL_RELAX_FACTOR = 0.01")):
    d = prepare_g(f"validate_{kw.lower()}", "1d", G_KEYWORDS + [(kw, val)])
    r = execute(a.rhf1d, d)
    ok = r.returncode == EXIT_FATAL and want in r.stdout
    ck.report(ok, f"validate SOLVER = FIXED_POINT with {kw} = {val}: "
                  f"{'stops at input' if ok else 'NOT STOPPED'} (status {r.returncode})",
              f"validation ({kw})")

# ne in the state requires the background in double (SOLVER.md Q11)
d = prepare_g("validate_ne", "1d", G_KEYWORDS + [("SOLVE_NE", "ITERATION"),
                                               ("BACKGROUND_DOUBLE", "FALSE")])
r = execute(a.rhf1d, d)
ok = r.returncode == EXIT_FATAL and "BACKGROUND_DOUBLE" in r.stdout and \
    "reading input file" not in r.stdout.split("TERMINATING")[0].split("readInput")[-1]
ck.report(ok, f"validate SOLVE_NE = ITERATION with BACKGROUND_DOUBLE = FALSE: "
              f"{'stops at input' if ok else 'NOT STOPPED'} (status {r.returncode})",
          "validation (BACKGROUND_DOUBLE)")

# Negative Ng orders and delays stop at input under every solver (B-45)
for kw, val, want in (("NG_ORDER", "-2", "Values of NG_ORDER and NG_DELAY"),
                      ("NG_DELAY", "-1", "Values of NG_ORDER and NG_DELAY"),
                      ("PRD_NG_ORDER", "-2", "Values of PRD_NG_ORDER and PRD_NG_DELAY"),
                      ("PRD_NG_DELAY", "-1", "Values of PRD_NG_ORDER and PRD_NG_DELAY")):
    d = prepare_g(f"validate_{kw.lower()}", "1d", [("SOLVER", "NG"), (kw, val)])
    r = execute(a.rhf1d, d)
    ok = r.returncode == EXIT_FATAL and want in r.stdout and \
        "readInput" in r.stdout.split(want)[0].rsplit("TERMINATING", 1)[-1]
    ck.report(ok, f"validate SOLVER = NG with {kw} = {val}: "
                  f"{'stops at input' if ok else 'NOT STOPPED'} (status {r.returncode})",
              f"validation ({kw} < 0)")

ck.finish("solver_g", work)
