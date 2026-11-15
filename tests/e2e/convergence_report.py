#!/usr/bin/env python3
"""End-to-end test of the convergence report of the G solvers (N-3).

usage: convergence_report.py --rhf1d PATH [--src DIR] [--work DIR]

Ca II with angle-averaged PRD (two PRD lines) on FAL-C, H_6 passive, as
solver_anderson.py:
 1. The per-iteration lines carry the residual of the rho blocks (" rho
    CA: delta = ...") and, with SOLVE_NE = ITERATION, that of ne, with
    all three G solvers.
 2. SOLVER = FIXED_POINT ends with " -- Estimated contraction ..." after
    the convergence message, which keeps its format; with fewer than four
    iterations the line says "no estimate". ANDERSON, NEWTON_KRYLOV, NG and
    FIXED_POINT with NG_ORDER > 0 print no such line (the estimate understates their error by up to 20
    times, docs/changes/N-3.md).
 3. The FIXED_POINT error estimate (residual / (1 - contraction)) at
    ITER_LIMIT = 1e-4 is within a factor 3 of the actual difference to an
    ANDERSON run at ITER_LIMIT = 1e-7 (the reference). The numbers are
    printed.
That the reporting leaves the results alone is checked by run_tests.sh
--exact against a build without it.
"""

import os
import re
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.keywords import H6_CAII
from rhtest.output import read_pops
from rhtest.run import execute, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d",))
work = WorkDir(a.src, a.work)
ck = Checker()

MAP = [("NG_ORDER", "0"), ("NG_DELAY", None), ("NG_PERIOD", None),
       ("PRD_NG_ORDER", "0"), ("PRD_ITER_LIMIT", "0")]


def run(name, keywords):
    d = prepare(work, name, "rhf1d",
                [("ATMOS_FILE", "../../Atmos/FALC_82.atmos"),
                 ("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
                 ("PRD_N_MAX_ITER", "2"), ("PRD_ANGLE_DEP", "PRD_ANGLE_AVER"),
                 ("N_THREADS", "1"), ("PRINT_CPU", "FALSE"),
                 ("STARTING_J", "NEW_J"), ("KURUCZ_DATA", "none"),
                 ("SOLVE_NE", "NONE"), ("COLL_RELAX_NSTEP", "0")] + keywords, H6_CAII)
    r = execute(a.rhf1d, d, timeout=900).require()
    return d, r.stdout


EST = re.compile(r" -- Estimated contraction ([0-9.]+) per iteration \(last (\d+) "
                 r"iterations\); the error of the populations may be up to about "
                 r"([0-9.E+-]+) \(residual / \(1 - contraction\)\)\n")
SOLVERS = {"FIXED_POINT": [("SOLVER", "FIXED_POINT")],
           "ANDERSON": [("SOLVER", "ANDERSON")],
           "NEWTON_KRYLOV": [("SOLVER", "NEWTON_KRYLOV")]}

# --- 1, 2. rho and ne lines, the contraction line after the convergence message
for name, kws in SOLVERS.items():
    d, out = run(f"lines_{name}", MAP + kws + [("SOLVE_NE", "ITERATION"),
                                               ("ITER_LIMIT", "1e-3"),
                                               ("N_MAX_ITER", "60")])
    # iterations that evaluated the residual: iterateAnderson prints a header
    # for a rejected trial point without one
    n_it = sum(1 for blk in re.split(r"\n -- Iteration +\d+\n", out)[1:]
               if re.search(r"^ ?\S+\s*, delta = ", blk, re.M))
    n_rho = len(re.findall(r"^ rho CA: delta = [0-9.E+-]+$", out, re.M))
    n_ne = len(re.findall(r"^ ne: delta = [0-9.E+-]+$", out, re.M))
    conv = re.search(r"\n -- Converged (in \d+ iterations|at base point \d+): max change "
                     r"[0-9.E+-]+ < ITER_LIMIT = [0-9.E+-]+\n(.*\n)", out)
    est = EST.match(conv.group(2)) if conv else None
    want_est = name == "FIXED_POINT"
    ok = (n_it > 3 and n_rho == n_it and n_ne == n_it and conv is not None and
          (est is not None) == want_est and
          ("Estimated contraction" in out) == want_est)
    ck.report(ok, f"{name}: {n_it} iterations, {n_rho} rho and {n_ne} ne lines; "
                  f"{est.group(0).strip() if est else 'no estimate line'}",
              f"report lines ({name})")

d, out = run("ng", MAP + [("SOLVER", "NG"), ("ITER_LIMIT", "1e-3"), ("N_MAX_ITER", "60")])
ok = "Converged in" in out and "Estimated contraction" not in out
ck.report(ok, f"NG: {'no estimate line' if ok else 'ESTIMATE LINE PRESENT or no convergence'}",
          "NG prints no estimate")

d, out = run("fp_ng", MAP + [("SOLVER", "FIXED_POINT"), ("NG_ORDER", "2"),
                            ("ITER_LIMIT", "1e-2"), ("N_MAX_ITER", "40")])
ok = re.search(r" -- Converged (in|at)", out) and "Estimated contraction" not in out
ck.report(bool(ok), f"FIXED_POINT with NG_ORDER = 2: {'no estimate line' if ok else 'WRONG'}",
          "FIXED_POINT with Ng prints no estimate")

d, out = run("short", MAP + [("SOLVER", "FIXED_POINT"), ("ITER_LIMIT", "1e-14"),
                             ("N_MAX_ITER", "2")])
ok = "NOT CONVERGED after 2 iterations" in out and \
    "Estimated contraction: no estimate (fewer than 4" in out
ck.report(ok, f"two iterations: {'no estimate' if ok else 'WRONG'}",
          "no estimate with two iterations")


# --- 3. The estimate against the actual difference to a tight run
def err(d1, d2):
    n1, n2 = read_pops(os.path.join(d1, "pops.CA.out")), \
        read_pops(os.path.join(d2, "pops.CA.out"))
    if len(n1) != len(n2):
        raise ValueError(f"pops.CA.out of {d1} and {d2} differ in size")
    return max(abs((x - y) / y) for x, y in zip(n1, n2) if y != 0)


ref, out_ref = run("ref_ANDERSON_1e-7", MAP + SOLVERS["ANDERSON"] +
                   [("N_MAX_ITER", "400"), ("ITER_LIMIT", "1e-7")])
d4, out4 = run("FIXED_POINT_1e-4", MAP + SOLVERS["FIXED_POINT"] +
               [("N_MAX_ITER", "400"), ("ITER_LIMIT", "1e-4")])
est = EST.search(out4)
converged = re.search(r" -- Converged (in|at)", out_ref) and est is not None and \
    re.search(r" -- Converged (in|at)", out4)
ck.check(converged, "FIXED_POINT 1e-4 and reference run converged, estimate line present")
if converged:
    lam, e_est, e_act = float(est.group(1)), float(est.group(3)), err(d4, ref)
    r = e_est / e_act if e_act > 0 else float("inf")
    ck.report(1.0 / 3.0 < r < 3.0,
              f"FIXED_POINT: contraction {lam:.3f}, estimated error {e_est:.2E}, "
              f"actual |x(1e-4) - x(ref)| {e_act:.2E}, ratio {r:.2f}",
              f"FIXED_POINT estimate {e_est:.2E} vs actual {e_act:.2E}")
else:
    ck.check(False, "FIXED_POINT estimate within a factor 3 of the actual error",
             "(not compared)")

ck.finish("convergence_report", work)
