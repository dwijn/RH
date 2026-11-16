#!/usr/bin/env python3
"""End-to-end test of SOLVER = NEWTON_KRYLOV (S-4, docs/SOLVER.md).

usage: solver_newton.py --rhf1d PATH [--src DIR] [--work DIR]

Ca II with angle-averaged PRD (two PRD lines) on FAL-C, H_6 passive, as
solver_anderson.py:
 1. With the defaults (LOG variables) and with RELATIVE variables, 30 G
    evaluations take Newton steps with GMRES products and an accepted
    line-search point, N_MAX_ITER caps the G evaluations, and two runs
    give bit-identical output.
 2. Conflicting keywords stop the run at input: a NEWTON_* keyword with
    SOLVER = ANDERSON, an ANDERSON_* keyword and NG_ORDER > 0 with
    SOLVER = NEWTON_KRYLOV, a value out of range, and a map keyword
    (PRD_ITER_LIMIT) with a message that names NEWTON_KRYLOV.
"""

import os
import re
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker, same_file
from rhtest.cli import parse
from rhtest.keywords import H6_CAII
from rhtest.run import EXIT_FATAL, execute, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d",))
work = WorkDir(a.src, a.work)
ck = Checker()

# NG_DELAY and NG_PERIOD, which SOLVER = NEWTON_KRYLOV refuses when set, are
# removed in case the keyword file sets them (the shipped one did before S-10)
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
    return d, execute(a.rhf1d, d)


SHORT = [("ITER_LIMIT", "1e-14")]
NEVAL = 30

# --- 1. Newton steps, the evaluation cap, reproducibility
for label, extra in (("defaults", []), ("RELATIVE", [("NEWTON_VARIABLES", "RELATIVE")])):
    pair = []
    for k in (1, 2):
        d, r = run(f"nk_{len(extra)}_{k}", MAP + [("SOLVER", "NEWTON_KRYLOV")] + extra +
                   SHORT + [("N_MAX_ITER", str(NEVAL))])
        r.require()
        pair.append((d, r.stdout))
    log = pair[0][1]
    m = re.search(r"-- G evaluations: (\d+)\n -- Newton: steps (\d+), G evaluations base (\d+),"
                  r" Krylov (\d+), line search (\d+)", log)
    n_eval, n_steps, n_base, n_kry, n_line = map(int, m.groups()) if m else (0,) * 5
    n_acc = len(re.findall(r"\|d\|_max [^;]+; lambda ", log))
    ok = (m is not None and n_eval <= NEVAL and n_eval == n_base + n_kry + n_line and
          n_steps > 0 and n_kry > 0 and n_acc > 0 and
          all(same_file(pair[0][0], pair[1][0], f)
              for f in ("pops.CA.out", "J.dat", "spectrum.out")))
    ck.report(ok, f"{label}: {n_eval} G evaluations (base {n_base}, Krylov {n_kry}, line search"
                  f" {n_line}), {n_steps} Newton steps, {n_acc} accepted; two runs "
                  f"{'bit-identical' if ok else 'DIFFER (or no Newton step, or cap exceeded)'}",
              f"newton ({label})")

# --- 2. Conflicting keywords stop the run at input
for label, kws, want in (
        ("NEWTON_NOISE with SOLVER = ANDERSON",
         MAP + [("SOLVER", "ANDERSON"), ("NEWTON_NOISE", "1e-8")], "NEWTON_NOISE"),
        ("ANDERSON_DEPTH with SOLVER = NEWTON_KRYLOV",
         MAP + [("SOLVER", "NEWTON_KRYLOV"), ("ANDERSON_DEPTH", "3")], "ANDERSON_DEPTH"),
        ("NG_ORDER = 2 with SOLVER = NEWTON_KRYLOV",
         MAP + [("SOLVER", "NEWTON_KRYLOV"), ("NG_ORDER", "2")], "NG_ORDER = 2"),
        ("NG_PERIOD with SOLVER = NEWTON_KRYLOV",
         MAP + [("SOLVER", "NEWTON_KRYLOV"), ("NG_PERIOD", "3")], "NG_PERIOD is set"),
        ("NEWTON_KRYLOV_MAX = 0",
         MAP + [("SOLVER", "NEWTON_KRYLOV"), ("NEWTON_KRYLOV_MAX", "0")], "out of range"),
        ("PRD_ITER_LIMIT = 1e-2 with SOLVER = NEWTON_KRYLOV",
         [("SOLVER", "NEWTON_KRYLOV"), ("NG_ORDER", "0"), ("PRD_NG_ORDER", "0"),
          ("PRD_ITER_LIMIT", "1e-2")], "SOLVER = NEWTON_KRYLOV: PRD_ITER_LIMIT")):
    d, r = run("validate", kws)
    ok = r.returncode == EXIT_FATAL and want in r.stdout
    ck.report(ok, f"validate {label}: {'stops at input' if ok else 'NOT STOPPED'}"
                  f" (status {r.returncode})", f"validate {label}")

ck.finish("solver_newton", work)
