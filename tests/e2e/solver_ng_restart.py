#!/usr/bin/env python3
"""End-to-end test of the legacy chemical equilibrium (N-36).

usage: solver_ng_restart.py --rhf1d PATH [--src DIR] [--work DIR]

SOLVER = NG with SOLVE_NE = ITERATION on FAL-C, Ng acceleration on:

 1. Fresh and restarted run reach the same converged state. H_6_CRD active,
    ITER_LIMIT 1e-8. A fresh run converges (about 95 iterations); a restart
    from its populations and mean intensities (STARTING_J = OLD_J, which
    makes the populations OLD_POPULATIONS) converges again (about 70
    iterations; the first Background() solves ne from the atmosphere file,
    so the restart starts a few per cent away). The two populations must
    agree to 2e-6 (relative, every level and depth point). They do when the
    neutral fraction of an active atom in the chemical equilibrium comes from
    its populations once those are initialized, in the iteration of both runs
    (2.6e-7, the convergence level); with the LTE fraction in a fresh run,
    as before N-36, they differ by 6.1e-6. This tests that both runs
    converge to the same point, not that the restart stays put (it does not
    in its first iterations, on either side of the change).
 2. The setup Background() runs before initSolution(), while the populations
    of an active atom are still zero (or, for H, the populations of the
    atmosphere file): it must not turn them into a NaN neutral fraction.
    Two runs of N_MAX_ITER = 2 must start iterating and finish (status 0 or
    1, not converged is fine) without a NaN or a signal: C I+II active with
    H_6 passive, and H_6_CRD active with HYDROGEN_LTE = TRUE (both crashed
    with SIGFPE in the first version of N-36).

Exit status 0 = pass.
"""
import os
import re
import shutil
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.output import read_pops
from rhtest.run import EXIT_OK, RAN, execute, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d",))
work = WorkDir(a.src, a.work)
ck = Checker()

TOL = 2e-6                      # fresh vs restart: 6.1e-6 before N-36, 2.6e-7 after


def run(name, init, startj, from_dir=None, atoms=None, extra=(), nmax="300",
        converge=True):
    d = prepare(work, name, "rhf1d",
                [("ATMOS_FILE", "../../Atmos/FALC_82.atmos"),
                 ("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
                 ("N_THREADS", "1"), ("PRINT_CPU", "FALSE"),
                 ("KURUCZ_DATA", "none"), ("PRD_N_MAX_ITER", "0"),
                 ("STARTING_J", startj), ("SOLVER", "NG"),
                 ("SOLVE_NE", "ITERATION"), ("BACKGROUND_DOUBLE", "TRUE"),
                 ("N_MAX_SCATTER", "50"), ("ITER_LIMIT", "1e-8"),
                 ("N_MAX_ITER", nmax), ("NG_ORDER", "2"), ("NG_DELAY", "10"),
                 ("NG_PERIOD", "3"), *extra],
                atoms or [("../../Atoms/H_6_CRD.atom", "ACTIVE", init, "pops.H.out")],
                copy=(), required=("molecules.input", "ray.input"))
    if from_dir:
        for f in ("pops.H.out", "J.dat"):
            shutil.copy(os.path.join(from_dir, f), d)
    r = execute(a.rhf1d, d)
    n = len(re.findall(r"^ -- Iteration +\d+", r.stdout, re.M))
    r.require(ok=(EXIT_OK,) if converge else RAN)      # 1: not converged, fine if not converge
    return d, n, r.stdout


d1, n1, _ = run("fresh", "ZERO_RADIATION", "NEW_J")
d2, n2, _ = run("restart", "OLD_POPULATIONS", "OLD_J", from_dir=d1)
p1, p2 = read_pops(os.path.join(d1, "pops.H.out")), read_pops(os.path.join(d2, "pops.H.out"))
if len(p1) != len(p2) or not p1:
    sys.exit(f"population files differ in size: {len(p1)} and {len(p2)} values")
err = max(abs(x - y) / max(abs(x), 1e-300) for x, y in zip(p1, p2))
ck.report(err < TOL, f"fresh run {n1} iterations, restart {n2} iterations; populations differ by "
                     f"{err:.2e} (limit {TOL:.0e})",
          f"fresh and restarted runs differ by {err:.2e}")

# --- 2. setup before initSolution: no NaN, the run iterates
for name, atoms, extra in (
        ("setup_C", [("../../Atoms/H_6.atom", "PASSIVE", "LTE_POPULATIONS"),
                     ("../../Atoms/C_I+II_9.atom", "ACTIVE", "ZERO_RADIATION", "pops.C.out")],
         ()),
        ("setup_HLTE", None, (("HYDROGEN_LTE", "TRUE"),))):
    d, n, out = run(name, "ZERO_RADIATION", "NEW_J", atoms=atoms, extra=extra,
                    nmax="2", converge=False)
    bad = n < 2 or re.search(r"\bnan\b|FATAL|TERMINATING", out, re.I) is not None
    ck.report(not bad, f"setup    {name:10s} {n} iterations: {'FAILED' if bad else 'ok'}",
              f"{name}: {n} iterations, NaN or error in the output")

ck.finish("solver_ng_restart", work)
