#!/usr/bin/env python3
"""End-to-end test of SOLVER = ANDERSON (S-3, docs/SOLVER.md).

usage: solver_anderson.py --rhf1d PATH [--src DIR] [--work DIR]

Ca II with angle-averaged PRD (two PRD lines) on FAL-C, H_6 passive, as
solver_g.py:
 1. ANDERSON_DEPTH = 0 repeats SOLVER = FIXED_POINT with NG_ORDER = 0 bit
    for bit (populations, J, spectrum) after 8 iterations: plain steps
    are copied, never transformed to the Anderson variables and back.
 2. With the defaults (depth 10, LOG variables) and with RELATIVE
    variables and PER_BLOCK coefficients, 15 iterations extrapolate, and
    two runs give bit-identical output.
 3. Conflicting keywords stop the run at input: an ANDERSON_* keyword with
    SOLVER = FIXED_POINT, NG_ORDER > 0 and NG_DELAY (set at all) with
    SOLVER = ANDERSON, and a map keyword (PRD_ITER_LIMIT) with a message
    that names ANDERSON; also an invalid ANDERSON_EXTRAPOLATE and
    ANDERSON_SAFEGUARD_WINDOW = 0.
 4. ANDERSON_EXTRAPOLATE (S-11): with depth 0, POPULATIONS and
    POPULATIONS_RHO repeat FIXED_POINT bit for bit as well; POPULATIONS
    extrapolates the Ca II populations only and POPULATIONS_RHO also the
    two rho blocks (start banner); a POPULATIONS_RHO run with PER_BLOCK
    coefficients (three solves) is bit-reproducible;
    and with ANDERSON_SAFEGUARD = 0.5 and ANDERSON_SAFEGUARD_WINDOW = 1,
    which force rejections, the carried
    blocks (rho, J) of every point equal those of G at the last accepted
    point bit for bit (the solver's own check, end of the run; it checks
    the packing and the copy of the carried blocks, not G). Without PRD
    lines (PRD_N_MAX_ITER = 0: all lines in CRD) POPULATIONS_RHO warns
    that it is the same as POPULATIONS, and runs.
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

# NG_DELAY and NG_PERIOD, which SOLVER = ANDERSON refuses when set, are
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

# --- 1. Depth 0 is FIXED_POINT without Ng
dirs = {}
for name, solver in (("fp", [("SOLVER", "FIXED_POINT")]),
                     ("aa0", [("SOLVER", "ANDERSON"), ("ANDERSON_DEPTH", "0")]),
                     ("aa0_pop", [("SOLVER", "ANDERSON"), ("ANDERSON_DEPTH", "0"),
                                  ("ANDERSON_EXTRAPOLATE", "POPULATIONS")]),
                     ("aa0_rho", [("SOLVER", "ANDERSON"), ("ANDERSON_DEPTH", "0"),
                                  ("ANDERSON_EXTRAPOLATE", "POPULATIONS_RHO")])):
    d, r = run(name, MAP + solver + SHORT + [("N_MAX_ITER", "8")])
    r.require()                               # 1: not converged, expected
    dirs[name] = d
for name in ("aa0", "aa0_pop", "aa0_rho"):
    for f in ("pops.CA.out", "J.dat", "spectrum.out"):
        ok = same_file(dirs["fp"], dirs[name], f)
        ck.report(ok, f"depth 0 {name:8s} {f:13s} "
                      f"{'bit-identical to FIXED_POINT' if ok else 'DIFFERS'}",
                  f"depth 0 {name} {f}")

# --- 2. Extrapolating runs are reproducible
for label, extra in (("defaults", []),
                     ("RELATIVE, PER_BLOCK", [("ANDERSON_VARIABLES", "RELATIVE"),
                                              ("ANDERSON_COEFFICIENTS", "PER_BLOCK")]),
                     ("POPULATIONS_RHO, RELATIVE, PER_BLOCK",
                      [("ANDERSON_EXTRAPOLATE", "POPULATIONS_RHO"),
                       ("ANDERSON_VARIABLES", "RELATIVE"),
                       ("ANDERSON_COEFFICIENTS", "PER_BLOCK")])):
    pair = []
    for k in (1, 2):
        d, r = run(f"aa_{len(extra)}_{k}", MAP + [("SOLVER", "ANDERSON")] + extra +
                   SHORT + [("N_MAX_ITER", "15")])
        r.require()
        pair.append((d, r.stdout))
    m = re.search(r"-- Anderson: extrapolated (\d+)", pair[0][1])
    n_ex = int(m.group(1)) if m else 0
    ok = n_ex > 0 and all(same_file(pair[0][0], pair[1][0], f)
                          for f in ("pops.CA.out", "J.dat", "spectrum.out"))
    ck.report(ok, f"{label}: {n_ex} extrapolations, two runs "
                  f"{'bit-identical' if ok else 'DIFFER (or no extrapolation)'}",
              f"reproducible ({label})")

# --- 3. Conflicting keywords stop the run at input
for label, kws, want in (
        ("ANDERSON_DEPTH with SOLVER = FIXED_POINT",
         MAP + [("SOLVER", "FIXED_POINT"), ("ANDERSON_DEPTH", "3")], "ANDERSON_DEPTH"),
        ("NG_ORDER = 2 with SOLVER = ANDERSON",
         MAP + [("SOLVER", "ANDERSON"), ("NG_ORDER", "2")], "NG_ORDER = 2"),
        ("NG_DELAY with SOLVER = ANDERSON",
         MAP + [("SOLVER", "ANDERSON"), ("NG_DELAY", "10")], "NG_DELAY is set"),
        ("PRD_ITER_LIMIT = 1e-2 with SOLVER = ANDERSON",
         [("SOLVER", "ANDERSON"), ("NG_ORDER", "0"), ("PRD_NG_ORDER", "0"),
          ("PRD_ITER_LIMIT", "1e-2")], "SOLVER = ANDERSON: PRD_ITER_LIMIT"),
        ("ANDERSON_EXTRAPOLATE with SOLVER = FIXED_POINT",
         MAP + [("SOLVER", "FIXED_POINT"), ("ANDERSON_EXTRAPOLATE", "POPULATIONS")],
         "ANDERSON_EXTRAPOLATE is set"),
        ("ANDERSON_EXTRAPOLATE = J",
         MAP + [("SOLVER", "ANDERSON"), ("ANDERSON_EXTRAPOLATE", "J")],
         "Invalid value for keyword ANDERSON_EXTRAPOLATE"),
        ("ANDERSON_SAFEGUARD_WINDOW = 0",
         MAP + [("SOLVER", "ANDERSON"), ("ANDERSON_SAFEGUARD_WINDOW", "0")],
         "ANDERSON_SAFEGUARD_WINDOW = 0 (1 .. 100)")):
    d, r = run("validate", kws)
    ok = r.returncode == EXIT_FATAL and want in r.stdout
    ck.report(ok, f"validate {label}: {'stops at input' if ok else 'NOT STOPPED'}"
                  f" (status {r.returncode})", f"validate {label}")

# --- 4. Reduced extrapolation: what is extrapolated, and the carried blocks
for setting, want in (("POPULATIONS", (492, 1)), ("POPULATIONS_RHO", (None, 3))):
    d, r = run(f"carry_{setting}", MAP + [("SOLVER", "ANDERSON"),
                                          ("ANDERSON_EXTRAPOLATE", setting),
                                          ("ANDERSON_SAFEGUARD", "0.5"),
                                          ("ANDERSON_SAFEGUARD_WINDOW", "1")] +
               SHORT + [("N_MAX_ITER", "15")])
    r.require()
    m = re.search(r"ANDERSON_EXTRAPOLATE = (\w+): (\d+) unknowns in (\d+) blocks", r.stdout)
    rej = re.search(r"rejected (\d+);", r.stdout)
    c = re.search(r"carried blocks equal to G at the last accepted point in (\d+) of"
                  r" (\d+) evaluations", r.stdout)
    ok_banner = (m is not None and m.group(1) == setting and int(m.group(3)) == want[1]
                 and (want[0] is None or int(m.group(2)) == want[0]))
    n_rej = int(rej.group(1)) if rej else 0
    n_eval = len(re.findall(r"-- Iteration +\d+", r.stdout))
    ok_carry = (c is not None and c.group(1) == c.group(2)
                and int(c.group(2)) == n_eval - 1 > 0)
    print(f"{setting}: banner {'ok' if ok_banner else 'WRONG'}"
          f" ({m.group(0) if m else 'missing'}); {n_rej} rejected; carried "
          f"{(c.group(1) + ' of ' + c.group(2)) if c else 'no check line'}")
    ck.check(ok_banner, f"{setting} banner")
    ck.check(n_rej > 0 and ok_carry, f"{setting} carried blocks (rejections {n_rej})")

d, r = run("rho_no_prd", MAP + [("SOLVER", "ANDERSON"),
                               ("ANDERSON_EXTRAPOLATE", "POPULATIONS_RHO"),
                               ("PRD_N_MAX_ITER", "0")] + SHORT + [("N_MAX_ITER", "4")])
ok = (r.ran and "no PRD lines" in r.stdout
      and "the same as POPULATIONS" in r.stdout
      and re.search(r"POPULATIONS_RHO: 492 unknowns in 1 blocks", r.stdout) is not None)
ck.report(ok, f"POPULATIONS_RHO without PRD lines: {'warns and runs' if ok else 'WRONG'}"
              f" (status {r.returncode})", "POPULATIONS_RHO without PRD lines")

ck.finish("solver_anderson", work)
