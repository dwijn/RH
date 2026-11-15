#!/usr/bin/env python3
"""End-to-end test of the convergence report of the G solvers (N-3).

usage: convergence_report.py --rhf1d PATH [--src DIR] [--work DIR]

Ca II with angle-averaged PRD (two PRD lines) on FAL-C, H_6 passive, as
solver_anderson.py:
 1. The per-iteration lines carry the residual of the rho blocks (" rho
    CA: delta = ...") and, with SOLVE_NE = ITERATION, that of ne, with
    all three G solvers.
 2. Each G solver ends with " -- Estimated contraction ..." after the
    convergence message, which keeps its format; with SOLVER = NG there is
    no such line; with fewer than three iterations the line says "no
    estimate".
 3. The error estimate (residual / (1 - contraction)) at ITER_LIMIT = 1e-4
    is of the size of the actual difference to a run of the same case at
    ITER_LIMIT = 1e-9 of ANDERSON (the reference). For FIXED_POINT without
    Ng the estimate is within a factor 3 of the difference: the contraction
    is that of G. For ANDERSON and NEWTON_KRYLOV it is an underestimate by
    up to an order of magnitude (the rate of the accelerated iteration is
    not that of G, and the residuals are not monotone), so only a loose
    check is made (factor 30). The numbers are printed.
That the reporting leaves the results alone is checked by run_tests.sh
--exact against a build without it.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--rhf1d", required=True)
p.add_argument("--src", default=os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
p.add_argument("--work", help="keep the runs here (default: a temporary dir)")
a = p.parse_args()
src = os.path.abspath(a.src)
work = os.path.abspath(a.work) if a.work else tempfile.mkdtemp(prefix="rh_e2e_")
os.makedirs(os.path.join(work, "runs"), exist_ok=True)
for d in ("Atoms", "Atmos", "Molecules"):          # run dirs sit two levels down
    link = os.path.join(work, d)
    if not os.path.lexists(link):
        os.symlink(os.path.join(src, d), link)
sys.path.insert(0, os.path.join(src, "tests"))
from analyze_run import read_pops  # noqa: E402

MAP = [("NG_ORDER", "0"), ("NG_DELAY", None), ("NG_PERIOD", None),
       ("PRD_NG_ORDER", "0"), ("PRD_ITER_LIMIT", "0")]


def setkw(text, key, value):
    """Set key = value in a keyword.input text; value None removes the key
    (comments it out), so that it is not set at all."""
    out, done = [], False
    for l in text.splitlines():
        s = l.strip()
        if not s.startswith("#") and "=" in s and s.split("=")[0].strip() == key:
            out.append(f"# {s}" if value is None else f"  {key} = {value}")
            done = True
        else:
            out.append(l)
    if not done and value is not None:
        out.append(f"  {key} = {value}")
    return "\n".join(out) + "\n"


def run(name, keywords):
    d = os.path.join(work, "runs", name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    kw = open(os.path.join(src, "rhf1d", "run", "keyword.input")).read()
    for k, v in [("ATMOS_FILE", "../../Atmos/FALC_82.atmos"),
                 ("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
                 ("PRD_N_MAX_ITER", "2"), ("PRD_ANGLE_DEP", "PRD_ANGLE_AVER"),
                 ("N_THREADS", "1"), ("PRINT_CPU", "FALSE"),
                 ("STARTING_J", "NEW_J"), ("KURUCZ_DATA", "none"),
                 ("SOLVE_NE", "NONE"), ("COLL_RELAX_NSTEP", "0")] + keywords:
        kw = setkw(kw, k, v)
    open(os.path.join(d, "keyword.input"), "w").write(kw)
    open(os.path.join(d, "atoms.input"), "w").write(
        "2\n  ../../Atoms/H_6.atom  PASSIVE  LTE_POPULATIONS\n"
        "  ../../Atoms/CaII.atom  ACTIVE  ZERO_RADIATION  pops.CA.out\n")
    for f in ("molecules.input", "ray.input"):
        if os.path.exists(os.path.join(src, "rhf1d", "run", f)):
            shutil.copy(os.path.join(src, "rhf1d", "run", f), d)
    try:
        r = subprocess.run([os.path.abspath(a.rhf1d)], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=900)
    except subprocess.TimeoutExpired:
        sys.exit(f"rhf1d timed out in {d}")
    open(os.path.join(d, "out.log"), "w").write(r.stdout)
    if r.returncode not in (0, 1):
        sys.exit(f"rhf1d failed in {d} (status {r.returncode})")
    return d, r.stdout


EST = re.compile(r" -- Estimated contraction ([0-9.]+) per iteration \(last (\d+) "
                 r"iterations\); the error of the populations may be up to about "
                 r"([0-9.E+-]+) \(residual / \(1 - contraction\)\)\n")
failed = []
SOLVERS = {"FIXED_POINT": [("SOLVER", "FIXED_POINT")],
           "ANDERSON": [("SOLVER", "ANDERSON")],
           "NEWTON_KRYLOV": [("SOLVER", "NEWTON_KRYLOV")]}

# --- 1, 2. rho and ne lines, the contraction line after the convergence message
for name, kws in SOLVERS.items():
    d, out = run(f"lines_{name}", MAP + kws + [("SOLVE_NE", "ITERATION"),
                                               ("ITER_LIMIT", "1e-3"),
                                               ("N_MAX_ITER", "60")])
    n_it = len(re.findall(r"\n -- Iteration +\d+\n", out))
    n_rho = len(re.findall(r"^ rho CA: delta = [0-9.E+-]+$", out, re.M))
    n_ne = len(re.findall(r"^ ne: delta = [0-9.E+-]+$", out, re.M))
    conv = re.search(r"\n -- Converged (in \d+ iterations|at base point \d+): max change "
                     r"[0-9.E+-]+ < ITER_LIMIT = [0-9.E+-]+\n(.*\n)", out)
    est = EST.match(conv.group(2)) if conv else None
    ok = n_it > 3 and n_rho == n_it and n_ne == n_it and est is not None
    print(f"{name}: {n_it} iterations, {n_rho} rho and {n_ne} ne lines; "
          f"{est.group(0).strip() if est else 'NO estimate line after the Converged line'}")
    if not ok:
        failed.append(f"report lines ({name})")

d, out = run("ng", MAP + [("SOLVER", "NG"), ("ITER_LIMIT", "1e-3"), ("N_MAX_ITER", "60")])
ok = "Converged in" in out and "Estimated contraction" not in out
print(f"NG: {'no estimate line' if ok else 'ESTIMATE LINE PRESENT or no convergence'}")
if not ok:
    failed.append("NG prints no estimate")

d, out = run("short", MAP + [("SOLVER", "FIXED_POINT"), ("ITER_LIMIT", "1e-14"),
                             ("N_MAX_ITER", "2")])
ok = "NOT CONVERGED after 2 iterations" in out and \
    "Estimated contraction: no estimate (fewer than 3" in out
print(f"two iterations: {'no estimate' if ok else 'WRONG'}")
if not ok:
    failed.append("no estimate with two iterations")


# --- 3. The estimate against the actual difference to a tight run
def err(d1, d2):
    n1, n2 = read_pops(os.path.join(d1, "pops.CA.out")), \
        read_pops(os.path.join(d2, "pops.CA.out"))
    nz = n2 != 0
    return float(np.max(np.abs((n1 - n2)[nz] / n2[nz])))


ref, out_ref = run("ref_ANDERSON_1e-9", MAP + SOLVERS["ANDERSON"] +
                   [("N_MAX_ITER", "400"), ("ITER_LIMIT", "1e-9")])
if not re.search(r" -- Converged (in|at)", out_ref):
    sys.exit(f"the reference run did not converge ({ref})")
for name, kws, factor in (("FIXED_POINT", SOLVERS["FIXED_POINT"], 3.0),
                          ("ANDERSON", SOLVERS["ANDERSON"], 30.0),
                          ("NEWTON_KRYLOV", SOLVERS["NEWTON_KRYLOV"], 30.0)):
    d4, out4 = run(f"{name}_1e-4", MAP + kws + [("N_MAX_ITER", "400"),
                                                 ("ITER_LIMIT", "1e-4")])
    est = EST.search(out4)
    if not re.search(r" -- Converged (in|at)", out4) or est is None:
        failed.append(f"{name} 1e-4 run")
        continue
    lam, e_est, e_act = float(est.group(1)), float(est.group(3)), err(d4, ref)
    r = e_est / e_act if e_act > 0 else float("inf")
    print(f"{name}: contraction {lam:.3f}, estimated error {e_est:.2E}, "
          f"actual |x(1e-4) - x(1e-9)| {e_act:.2E}, ratio {r:.2f}")
    if not (1.0 / factor < r < factor):
        failed.append(f"{name} estimate {e_est:.2E} vs actual {e_act:.2E}")

if failed:
    sys.exit("FAILED: " + "; ".join(failed) + f"\n runs in {work}")
if not a.work:
    shutil.rmtree(work)
print("convergence_report: all checks passed")
