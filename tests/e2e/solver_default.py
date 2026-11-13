#!/usr/bin/env python3
"""End-to-end test of the default solver (S-10, docs/SOLVER.md, Anderson as
the default solver).

usage: solver_default.py --rhf1d PATH --rhsphere PATH --rhsc2d PATH
                         [--src DIR] [--work DIR]

Each program runs its shipped example (run/) for one or two iterations:
 1. rhf1d, rhsphere and rhsc2d iterate with ANDERSON when the keyword file
    does not set SOLVER, and say so (" -- SOLVER = ANDERSON (the default)").
    rhf1d and rhsphere extrapolate the whole state (ANDERSON_EXTRAPOLATE =
    ALL), rhsc2d only the populations and rho (POPULATIONS_RHO, S-11).
 2. A 1-D keyword file written for the old default (Ng keywords, no SOLVER)
    stops at input (status 70) with messages that name the default and
    SOLVER = NG as a fix; with SOLVER = NG added it runs, and the banner
    does not call NG the default.
Exit status 0 = pass.
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

p = argparse.ArgumentParser()
p.add_argument("--rhf1d", required=True)
p.add_argument("--rhsphere", required=True)
p.add_argument("--rhsc2d", required=True)
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


def setkw(text, key, value):
    out, done = [], False
    for l in text.splitlines():
        s = l.strip()
        if not s.startswith("#") and "=" in s and s.split("=")[0].strip() == key:
            out.append(f"  {key} = {value}")
            done = True
        else:
            out.append(l)
    if not done:
        out.append(f"  {key} = {value}")
    return "\n".join(out) + "\n"


def run(name, program, rundir, keywords):
    d = os.path.join(work, "runs", name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    shutil.copytree(os.path.join(src, rundir, "run"), d)
    kw = open(os.path.join(d, "keyword.input")).read()
    for k, v in [("N_THREADS", "1"), ("PRINT_CPU", "FALSE")] + keywords:
        kw = setkw(kw, k, v)
    open(os.path.join(d, "keyword.input"), "w").write(kw)
    try:
        r = subprocess.run([os.path.abspath(program)], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        sys.exit(f"{program} timed out in {d}")
    open(os.path.join(d, "out.log"), "w").write(r.stdout)
    return r


failed = []


def check(label, ok, r):
    print(f"{label}: {'ok' if ok else 'FAILED'} (status {r.returncode})")
    if not ok:
        failed.append(label)


# --- 1. The default of each geometry
for name, program, rundir, want in (
        ("rhf1d", a.rhf1d, "rhf1d", "ALL"),
        ("rhsphere", a.rhsphere, "rhsphere", "ALL"),
        ("rhsc2d", a.rhsc2d, "rhsc2d", "POPULATIONS_RHO")):
    r = run(f"default_{name}", program, rundir, [("N_MAX_ITER", "2")])
    ok = r.returncode in (0, 1) and " -- SOLVER = ANDERSON (the default)\n" in r.stdout \
        and "-- Anderson: extrapolated" in r.stdout \
        and f"ANDERSON_EXTRAPOLATE = {want}:" in r.stdout
    check(f"{name} without SOLVER iterates with ANDERSON, {want} (the default)", ok, r)

# --- 2. A keyword file written for the old default
OLD = [("NG_ORDER", "2"), ("NG_DELAY", "10"), ("NG_PERIOD", "3"),
       ("PRD_ITER_LIMIT", "1.0E-2")]
r = run("old_file", a.rhf1d, "rhf1d", OLD)
ok = r.returncode == 70 and "SOLVER = ANDERSON (the default)" in r.stdout and \
    "Set SOLVER = NG" in r.stdout and "4 keyword settings conflict" in r.stdout
check("rhf1d, old keyword file without SOLVER stops at input", ok, r)
r = run("old_file_ng", a.rhf1d, "rhf1d", OLD + [("SOLVER", "NG"), ("N_MAX_ITER", "2")])
ok = r.returncode in (0, 1) and " -- SOLVER = NG\n" in r.stdout
check("rhf1d, old keyword file with SOLVER = NG runs with NG", ok, r)

if failed:
    sys.exit("FAILED: " + "; ".join(failed) + f"\n runs in {work}")
if not a.work:
    shutil.rmtree(work)
print("solver_default: all checks passed")
