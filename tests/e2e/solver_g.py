#!/usr/bin/env python3
"""End-to-end test of the map G and SOLVER = FIXED_POINT (S-2, docs/SOLVER.md).

usage: solver_g.py --rhf1d PATH --gcheck1d PATH [--gcheck2d PATH]
                   [--src DIR] [--work DIR]

Ca II with angle-averaged PRD (two PRD lines) on FAL-C, H_6 passive:

 1. history: SOLVER = FIXED_POINT with NG_ORDER = 0 reproduces the default
    solver with NG_ORDER = PRD_NG_ORDER = PRD_ITER_LIMIT = 0 bit for bit
    (populations, J, spectrum) after 5 iterations.
 2. gcheck (rhf1d, and rhsc2d at 2 threads with --gcheck2d): G(x_a) twice
    gives bit-identical results, and G(x_b) between two calls of G(x_a)
    changes nothing.
 3. validation: SOLVER = FIXED_POINT with the default PRD_ITER_LIMIT stops
    at input with a message that names PRD_ITER_LIMIT.

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
p.add_argument("--gcheck1d", required=True)
p.add_argument("--gcheck2d")
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

G_KEYWORDS = [("SOLVER", "FIXED_POINT"), ("NG_ORDER", "0"),
              ("PRD_NG_ORDER", "0"), ("PRD_ITER_LIMIT", "0")]


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


def prepare(name, geom, keywords):
    rundir = {"1d": "rhf1d", "2d": "rhsc2d"}[geom]
    atmos = {"1d": "../../Atmos/FALC_82.atmos", "2d": "../../Atmos/FALC_03x82.atmos"}[geom]
    d = os.path.join(work, "runs", name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    kw = open(os.path.join(src, rundir, "run", "keyword.input")).read()
    for k, v in [("ATMOS_FILE", atmos), ("STOKES_MODE", "NO_STOKES"),
                 ("STOKES_INPUT", "none"), ("PRD_N_MAX_ITER", "2"),
                 ("PRD_ANGLE_DEP", "PRD_ANGLE_AVER"), ("N_THREADS", "1"),
                 ("PRINT_CPU", "FALSE"), ("STARTING_J", "NEW_J"),
                 ("KURUCZ_DATA", "none"), ("SOLVE_NE", "NONE"),
                 ("COLL_RELAX_NSTEP", "0")] + keywords:
        kw = setkw(kw, k, v)
    open(os.path.join(d, "keyword.input"), "w").write(kw)
    open(os.path.join(d, "atoms.input"), "w").write(
        "2\n  ../../Atoms/H_6.atom  PASSIVE  LTE_POPULATIONS\n"
        "  ../../Atoms/CaII.atom  ACTIVE  ZERO_RADIATION  pops.CA.out\n")
    for f in ("molecules.input", "ray.input"):
        if os.path.exists(os.path.join(src, rundir, "run", f)):
            shutil.copy(os.path.join(src, rundir, "run", f), d)
    return d


def execute(program, d, args=()):
    try:
        r = subprocess.run([program, *args], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        sys.exit(f"{program} timed out in {d}")
    open(os.path.join(d, "out.log"), "w").write(r.stdout)
    return r


failed = []

# --- 1. Iteration history of FIXED_POINT = default solver without Ng
runs = {}
for name, extra in (("hist_ng", [("SOLVER", "NG")] + G_KEYWORDS[1:]),
                    ("hist_fp", G_KEYWORDS)):
    d = prepare(name, "1d", extra + [("N_MAX_ITER", "5"), ("ITER_LIMIT", "1e-14")])
    r = execute(os.path.abspath(a.rhf1d), d)
    if r.returncode not in (0, 1):            # 1: not converged, expected
        sys.exit(f"rhf1d failed in {d} (status {r.returncode})")
    runs[name] = d
for f in ("pops.CA.out", "J.dat", "spectrum.out"):
    same = open(os.path.join(runs["hist_ng"], f), "rb").read() == \
        open(os.path.join(runs["hist_fp"], f), "rb").read()
    print(f"history  {f:13s} {'bit-identical' if same else 'DIFFERS'}")
    if not same:
        failed.append(f"history {f}")

# --- 2. gcheck: twice and history independence
checks = [("gcheck_1d", "1d", a.gcheck1d, "1")]
if a.gcheck2d:
    checks.append(("gcheck_2d", "2d", a.gcheck2d, "2"))
for name, geom, program, threads in checks:
    d = prepare(name, geom, G_KEYWORDS + [("N_THREADS", threads)])
    r = execute(os.path.abspath(program), d, ["--steps", "3", "--lag", "1"])
    lines = [l for l in r.stdout.splitlines() if l.startswith("gcheck: twice")
             or l.startswith("gcheck: history")]
    if r.returncode != 0 or len(lines) != 2:
        tail = "\n".join(r.stdout.splitlines()[-30:])
        sys.exit(f"{program} failed in {d} (status {r.returncode}):\n{tail}")
    for l in lines:
        print(f"{name} ({threads} thread{'s' if threads != '1' else ''}): {l[8:]}")
        if "bit-identical" not in l:
            failed.append(f"{name}: {l}")

# --- 3. Conflicting keywords stop the run at input
d = prepare("validate", "1d", [("SOLVER", "FIXED_POINT")])
r = execute(os.path.abspath(a.rhf1d), d)
ok = r.returncode == 70 and "PRD_ITER_LIMIT" in r.stdout and \
    "reading input file" not in r.stdout.split("TERMINATING")[0].split("readInput")[-1]
print(f"validate SOLVER = FIXED_POINT with PRD_ITER_LIMIT = 1e-2: "
      f"{'stops at input' if ok else 'NOT STOPPED'} (status {r.returncode})")
if not ok:
    failed.append("validation")

if failed:
    sys.exit("FAILED: " + "; ".join(failed))
print("PASS")
