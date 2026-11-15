#!/usr/bin/env python3
"""End-to-end test: the legacy iteration does not depend on how it was started (N-36).

usage: solver_ng_restart.py --rhf1d PATH [--src DIR] [--work DIR]

SOLVER = NG with SOLVE_NE = ITERATION on FAL-C with H_6_CRD active (all
molecules in the chemical equilibrium), Ng acceleration on, ITER_LIMIT 1e-8:

 1. A fresh run converges (about 95 iterations).
 2. A restart from its populations and mean intensities (STARTING_J = OLD_J,
    which makes the populations OLD_POPULATIONS), same settings, converges
    again (about 70 iterations: the first Background() solves ne from the
    atmosphere file, so the restart starts a few per cent away).
 3. The two converged populations agree to 2e-6 (relative, every level and
    depth point). Both are fixed points of the same equations only when the
    neutral fraction of an active atom in the chemical equilibrium comes
    from its populations in both runs: with the LTE fraction in a fresh
    run, as before N-36, they differ by 6e-6 (H2, H-, ne moved by 1e-6 in
    the model they solve; with the fix 3e-7, the convergence level).

Exit status 0 = pass.
"""
import argparse
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile

p = argparse.ArgumentParser()
p.add_argument("--rhf1d", required=True)
p.add_argument("--src", default=os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
p.add_argument("--work", help="run directory, kept (default: a temporary dir under "
                               "$TMPDIR, removed when the test passes, kept and named "
                               "on failure)")
a = p.parse_args()
src = os.path.abspath(a.src)
work = os.path.abspath(a.work) if a.work else tempfile.mkdtemp(prefix="rh_e2e_")

TOL = 2e-6                      # fresh vs restart: 6.1e-6 before N-36, 2.6e-7 after


def _name_work(*exc):                  # an uncaught exception keeps the runs: say where
    sys.__excepthook__(*exc)
    print(f"runs kept in {work}", file=sys.stderr)


sys.excepthook = _name_work
os.makedirs(os.path.join(work, "runs"), exist_ok=True)
for d in ("Atoms", "Atmos", "Molecules"):          # run dirs sit two levels down
    link = os.path.join(work, d)
    if not os.path.lexists(link):
        os.symlink(os.path.join(src, d), link)


def read_pops(path):
    """The populations n of a pops.*.out file as a flat list (struct only:
    the CI build jobs have no numpy)."""
    b = open(path, "rb").read()
    for off in range(4, min(len(b) - 8, 512), 4):
        nl, ns = struct.unpack(">2i", b[off:off + 8])
        if 0 < nl < 100000 and 0 < ns and off + 8 + 16 * nl * ns == len(b):
            return struct.unpack(f">{nl * ns}d", b[off + 8:off + 8 + 8 * nl * ns])
    raise ValueError(f"cannot parse {path}")


def setkw(text, key, value):
    """Set key = value in a keyword.input text."""
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


def run(name, init, startj, from_dir=None):
    d = os.path.join(work, "runs", name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    kw = open(os.path.join(src, "rhf1d", "run", "keyword.input")).read()
    for k, v in [("ATMOS_FILE", "../../Atmos/FALC_82.atmos"),
                 ("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
                 ("N_THREADS", "1"), ("PRINT_CPU", "FALSE"),
                 ("KURUCZ_DATA", "none"), ("PRD_N_MAX_ITER", "0"),
                 ("STARTING_J", startj), ("SOLVER", "NG"),
                 ("SOLVE_NE", "ITERATION"), ("BACKGROUND_DOUBLE", "TRUE"),
                 ("N_MAX_SCATTER", "50"), ("ITER_LIMIT", "1e-8"),
                 ("N_MAX_ITER", "300"), ("NG_ORDER", "2"), ("NG_DELAY", "10"),
                 ("NG_PERIOD", "3")]:
        kw = setkw(kw, k, v)
    open(os.path.join(d, "keyword.input"), "w").write(kw)
    open(os.path.join(d, "atoms.input"), "w").write(
        f"1\n  ../../Atoms/H_6_CRD.atom  ACTIVE  {init}  pops.H.out\n")
    for f in ("molecules.input", "ray.input"):
        shutil.copy(os.path.join(src, "rhf1d", "run", f), d)
    if from_dir:
        for f in ("pops.H.out", "J.dat"):
            shutil.copy(os.path.join(from_dir, f), d)
    try:
        r = subprocess.run([os.path.abspath(a.rhf1d)], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        sys.exit(f"rhf1d timed out in {d}")
    open(os.path.join(d, "out.log"), "w").write(r.stdout)
    n = len(re.findall(r"^ -- Iteration +\d+", r.stdout, re.M))
    if r.returncode != 0:
        sys.exit(f"rhf1d did not converge in {d} (status {r.returncode}, "
                 f"{n} iterations)")
    return d, n


failed = []
d1, n1 = run("fresh", "ZERO_RADIATION", "NEW_J")
d2, n2 = run("restart", "OLD_POPULATIONS", "OLD_J", from_dir=d1)
p1, p2 = read_pops(os.path.join(d1, "pops.H.out")), read_pops(os.path.join(d2, "pops.H.out"))
err = max(abs(x - y) / abs(x) for x, y in zip(p1, p2))
print(f"fresh run {n1} iterations, restart {n2} iterations; populations differ by "
      f"{err:.2e} (limit {TOL:.0e})")
if not err < TOL:
    failed.append(f"fresh and restarted runs differ by {err:.2e}")

if failed:
    sys.exit("FAILED: " + "; ".join(failed) + f"\n runs in {work}")
if not a.work:
    shutil.rmtree(work)
print("PASS")
