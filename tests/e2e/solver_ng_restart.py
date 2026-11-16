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


def run(name, init, startj, from_dir=None, atoms=None, extra=(), nmax="300",
        converge=True):
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
                 ("N_MAX_ITER", nmax), ("NG_ORDER", "2"), ("NG_DELAY", "10"),
                 ("NG_PERIOD", "3"), *extra]:
        kw = setkw(kw, k, v)
    open(os.path.join(d, "keyword.input"), "w").write(kw)
    open(os.path.join(d, "atoms.input"), "w").write(
        atoms or f"1\n  ../../Atoms/H_6_CRD.atom  ACTIVE  {init}  pops.H.out\n")
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
    if r.returncode != 0 and (converge or r.returncode != 1):
        sys.exit(f"rhf1d failed in {d} (status {r.returncode}, {n} iterations)")
    return d, n, r.stdout


failed = []
d1, n1, _ = run("fresh", "ZERO_RADIATION", "NEW_J")
d2, n2, _ = run("restart", "OLD_POPULATIONS", "OLD_J", from_dir=d1)
p1, p2 = read_pops(os.path.join(d1, "pops.H.out")), read_pops(os.path.join(d2, "pops.H.out"))
if len(p1) != len(p2) or not p1:
    sys.exit(f"population files differ in size: {len(p1)} and {len(p2)} values")
err = max(abs(x - y) / max(abs(x), 1e-300) for x, y in zip(p1, p2))
print(f"fresh run {n1} iterations, restart {n2} iterations; populations differ by "
      f"{err:.2e} (limit {TOL:.0e})")
if not err < TOL:
    failed.append(f"fresh and restarted runs differ by {err:.2e}")

# --- 2. setup before initSolution: no NaN, the run iterates
for name, atoms, extra in (
        ("setup_C", "2\n  ../../Atoms/H_6.atom  PASSIVE  LTE_POPULATIONS\n"
                    "  ../../Atoms/C_I+II_9.atom  ACTIVE  ZERO_RADIATION  pops.C.out\n", ()),
        ("setup_HLTE", None, (("HYDROGEN_LTE", "TRUE"),))):
    d, n, out = run(name, "ZERO_RADIATION", "NEW_J", atoms=atoms, extra=extra,
                    nmax="2", converge=False)
    bad = n < 2 or re.search(r"\bnan\b|FATAL|TERMINATING", out, re.I) is not None
    print(f"setup    {name:10s} {n} iterations: {'FAILED' if bad else 'ok'}")
    if bad:
        failed.append(f"{name}: {n} iterations, NaN or error in the output")

if failed:
    sys.exit("FAILED: " + "; ".join(failed) + f"\n runs in {work}")
if not a.work:
    shutil.rmtree(work)
print("PASS")
