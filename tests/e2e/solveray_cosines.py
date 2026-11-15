#!/usr/bin/env python3
"""End-to-end test of the direction cosines of solveray, B-34.

usage: solveray_cosines.py --rhsc2d P --solveray2d P --solveray3d P
                           [--src DIR] [--work DIR] [-v]

solveray reads two cosines from ray.input and derives the third from
mux^2 + muy^2 + muz^2 = 1. A sum of squares just over 1 (an input rounded
up, or FMA contraction) made the root a NaN and the run completed:
polarized 2-D results were NaN (muy only enters Bproject), in 3-D muz was
NaN (a crash here). 2-D runs a short FAL-C case with a magnetic field
(make_field.py), then solveray for:

 2-D: mux = 0.8541662602, muz = 0.52, whose squares sum to 1 + 6e-11 in
      double precision: finite I, Q, U, V, and equal to the run with
      mux = 0.8541662601 to 2e-4 of the peak (muy = 0 against muy = 5e-6,
      so not to rounding); mux = muz = 0.9 is an error that names the
      sum. A valid ray gives finite output.
 3-D: mux = 0.8660254038, muy = 0.5 (sum 1 + 3e-11) and mux = 0.6,
      muy = 0.8 (muz = 0) are errors, since muz must be positive (the
      stencil distances divide by it). These stop before the results of a
      main run are read, so they run in a copy of rhsc3d/run (a 3-D main
      run would take ~40 s).

Python 3 standard library only. Exit status 0 = pass.
"""
import argparse
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile

here = os.path.dirname(os.path.abspath(__file__))
p = argparse.ArgumentParser()
for k in ("rhsc2d", "solveray2d", "solveray3d"):
    p.add_argument("--" + k, required=True)
p.add_argument("--src", default=os.path.dirname(os.path.dirname(here)),
               help="RH source tree (Atoms/, Atmos/, Molecules/, rhsc*/run/)")
p.add_argument("--work", help="run directory, kept (default: a temporary dir under "
                               "$TMPDIR, removed when the test passes, kept and named "
                               "on failure)")
p.add_argument("-v", "--verbose", action="store_true")
a = p.parse_args()
src = os.path.abspath(a.src)
work = os.path.abspath(a.work) if a.work else tempfile.mkdtemp(prefix="rh_e2e_")


def _name_work(*exc):                  # an uncaught exception keeps the runs: say where
    sys.__excepthook__(*exc)
    print(f"runs kept in {work}", file=sys.stderr)


sys.excepthook = _name_work
os.makedirs(os.path.join(work, "runs"), exist_ok=True)
for d in ("Atoms", "Atmos", "Molecules"):          # run dirs sit two levels down
    link = os.path.join(work, d)
    if not os.path.lexists(link):
        os.symlink(os.path.join(src, d), link)

failed = nchecks = 0


def check(ok, what, detail=""):
    global failed, nchecks
    nchecks += 1
    failed += not ok
    if a.verbose or not ok:
        print(f"{'PASS' if ok else 'FAIL'}  {what} {detail}")


def execute(prog, d, timeout=600):
    r = subprocess.run([os.path.abspath(prog)], cwd=d, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, timeout=timeout)
    return r.returncode, r.stdout


def main_run(geo, rhsc):
    """A short run whose directory solveray can start from."""
    d = os.path.join(work, "runs", geo)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(os.path.join(src, f"rhsc{geo}", "run"), d)
    atmos = {"2d": "FALC_03x82", "3d": "FALC_04x04x82"}[geo]
    bfile = os.path.join(work, f"field_{geo}.B")
    subprocess.run([sys.executable, "-I", os.path.join(src, "tests", "make_field.py"),
                    geo, os.path.join(src, "Atmos", atmos + ".atmos"), bfile],
                   check=True, stdout=subprocess.DEVNULL)
    with open(os.path.join(d, "keyword.input"), "a") as f:
        f.write(f"  STOKES_INPUT = {bfile}\n  STOKES_MODE = FIELD_FREE\n"
                "  N_MAX_ITER = 2\n  N_THREADS = 1\n")
    rc, out = execute(rhsc, d)
    if rc not in (0, 1) or "TERMINATING_ERROR" in out:    # 1: not converged
        sys.exit(f"rhsc{geo} failed in {d} (status {rc}):\n" + out[-2000:])
    return d


def solveray(prog, rundir, name, c1, c2):
    d = os.path.join(work, "runs", name)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(rundir, d)
    with open(os.path.join(d, "ray.input"), "w") as f:
        f.write(f"{c1} {c2}\n0\n")
    rc, out = execute(prog, d)
    files = [x for x in os.listdir(d) if x.startswith("spectrum_")]
    return rc, out, (os.path.join(d, files[0]) if files else None)


def read_ray(path):
    """I, Q, U, V of a solveray output: mux, mu (2 doubles), I, the int
    Nspect = 0 of ray.input, then Q, U, V, all big-endian XDR."""
    b = open(path, "rb").read()
    M = (len(b) - 20) // 32
    assert 32 * M + 20 == len(b), "unexpected solveray output layout"
    off = 16 + 8 * M + 4
    return [struct.unpack(f">{M}d", b[o:o + 8 * M])
            for o in (16, off, off + 8 * M, off + 16 * M)]


def finite(v):
    return all(math.isfinite(x) for s in v for x in s)


def maxdiff(u, v):
    """Largest difference of any of I, Q, U, V, each relative to its own peak."""
    return max(max(abs(x - y) for x, y in zip(s, t)) / max(max(map(abs, t)), 1e-300)
               for s, t in zip(u, v))


d2 = main_run("2d", a.rhsc2d)
rc, out, f_up = solveray(a.solveray2d, d2, "ray2d_up", "0.8541662602", "0.52")
check(rc == 0 and f_up and finite(read_ray(f_up)),
      "2-D: sum of squares 1 + 6e-11 gives finite I, Q, U, V", f"(status {rc})")
rc, out, f_dn = solveray(a.solveray2d, d2, "ray2d_dn", "0.8541662601", "0.52")
if f_up and f_dn:
    err = maxdiff(read_ray(f_up), read_ray(f_dn))
    check(err < 2e-4, "2-D: clamped ray equals the ray with the cosine rounded down",
          f"({err:.1e})")
rc, out, f = solveray(a.solveray2d, d2, "ray2d_bad", "0.9", "0.9")
check(rc != 0 and "mux^2 + muz^2 = 1.62" in out and f is None,
      "2-D: mux = muz = 0.9 is an error naming the sum", f"(status {rc})")
rc, out, f = solveray(a.solveray2d, d2, "ray2d_ok", "0.5", "0.5")
check(rc == 0 and f and finite(read_ray(f)), "2-D: valid ray is finite", f"(status {rc})")

d3 = os.path.join(work, "runs", "3d")
shutil.rmtree(d3, ignore_errors=True)
shutil.copytree(os.path.join(src, "rhsc3d", "run"), d3)
for name, c1, c2 in (("up", "0.8660254038", "0.5"), ("flat", "0.6", "0.8"),
                     ("bad", "0.9", "0.9")):
    rc, out, f = solveray(a.solveray3d, d3, "ray3d_" + name, c1, c2)
    check(rc != 0 and "no upward component" in out and f is None,
          f"3-D: mux = {c1}, muy = {c2} is an error (muz must be positive)", f"(status {rc})")

print(f"solveray_cosines: {nchecks} checks, {failed} failed")
if failed:
    print(f"runs kept in {work}")
elif not a.work:
    shutil.rmtree(work)
sys.exit(1 if failed else 0)
