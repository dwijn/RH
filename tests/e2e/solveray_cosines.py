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
      main run are read, so they run in a copy of examples/rhsc3d (a 3-D main
      run would take ~40 s).

Python 3 standard library only. Exit status 0 = pass.
"""
import math
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker, maxabs, maxdiff
from rhtest.cli import parse
from rhtest.output import read_solveray
from rhtest.run import copy_dir, execute, make_field, template
from rhtest.workdir import WorkDir

a = parse(required=("rhsc2d", "solveray2d", "solveray3d"), verbose=True,
          src_help="RH source tree (Atoms/, Atmos/, Molecules/, examples/)")
src = a.src
work = WorkDir(src, a.work)
ck = Checker(a.verbose)


def main_run(geo, rhsc):
    """A short run whose directory solveray can start from."""
    d = copy_dir(work, geo, template(src, f"rhsc{geo}"))
    atmos = {"2d": "FALC_03x82", "3d": "FALC_04x04x82"}[geo]
    bfile = work.join(f"field_{geo}.B")
    make_field(src, geo, os.path.join(src, "Atmos", atmos + ".atmos"), bfile)
    with open(os.path.join(d, "keyword.input"), "a") as f:
        f.write(f"  STOKES_INPUT = {bfile}\n  STOKES_MODE = FIELD_FREE\n"
                "  N_MAX_ITER = 2\n  N_THREADS = 1\n")
    execute(rhsc, d).require(bad_text="TERMINATING_ERROR")    # 1: not converged
    return d


def solveray(prog, rundir, name, c1, c2):
    d = copy_dir(work, name, rundir)
    with open(os.path.join(d, "ray.input"), "w") as f:
        f.write(f"{c1} {c2}\n0\n")
    r = execute(prog, d)
    rc, out = r.returncode, r.stdout
    files = [x for x in os.listdir(d) if x.startswith("spectrum_")]
    return rc, out, (os.path.join(d, files[0]) if files else None)


def finite(v):
    return all(math.isfinite(x) for s in v for x in s)


def ray_diff(u, v):
    """Largest difference of any of I, Q, U, V, each relative to its own peak."""
    return max(maxdiff(s, t) / max(maxabs(t), 1e-300) for s, t in zip(u, v))


d2 = main_run("2d", a.rhsc2d)
rc, out, f_up = solveray(a.solveray2d, d2, "ray2d_up", "0.8541662602", "0.52")
ck.check(rc == 0 and f_up and finite(read_solveray(f_up)),
         "2-D: sum of squares 1 + 6e-11 gives finite I, Q, U, V", f"(status {rc})")
rc, out, f_dn = solveray(a.solveray2d, d2, "ray2d_dn", "0.8541662601", "0.52")
if f_up and f_dn:
    err = ray_diff(read_solveray(f_up), read_solveray(f_dn))
    ck.check(err < 2e-4, "2-D: clamped ray equals the ray with the cosine rounded down",
             f"({err:.1e})")
rc, out, f = solveray(a.solveray2d, d2, "ray2d_bad", "0.9", "0.9")
ck.check(rc != 0 and "mux^2 + muz^2 = 1.62" in out and f is None,
         "2-D: mux = muz = 0.9 is an error naming the sum", f"(status {rc})")
rc, out, f = solveray(a.solveray2d, d2, "ray2d_ok", "0.5", "0.5")
ck.check(rc == 0 and f and finite(read_solveray(f)), "2-D: valid ray is finite", f"(status {rc})")

d3 = copy_dir(work, "3d", template(src, "rhsc3d"))
for name, c1, c2 in (("up", "0.8660254038", "0.5"), ("flat", "0.6", "0.8"),
                     ("bad", "0.9", "0.9")):
    rc, out, f = solveray(a.solveray3d, d3, "ray3d_" + name, c1, c2)
    ck.check(rc != 0 and "no upward component" in out and f is None,
             f"3-D: mux = {c1}, muy = {c2} is an error (muz must be positive)", f"(status {rc})")

ck.finish("solveray_cosines", work)
