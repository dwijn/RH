#!/usr/bin/env python3
"""End-to-end test of hybrid PRD (PRD_ANGLE_APPROX), B-31.

usage: prd_hybrid.py --rhf1d PATH [--rhsc2d PATH] [--src DIR] [--work DIR] [-v]

In a static atmosphere the gas frame is the observer's frame, so the
gas-frame mean intensity Jgas of the hybrid approximation equals J, and
hybrid PRD gives the same answer as angle-averaged PRD. Ca II with PRD,
three iterations, rhf1d on FAL-C with one ray and with three rays, and
(with --rhsc2d) rhsc2d on three FAL-C columns with one ray and with six:

 1. Hybrid run: Jgas = J at every wavelength and depth.
 2. J of the hybrid run equals J of the angle-averaged run.

Both hold to rounding. Before the fix of B-31, with one ray, Jgas grew by
J with every formal solution at the wavelengths without PRD line opacity
and stayed zero at those of the PRD lines; in 2-D and 3-D it was 2 J at
the wavelengths without PRD line opacity. Exit status 0 = pass.
"""
import argparse
import os
import shutil
import struct
import subprocess
import sys
import tempfile

TOL = 1e-12

p = argparse.ArgumentParser()
p.add_argument("--rhf1d", required=True)
p.add_argument("--rhsc2d", help="also test rhsc2d")
p.add_argument("--src", default=os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))),
    help="RH source tree (Atoms/, Atmos/, Molecules/, rhf1d/run/)")
p.add_argument("--work", help="keep the runs here (default: a temporary dir)")
p.add_argument("-v", "--verbose", action="store_true")
a = p.parse_args()
rhf1d = os.path.abspath(a.rhf1d)
rhsc2d = os.path.abspath(a.rhsc2d) if a.rhsc2d else None
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


def read_doubles(path):
    raw = open(path, "rb").read()
    return struct.unpack(f"={len(raw) // 8}d", raw)


def run(name, prog, angles, mode):
    """prog: (program, run dir in the source tree, atmosphere); angles:
    keyword and value. Returns J and, for the hybrid mode, Jgas."""
    program, rundir, atmos = prog
    d = os.path.join(work, "runs", name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    kw = open(os.path.join(src, rundir, "run", "keyword.input")).read()
    for k, v in [("ATMOS_FILE", atmos),
                 ("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
                 angles, ("N_MAX_ITER", "3"),
                 ("PRD_N_MAX_ITER", "2"), ("PRD_ANGLE_DEP", mode),
                 ("N_THREADS", "1"), ("PRINT_CPU", "FALSE"),
                 ("STARTING_J", "NEW_J"), ("KURUCZ_DATA", "none"),
                 ("SOLVE_NE", "NONE")]:
        kw = setkw(kw, k, v)
    open(os.path.join(d, "keyword.input"), "w").write(kw)
    open(os.path.join(d, "atoms.input"), "w").write(
        "2\n  ../../Atoms/H_6.atom  PASSIVE  LTE_POPULATIONS\n"
        "  ../../Atoms/CaII.atom  ACTIVE  ZERO_RADIATION  pops.CA.out\n")
    for f in ("molecules.input", "ray.input"):
        if os.path.exists(os.path.join(src, rundir, "run", f)):
            shutil.copy(os.path.join(src, rundir, "run", f), d)
    try:
        r = subprocess.run([program], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=600)
    except subprocess.TimeoutExpired:
        sys.exit(f"{program} timed out in {d}")
    open(os.path.join(d, "out.log"), "w").write(r.stdout)
    # --- Three iterations do not converge: exit status 1 is expected
    if r.returncode not in (0, 1) or "ERROR" in r.stdout:
        sys.exit(f"{program} failed in {d} (status {r.returncode}); see out.log")
    J = read_doubles(os.path.join(d, "J.dat"))
    Jgas = read_doubles(os.path.join(d, "Jgas.dat")) if mode == "PRD_ANGLE_APPROX" else None
    return J, Jgas


def maxrel(x, y):
    return max(abs(u - v) / max(abs(u), abs(v), 1e-300) for u, v in zip(x, y))


cases = [("1d", (rhf1d, "rhf1d", "../../Atmos/FALC_82.atmos"), "1 ray", ("NRAYS", "1")),
         ("1d", (rhf1d, "rhf1d", "../../Atmos/FALC_82.atmos"), "3 rays", ("NRAYS", "3"))]
if rhsc2d:
    prog2 = (rhsc2d, "rhsc2d", "../../Atmos/FALC_03x82.atmos")
    cases += [("2d", prog2, "1 ray", ("ANGLE_SET", "SET_VERTICAL")),
              ("2d", prog2, "6 rays", ("ANGLE_SET", "SET_GL_3X2"))]
failed = nchecks = 0
for geo, prog, label, angles in cases:
    tag = f"{geo}_{label.replace(' ', '')}"
    J_av, _ = run(f"aver_{tag}", prog, angles, "PRD_ANGLE_AVER")
    J_hy, Jgas = run(f"hyb_{tag}", prog, angles, "PRD_ANGLE_APPROX")
    checks = [(f"{geo}, {label}: hybrid Jgas vs J", maxrel(Jgas, J_hy)),
              (f"{geo}, {label}: hybrid J vs angle-averaged J", maxrel(J_hy, J_av))]
    for what, err in checks:
        ok = len(J_hy) == len(J_av) == len(Jgas) and err <= TOL
        failed += not ok
        nchecks += 1
        if a.verbose or not ok:
            print(f"{'PASS' if ok else 'FAIL'}  {what}: {err:.2e} (tol {TOL:g})")

print(f"prd_hybrid: {nchecks} checks, {failed} failed")
if not a.work:
    shutil.rmtree(work)
sys.exit(1 if failed else 0)
