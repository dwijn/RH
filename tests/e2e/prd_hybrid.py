#!/usr/bin/env python3
"""End-to-end test of hybrid PRD (PRD_ANGLE_APPROX), B-31.

usage: prd_hybrid.py --rhf1d PATH [--rhsc2d PATH] [--rhsc3d PATH] [--src DIR] [--work DIR] [-v]

In a static atmosphere the gas frame is the observer's frame, so the
gas-frame mean intensity Jgas of the hybrid approximation equals J, and
hybrid PRD gives the same answer as angle-averaged PRD. Ca II with PRD,
three iterations, rhf1d on FAL-C with one ray and with three rays,
(with --rhsc2d) rhsc2d on three FAL-C columns with one ray and with six,
and (with --rhsc3d) rhsc3d on 4x4 FAL-C columns with one ray and with six:

 1. Hybrid run: Jgas = J at every wavelength and depth.
 2. J of the hybrid run equals J of the angle-averaged run.

1 holds to 1e-12. In 2 the two runs sum the same intensities in a
different order, and three iterations of the PRD and statistical
equilibrium amplify the rounding (2e-12 on macOS/arm64, where the
compiler fuses multiply-adds), so 2 is checked to 1e-9; the bug gave
differences of order 1. Before the fix of B-31, with one ray, Jgas grew by
J with every formal solution at the wavelengths without PRD line opacity
and stayed zero at those of the PRD lines; in 2-D and 3-D it was 2 J at
the wavelengths without PRD line opacity.

Moving atmosphere (P-17): rhf1d on FAL-C with an outflow of 10 km/s above
the lower chromosphere (tests/atmos/FALC_v10.atmos, 29 depths at rest),
three rays, hybrid PRD with the interpolation tables and with
PRDH_LIMIT_MEM = TRUE, which computes the same weights when needed:

 3. At the depths at rest the gas frame is the observer's frame:
    Jgas = J at every wavelength, to 1e-12. Before P-17 the table gave
    2 J at the last wavelength of the grid and PRDH_LIMIT_MEM gave Jgas
    without the intensity of the first wavelength there.
 4. The PRDH_LIMIT_MEM run equals the table run to the bit (J, Jgas,
    populations).

Exit status 0 = pass.
"""
import argparse
import os
import shutil
import struct
import subprocess
import sys
import tempfile

TOL_JGAS = 1e-12   # Jgas vs J of the same run
TOL_RUNS = 1e-9    # J of the hybrid run vs J of the angle-averaged run

p = argparse.ArgumentParser()
p.add_argument("--rhf1d", required=True)
p.add_argument("--rhsc2d", help="also test rhsc2d")
p.add_argument("--rhsc3d", help="also test rhsc3d (B-7)")
p.add_argument("--src", default=os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))),
    help="RH source tree (Atoms/, Atmos/, Molecules/, rhf1d/run/)")
p.add_argument("--work", help="keep the runs here (default: a temporary dir)")
p.add_argument("-v", "--verbose", action="store_true")
a = p.parse_args()
rhf1d = os.path.abspath(a.rhf1d)
rhsc2d = os.path.abspath(a.rhsc2d) if a.rhsc2d else None
rhsc3d = os.path.abspath(a.rhsc3d) if a.rhsc3d else None
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


def run(name, prog, angles, mode, extra=()):
    """prog: (program, run dir in the source tree, atmosphere); angles:
    keyword and value; extra: more keywords. Returns J and, for the
    hybrid mode, Jgas."""
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
                 ("SOLVE_NE", "NONE")] + list(extra):
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
        tail = "\n".join(r.stdout.splitlines()[-40:])   # FP trap backtrace (B-33)
        sys.exit(f"{program} failed in {d} (status {r.returncode}); end of out.log:\n{tail}")
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
if rhsc3d:
    prog3 = (rhsc3d, "rhsc3d", "../../Atmos/FALC_04x04x82.atmos")
    cases += [("3d", prog3, "1 ray", ("ANGLE_SET", "SET_VERTICAL")),
              ("3d", prog3, "6 rays", ("ANGLE_SET", "SET_GL_3X2"))]
failed = nchecks = 0
for geo, prog, label, angles in cases:
    tag = f"{geo}_{label.replace(' ', '')}"
    J_av, _ = run(f"aver_{tag}", prog, angles, "PRD_ANGLE_AVER")
    J_hy, Jgas = run(f"hyb_{tag}", prog, angles, "PRD_ANGLE_APPROX")
    checks = [(f"{geo}, {label}: hybrid Jgas vs J", maxrel(Jgas, J_hy), TOL_JGAS),
              (f"{geo}, {label}: hybrid J vs angle-averaged J", maxrel(J_hy, J_av), TOL_RUNS)]
    for what, err, tol in checks:
        ok = len(J_hy) == len(J_av) == len(Jgas) and err <= tol
        failed += not ok
        nchecks += 1
        if a.verbose or not ok:
            print(f"{'PASS' if ok else 'FAIL'}  {what}: {err:.2e} (tol {tol:g})")

# --- Moving atmosphere: Jgas = J at rest, PRDH_LIMIT_MEM equals the table
moving = os.path.join(src, "tests", "atmos", "FALC_v10.atmos")
v = []
for l in open(moving):
    f = l.split()
    if len(f) == 5 and not l.lstrip().startswith("*"):
        v.append(float(f[3]))
rest = [k for k, x in enumerate(v) if x == 0.0]
Nspace = len(v)
prog1 = (rhf1d, "rhf1d", moving)
J_tab, Jgas_tab = run("hyb_1d_v10_table", prog1, ("NRAYS", "3"), "PRD_ANGLE_APPROX")
J_mem, Jgas_mem = run("hyb_1d_v10_limitmem", prog1, ("NRAYS", "3"), "PRD_ANGLE_APPROX",
                      [("PRDH_LIMIT_MEM", "TRUE")])
d_tab = os.path.join(work, "runs", "hyb_1d_v10_table")
d_mem = os.path.join(work, "runs", "hyb_1d_v10_limitmem")
for label, J, Jgas in (("table", J_tab, Jgas_tab), ("PRDH_LIMIT_MEM", J_mem, Jgas_mem)):
    Nspect = len(J) // Nspace
    err, where = 0.0, None
    for la in range(Nspect):
        for k in rest:
            i = la * Nspace + k
            e = abs(Jgas[i] - J[i]) / max(abs(J[i]), 1e-300)
            if e > err:
                err, where = e, (la, k)
    ok = len(rest) > 0 and len(J) == Nspect * Nspace and err <= TOL_JGAS
    failed += not ok
    nchecks += 1
    if a.verbose or not ok:
        print(f"{'PASS' if ok else 'FAIL'}  1d moving, {label}: Jgas vs J at the {len(rest)} "
              f"depths at rest: {err:.2e} (tol {TOL_JGAS:g}), worst at wavelength, depth {where} "
              f"of {Nspect}")
same = (J_tab == J_mem and Jgas_tab == Jgas_mem and
        open(os.path.join(d_tab, "pops.CA.out"), "rb").read() ==
        open(os.path.join(d_mem, "pops.CA.out"), "rb").read())
failed += not same
nchecks += 1
if a.verbose or not same:
    print(f"{'PASS' if same else 'FAIL'}  1d moving: PRDH_LIMIT_MEM run equals the table run "
          f"to the bit (J {maxrel(J_tab, J_mem):.2e}, Jgas {maxrel(Jgas_tab, Jgas_mem):.2e})")

print(f"prd_hybrid: {nchecks} checks, {failed} failed")
if not a.work:
    shutil.rmtree(work)
sys.exit(1 if failed else 0)
