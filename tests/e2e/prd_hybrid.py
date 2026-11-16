#!/usr/bin/env python3
"""End-to-end test of hybrid PRD (PRD_ANGLE_APPROX), B-31.

usage: prd_hybrid.py --rhf1d PATH [--rhsc2d PATH] [--rhsc3d PATH] [--src DIR] [--work DIR] [-v]

In a static atmosphere the gas frame is the observer's frame, so the
gas-frame mean intensity Jgas of the hybrid approximation equals J, and
hybrid PRD gives the same answer as angle-averaged PRD. Ca II with PRD,
rhf1d on FAL-C with one ray and with three rays, (with --rhsc2d) rhsc2d on
three FAL-C columns with one ray and with six, three iterations each, and
(with --rhsc3d) rhsc3d on 4x4 FAL-C columns with six rays, one iteration
(the identity holds after any number of iterations, and 3-D is slow):

 1. Hybrid run: Jgas = J at every wavelength and depth.
 2. J of the hybrid run equals J of the angle-averaged run.

1 holds to 1e-12. In 2 the two runs sum the same intensities in a
different order, and the iterations of the PRD and statistical
equilibrium (three; one in 3-D, where the difference is 0) amplify the rounding (2e-12 on macOS/arm64, where the
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
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker, maxrel, same_file
from rhtest.cli import parse
from rhtest.keywords import H6_CAII
from rhtest.output import read_doubles
from rhtest.run import execute, prepare
from rhtest.workdir import WorkDir

TOL_JGAS = 1e-12   # Jgas vs J of the same run
TOL_RUNS = 1e-9    # J of the hybrid run vs J of the angle-averaged run

a = parse(required=("rhf1d",), optional=("rhsc2d", "rhsc3d"), verbose=True)
rhf1d, rhsc2d, rhsc3d, src = a.rhf1d, a.rhsc2d, a.rhsc3d, a.src
work = WorkDir(src, a.work)
ck = Checker(a.verbose)


def run(name, prog, angles, mode, extra=(), iters=3):
    """prog: (program, run dir in the source tree, atmosphere); angles:
    keyword and value; extra: more keywords. Returns the run directory,
    J and, for the hybrid mode, Jgas."""
    program, rundir, atmos = prog
    d = prepare(work, name, rundir,
                [("ATMOS_FILE", atmos),
                 ("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
                 angles, ("N_MAX_ITER", str(iters)),
                 ("PRD_N_MAX_ITER", "2"), ("PRD_ANGLE_DEP", mode),
                 ("N_THREADS", "1"), ("PRINT_CPU", "FALSE"),
                 ("STARTING_J", "NEW_J"), ("KURUCZ_DATA", "none"),
                 ("SOLVE_NE", "NONE")] + list(extra), H6_CAII)
    # Three iterations (or one) do not converge: exit status 1 is expected
    execute(program, d).require(bad_text="ERROR")
    J = read_doubles(os.path.join(d, "J.dat"))
    Jgas = read_doubles(os.path.join(d, "Jgas.dat")) if mode == "PRD_ANGLE_APPROX" else None
    return d, J, Jgas


cases = [("1d", (rhf1d, "rhf1d", "../../Atmos/FALC_82.atmos"), "1 ray", ("NRAYS", "1")),
         ("1d", (rhf1d, "rhf1d", "../../Atmos/FALC_82.atmos"), "3 rays", ("NRAYS", "3"))]
if rhsc2d:
    prog2 = (rhsc2d, "rhsc2d", "../../Atmos/FALC_03x82.atmos")
    cases += [("2d", prog2, "1 ray", ("ANGLE_SET", "SET_VERTICAL")),
              ("2d", prog2, "6 rays", ("ANGLE_SET", "SET_GL_3X2"))]
if rhsc3d:
    prog3 = (rhsc3d, "rhsc3d", "../../Atmos/FALC_04x04x82.atmos")
    cases += [("3d", prog3, "4 rays", ("ANGLE_SET", "SET_GL_1X1"), 1)]
for geo, prog, label, angles, *rest in cases:
    iters = rest[0] if rest else 3
    tag = f"{geo}_{label.replace(' ', '')}"
    _, J_av, _ = run(f"aver_{tag}", prog, angles, "PRD_ANGLE_AVER", iters=iters)
    _, J_hy, Jgas = run(f"hyb_{tag}", prog, angles, "PRD_ANGLE_APPROX", iters=iters)
    checks = [(f"{geo}, {label}: hybrid Jgas vs J", maxrel(Jgas, J_hy), TOL_JGAS),
              (f"{geo}, {label}: hybrid J vs angle-averaged J", maxrel(J_hy, J_av), TOL_RUNS)]
    for what, err, tol in checks:
        ok = len(J_hy) == len(J_av) == len(Jgas) and err <= tol
        ck.check(ok, what, f"{err:.2e} (tol {tol:g})")

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
d_tab, J_tab, Jgas_tab = run("hyb_1d_v10_table", prog1, ("NRAYS", "3"), "PRD_ANGLE_APPROX")
d_mem, J_mem, Jgas_mem = run("hyb_1d_v10_limitmem", prog1, ("NRAYS", "3"), "PRD_ANGLE_APPROX",
                      [("PRDH_LIMIT_MEM", "TRUE")])
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
    ck.check(ok, f"1d moving, {label}: Jgas vs J at the {len(rest)} depths at rest",
             f"{err:.2e} (tol {TOL_JGAS:g}), worst at wavelength, depth {where} of {Nspect}")
same = J_tab == J_mem and Jgas_tab == Jgas_mem and same_file(d_tab, d_mem, "pops.CA.out")
ck.check(same, "1d moving: PRDH_LIMIT_MEM run equals the table run to the bit",
         f"(J {maxrel(J_tab, J_mem):.2e}, Jgas {maxrel(Jgas_tab, Jgas_mem):.2e})")

ck.finish("prd_hybrid", work)
