#!/usr/bin/env python3
"""End-to-end test of ANGLE_SET = SET_EDDINGTON in rhsc2d, B-57.

usage: eddington_2d.py --rhf1d P --rhsc2d P [--solveray2d P] [--chartflow2d P]
                       [--rhsc3d P] [--src DIR] [--work DIR] [-v]

In 2-D, SET_EDDINGTON (the old EDDINGTON keyword, lost in upstream df4c38f)
solves every column on its own with the Feautrier scheme: one ray at
mu = 1/sqrt(3), weight 1, no horizontal transfer, no stencils. Ca II (with
PRD) on FAL-C, three columns:

 1. The run converges, and every column of a horizontally uniform atmosphere
    gives the same J and emergent I, exactly.
 2. It agrees with rhf1d, ANGLE_SET = SET_EDDINGTON, on the same atmosphere
    (one ray at mu = 1/sqrt(3) as well) to the difference of the two formal
    solvers (Feautrier against DELO-Bezier on 82 depths, not round-off).
    Measured (2-D to ITER_LIMIT 1e-4, 1-D to 1e-6): emergent I differs by
    0.4% (median over the 287 wavelengths) and by 1.8% at most from 300 nm
    on (in the Ca II K core); below 150 nm, where the optical depth steps of
    the grid are large and I is tiny, by 6% (median); at 125 nm, a deep
    minimum of I, by 50% and sensitive to the convergence, so no worst case
    is bounded there. J relative to its peak over depth: 0.3% at most from
    300 nm on, 11% at most below 150 nm. For scale: a wrong mu (1 for 1/sqrt(3))
    moves I by 12% (median, from 300 nm). The bounds have a
    factor of 1.5 to 3 headroom.
 3. Column independence: in a run with different columns (ne scaled
    by 1, 1.3, 0.7), each column equals the run with that column's atmosphere
    in every column, exactly (10 iterations of SOLVER = NG without acceleration, which
    couples the points, as does the default SOLVER = ANDERSON).
 4. solveray of a Eddington run uses the short characteristics of its own
    ray (finite I, the same in every column).

    LIMIT_MEMORY = TRUE (J in a buffer, not in spectrum.J) gives the same
    J.dat and spectrum.out as FALSE, to the bit (it crashed when the Eddington
    branch wrote spectrum.J directly).
 4b. solveray is not refused when a magnetic field and SET_EDDINGTON are in its
    keyword.input: it solves a polarized ray with short characteristics.

Not available with SET_EDDINGTON (needs angle-dependent opacity, for which
there is no solver), each stops with a message and status 70:
 5. at readInput: 3-D geometry; in checkEddington() (rhsc2d, gcheck, after
    Background()): BACKGROUND_POLARIZATION, a magnetic field, PRD_ANGLE_DEP
    and PRD_ANGLE_APPROX with a PRD line, a moving atmosphere with active
    lines or with background lines; Formal() is the backstop.
    chartflow stops too. Before B-57 all of these gave "Did not recognize
    angleset" (or nothing: input.Eddington was never set).

Python 3 standard library only. Exit status 0 = pass.
"""
import math
import os
import statistics
import struct
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker, maxrel, same_file
from rhtest.cli import parse
from rhtest.keywords import H6_CAII, edit_file
from rhtest.output import read_bytes, read_doubles, read_solveray_intensity, read_spectrum
from rhtest.run import EXIT_FATAL, copy_dir, execute, make_field, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d", "rhsc2d"),
          optional=("solveray2d", "chartflow2d", "rhsc3d"), verbose=True,
          src_help="RH source tree (Atoms/, Atmos/, Molecules/, rh*/run/)")
src = a.src
work = WorkDir(src, a.work)
ck = Checker(a.verbose)

CAII = H6_CAII
BACKGROUND_ONLY = H6_CAII[:1]


def run(name, program, rundir, kws, atoms=CAII, files=(), atmos_bytes=None,
        written=None, timeout=300):
    """Run program in a fresh copy of the run files of rundir with the
    keywords kws, the files `written` ({name: text}) added. Returns the
    Result."""
    base = [("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
            ("N_THREADS", "1"), ("PRINT_CPU", "FALSE"), ("STARTING_J", "NEW_J"),
            ("KURUCZ_DATA", "none"), ("SOLVE_NE", "NONE"), ("PRD_N_MAX_ITER", "2"),
            ("PRD_ANGLE_DEP", "PRD_ANGLE_AVER"), ("LAMBDA_REF", "500.0")]
    d = prepare(work, name, rundir, base + list(kws), atoms, copy=(),
                required=("molecules.input", "ray.input") + tuple(files))
    for f, text in (written or {}).items():
        open(os.path.join(d, f), "w").write(text)
    if atmos_bytes is not None:
        open(os.path.join(d, "atmos.xdr"), "wb").write(atmos_bytes)
    return execute(program, d, timeout=timeout)


def must_run(*args, **kw):
    return run(*args, **kw).require(ok=None, bad_text="ERROR")    # any status, FP trap (B-33)


def read_I(d, Nx):
    """Wavelengths and the emergent I per wavelength and column (one ray)."""
    return read_spectrum(os.path.join(d, "spectrum.out"), Nx)


# --- 2-D atmosphere files: the layout read by rhsc2d/readatmos_xdr.c
#     (see tests/make_tiled_2d_atmos.py): Nx, Nz, NH, three boundary
#     codes, dx[Nx], z[Nz], then T, ne, vturb, vx, vz (each [Nz][Nx]) and
#     nH[NH][Nz][Nx]; big-endian

def read_atmos2d(path):
    b = read_bytes(path)
    Nx, Nz, NH = struct.unpack(">3i", b[:12])
    N = Nx * Nz
    pos = 24 + 8 * Nx + 8 * Nz
    head = b[:pos]
    fields = [list(struct.unpack(f">{N}d", b[pos + 8 * N * i:pos + 8 * N * (i + 1)]))
              for i in range(5)]
    tail = b[pos + 40 * N:]
    if len(tail) != 8 * NH * N:
        raise ValueError(f"{path}: unexpected layout of the 2-D atmosphere")
    return Nx, Nz, head, fields, tail


def write_atmos2d(Nx, Nz, head, fields, tail):
    return head + b"".join(struct.pack(f">{Nx * Nz}d", *f) for f in fields) + tail


Nx, Nz, head, fields, tail = read_atmos2d(os.path.join(src, "Atmos", "FALC_03x82.atmos"))
rmean = lambda xs: statistics.median(xs)

# ===================================================================
# 1, 2: uniform atmosphere, convergence, columns, rhf1d

EDD = [("ANGLE_SET", "SET_EDDINGTON"), ("ITER_LIMIT", "1.0E-4"), ("N_MAX_ITER", "300")]
r2 = must_run("e2d", a.rhsc2d, "rhsc2d", EDD + [
    ("ATMOS_FILE", "../../Atmos/FALC_03x82.atmos")])
d2 = r2.dir
ck.check(r2.returncode == 0 and "Converged in" in r2.stdout, "2-D SET_EDDINGTON run converges",
         f"(status {r2.returncode})")
r1 = must_run("e1d", a.rhf1d, "rhf1d", [
    ("ANGLE_SET", "SET_EDDINGTON"), ("NRAYS", "1"), ("ITER_LIMIT", "1.0E-6"),
    ("N_MAX_ITER", "300"), ("ATMOS_FILE", "../../Atmos/FALC_82.atmos")])
d1 = r1.dir
ck.check(r1.returncode == 0 and "Converged in" in r1.stdout, "1-D SET_EDDINGTON run converges",
         f"(status {r1.returncode})")

lam1, I1 = read_I(d1, 1)
lam2, I2 = read_I(d2, Nx)
J1 = read_doubles(os.path.join(d1, "J.dat"))
J2 = read_doubles(os.path.join(d2, "J.dat"))
ns = len(lam1)
ck.check(len(lam2) == ns and len(J1) == ns * Nz and len(J2) == ns * Nz * Nx,
         "1-D and 2-D runs have the same wavelengths and shapes")
if ck.failed:
    ck.finish("eddington_2d", work)

# --- every column equal
spreadI = max((max(row) - min(row)) / max(abs(x) for x in row) for row in I2)
spreadJ = max(maxrel(J2[la * Nz * Nx:(la + 1) * Nz * Nx:Nx],
                     J2[la * Nz * Nx + c:(la + 1) * Nz * Nx:Nx])
              for la in range(ns) for c in range(1, Nx))
ck.check(spreadI == 0 and spreadJ == 0,
         "uniform atmosphere: every column gives the same I and J",
         f"(I {spreadI:.1e}, J {spreadJ:.1e})")

# --- against rhf1d
Ie = [abs(row2[0] - row1[0]) / abs(row1[0]) for row1, row2 in zip(I1, I2)]
Jpk = []
for la in range(ns):
    j1 = J1[la * Nz:(la + 1) * Nz]
    j2 = J2[la * Nz * Nx:(la + 1) * Nz * Nx:Nx]
    Jpk.append(max(abs(u - v) for u, v in zip(j1, j2)) / max(abs(u) for u in j1))
vis = [la for la in range(ns) if lam1[la] >= 300.0]
euv = [la for la in range(ns) if lam1[la] < 150.0]
stats = [("I, median over wavelengths", rmean(Ie), 1.0e-2),
         ("I, worst from 300 nm", max(Ie[la] for la in vis), 3.0e-2),
         ("I, median below 150 nm", rmean([Ie[la] for la in euv]), 0.15),
         ("J (of its peak), worst from 300 nm", max(Jpk[la] for la in vis), 1.0e-2),
         ("J (of its peak), worst below 150 nm", max(Jpk[la] for la in euv), 0.15)]
for what, val, tol in stats:
    ck.check(val <= tol, f"2-D Feautrier against 1-D SET_EDDINGTON: {what}",
             f"({val:.2e}, bound {tol:g})")

# ===================================================================
# 3: column independence

ITERS = [("ANGLE_SET", "SET_EDDINGTON"), ("SOLVER", "NG"), ("N_MAX_ITER", "10"), ("ITER_LIMIT", "1.0E-30"),
         ("NG_ORDER", "0"), ("PRD_NG_ORDER", "0"), ("PRD_ITER_LIMIT", "0"),
         ("N_MAX_SCATTER", "2"), ("ATMOS_FILE", "atmos.xdr")]
scale = [1.0, 1.3, 0.7]


def atmos_ne(factors):
    """The atmosphere with ne of column l multiplied by factors[l]. (Not T:
    the wavelength grids of the lines follow the largest Doppler width of
    the atmosphere, so every run would have its own grid.)"""
    f = [list(x) for x in fields]
    f[1] = [fields[1][i] * factors[i % Nx] for i in range(Nx * Nz)]
    return write_atmos2d(Nx, Nz, head, f, tail)


def run_cols(name, factors):
    d = must_run(name, a.rhsc2d, "rhsc2d", ITERS, atmos_bytes=atmos_ne(factors)).dir
    return read_doubles(os.path.join(d, "J.dat")), read_I(d, Nx)[1]


Jv, Iv = run_cols("cols_varied", scale)
worst = 0.0
for c in range(Nx):
    Ju, Iu = run_cols(f"cols_uniform{c}", [scale[c]] * Nx)
    for la in range(ns):
        worst = max(worst,
                    maxrel(Jv[la * Nz * Nx + c:(la + 1) * Nz * Nx:Nx],
                           Ju[la * Nz * Nx + c:(la + 1) * Nz * Nx:Nx]),
                    abs(Iv[la][c] - Iu[la][c]) / abs(Iu[la][c]))
ck.check(worst == 0, "columns are independent: varied atmosphere against "
         "the uniform atmosphere of each column", f"(worst relative difference {worst:.1e})")

# ===================================================================
# 4: solveray after an Eddington run

# --- LIMIT_MEMORY = TRUE gives the same J and I (CRD Ca II: PRD cannot
#     limit memory)
CAII_CRD = [CAII[0], (CAII[1][0].replace("CaII.atom", "CaII_CRD.atom"),) + CAII[1][1:]]
LM = [("ANGLE_SET", "SET_EDDINGTON"), ("SOLVER", "NG"), ("N_MAX_ITER", "3"),
      ("NG_ORDER", "0"), ("PRD_NG_ORDER", "0"), ("PRD_ITER_LIMIT", "0"),
      ("ATMOS_FILE", "../../Atmos/FALC_03x82.atmos")]
dm = must_run("lm_false", a.rhsc2d, "rhsc2d", LM + [("LIMIT_MEMORY", "FALSE")],
              atoms=CAII_CRD).dir
rt = run("lm_true", a.rhsc2d, "rhsc2d", LM + [("LIMIT_MEMORY", "TRUE")], atoms=CAII_CRD)
dt = rt.dir
same = rt.ran and all(
    os.path.exists(os.path.join(dt, f)) and same_file(dm, dt, f)
    for f in ("J.dat", "spectrum.out"))
ck.check(same, "LIMIT_MEMORY = TRUE gives the same J.dat and spectrum.out as FALSE, to the bit",
         f"(status {rt.returncode})")

# --- a magnetic field with SET_EDDINGTON in solveray's keyword.input: not
#     refused (solveray solves its ray with short characteristics)
bfile = work.join("field_2d.B")
make_field(src, "2d", os.path.join(src, "Atmos", "FALC_03x82.atmos"), bfile)
if a.solveray2d:
    db = must_run("ray_B_main", a.rhsc2d, "rhsc2d",
                  [("ANGLE_SET", "SET_VERTICAL"), ("SOLVER", "NG"),
                   ("N_MAX_ITER", "2"), ("NG_ORDER", "0"),
                   ("PRD_NG_ORDER", "0"), ("PRD_ITER_LIMIT", "0"),
                   ("ATMOS_FILE", "../../Atmos/FALC_03x82.atmos"),
                   ("STOKES_INPUT", bfile), ("STOKES_MODE", "FIELD_FREE")]).dir
    edit_file(os.path.join(db, "keyword.input"), [("ANGLE_SET", "SET_EDDINGTON")])
    open(os.path.join(db, "ray.input"), "w").write("0.5 0.7\n0\n")
    r = execute(a.solveray2d, db, timeout=300, log="solveray.log")
    ck.check(r.returncode == 0 and "ERROR" not in r.stdout and
             any(x.startswith("spectrum_") for x in os.listdir(db)),
             "solveray with SET_EDDINGTON and a magnetic field in keyword.input runs",
             f"(status {r.returncode})")

if a.solveray2d:
    d = copy_dir(work, "ray", d2)
    open(os.path.join(d, "ray.input"), "w").write("0.0 1.0\n0\n")
    r = execute(a.solveray2d, d, timeout=300, log=None)
    files = [x for x in os.listdir(d) if x.startswith("spectrum_")]
    ok = r.returncode == 0 and "ERROR" not in r.stdout and len(files) == 1
    if ok:
        M = ns * Nx
        try:
            vals = read_solveray_intensity(os.path.join(d, files[0]), M)
        except ValueError as exc:          # wrong size: a failed check, not the end
            print(exc)
            ok = False
        else:
            rows = [vals[i * Nx:(i + 1) * Nx] for i in range(ns)]
            ok = (all(math.isfinite(v) for v in vals)
                  and all(max(row) - min(row) <= 1e-12 * max(abs(x) for x in row)
                          for row in rows))
    ck.check(ok, "solveray after a SET_EDDINGTON run: finite, the same in every column",
             f"(status {r.returncode})")

# ===================================================================
# 5: what is refused

moving = [list(x) for x in fields]
moving[4] = [5.0] * (Nx * Nz)                  # vz in km/s: a moving atmosphere
moving_bytes = write_atmos2d(Nx, Nz, head, moving, tail)

# the Kurucz line list of the background-lines case: Fe I 630.25 nm
FE_LIST = {"kurucz.input": "../../Atoms/Kurucz/gf06302_1564.kur\n"}

EDDQ = [("ANGLE_SET", "SET_EDDINGTON"), ("N_MAX_ITER", "2")]
refusals = [
    ("background polarization", a.rhsc2d, "rhsc2d", EDDQ +
     [("BACKGROUND_POLARIZATION", "TRUE")], {}, "cannot treat background polarization"),
    ("magnetic field", a.rhsc2d, "rhsc2d", EDDQ +
     [("STOKES_INPUT", bfile),
      ("STOKES_MODE", "FULL_STOKES")],
     {}, "cannot treat a magnetic field"),
    ("PRD_ANGLE_APPROX with a PRD line", a.rhsc2d, "rhsc2d", EDDQ +
     [("PRD_ANGLE_DEP", "PRD_ANGLE_APPROX")], {},
     "cannot treat angle-dependent PRD"),
    ("PRD_ANGLE_DEP with a PRD line", a.rhsc2d, "rhsc2d", EDDQ +
     [("PRD_ANGLE_DEP", "PRD_ANGLE_DEP")], {}, "cannot treat angle-dependent PRD"),
    ("moving atmosphere, active lines", a.rhsc2d, "rhsc2d", EDDQ +
     [("ATMOS_FILE", "atmos.xdr")], {"atmos_bytes": moving_bytes},
     "cannot treat a moving atmosphere"),
    ("moving atmosphere, background lines", a.rhsc2d, "rhsc2d", EDDQ +
     [("ATMOS_FILE", "atmos.xdr"), ("KURUCZ_DATA", "kurucz.input"),
      ("LAMBDA_REF", "630.42")],
     {"atmos_bytes": moving_bytes, "atoms": BACKGROUND_ONLY, "written": FE_LIST},
     "cannot treat a moving atmosphere with background lines"),
]
if a.rhsc3d:
    refusals.append(("3-D geometry", a.rhsc3d, "rhsc3d", [("ANGLE_SET", "SET_EDDINGTON"),
                     ("N_MAX_ITER", "2")], {}, "not available in 3-D"))
for what, prog, rundir, kws, extra, msg in refusals:
    r = run("refuse_" + what.replace(" ", "_").replace(",", ""), prog, rundir, kws, **extra)
    rc, out = r.returncode, r.stdout
    flat = " ".join(out.split())
    ck.check(rc == EXIT_FATAL and " ".join(msg.split()) in flat and "Did not recognize" not in out,
             f"SET_EDDINGTON refused: {what}", f"(status {rc}, message '{msg}')")

if a.chartflow2d:
    r = run("refuse_chartflow", a.chartflow2d, "rhsc2d", EDDQ, files=("flow.input",))
    rc = r.returncode
    ck.check(rc == EXIT_FATAL and "no short characteristics to analyze" in " ".join(r.stdout.split()),
             "SET_EDDINGTON refused: chartflow", f"(status {rc})")

ck.finish("eddington_2d", work)
