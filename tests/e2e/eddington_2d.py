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
import argparse
import math
import os
import shutil
import statistics
import struct
import subprocess
import sys
import tempfile

here = os.path.dirname(os.path.abspath(__file__))
p = argparse.ArgumentParser()
p.add_argument("--rhf1d", required=True)
p.add_argument("--rhsc2d", required=True)
p.add_argument("--solveray2d")
p.add_argument("--chartflow2d")
p.add_argument("--rhsc3d")
p.add_argument("--src", default=os.path.dirname(os.path.dirname(here)),
               help="RH source tree (Atoms/, Atmos/, Molecules/, rh*/run/)")
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


CAII = ("2\n  ../../Atoms/H_6.atom  PASSIVE  LTE_POPULATIONS\n"
        "  ../../Atoms/CaII.atom  ACTIVE  ZERO_RADIATION  pops.CA.out\n")
BACKGROUND_ONLY = "1\n  ../../Atoms/H_6.atom  PASSIVE  LTE_POPULATIONS\n"


def run(name, program, rundir, kws, atoms=CAII, files=(), atmos_bytes=None,
        written=None, timeout=300):
    """Run program in a fresh copy of the run files of rundir with the
    keywords kws, the files `written` ({name: text}) added. Returns (status,
    output, directory)."""
    d = os.path.join(work, "runs", name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    kw = open(os.path.join(src, rundir, "run", "keyword.input")).read()
    base = [("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
            ("N_THREADS", "1"), ("PRINT_CPU", "FALSE"), ("STARTING_J", "NEW_J"),
            ("KURUCZ_DATA", "none"), ("SOLVE_NE", "NONE"), ("PRD_N_MAX_ITER", "2"),
            ("PRD_ANGLE_DEP", "PRD_ANGLE_AVER"), ("LAMBDA_REF", "500.0")]
    for k, v in base + list(kws):
        kw = setkw(kw, k, v)
    open(os.path.join(d, "keyword.input"), "w").write(kw)
    open(os.path.join(d, "atoms.input"), "w").write(atoms)
    for f in ("molecules.input", "ray.input") + tuple(files):
        shutil.copy(os.path.join(src, rundir, "run", f), d)
    for f, text in (written or {}).items():
        open(os.path.join(d, f), "w").write(text)
    if atmos_bytes is not None:
        open(os.path.join(d, "atmos.xdr"), "wb").write(atmos_bytes)
    try:
        r = subprocess.run([os.path.abspath(program)], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        sys.exit(f"{program} timed out in {d}")
    open(os.path.join(d, "out.log"), "w").write(r.stdout)
    return r.returncode, r.stdout, d


def must_run(*args, **kw):
    rc, out, d = run(*args, **kw)
    if "ERROR" in out:
        tail = "\n".join(out.splitlines()[-40:])           # FP trap backtrace (B-33)
        sys.exit(f"{args[1]} failed in {d} (status {rc}); end of out.log:\n{tail}")
    return rc, out, d


def read_doubles(path):
    raw = open(path, "rb").read()
    return struct.unpack(f"={len(raw) // 8}d", raw)


def read_spectrum(d, Nx):
    """Wavelengths and the emergent I per wavelength and column (one ray)."""
    raw = open(os.path.join(d, "spectrum.out"), "rb").read()
    ns = struct.unpack(">i", raw[:4])[0]
    lam = struct.unpack(f">{ns}d", raw[4:4 + 8 * ns])
    pos = 4 + 8 * ns
    I = struct.unpack(f">{ns * Nx}d", raw[pos:pos + 8 * ns * Nx])
    return lam, [I[i * Nx:(i + 1) * Nx] for i in range(ns)]


# --- 2-D atmosphere files: the layout read by rhsc2d/readatmos_xdr.c
#     (see tests/make_tiled_2d_atmos.py): Nx, Nz, NH, three boundary
#     codes, dx[Nx], z[Nz], then T, ne, vturb, vx, vz (each [Nz][Nx]) and
#     nH[NH][Nz][Nx]; big-endian

def read_atmos2d(path):
    b = open(path, "rb").read()
    Nx, Nz, NH = struct.unpack(">3i", b[:12])
    N = Nx * Nz
    pos = 24 + 8 * Nx + 8 * Nz
    head = b[:pos]
    fields = [list(struct.unpack(f">{N}d", b[pos + 8 * N * i:pos + 8 * N * (i + 1)]))
              for i in range(5)]
    tail = b[pos + 40 * N:]
    assert len(tail) == 8 * NH * N, "unexpected layout of the 2-D atmosphere"
    return Nx, Nz, head, fields, tail


def write_atmos2d(Nx, Nz, head, fields, tail):
    return head + b"".join(struct.pack(f">{Nx * Nz}d", *f) for f in fields) + tail


def maxrel(x, y):
    return max(abs(u - v) / max(abs(u), abs(v), 1e-300) for u, v in zip(x, y))


Nx, Nz, head, fields, tail = read_atmos2d(os.path.join(src, "Atmos", "FALC_03x82.atmos"))
rmean = lambda xs: statistics.median(xs)

# ===================================================================
# 1, 2: uniform atmosphere, convergence, columns, rhf1d

EDD = [("ANGLE_SET", "SET_EDDINGTON"), ("ITER_LIMIT", "1.0E-4"), ("N_MAX_ITER", "300")]
rc2, out2, d2 = must_run("e2d", a.rhsc2d, "rhsc2d", EDD + [
    ("ATMOS_FILE", "../../Atmos/FALC_03x82.atmos")])
check(rc2 == 0 and "Converged in" in out2, "2-D SET_EDDINGTON run converges",
      f"(status {rc2})")
rc1, out1, d1 = must_run("e1d", a.rhf1d, "rhf1d", [
    ("ANGLE_SET", "SET_EDDINGTON"), ("NRAYS", "1"), ("ITER_LIMIT", "1.0E-6"),
    ("N_MAX_ITER", "300"), ("ATMOS_FILE", "../../Atmos/FALC_82.atmos")])
check(rc1 == 0 and "Converged in" in out1, "1-D SET_EDDINGTON run converges",
      f"(status {rc1})")

lam1, I1 = read_spectrum(d1, 1)
lam2, I2 = read_spectrum(d2, Nx)
J1 = read_doubles(os.path.join(d1, "J.dat"))
J2 = read_doubles(os.path.join(d2, "J.dat"))
ns = len(lam1)
check(len(lam2) == ns and len(J1) == ns * Nz and len(J2) == ns * Nz * Nx,
      "1-D and 2-D runs have the same wavelengths and shapes")
if failed:
    sys.exit(1)

# --- every column equal
spreadI = max((max(row) - min(row)) / max(abs(x) for x in row) for row in I2)
spreadJ = max(maxrel(J2[la * Nz * Nx:(la + 1) * Nz * Nx:Nx],
                     J2[la * Nz * Nx + c:(la + 1) * Nz * Nx:Nx])
              for la in range(ns) for c in range(1, Nx))
check(spreadI == 0 and spreadJ == 0,
      "uniform atmosphere: every column gives the same I and J",
      f"(I {spreadI:.1e}, J {spreadJ:.1e})")

# --- against rhf1d
Ie = [abs(r2[0] - r1[0]) / abs(r1[0]) for r1, r2 in zip(I1, I2)]
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
    check(val <= tol, f"2-D Feautrier against 1-D SET_EDDINGTON: {what}",
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
    rc, out, d = must_run(name, a.rhsc2d, "rhsc2d", ITERS,
                          atmos_bytes=atmos_ne(factors))
    return read_doubles(os.path.join(d, "J.dat")), read_spectrum(d, Nx)[1]


Jv, Iv = run_cols("cols_varied", scale)
worst = 0.0
for c in range(Nx):
    Ju, Iu = run_cols(f"cols_uniform{c}", [scale[c]] * Nx)
    for la in range(ns):
        worst = max(worst,
                    maxrel(Jv[la * Nz * Nx + c:(la + 1) * Nz * Nx:Nx],
                           Ju[la * Nz * Nx + c:(la + 1) * Nz * Nx:Nx]),
                    abs(Iv[la][c] - Iu[la][c]) / abs(Iu[la][c]))
check(worst == 0, "columns are independent: varied atmosphere against "
      "the uniform atmosphere of each column", f"(worst relative difference {worst:.1e})")

# ===================================================================
# 4: solveray after an Eddington run

# --- LIMIT_MEMORY = TRUE gives the same J and I (CRD Ca II: PRD cannot
#     limit memory)
CAII_CRD = CAII.replace("CaII.atom", "CaII_CRD.atom")
LM = [("ANGLE_SET", "SET_EDDINGTON"), ("SOLVER", "NG"), ("N_MAX_ITER", "3"),
      ("NG_ORDER", "0"), ("PRD_NG_ORDER", "0"), ("PRD_ITER_LIMIT", "0"),
      ("ATMOS_FILE", "../../Atmos/FALC_03x82.atmos")]
rcm, outm, dm = must_run("lm_false", a.rhsc2d, "rhsc2d", LM + [("LIMIT_MEMORY", "FALSE")],
                         atoms=CAII_CRD)
rct, outt, dt = run("lm_true", a.rhsc2d, "rhsc2d", LM + [("LIMIT_MEMORY", "TRUE")],
                    atoms=CAII_CRD)
same = rct in (0, 1) and all(
    os.path.exists(os.path.join(dt, f)) and
    open(os.path.join(dm, f), "rb").read() == open(os.path.join(dt, f), "rb").read()
    for f in ("J.dat", "spectrum.out"))
check(same, "LIMIT_MEMORY = TRUE gives the same J.dat and spectrum.out as FALSE, to the bit",
      f"(status {rct})")

# --- a magnetic field with SET_EDDINGTON in solveray's keyword.input: not
#     refused (solveray solves its ray with short characteristics)
bfile = os.path.join(work, "field_2d.B")
subprocess.run([sys.executable, "-I", os.path.join(src, "tests", "make_field.py"),
                "2d", os.path.join(src, "Atmos", "FALC_03x82.atmos"), bfile],
               check=True, stdout=subprocess.DEVNULL)
if a.solveray2d:
    rcb, outb, db = must_run("ray_B_main", a.rhsc2d, "rhsc2d",
                             [("ANGLE_SET", "SET_VERTICAL"), ("SOLVER", "NG"),
                              ("N_MAX_ITER", "2"), ("NG_ORDER", "0"),
                              ("PRD_NG_ORDER", "0"), ("PRD_ITER_LIMIT", "0"),
                              ("ATMOS_FILE", "../../Atmos/FALC_03x82.atmos"),
                              ("STOKES_INPUT", bfile), ("STOKES_MODE", "FIELD_FREE")])
    kw = open(os.path.join(db, "keyword.input")).read()
    open(os.path.join(db, "keyword.input"), "w").write(
        setkw(kw, "ANGLE_SET", "SET_EDDINGTON"))
    open(os.path.join(db, "ray.input"), "w").write("0.5 0.7\n0\n")
    r = subprocess.run([os.path.abspath(a.solveray2d)], cwd=db, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, timeout=300)
    check(r.returncode == 0 and "ERROR" not in r.stdout and
          any(x.startswith("spectrum_") for x in os.listdir(db)),
          "solveray with SET_EDDINGTON and a magnetic field in keyword.input runs",
          f"(status {r.returncode})")

if a.solveray2d:
    d = os.path.join(work, "runs", "ray")
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(d2, d)
    open(os.path.join(d, "ray.input"), "w").write("0.0 1.0\n0\n")
    r = subprocess.run([os.path.abspath(a.solveray2d)], cwd=d, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True, timeout=300)
    files = [x for x in os.listdir(d) if x.startswith("spectrum_")]
    ok = r.returncode == 0 and "ERROR" not in r.stdout and len(files) == 1
    if ok:
        raw = open(os.path.join(d, files[0]), "rb").read()
        M = ns * Nx
        # mux, muz (2 doubles), I, then the int Nspect = 0
        vals = struct.unpack(f">{M}d", raw[16:16 + 8 * M])
        rows = [vals[i * Nx:(i + 1) * Nx] for i in range(ns)]
        ok = (len(raw) == 16 + 8 * M + 4 and all(math.isfinite(v) for v in vals)
              and all(max(r) - min(r) <= 1e-12 * max(abs(x) for x in r) for r in rows))
    check(ok, "solveray after a SET_EDDINGTON run: finite, the same in every column",
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
    rc, out, d = run("refuse_" + what.replace(" ", "_").replace(",", ""), prog, rundir,
                     kws, **extra)
    flat = " ".join(out.split())
    check(rc == 70 and " ".join(msg.split()) in flat and "Did not recognize" not in out,
          f"SET_EDDINGTON refused: {what}", f"(status {rc}, message '{msg}')")

if a.chartflow2d:
    rc, out, d = run("refuse_chartflow", a.chartflow2d, "rhsc2d", EDDQ,
                     files=("flow.input",))
    check(rc == 70 and "no short characteristics to analyze" in " ".join(out.split()),
          "SET_EDDINGTON refused: chartflow", f"(status {rc})")

print(f"eddington_2d: {nchecks} checks, {failed} failed")
if failed:
    print(f"runs kept in {work}")
elif not a.work:
    shutil.rmtree(work)
sys.exit(1 if failed else 0)
