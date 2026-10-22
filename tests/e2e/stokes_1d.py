#!/usr/bin/env python3
"""End-to-end Stokes tests of rhf1d with a fictional atom in LTE.

usage: stokes_1d.py --rhf1d PATH [--src DIR] [--work DIR] [-v]

The fictional atom (TI_fict.atom: two lines with anomalous Zeeman patterns,
g_eff = 4/3 and 7/6) is active with LTE populations and N_MAX_ITER = 0, so
rhf1d computes one formal solution of the full Stokes vector through the
whole program: atom and field input, Bproject, profiles, opacities and the
polarized solver. N_MAX_SCATTER = 1 makes that solution use the initial
J = 0 for background scattering in every run, so all runs solve the same
linear problem:

 1. No field: Q = U = V = 0 exactly.
 2. Mirror y -> -y, the plane of the rays (x-z): the field, an axial
    vector, becomes (-Bx, By, -Bz), i.e. gamma -> pi - gamma, chi -> pi - chi,
    and the Stokes vector (I, Q, -U, -V). Exact (parity), with and without
    magneto-optical effects, any depth dependence, every ray.
 3. Without magneto-optical effects, B -> -B (gamma -> pi - gamma,
    chi -> chi + pi): (I, Q, U, -V), exact. With them it is not a symmetry:
    Faraday rotation is odd in B.
 4. Vertical ray, field azimuth turned by d: (Q, U) turn by 2 d. Exact
    for the transfer equation and, since the polarized Bezier solver
    limits the derivatives of (S_Q, S_U), (eta_Q, eta_U) and (rho_Q, rho_U)
    as vectors (N-25), also for the solver, at every depth resolution.
    (Limited component by component, the deviation was 7e-3 of max |Q|,
    |U| on FAL-C and vanished only with refinement.)
 5. Weak field. To first order in B, for a field that is constant along
    the ray and LTE (RH's sign of V, see docs/unit-tests, section 8):
        V = -g_eff v_B cos(theta) dI/dv,
    v the wavelength in Doppler units. A uniform vertical velocity vz
    shifts the profiles by -mu vz (to the blue for an upflow), so
        V = -g_eff u_B (cos(theta) / mu) dI/dvz,  u_B = e lambda0 B / (4 pi m_e).
    dI/dvz comes from runs with vz = +-dv, +-2dv (fourth order); B = 1 G
    makes the O(B^2) terms negligible. V and dI/dvz reach the solver in
    different ways (an extra source term vs. shifted opacities), so they
    agree to discretization error: checked on a smooth grey atmosphere at
    40, 80, 160 and 320 depth points, where the difference must converge.

Tests 1-3 use FAL-C with a depth-dependent field, 4 FAL-C and the grey
atmosphere, 5 the grey atmosphere. Exit status 0 = pass.
"""
import argparse
import math
import os
import shutil
import struct
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

p = argparse.ArgumentParser()
p.add_argument("--rhf1d", required=True)
p.add_argument("--src", default=os.path.dirname(os.path.dirname(HERE)),
               help="RH source tree (Atoms/, Atmos/, Molecules/, rhf1d/run/)")
p.add_argument("--work", help="keep the runs here (default: a temporary dir)")
p.add_argument("-v", "--verbose", action="store_true")
a = p.parse_args()
rhf1d = os.path.abspath(a.rhf1d)
src = os.path.abspath(a.src)

work = os.path.abspath(a.work) if a.work else tempfile.mkdtemp(prefix="rh_e2e_")
os.makedirs(os.path.join(work, "runs"), exist_ok=True)
for d in ("Atoms", "Atmos", "Molecules"):          # run dirs sit two levels down
    link = os.path.join(work, d)
    if not os.path.lexists(link):
        os.symlink(os.path.join(src, d), link)

# --- Physical input ---------------------------------------------------------

E_CHARGE, M_ELECTRON = 1.602176634e-19, 9.1093837015e-31
LINES = [  # vacuum lambda0 [nm] from the level energies, g_eff from LS coupling
    (1e7 / 17000.0, None, (0.5, 0, 0.5), (0.5, 1, 0.5)),
    (1e7 / 17200.0, None, (0.5, 0, 0.5), (0.5, 1, 1.5))]


def lande(S, L, J):
    return 1.0 + (J * (J + 1) + S * (S + 1) - L * (L + 1)) / (2 * J * (J + 1))


def g_eff(lo, up):
    (Sl, Ll, Jl), (Su, Lu, Ju) = lo, up
    gl, gu = lande(Sl, Ll, Jl), lande(Su, Lu, Ju)
    return 0.5 * (gu + gl) + 0.25 * (gu - gl) * (Ju * (Ju + 1) - Jl * (Jl + 1))


LINES = [(lam, g_eff(lo, up)) for lam, _, lo, up in LINES]
assert abs(LINES[0][1] - 4 / 3) < 1e-15 and abs(LINES[1][1] - 7 / 6) < 1e-15


def gauss_legendre_01(n):
    """Nodes and weights on [0, 1] (Newton on P_n)."""
    out = []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-16:
                break
        out.append(0.5 * (1 - x))
    return sorted(out)


NRAYS = 3
MU = gauss_legendre_01(NRAYS)        # rhf1d's rays: muz, mux = sqrt(1 - muz^2), muy = 0


def write_atmos_text(path, logm, T, ne, v, vturb, ident):
    """1-D atmosphere on a column mass scale; no hydrogen populations (RH
    computes them in LTE, with ne from SOLVE_NE)."""
    with open(path, "w") as f:
        f.write(f"* {ident}\n  {ident}\n  Mass scale\n* log g [cm s^-2]\n  4.44\n"
                f"* Ndep\n  {len(logm)}\n*\n* lg m  T  ne [cm^-3]  V  Vturb [km/s]\n")
        for row in zip(logm, T, ne, v, vturb):
            f.write("  " + "  ".join(f"{x:.12E}" for x in row) + "\n")


def grey(n, vz=0.0):
    """Smooth grey atmosphere: T^4 = 3/4 Teff^4 (tau + 2/3), tau = m / (1 g cm^-2),
    log m from -5 to 1 at n points."""
    logm = [-5.0 + 6.0 * k / (n - 1) for k in range(n)]
    T = [(0.75 * 5778.0 ** 4 * (10 ** lm + 2.0 / 3.0)) ** 0.25 for lm in logm]
    ne = [1e12 * 10 ** (0.8 * lm) for lm in logm]      # first guess only
    return logm, T, ne, [vz] * n, [1.0] * n


def falc_v(vz):
    """FAL-C (Atmos/FALC_82.atmos) as is, with uniform V = vz [km/s]."""
    out, ndep, row = [], None, 0
    for l in open(os.path.join(src, "Atmos", "FALC_82.atmos")).read().splitlines():
        w = l.split()
        if ndep is None or not w or l.lstrip().startswith("*") or row >= ndep:
            if ndep is None and w and not l.lstrip().startswith("*") and len(w) == 1:
                try:
                    ndep = int(w[0])
                except ValueError:
                    pass
            out.append(l)
            continue
        w[3] = f"{vz:.12E}"
        out.append("  " + "  ".join(w))
        row += 1
    return "\n".join(out) + "\n"


def write_B(path, B, gamma, chi):
    with open(path, "wb") as f:
        for q in (B, gamma, chi):
            f.write(struct.pack(f">{len(q)}d", *q))


# --- Running rhf1d ----------------------------------------------------------

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


def run(name, atmos_text, ndep, field, vertical=False, solve_ne=False, mo="TRUE"):
    """field: (B [T], gamma, chi) per depth, or None for B = 0. Returns
    (lambda, [I, Q, U, V]) with each Stokes array [nspect][nray]."""
    d = os.path.join(work, "runs", name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    open(os.path.join(d, "atmos.atmos"), "w").write(atmos_text)
    B, gam, chi = field if field else ([0.0] * ndep,) * 3
    write_B(os.path.join(d, "B.dat"), B, gam, chi)
    shutil.copy(os.path.join(HERE, "TI_fict.atom"), d)
    kw = open(os.path.join(src, "rhf1d", "run", "keyword.input")).read()
    for k, v in [("ATMOS_FILE", "atmos.atmos"), ("STOKES_MODE", "FULL_STOKES"),
                 ("STOKES_INPUT", "B.dat"), ("N_MAX_ITER", "0"),
                 ("N_MAX_SCATTER", "1"), ("N_THREADS", "1"),
                 ("PRINT_CPU", "FALSE"), ("VACUUM_TO_AIR", "FALSE"),
                 ("MAGNETO_OPTICAL", mo), ("BACKGROUND_POLARIZATION", "FALSE"),
                 ("STARTING_J", "NEW_J"), ("KURUCZ_DATA", "none"),
                 ("SOLVE_NE", "ONCE" if solve_ne else "NONE")]:
        kw = setkw(kw, k, v)
    kw = setkw(kw, "ANGLE_SET", "SET_VERTICAL") if vertical else setkw(kw, "NRAYS", str(NRAYS))
    open(os.path.join(d, "keyword.input"), "w").write(kw)
    open(os.path.join(d, "atoms.input"), "w").write(
        "2\n  ../../Atoms/H_6.atom  PASSIVE  LTE_POPULATIONS\n"
        "  TI_fict.atom  ACTIVE  LTE_POPULATIONS  pops.TI.out\n")
    shutil.copy(os.path.join(src, "rhf1d", "run", "molecules.input"), d)
    r = subprocess.run([rhf1d], cwd=d, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, text=True)
    open(os.path.join(d, "out.log"), "w").write(r.stdout)
    if r.returncode != 0 or "ERROR" in r.stdout:
        tail = "\n".join(r.stdout.splitlines()[-40:])   # FP trap backtrace (B-33)
        sys.exit(f"rhf1d failed in {d} (status {r.returncode}); end of out.log:\n{tail}")
    raw = open(os.path.join(d, "spectrum.out"), "rb").read()
    ns = struct.unpack(">i", raw[:4])[0]
    lam = struct.unpack(f">{ns}d", raw[4:4 + 8 * ns])
    nr = 1 if vertical else NRAYS
    size, pos = ns * nr, 4 + 8 * ns
    I = struct.unpack(f">{size}d", raw[pos:pos + 8 * size]); pos += 8 * size + 12
    QUV = [struct.unpack(f">{size}d", raw[pos + 8 * size * i:pos + 8 * size * (i + 1)])
           for i in range(3)]
    if pos + 24 * size != len(raw):
        sys.exit(f"{d}/spectrum.out: unexpected size")
    shape = lambda q: [q[l * nr:(l + 1) * nr] for l in range(ns)]
    return lam, [shape(q) for q in [I] + QUV]


failures = []


def check(name, value, tol, at_least=False):
    """value <= tol, or value >= tol with at_least."""
    ok = value >= tol if at_least else value <= tol
    if not ok:
        failures.append(name)
    if a.verbose or not ok:
        print(f"{'ok  ' if ok else 'FAIL'} {name}: {value:.3e} "
              f"({'at least' if at_least else 'tol'} {tol:.1e})")


def maxabs(q):
    return max(abs(x) for row in q for x in row)


def maxdiff(p, q, sign=1.0):
    return max(abs(x - sign * y) for rp, rq in zip(p, q) for x, y in zip(rp, rq))


def field(n):
    """Depth-dependent field: 200 G at the top to 2 kG at the bottom, gamma
    from 70 to 20 deg, chi from 100 to 10 deg."""
    t = [k / (n - 1) for k in range(n)]
    return ([1e-4 * (200.0 + 1800.0 * s * s) for s in t],
            [math.radians(70.0 - 50.0 * s) for s in t],
            [math.radians(100.0 - 90.0 * s) for s in t])


def rotate(sp, turn):
    """Q, U of sp turned by 2 turn."""
    c2, s2 = math.cos(2 * turn), math.sin(2 * turn)
    Q = [[c2 * q - s2 * u for q, u in zip(rq, ru)] for rq, ru in zip(sp[1], sp[2])]
    U = [[s2 * q + c2 * u for q, u in zip(rq, ru)] for rq, ru in zip(sp[1], sp[2])]
    return Q, U


# --- 1-3: exact symmetries on FAL-C ------------------------------------------

ndep = 82
B0, G0, C0 = field(ndep)
atm = falc_v(0.0)

lam, free = run("free", atm, ndep, None)
for i, s in enumerate("QUV"):
    check(f"no field: max|{s}| / max I", maxabs(free[i + 1]) / maxabs(free[0]), 0.0)

for mo in ("TRUE", "FALSE"):
    _, base = run(f"base_mo{mo}", atm, ndep, (B0, G0, C0), mo=mo)
    scale = [maxabs(q) for q in base]
    if min(scale[1:]) < 1e-3 * scale[0]:
        sys.exit("the field gives too little polarization for a meaningful test")
    _, mir = run(f"mirror_mo{mo}", atm, ndep,
                 (B0, [math.pi - g for g in G0], [math.pi - c for c in C0]), mo=mo)
    for i, s in enumerate("IQUV"):
        check(f"mirror y -> -y, MO {mo}: {s}",
              maxdiff(mir[i], base[i], -1.0 if s in "UV" else 1.0) / scale[i], 1e-12)
    if mo == "FALSE":
        _, neg = run("minusB", atm, ndep,
                     (B0, [math.pi - g for g in G0], [c + math.pi for c in C0]), mo=mo)
        for i, s in enumerate("IQUV"):
            check(f"B -> -B, no MO: {s}",
                  maxdiff(neg[i], base[i], -1.0 if s == "V" else 1.0) / scale[i], 1e-12)

# --- 4 on FAL-C: vertical ray, field azimuth turned by 35 deg ----------------

turn = math.radians(35.0)                       # also used on the grey atmosphere
_, vb = run("rot_falc", atm, ndep, (B0, G0, C0), vertical=True)
_, vt = run("rot_falc_turned", atm, ndep, (B0, G0, [c + turn for c in C0]), vertical=True)
Qr, Ur = rotate(vb, turn)
check("frame rotation, FAL-C: I", maxdiff(vt[0], vb[0]) / maxabs(vb[0]), 1e-12)
check("frame rotation, FAL-C: V", maxdiff(vt[3], vb[3]) / maxabs(vb[3]), 1e-12)
check("frame rotation, FAL-C: Q, U",
      max(maxdiff(vt[1], Qr), maxdiff(vt[2], Ur)) / max(maxabs(vb[1]), maxabs(vb[2])), 1e-12)

# --- 4, 5: frame rotation and weak field on smooth grey atmospheres ---------

Bw, gw, cw = 1e-4, math.radians(30.0), math.radians(20.0)   # 1 G
dv = 0.02                                                   # km/s
errs, rots = [], []
for n in (40, 80, 160, 320):
    def text(vz):
        path = os.path.join(work, "tmp.atmos")
        write_atmos_text(path, *grey(n, vz), f"GREY_{n}")
        return open(path).read()

    Bn, Gn, Cn = field(n)
    _, vb = run(f"rot_{n}", text(0.0), n, (Bn, Gn, Cn), vertical=True, solve_ne=True)
    _, vt = run(f"rot_{n}_turned", text(0.0), n, (Bn, Gn, [c + turn for c in Cn]),
                vertical=True, solve_ne=True)
    Qr, Ur = rotate(vb, turn)
    rot = max(maxdiff(vt[0], vb[0]) / maxabs(vb[0]), maxdiff(vt[3], vb[3]) / maxabs(vb[3]),
              max(maxdiff(vt[1], Qr), maxdiff(vt[2], Ur)) / max(maxabs(vb[1]), maxabs(vb[2])))
    rots.append(rot)
    if a.verbose:
        print(f"     N = {n:3d}: vertical ray, chi + 35 deg: largest deviation {rot:.2e}")

    lam, sp = run(f"wf_{n}", text(0.0), n, ([Bw] * n, [gw] * n, [cw] * n), solve_ne=True)
    I = {s: run(f"wf_{n}_v{s:+g}", text(s * dv), n, None, solve_ne=True)[1][0]
         for s in (-2, -1, 1, 2)}
    err = 0.0
    for mu in range(NRAYS):
        muz = MU[mu]
        cos_th = math.sqrt(1 - muz * muz) * math.sin(gw) * math.cos(cw) + muz * math.cos(gw)
        for lam0, g in LINES:
            idx = [l for l in range(len(lam)) if abs(lam[l] - lam0) < 1.0]
            uB = E_CHARGE * lam0 * 1e-9 * Bw / (4 * math.pi * M_ELECTRON)    # m/s
            V = [sp[3][l][mu] for l in idx]
            dIdv = [(8 * (I[1][l][mu] - I[-1][l][mu]) - (I[2][l][mu] - I[-2][l][mu]))
                    / (12 * dv * 1e3) for l in idx]
            Vwf = [-g * uB * cos_th / muz * d for d in dIdv]
            e = max(abs(x - y) for x, y in zip(V, Vwf)) / max(abs(x) for x in V)
            err = max(err, e)
            if a.verbose:
                print(f"     N = {n:3d}  mu = {muz:.3f}  {lam0:.2f} nm: "
                      f"max|V - V_wf| / max|V| = {e:.2e}")
    errs.append(err)

check("weak field, N = 320: max|V - V_wf| / max|V|", errs[-1], 1e-5)
for n, r in zip((40, 80, 160, 320), rots):
    check(f"frame rotation, N = {n}: largest deviation", r, 1e-12)
for n, e0, e1 in zip((80, 160, 320), errs, errs[1:]):
    check(f"weak field, error ratio N = {n // 2} -> {n}", e0 / e1, 8.0, at_least=True)

if not a.work:
    shutil.rmtree(work)
print(f"{'FAILED: ' + ', '.join(failures) if failures else 'all passed'}")
sys.exit(1 if failures else 0)
