#!/usr/bin/env python3
"""End-to-end test of the final state of SOLVER = NG with HYDROSTATIC = TRUE (N-16).

usage: hydrostatic_final.py --rhf1d PATH [--src DIR] [--work DIR]

H_6 and Ca II active on FAL-C, SOLVE_NE = ITERATION. The hydrostatic step
(Hydrostatic) rescales only hydrogen to the new nHtot, so a state that is
written right after it has the other active atoms off by the last change of
nHtot. The output of a run must be consistent: sum of the populations of
Ca II = abundance(Ca) * nHtot, with nHtot the sum of the H populations in
atmos.out (at the depth points without H2) and the abundance from its element
list.
 1. A run that stops at N_MAX_ITER = 3 (exit status 1, NOT CONVERGED) is
    consistent to 1e-5. Without the fix the last hydrostatic step leaves
    it off by 0.1.
 2. A run that converges (ITER_LIMIT = 5e-2) is consistent to 1e-5.
"""

import argparse
import os
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
    """(Nlevel, Nspace, NLTE populations as a flat list) of a pops.*.out file
    (struct only, as the CI build jobs have no numpy)."""
    b = open(path, "rb").read()
    for off in range(4, min(len(b) - 8, 512), 4):
        nl, ns = struct.unpack(">2i", b[off:off + 8])
        if 0 < nl < 100000 and 0 < ns and off + 8 + 16 * nl * ns == len(b):
            return nl, ns, struct.unpack(f">{nl * ns}d", b[off + 8:off + 8 + 8 * nl * ns])
    raise ValueError(f"cannot parse {path}")


def read_atmos(path, nspace):
    """nHtot (sum of the H populations) and the element abundances of atmos.out
    (writeAtmos: ints Nlevel, Nelem, bool moving; T, ne, vturb; H populations;
    counted atmosphere ID; Nelem x (counted ID, weight, abundance))."""
    b = open(path, "rb").read()
    nlev, nelem, _ = struct.unpack(">3i", b[:12])
    off = 12 + 3 * 8 * nspace
    h = struct.unpack(f">{nlev * nspace}d", b[off:off + 8 * nlev * nspace])
    off += 8 * nlev * nspace

    def cstring(off):                      # xdr_short (4 bytes), then xdr_string
        n = struct.unpack(">i", b[off:off + 4])[0]     # (length, padded data)
        if n == 0:
            return "", off + 4
        m = struct.unpack(">i", b[off + 4:off + 8])[0]
        s = b[off + 8:off + 8 + m].decode()
        return s, off + 8 + 4 * ((m + 3) // 4)
    _, off = cstring(off)
    ab = {}
    for _ in range(nelem):
        el, off = cstring(off)
        w, abund = struct.unpack(">2d", b[off:off + 16])
        off += 16
        ab[el.strip().upper()] = abund
    nh = [sum(h[l * nspace + k] for l in range(nlev)) for k in range(nspace)]
    return nh, ab


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


def run(name, keywords):
    d = os.path.join(work, "runs", name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    kw = open(os.path.join(src, "rhf1d", "run", "keyword.input")).read()
    for k, v in [("ATMOS_FILE", "../../Atmos/FALC_82.atmos"),
                 ("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
                 ("N_THREADS", "1"), ("PRINT_CPU", "FALSE"), ("STARTING_J", "NEW_J"),
                 ("KURUCZ_DATA", "none"), ("SOLVE_NE", "ITERATION"),
                 ("COLL_RELAX_NSTEP", "0"), ("SOLVER", "NG"),
                 ("HYDROSTATIC", "TRUE")] + keywords:
        kw = setkw(kw, k, v)
    open(os.path.join(d, "keyword.input"), "w").write(kw)
    open(os.path.join(d, "atoms.input"), "w").write(
        "2\n  ../../Atoms/H_6.atom  ACTIVE  ZERO_RADIATION  pops.H.out\n"
        "  ../../Atoms/CaII.atom  ACTIVE  ZERO_RADIATION  pops.CA.out\n")
    for f in ("molecules.input", "ray.input"):
        if os.path.exists(os.path.join(src, "rhf1d", "run", f)):
            shutil.copy(os.path.join(src, "rhf1d", "run", f), d)
    try:
        r = subprocess.run([os.path.abspath(a.rhf1d)], cwd=d, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, timeout=900)
    except subprocess.TimeoutExpired:
        sys.exit(f"rhf1d timed out in {d}")
    open(os.path.join(d, "out.log"), "w").write(r.stdout)
    if r.returncode not in (0, 1):
        sys.exit(f"rhf1d failed in {d} (status {r.returncode})")
    return d, r.stdout, r.returncode


def inconsistency(d):
    """max over the depth points without H2 of |sum_i n_CaII,i / (abund_Ca nHtot) - 1|.
    nHtot is the sum of the H populations of atmos.out only where no H2 is
    present: the depth points with T outside [Tmin, Tmax] of H2.molecule, the
    condition of the hydrostatic step (molecules.out: Nmolecule, counted ID of
    the first molecule, Nv, NJ, Ediss, Tmin, Tmax)."""
    nl, ns, n = read_pops(os.path.join(d, "pops.CA.out"))
    nh, ab = read_atmos(os.path.join(d, "atmos.out"), ns)
    mol = open(os.path.join(d, "molecules.out"), "rb").read()
    tmin, tmax = struct.unpack(">2d", mol[32:48])
    atm = open(os.path.join(d, "atmos.out"), "rb").read()
    T = struct.unpack(f">{ns}d", atm[12:12 + 8 * ns])
    ks = [k for k in range(ns) if not tmin <= T[k] <= tmax]
    return max(abs(sum(n[i * ns + k] for i in range(nl)) / (ab["CA"] * nh[k]) - 1.0)
               for k in ks)


failed = []
TOL = 1e-5      # the starting state is off by 1.1e-6 (N_MAX_ITER = 1), unrelated to the hydrostatic step

d, out, rc = run("unconverged", [("N_MAX_ITER", "3"), ("ITER_LIMIT", "1e-3")])
c = inconsistency(d)
ok = rc == 1 and "NOT CONVERGED after 3 iterations" in out and c < TOL
print(f"N_MAX_ITER = 3: status {rc}, Ca II vs nHtot of the output {c:.2E}"
      f" ({'ok' if ok else 'FAIL'}, tolerance {TOL:.0E})")
if not ok:
    failed.append(f"unconverged run inconsistent by {c:.2E}")

d, out, rc = run("converged", [("N_MAX_ITER", "100"), ("ITER_LIMIT", "5e-2")])
c = inconsistency(d)
ok = rc == 0 and "Converged in" in out and c < TOL
print(f"converged: status {rc}, Ca II vs nHtot of the output {c:.2E}"
      f" ({'ok' if ok else 'FAIL'})")
if not ok:
    failed.append(f"converged run: status {rc}, inconsistency {c:.2E}")

if failed:
    sys.exit("FAILED: " + "; ".join(failed) + f"\n runs in {work}")
if not a.work:
    shutil.rmtree(work)
print("hydrostatic_final: all checks passed")
