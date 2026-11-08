#!/usr/bin/env python3
"""Residual of statistical equilibrium in a finished RH run (any geometry).

usage: se_residual.py RUNDIR [-v] [--json]

For every active atom with atom.X.out, pops.X.out, radrate.X.out and
collrate.X.out in RUNDIR, evaluates at each spatial point k and level c the
SE residual F_SE (S-1 in docs/SOLVER.md)

    F_c = sum_b n_b P_bc - n_c sum_b P_cb,

with P the collisional rates plus the radiative rates RH wrote for the
final formal solution, and from it

    r_c = |F_c| / (n_c sum_b P_cb)       relative imbalance of level c: the
                                         relative change one Jacobi step
                                         n_c <- in_c / out-rate_c would make
    s_c = r_c n_c / n_tot                the same in units of the total
                                         population n_tot = sum_c n*_c

and reports (the norm is still open in SOLVER.md, so all of them):

    max_rel           max r_c over all levels and points
    max_rel@1e-10, max_rel@1e-06, max_rel@1e-03
                      max r_c over levels with n_c / n_tot >= threshold
    rms               RMS of s_c over all levels and points
    cons              conservation row: max over points of
                      |sum_c n_c - n_tot| / n_tot

After convergence they should be of the order of ITER_LIMIT or less: the
iteration may take any path, but its end point must satisfy the plain rate
equations with the radiation field the populations produce (B-6).

What the rates are: those of the final sweep after Iterate() (rhf1d.c and
friends), a formal solution with the populations fixed and no operator,
with the background-scattering J iterated up to N_MAX_SCATTER times (to
ITER_LIMIT) and the PRD rho carried over from the last iteration (not
recomputed). An unconverged rho or background J therefore shows up in
F_SE; ne is taken as given.

The run needs RADRATE_OUTPUT and COLLRATE_OUTPUT set, and no polarized
final solution (STOKES_MODE = NO_STOKES): the polarized formal solution
after convergence changes J, so its rates do not belong to the populations.
Fixed transitions are not written by RH and are not included (warning).
-v prints the largest residual of every level; --json prints one JSON
line {ID: metrics, ..., "all": metrics over all atoms}.
"""
import glob
import json
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_run import read_pops  # noqa: E402

THRESHOLDS = (1e-10, 1e-6, 1e-3)


class XDR:
    def __init__(self, data):
        self.b, self.p = data, 0

    def int(self):
        v = struct.unpack(">i", self.b[self.p:self.p + 4])[0]
        self.p += 4
        return v

    def double(self):
        v = struct.unpack(">d", self.b[self.p:self.p + 8])[0]
        self.p += 8
        return v

    def counted_string(self):
        n = self.int()                      # xdr_short, 4 bytes
        if n == 0:
            return ""
        m = self.int()                      # xdr_string's own length
        s = self.b[self.p:self.p + m].decode()
        self.p += (m + 3) // 4 * 4
        return s


def read_atom(path):
    x = XDR(open(path, "rb").read())
    x.int()                                 # active
    Nlevel, Nline, Ncont, Nfixed = (x.int() for _ in range(4))
    x.double(), x.double()                  # abundance, weight
    labels = [x.counted_string().strip() for _ in range(Nlevel)]
    x.p += 8 * Nlevel * 2 + 4 * Nlevel      # g, E, stage
    trans, lines = [], []
    for _ in range(Nline):
        x.int()
        i, j = x.int(), x.int()
        Nlambda, _ = x.int(), x.int()
        lambda0 = x.double()
        x.int(), x.double()
        trans.append((i, j))
        lines.append((lambda0, Nlambda))
    for _ in range(Ncont):
        x.int()
        i, j = x.int(), x.int()
        x.int(), x.int(), x.double(), x.int(), x.double()
        trans.append((i, j))
    return Nlevel, trans, Nfixed, labels, lines


def residual(rundir, ID, verbose, quiet=False):
    """Metrics of one atom; also returns the arrays r, s for pooling."""
    Nlevel, trans, Nfixed, labels, _ = read_atom(os.path.join(rundir, f"atom.{ID}.out"))
    n = read_pops(os.path.join(rundir, f"pops.{ID}.out"))
    nstar = read_pops(os.path.join(rundir, f"pops.{ID}.out"), lte=True)
    Nspace = n.shape[1]
    R = np.fromfile(os.path.join(rundir, f"radrate.{ID}.out"), ">f8")
    C = np.fromfile(os.path.join(rundir, f"collrate.{ID}.out"), ">f8")
    if R.size != 2 * len(trans) * Nspace or C.size != Nlevel**2 * Nspace:
        sys.exit(f"{ID}: rate files do not match the atom ({R.size}, {C.size})")
    R = R.reshape(len(trans), 2, Nspace)
    C = C.reshape(Nlevel, Nlevel, Nspace)   # C[i][j] = rate j -> i

    inflow = np.zeros((Nlevel, Nspace))
    outflow = np.zeros((Nlevel, Nspace))
    for c in range(Nlevel):
        for b in range(Nlevel):
            if b != c:
                inflow[c] += n[b] * C[c, b]
                outflow[c] += n[c] * C[b, c]
    for (i, j), (Rij, Rji) in zip(trans, R):
        up, down = n[i] * Rij, n[j] * Rji
        outflow[i] += up
        inflow[j] += up
        outflow[j] += down
        inflow[i] += down
    # --- levels with no outflow (population underflowed to 0) are skipped
    r = np.where(outflow > 0, np.abs(inflow - outflow) / np.where(outflow > 0, outflow, 1), 0.0)
    ntot = nstar.sum(axis=0)
    frac = n / ntot
    s = r * frac
    m = metrics(r, s, frac)
    m["cons"] = float(np.max(np.abs(n.sum(axis=0) - ntot) / ntot))
    if not quiet:
        c, k = np.unravel_index(np.argmax(r), r.shape)
        note = f"  ({Nfixed} fixed transitions not included)" if Nfixed else ""
        print(f"{ID}: max residual {r.max():.3e} at level {c} ({labels[c]}), "
              f"point {k}{note}")
        print("    " + ", ".join(f"{key} {v:.3e}" for key, v in m.items()))
        if verbose:
            for c in range(Nlevel):
                print(f"  level {c:3d} {labels[c]:30s} {r[c].max():.3e}")
    return m, r, s, frac


def metrics(r, s, frac):
    m = {"max_rel": float(r.max())}
    for t in THRESHOLDS:
        sel = frac >= t
        m[f"max_rel@{t:.0e}"] = float(r[sel].max()) if sel.any() else 0.0
    m["rms"] = float(np.sqrt(np.mean(s**2)))
    return m


def run_metrics(rundir, verbose=False, quiet=True):
    """{ID: metrics, "all": pooled metrics}, or {} if no atom has rates."""
    out, pool = {}, []
    for f in sorted(glob.glob(os.path.join(rundir, "atom.*.out"))):
        ID = os.path.basename(f)[5:-4]
        if all(os.path.exists(os.path.join(rundir, f"{p}.{ID}.out"))
               for p in ("pops", "radrate", "collrate")):
            m, r, s, frac = residual(rundir, ID, verbose, quiet)
            out[ID] = m
            pool.append((r.ravel(), s.ravel(), frac.ravel()))
    if pool:
        r, s, frac = (np.concatenate(a) for a in zip(*pool))
        out["all"] = metrics(r, s, frac)
        out["all"]["cons"] = max(m["cons"] for k, m in out.items() if k != "all")
    return out


def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        sys.exit(__doc__)
    rundir, verbose = sys.argv[1], "-v" in sys.argv[2:]
    kw = os.path.join(rundir, "keyword.input")
    if os.path.exists(kw):
        for l in open(kw):
            f = l.split("#")[0].split("=")
            if len(f) == 2 and f[0].strip() == "STOKES_MODE" and f[1].strip() != "NO_STOKES":
                print(f"warning: STOKES_MODE = {f[1].strip()}: the polarized final "
                      "solution changes J, the rates do not belong to the populations",
                      file=sys.stderr)
    as_json = "--json" in sys.argv[2:]
    out = run_metrics(rundir, verbose, quiet=as_json)
    if not out:
        sys.exit(f"no atom with pops, radrate and collrate output in {rundir}")
    if as_json:
        print(json.dumps(out))
    elif len(out) > 2:
        print("all: " + ", ".join(f"{key} {v:.3e}" for key, v in out["all"].items()))


if __name__ == "__main__":
    main()
