#!/usr/bin/env python3
"""Residual of statistical equilibrium in a finished RH run (any geometry).

usage: se_residual.py RUNDIR [-v]

For every active atom with atom.X.out, pops.X.out, radrate.X.out and
collrate.X.out in RUNDIR, evaluates at each spatial point and level c

    r_c = | sum_b n_b P_bc - n_c sum_b P_cb | / (n_c sum_b P_cb),

with P the collisional rates plus the radiative rates RH wrote for the
final formal solution, and prints the largest r_c. After convergence it
should be of the order of ITER_LIMIT or less: the iteration may take any
path, but its end point must satisfy the plain rate equations with the
radiation field the populations produce (B-6).

The run needs RADRATE_OUTPUT and COLLRATE_OUTPUT set, and no polarized
final solution (STOKES_MODE = NO_STOKES): the polarized formal solution
after convergence changes J, so its rates do not belong to the populations.
Fixed transitions are not written by RH and are not included (warning).
-v prints the largest residual of every level.
"""
import glob
import os
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analyze_run import read_pops  # noqa: E402


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
    trans = []
    for _ in range(Nline):
        x.int()
        i, j = x.int(), x.int()
        x.int(), x.int(), x.double(), x.int(), x.double()
        trans.append((i, j))
    for _ in range(Ncont):
        x.int()
        i, j = x.int(), x.int()
        x.int(), x.int(), x.double(), x.int(), x.double()
        trans.append((i, j))
    return Nlevel, trans, Nfixed, labels


def residual(rundir, ID, verbose):
    Nlevel, trans, Nfixed, labels = read_atom(os.path.join(rundir, f"atom.{ID}.out"))
    n = read_pops(os.path.join(rundir, f"pops.{ID}.out"))
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
    c, k = np.unravel_index(np.argmax(r), r.shape)
    note = f"  ({Nfixed} fixed transitions not included)" if Nfixed else ""
    print(f"{ID}: max residual {r.max():.3e} at level {c} ({labels[c]}), "
          f"point {k}{note}")
    if verbose:
        for c in range(Nlevel):
            print(f"  level {c:3d} {labels[c]:30s} {r[c].max():.3e}")
    return r.max()


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
    found = False
    for f in sorted(glob.glob(os.path.join(rundir, "atom.*.out"))):
        ID = os.path.basename(f)[5:-4]
        if all(os.path.exists(os.path.join(rundir, f"{p}.{ID}.out"))
               for p in ("pops", "radrate", "collrate")):
            residual(rundir, ID, verbose)
            found = True
    if not found:
        sys.exit(f"no atom with pops, radrate and collrate output in {rundir}")


if __name__ == "__main__":
    main()
