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
 3. A run stopped with STOP_RH during iteration 3 ends after 3 iterations
    and is consistent to 1e-5. With STOP_RH present at the start no iteration
    runs ("after 0 iterations").
"""

import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.output import read_atmos_out, read_first_molecule, read_pops_array
from rhtest.run import EXIT_NOT_CONVERGED, EXIT_OK, execute_streaming, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d",))
work = WorkDir(a.src, a.work)
ck = Checker()


def run(name, keywords, stop_at=None):
    d = prepare(work, name, "rhf1d",
                [("ATMOS_FILE", "../../Atmos/FALC_82.atmos"),
                 ("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"),
                 ("N_THREADS", "1"), ("PRINT_CPU", "FALSE"), ("STARTING_J", "NEW_J"),
                 ("KURUCZ_DATA", "none"), ("SOLVE_NE", "ITERATION"),
                 ("COLL_RELAX_NSTEP", "0"), ("SOLVER", "NG"),
                 ("HYDROSTATIC", "TRUE")] + keywords,
                [("../../Atoms/H_6.atom", "ACTIVE", "ZERO_RADIATION", "pops.H.out"),
                 ("../../Atoms/CaII.atom", "ACTIVE", "ZERO_RADIATION", "pops.CA.out")])
    # The messages go to stderr (unbuffered), so they arrive as they are
    # written. With stop_at = k, STOP_RH is created when the line of
    # iteration k with its H residual appears: three PRD sweeps (each of
    # the order of 0.1 s) come before the end of that iteration, where the
    # stop request is read, so it is read at the end of iteration k.
    if stop_at == 0:                       # requested before the first iteration
        open(os.path.join(d, "STOP_RH"), "w").close()
    seen = []                              # [True] after the header of iteration stop_at

    def on_line(line):
        if stop_at is None:
            return
        if line.strip() == f"-- Iteration {stop_at:3d}":
            seen[:] = [True]
        elif seen and line.startswith(" H , delta"):
            open(os.path.join(d, "STOP_RH"), "w").close()
            seen[:] = []

    r = execute_streaming(a.rhf1d, d, on_line, timeout=900).require()
    return d, r.stdout, r.returncode


def inconsistency(d):
    """max over the depth points without H2 of |sum_i n_CaII,i / (abund_Ca nHtot) - 1|.
    nHtot is the sum of the H populations of atmos.out only where no H2 is
    present: the depth points with T outside [Tmin, Tmax] of H2.molecule, the
    condition of the hydrostatic step (molecules.out; the first molecule must
    be H2)."""
    nl, ns, n = read_pops_array(os.path.join(d, "pops.CA.out"))
    T, nh, ab = read_atmos_out(os.path.join(d, "atmos.out"), ns)
    ident, tmin, tmax = read_first_molecule(os.path.join(d, "molecules.out"))
    if ident != "H2":
        sys.exit(f"first molecule of {d}/molecules.out is {ident!r}, not H2")
    ks = [k for k in range(ns) if not tmin <= T[k] <= tmax]
    if len(ks) < 5:
        sys.exit(f"only {len(ks)} depth points outside the H2 range [{tmin}, {tmax}] K:"
                 f" the test cannot see the problem")
    return max(abs(sum(n[i * ns + k] for i in range(nl)) / (ab["CA"] * nh[k]) - 1.0)
               for k in ks)


TOL = 1e-5      # the starting state is off by 1.1e-6 (N_MAX_ITER = 1), unrelated to the hydrostatic step

d, out, rc = run("unconverged", [("N_MAX_ITER", "3"), ("ITER_LIMIT", "1e-3")])
c = inconsistency(d)
ok = rc == EXIT_NOT_CONVERGED and "NOT CONVERGED after 3 iterations" in out and c < TOL
ck.report(ok, f"N_MAX_ITER = 3: status {rc}, Ca II vs nHtot of the output {c:.2E}"
              f" ({'ok' if ok else 'FAIL'}, tolerance {TOL:.0E})",
          f"unconverged run inconsistent by {c:.2E}")

d, out, rc = run("converged", [("N_MAX_ITER", "100"), ("ITER_LIMIT", "5e-2")])
c = inconsistency(d)
ok = rc == EXIT_OK and "Converged in" in out and c < TOL
ck.report(ok, f"converged: status {rc}, Ca II vs nHtot of the output {c:.2E}"
              f" ({'ok' if ok else 'FAIL'})",
          f"converged run: status {rc}, inconsistency {c:.2E}")

# STOP_RH created during iteration 3: the run ends after 3 iterations, like
# N_MAX_ITER = 3, and must be as consistent
d, out, rc = run("stopped", [("N_MAX_ITER", "100"), ("ITER_LIMIT", "1e-3")], stop_at=3)
c = inconsistency(d)
ok = rc == EXIT_NOT_CONVERGED and "Stopped on request (STOP_RH) after 3 iterations" in out \
    and c < TOL
ck.report(ok, f"STOP_RH in iteration 3: status {rc}, Ca II vs nHtot of the output {c:.2E}"
              f" ({'ok' if ok else 'FAIL'})",
          f"stopped run: status {rc}, inconsistency {c:.2E}")

# STOP_RH present at the start: no iteration runs
d, out, rc = run("stopped0", [("N_MAX_ITER", "100"), ("ITER_LIMIT", "1e-3")], stop_at=0)
ok = rc == EXIT_NOT_CONVERGED and "Stopped on request (STOP_RH) after 0 iterations" in out
ck.report(ok, f"STOP_RH at the start: status {rc} ({'ok' if ok else 'FAIL'})",
          f"stop at the start: status {rc}")

ck.finish("hydrostatic_final", work)
