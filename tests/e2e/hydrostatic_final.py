#!/usr/bin/env python3
"""End-to-end test of the final state of SOLVER = NG with HYDROSTATIC = TRUE (N-16).

usage: hydrostatic_final.py --rhf1d PATH [--src DIR] [--work DIR]

H_6 and Ca II active on FAL-C, SOLVE_NE = ITERATION. The hydrostatic step
(Hydrostatic) rescales only hydrogen to the new nHtot, so a state that is
written right after it has the other active atoms off by the last change of
nHtot. The output of a run must be consistent at every depth point: sum of
the populations of Ca II = abundance(Ca) * nHtot, with the abundance from the
element list of atmos.out and nHtot the hydrogen nuclei: the sum of the H
populations in atmos.out (atoms and protons) plus twice the H2 density of
molecules.out. Up to 6.2e-4 of nHtot is in H2 on FAL-C; the hydrogen in H-
and other molecules (1e-6) is left out. Before N-44 the hydrostatic step
counted H2 once, which the converged case sees (6.2e-4).
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
    seen = False                           # the header of iteration stop_at has appeared

    def on_line(line):
        nonlocal seen
        if stop_at is None:
            return
        if line.strip() == f"-- Iteration {stop_at:3d}":
            seen = True
        elif seen and line.startswith(" H , delta"):
            open(os.path.join(d, "STOP_RH"), "w").close()
            seen = False

    r = execute_streaming(a.rhf1d, d, timeout=900, on_line=on_line).require()
    return d, r.stdout, r.returncode


def inconsistency(d):
    """max over the depth points of |sum_i n_CaII,i / (abund_Ca nHtot) - 1|,
    nHtot = sum of the H populations of atmos.out + 2 n_H2 (molecules.out;
    the first molecule must be H2)."""
    nl, ns, n = read_pops_array(os.path.join(d, "pops.CA.out"))
    T, nh, ab = read_atmos_out(os.path.join(d, "atmos.out"), ns)
    path = os.path.join(d, "molecules.out")
    try:
        ident, tmin, tmax, nH2 = read_first_molecule(path, ns)
        problem = None if ident == "H2" else f"it is {ident!r}"
    except ValueError as e:
        problem = str(e)
    if problem:
        sys.exit(f"{path}: the first molecule must be H2 ({problem})")
    if max(nH2[k] / nh[k] for k in range(ns)) < 1e-4:
        sys.exit("less than 1e-4 of the hydrogen in H2 at every depth point:"
                 " the test cannot see N-44")
    return max(abs(sum(n[i * ns + k] for i in range(nl)) / (ab["CA"] * (nh[k] + 2.0 * nH2[k]))
                   - 1.0) for k in range(ns))


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
