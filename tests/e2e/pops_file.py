#!/usr/bin/env python3
"""End-to-end test of reading populations files with OLD_POPULATIONS (B-60, B-61).

usage: pops_file.py --rhf1d PATH [--src DIR] [--work DIR] [-v]

FAL-C (rhf1d/run), Ca II active, 1 thread, ITER_LIMIT 1e-4.
Run A (SOLVE_NE = NONE, ZERO_RADIATION) writes pops.CA.out.

B-60: readPopulations read the LTE populations n* from the file as well,
after SetLTEQuantities had computed them for the current atmosphere and
electron density, and an active atom kept them for the whole run. Now only
n is read:

 1. A2 (SOLVE_NE = NONE, OLD_POPULATIONS from A): n* as in A (the same
    electron density), n as in A to the convergence tolerance.
 2. B (SOLVE_NE = ONCE, OLD_POPULATIONS from A) and C (SOLVE_NE = ONCE,
    ZERO_RADIATION): the electron density differs from A by -15 % to +9 %;
    B has the n* of C (exactly), not those of A, and converges to the
    populations of C. (Before: B kept A's n* and its populations differed
    from C's by up to 15 %.)

B-61: a populations file of another atmosphere stopped nothing (a warning
when the IDs differed, nothing when the current ID was part of the file's);
another number of levels or depths gave "Unable to read from input file".
Now each stops at input with a message and status 70:

 3. ID of the atmosphere changed to FALC_XX, and to FALC_8 (part of the ID
    of the file, FALC_82).
 4. The populations of Ca II read for He (another number of levels).

Python 3 standard library only. Exit status 0 = pass.
"""
import os
import shutil
import struct
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker, maxrel
from rhtest.cli import parse
from rhtest.output import read_bytes, xdr_counted_string
from rhtest.run import EXIT_FATAL, RAN, execute, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d",), verbose=True)
work = WorkDir(a.src, a.work)
ck = Checker(a.verbose)

H6 = ("../../Atoms/H_6.atom", "PASSIVE", "LTE_POPULATIONS")
CA = "../../Atoms/CaII.atom"
KWS = [("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"), ("N_THREADS", "1"),
       ("PRINT_CPU", "FALSE"), ("STARTING_J", "NEW_J"), ("KURUCZ_DATA", "none"),
       ("N_MAX_ITER", "200"), ("ITER_LIMIT", "1.0E-4")]


def read_pops(path):
    """(n, n*) of a pops file, each a list of Nlevel lists of Nspace values."""
    b = read_bytes(path)
    _, off = xdr_counted_string(b, 0)
    nl, ns = struct.unpack(">2i", b[off:off + 8])
    v = struct.unpack(f">{2 * nl * ns}d", b[off + 8:off + 8 + 16 * nl * ns])
    rows = [v[i * ns:(i + 1) * ns] for i in range(2 * nl)]
    return rows[:nl], rows[nl:]


def run(name, solve_ne, atoms, files=None, kws=()):
    d = prepare(work, name, "rhf1d", KWS + [("SOLVE_NE", solve_ne)] + list(kws), atoms,
                copy=(), required=("molecules.input", "ray.input"))
    for f, path in (files or {}).items():
        shutil.copy(path, os.path.join(d, f))
    return execute(a.rhf1d, d, timeout=900)


def refused(name, r, text):
    flat = " ".join(r.stdout.split())
    ck.check(r.returncode == EXIT_FATAL and text in flat, f"{name}: stops ({text!r})",
             f"(status {r.returncode})" if r.returncode != EXIT_FATAL else "")


A = run("A", "NONE", [H6, (CA, "ACTIVE", "ZERO_RADIATION")]).require(ok=RAN)
pops_A = os.path.join(A.dir, "pops.CA.out")
old = [H6, (CA, "ACTIVE", "OLD_POPULATIONS", "pops.CA.in")]
A2 = run("A2", "NONE", old, {"pops.CA.in": pops_A}).require(ok=RAN)
B = run("B", "ONCE", old, {"pops.CA.in": pops_A}).require(ok=RAN)
C = run("C", "ONCE", [H6, (CA, "ACTIVE", "ZERO_RADIATION")]).require(ok=RAN)
P = {r.dir: read_pops(os.path.join(r.dir, "pops.CA.out")) for r in (A, A2, B, C)}
nA, sA = P[A.dir]
nA2, sA2 = P[A2.dir]
nB, sB = P[B.dir]
nC, sC = P[C.dir]

# --- 1, 2: n* of the current atmosphere and electron density --     ----

ck.check(sA2 == sA, "A2 (same electron density): n* as in A")
ck.value("A2: n against A", maxrel(nA2, nA), 5e-3)
ck.check(sB == sC, "B (SOLVE_NE = ONCE, read from A): n* of C, computed for this run")
ck.value("B: n* differs from A's (the electron density changed)", maxrel(sB, sA), 0.05,
         at_least=True)
ck.value("B: n against C (same problem, other start)", maxrel(nB, nC), 5e-3)

# --- 3: another atmosphere --                                      ------

with open(os.path.join(a.src, "Atmos", "FALC_82.atmos")) as f:
    atmos = f.read()
for new_id in ("FALC_XX", "FALC_8"):
    d = prepare(work, f"atmos_{new_id}", "rhf1d",
                KWS + [("SOLVE_NE", "NONE"), ("ATMOS_FILE", "falc.atmos")], old,
                copy=(), required=("molecules.input", "ray.input"))
    with open(os.path.join(d, "falc.atmos"), "w") as f:
        f.write(atmos.replace("\n  FALC_82\n", f"\n  {new_id}\n", 1))
    shutil.copy(pops_A, os.path.join(d, "pops.CA.in"))
    refused(f"atmosphere ID {new_id}", execute(a.rhf1d, d, timeout=300),
            f'Populations were computed with atmosphere "FALC_82", not with the current one, '
            f'"{new_id}"')

# --- 4: another model --                                           ------

r = run("he_from_ca", "NONE",
        [H6, (CA, "ACTIVE", "ZERO_RADIATION"),
         ("../../Atoms/He.atom", "PASSIVE", "OLD_POPULATIONS", "pops.CA.in")],
        {"pops.CA.in": pops_A})
refused("Ca II populations for He", r, "Populations have 6 levels at")

ck.finish("pops_file", work)
