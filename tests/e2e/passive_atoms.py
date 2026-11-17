#!/usr/bin/env python3
"""End-to-end test of the initial solutions of PASSIVE atoms in atoms.input (B-59).

usage: passive_atoms.py --rhf1d PATH [--src DIR] [--work DIR] [-v]

A PASSIVE atom keeps the populations it starts with: LTE (LTE_POPULATIONS),
or populations read from a file (OLD_POPULATIONS). Passive hydrogen takes its
populations from the atmosphere (LTE with HYDROGEN_LTE). Before B-59 RH
accepted input that it then ignored, without a word: ZERO_RADIATION gave
LTE, a populations file without OLD_POPULATIONS was not read, and passive
hydrogen read OLD_POPULATIONS and threw them away. Each now stops at input
with a message and status 70:

 1. PASSIVE atom with ZERO_RADIATION, without and with a populations file.
 2. PASSIVE atom with LTE_POPULATIONS and a populations file.
 3. PASSIVE hydrogen with OLD_POPULATIONS.

What stays valid:

 4. PASSIVE atom with LTE_POPULATIONS (no file): runs.
 5. ACTIVE atom with LTE_POPULATIONS and a file (read with STARTING_J =
    OLD_J): runs.
 6. PASSIVE atom with OLD_POPULATIONS: He solved ACTIVE writes pops.HE.out;
    He PASSIVE reads it ("NLTE populations for atom HE"), and the spectrum
    differs from He PASSIVE in LTE (the He I and He II continua below
    50.4 nm), so the populations read are the ones used.

Python 3 standard library only. Exit status 0 = pass.
"""
import os
import shutil
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.output import read_bytes
from rhtest.run import EXIT_FATAL, RAN, execute, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d",), verbose=True)
work = WorkDir(a.src, a.work)
ck = Checker(a.verbose)

H6 = ("../../Atoms/H_6.atom", "PASSIVE", "LTE_POPULATIONS")
CAII = ("../../Atoms/CaII.atom", "ACTIVE", "ZERO_RADIATION", "pops.CA.out")
HE = "../../Atoms/He.atom"
KWS = [("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"), ("N_THREADS", "1"),
       ("PRINT_CPU", "FALSE"), ("STARTING_J", "NEW_J"), ("KURUCZ_DATA", "none"),
       ("N_MAX_ITER", "2"), ("SOLVE_NE", "NONE")]


def run(name, atoms, files=None):
    """rhf1d on the rhf1d example with the atoms of atoms.input; files
    {name: path} are copied into the run directory first."""
    d = prepare(work, name, "rhf1d", KWS, atoms, copy=(),
                required=("molecules.input", "ray.input"))
    for f, path in (files or {}).items():
        shutil.copy(path, os.path.join(d, f))
    return execute(a.rhf1d, d, timeout=600)


def refused(name, atoms, text, files=None):
    r = run(name, atoms, files)
    flat = " ".join(r.stdout.split())
    ck.check(r.returncode == EXIT_FATAL and text in flat,
             f"{name}: stops at input ({text!r})",
             f"(status {r.returncode})" if r.returncode != EXIT_FATAL else "")


# --- 1-3: refused --                                         ------------

refused("passive_zero_radiation", [H6, CAII, (HE, "PASSIVE", "ZERO_RADIATION")],
        "ZERO_RADIATION is for ACTIVE atoms only (PASSIVE atom HE")
refused("passive_zero_radiation_file", [H6, CAII, (HE, "PASSIVE", "ZERO_RADIATION", "pops.HE.out")],
        "ZERO_RADIATION is for ACTIVE atoms only (PASSIVE atom HE")
refused("passive_lte_file", [H6, CAII, (HE, "PASSIVE", "LTE_POPULATIONS", "pops.HE.out")],
        "Populations file pops.HE.out of PASSIVE atom HE would not be read")

# --- 4-6: valid --                                           ------------

lte = run("passive_lte", [H6, CAII, (HE, "PASSIVE", "LTE_POPULATIONS")])
ck.check(lte.ran, "passive He with LTE_POPULATIONS runs", f"(status {lte.returncode})")

act = run("active_lte_file", [H6, ("../../Atoms/CaII.atom", "ACTIVE", "LTE_POPULATIONS", "pops.CA.out")])
ck.check(act.ran, "active Ca II with LTE_POPULATIONS and a file runs", f"(status {act.returncode})")

he = run("active_he", [H6, CAII, (HE, "ACTIVE", "ZERO_RADIATION", "pops.HE.out")])
pops = os.path.join(he.dir, "pops.HE.out")
ck.check(he.ran and os.path.exists(pops), "active He writes pops.HE.out", f"(status {he.returncode})")

if os.path.exists(pops):
    # passive hydrogen cannot read populations, whatever the file
    refused("passive_h_old", [("../../Atoms/H_6.atom", "PASSIVE", "OLD_POPULATIONS", "pops.HE.out"),
                              CAII], "PASSIVE hydrogen cannot read OLD_POPULATIONS",
            files={"pops.HE.out": pops})
    old = run("passive_old", [H6, CAII, (HE, "PASSIVE", "OLD_POPULATIONS", "pops.HE.out")],
              files={"pops.HE.out": pops})
    ck.check(old.ran and "NLTE populations for atom HE" in old.stdout,
             "passive He reads OLD_POPULATIONS", f"(status {old.returncode})")
    if old.ran and lte.ran:
        s_old = read_bytes(os.path.join(old.dir, "spectrum.out"))
        s_lte = read_bytes(os.path.join(lte.dir, "spectrum.out"))
        ck.check(len(s_old) == len(s_lte) and s_old != s_lte,
                 "the populations read change the spectrum (not LTE)")

ck.finish("passive_atoms", work)
