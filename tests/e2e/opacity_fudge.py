#!/usr/bin/env python3
"""End-to-end test of the opacity fudge table check (B-65).

usage: opacity_fudge.py --rhf1d PATH [--src DIR] [--work DIR] [-v]

OPACITY_FUDGE names a table of wavelength and three values, each a factor
minus 1 for the H-minus, scattering and metal background opacities,
interpolated with Linear. It was read without a check:

 1. one wavelength: Linear read one element past the table;
 2. a value below -1: a negative factor, a negative opacity;
 3. wavelengths out of order or repeated: Linear takes the direction from
    the first two, and divides by zero on a repeat;
 4. a line with fewer than four numbers: the rest of the table uninitialized.

Each now stops at input with a message and status 70. The shipped table
(Atmos/opacity_fudge.input) still runs.

Python 3 standard library only. Exit status 0 = pass.
"""
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.run import EXIT_FATAL, EXIT_OK, execute, prepare
from rhtest.workdir import WorkDir

a = parse(required=("rhf1d",), verbose=True)
work = WorkDir(a.src, a.work)
ck = Checker(a.verbose)

KWS = [("STOKES_MODE", "NO_STOKES"), ("STOKES_INPUT", "none"), ("N_THREADS", "1"),
       ("PRINT_CPU", "FALSE"), ("KURUCZ_DATA", "none"), ("N_MAX_ITER", "1"),
       ("OPACITY_FUDGE", "fudge.input")]
ATOMS = [("../../Atoms/H_6.atom", "PASSIVE", "LTE_POPULATIONS")]


def case(name, table, text=None):
    """Run with fudge.input = table; text = None: must run, else must stop with it."""
    d = prepare(work, name, "rhf1d", KWS, ATOMS, copy=("molecules.input", "ray.input"))
    with open(os.path.join(d, "fudge.input"), "w") as f:
        f.write(table)
    r = execute(a.rhf1d, d, timeout=300)
    flat = " ".join(r.stdout.split())
    if text is None:
        ck.check(r.returncode in (EXIT_OK, 1) and "Opacity fudge file" not in flat,
                 f"{name}: runs", f"(status {r.returncode})")
    else:
        ck.check(r.returncode == EXIT_FATAL and text in flat, f"{name}: stops ({text!r})",
                 f"(status {r.returncode})" if r.returncode != EXIT_FATAL else "")


with open(os.path.join(a.src, "Atmos", "opacity_fudge.input")) as f:
    case("shipped", f.read())

case("one_wavelength", "1\n500.0  0.0 0.0 0.0\n", "the number of wavelengths, at least 2")
case("below_minus_one", "2\n300.0  0.0 -1.5 0.0\n800.0  0.0 0.0 0.0\n",
     "value -1.5 below -1")
case("out_of_order", "3\n300.0 0 0 0\n500.0 0 0 0\n400.0 0 0 0\n",
     "wavelengths must be strictly increasing or decreasing (400 nm after 500 nm)")
case("repeated", "2\n300.0 0 0 0\n300.0 0 0 0\n",
     "wavelengths must be strictly increasing or decreasing (300 nm repeated)")
case("short_line", "2\n300.0 0 0\n800.0 0 0 0\n", "wavelength 1 of 2: need four numbers")
case("zero_factor", "2\n300.0 -1.0 0.0 0.0\n800.0 -1.0 0.0 0.0\n")

ck.finish("opacity_fudge", work)
