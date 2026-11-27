#!/usr/bin/env python3
"""End-to-end test of tools/waveinfo (B-75, B-76).

usage: waveinfo.py --waveinfo PATH [--src DIR] [--work DIR] [-v]

A wavelength file (XDR) holds the number of wavelengths, an int, and then
the wavelengths, doubles in nm. waveinfo never checked its reads: an empty
file gave "Found 0 wavelengths, ranging from 0.000 to 0.000 [nm]" from an
uninitialized count, reading outside the array, and a file with a count of
5 and two values gave "Found 5 wavelengths". A file that did not open
ended with status 0 and a file that was read with status 1. Now:

 1. A file of 7 wavelengths: their number and range, status 0.
 2. A count of 0: "Found no wavelengths", status 0.
 3. An empty file: "no number of wavelengths", status 1.
 4. A count of 5 and two values: "it ends before its 5 wavelengths",
    status 1.
 5. A negative count: "a negative number of wavelengths", status 1.
 6. A file that does not exist: "Error opening file", status 1.

With -full, waveinfo lists the wavelengths five to a row. The last row
was always the last five, so a number that is not a multiple of 5 repeated
up to four of them (B-76: 7 listed as 1-5, 3-7). Now each comes once:

 7. -full on 7, 10 and 12 wavelengths: rows of five, the last one shorter.

Python 3 standard library only. Exit status 0 = pass.
"""
import os
import struct
import subprocess
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rhtest.check import Checker
from rhtest.cli import parse
from rhtest.workdir import WorkDir

a = parse(required=("waveinfo",), verbose=True)
work = WorkDir(a.src, a.work)
ck = Checker(a.verbose)


def wavefile(name, count, values):
    path = work.join(name)
    with open(path, "wb") as f:
        if count is not None:
            f.write(struct.pack(">i", count))
        f.write(struct.pack(f">{len(values)}d", *values))
    return path


def waveinfo(path, status, text):
    r = subprocess.run([a.waveinfo, path], capture_output=True, text=True, timeout=60)
    out = " ".join((r.stdout + r.stderr).split())
    ck.check(r.returncode == status and text in out,
             f"{os.path.basename(path)}: status {status}, {text!r}",
             f"(status {r.returncode}: {out!r})")


SEVEN = [121.567, 279.635, 393.366, 396.847, 588.995, 656.280, 854.209]
waveinfo(wavefile("seven.wave", 7, SEVEN), 0,
         "Found 7 wavelengths, ranging from 121.567 to 854.209 [nm] (in vacuum)")
waveinfo(wavefile("zero.wave", 0, []), 0, "Found no wavelengths")
waveinfo(wavefile("empty.wave", None, []), 1, "no number of wavelengths")
waveinfo(wavefile("short.wave", 5, SEVEN[:2]), 1, "it ends before its 5 wavelengths")
waveinfo(wavefile("negative.wave", -3, []), 1, "a negative number of wavelengths (-3)")
waveinfo(work.join("no_such.wave"), 1, "Error opening file")


def listing(n):
    """-full on n wavelengths 100, 101, ...: stdout must be the rows of
    five, each wavelength once, then an empty line."""
    values = [100.0 + k for k in range(n)]
    path = wavefile(f"full{n}.wave", n, values)
    r = subprocess.run([a.waveinfo, "-full", path], capture_output=True, text=True,
                       timeout=60)
    rows = ["   ".join(f"{v:9.3f}" for v in values[k:k + 5]) for k in range(0, n, 5)]
    want = "\n".join(rows) + "\n\n"
    ck.check(r.returncode == 0 and r.stdout == want,
             f"-full on {n} wavelengths: rows of five, each once",
             f"(status {r.returncode}: {r.stdout!r})")


for n in (7, 10, 12):
    listing(n)

ck.finish("waveinfo", work)
