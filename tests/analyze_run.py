#!/usr/bin/env python3
"""Summarize RH runs for the stress / convergence suite.

usage: analyze_run.py summary RUNDIR ITER_LIMIT
       analyze_run.py pops RUNDIR_A RUNDIR_B

summary prints one JSON line with
    iters       number of main iterations performed
    delta       max population change in the last iteration (all atoms)
    converged   delta < ITER_LIMIT
    trend       'monotone', 'oscillating' or 'stalled' over the last iterations
    wall, rss   wall-clock seconds and peak resident memory [MB]
    status      'ok', or the first error found in the log

pops prints the largest relative difference of the NLTE populations of every
pops.*.out file present in both run dirs, using only entries above 1e-10 of
the level-summed population at that point, so that tiny populations of high
levels do not dominate.
"""
import glob
import json
import math
import os
import re
import struct
import sys

import numpy as np


def summary(rundir, iter_limit):
    log = open(os.path.join(rundir, "log.log"), errors="replace").read()
    status = "ok"
    m = re.search(r"(TERMINATING_ERROR|AddressSanitizer|runtime error|"
                  r"Segmentation fault)[^\n]*\n?[^\n]*", log)
    if m:
        status = " ".join(m.group(0).split())[:160]
    elif re.search(r"(^|[^A-Za-z])-?nan([^A-Za-z]|$)", log, re.I):
        status = "NaN in log"
    # Per iteration, the largest of the per-atom/molecule deltas
    hist = []
    for block in re.split(r"\n -- Iteration\s+\d+\n", log)[1:]:
        d = [float(x) for x in re.findall(
            r"^ ?\S+\s*, delta = ([-+0-9.Ee]+|nan|NaN)", block, re.M)]
        if d:
            hist.append(max(d))
    iters = len(hist)
    delta = hist[-1] if hist else float("nan")
    tail = hist[-6:]
    trend = "n/a"
    if len(tail) >= 4:
        ups = sum(b > a for a, b in zip(tail, tail[1:]))
        ratio = tail[-1] / tail[0] if tail[0] > 0 else float("nan")
        if ups >= 2:
            trend = "oscillating"
        elif ratio > 0.5:
            trend = "stalled"
        else:
            trend = "monotone"
    wall = rss = None
    tf = os.path.join(rundir, "time.txt")
    if os.path.exists(tf):
        t = open(tf).read()
        m = re.search(r"Elapsed \(wall clock\) time.*: (?:(\d+):)?(\d+):([\d.]+)", t)
        if m:
            wall = int(m.group(1) or 0) * 3600 + int(m.group(2)) * 60 + float(m.group(3))
        m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
        if m:
            rss = int(m.group(1)) / 1024
        # Exit status 1-63 means the run completed (bits, e.g. 1 = not
        # converged, which "converged" already reports); 64 and up is a
        # failure, as is a signal
        m = re.search(r"Exit status: (\d+)", t)
        if (re.search(r"Command terminated by signal", t) or
                (m and int(m.group(1)) >= 64)) and status == "ok":
            status = "nonzero exit: " + " ".join(re.findall(r"(?:signal \d+|Exit status: \d+)", t))
    print(json.dumps(dict(iters=iters, delta=delta,
                          converged=bool(delta < iter_limit) if not math.isnan(delta) else False,
                          trend=trend, wall=wall, rss=rss, status=status,
                          history=[float(f"{h:.3e}") for h in hist])))


def read_pops(path):
    b = open(path, "rb").read()
    for off in range(4, min(len(b) - 8, 512), 4):
        nl, ns = struct.unpack(">2i", b[off:off + 8])
        if 0 < nl < 100000 and 0 < ns and off + 8 + 16 * nl * ns == len(b):
            return np.frombuffer(b, ">f8", nl * ns, off + 8).reshape(nl, ns)
    raise ValueError(f"cannot parse {path}")


def pops(a, b):
    out = {}
    for fa in sorted(glob.glob(os.path.join(a, "pops.*.out"))):
        fb = os.path.join(b, os.path.basename(fa))
        if not os.path.exists(fb):
            continue
        try:
            na, nb = read_pops(fa), read_pops(fb)
        except ValueError as e:
            out[os.path.basename(fa)] = str(e); continue
        if na.shape != nb.shape:
            out[os.path.basename(fa)] = "shape mismatch"; continue
        ref = np.abs(na).sum(axis=0)
        mask = np.abs(na) > 1e-10 * ref
        rel = np.abs(na - nb) / np.maximum(np.abs(na), 1e-300)
        out[os.path.basename(fa)] = float(rel[mask].max()) if mask.any() else 0.0
    print(json.dumps(out))


if __name__ == "__main__":
    if sys.argv[1] == "summary":
        summary(sys.argv[2], float(sys.argv[3]))
    elif sys.argv[1] == "pops":
        pops(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
