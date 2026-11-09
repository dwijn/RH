#!/usr/bin/env python3
"""Summarize RH runs for the stress / convergence suite.

usage: analyze_run.py summary RUNDIR ITER_LIMIT
       analyze_run.py pops RUNDIR_A RUNDIR_B
       analyze_run.py compare RUNDIR REFDIR

summary prints one JSON line with
    iters       number of main iterations performed
    delta       max population change in the last iteration (all atoms)
    converged   delta < ITER_LIMIT
    trend       'monotone', 'oscillating' or 'stalled' over the last iterations
    wall, rss   wall-clock seconds and peak resident memory [MB]
    status      'ok', or the first error found in the log
    work        the work counters RH prints at the end (SOLVER.md 4.1):
                nspect, wavelength solves and sweeps per phase (initial,
                iterate, final) and kind (operator, plain, PRD), fs = all
                wavelength solves / nspect (formal-solution equivalents),
                calls and seconds of Background() and Hydrostatic(),
                seconds in Iterate(), threads
    fse         SE residual metrics over all active atoms
                (se_residual.py; needs RADRATE_OUTPUT and COLLRATE_OUTPUT)

pops prints the largest relative difference of the NLTE populations of every
pops.*.out file present in both run dirs, using only entries above 1e-10 of
the level-summed population at that point, so that tiny populations of high
levels do not dominate.

compare prints one JSON line with the error of RUNDIR against a reference
solution REFDIR (same case, S-1 in docs/SOLVER.md):
    pops        max relative population error over all active atoms,
                restricted to levels with n / n_tot >= 1e-10, 1e-6, 1e-3
                (n_tot = sum of the LTE populations at that point)
    spec_wl     max over wavelengths of max|dI| / max|I_ref|, both over
                rays and columns at that wavelength (compare_spectra.py);
                Stokes I only when the runs also wrote Q, U and V
    spec_line   max over active lines of max|dI| / max|I_ref| over the
                line's window: the line's Nlambda wavelengths of the
                spectrum nearest to lambda0 (an approximation of its own
                grid, which RH does not write)
    worst_line  wavelength of the line that sets spec_line
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
    out = dict(iters=iters, delta=delta,
               converged=bool(delta < iter_limit) if not math.isnan(delta) else False,
               trend=trend, wall=wall, rss=rss, status=status,
               history=[float(f"{h:.3e}") for h in hist])
    w = work(log)
    if w:
        out["work"] = w
    try:
        from se_residual import run_metrics
        fse = run_metrics(rundir)
        if fse:
            out["fse"] = fse["all"]
    except (OSError, ValueError, SystemExit) as e:
        out["fse"] = f"unreadable: {e}"
    print(json.dumps(out))


def work(log):
    """The work block printed by printWork() (work.c), or None."""
    m = re.search(r"-- Work \(wavelength solves; one formal solution = (\d+)\)\n"
                  r"((?:\s+\w+(?:\s+\w+ \d+ \(\d+ sweeps\))+\n)+)"
                  r"\s+Background (\d+) calls ([\d.]+) s;\s+Hydrostatic (\d+) calls ([\d.]+) s;"
                  r"\s+Iterate ([\d.]+) s;\s+threads (\d+)", log)
    if not m:
        return None
    w = dict(nspect=int(m.group(1)))
    total = 0
    for line in m.group(2).splitlines():
        f = line.split()
        if not f:
            continue
        w[f[0]] = {k: int(v) for k, v in re.findall(r"(\w+) (\d+) \(", line)}
        w[f[0]].update({k + "_sweeps": int(v) for k, v in re.findall(r"(\w+) \d+ \((\d+) sweeps", line)})
        total += sum(int(v) for _, v in re.findall(r"(\w+) (\d+) \(", line))
    w["fs"] = total / w["nspect"]
    w.update(background_calls=int(m.group(3)), background_s=float(m.group(4)),
             hydrostatic_calls=int(m.group(5)), hydrostatic_s=float(m.group(6)),
             iterate_s=float(m.group(7)), threads=int(m.group(8)))
    return w


def read_pops(path, lte=False):
    """NLTE populations n[level, point] of a pops.*.out file, or the LTE
    populations n* (written after them) with lte=True."""
    b = open(path, "rb").read()
    for off in range(4, min(len(b) - 8, 512), 4):
        nl, ns = struct.unpack(">2i", b[off:off + 8])
        if 0 < nl < 100000 and 0 < ns and off + 8 + 16 * nl * ns == len(b):
            return np.frombuffer(b, ">f8", nl * ns,
                                 off + 8 + (8 * nl * ns if lte else 0)).reshape(nl, ns)
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


def compare(run, ref):
    from compare_spectra import load
    from se_residual import THRESHOLDS, read_atom
    out = {"pops": {}}
    worst = {t: 0.0 for t in THRESHOLDS}
    lines = []
    for fr in sorted(glob.glob(os.path.join(ref, "pops.*.out"))):
        f = os.path.join(run, os.path.basename(fr))
        if not os.path.exists(f):
            out["pops"] = f"missing {os.path.basename(fr)}"; break
        nr, n = read_pops(fr), read_pops(f)
        if nr.shape != n.shape:
            out["pops"] = f"shape mismatch {os.path.basename(fr)}"; break
        frac = nr / read_pops(fr, lte=True).sum(axis=0)
        rel = np.abs(n - nr) / np.maximum(np.abs(nr), 1e-300)
        for t in THRESHOLDS:
            sel = frac >= t
            if sel.any():
                worst[t] = max(worst[t], float(rel[sel].max()))
        ID = os.path.basename(fr)[5:-4]
        fa = os.path.join(ref, f"atom.{ID}.out")
        if os.path.exists(fa):
            lines += read_atom(fa)[4]
    else:
        out["pops"] = {f"{t:.0e}": v for t, v in worst.items()}
    try:
        _, lam, Ir, ns = load(os.path.join(ref, "spectrum.out"))
        _, lam2, I, ns2 = load(os.path.join(run, "spectrum.out"))
    except (OSError, ValueError) as e:
        out["spec"] = f"unreadable: {e}"
        print(json.dumps(out)); return
    if ns != ns2 or Ir["I"].size != I["I"].size or not np.array_equal(lam, lam2):
        out["spec"] = "wavelength grids differ"
        print(json.dumps(out)); return
    # I[Nspect][rays, columns] (initial_xdr.c): wavelength slowest (Q-22)
    Ir, I = Ir["I"].reshape(ns, -1), I["I"].reshape(ns, -1)
    d = np.abs(I - Ir).max(axis=1)
    pk = np.abs(Ir).max(axis=1)
    out["spec_wl"] = float((d / np.where(pk > 0, pk, 1)).max())
    best, where = 0.0, None
    for lambda0, Nlambda in lines:
        i = np.argsort(np.abs(lam - lambda0))[:Nlambda]
        e = float(d[i].max() / max(pk[i].max(), 1e-300))
        if e > best:
            best, where = e, lambda0
    if lines:
        out["spec_line"], out["worst_line"] = best, where
    print(json.dumps(out))


if __name__ == "__main__":
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    if len(sys.argv) == 4 and sys.argv[1] == "summary":
        summary(sys.argv[2], float(sys.argv[3]))
    elif len(sys.argv) == 4 and sys.argv[1] == "pops":
        pops(sys.argv[2], sys.argv[3])
    elif len(sys.argv) == 4 and sys.argv[1] == "compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
