#!/usr/bin/env python3
"""Compare two RH spectrum.out files.

usage: compare_spectra.py A B [--stokes] [--tol 1e-8] [--exact]

spectrum.out layout (writespect_xdr.c), all big-endian XDR:
    int     Nspect
    double  lambda[Nspect]
    double  I[Nspect * Nrays * (Nx * Ny)]
    bool(4 bytes) vacuum_to_air, double vacuum_to_air_limit
    [Stokes_Q, Stokes_U, Stokes_V, each the same size as I, only if a
     magnetic field was read (use --stokes)]

Intensity is compared per wavelength, relative to that wavelength's peak
|I|. Stokes Q, U, V are compared relative to the peak of the whole array,
since they can be zero at many wavelengths.

--exact demands byte-identical files (use for pure memory-management
changes). Otherwise the largest difference must be <= tol. Floating point
sums are not associative, so different thread counts and re-ordered
accumulations legitimately differ in the last bits.

Exit status: 0 = match, 1 = mismatch, 2 = unreadable/incompatible files.
"""
import argparse
import struct
import sys

import numpy as np


def load(path, stokes):
    raw = open(path, "rb").read()
    ns = struct.unpack(">i", raw[:4])[0]
    lam = np.frombuffer(raw, dtype=">f8", count=ns, offset=4)
    body = len(raw) - 4 - 8 * ns - 12
    narr = 4 if stokes else 1
    if ns <= 0 or body <= 0 or body % (8 * narr) or (body // (8 * narr)) % ns:
        raise ValueError(f"{path}: size {len(raw)} does not fit the "
                         f"spectrum.out layout for Nspect={ns}, stokes={stokes}")
    n_i = body // (8 * narr)
    off = 4 + 8 * ns
    arrays = {}
    for k, name in enumerate("IQUV"[:narr]):
        start = off if k == 0 else off + 8 * n_i + 12 + 8 * n_i * (k - 1)
        arrays[name] = np.frombuffer(raw, dtype=">f8", count=n_i,
                                     offset=start).astype(np.float64)
    return raw, lam.astype(np.float64), arrays, ns


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--stokes", action="store_true")
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--exact", action="store_true")
    args = ap.parse_args()

    try:
        raw_a, lam_a, arr_a, ns = load(args.a, args.stokes)
        raw_b, lam_b, arr_b, ns_b = load(args.b, args.stokes)
    except (OSError, ValueError, struct.error) as e:
        print(f"UNREADABLE: {e}")
        return 2

    if ns != ns_b or not np.array_equal(lam_a, lam_b):
        print(f"WAVELENGTH GRID DIFFERS: Nspect {ns} vs {ns_b}")
        return 1
    if len(raw_a) != len(raw_b):
        print(f"SIZE DIFFERS: {len(raw_a)} vs {len(raw_b)} bytes")
        return 1
    if args.exact:
        if raw_a == raw_b:
            print("byte-identical")
            return 0
        n = sum(x != y for x, y in zip(raw_a, raw_b))
        print(f"NOT BYTE-IDENTICAL: {n} differing bytes")
        return 1

    worst = {}
    for name in arr_a:
        x, y = arr_a[name], arr_b[name]
        if not (np.isfinite(x).all() and np.isfinite(y).all()):
            print(f"NON-FINITE VALUES in {name}")
            return 1
        d = np.abs(x - y)
        if name == "I":
            rows = d.reshape(ns, -1).max(axis=1)
            peak = np.abs(x).reshape(ns, -1).max(axis=1)
            rel = np.where(peak > 0, rows / np.where(peak > 0, peak, 1), rows)
            i = int(np.argmax(rel))
            worst[name] = (rel[i], f"lambda={lam_a[i]:.4f}nm (index {i})")
        else:
            peak = np.abs(x).max()
            rel = d.max() / peak if peak > 0 else d.max()
            worst[name] = (rel, "")

    bad = [n for n, (r, _) in worst.items() if r > args.tol]
    summary = "  ".join(f"{n}={r:.1e}" for n, (r, _) in worst.items())
    where = worst["I"][1]
    if bad:
        print(f"MISMATCH (tol {args.tol:.0e}): {summary}  worst I at {where}")
        return 1
    print(f"match (tol {args.tol:.0e}): {summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
