"""Comparison of two spectrum.out files (what tests/compare_spectra.py
prints and returns, without numpy)."""
import math
import struct

from . import output


def _differing_bytes(a, b):
    n, step = 0, 1 << 16
    for i in range(0, len(a), step):
        x, y = a[i:i + step], b[i:i + step]
        if x != y:
            n += sum(u != v for u, v in zip(x, y))
    return n


def compare_spectra(a, b, tol=1e-8, exact=False, stokes=None):
    """(status, message) of the comparison of the spectrum.out files a and b
    (status 0 = match, 1 = mismatch, 2 = unreadable or incompatible files),
    as compare_spectra.py A B [--stokes|--no-stokes] --tol TOL [--exact]
    reports it.

    Intensity is compared per wavelength, relative to that wavelength's peak
    |I|; Stokes Q, U, V relative to the peak of the whole array (they can be
    zero at many wavelengths). exact demands byte-identical files; else the
    largest difference must be <= tol. stokes=True or False states the
    layout and fails if the file has the other one, None finds it."""
    try:
        raw_a, lam_a, arr_a, ns = output.read_spectrum_file(a, stokes)
        raw_b, lam_b, arr_b, ns_b = output.read_spectrum_file(b, stokes)
    except (OSError, ValueError, struct.error) as e:
        return 2, f"UNREADABLE: {e}"

    if ns != ns_b or lam_a != lam_b:
        return 1, f"WAVELENGTH GRID DIFFERS: Nspect {ns} vs {ns_b}"
    if len(raw_a) != len(raw_b):
        return 1, f"SIZE DIFFERS: {len(raw_a)} vs {len(raw_b)} bytes"
    if exact:
        if raw_a == raw_b:
            return 0, "byte-identical"
        return 1, f"NOT BYTE-IDENTICAL: {_differing_bytes(raw_a, raw_b)} differing bytes"

    worst = {}
    for name in arr_a:
        x, y = arr_a[name], arr_b[name]
        if not (all(map(math.isfinite, x)) and all(map(math.isfinite, y))):
            return 1, f"NON-FINITE VALUES in {name}"
        if name == "I":
            n = len(x) // ns
            rel = []
            for i in range(ns):
                xs, ys = x[i * n:(i + 1) * n], y[i * n:(i + 1) * n]
                row = max(abs(u - v) for u, v in zip(xs, ys))
                peak = max(abs(u) for u in xs)
                rel.append(row / peak if peak > 0 else row)
            i = max(range(ns), key=rel.__getitem__)       # the first maximum
            worst[name] = (rel[i], f"lambda={lam_a[i]:.4f}nm (index {i})")
        else:
            d = max(abs(u - v) for u, v in zip(x, y))
            peak = max(abs(u) for u in x)
            worst[name] = (d / peak if peak > 0 else d, "")

    bad = [n for n, (r, _) in worst.items() if r > tol]
    summary = "  ".join(f"{n}={r:.1e}" for n, (r, _) in worst.items())
    where = worst["I"][1]
    if bad:
        return 1, f"MISMATCH (tol {tol:.0e}): {summary}  worst I at {where}"
    return 0, f"match (tol {tol:.0e}): {summary}"
