#!/usr/bin/env python3
"""Compare two RH spectrum.out files.

usage: compare_spectra.py A B [--stokes | --no-stokes] [--tol 1e-8] [--exact]

spectrum.out layout (writespect_xdr.c), all big-endian XDR:
    int     Nspect
    double  lambda[Nspect]
    double  I[Nspect * Nrays * (Nx * Ny)]
    bool(4 bytes) vacuum_to_air, double vacuum_to_air_limit
    [Stokes_Q, Stokes_U, Stokes_V, each the same size as I, only if a
     magnetic field was read or BACKGROUND_POLARIZATION = TRUE]

I is stored with the wavelength as the slowest index, I[Nspect][rays...].
The file size does not tell the two layouts apart: Q, U and V read as part
of I still make a multiple of Nspect values (Q-21). The layout is found
from the vacuum_to_air record instead, which must follow I: a bool (0 or
1) and VACUUM_TO_AIR_LIMIT (spectrum.h, 200.0 here and upstream).
--stokes or --no-stokes states the layout and fails if the file has the
other one.

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
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rhtest import compare, output      # noqa: E402


def load(path, stokes=None):
    """Read spectrum.out (rhtest.output.read_spectrum_file finds the layout;
    stokes=True or False demands one), for analyze_run.py: the raw bytes,
    lambda, the arrays I (and Q, U, V) as numpy arrays and Nspect."""
    import numpy as np
    raw, lam, arrays, ns = output.read_spectrum_file(path, stokes)
    return (raw, np.frombuffer(lam, dtype=np.float64),
            {k: np.frombuffer(v, dtype=np.float64) for k, v in arrays.items()}, ns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    pol = ap.add_mutually_exclusive_group()
    pol.add_argument("--stokes", dest="stokes", action="store_true",
                     default=None)
    pol.add_argument("--no-stokes", dest="stokes", action="store_false")
    ap.add_argument("--tol", type=float, default=1e-8)
    ap.add_argument("--exact", action="store_true")
    args = ap.parse_args()
    # The comparison is rhtest.compare.compare_spectra, which run_tests.py and
    # run_stress.py use as well
    status, message = compare.compare_spectra(args.a, args.b, tol=args.tol, exact=args.exact,
                                              stokes=args.stokes)
    print(message)
    return status


if __name__ == "__main__":
    sys.exit(main())
