#!/usr/bin/env python3
"""Build Atoms/Hminus_bf.ecsv and the fit constants of Hminus_bf (N-35).

The H- bound-free cross-section in RH has three parts, all derived from the
merged and scaled single-detachment cross-section of McLaughlin, Stancil,
Sadeghpour & Forrey (2017, J. Phys. B 50, 114001), supplementary file
JPBaa6c1fSuppdata1.dat (photon energy [eV], cross-section [Mb]):

  lambda < lambda_min of the table:
      sigma = A E^p (1 + (Ec/E)^q)^(-s)          (E in eV, sigma in Mb)
      (A, p) fitted to the file at 50-1000 eV, where the file is a power law;
      (Ec, q, s) fitted from 85 nm to 45 eV plus all E >= 50 eV, with
      25.5-28.5 nm left out (a feature of the R-matrix data).
  table, lambda_min <= lambda <= lambda_max (80-125 nm):
      linear interpolation in non-uniform knots, placed greedily so that the
      interpolation stays within TOL of the 0.1 nm bin averages (empty bins
      bridged); the end knots take the values of the fits.
  lambda > lambda_max of the table, up to the threshold lambda0:
      sigma = lambda^3 x^(3/2) sum_n C_n x^(n/2)  (lambda in micron, Mb),
      x = 1/lambda - 1/lambda0, the form of John (1988, A&A 193, 189),
      refitted to the file at lambda > 120 nm.

Usage:
    python3 hminus_bf_table.py McL_FILE OUT.ecsv

Needs numpy, scipy and astropy. The McLaughlin file itself is not part of RH:
the authors retain its copyright. Do not commit it.
"""
import sys

import numpy as np
from astropy.table import Table
import astropy.units as u
from scipy.optimize import least_squares

HC = 1239.84198           # h c [eV nm]
E0 = 0.754204             # threshold adopted by McLaughlin et al. (Drake 1988) [eV]
L0 = HC / E0 / 1000       # threshold wavelength [micron]

N_JOHN = 6                # terms in the John form
LAM_JOHN_FIT = 120.0      # [nm] John form fitted above this
LAM_SHORT_FIT = 85.0      # [nm] short-wavelength form fitted below this
E_RMPS_MAX = 45.0         # [eV] ... from LAM_SHORT_FIT up to this energy,
E_POWER = (50.0, 1000.0)  # [eV] ... plus the power-law range
EXCLUDE = (25.5, 28.5)    # [nm] left out of the short-wavelength fit
TAB_LO, TAB_HI = 80.0, 125.0   # [nm] table range
BIN = 0.1                 # [nm] bin width for the table averages
TOL = 0.01                # relative tolerance of the knot placement


def load(path):
    E, s = np.loadtxt(path, comments='#', unpack=True)
    o = np.argsort(E)                  # one row is out of order
    return E[o], s[o]


def bin_loglam(E, s, lmin, lmax, nbins):
    """Means over log-lambda bins: equal weight per d ln lambda in the fits."""
    lam = HC / E
    keep = (lam >= lmin) & (E > E0 + 1e-4)
    edges = np.geomspace(lmin, lmax, nbins)
    idx = np.digitize(lam, edges)
    out_l, out_s = [], []
    for k in range(1, len(edges)):
        m = keep & (idx == k)
        if m.any():
            out_l.append(np.exp(np.log(lam[m]).mean()))
            out_s.append(s[m].mean())
    return np.array(out_l), np.array(out_s)


def fit_powerlaw(E, s):
    m = (E >= E_POWER[0]) & (E <= E_POWER[1])
    p, lnA = np.polyfit(np.log(E[m]), np.log(s[m]), 1)
    return np.exp(lnA), p


def sigma_short(lam_nm, A, p, Ec, q, s):
    E = HC / np.asarray(lam_nm, float)
    return A * E**p * (1 + (Ec / E)**q)**(-s)


def fit_short(E, s, A, p):
    lam, sig = bin_loglam(E, s, HC / 1000, 95.0, 600)
    Eb = HC / lam
    m = ((Eb >= HC / LAM_SHORT_FIT) & (Eb <= E_RMPS_MAX)) | (Eb >= E_POWER[0])
    m &= ~((lam > EXCLUDE[0]) & (lam < EXCLUDE[1]))
    r = least_squares(lambda x: np.log(sigma_short(lam[m], A, p, *x) / sig[m]),
                      [30.0, 2.0, 1.0])
    return r.x


def sigma_long(lam_nm, C):
    lu = np.asarray(lam_nm, float) / 1e3
    x = np.clip(1 / lu - 1 / L0, 0, None)
    return lu**3 * x**1.5 * sum(c * x**(n / 2) for n, c in enumerate(C))


def fit_long(E, s):
    lam, sig = bin_loglam(E, s, 19.0, L0 * 1e3, 1500)
    k = np.sqrt(1e3 / lam - 1 / L0)
    f = sig / ((lam / 1e3)**3 * k**3)
    m = lam > LAM_JOHN_FIT
    A = np.array([k[m]**n for n in range(N_JOHN)]).T
    w = 1 / f[m]
    c, *_ = np.linalg.lstsq(A * w[:, None], f[m] * w, rcond=None)
    return c


def bin_average(E, s, lo, hi, width):
    lam = HC / E
    edges = np.arange(lo, hi + width / 2, width)
    idx = np.digitize(lam, edges)
    vals = np.array([s[idx == k].mean() if np.any(idx == k) else np.nan
                     for k in range(1, len(edges))])
    return 0.5 * (edges[1:] + edges[:-1]), vals


def greedy_knots(x, y, tol):
    """Fewest knots whose linear interpolation stays within tol of all y."""
    knots, i = [0], 0
    while i < len(x) - 1:
        j = i + 1
        while j + 1 < len(x):
            jj = j + 1
            yi = np.interp(x[i:jj + 1], [x[i], x[jj]], [y[i], y[jj]])
            if np.max(np.abs(yi / y[i:jj + 1] - 1)) > tol:
                break
            j = jj
        knots.append(j)
        i = j
    return np.array(knots)


def main(mcl_file, out):
    E, s = load(mcl_file)
    A, p = fit_powerlaw(E, s)
    Ec, q, ss = fit_short(E, s, A, p)
    C = fit_long(E, s)
    short = (A, p, Ec, q, ss)

    xc, v = bin_average(E, s, TAB_LO, TAB_HI, BIN)
    ok = ~np.isnan(v)
    ok[[0, -1]] = False       # bin centres next to the pinned ends
    x = np.r_[TAB_LO, xc[ok], TAB_HI]
    y = np.r_[sigma_short(TAB_LO, *short), v[ok], sigma_long(TAB_HI, C)]
    kn = greedy_knots(x, y, TOL)
    lam_t, sig_t = x[kn], y[kn]

    t = Table()
    t['lambda'] = np.round(lam_t, 4) * u.nm
    t['sigma'] = np.array([float(f'{v:.6g}') for v in sig_t]) * u.Mbarn
    t['lambda'].description = 'vacuum wavelength'
    t['sigma'].description = 'H- bound-free (photodetachment) cross-section'
    t.meta['title'] = 'H- bound-free cross-section, 80-125 nm, for RH (Hminus_bf)'
    t.meta['source'] = ('McLaughlin B M, Stancil P C, Sadeghpour H R and Forrey R C 2017, '
                        'J. Phys. B: At. Mol. Opt. Phys. 50, 114001, supplementary data '
                        'JPBaa6c1fSuppdata1.dat (merged and scaled single detachment)')
    t.meta['doi'] = '10.1088/1361-6455/aa6c1f'
    t.meta['derivation'] = (f'Means of the source over {BIN} nm bins in wavelength (empty bins '
                            f'bridged linearly); knots placed so that linear interpolation '
                            f'stays within {TOL:.0%} of the bin means; the end knots take the '
                            'values of the fits RH uses outside the table.')
    t.meta['generator'] = 'docs/opacity/hminus_bf_table.py'
    t.meta['rh_reader'] = ('RH reads only this form: space-delimited, the columns lambda [nm] '
                           'and sigma [Mbarn] in this order, lambda strictly increasing. '
                           'It applies its fits below the first and above the last knot; '
                           'the table must cover 85-120 nm.')
    t.write(out, format='ascii.ecsv', overwrite=True)

    cstr = lambda a: ', '.join(f'{v:.10g}' for v in a)
    print(f'table: {len(lam_t)} knots, {lam_t[0]}-{lam_t[-1]} nm, '
          f'spacing {np.diff(lam_t).min():.2f}-{np.diff(lam_t).max():.2f} nm -> {out}')
    print(f'short: A = {A:.10g}, p = {p:.10g}, Ec = {Ec:.10g}, q = {q:.10g}, s = {ss:.10g}')
    print(f'long:  lambda0 = {L0:.10g} micron, C = {{{cstr(C)}}}')
    lam_c, sig_c = bin_loglam(E, s, 19.0, L0 * 1e3, 1500)
    for lo, hi, fn in ((30, 80, lambda l: sigma_short(l, *short)),
                       (125, 1500, lambda l: sigma_long(l, C))):
        m = (lam_c >= lo) & (lam_c < hi)
        print(f'  max |fit/McL - 1| {lo}-{hi} nm: {np.abs(fn(lam_c[m]) / sig_c[m] - 1).max():.4f}')


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
