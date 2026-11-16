#!/usr/bin/env python3
"""H- opacity checks for N-35 (bound-free) and N-45 (free-free) (numpy only).

Evidence: claude/RH/evidence/st1-n35-t2 (README) on the HAO file system.

Reproduces, from the coefficients as printed in John (1988, A&A 193, 189):
  bf: Geltman (1962) table (natural cubic spline, as in hydrogen.c) against the
      John Table 2 fit; behaviour of the fit below 125 nm; optional comparison
      with McLaughlin et al. (2017) supplementary data (argument: path to
      JPBaa6c1fSuppdata1.dat).
  ff: Stilley & Callaway (1970) table (hydrogen.c) against John Tables 3a/3b;
      3b below 182.3 nm; 3a below theta = 0.5; ff/bf ratio.
  Bell & Berrington (1987) Table 1 at theta = 0.5, 1, 2, 3.6 (read from the paper).
  John Table 1 (total) reproduced from Eq. (3)-(6): its bound-free part is 10x too high.
  LTE order-of-magnitude estimate of what dominates below 125 nm.
Usage: python3 -I hminus_check.py [JPBaa6c1fSuppdata1.dat]
"""
import sys
import numpy as np

HC_EVNM = 1239.84198          # h c [eV nm]
K_EV = 8.617333e-5            # k [eV/K]

# ---------------- bound-free ----------------
LAM_G = np.array([0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700,
                  750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350,
                  1400, 1450, 1500, 1550, 1600, 1641.9])
SIG_G = np.array([0, .15, .33, .57, .85, 1.17, 1.52, 1.89, 2.23, 2.55, 2.84, 3.11, 3.35, 3.56,
                  3.71, 3.83, 3.92, 3.95, 3.93, 3.85, 3.73, 3.58, 3.38, 3.14, 2.85, 2.54, 2.20,
                  1.83, 1.46, 1.06, .71, .40, .17, 0]) * 10.0      # 1e-21 m^2 -> Mb


def natural_spline(x, y):
    n = len(x)
    h = np.diff(x)
    A = np.zeros((n, n)); r = np.zeros(n)
    A[0, 0] = A[-1, -1] = 1
    for i in range(1, n - 1):
        A[i, i - 1:i + 2] = h[i - 1], 2 * (h[i - 1] + h[i]), h[i]
        r[i] = 6 * ((y[i + 1] - y[i]) / h[i] - (y[i] - y[i - 1]) / h[i - 1])
    M = np.linalg.solve(A, r)

    def f(t):
        t = np.atleast_1d(np.asarray(t, float))
        i = np.clip(np.searchsorted(x, t) - 1, 0, n - 2)
        a, b = x[i + 1] - t, t - x[i]
        return (M[i] * a**3 + M[i + 1] * b**3) / (6 * h[i]) + \
            (y[i] / h[i] - M[i] * h[i] / 6) * a + (y[i + 1] / h[i] - M[i + 1] * h[i] / 6) * b
    return f


geltman = natural_spline(LAM_G, SIG_G)
C_BF = [152.519, 49.534, -118.858, 92.536, -34.194, 4.982]   # John Table 2
LAM0 = 1.6419                                                # micron


def john_bf(lam_nm):
    """sigma [Mb] = lam^3 x^1.5 sum C_n x^((n-1)/2), lam in micron."""
    l = np.asarray(lam_nm, float) / 1000
    x = 1 / l - 1 / LAM0
    return l**3 * x**1.5 * sum(c * x**(n / 2) for n, c in enumerate(C_BF))


# ---------------- free-free ----------------
LAM_SC = np.array([0.0, 303.8, 455.6, 506.3, 569.5, 650.9, 759.4, 911.3, 1013.0, 1139.0,
                   1302.0, 1519.0, 1823.0, 2278.0, 3038.0, 4556.0, 9113.0])
TH_SC = np.arange(0.5, 2.01, 0.1)
K_SC = np.array([
    [0] * 16,
    [3.44e-02, 4.18e-02, 4.91e-02, 5.65e-02, 6.39e-02, 7.13e-02, 7.87e-02, 8.62e-02, 9.36e-02, 1.01e-01, 1.08e-01, 1.16e-01, 1.23e-01, 1.30e-01, 1.38e-01, 1.45e-01],
    [7.80e-02, 9.41e-02, 1.10e-01, 1.25e-01, 1.40e-01, 1.56e-01, 1.71e-01, 1.86e-01, 2.01e-01, 2.16e-01, 2.31e-01, 2.45e-01, 2.60e-01, 2.75e-01, 2.89e-01, 3.03e-01],
    [9.59e-02, 1.16e-01, 1.35e-01, 1.53e-01, 1.72e-01, 1.90e-01, 2.08e-01, 2.25e-01, 2.43e-01, 2.61e-01, 2.78e-01, 2.96e-01, 3.13e-01, 3.30e-01, 3.47e-01, 3.64e-01],
    [1.21e-01, 1.45e-01, 1.69e-01, 1.92e-01, 2.14e-01, 2.36e-01, 2.58e-01, 2.80e-01, 3.01e-01, 3.22e-01, 3.43e-01, 3.64e-01, 3.85e-01, 4.06e-01, 4.26e-01, 4.46e-01],
    [1.56e-01, 1.88e-01, 2.18e-01, 2.47e-01, 2.76e-01, 3.03e-01, 3.31e-01, 3.57e-01, 3.84e-01, 4.10e-01, 4.36e-01, 4.62e-01, 4.87e-01, 5.12e-01, 5.37e-01, 5.62e-01],
    [2.10e-01, 2.53e-01, 2.93e-01, 3.32e-01, 3.69e-01, 4.06e-01, 4.41e-01, 4.75e-01, 5.09e-01, 5.43e-01, 5.76e-01, 6.08e-01, 6.40e-01, 6.72e-01, 7.03e-01, 7.34e-01],
    [2.98e-01, 3.59e-01, 4.16e-01, 4.70e-01, 5.22e-01, 5.73e-01, 6.21e-01, 6.68e-01, 7.15e-01, 7.60e-01, 8.04e-01, 8.47e-01, 8.90e-01, 9.32e-01, 9.73e-01, 1.01e+00],
    [3.65e-01, 4.39e-01, 5.09e-01, 5.75e-01, 6.39e-01, 7.00e-01, 7.58e-01, 8.15e-01, 8.71e-01, 9.25e-01, 9.77e-01, 1.03e+00, 1.08e+00, 1.13e+00, 1.18e+00, 1.23e+00],
    [4.58e-01, 5.50e-01, 6.37e-01, 7.21e-01, 8.00e-01, 8.76e-01, 9.49e-01, 1.02e+00, 1.09e+00, 1.15e+00, 1.22e+00, 1.28e+00, 1.34e+00, 1.40e+00, 1.46e+00, 1.52e+00],
    [5.92e-01, 7.11e-01, 8.24e-01, 9.31e-01, 1.03e+00, 1.13e+00, 1.23e+00, 1.32e+00, 1.40e+00, 1.49e+00, 1.57e+00, 1.65e+00, 1.73e+00, 1.80e+00, 1.88e+00, 1.95e+00],
    [7.98e-01, 9.58e-01, 1.11e+00, 1.25e+00, 1.39e+00, 1.52e+00, 1.65e+00, 1.77e+00, 1.89e+00, 2.00e+00, 2.11e+00, 2.21e+00, 2.32e+00, 2.42e+00, 2.51e+00, 2.61e+00],
    [1.14e+00, 1.36e+00, 1.58e+00, 1.78e+00, 1.98e+00, 2.17e+00, 2.34e+00, 2.52e+00, 2.68e+00, 2.84e+00, 3.00e+00, 3.15e+00, 3.29e+00, 3.43e+00, 3.57e+00, 3.70e+00],
    [1.77e+00, 2.11e+00, 2.44e+00, 2.75e+00, 3.05e+00, 3.34e+00, 3.62e+00, 3.89e+00, 4.14e+00, 4.39e+00, 4.63e+00, 4.86e+00, 5.08e+00, 5.30e+00, 5.51e+00, 5.71e+00],
    [3.10e+00, 3.71e+00, 4.29e+00, 4.84e+00, 5.37e+00, 5.87e+00, 6.36e+00, 6.83e+00, 7.28e+00, 7.72e+00, 8.14e+00, 8.55e+00, 8.95e+00, 9.33e+00, 9.71e+00, 1.01e+01],
    [6.92e+00, 8.27e+00, 9.56e+00, 1.08e+01, 1.19e+01, 1.31e+01, 1.42e+01, 1.52e+01, 1.62e+01, 1.72e+01, 1.82e+01, 1.91e+01, 2.00e+01, 2.09e+01, 2.17e+01, 2.25e+01],
    [2.75e+01, 3.29e+01, 3.80e+01, 4.28e+01, 4.75e+01, 5.19e+01, 5.62e+01, 6.04e+01, 6.45e+01, 6.84e+01, 7.23e+01, 7.60e+01, 7.97e+01, 8.32e+01, 8.67e+01, 9.01e+01]])


def rh_ff(lam_nm, theta):
    """What hydrogen.c does: John 3a (no clamp) from 9113 nm, else sc_ff."""
    if lam_nm >= 9113.0:
        return john_ff_table(lam_nm, theta, T3A)
    return sc_ff(lam_nm, theta)


def sc_ff(lam_nm, theta):
    """S&C as in hydrogen.c: bilinear, theta clamped; units 1e-26 cm^4/dyn."""
    th = np.clip(theta, 0.5, 2.0)
    i = np.clip(np.searchsorted(LAM_SC, lam_nm) - 1, 0, len(LAM_SC) - 2)
    t = (lam_nm - LAM_SC[i]) / (LAM_SC[i + 1] - LAM_SC[i])
    row = (1 - t) * K_SC[i] + t * K_SC[i + 1]
    return np.interp(th, TH_SC, row)


# John Table 3a (lambda > 0.3645 um) and 3b (0.1823-0.3645 um); rows n=1..6, cols A..F
T3A = np.array([[0, 0, 0, 0, 0, 0],
                [2483.3460, 285.8270, -2054.2910, 2827.7760, -1341.5370, 208.9520],
                [-3449.8890, -1158.3820, 8746.5230, -11485.6320, 5303.6090, -812.9390],
                [2200.0400, 2427.7190, -13651.1050, 16755.5240, -7510.4940, 1132.7380],
                [-696.2710, -1841.4000, 8624.9700, -10051.5300, 4400.0670, -655.0200],
                [88.2830, 444.5170, -1863.8640, 2095.2880, -901.7880, 132.9850]])
T3B = np.array([[518.1021, -734.8666, 1021.1775, -479.0721, 93.1373, -6.4285],
                [473.2636, 1443.4137, -1977.3395, 922.3575, -178.9275, 12.3600],
                [-482.2089, -737.1616, 1096.8827, -521.1341, 101.7963, -7.0571],
                [115.5291, 169.6374, -245.6490, 114.2430, -21.9972, 1.5097],
                [0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])


def john_ff_table(lam_nm, theta, T):
    l = lam_nm / 1000
    pw = np.array([l**2, 1, 1 / l, 1 / l**2, 1 / l**3, 1 / l**4])
    n = np.arange(1, 7)
    # Eq. (6): 1e-29 cm^4/dyn -> report in 1e-26 cm^4/dyn
    return 1e-3 * np.sum(theta**((n + 1) / 2) * (T @ pw))


def john_ff(lam_nm, theta, tail=True):
    """3a above 364.5 nm, 3b from 182.3 nm, k ~ lambda^2 below (if tail)."""
    if lam_nm >= 364.5:
        return john_ff_table(lam_nm, theta, T3A)
    if lam_nm >= 182.3 or not tail:
        return john_ff_table(lam_nm, theta, T3B)
    return john_ff_table(182.3, theta, T3B) * (lam_nm / 182.3)**2


def john_bf_kappa(lam_nm, theta):
    """Eq. (3): k_bf per H atom per unit Pe, in 1e-26 cm^4/dyn."""
    T = 5040.0 / theta
    lA = lam_nm * 10
    a = 1.439e8
    sig = john_bf(lam_nm) * 1e-18 if lam_nm < 1641.9 else 0.0
    return 0.750 * T**-2.5 * np.exp(a / (LAM0 * 1e4) / T) * (1 - np.exp(-a / (lA * T))) * sig / 1e-26


def hdr(s):
    print('\n=== ' + s)


def main():
    hdr('bf: Geltman spline vs John Table 2 [Mb]')
    for l in [30, 50, 100, 106, 110, 125, 150, 170, 200, 300, 500, 800, 1000, 1400, 1500, 1600, 1620, 1635]:
        g, j = geltman(l)[0], john_bf(l)
        print(f'  {l:6.0f} nm  Geltman {g:7.3f}  John {j:7.3f}  G/J {g / j:6.3f}')
    l = np.linspace(60, 160, 100001)
    s = john_bf(l)
    k = np.argmin(s)
    print(f'  John minimum: {s[k]:.3f} Mb at {l[k]:.1f} nm')
    e = 1e-4
    slope = (np.log(john_bf(125 * (1 + e))) - np.log(john_bf(125 * (1 - e)))) / (2 * e)
    print(f'  d ln sigma / d ln lambda at 125 nm: {slope:+.3f}')
    print('  lambda*sigma [nm Mb] at 0.01, 1, 10 nm: '
          + ', '.join(f'{x * john_bf(x):.0f}' for x in (0.01, 1, 10))
          + f'  -> 1000 C6 = {C_BF[5] * 1000:.0f} (sigma -> C6/lambda[um] as lambda -> 0)')
    print(f'  step at 125 nm: John/Geltman = {john_bf(125) / geltman(125)[0]:.3f}, '
          f'1 - G/J = {1 - geltman(125)[0] / john_bf(125):.3f}')

    hdr('John Table 1 (total, 1e-26 cm^4/dyn) reproduced from Tables 2, 3a, 3b')
    T1 = {0.2: [1.676, 3.135, 9.077, 22.41, 50.00, 104.0, 205.6, 390.8, 720.0],
          0.3: [2.850, 5.355, 15.54, 38.38, 85.65, 178.2, 352.2, 669.3, 1233.],
          0.5: [4.925, 9.423, 27.83, 69.12, 154.5, 321.7, 635.9, 1208., 2226.],
          1.0: [5.400, 10.64, 33.11, 85.18, 194.5, 410.0, 816.4, 1558., 2877.],
          1.6: [0.998, 1.321, 2.318, 4.244, 8.038, 15.37, 29.17, 54.38, 99.32],
          2.0: [1.335, 1.610, 2.124, 2.610, 3.029, 3.466, 3.862, 4.237, 4.584],
          3.0: [2.966, 3.562, 4.685, 5.720, 6.694, 7.605, 8.463, 9.276, 10.06]}
    ths = [0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    worst = 0
    for lu, row in T1.items():
        r = [(john_bf_kappa(lu * 1000, th) + john_ff(lu * 1000, th)) / v for th, v in zip(ths, row)]
        worst = max(worst, max(abs(np.array(r) - 1)))
        print(f'  {lu:4.1f} um  computed/Table1: ' + ' '.join(f'{x:.4f}' for x in r))
    print(f'  max |ratio-1| = {worst:.2e}')

    print('  NOTE: Table 1 is reproduced by 10*k_bf + k_ff (see below); Eq. (3) constant 0.750 is the')
    print('  exact Saha factor, and k_ff matches Bell & Berrington, so Table 1 bf is 10x too high (misprint).')
    for lu, row in T1.items():
        if lu > 1.6:
            continue
        r = [(10 * john_bf_kappa(lu * 1000, th) + john_ff(lu * 1000, th)) / v for th, v in zip(ths, row)]
        print(f'  {lu:4.1f} um  (10 bf + ff)/Table1: ' + ' '.join(f'{x:.4f}' for x in r))
    me, h, kB = 9.1093837e-28, 6.62607015e-27, 1.380649e-16
    print(f'  Saha constant (h^2/2 pi m)^1.5 k^-2.5 / 4 = {(h**2 / (2 * np.pi * me))**1.5 * kB**-2.5 / 4:.4f}')

    hdr('ff: John 3a/3b vs Bell & Berrington (1987) Table 1 (read from the scan), 1e-26 cm^4/dyn')
    BB_LAM = [151890, 113918, 91134, 45567, 30378, 22784, 18227, 15189, 13019, 11392, 10126, 9113,
              7595, 6510, 5696, 5063, 4557, 3645, 3038, 2604, 2278, 1823]
    BB = {0.5: [75.1, 42.3, 27.0, 6.79, 3.04, 1.73, 1.11, .781, .579, .448, .358, .293, .208, .154,
                .121, .0965, .0791, .0520, .0364, .0277, .0228, .0178],
          1.0: [144, 80.8, 51.9, 13.1, 5.86, 3.37, 2.17, 1.52, 1.13, .871, .702, .576, .409, .309,
                .241, .195, .161, .108, .0789, .0615, .0499, .0402],
          2.0: [253, 142, 91.4, 23.0, 10.3, 5.95, 3.80, 2.69, 2.02, 1.57, 1.28, 1.06, .774, .597,
                .475, .390, .327, .225, .169, .136, .110, .0896],
          3.6: [388, 219, 140, 35.0, 15.6, 9.06, 5.97, 4.31, 3.31, 2.60, 2.17, 1.83, 1.36, 1.07,
                .861, .711, .602, .418, .318, .262, .211, .172]}
    for th, row in BB.items():
        r = np.array([john_ff(l / 10, th) / v for l, v in zip(BB_LAM, row)])
        s = np.array([rh_ff(l / 10, th) / v for l, v in zip(BB_LAM, row)])
        print(f'  theta {th}: John/BB {r.min():.3f}..{r.max():.3f} (worst at {BB_LAM[np.argmax(abs(r - 1))]} A); '
              f'S&C(RH)/BB {s.min():.3f}..{s.max():.3f} (worst at {BB_LAM[np.argmax(abs(s - 1))]} A)')

    print('  RH jump at 9113 nm (3a / S&C just below), theta 0.5, 1, 2, 2.5, 3, 3.6: ' + ' '.join(
        f'{john_ff_table(9113.0, th, T3A) / sc_ff(9112.999, th):.3f}' for th in (0.5, 1.0, 2.0, 2.5, 3.0, 3.6)))

    hdr('ff: S&C (hydrogen.c) / John 3a-3b, theta 0.5..2.0')
    for l in [303.8, 364.5, 455.6, 506.3, 650.9, 911.3, 1823.0, 4556.0, 9113.0]:
        r = [sc_ff(l, th) / john_ff(l, th) for th in TH_SC]
        print(f'  {l:7.1f} nm  min {min(r):.3f} max {max(r):.3f}')
    print('  3a/3b at 364.5 nm: ' + ' '.join(
        f'{john_ff_table(364.5, th, T3A) / john_ff_table(364.5, th, T3B):.4f}' for th in (0.5, 1.0, 2.0, 3.6)))
    print('  S&C-path at 182.3, 250 nm (linear to 0) / John 3b: ' + ', '.join(
        f'{l}: {min(sc_ff(l, th) / john_ff(l, th) for th in TH_SC):.2f}-'
        f'{max(sc_ff(l, th) / john_ff(l, th) for th in TH_SC):.2f}' for l in (182.3, 250.0)))
    print('  S&C theta clamp: theta=3.0 (1680 K) S&C uses theta=2: S&C/John at 911, 4556 nm: '
          + ', '.join(f'{sc_ff(l, 3.0) / john_ff(l, 3.0):.3f}' for l in (911.3, 4556.0)))

    hdr('ff: Table 3b below 182.3 nm (no tail): k/lambda^2 [1e-26 cm^4/dyn / um^2]')
    for th in (0.5, 1.0, 2.0):
        vals = [(l, john_ff_table(l, th, T3B) / (l / 1000)**2) for l in (182.3, 160, 150, 140, 120, 100)]
        print(f'  theta {th}: ' + ' '.join(f'{l:.0f}:{v:.3g}' for l, v in vals))
    l = np.linspace(100, 182.3, 8231)
    for th in (0.5, 0.6, 0.8):
        v = np.array([john_ff_table(x, th, T3B) for x in l])
        neg = l[v < 0]
        print(f'  theta {th}: 3b negative below {neg.max():.1f} nm' if neg.size else f'  theta {th}: positive down to 100 nm')

    hdr('ff/bf (John, per H atom per Pe), 150-300 nm')
    for l in (150, 182.3, 200, 250, 300):
        print(f'  {l:6.1f} nm: ' + ' '.join(f'th{th}:{john_ff(l, th) / john_bf_kappa(l, th):.3f}' for th in (0.6, 0.8, 1.0, 1.2)))

    hdr('ff: Table 3a below theta = 0.5 (extrapolation), lambda = 1 um, 10 um, 1 mm')
    for l in (1000.0, 1e4, 1e6):
        v = [john_ff_table(l, th, T3A) / (l / 1000)**2 for th in (0.5, 0.4, 0.3, 0.2, 0.15, 0.1)]
        print(f'  {l:9.0f} nm k/lambda^2: ' + ' '.join(f'{x:.4g}' for x in v))
    th = np.linspace(0.05, 0.5, 451)
    for l in (1000.0, 1e4, 1e6):
        v = np.array([john_ff_table(l, t, T3A) for t in th])
        print(f'  {l:9.0f} nm: positive for theta >= {th[np.argmax(v > 0)] if (v > 0).any() else np.nan:.2f}; '
              f'monotonic in theta: {bool(np.all(np.diff(v) > 0))}')

    hdr('LTE order of magnitude, T = 6400 K, ne = 7e13 cm^-3')
    T, ne = 6400.0, 7e13
    me, h, kB = 9.1093837e-28, 6.62607015e-27, 1.380649e-16
    lam_th = (h**2 / (2 * np.pi * me * kB * T))**1.5
    hm = ne / 4 * lam_th * np.exp(0.754204 / (K_EV * T))
    print(f'  n(H-)/n(H I) = {hm:.3g}')
    # C I/C II: U(CI)~9.3 (3P, 1D), U(CII)~5.9; chi = 11.2603 eV
    ratio_cii = 2 * 5.9 / 9.3 / ne / lam_th * np.exp(-11.2603 / (K_EV * T))
    fCI = 1 / (1 + ratio_cii)
    AC = 2.7e-4
    f3P = 9 / 9.3
    f1D = 5 / 9.3 * np.exp(-1.2637 / (K_EV * T))
    print(f'  C I fraction {fCI:.3f}; ground 3P share {f3P:.3f}, 1D share {f1D:.3g}')
    sigC = 1e-17
    sH = john_bf(106)
    print(f'  C I 3P bf / H- bf (sigma_C = 1e-17, sigma_H- = {sH:.2f} Mb at 106 nm): '
          f'{AC * fCI * f3P * sigC / (hm * sH * 1e-18):.3g}')
    sH = john_bf(120)
    print(f'  C I 1D bf / H- bf at 120 nm (sigma_C = 1e-17 assumed): '
          f'{AC * fCI * f1D * sigC / (hm * sH * 1e-18):.3g}')
    # Rayleigh, Lyman-alpha wing, classical oscillator, f = 0.4162
    w = 121.567 / 125.0
    sR = 6.652e-25 * 0.4162**2 * w**4 / (1 - w**2)**2
    print(f'  H Rayleigh (one-oscillator Lya) at 125 nm: {sR:.3g} cm^2; / H- bf: {sR / (hm * john_bf(125) * 1e-18):.3g}')

    if len(sys.argv) > 1:
        hdr('McLaughlin et al. 2017 (supplementary data 1)')
        E, sig = np.loadtxt(sys.argv[1], comments='#', unpack=True)
        lam = HC_EVNM / E
        print(f'  {len(E)} points, E {E.min():.6g}-{E.max():.6g} eV (lambda {lam.min():.3g}-{lam.max():.6g} nm); '
              f'monotonic E: {bool(np.all(np.diff(E) > 0))}')
        dE = np.diff(E)
        for a, b in ((0.75, 1), (1, 3), (3, 9), (9, 10.9), (10.9, 14.4), (14.4, 50), (50, 1e5)):
            m = (E[:-1] >= a) & (E[:-1] < b)
            if m.any():
                print(f'    {a:g}-{b:g} eV: {m.sum()} pts, dE {dE[m].min():.3g}-{dE[m].max():.3g} eV')
        o = np.argsort(lam)
        lo, so = lam[o], sig[o]

        def mcl(x):
            return np.interp(x, lo, so)

        def binavg(a, b):
            m = (lo >= a) & (lo <= b)
            return np.trapezoid(so[m], lo[m]) / (lo[m][-1] - lo[m][0])
        print('  pointwise McLaughlin / John and / Geltman:')
        for l in (125, 130, 150, 170, 200, 250, 300, 400, 500, 800, 1000, 1200, 1400, 1500, 1600, 1620, 1635):
            print(f'    {l:5d} nm  McL {mcl(l):7.3f}  /John {mcl(l) / john_bf(l):6.3f}  /Geltman {mcl(l) / geltman(l)[0]:6.3f}')
        ll = np.linspace(125, 1630, 20000)
        r = mcl(ll) / john_bf(ll)
        for a, b in ((125, 170), (170, 300), (300, 1500), (1500, 1630), (170, 1630)):
            m = (ll >= a) & (ll <= b)
            print(f'    McL/John {a}-{b} nm: {r[m].min():.3f} .. {r[m].max():.3f}')
        print('  below 125 nm, bin averages over 5-nm bins [Mb]: McL, Geltman, Geltman*1.24, John')
        for a in (30, 50, 70, 85, 90, 95, 100, 105, 110, 115, 120):
            b = a + 5
            c = 0.5 * (a + b)
            print(f'    {a:3d}-{b:3d} nm: McL {binavg(a, b):7.3f}  G {geltman(c)[0]:6.3f}  '
                  f'G*1.24 {1.24 * geltman(c)[0]:6.3f}  John {john_bf(c):7.3f}')
        print(f'  jump at 125 nm if McL below / John above: McL(125)/John(125) = {mcl(125) / john_bf(125):.3f}')
        m = (lo > 86) & (lo < 114)
        print(f'  resonance region 86-114 nm: McL min {so[m].min():.3g} max {so[m].max():.3g} Mb')


if __name__ == '__main__':
    main()
