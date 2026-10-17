#!/usr/bin/env python3
"""Coefficients of Weideman's (1994) rational approximation of the
Faddeeva function, as used in voigt.c:

    w(z) = 2 p(Z) / (L - iz)^2 + 1 / (sqrt(pi) (L - iz)),
    Z = (L + iz) / (L - iz),  L = sqrt(N / sqrt(2)),

where p is the polynomial of degree N - 1 whose coefficients are the
Fourier coefficients of exp(-t^2) (L^2 + t^2) on t = L tan(theta / 2)
(J. A. C. Weideman, SIAM J. Numer. Anal. 31, 1497, 1994). Computed with
mpmath at 40 digits; prints the C initializer, highest power first.

    python3 tests/unit/make_weideman_coefficients.py [N]   (default 32)
"""
import sys

import mpmath

mpmath.mp.dps = 40


def coefficients(N):
    M = 2 * N
    L = mpmath.sqrt(N / mpmath.sqrt(2))
    f = [mpmath.mpf(0)]
    for k in range(-M + 1, M):
        t = L * mpmath.tan(k * mpmath.pi / (2 * M))
        f.append(mpmath.exp(-t * t) * (L * L + t * t))
    n = len(f)
    fs = f[n // 2:] + f[:n // 2]          # fftshift
    a = [mpmath.fsum(fs[m] * mpmath.cos(2 * mpmath.pi * j * m / n)
                     for m in range(n)) / (2 * M) for j in range(n)]
    return L, a[1:N + 1][::-1]


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 32
    L, a = coefficients(N)
    print("#define WEIDEMAN_N %d" % N)
    print("#define WEIDEMAN_L %s\n" % mpmath.nstr(L, 21))
    print("static const double weideman_a[WEIDEMAN_N] = {")
    for i in range(0, N, 3):
        print("\t" + ", ".join(mpmath.nstr(x, 17, min_fixed=1, max_fixed=0)
                               for x in a[i:i + 3]) + ",")
    print("};")


if __name__ == "__main__":
    main()
