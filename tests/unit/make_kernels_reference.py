#!/usr/bin/env python3
"""Reference values for the kernel unit tests (test_interp, test_functions).

All inputs are IEEE doubles (python floats), and everything is computed from
those exact values with mpmath at 50 significant digits and written as a C
header, so that the tests themselves need no Python. Regenerate with

    python3 tests/unit/make_kernels_reference.py > tests/unit/kernels_reference.h

Contents:
  -- E1(x), E2(x)                    (mpmath e1, expint)
  -- Planck function B_nu(lambda, T) with the constants of constant.h
  -- w_k = gamma(k+1, dtau), the incomplete gamma functions that w2/w3
     return, and U_k = exp(-dtau) int_0^dtau t^k exp(t) dt (U3)
  -- the natural cubic spline and the exponential spline (Stoer & Bulirsch
     tension form, uniform tension sigma = tension (N-1) / |x_N-1 - x_0|)
     through a table of exp(sin x) on a non-uniform grid, obtained here by
     solving the continuity conditions of the first derivative directly
     (not by transcribing the tridiagonal recurrence of the C code).
"""
import mpmath as mp

mp.mp.dps = 50

# --- constants of constant.h (CODATA 1986 values, used by Planck) ---
CLIGHT = mp.mpf("2.99792458E+08")
HPLANCK = mp.mpf("6.6260755E-34")
KBOLTZMANN = mp.mpf("1.380658E-23")
NM_TO_M = mp.mpf("1.0E-09")

EXP_X = ([1e-6, 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0,
          1.0000001, 1.1, 1.3, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 7.0, 10.0,
          15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 79.0, 80.0, 81.0,
          100.0, 200.0, 500.0, 600.0])

PLANCK_LAMBDA = [100.0, 500.0, 1600.0, 5000.0, 1e5, 1e7, 1e9]
PLANCK_T = [3000.0, 5000.0, 6000.0, 1e4, 3e4, 1e5, 1e6]

DTAU = [0.0, 1e-9, 1e-6, 1e-5, 1e-4, 1.9e-4, 2e-4, 2.1e-4, 4.9e-4, 5e-4,
        5.1e-4, 1e-3, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 49.0, 50.0,
        51.0, 100.0, 199.0, 200.0, 201.0, 500.0, 1000.0]

SPLINE_N = 17
SPLINE_TENSIONS = [0.5, 3.0]


def c_double(x):
    if x == 0:
        return "0.0"
    if abs(x) < 1e-300:
        return "0.0"
    return mp.nstr(x, 17, min_fixed=1, max_fixed=0)


def array(name, vals, ctype="double"):
    body = "".join("\t%s,\n" % c_double(v) for v in vals)
    return "static const %s %s[%d] = {\n%s};\n" % (ctype, name, len(vals), body)


def spline_grid():
    x = [6.0 * (j / (SPLINE_N - 1.0)) ** 1.3 for j in range(SPLINE_N)]
    y = [float(mp.exp(mp.sin(mp.mpf(xj)))) for xj in x]
    # evaluation points strictly inside, off the nodes: three per interval
    xe = []
    for j in range(SPLINE_N - 1):
        for f in (0.15, 0.5, 0.8):
            xe.append(x[j] + f * (x[j + 1] - x[j]))
    return x, y, xe


def natural_spline(x, y, xe):
    """Natural cubic spline: second derivative M_j, continuity of s'."""
    n = len(x)
    xm = [mp.mpf(v) for v in x]
    ym = [mp.mpf(v) for v in y]
    h = [xm[j + 1] - xm[j] for j in range(n - 1)]
    A = mp.zeros(n - 2, n - 2)
    b = mp.zeros(n - 2, 1)
    for i in range(1, n - 1):
        # h_{i-1} M_{i-1} + 2 (h_{i-1} + h_i) M_i + h_i M_{i+1} = 6 (D_i - D_{i-1})
        r = i - 1
        A[r, r] = 2 * (h[i - 1] + h[i])
        if r > 0:
            A[r, r - 1] = h[i - 1]
        if r < n - 3:
            A[r, r + 1] = h[i]
        b[r] = 6 * ((ym[i + 1] - ym[i]) / h[i] - (ym[i] - ym[i - 1]) / h[i - 1])
    Mi = mp.lu_solve(A, b)
    M = [mp.mpf(0)] + [Mi[k] for k in range(n - 2)] + [mp.mpf(0)]
    out = []
    for xv in xe:
        xv = mp.mpf(xv)
        j = max(k for k in range(n - 1) if xm[k] <= xv)
        t = (xv - xm[j]) / h[j]
        out.append((1 - t) * ym[j] + t * ym[j + 1] +
                   h[j] ** 2 / 6 * ((1 - t) * ((1 - t) ** 2 - 1) * M[j] +
                                    t * (t ** 2 - 1) * M[j + 1]))
    return out


def exp_spline(x, y, xe, tension):
    n = len(x)
    xm = [mp.mpf(v) for v in x]
    ym = [mp.mpf(v) for v in y]
    sigma = mp.mpf(tension) * (n - 1) / abs(xm[-1] - xm[0])
    h = [xm[j + 1] - xm[j] for j in range(n - 1)]
    sh = [mp.sinh(sigma * hj) for hj in h]
    ch = [mp.cosh(sigma * hj) for hj in h]
    # s'(x_i-) from interval i-1 and s'(x_i+) from interval i, M_0 = M_{n-1} = 0
    #   s'_j(x_j)   = D_j + (-sigma ch_j/sh_j + 1/h_j) M_j + (sigma/sh_j - 1/h_j) M_{j+1}
    #   s'_j(x_j+1) = D_j + (-sigma/sh_j + 1/h_j) M_j + (sigma ch_j/sh_j - 1/h_j) M_{j+1}
    A = mp.zeros(n - 2, n - 2)
    b = mp.zeros(n - 2, 1)
    for i in range(1, n - 1):
        r = i - 1
        Dm = (ym[i] - ym[i - 1]) / h[i - 1]
        Dp = (ym[i + 1] - ym[i]) / h[i]
        # left (interval i-1 at its right end), right (interval i at its left end)
        # s'_left = Dm + (-sigma/sh + 1/h) M_{i-1} + (sigma ch/sh - 1/h) M_i
        # s'_right = Dp + (-sigma ch/sh + 1/h) M_i + (sigma/sh - 1/h) M_{i+1}
        cL_prev = -sigma / sh[i - 1] + 1 / h[i - 1]
        cL_self = sigma * ch[i - 1] / sh[i - 1] - 1 / h[i - 1]
        cR_self = -sigma * ch[i] / sh[i] + 1 / h[i]
        cR_next = sigma / sh[i] - 1 / h[i]
        # left - right = 0
        A[r, r] = cL_self - cR_self
        if r > 0:
            A[r, r - 1] = cL_prev
        if r < n - 3:
            A[r, r + 1] = -cR_next
        b[r] = Dp - Dm
    Mi = mp.lu_solve(A, b)
    M = [mp.mpf(0)] + [Mi[k] for k in range(n - 2)] + [mp.mpf(0)]
    out = []
    for xv in xe:
        xv = mp.mpf(xv)
        j = max(k for k in range(n - 1) if xm[k] <= xv)
        t = (xv - xm[j]) / h[j]
        out.append((1 - t) * ym[j] + t * ym[j + 1] +
                   (mp.sinh(sigma * (xm[j + 1] - xv)) / sh[j] - (1 - t)) * M[j] +
                   (mp.sinh(sigma * (xv - xm[j])) / sh[j] - t) * M[j + 1])
    return out


def w_ref(dt):
    dt = mp.mpf(dt)
    return [mp.gammainc(k + 1, 0, dt) for k in range(3)]


def u_ref(dt):
    dt = mp.mpf(dt)
    # exp(-dt) int_0^dt t^k exp(t) dt, by quadrature (independent of the
    # closed forms used by the C code)
    return [mp.exp(-dt) * mp.quad(lambda t, k=k: t ** k * mp.exp(t), [0, dt])
            if dt > 0 else mp.mpf(0) for k in range(3)]


def main():
    print("/* Generated by tests/unit/make_kernels_reference.py -- do not edit."
          "\n   mpmath at 50 digits from the exact double inputs. */\n")

    print("#define EXP_REF_N %d" % len(EXP_X))
    print(array("exp_ref_x", EXP_X))
    print(array("exp_ref_E1", [mp.e1(mp.mpf(v)) for v in EXP_X]))
    print(array("exp_ref_E2", [mp.expint(2, mp.mpf(v)) for v in EXP_X]))

    pl = [(lam, T) for lam in PLANCK_LAMBDA for T in PLANCK_T
          if HPLANCK * CLIGHT / (KBOLTZMANN * NM_TO_M * lam * T) <= 149]
    print("#define PLANCK_REF_N %d" % len(pl))
    print(array("planck_ref_lambda", [p[0] for p in pl]))
    print(array("planck_ref_T", [p[1] for p in pl]))
    vals = []
    for lam, T in pl:
        nu_c = 1 / (NM_TO_M * mp.mpf(lam))   # nu / c
        x = HPLANCK * CLIGHT * nu_c / (KBOLTZMANN * mp.mpf(T))
        vals.append(2 * HPLANCK * CLIGHT * nu_c ** 3 / mp.expm1(x))
    print(array("planck_ref_B", vals))

    print("#define DTAU_REF_N %d" % len(DTAU))
    print(array("dtau_ref", DTAU))
    wr = [w_ref(d) for d in DTAU]
    ur = [u_ref(d) for d in DTAU]
    for k in range(3):
        print(array("w_ref_%d" % k, [w[k] for w in wr]))
    for k in range(3):
        print(array("u_ref_%d" % k, [u[k] for u in ur]))

    x, y, xe = spline_grid()
    print("#define SPLINE_REF_N %d" % len(x))
    print("#define SPLINE_REF_NE %d" % len(xe))
    print(array("spline_ref_x", x))
    print(array("spline_ref_y", y))
    print(array("spline_ref_xe", xe))
    print(array("spline_ref_natural", natural_spline(x, y, xe)))
    print("#define EXPSPLINE_REF_NT %d" % len(SPLINE_TENSIONS))
    print(array("expspline_ref_tension", SPLINE_TENSIONS))
    print("static const double expspline_ref_y[EXPSPLINE_REF_NT][SPLINE_REF_NE] = {")
    for s in SPLINE_TENSIONS:
        print("\t{\n" + "".join("\t\t%s,\n" % c_double(v) for v in exp_spline(x, y, xe, s)) + "\t},")
    print("};")


if __name__ == "__main__":
    main()
