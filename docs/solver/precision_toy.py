#!/usr/bin/env python3
"""Toy model for SOLVER.md section 3.7: precision of G and Newton-Krylov.

Fixed-point problem x = G(x), G(x) = M x + b + c sin(x), N = 200, with the
slowest mode of M at 0.995 (kappa(I - M) ~ 200) and c = 0.002, so that G is
a contraction with a unique fixed point. Inexact Newton on R(x) = G(x) - x
with finite-difference Jacobian products and restarted GMRES(40).

The variants differ in the precision of the R used for the Newton
right-hand side and of the R evaluations inside the Jacobian products:

  dbl  everything double
  qbg  the x-dependent term computed in double and rounded to float, as
       RH's background opacities are when Background() runs inside G
  flt  everything float32

Printed: the true residual ||R(x)|| (double) after every second Newton step,
and the number of R evaluations (the unit of work).

Usage: OMP_NUM_THREADS=1 python3 precision_toy.py   (numpy, scipy; ~5 s)
"""
import inspect

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres

rng = np.random.default_rng(1)
N = 200
Q, _ = np.linalg.qr(rng.standard_normal((N, N)))
lam = rng.uniform(0, 0.9, N)
lam[0] = 0.995
M = (Q * lam) @ Q.T
b = rng.standard_normal(N) * 0.1
c = 0.002       # rho(M) + c < 1: contraction, unique fixed point

# scipy renamed gmres(tol=...) to rtol in 1.12
TOL = "rtol" if "rtol" in inspect.signature(gmres).parameters else "tol"


def G(x, mode):
    if mode == "dbl":
        return M @ x + b + c * np.sin(x)
    if mode == "qbg":
        return M @ x + b + (c * np.sin(x)).astype(np.float32).astype(np.float64)
    if mode == "flt":
        x32 = x.astype(np.float32)
        return (M.astype(np.float32) @ x32 + b.astype(np.float32)
                + np.float32(c) * np.sin(x32)).astype(np.float64)
    raise ValueError(mode)


def R(x, mode):
    return G(x, mode) - x


def jfnk(rhs_mode, jv_mode, eps_rel, nnewton=12, eta=1e-3):
    x = np.zeros(N)
    nR = 0
    hist = []
    for _ in range(nnewton):
        r = R(x, rhs_mode)
        nR += 1
        hist.append(np.linalg.norm(R(x, "dbl")))    # true residual, not counted
        h = eps_rel * (1 + np.linalg.norm(x))
        r0 = R(x, jv_mode)                          # base of the differences
        nR += 1

        def mv(v):
            nonlocal nR
            nv = np.linalg.norm(v)
            if nv == 0:
                return np.zeros_like(v)
            nR += 1
            return (R(x + h * v / nv, jv_mode) - r0) / h * nv

        A = LinearOperator((N, N), matvec=mv)
        d, _ = gmres(A, -r, restart=40, maxiter=20, atol=0, **{TOL: eta})
        x = x + d
    hist.append(np.linalg.norm(R(x, "dbl")))
    return hist, nR


CASES = [
    ("all double, eps=1.5e-8", "dbl", "dbl", 1.5e-8),
    ("float bg, eps=1.5e-8", "qbg", "qbg", 1.5e-8),
    ("float bg, eps=1e-4 (from the noise)", "qbg", "qbg", 1e-4),
    ("R double, Jv all float, eps=1e-3", "dbl", "flt", 1e-3),
    ("all float, eps=1e-3", "flt", "flt", 1e-3),
]

for eta, nnewton in ((1e-3, 12), (1e-1, 16)):
    print(f"forcing term eta = {eta:g}: true ||R|| after Newton steps 0, 2, 4, ...; R evaluations")
    for name, rm, jm, e in CASES:
        h, n = jfnk(rm, jm, e, nnewton=nnewton, eta=eta)
        print(f"  {name:38s}", " ".join(f"{v:8.1e}" for v in h[::2]), f"  nR={n}")
