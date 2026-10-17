ids: Q-14, B-5, N-17, N-18
date: 2026-10-17
title: 1-D formal solvers tested on analytic slabs; Bezier and linear weights accurate to rounding

- `test_formal_1d` runs the 1-D solvers on slabs whose source function is a
  polynomial in optical depth, where the intensities in both directions are
  known in closed form. It checks exactness for a linear S, the order of
  convergence as the grid is refined (constant and exponentially varying
  opacity), the diagonal operator against dI/dS, and the Bezier weights
  against mpmath. Measured at 16 points per decade:

  | Solver | order, emergent I | error, emergent I | Psi |
  |---|---|---|---|
  | BEZIER3 (default) | 4 | 9e-7 | approximate, 28% from dI/dS (N-19) |
  | LINEAR | 2 | 1e-3 | exact |
  | FEAUTRIER | 2 | 4e-4 | exact |
  | FEAUTRIER_HERMITE | 4 (2 with varying opacity) | 4e-7 | not implemented |

- B-5: the linear short-characteristics solver computed the last point of
  every ray with the weights of the previous interval (`w2` not called
  again). Its emergent intensity was off by up to 1e-5 relative for a
  linear S, which it should reproduce exactly; now 1e-15.
- N-17: the cubic Bezier weights used a Taylor series to dt³ below
  dt = 0.05 (up to 1.2e-5 relative error) and closed forms that lose
  digits to cancellation above. Now an 18-term series below dt = 1 and the
  closed forms above: within 6.5e-15 of mpmath for dt from 1e-8 to 1e3.
  Bezier3 is the default solver in all geometries, so this changes most
  results by up to about 1e-6.
- N-18: `w2` and `w3` sum the incomplete gamma series (all terms positive)
  below dtau = 1 instead of a two-term series below 5e-4 and the closed
  forms above.
