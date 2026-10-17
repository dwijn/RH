ids: Q-14, B-5, N-17, N-18
date: 2026-10-17
title: 1-D formal solvers tested on analytic slabs; accurate Bezier and linear weights

- `test_formal_1d` runs the 1-D solvers on slabs whose source function is a
  polynomial in optical depth, where the intensities in both directions are
  known in closed form. It checks exactness for a linear S, the order of
  convergence as the grid is refined (constant and exponentially varying
  opacity), the diagonal operator against dI/dS, and the Bezier weights
  against mpmath. Measured at 16 points per decade:

  | Solver | order, emergent I | error, emergent I | Psi |
  |---|---|---|---|
  | BEZIER3 (default; 1-D only with velocities or Stokes, see N-23) | 4 | 9e-7 | approximate, 28% from dI/dS (N-19) |
  | LINEAR | 2 | 1e-3 | exact |
  | FEAUTRIER | 2 | 4e-4 | exact |
  | FEAUTRIER_HERMITE | 4 (2 with varying opacity) | 4e-7 | not implemented |

- B-5: the linear short-characteristics solver computed the last point of
  every ray with the weights of the previous interval (`w2` not called
  again). Its emergent intensity was off by up to 1e-5 relative for a
  linear S, which it should reproduce exactly; now 1e-15.
- N-17: the cubic Bezier weights used a Taylor series to dt³ below
  dt = 0.05 (up to 1.2e-5 relative error) and closed forms that lose
  digits to cancellation above. Now the series below |dt| = 0.2, with as
  many terms as a relative error of 1e-12 needs (2 terms below 2e-6, up
  to 9 near 0.2), the closed forms up to 40 and their expansion in 1/dt
  above (finite for infinite dt): within 1.1e-12 of mpmath for dt from
  1e-8 to 1e3. Bezier3 is the default solver in all geometries and for
  Stokes, so this changes those results by up to about 7e-7. Cost: in
  the 2-D memory case half of the 7.5e8 calls have dt < 1e-3; with the
  terms matched to dt, run time is +3% (2-D) and +1% (3-D) at 1 thread;
  a fixed 18-term series cost +15%.
- N-18: `w2` and `w3` sum the incomplete gamma series (all terms positive)
  below |dtau| = 1 instead of a two-term series below 5e-4 and the closed
  forms above: within 3.5e-16 of mpmath (was up to 2e-6).
