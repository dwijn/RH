ids: Q-14, B-27, B-28
date: 2026-10-19
title: Zeeman splitting and polarized profiles tested from first principles; NaN profiles for J = 0 -> 0 lines and fields along a ray

- `test_zeeman` (new) builds fictional lines: the normal triplet, a
  triplet with an explicit effective Lande factor, the anomalous
  patterns of 2S1/2-2P1/2 and 2S1/2-2P3/2 (as Na D1, D2), 5P2-5D2
  (dJ = 0), 5P1-5D0 (upper J smaller), the null line 5D0-5F1
  (g(5F1) = 0) and 1S0-3P0 (J = 0 -> 0). It checks the components,
  shifts and strengths (3j symbols), g from LS coupling and g_eff (its
  formula and the barycentre of the sigma components), and the seven
  profiles from `Profile()` for seven fields (0 to 0.5 T, toward and
  away from the observer, transverse, along an oblique ray), two rays
  and both directions, against the propagation matrix derived from the
  circular and linear normal modes (docs/unit-tests, section 8).
  All agree to 1e-16 once Stokes V changes sign: RH is internally
  consistent, with positive V meaning negative helicity (opposite to
  IAU 1973; positive blue lobe for a field toward the observer).
- Milne-Eddington: RH's profiles of one point turned into a slab with
  constant field and line-to-continuum ratio and S linear in tau_c,
  solved with the polarized solver, agree with the Unno-Rachkovsky
  solution to 1e-13 for every Stokes parameter across the line; in a
  weak field V = -g_eff vB cos(theta) dI/dv holds to 4-6e-4 of max |V|
  (the O(vB^2) terms).
- B-27: a J = 0 -> 0 line was marked polarizable and `Zeeman()`
  normalized its only component (strength 0) by zero: NaN profiles.
  The decision moved to `zeemanPolarizable()` (zeeman.c), which rejects
  J = 0 -> 0; `Zeeman()` no longer divides by an empty sum.
- B-28: `Bproject()` (1-D, 2-D, 3-D) divided by 1 - cos^2(gamma) for the
  field azimuth, 0/0 for a field exactly along an oblique ray: NaN in
  the Q and U profiles of every polarizable line there. Now chi = 0
  there (its terms carry sin^2(gamma) = 0). Results elsewhere are
  bit-identical.
- The unit-test report gains section 8: the Zeeman pattern, the
  propagation matrix from the normal modes, RH's sign of V,
  Milne-Eddington and the weak-field limit.
