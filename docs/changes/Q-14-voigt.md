ids: Q-14, B-20
date: 2026-10-17
title: Unit tests, starting with the Voigt functions; VoigtHui finite at a = v = 0

- New test tier: unit tests in `tests/unit/`, small programs linked against
  `librh` that check single routines on problems with known answers. They
  run with `ctest` (label `unit`), in CI on every build and with
  AddressSanitizer and UBSan in the sanitizer job. The code-format check
  now covers `tests/` as well (except generated tables).
- `test_voigt` checks the Voigt function H and the Faraday-Voigt function
  F of every algorithm against the Faddeeva function w(v + ia) computed
  with mpmath at 40 digits (`tests/unit/make_voigt_reference.py`, 18 damping
  parameters from 0 to 100 times 39 frequencies from 0 to 10⁴), plus
  symmetry in v and finite values. Measured accuracy (gcc, x86-64), as
  relative error where H ≥ 10⁻¹² H(a, 0), and as error relative to line
  centre:

  | Algorithm | Used for | rel. error | rel. to peak | F, rel. to max F | ns/call |
  |---|---|---|---|---|---|
  | ARMSTRONG | all line profiles without Stokes | 1.2e-4 (v ≈ 5, small a) | 1.3e-6 | — | 97 |
  | RYBICKI | `giigen.c`, `solve1d`, `solve2d` | 1.8e-8 | 7.6e-11 | — | 192 |
  | HUI_ETAL | angle-dependent PRD (`rii.c`, `riiplane.c`) | no Gaussian wings for small a | 3.9e-6 | 6.4e-6 | 66 |
  | HUMLICEK | profiles with Stokes (H and F) | 8.4e-5 | 2.9e-5 | 4.9e-5 | 46 |

  The test's bounds are about twice these, so that a real loss of accuracy
  fails. These numbers are the input for Q-11 (which algorithm to keep).
- B-20: `VoigtHui()` returned NaN at a = v = 0 (0/0: numerator and
  denominator were each multiplied by z once too often). The polynomials
  are now evaluated by plain Horner; other values change only by rounding.
