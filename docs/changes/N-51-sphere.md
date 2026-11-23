ids: N-51
date: 2026-11-23
title: Background scattering iterated at each wavelength in spherical symmetry (rhsphere)

The background scattering iteration of N-51 (with the stopping test of N-66 and the default of N-65) now also runs in `rhsphere`. Until now `rhsphere` kept the plain lambda iteration, one step per formal solution, which needs of order 1/ε sweeps where a scattering line is optically thick. This closes N-51; 1-D polarized runs (N-67) and 2-D and 3-D (N-68) are on the wishlist.

**The change.**
- **The sweep.** `rhsphere` supplies the formal solution of the iteration: Feautrier along every ray, J from P with the ray weights of each shell, and the tridiagonal of dJ/du for the preconditioned step (`BACKGROUND_OPERATOR = TRIDIAGONAL`, the default).
- **The operator.** `Feautrier` (`rhsphere/feautrier.c`) returns, on request, the tridiagonal of dP/dS besides its diagonal Psi. Column j of the inverse of the Feautrier matrix falls off from its diagonal as P_k = P_k+1 / (1 + F_k) above j and P_k = P_k−1 / (1 + G_k) below it, with F and G of the forward and backward elimination, so the off-diagonals cost no more than Psi.
- **Positivity.** The Feautrier matrix (`STANDARD`) is an M-matrix with row sums ≥ 1, so dP/dS ≥ 0 with row sums ≤ 1, and the ray weights of each shell add up to 1 (to rounding, in the test below). ‖Λσ‖∞ is then at most the largest albedo, which is what the error bound of the stopping test (N-66) needs, and no wavelength needs the linear fallback of N-52.
- **The fallback hook.** Whether iterated wavelengths switch to the linear solver was decided in `bgscatter.c` from `S_INTERPOLATION`, which `rhsphere` reads but does not use. Each geometry now says it (`bgScatterLinearFallback`): `rhf1d` with `S_BEZIER3` only.
- `rhsphere/run/keyword.input` documents the keywords.

**Results.** FAL-C in spherical symmetry (`FALC_82_s`, 88 rays), the example's atoms (H₆ and Ca II active, PRD), the Kurucz lines of 256–260 nm with `RLK_SCATTER` (7292 lines) and the wavelength table 256.5–259.5 nm of the N-65 reproducer (6000 points), `ITER_LIMIT` 10⁻³, 8 threads. Errors are of the emergent intensity of the central ray, against the new code at `ITER_LIMIT` 10⁻⁶ and `BACKGROUND_ITER_LIMIT` 10⁻⁹.

| | iterations | background solves | wall | largest error | median | wavelengths off by > 10 % |
|---|---|---|---|---|---|---|
| before | 26 | — | 0:48 | 0.90 (I 10 times too small) | 1.7×10⁻³ | 24 % |
| now | 15 | 476k | 0:49 | 3.5×10⁻⁵ | 2×10⁻⁶ | 0 |
| before, `ITER_LIMIT` 10⁻⁶ | 225 (reported converged) | — | 5:34 | 0.64 | 3×10⁻⁷ | 4 % |

**Check against the lambda iteration.** Part of the window (257.70–257.80 nm, 200 points), all atoms passive (LTE), so that only the background scattering is iterated:
- The iteration (`BACKGROUND_ITER_LIMIT` 10⁻¹⁰) against the plain lambda iteration after N sweeps (`BACKGROUND_N_MAX_ITER` 0): J/J_iterated − 1 is −0.87 at worst after 1000 sweeps, −0.56 after 4000 and −0.084 after 16000, and never above 7×10⁻¹² (the lambda iteration of a positive operator rises monotonically to its fixed point).
- The lambda iteration of the code before this change gives the same J as the new code with `BACKGROUND_N_MAX_ITER` 0, bit for bit.

**`run_tests.py`** against B-73 (the base): 37 of 38 cases identical; `sph/small` (the shipped spherical example, stopped after 4 iterations, unconverged) changes by up to 2.1×10⁻² (I at 95.05 nm), because every wavelength is now iterated. Run to convergence (`ITER_LIMIT` 10⁻⁶, 28 iterations before and after), the example agrees with a run at `BACKGROUND_ITER_LIMIT` 10⁻¹⁰ to 6.8×10⁻⁸ before and 2.6×10⁻⁹ now: it has no strongly scattering wavelengths. The 1-D, 2-D and 3-D cases are byte-identical.

Tests:
- `test_formal_sphere`, `check_operator`: along every ray, Psi and the tridiagonal from `Feautrier` against the columns of dP/dS from perturbations of S, with the opacity varying over three decades: equal to 3×10⁻¹⁵ (the test allows 10⁻¹²); every entry of dP/dS ≥ 0, row sums ≤ 1, ray weights of each shell adding up to 1.

Evidence: `/hao/solar1/dwijn/claude/RH/evidence/n51-sphere` (README, run setup, comparison scripts).
