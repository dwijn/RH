ids: N-53
date: 2026-11-20
title: H2- and H2+ free-free: no more zero opacity beyond the last tabulated wavelength

`H2minus_ff` (Bell 1980 table, 19 wavelengths 0 and 350.5–15188.3 nm × 8 values of θ = 5040/T, θ clamped to [0.5, 3.6]) and `H2plus_ff` (Bates 1952 table, 15 wavelengths 0 and 384.6–10000 nm × 10 temperatures, T clamped to 2500–12000 K) returned FALSE for λ ≥ the last tabulated wavelength, that is, no opacity at all, although both processes are nonzero there. This is part (c) of N-53; the agreed policy is that a table outside its range must not return 0 for a physically nonzero quantity, that it holds the end value or continues a power law where the physics gives one, and that a message is written once per table per run.

**Check of the power law.** The exponent of κ(λ) fitted between the last wavelengths of the H2- table (9113.0, 11391.3, 15188.3 nm) is, per θ = 0.5, 0.8, 1.0, 1.2, 1.6, 2.0, 2.8, 3.6, between the last two columns: 1.998, 1.995, 1.989, 1.997, 1.994, 2.004, 2.001, 1.987 (the interval before: 1.999, 1.997, 1.998, 1.990, 1.985, 1.991, 1.988, 1.985). The table follows λ² to 1.5 % over the whole θ range, as free-free absorption of a slowly varying cross-section does at long wavelength. For H2+ the table is flat or falling at the long end (κ at 10000 nm over κ at 5000 nm is 0.86 at 2500 K and 0.79 at 12000 K), so no power law is available.

**Fix.**
- `H2minus_ff`: for λ > 15188.3 nm, κ(λ, θ) = κ(15188.3 nm, θ)·(λ/15188.3 nm)². The same factor holds for every θ, so it is applied to the last-row interpolation in θ.
- `H2plus_ff`: for λ > 10000 nm, the last column is held (the same values as at 10000 nm).
- Both routines write one WARNING per run (an `atomic_flag`, reset by `hydrogenPrepare`, so the routines stay callable from several threads) the first time a wavelength beyond the table is used. Runs that stay inside the tables write nothing.
- The λ = 0 zero row at the short end is left as it is: below 350.5 nm (H2-) and 384.6 nm (H2+), κ is interpolated linearly down to 0 at 0 nm. The zero row only closes the table at the origin; it is not a fit. For H2- the first two columns give an exponent of 2.0 (κ(414.2)/κ(350.5) at θ = 0.5 is 1.40 for a wavelength ratio of 1.18), so at 300 nm linear interpolation gives 86 % of the 350.5 nm value where λ² would give 73 %. It is not changed here.

**Measured.** `tests/unit/test_h2ff.c` (new, no numpy): the table values at the last node, κ ∝ λ² of H2- at 20, 50, 100 µm, the held H2+ value, continuity at the last node (1e-9 relative below it), and the message counts (one per routine). Opacity χ (m⁻¹) before (always 0) and after, at cool (T = 4000 K, n_e = 10¹⁷, n_H = 10²³, n_H2 = 10²¹, n_p = 10¹² m⁻³) and warm (T = 6500 K, n_e = n_p = 10¹⁹, n_H = 10²³, n_H2 = 10¹⁸ m⁻³) points, with H- free-free for scale:

| λ (µm) | depth | H2- ff | H2+ ff | H- ff |
|---|---|---|---|---|
| 20 | cool | 1.05e-8 | 9.0e-15 | 1.68e-6 |
| 20 | warm | 1.41e-9 | 5.6e-8 | 1.79e-4 |
| 50 | cool | 6.55e-8 | 9.0e-15 | 1.05e-5 |
| 50 | warm | 8.81e-9 | 5.6e-8 | 1.12e-3 |
| 100 | cool | 2.62e-7 | 9.0e-15 | 4.21e-5 |
| 100 | warm | 3.53e-8 | 5.6e-8 | 4.47e-3 |

H2- ff adds 0.6 % to the H- free-free opacity at the cool point and 0.8 % at the warm point at these wavelengths (it is a trace molecule there); H2+ ff is 10⁻³ of H- at the warm point and nothing at the cool point.

**Results that change.** Only runs with wavelengths beyond 10 µm (H2+) or 15.2 µm (H2-). No case of `tests/run_tests.py` reaches them (no message in the 53 run logs); its 14 candidate-versus-baseline comparisons (1-D, 2-D, 3-D, spherical, Stokes) are identical to 0.0e+00, and the multi-thread comparisons are at 10⁻¹² to 10⁻¹³ as for unmodified code. Results at shorter wavelengths are bit-identical.
