ids: N-53
date: 2026-11-21
title: H collision tables (Johnson 1972) extended from 30 000 K to 10⁶ K

Part (d) of N-53, the data. The electron collision tables (`CE`, `CI`) of every shipped H atom end at 30 000 K. Above that, RH holds the table's end value. For `CE` that means the effective collision strength Ω grows in proportion to T, about 3.3× by 10⁵ K, where it should grow roughly logarithmically. FAL-type atmospheres reach 10⁵ K, and in their transition region the hydrogen rates came from that extrapolation.

**Where the tables come from.** All the shipped H tables are Johnson (1972, ApJ 174, 227), evaluated by `tools/make_h.c` at 3000–30 000 K. Johnson's formulas hold at any temperature and tend to the Bethe limit at high temperature: Lyα's Ω at 10⁸ K is 29.1 against the Bethe 29.9 (`claude/RH/evidence/n53d-johnson`). Extending them needs no new source. Better data near threshold (R-matrix, n ≤ 5) is F-12.

**The change.**
- `tools/make_h.c` has 15 temperatures, 3000 K to 10⁶ K, keeping the old six. Above 30 000 K: 5×10⁴, 7×10⁴, 10⁵, 1.5×10⁵, 2×10⁵, 3×10⁵, 5×10⁵, 7×10⁵, 10⁶ K. With these nodes, RH's PCHIP interpolation in log T and log C is within 0.64 % of Johnson above 30 000 K (median 0.09 %, all 190 records of 20-level H). An 11-node grid was within 1.2 %.
- `make_h` also goes up to 20 levels (for H_20), and labels the continuum H II for any number of levels (it took a bound level's label unless there were 11).
- Every shipped H atom gets its `TEMP` line and all its `CE` and `CI` records from `make_h 20`: H, H_3, H_4, H_6, H_6_CRD, H_6_noLybeta, H_9, H_10, H_20, H_CRD. On the old grid the new values are the old ones, except 12 records of H, H_10 and H_CRD that differ in the last printed digit (at most 9×10⁻⁴). Those files came from an older version of the tool. The script is `claude/RH/evidence/n53d-johnson/regen_htables.py`.
- No code change: this pull request changes data and the tool that makes it. Above 10⁶ K the end value is still held. The other parts of (d), a report when an atmosphere leaves a table's range and the rule for `CE` beyond it, are not done here.

**Measured** (FAL-type models to 10⁵ K; old against new data with the same build; `claude/RH/evidence/n53d-johnson/measure`):

| case | H ground level, where it holds ≥ 10⁻³ of H (≥ 10⁻⁶) | Ca II / Mg II | Lyα core, limb (µ = 0.05) / µ = 0.95 | Hα |
|---|---|---|---|---|
| FAL-C, H_6 + Ca II + Mg II | 5.7×10⁻² (0.25) | 3.9×10⁻⁴ / 7.0×10⁻⁴ | 0.16 / 0.033 | ≤ 1×10⁻⁴ |
| FAL-C, H_20 | 4.7×10⁻² (0.22) | — | 0.15 / 0.031 | ≤ 8×10⁻⁵ |
| FAL-A, H + Ca II | 5.6×10⁻² (0.24) | 2.4×10⁻⁴ | 0.14 / 0.023 | ≤ 3×10⁻⁵ |
| FAL-XCO, H_6_CRD + Ca II | 5.7×10⁻² (0.25) | 9×10⁻⁵ | 0.10 / 0.015 | ≤ 2×10⁻⁵ |
| FAL-C shock, H_6_CRD, `SOLVE_NE = ITERATION` | 5.3×10⁻² (0.24) | — | 0.074 / 0.015 | ≤ 1.2×10⁻⁴ |
| T5750 (to 10⁴ K), H_6 | 0 | — | 0 | 0 |

- All runs converge in the same number of iterations with both data sets.
- The changes are where T > 30 000 K: the H ground level at 36 600 K, where it holds 10⁻³ of H, and at 10⁵ K, where it holds 2×10⁻⁵. The Lyα and Lyβ core intensities follow.
- The T5750 model never exceeds 10⁴ K and is bit-identical, as it should be with PCHIP: the new nodes do not reach below 30 000 K.

**With upstream RH** (c95bdf6, its own H atoms regenerated the same way, rhf1d run by hand): the same changes within a few percent, e.g. FAL-C H_6 + Ca II + Mg II 5.7×10⁻² in the H ground level and 0.16 in the Lyα core. Two differences:
- Upstream interpolates the tables with a global cubic spline, so the new nodes also move its rates below 30 000 K, by about 10⁻⁸ in T5750.
- Upstream stops with a segmentation fault in the shock case with `SOLVE_NE = ITERATION`, with either data set (probably B-1, fixed here; not checked).

**Results.** Change for every input with an active H atom and temperatures above 30 000 K. `run_tests.py` takes the atomic data from the candidate tree for both builds, so it cannot see a data change: `--exact` against `modernize` a34785a passes. Its H_6 cases (`1d/twoactive`, `sph/small`) change with the data, as above.
