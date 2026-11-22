ids: N-53
date: 2026-11-21
title: One warning per atom and table type when collision tables are used outside their temperature range

The last part of N-53 (d). The collision tables of `CollisionRate` (`OMEGA`, `CE`, `CI`, `CP`, `CH0`, `CH+`, `CH`) are PCHIP in log T and log of the value (N-32). Outside their `TEMP` grid they hold the end value, and nothing said so: `rateCheckOutside` reports only the analytic formulas.

**The change.**
- For each table type, `CollisionRate` counts the records of the atom whose grid does not cover the atmosphere.
- It reports one WARNING per atom and type, once per run. The report is kept with the atom's other once-per-run reports, so later calls (a new n_e) do not repeat it.
- The warning gives:
  - how many tables of that type are outside their range, and how many below and above;
  - the table-locations (record × depth point) affected;
  - the span of the grids, and the temperatures of the atmosphere outside them;
  - what is held. For `CE` that is C, so the effective collision strength grows as T and the rate as √T.

For example, with FAL-C and the shipped Ca II:

  `Atom CA: 10 of 10 OMEGA tables are used outside their temperature range (0 below, 10 above; 120 table-locations of 820; table grids span T = 3000 to 38000 K, atmosphere outside them at T = 41180 to 100000 K); the collision strength is held`

No rate changes. The rule beyond a table stays as it is. The shipped H tables reach 10⁶ K with #544, so with them only coronal temperatures leave the H tables.

**Results.** Unchanged: `run_tests.py --exact` against `modernize` e5a8eaa (60 of 60). Its runs show the new warnings for Ca II (`OMEGA`, to 38 000 K), and for H (`CE`, `CI`, to 30 000 K) until #544.

Tests: new unit test `test_collision_tables`, which runs `CollisionRate` on synthetic atom files.
- Inside its range, a `CE` table of a power-law Ω gives the `OMEGA` table's rate.
- Outside it, every type holds its end value (`CE` with the rate going as √T, `OMEGA` as T^-1/2).
- There is one message per type, with the counts and ranges, none when the tables cover the atmosphere, no repeat on the next call, and the message again after `freeRateReports`.

With the report switched off, 12 of its 220 checks fail.
