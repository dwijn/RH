ids: N-53
date: 2026-11-20
title: Molecules no longer vanish outside the temperature range of their fits (partition function held, ln K extrapolated in 5040/T)

`partfunction` and `equilconstant` (`chemequil.c`) returned 0 for T outside [Tmin, Tmax], the range of the fit read from the molecule file. All 19 shipped molecules have Tmin = 1000 K (H2 and most others Tmax = 9000 K), so below 1000 K and above Tmax every molecule had K = 0 and was silently absent. `Hydrostatic` (`rhf1d/hydrostat.c`) also dropped H2 from its Newton system there, and `readMolecule` stored pf = 0, which `opacity.c` and `initial_xdr.c` divide by for an active molecule.

**The fix** (the agreed N-53 policy: tables and fits outside their range must not return 0 for a physically nonzero quantity).
- **Partition function:** evaluated at T clamped to [Tmin, Tmax], so it is held at its edge value.
- **Equilibrium constant:** ln K is extrapolated linearly in θ = 5040/T from the nearer edge, with the value and the θ-derivative of ln K at that edge (van 't Hoff: ln K is close to linear in 1/T, the dissociation energy is the slope). It is done for all fit types (`KURUCZ_70`, `KURUCZ_85`, `SAUVAL_TATUM_84`/`IRWIN_81`, `TSUJI_73`) through a log form of each fit, `lnEquilConstant`. The slope is a centred difference over ±0.1 % of θ_edge: the fits are smooth (polynomials in T, ln T or log θ), so the truncation error is negligible and no per-fit analytic derivative has to be kept in step with the fit. ln K is capped at 700 to stay below the double range (reached below about 170 K for CO). Inside the range the old code runs unchanged.
- **Message:** `readMolecule` prints one WARNING per molecule when any depth point is outside [Tmin, Tmax] (number of points, atmosphere T range), not per depth point.
- **`hydrostat.c`:** H2 is always in the system (5 equations). The lower cut is gone, and so is the upper one: with the extrapolated K above 9000 K (FAL-C, 33 points above Tmax up to 10⁵ K, nH2/nH down to 3×10⁻⁸) the Newton iteration converges at every point, no "not converged" warning; run with active H, `SOLVE_NE = ITERATION`.
- **Users of pf = 0 / K = 0 as an "absent" flag:** none found. `Phi[i]` (`chemequil.c`) is only used as K; `molecule->pf` is read in `opacity.c` (guarded by n > 0, now pf > 0 always), `initial_xdr.c` and `ltepops.c`; `Tmin`/`Tmax` are only written to `molecules.out` and read by the Python reader.

**Accuracy of the extrapolation** (`test_molecule_range --table`: Tmin raised by dT, K at the original Tmin extrapolated from there against the fit, ratio):

| | Tmin+200 K | Tmin+500 K | Tmax-200 K | Tmax-500 K |
|---|---|---|---|---|
| H2 | 1.0015 | 1.0044 | 1.0007 | 1.0024 |
| CO | 0.9972 | 0.9776 | 0.9928 | 0.9612 |
| OH | 1.0019 | 1.0036 | 0.9998 | 0.9986 |
| CH | 0.9972 | 0.9761 | 0.9996 | 0.9977 |
| H2O | 1.0015 | 0.9957 | 1.0034 | 1.0141 |
| H2+ | 0.9676 | 0.8360 | 1.0342 | 1.2061 |

Within 0.5 % over 200 K for all but H2+ (KURUCZ_70, strongly curved, 3 %). In the full code (FAL-C, `HYDROSTATIC = TRUE`, H2 range cut to [5200, 9000] K so that the extrapolation covers 4500 to 5200 K): n(H2) is within 0.9 % of the run with the real fit (max over 4500 to 9000 K, at the point farthest below the edge), CO, OH, CH within 3×10⁻⁶, ne and nHtot within 1.3×10⁻⁶ and 1.8×10⁻⁶; with the old cut, H2 is 0 below 5200 K and ne and nHtot are off by 4.5×10⁻⁴ and 6.9×10⁻⁴, CO by 1.1×10⁻³.

**Results.** FAL-C and the shipped regression cases (T from 4170 K up to 10⁵ K) change only at the two deepest FAL-C points (9140 and 9400 K, above Tmax of H2 and most others), where H2 now exists (nH2/nH 2×10⁻⁵): the hydrogen atoms there change by 4×10⁻⁵ and the emergent intensity by at most 7×10⁻⁹ (at 90 nm), `run_tests.py` (tolerance mode against a build of `modernize`, 4b665de): 20 of 26 comparisons pass; the six that differ are the cases that contain FAL-C with its deepest points, the five 1-D cases (`nostokes`, `fieldfree`, `fullstokes`, `crd`: 6.8×10⁻⁹ at 90 nm, `twoactive`: 7.0×10⁻⁹ at 177.6 nm) and `sph/small` (9.8×10⁻⁹ at 177.6 nm), all above the 10⁻⁹ tolerance; Q, U, V differ by at most 2.5×10⁻¹²; the 2-D and 3-D cases differ too but within their tolerances (at most 3×10⁻¹⁰ in I, 1.7×10⁻¹⁰ in V; tolerances 10⁻⁷ and 10⁻⁹), and the memory case is unchanged (103 MB). The change is therefore a change of results: label `changes-results`. Atmospheres with T < 1000 K (or between Tmax and the hot end) now have molecules, as the equilibrium says.

**Cool pocket test** (`FALC_cool800.atmos`, `FALC_cool1500.atmos`: FAL-C with T lowered by a Gaussian in log column mass around the temperature minimum to Tmin = 800 and 1500 K; made by `mk_cool.py`; molecule list of the shipped `molecules.input`). At 800 K the old code gave n(H2) = n(CO) = n(OH) = 0 (all hydrogen atomic at 800 K); now 2 n(H2)/nHtot = 1.00, n(CO) = 5.0×10¹⁷ m⁻³ (all carbon), n(H2O) = 4.3×10¹⁷ m⁻³, n(OH) = 4×10⁴ m⁻³.

**Known limit this exposes (not changed here).** `ChemicalEquilibrium` is called with `N_MAX_CHEM_ITER` = 10 Newton iterations from the "everything dissociated" start. Below about 2000 K this is not enough, even inside the fit range and with the old code: at 1500 K the old code ends with "Iteration not converged" and negative number densities (CH −8×10¹⁵ m⁻³, CN, H2O, C2 negative), and the same now happens for T < 1000 K. With `N_MAX_CHEM_ITER` = 100 the cool-pocket atmospheres (800 K and 1500 K) converge at every point (tried 30: not enough, 100 and 300: all converge). This needs its own fix before cool atmospheres give usable molecular densities.

Tests: new unit test `test_molecule_range` (all fit types, continuity of K and of dlnK/dθ at both edges, straight line in θ outside, pf held, K positive and finite).
