# RH modernization log

A working record of the effort to modernize RH: issues found, baseline
measurements, and the changes made, each with how it was verified against
upstream. Newest changes go at the top of the changelog.

- **Goals:** modern code quality and one consistent style; speed and
  efficiency; accuracy and robust convergence in strongly NLTE problems; a
  smaller memory footprint. Speed and memory matter most for 2-D/3-D runs with
  several many-level active atoms; convergence matters most where NLTE effects
  are strong.
- **Not a goal:** keeping the diff against upstream small.
- **Reference:** upstream is `../RH-upstream` (commit `c95bdf6`). It is built
  with the Fortran Voigt routines (`HAVE_F90`); a build without them is used
  when Stokes results must match bit for bit (see C-1 in the changelog).
- **Companion document:** the *RH code map* artifact (how the code is put
  together, which parts are threaded) from an earlier session:
  <https://claude.ai/artifact/UXA2hsvfBEMVFQEzBCLPFW>. Changes below note which
  of its sections they make out of date.

## At a glance

*Status on 2026-11-10.* The work is on the fork
[dwijn/RH](https://github.com/dwijn/RH), integration branch `modernize`.
Every change is a pull request, verified against upstream RH (`c95bdf6`)
or against its base commit. 151 issues are logged; 75 are fixed and 2
closed without a fix.
Plan steps 1–3 (crash fixes, build and hygiene, tests) are done, and
step 4 (the solver) is under way. Steps 5–7 follow: speed, memory and
structure.

- **Wrong answers that looked right, also in upstream RH.**
  - 52 of the bugs found also exist upstream (§3.9), and 42 of them are
    fixed here.
  - **Formal solvers:**
    - 2-D and 3-D long characteristics stopped one section short (B-24):
      7–11 % off on a slab with an exact solution, and up to 25 % in the
      shipped examples.
    - The polarized DELO-Bezier solvers had two wrong terms (B-25):
      Stokes V was 38 % off for Ca II on FAL-C.
  - **Statistical equilibrium:** converged populations did not satisfy
    statistical equilibrium where transitions overlap (B-6).
  - **Chemical equilibrium:** the stored H⁻ density was too large by the
    inverse neutral fraction of hydrogen (B-40), up to 10¹¹ at the top of
    FAL-C. With the fix, continua below 200 nm change by up to 9–15 %
    on inclined rays, and the limb rays of the spherical example by up
    to 17 %.
  - **Atomic data:**
    - Two model atoms had wrong collision strengths (B-42 for C II,
      B-43 for Ca II).
    - Tabulated rates are now interpolated monotonically in log-log, so
      they no longer go negative (N-32).
  - The unit tests check against problems with known answers (exact
    slabs, the analytic polarized solution, first-principles Zeeman
    patterns), which is how most of these were found.
- **Convergence: a new solver, `ANDERSON`, the default in 1-D since
  S-10.** RH's Ng is replaced by stabilized Anderson acceleration of an
  explicit, deterministic MALI map G (S-2, S-3; design and reasoning in
  [`SOLVER.md`](https://github.com/dwijn/RH/blob/modernize/docs/SOLVER.md)).
  - On the 1-D test set (§2.5) it converges every case. That includes
    the six where Ng fails, three of them with a crash.
  - Applied to the same map, it needs 26–77 % less work than Ng on the
    easy cases and 3× less on the slow one.
  - On the shipped examples it never needs more iterations than Ng. At
    `ITER_LIMIT` 10⁻⁵, `rhsphere` needs half as many.
  - Rate equations are solved by GTH elimination (N-31). He_large
    stalled at 2×10⁻³ because of rounding in the LU solve. Its plain
    iteration now converges, and so does `ANDERSON`; with Ng it still
    depends on the extrapolation (N-28).
  - Runs report whether they converged, also in the exit status (N-1).
  - Every keyword that is set is checked against the solver, and a
    conflict stops the run with a clear message.
  - A Newton–Krylov solver (S-4) is also built: G is bit-reproducible at
    a fixed thread count, so its derivatives can be taken by finite
    differences.
- **Speed:**
  - PRD scattering integrals are 4× faster, and whole PRD runs 1.3–1.4×
    (P-8).
  - The vectorized per-transition loops and `-O3` give 1.15–1.23× (P-7).
  - Threading gave 1.7–1.8× on 2-D and 3-D cases at 16 threads in the
    perf set of 2026-10-16 (§2.3). It costs 8–15 % more memory for
    per-thread copies (M-3).
- **Memory:**
  - Line profiles and the background are stored in single precision:
    the 2-D perf case peaks at 10.7 GB instead of 17.5 GB (M-1, M-7).
  - Realistic 3-D runs are still out of reach, at tens of TB. The
    causes are line profiles stored per ray, dense rate and collision
    matrices, per-thread copies and PRD weights (§3.3, plan step 6).
  - Anderson's history can cover only the populations and the PRD ρ
    (S-11), which makes it the default in 2-D and 3-D at 1.65× Ng's
    memory on the 2-D perf case instead of 3.8×.
- **Code and testing:**
  - CMake builds. The Fortran is gone, and one C Voigt function replaces
    five.
  - Zero compiler warnings, enforced as errors in CI, and one code style
    (clang-format).
  - CI on Linux and macOS with gcc and clang, also under sanitizers. 27
    tests run in ctest: 17 unit, 6 end-to-end and 4 that start each
    program.
  - Every pull request is compared with its base commit, and drift from
    upstream is reported.
  - A stress suite covers convergence and performance, including cuts of
    a Bifrost snapshot.
- **For codes that call RH as a library (STIC):**
  - 1-D columns get the robust solver by default.
  - Each run says whether it converged.
  - The map G is deterministic, which is what derivatives and
    Newton–Krylov need.
  - The next structural step is that errors return to the caller
    instead of calling `exit()` (Q-7, high severity, first in plan
    step 7).
- **Next:**
  - Write up the S-11 measurements; settle the PRD sub-iteration
    damping (N-38).
  - The rest of step 4, then speed (step 5) and memory (step 6).
  - Report the upstream bugs of §3.9 to upstream RH.

Issue IDs: **B** bug or crash, **N** numerics/convergence, **P** speed,
**M** memory, **T** threading, **Q** code quality / build, **F** feature,
**S** solver redesign (design and work items in [`SOLVER.md`](https://github.com/dwijn/RH/blob/modernize/docs/SOLVER.md)).
Severity: high, medium, low. Status: open, fixed (with changelog entry),
wontfix.

---

## 1. Test infrastructure

| Tool | Purpose |
|---|---|
| `tests/run_tests.sh` | Regression: candidate vs baseline at 4 iterations on 1-D (3 Stokes modes, CRD, two active atoms), 2-D, 3-D, spherical; thread-count consistency (2, 3, 8 threads vs 1); optional ASan build; peak memory on a large 2-D grid; Kurucz blank-line crash. `--exact` for byte-identical changes. |
| `tests/run_stress.sh` | Stress/convergence/performance: runs NLTE problems **to convergence** on one or more builds and reports iterations, final change, trend (monotone / oscillating / stalled), wall time, peak memory, and the population and spectrum differences between builds. `--set conv` (24 cases) or `--set perf`. Run with `--work /hao/solar1/dwijn/rh-stress` (large outputs). |
| `tests/analyze_run.py` | Parses a run's log and `time.txt` (iterations, convergence, trend, status) and compares `pops.*.out` between runs. |
| `tests/make_perturbed_atmos.py` | Builds inhomogeneous 2-D/3-D atmospheres from FAL-C: transition region shifted by up to ±250 km, chromospheric T modulated by ±35 %, velocities up to ±8 km/s (vertical) and ±4 km/s (horizontal). Deterministic and periodic. |
| `tests/atmos/FALC_v10.atmos`, `FALC_shock.atmos` | 1-D FAL-C with a ±10 km/s velocity field, and with a 2.2× temperature spike plus a velocity jump. |

Gaps in the test suite before this effort: no case ran to convergence (all
stopped at 4 iterations), no strongly NLTE or many-level atom, no
inhomogeneous multi-D atmosphere (the shipped 2-D/3-D models are FAL-C
repeated in every column), no timing case. No unit tests, no CI.

Notes:
- 1-D stress cases run on 1 thread, because thread count changes summation
  order and some problems are sensitive enough that this changes the
  converged answer (N-3).
- The pre-existing `tests/make_tiled_2d_atmos.py` reads the 2-D file as if x
  were the slow index; RH stores height-major with x fastest (`k*Nx + l`).
  Its tiled atmosphere is scrambled. Harmless for the memory check it serves,
  but should be fixed (Q-12).

---

## 2. Baseline results

### 2.1 Regression suite (`run_tests.sh`, WIP tree vs upstream)

32 of 34 checks pass. The two failures are `1d/fieldfree` and `1d/fullstokes`
(Q/U differ by 1.5e-6, I by 4e-8). **Cause verified:** upstream uses the
Fortran Humlicek/Hui Voigt routines, the WIP tree the C ones. Upstream built
without `HAVE_F90` is bit-identical to the WIP tree in both cases. All other
cases are bit-identical to upstream at 1 thread; thread-count noise is 1e-12
to 5e-9; peak memory on the tiled 2-D grid is 208 MB vs 206 MB.

### 2.2 Convergence (`run_stress.sh --set conv`)

1-D cases at 1 thread, `ITER_LIMIT = 1e-3`, `N_MAX_ITER = 200` unless noted.
**The WIP tree is bit-identical to upstream in every case that produces
output** (populations and spectrum differences exactly 0), and takes the same
wall time at 1 thread: its speedups come from threading only.

| Case | Iterations | Final change | Outcome (upstream = WIP) |
|---|---|---|---|
| FALC H_6 + CaII + MgII | 22 | 7.5e-4 | converged |
| FALA H (9 lev) + CaII | 26 | 8.7e-4 | converged |
| FALC H_20 | 16 | 7.5e-4 | converged |
| FALC MgI_66 | 16 | 4.6e-4 | converged (38 s) |
| FALC Fe33 | 20 | 8.9e-4 | converged |
| FALXCO H_6_CRD + CaII | 19 | 3.0e-4 | converged |
| FALC_v10 CaII, Bezier | 14 | 7.4e-4 | converged |
| FALC H_6 + CaII + MgII, no Ng | 50 | 9.8e-4 | converged, stalling (Ng saves 56 %) |
| FALC He_large | 108 | 6.9e-4 | converged slowly; 48 at 4 threads with a 1.4 % different answer (N-3) |
| FALXCO H_6 + CaII (PRD) | 200 | 0.27 | **not converged** (N-4: same atoms in CRD converge in 19) |
| T5750 H_6 | 200 | 4.2e-3 | **not converged**, oscillating |
| FALC H_6, no Ng, limit 1e-8 | 300 | 5.4e-4 | **stalls** (PRD limit cycle, N-4) |
| FALC CI–VI Judge (115 lev) | 150 | 1.0 | **diverges**, SIGFPE after 15 min (N-7) |
| FALC_v10 CaII, parabolic | 3 | 1e13 | **SIGFPE** (N-5) |
| FALC_shock H_6 | 14 | 6e16 | **SIGFPE** with Ng (N-2) |
| FALC CaII_large | 0 | — | **SIGFPE** at setup (B-4) |
| FALC H_6 + CaII, `SOLVE_NE` | 1 | — | **segfault** (B-1) |
| FALC CaII, `PRD_NG_ORDER = 2` | 1 | — | **SIGFPE** (B-2) |
| FALXCO H_6, `COLL_RELAX_*` | 1 | 1.4e-12 | **reports converged with LTE populations** (B-3) |

None of the non-converged runs says so: they exit with status 0 and write
normal output (N-1).

Multi-D, 8 threads:

| Case | Upstream | WIP |
|---|---|---|
| 2-D perturbed FAL-C 32×82, H_6 + CaII + MgII (PRD), Bezier | converged in 76 iterations, 1318 s, 1.20 GB | 76 iterations, **644 s (2.0×)**, 1.28 GB; populations agree to 4e-8 (summation order) |
| 3-D perturbed FAL-C 8×8×82, H_6 + CaII | diverges, SIGFPE at iteration 20 (785 s) | same (578 s) — see N-10 |

First run (12 cases, 1-D at 4 threads), WIP vs upstream:

| Case | Upstream | WIP |
|---|---|---|
| FALC H_6+CaII+MgII, FALA H+CaII, H_20, MgI_66, Fe33 | converge in 16–26 iterations | same iterations; 2.4–3.6× faster wall time; peak memory +5–30 % |
| FALC He_large | 48 iterations | 75 iterations, populations differ by 1.1 % (see N-3: thread-order sensitivity, not a WIP bug; bit-identical at 1 thread) |
| FALXCO H_6+CaII | no convergence in 200 | same |
| CaII_large; SOLVE_NE=ITERATION | crash | crash |
| 3-D perturbed 8×8, H_6+CaII | diverges from iteration 1, SIGFPE at 20 (13 min) | same |

After the fix of B-29 (#187), `run_stress.sh --set conv` on `modernize`
(2026-10-20, with the Bifrost cuts), changes from the tables above:

| Case | Before | After #187 |
|---|---|---|
| FALXCO H_6 + CaII (PRD) | not converged after 200 (0.27) | converged in 109 |
| FALC He_large | 108 | 39 |
| FALC CI–VI Judge | SIGFPE after 15 min | no crash; not converged after 200 (delta 13) |
| FALXCO H_6, `COLL_RELAX_*` | reported converged with LTE populations (B-3) | not converged after 200 (7.8e-2) |
| Bifrost column 260, CaII PRD | SIGFPE at iteration 3 | converged in 30 |
| Bifrost 2-D slice 126×89, CaII angle-averaged PRD | SIGFPE at iteration 2 | limit cycle, 1.3e-2 after 200 (N-4) |
| Bifrost 3-D 24×24×89, CaII hybrid PRD | crash | converged in 43 (1 h 24 min, 35 GB, 16 threads) |
| 3-D perturbed 8×8×82, H_6 + CaII | SIGFPE at iteration 20 | SIGFPE at iteration 5 (N-10) |

FALC_shock H_6 still dies with SIGFPE with Ng (N-2), FALC H_6 at 1e-8
still stalls (N-4), and CaII_large now stops with a clear error about its
collision data (B-4).

### 2.3 Performance

Measured 2-D profile (ns_y27_63x64, Nspace 4032, Nspect 3151, 12 rays; H 9
levels, CaII 6, MgII 11; PRD angle-averaged), per main iteration:

| Part | Share / time |
|---|---|
| Short-characteristics solve (`Quadr` 20 %, `Piecewise_2D` 14 %) | ~36 % |
| Per-transition work (`Opacity` 17 %, `addtoGamma` 8 %, `addtoRates` 8 %, `addtoCoupling` 5.5 %) | ~38 % |
| PRD scattering integral (`Linear`, `Hunt` 2e9 calls, `prdIntegrate`, `GII`) | ~12 % |
| Voigt profiles (one-time) | ~7.5 % |
| Statistical equilibrium | 0.02 s |
| "Spectrum & Operator" pass, -O2 / -O3 -march=native | 41.3 s / 35.6 s |
| Peak memory | 2.3 GB, of which 1.3 GB line profiles stored per ray |

Earlier profile (code map, Bifrost 2-D slice, Ca II PRD, 250 s at 1 thread):
~40 % serial (PRD scatter integral, output writing, profiles), so speedup
capped near 2.5×. With the WIP threading the same run went 124 s → 69 s at 8
threads; writing `opacity.out` is 26 s of the 69 s.

Perf set (`run_stress.sh --set perf`, 2026-10-15/16, idle 48-core node,
16 threads, fixed iteration count). Upstream is c95bdf6 built without
Fortran; modernize is 2868d18. Atoms H_6, CaII, MgII (PRD) and Fe33 active,
on perturbed FAL-C (82 depths):

| Case | Upstream | modernize | Speedup | Peak RSS upstream / modernize |
|---|---|---|---|---|
| 2-D 128×82, 5 iterations | 1776 s | 1005 s | 1.77× | 15.1 / 17.4 GB (+15 %) |
| 3-D 24×24×82, 3 iterations | 13350 s | 7865 s | 1.70× | 127 / 137 GB (+8 %) |
| 3-D 16×16×82 MgI_66, 3 iterations | — | — | — | not run (time cap) |

- Both trees follow identical iteration histories, so the speedup is pure
  threading.
- The extra memory in modernize is the per-thread Γ and rate copies (M-3).
- Output is a large serial share in 3-D. Writing took about 21 minutes,
  mostly `opacity.out` at 109 GB of the 136 GB run directory. That is 9 % of
  the upstream wall time and 16 % of modernize's (P-10, M-7).
- The 3-D case is far from converging, with population changes of 1.6e8,
  5.6e5 and 3.1e4 over its three iterations (N-10). The timing is still
  representative per iteration.

### 2.4 Mixed precision (feasibility study, 2026-11-06)

Can RH run (mostly) in single precision, to save memory and time and to run on
consumer GPUs? Branch `exp/mixed-precision` (local, not for merging) adds a
run-time switch `RH_MIXPREC` that puts selected stages in float: line profiles,
background opacities, populations, active opacities, J, rate and Γ terms and
their accumulators, the SE solve, the 1-D formal solver (real float kernels),
and the PRD scattering integrals. Report: [https://claude.ai/artifact/UxvvK4AzCyLc1tcFxBAuD6](https://claude.ai/artifact/UxvvK4AzCyLc1tcFxBAuD6).

- A naive float port fails on strongly scattering lines. J − S and Ieff = I − Ψη
  cancel, and rounding u gives population errors ~u/ε. Hydrogen stalls at a
  population change of 1e-2 with errors of 1–3 %; a two-level atom has 0.15 %
  at ε = 1e-4 and NaN at ε = 1e-8. Double has the same law (N-30).
- With four changes a float RH converges in the double iteration count on
  every 1-D test case and agrees with double to ≤ 6e-6:
  1. a hybrid formal solver that carries D = I − S where |I − S| ≤ S/2 and I
     elsewhere (D alone loses the EUV, by up to 16 %);
  2. rate terms formed from D and 1 − Ψχ, and every S and η computed once and
     shared by the formal solution and the rate equations (a one-ulp mismatch
     stalls hydrogen at 6e-3);
  3. intensities scaled per wavelength (SI UV intensities underflow float),
     cell widths instead of heights;
  4. the SE solve and the rate accumulators in double. Float LU fails even
     with refinement and equilibration (Fe33, the Bifrost column).
- Line profiles are 86 % of the memory of a 2-D run. Float storage halves
  them, with populations changed by ≤ 3e-6 (M-1).
- PRD integrals are safe in float in all three PRD modes (≤ 2e-5). Their cost
  is the re-interpolation of J on every call. Folding it into the weights once
  is precision-independent, exact to 1e-9, and about 10× faster (P-8, M-4).
  Angle-dependent PRD is serial and I/O-bound (P-14).
- Speed: float gains 1.0–1.3× in RH's present scalar loops and 2–2.6× in
  vectorized ones. The per-transition loops do not vectorize today (P-7).
  Restructuring matters more than precision on a CPU; on a GPU float is
  required (F-4).
- Full Stokes (report §5, local branch `exp/float-stokes`): the DELO-Bezier
  solver in float loses up to 2e-5 of I_max on ~1000-point rays, because the
  attenuation of optically thin steps is rounded. Increments with a float-float
  (Kahan) sum of I bring it to ≤ 4e-7, with Q, U, V errors proportional to the
  signal. Batched over 16 wavelengths the float kernel is 4.4× faster than the
  scalar double one.
- Step 1 (float profiles and background file) is done: #230.
- Step 2 (PRD interpolation folded into the weights, angle-dependent PRD
  threaded) is done: #232.
- Step 3 (per-transition loops vectorized, -O3 by default) is done: #238.
  Whole runs 1.15–1.23× faster, results byte-identical.

### 2.5 Solver comparison (S-5, 2026-11-09)

The experiment of SOLVER.md S-5, with the criteria and the calibration
fixed before the runs: the 1-D case set at 1 thread, `ITER_LIMIT` 10⁻⁶,
`N_MAX_SCATTER` 50, `N_MAX_ITER` 3000, modernize 0ec2849 (with S-3, #284),
errors against references converged to 10⁻¹⁰ at f5e9847 (the ne cases with
ne in the state, N-36). Solvers:

- default Ng: the shipped `NG_ORDER` 2, `NG_DELAY` 10, `NG_PERIOD` 3 and the
  PRD early exit of the stress keyword files;
- plain G: `FIXED_POINT`, no Ng;
- Ng on G: `FIXED_POINT` with Ng 2/10/3;
- Anderson: `ANDERSON` at depth m = 3, 5, 10 (LOG variables, GLOBAL
  coefficients, safeguard 1, delay 2).

Work in formal-solution equivalents (all wavelength solves / Nspect;
`Background()` calls matter only for solvene, where they track the
iterations), with G evaluations or iterations in parentheses:

| Case | Role | Default Ng | Plain G | Ng on G | Anderson m=3 | m=5 | m=10 |
|---|---|---|---|---|---|---|---|
| `falc_h6_ca_mg` | easy | 163 (74) | 435 (172) | 213 (75) | 165 (54) | 144 (45) | 130 (39) |
| `falc_h20` | easy | 54 (44) | 145 (135) | 54 (44) | 67 (57) | 49 (39) | 40 (30) |
| `falc_mgi66` | easy | 112 (24) | 112 (47) | 112 (24) | 89 (28) | 83 (25) | 82 (24) |
| `falc_solvene` | easy, ne | 148 (69) | 517 (243) | 519 (244) | 174 (67) | 147 (54) | 121 (40) |
| `falxco_h6_ca` | slow | 619 (421) | 501 (242) | 405 (192) | 193 (84) | 140 (57) | 131 (52) |
| `t5750_h6` | failing | not conv. (3000) | 2879 (1517) | 2330 (1225) | 524 (264) | 432 (215) | 300 (145) |
| `shock_h6` | failing | SIGFPE at 10 | 517 (249) | SIGFPE at 9 | 556 (270) | 580 (282) | 320 (146) |
| `falc_helarge` | failing | SIGFPE at 14 | 416 (365) | SIGFPE at 14 | 144 (93) | 143 (92) | 139 (88) |
| `falc_civi` | failing | SIGFPE at 12 | 256 (179) | SIGFPE at 12 | 129 (61) | 133 (55) | 151 (51) |
| `bf_c100_prd` | failing | not conv. (3000) | 441 (211) | 131 (58) | 99 (42) | 89 (37) | 77 (31) |
| `shock_h6crd_ne` | ne, report | not conv. (3000) | NaN at 17 | NaN at 33 | 146 (95) | 118 (65) | 109 (57) |

Error against the reference (populations at levels above 10⁻⁶ of the
total / emergent spectrum, relative to the peak) and SE residual ‖F_SE‖
(max, levels above 10⁻⁶), for the converged runs:

| Case | Ng on G | Anderson m=5 | Anderson m=10 |
|---|---|---|---|
| `falc_h6_ca_mg` | 3e-06 / 5e-06 / 2e-08 | 5e-06 / 6e-06 / 3e-07 | 4e-06 / 5e-06 / 3e-07 |
| `falc_h20` | 6e-06 / 7e-06 / 2e-07 | 7e-06 / 9e-06 / 4e-07 | 4e-06 / 6e-06 / 4e-07 |
| `falc_mgi66` | 5e-06 / 4e-05 / 4e-06 | 6e-07 / 4e-06 / 2e-07 | 3e-07 / 3e-06 / 3e-07 |
| `falc_solvene` | 5e-06 / 1e-05 / 2e-06 | 6e-06 / 1e-05 / 4e-07 | 1e-06 / 2e-06 / 5e-07 |
| `falxco_h6_ca` | 6e-05 / 1e-07 / 3e-07 | 1e-05 / 7e-07 / 1e-07 | 1e-05 / 1e-06 / 8e-08 |
| `t5750_h6` | 5e-05 / 5e-05 / 1e-07 | 5e-05 / 3e-05 / 5e-08 | 5e-05 / 3e-05 / 7e-08 |
| `shock_h6` | — | 3e-06 / 3e-06 / 3e-07 | 5e-06 / 7e-06 / 2e-07 |
| `falc_helarge` | — | 2e-06 / 3e-06 / 3e-07 | 3e-06 / 2e-06 / 2e-07 |
| `falc_civi` | — | 6e-07 / 5e-06 / 1e-07 | 6e-06 / 6e-06 / 3e-07 |
| `bf_c100_prd` | 3e-06 / 2e-06 / 4e-07 | 6e-06 / 4e-06 / 4e-07 | 3e-06 / 1e-06 / 3e-07 |
| `shock_h6crd_ne` | — | 4e+01 / 3e+00 / 5e-07 | 4e+01 / 3e+00 / 7e-07 |

All converged runs are within 10⁻⁴ of their reference. `shock_h6crd_ne`
has no converged reference (plain G diverges), so its error column is
meaningless. Anderson's solutions satisfy SE to 5–7×10⁻⁷ with positive
populations, and m = 3, 5 and 10 agree to 10⁻⁵. The default Ng's solvene
run, judged against the default-solver reference, is within
8.8×10⁻⁶ / 1.9×10⁻⁵.

Ablations of Anderson m = 5 (work, with iterations and rejected points):

| Case | m=5 (LOG, GLOBAL, c = 1) | RELATIVE | PER_BLOCK | safeguard off |
|---|---|---|---|---|
| `falc_h6_ca_mg` | 144 (45), 0 rej. | 174 (58), 0 rej. | 300 (113), 48 rej. | 144 (45), 0 rej. |
| `falxco_h6_ca` | 140 (57), 0 rej. | 617 (301), 106 rej. | 388 (184), 82 rej. | 140 (57), 0 rej. |
| `shock_h6` | 580 (282), 78 rej. | 442 (210), 80 rej. | 600 (293), 138 rej. | 1721 (881), 0 rej. |
| `falc_helarge` | 143 (92), 3 rej. | 254 (203), 51 rej. | 190 (139), 23 rej. | 134 (83), 0 rej. |

Anderson's counters (rejected by the safeguard / guard damped / guard plain
step / rate equations solved with LU instead of GTH):

| Case | m=5: rejected / guard damped / guard plain / LU | m=10 |
|---|---|---|
| `falc_h6_ca_mg` | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| `falc_h20` | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| `falc_mgi66` | 0 / 17 / 6 / 0 | 0 / 17 / 6 / 0 |
| `falc_solvene` | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| `falxco_h6_ca` | 0 / 1 / 0 / 0 | 1 / 1 / 0 / 0 |
| `t5750_h6` | 27 / 18 / 5 / 0 | 18 / 9 / 2 / 0 |
| `shock_h6` | 78 / 0 / 0 / 1 | 31 / 0 / 0 / 1 |
| `falc_helarge` | 3 / 105 / 44 / 0 | 3 / 105 / 44 / 0 |
| `falc_civi` | 0 / 36 / 12 / 0 | 1 / 31 / 11 / 0 |
| `bf_c100_prd` | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| `shock_h6crd_ne` | 1 / 13 / 5 / 49 | 1 / 13 / 5 / 49 |

Change per evaluation on the slow case `falxco_h6_ca`, every 10th
evaluation:

| Solver | Evaluations | Change at 1, 11, 21, … |
|---|---|---|
| plain G | 242 | 4e0 2e-1 9e-2 6e-2 4e-2 3e-2 2e-2 9e-3 6e-3 3e-3 2e-3 1e-3 6e-4 3e-4 2e-4 1e-4 7e-5 4e-5 2e-5 1e-5 9e-6 5e-6 3e-6 2e-6 1e-6 |
| Ng on G | 192 | 4e0 1e-1 3e-1 9e-2 3e-2 8e-2 2e-2 2e-2 5e-3 2e-3 9e-4 2e-3 5e-4 1e-4 2e-4 6e-5 2e-5 3e-5 2e-5 3e-6 |
| Anderson m = 5 | 57 | 4e0 7e-2 1e-2 1e-3 9e-5 8e-6 |
| Anderson m = 10 | 52 | 4e0 7e-2 3e-2 1e-3 5e-5 1e-6 |

Wall times were measured with up to 27 runs at once on 48 cores and are
not reported; the work counts are the measure. Data, tables of every run
and every metric: `claude-tmp/s1/s5` (`table.md`). Analysis and decisions:
SOLVER.md S-5, Results.

---

## 3. Issues

### 3.1 Bugs and crashes

| ID | Sev | Issue | Where | Evidence | Status |
|---|---|---|---|---|---|
| [B-1](https://github.com/dwijn/RH/issues/2) | high | `SOLVE_NE = ITERATION` aborts on the second `Background()` call ("Reached end of datafile"): collision data re-read from an open file without rewinding; `offset_coll` saved but never used; `atom->C` reallocated (leak) | `ltepops.c:239`, `readatom.c:568-569`, `collision.c:482` | crashes upstream and WIP | fixed |
| [B-2](https://github.com/dwijn/RH/issues/3) | high | `PRD_NG_ORDER > 0` with default `PRD_NG_PERIOD = 0` divides by zero | `readinput.c:99`, `accelerate.c:88` | SIGFPE | fixed |
| [B-3](https://github.com/dwijn/RH/issues/4) | high | Collisional relaxation (`COLL_RELAX_NSTEP/FACTOR`) stops after one iteration with LTE populations: convergence test fires during the relaxation phase while collisions are scaled by up to 1e38 | `iterate.c:132`, `fillgamma.c:59` | FALXCO H_6: 1 iteration, delta 1.4e-12, LTE result | fixed |
| [B-4](https://github.com/dwijn/RH/issues/5) | high | `CaII_large` crashes at setup (SIGFPE in `splineCoef` from `SortLambda`) | `sortlambda.c:387` | upstream and WIP | fixed |
| [B-5](https://github.com/dwijn/RH/issues/6) | medium | Linear short-characteristics end point uses stale weights (`w2` not recomputed) | `rhf1d/piecewise_1D.c:126` | flux changes up to 1e-3 when fixed | fixed |
| [B-6](https://github.com/dwijn/RH/issues/7) | high | The preconditioned rate equations converge to the wrong populations wherever active transitions overlap other than by sharing an upper level. `Ieff` subtracted Psi times the emissivity of all active transitions, but Gamma put back only part of it: pairs of levels without a transition were missing, one of the two chi terms was missing on connected pairs, and the term for a lower level that is also an upper level was added once per transition sharing it (the multi-ion double counting first suspected). The converged populations then violate statistical equilibrium with their own radiation field. Lightweaver's form avoids the double count but has the same fixed-point error | `fillgamma.c` addtoGamma, `rhsphere/formal.c` addtoGamma_sphere | residual of the plain rate equations after convergence to 1e-10 (FAL-C, CRD, no polarization): 1e-7 Ca II (H and K wings), 7e-4 H_6, 3e-2 MgI+II; with the fix 1e-11 to 1e-9. MgI+II 4s and 3d populations move 2-4 % (median), up to 22 % (34 % with PRD h and k). `tests/unit/test_coupling` | fixed |
| [B-7](https://github.com/dwijn/RH/issues/8) | medium | 3-D: `Jgas` accumulated with `wmu` inside the to_obs loop, looks double counted | `rhsc3d/formal.c:370-374` | code reading | open |
| [B-8](https://github.com/dwijn/RH/issues/9) | medium | `calcDistance()` falls off the end without a return value | `rhsc3d/geometry.c:265` | -Wreturn-type | fixed |
| [B-9](https://github.com/dwijn/RH/issues/10) | medium | `pread`/`pwrite` results compared as `int` against `size_t`; -1 error mishandled | `readj.c:58,82,146,171,200,229,413,442` | -Wsign-compare | fixed |
| [B-10](https://github.com/dwijn/RH/issues/11) | medium | 32-bit overflow of `Nlambda*Nspace` sizes (`NgInit(int N)`, `Nlamu`) once > 2³¹ (e.g. >163 PRD wavelengths at 13 M points) | `accelerate.c:34`, `redistribute.c` | code reading | open |
| [B-11](https://github.com/dwijn/RH/issues/12) | low | Leaks: `backgrflags`/`backgrrecno` re-malloc'd every `Background()`; `spectrum.I` allocated twice in 3-D with Stokes | `background.c:300,312`, `initial_xdr.c:386,392` | code reading | fixed |
| [B-12](https://github.com/dwijn/RH/issues/13) | low | AR85-CHP collisions hard-code `atmos.H->n[5]` as protons (wrong unless H has 6 levels; no shipped atom uses it) | `collision.c:883` | code reading | open |
| [B-13](https://github.com/dwijn/RH/issues/14) | low | `ERROR_LEVEL_1` for an unknown solver leaves I and Psi undefined and continues. **Not a bug:** `ERROR_LEVEL_1` terminates the run (see Q-16), and the keyword parser already rejects unknown solvers | `rhf1d/formal.c:209-231` | code reading | wontfix |
| [B-14](https://github.com/dwijn/RH/issues/15) | low | `error.c` `sprintf`s a 512-byte message plus routine name into a 512-byte buffer | `error.c:51` | code reading | open |
| [B-15](https://github.com/dwijn/RH/issues/80) | low | `Atoms/CaII_large.atom` has 13 `OMEGA` collision records equal to `Inf`; since B-4 the run stops at setup with a clear error instead of a SIGFPE | `Atoms/CaII_large.atom` | `1d_falc_cail`: "collision record OMEGA 11 10 has a non-finite value (inf)" | fixed |
| [B-16](https://github.com/dwijn/RH/issues/93) | medium | Collision data must follow a `TEMP` line, but this is not enforced: a tabulated keyword (OMEGA, CE, CI, CP, CH0, CH+, CH) before any `TEMP` reads an uninitialized item count and a NULL temperature grid (crash); after a fixed-size keyword (AR85-*, BURGESS, SHULL82, BADNELL, SUMMERS) without a new `TEMP` it inherits that keyword's item count, so a line with that many numbers gives silently wrong rates | `collision.c:469` | -Wmaybe-uninitialized triage; no shipped atom file is affected | fixed |
| [B-17](https://github.com/dwijn/RH/issues/94) | low | `writeMULTIatmos()` writes an uninitialized `height` for a GEOMETRIC depth scale (no branch for it); latent, both callers use other scales | `tools/writemulti.c:89` | -Wmaybe-uninitialized triage | fixed |
| [B-18](https://github.com/dwijn/RH/issues/95) | low | `solve1d` reads the boundary-condition codes `btop`/`bbot` from stdin without validating them; an out-of-range code now gives a zero boundary intensity | `rhf1d/solve1d.c` | -Wmaybe-uninitialized triage | fixed |
| [B-19](https://github.com/dwijn/RH/issues/113) | low | `StopRequested()` opens `STOP_RH` with `fopen` to test for it and never closes the file | `stopreq.c:39` | code reading | fixed |
| [B-20](https://github.com/dwijn/RH/issues/128) | low | `VoigtHui()` returns NaN at a = v = 0: the Horner loops multiply numerator and denominator by z once too often, so the ratio is 0/0 at z = 0. Used for angle-dependent PRD (`rii.c`, `riiplane.c`), where a = 0 occurs only with zero damping | `voigt.c:324` | `test_voigt` (Q-14) | fixed |
| [B-21](https://github.com/dwijn/RH/issues/135) | low | `BiLinear()` reads one element past `a_table` and `f` when `a` (or `b`) equals the last table node, extrapolates below the table but not above, and keeps its hunt indices in function statics; only caller is the Voigt `LOOKUP` table | `linear.c:56,75` | `test_interp` (Q-14) | fixed |
| [B-22](https://github.com/dwijn/RH/issues/136) | low | `w3js()` accepts non-integer j−m and j1+j2+j3 and returns nonzero values (`w3js(0.5,0.5,0.5,0.5,−0.5,0)` = 1); `fact()` returns 1 for non-integer arguments; the symbol loses accuracy above j ≈ 40 (1e-9) and returns ±inf for j = 50–56 and 0 from 58 (factorial overflow); only caller is `paschen.c` at small J | `wigner.c:49,167` | `test_functions` (Q-14) | fixed |
| [B-23](https://github.com/dwijn/RH/issues/137) | low | `Hunt()` returns one index too low when hunting downward onto a value equal to a table node, disagreeing with `Locate()`; harmless for linear and spline interpolation | `hunt.c:40,55` | `test_interp` (Q-14) | fixed |
| [B-24](https://github.com/dwijn/RH/issues/149) | high | 2-D and 3-D long characteristics (periodic boundaries, used where a short characteristic crosses column lines or vertical planes before the next layer) stop one section short: `SolveLong`/`SolveLongStokes` treat the last crossing as the local point and return the intensity one crossing upwind of the true upwind point. Linear S on a uniform slab: emergent I 7% (2-D) and 11% (3-D) wrong, identical columns differ by 15%; the shipped uniform 2-D/3-D examples have column spreads of 2.4e-3 and 1e-2 (2e-4 and 1.9e-3 when fixed) | `rhsc2d/solvelong.c:60,140`, `rhsc3d/solvelong.c:58,138` | `test_formal_2d`, `test_formal_3d` (Q-14) | fixed |
| [B-25](https://github.com/dwijn/RH/issues/155) | high | Polarized DELO-Bezier solvers wrong in two control-point terms (1-D, 2-D, 3-D): (i) the derivative terms are paired with the wrong weights, γ S'(c) − θ S'(u) instead of γ S'(u) − θ S'(c), which makes them second order (with K' = 0 they differ from the unpolarized Bezier solver by 2.4% at 81 points); (ii) the sign of dK'/dτ at the central point is wrong, which makes them first order when K varies with depth. Converged Ca II full Stokes on FALC (82 points): Stokes V 38% from the converged answer (3% when fixed), the 854.2 nm V amplitude 21% too large, 1.9x more iterations | `rhf1d/bezier_1D.c:228,251`, `rhsc2d/bezier_2D.c:254,277`, `rhsc3d/bezier_3D.c:151,173` | `test_formal_stokes` (Q-14) | fixed |
| [B-26](https://github.com/dwijn/RH/issues/156) | medium | 2-D and 3-D polarized solvers use a stale or uninitialized source vector in the last layer of each ray: 2-D takes the value of the last point of the layer below, 3-D reads an uninitialized local variable. The emergent Stokes vector (top layer of the outgoing rays) is computed with the wrong source vector; on a test slab the bottom layer of I- was 38% off | `rhsc2d/bezier_2D.c:188,230`, `rhsc3d/bezier_3D.c:85,127` | `test_formal_2d`, `test_formal_3d` (Q-14) | fixed |
| [B-27](https://github.com/dwijn/RH/issues/166) | low | A J = 0 -> 0 line is marked polarizable (`readAtom` only requires abs(Ju - Jl) <= 1) and `Zeeman()` normalizes its single component, of strength 0, by a zero sum: every profile of the line is NaN in a full-Stokes run. Such a line has no single-photon strength at all, but a mislabeled atom file poisons the whole run | `readatom.c:339`, `zeeman.c:290` | `test_zeeman` (Q-14) | fixed |
| [B-28](https://github.com/dwijn/RH/issues/167) | low | `Bproject()` (1-D, 2-D, 3-D) divides by the squared transverse field, 1 - cos^2(gamma), to get cos 2chi and sin 2chi: a field exactly along an oblique ray gives 0/0, and phi_Q, phi_U, psi_Q, psi_U become NaN for every polarizable line at that point (the terms carry a factor sin^2(gamma) = 0). A uniform inclined field along one of the rays triggers it | `rhf1d/project.c:70`, `rhsc2d/project.c:71`, `rhsc3d/project.c:76` | `test_zeeman` (Q-14) | fixed |
| [B-29](https://github.com/dwijn/RH/issues/172) | high | Ca II with PRD in realistic (Bifrost) atmospheres: the PRD sub-iterations diverge (delta 1e4 to 1e31 within a few iterations), the populations go negative, the total opacity passes through zero, S = eta/chi overflows and the run dies with a floating-point exception (in the Bezier solver, `Accelerate`, `addtoRates` or the LU solve, wherever the overflow first lands). Nothing guards against negative populations or opacities. Velocities are the main trigger: in 1-D, 35 of 126 columns of a Bifrost slice (en024048_hion, y = 0, every 4th column) crash; with the velocities zeroed 33 of those converge, with Ca II in CRD 34 (the other one needs more than 200 iterations). The 2-D slice (126 x 89) crashes with angle-averaged PRD, with `PRD_ANGLE_APPROX`, and with vx = 0; static it converges in 44 iterations, in CRD in 25. The 3-D cut (24 x 24 x 89) crashes the same way. Upstream RH crashes the same way on the 2-D slice | PRD iteration (`scatter.c`); no positivity guard after the statistical equilibrium or in `Opacity` | `tests/make_bifrost_atmos.py DIR OUT --1d --x0 260 --zstride 2`, rhf1d with the example CaII.atom: delta 2e31 at iteration 3 | fixed |
| [B-30](https://github.com/dwijn/RH/issues/173) | low | `OH_bf_opac()` and `CH_bf_opac()` pass uninitialised start indices to `Hunt()`. valgrind reports conditional jumps on and uses of uninitialised values in every run with molecules, and the tainted opacities reach the background scratch file (`writeBackground`: pwrite of uninitialised bytes). `Hunt` brackets correctly from any start, so results are unaffected (byte-identical with the indices set to 0); undefined behaviour all the same | `ohchbf.c:54,384` | valgrind on any rhf1d run | fixed |
| [B-31](https://github.com/dwijn/RH/issues/183) | low | Hybrid PRD (`PRD_ANGLE_APPROX`) with a single ray (`Nrays = 1`): `Jgas` is zeroed before each formal solution only when Nrays > 1 (`iterate.c`), but the angle-independent branch of `Formal` adds to it regardless, so it accumulates over iterations and the PRD scattering uses an inflated gas-frame J. Found in the review of #182 by reading the code; not reproduced | `iterate.c:258-266`, `rhf1d/formal.c`, `rhsc2d/formal.c`, `rhsc3d/formal.c` (angle-independent branch) | reading | fixed |
| [B-32](https://github.com/dwijn/RH/issues/184) | medium | The coefficient table of hybrid PRD (`spectrum.nc`, `iprdh`, `cprdh`, built in `initSolution`) is indexed and sized with `int` arithmetic: `2*Nrays*Nspect*Nspace` and the cumulative coefficient counts overflow for large 3-D grids (e.g. 24 rays x 2000 wavelengths x 256 x 256 x 82), giving a too-small allocation and corrupt indices. Found in the review of #182 | `initial_xdr.c:180-230`, `spectrum.h` | reading | fixed |
| [B-33](https://github.com/dwijn/RH/issues/204) | low | Intermittent SIGFPE in `e2e_stokes_1d` (rhf1d, run `rot_40_turned`, an LTE full-Stokes run on a grey atmosphere with `SOLVE_NE = ONCE`, one thread) in the `build (ubuntu-24.04, clang)` CI job of #202, 2026-10-21; the re-run passed. Not reproduced: 60/60 local runs with clang (LLVM 14) pass, valgrind `--track-origins` over all runs of the test reports only an inlined-`strcmp` false positive (`collision.c:565`); no earlier occurrence in the previous 40 failed CI runs. Nondeterminism in a single-threaded run points to an uninitialized read or to the trap state; the e2e script does not print `out.log`, so the fpehandler backtrace was lost. Second occurrence in the `build (ubuntu-24.04, gcc)` job of #211, 2026-10-21 (same run; re-run passed); not reproduced locally with gcc 8.5 (30/30) or gcc 13.4 (20/20), with `MALLOC_PERTURB_`, or under valgrind. Suspect: runner hardware (the Ubuntu runners vary between Intel Xeon 8573C, AMD EPYC 9V74 and EPYC 7763, the last without AVX-512; glibc 2.39 picks libm paths by CPU). #212 adds diagnostics: the e2e scripts print the end of `out.log` on failure and CI prints the runner CPU, flags and libc. **Cause (2026-11-06):** `printTotalCPU()` prints user and system time as percentages of their sum, taken from `times()` in clock ticks of 10 ms; when the whole run uses less than one tick the sum is 0 and 0/0 raises FE_INVALID, which `SetFPEtraps()` turns into SIGFPE. `rot_40_turned` (LTE, `N_MAX_ITER = 0`) is at that threshold on the fastest runners: with ubuntu-26.04, 3 of 4 jobs on AMD EPYC 9V45 (Zen 5) failed, 0 of 14 on other CPUs. Reproduced on sundog2 by preloading a `times()` that reports zero ticks; `addr2line` gives `printTotalCPU` (`getcpu.c:101`) called from `main` | `getcpu.c:101` | CI #202 (run 39647896365, attempt 1), CI #211 (run 39696565002, attempt 1) | fixed |
| [B-34](https://github.com/dwijn/RH/issues/215) | medium | The 2-D and 3-D `solveray` compute the third direction cosine as `sqrt(1 - (a^2 + b^2))` without a clamp or a consistency check: when the given cosines' squares sum to just over 1 (one rounding step, e.g. `mux = sqrt(1 - muz^2)` rounded up in `ray.input`, or an FMA-contracted sum on arm64), it is NaN and the run completes normally. 2-D: `muy` only enters `Bproject()`, so every polarized result is NaN (I, Q, U, V at all line wavelengths, all columns) while the continuum stays finite; 3-D: `muz` itself is NaN (not tested). The input checks test each cosine separately, so `mux = muz = 0.9` also passes; the 3-D `muy` error message names the wrong interval. Suggested fix: error out when the sum exceeds 1 by more than 1e-12, otherwise `sqrt(fmax(0, ...))` | `rhsc2d/solveray.c:98-123`, `rhsc3d/solveray.c:95-119` | 2-D MURaM slab, Fe I 617.3 nm full Stokes, macOS arm64, modernize 810afd7 and upstream: NaN at mu = 0.52, 0.50, 0.49, 0.48, 0.45, 0.447 (with `mux` to 10 decimals) and at 0.44 with a sum of exactly 1.0 in Python; fine with `mux` rounded down | open |
| [B-35](https://github.com/dwijn/RH/issues/221) | low | `rhanalyze` (`python/rhanalyze/rhopacity.py`) returns the wrong record as the background scattering opacity at polarized wavelengths of a Stokes run: after the emissivity block (4 records: I, Q, U, V) it skips 3 records instead of 4, so `scatt` is the Stokes-V emissivity, and the source function it builds, `S = (eta_as + eta_c + J*scatt)/chi`, is wrong there. Unpolarized wavelengths and the IDL reader (`readopacity.pro`, which reads the records in sequence) are unaffected. Found while checking every reader of `background.dat` for P-13 | `python/rhanalyze/rhopacity.py:178` | 1-D Ca II example in full Stokes with `gf0400.10`: at 91 of the 94 polarized wavelengths `scatt` is negative; fixed, it equals the last record of each block and agrees with the unpolarized neighbours to 5 % | fixed |
| [B-36](https://github.com/dwijn/RH/issues/235) | low | `rowcol` (density-dependent dielectronic recombination, collision keyword `SUMMERS`) covers only the first five of the six rows of its periodic-table table: for a recombined ion with 55 to 86 electrons the row is never set, and the density factor is computed from an uninitialized value. Only `C_I+II_9.atom` uses `SUMMERS` among the shipped atoms (at most 6 electrons), so no shipped model is affected. GCC at -O3 flags it as maybe-uninitialized | `collision.c:84` | Found when -O3 became the default (P-7): the loop runs `i < 5` over the seven `istart` entries (rows 1-6) | fixed |
| [B-37](https://github.com/dwijn/RH/issues/241) | low | `chartflow` (2-D flow maps): at wavelengths without angle-dependent opacity the total opacity and source function add `chi_as` and `eta_as`, arrays that are allocated and never written, instead of the active opacity `as->chi` and `as->eta` that the angle-dependent branch uses. The result depends on what `malloc` returns; fresh pages are zero, so the active-line opacity and emissivity are silently left out. GCC 15 at -O3 reports the uninitialized read | `rhsc2d/chartflow.c:227` | Found by CI (GCC 15, -O3, -Werror) on #238; moving 2-D small case and static `FALC_03x82`, 41 wavelengths each: flux files unchanged by the fix | fixed |
| [B-38](https://github.com/dwijn/RH/issues/243) | low | Active molecules: `initGammaMolecule` resets only the Δv = 1 entries of `Gamma` (the fixed collisional rates), while `addtoGamma` and the thread merge add the radiative rates of every vibration-rotation line with `+=`. The entries of lines with Δv ≥ 2 (overtone bands, e.g. the CO first overtone) are never reset, so their radiative rates accumulate from one iteration to the next and the populations solve a rate matrix that grows with the iteration count. The atom version (`initGammaAtom`) sets every entry. The shipped CO line lists contain only Δv = 1 lines, so no shipped case is affected; any list with overtone lines is. Also a hidden state of the iteration that S-2 must remove | `fillgamma.c:70-88` | Code reading (SOLVER.md audit, 2026-11-07): the only writes to `molecule->Gamma` are `initGammaMolecule`, `addtoGamma` (`fillgamma.c:323`) and the merge (`fillgamma.c:171`); `statEquilMolecule` works on a copy and rebuilds only the diagonal. The GOORVITCH94 and Kurucz molecular line formats accept any vi, vj (`readmolecule.c:668-837`); `Molecules/CO/*_dv=1_*` hold Δv = 1 only. Same code upstream | fixed |
| [B-39](https://github.com/dwijn/RH/issues/244) | medium | `Background()` is not safe to call more than once: every call opens `background.dat` again (`O_TRUNC`) without closing the previous descriptor, reads the Kurucz line lists into a new array (`readKuruczLines` sets `atmos.rlk_lines = NULL` when `Nrlk == 0`, dropping the old one), and re-reads the Barklem tables and rewrites `metals.out`, `molecules.out` and `brs.out`; B-11 is the same for `backgrflags`/`backgrrecno`. With `SOLVE_NE = ITERATION` that is one descriptor and one line list per iteration. S-2 with ne in the state (S-2b in SOLVER.md) calls `Background()` inside G, once per Jacobian-vector product, thousands of times: the common descriptor limit of 1024 ends the run, and large line lists are re-read every time. Not reachable upstream, where the second call aborts (B-1) | `background.c:723-760`, `kurucz.c:172` | strace of `1d_falc_solvene` (b5aa20c, 27 iterations): 27 `openat` calls of `background.dat` with `O_TRUNC` returning descriptors 5, 12, 13, …, 37, none closed; `kurucz.input`, the three Barklem tables and the `.out` files are opened 27 times each | fixed |
| [B-40](https://github.com/dwijn/RH/issues/247) | low | The H⁻ density stored by `ChemicalEquilibrium` (`atmos.nHmin`, used for the H⁻ bound-free opacity and written to `molecules.out`) is `ne * n_H * Phi_H-` with `n_H` all atomic hydrogen, neutral atoms plus protons, while the hydrogen conservation equation solved a few lines earlier uses the neutral part, `ne * fn0 * n_H * Phi_H-`. The stored value is therefore too large by 1/fn0 and inconsistent with the solution of the chemical equilibrium. Where hydrogen is mostly neutral (photospheres of cool stars) the difference is negligible; it grows where H⁻ still matters and hydrogen is partly ionized (hotter stars, deep layers) | `chemequil.c:336` | Code reading (SOLVER.md audit, 2026-11-07), confirmed by a build with `fn0[0]` in the stored value: on FALC_82 the ratio fixed/current is the neutral fraction, 1e-11 at the top, 0.72 at depth 42, 1 through the photosphere and 0.97 at the bottom (T5750: 0.95 at the bottom). The layers with fn0 < 0.9 hold 3e-9 of the H⁻, so the emergent spectra of `1d_falc_h6_ca_mg`, `1d_fala_h_ca` and `1d_t5750_h6` change by at most 5e-4 (median 3e-5), the level of convergence noise at ITER_LIMIT = 1e-3. Same code upstream | fixed |
| [B-41](https://github.com/dwijn/RH/issues/255) | low | A terminating error message overflows a stack buffer when it is long: `ErrorExit` formats it with `sprintf` into `MAX_MESSAGE_LENGTH` (512) bytes, about 50 fixed characters plus the routine name and the message, so a routine name plus message longer than about 460 characters writes past the buffer. No message in the current code comes close (the longest found, with a 160-character keyword value, is ~260 characters), so the defect is latent; it surfaced while writing the S-2 keyword-conflict messages, which are now one warning per conflict plus a short error. `messageStr`, into which callers `sprintf`, is also 512 bytes | `error.c:46-50` | Code reading (S-2, 2026-11-08) and a scan of every `sprintf(messageStr, ...)` and `sprintf(errorStr, ...)` format for unbounded `%s` arguments. Same code upstream | open |
| [B-42](https://github.com/dwijn/RH/issues/258) | medium | `Atoms/C_I+II_9.atom` has two corrupt C II fine-structure collision strengths: in `OMEGA 4 3` (2p ²P°1/2–²P°3/2, the [C II] 158 µm line) and `OMEGA 7 6` (2s2p² ²D5/2–²D3/2) the value at 9000 K is `9.000E+03`, the node temperature copied into the data column, between neighbours of 2.08/2.12 and 1.91/1.90. The interpolant (natural cubic spline in T) passes through 9000 at 9 kK and, since a spline is not local, rings across the table: it is negative in six intervals between 3 and 15 kK (down to −1228) and more than 10% off the spline through a corrected value from 1.1 to 19 kK, so these collisional rates are wrong through the whole chromospheric range, up to 4000 times too large or negative. `OMEGA 4 3` also has 2.656 at 18 kK between 2.247 and 2.270, probably a typo for about 2.26 | `Atoms/C_I+II_9.atom:364-365` | Scan of the collision tables of all 88 atom files for records whose interpolant (`collision.c:729`) goes negative, and for isolated spikes, 2026-11-08: three bad tables (B-42, B-43, B-44), every other hit is a copy of one of them. Not used by any test or example. Same file upstream | fixed |
| [B-43](https://github.com/dwijn/RH/issues/259) | medium | The Ca II collision strengths in `CaII.atom` and its copies (`CaII_CRD`, `CaII_XRD`, `CaII_comp`, `CaII_fixed`, `CaI+II`) disagree with the R-matrix values of Meléndez, Bautista & Badnell (2007, A&A 469, 1203). 3d ²D3/2–²D5/2: 1.78, 2.87, 4.69, 12.66 at 3, 5, 7, 15 kK against 15.6, 15.3, 15.2, 14.4, then 3.19 and 7.93 at 50 and 100 kK (Meléndez stops at 38 kK, where it is 11.6). 4p ²P°1/2–²P°3/2: 0.85 to 4.0 at 3 to 15 kK against a nearly flat 5.3 to 5.5. 4s–3d: half of Meléndez's. 4s–4p and 3d–4p: within 25%. The file gives no source. The non-monotone 3d–3d row makes the interpolated Ω negative between 50 and 100 kK (−0.60 at 67 kK), so the collisional rates between the two 3d levels are negative at the top of FAL-C (k = 5–7, 60–71 kK). The larger models (`CaII_large`, `CaII_medium`, `CaII_med_oldbf`) give 3d–3d a constant 2102 (van Regemorter), 140 times Meléndez's. Effect on results not measured yet | `Atoms/CaII.atom:154-161` and copies | `1d_falc_h6_ca_mg` with each rate-matrix entry split into its collisional, spontaneous, Ieff and preconditioning terms (GTH test of SOLVER.md §9 Q12, 2026-11-08): C(1→2) = −50 s⁻¹ and C(2→1) = −34 s⁻¹ at k = 5, the only negative off-diagonal entries in H, Ca II and Mg II over 40 iterations. Meléndez's Υ from VizieR J/A+A/469/1203. `CaII.atom` is active in the default `rhf1d` run and in 16 stress cases. Same files upstream | fixed |
| [B-44](https://github.com/dwijn/RH/issues/260) | low | Nine of the 45 Mg II collision strengths in `MgII.atom` (and the copies in `MgII_CRD`, `MgI+II`, `MgI_b+MgII`) disagree with Sigut & Pradhan (1995, J. Phys. B 28, 4879), the source the file names. 3d–4f (`OMEGA 9 4`, `OMEGA 9 5`): ten times too small at 1–10 kK and a hundred times at 20 and 30 kK, where the values also drop tenfold (3.88, 4.25, 4.57, 5.23, 0.624, 0.708 for ²D5/2 against 38.8, 42.5, 45.7, 52.3, 62.4, 70.8); the interpolant went negative between 20 and 30 kK (−0.06 at 24 kK). The tenfold drop looked like an exponent typo for 6.24 and 7.08, but Table 5 gives 3d–4f as 64.7 to 118, so the whole row is a factor 10 low. 4s–3p, 5s–3p and 5s–3d: both J components carry the full LS value instead of splitting it in proportion to g(J), so they are 3 and 1.5 (3p1/2, 3p3/2) and 1.67 and 2.5 (3d5/2, 3d3/2) times too large. 3s–3p3/2 at 10 kK: 12.3 for 11.3. The other 36 rows agree within 3%. Fixed by rebuilding all 45 rows from S&P: fine-structure pairs from Table 6, a J level to an unsplit term from the sum of the Table 6 components where it has them, otherwise the Table 5 LS value times g(J)/g(term) | `Atoms/MgII.atom:286-344` and copies | Same scan as B-42; S&P Tables 5 and 6 transcribed from the paper and checked against it, 2026-11-12. Effect: `1d_falc_h6_ca_mg` and `1d_falc_nong` converge in the same number of iterations (20, 72), Mg II populations change by up to 4.7%, H and Ca by at most 4e-4 (convergence noise). On FAL-C with a 330-line Mg I+II atom and NLTE Fe, the Mg II k core irradiance drops 4.5%, h 0.3%, the subordinate lines < 0.4%. `MgII.atom` is active in 5 stress cases. Same files upstream | fixed |
| [B-45](https://github.com/dwijn/RH/issues/298) | medium | A negative `NG_ORDER` (or `PRD_NG_ORDER`) is accepted and crashes the run: with -2, `Accelerate` takes `count % (Norder + 2)`, a modulo by zero (SIGFPE at the first iteration); other negative values give `matrix_double` a negative or zero row count, and `PRD_NG_ORDER = -2` corrupts the heap (abort, "corrupted double-linked list", in the first iteration). Upstream the same input ends in a terminating error from `readBackground` ("Error reading file"), a symptom of the corrupted allocation. readInput checks `NG_PERIOD` (B-2) but not the sign of the orders or delays | `readinput.c` (Ng checks), `accelerate.c:74` (`NgInit`, `Accelerate`) | `rhf1d` example with `NG_ORDER = -2`, `NG_DELAY = 0`: SIGFPE here (status 136), "Error reading file" upstream (c95bdf6), where `NG_ORDER = 2` runs; `PRD_NG_ORDER = -2`: abort (status 134); found in review of #296, 2026-11-10 | fixed |
| [B-46](https://github.com/dwijn/RH/issues/299) | low | Under the G solvers (`FIXED_POINT`, `ANDERSON`, `NEWTON_KRYLOV`) readInput refuses `PRD_N_MAX_ITER` < 1 ("the map needs at least one PRD sub-iteration"), but `PRD_N_MAX_ITER` < 1 is RH's documented switch for treating every line in CRD: `readAtom` then marks no line PRD, and G has no PRD sub-iterations. Valid input stops the run | `readinput.c` (G solver checks, S-2) | `readatom.c:229`; found while implementing S-10 step 1, 2026-11-10. Exists only here (S-2) | fixed |
| [B-47](https://github.com/dwijn/RH/issues/313) | low | `tests/run_stress.sh` runs nothing and exits 0 when its selection is empty, and says nothing about why: `--only` naming a case outside `--set` (default `conv`, so `--only 2d_bf_s1_perf` alone selects nothing), a bifrost case without `--bifrost`/`RH_BIFROST_DIR`, or a misspelled `--set` value. The report table comes out empty | `tests/run_stress.sh:149` (`selected`) | A relaunch of the S-11 perf runs on 2026-11-11 started without `--set perf` and ended at once with an empty table, costing a restart cycle. Exists only here (run_stress.sh is not upstream) | open |
| [B-48](https://github.com/dwijn/RH/issues/316) | medium | `rhanalyze` imports `xdrlib`, deprecated in Python 3.11 and removed in 3.13. On 3.13 and later it switches to `mda_xdrlib`, an outside package (one or two maintainers) that RH declares nowhere, so on a fresh 3.13+ install every reader module fails at import. The switch also tests only the minor version (`sys.version_info[1] < 13`, wrong from 4.0 on). CI runs 3.12, so the 3.13+ path is never exercised | `python/rhanalyze/` (import block of 7 modules: rhatmos, rhatom, rhgeometry, rhinputs, rhopacity, rhspectrum, rhtools) | Found while planning the move to Python 3.14, 2026-11-12. The API used is small: `Unpacker` (unpack_int/float/double/string/fstring/farray, get_position, reset, done, get_buffer) and `Packer` (pack_int/float/double/farray, get_buffer). Same code upstream (c6b5e2a) | open |
| [B-49](https://github.com/dwijn/RH/issues/322) | medium | `rhanalyze.rhout` cannot read any 3-D output: `rhatmos` matches the geometry type against `"THREEE_D_PLANE"` (three E's) in `molecule.read`, `molecule.read_populations` and `atmos.read_Hminus`, so `dim` is never set and reading `molecules.out` (written by every 3-D run) stops with `UnboundLocalError: ... 'dim'` | `python/rhanalyze/rhatmos.py:40,66,176` | `rhout` on the regression output `3d_small` (test harness, 2026-11-12): UnboundLocalError with the old (`xdrlib`) and new (B-48) readers alike; with the spelling corrected the whole output reads (13 molecules, `Hminus` 4×4×82, `I` 4×4×24×286). Same code upstream (c95bdf6) | open |
| [B-50](https://github.com/dwijn/RH/issues/323) | medium | `input_atmos.write(outfile, Bfile)` ignores its `Bfile` argument and writes the magnetic field to `self.Bfile`, the file the atmosphere was read from: the input B file is overwritten, and the requested one is never written. An atmosphere read without a B file writes none, whatever `Bfile` says | `python/rhanalyze/rhatmos.py:456-466` (`input_atmos.write`) | `input_atmos("ONE_D_PLANE", FALC_82.atmos, Bfile=in.B)`, `B` doubled, `write(outfile=out.atmos, Bfile=out.B)`: `in.B` changed (sha256 0923ac0f… → 2d8675b2…), `out.B` missing (2026-11-12). Found in the review of #318. Same code upstream (c95bdf6) | open |
| [B-51](https://github.com/dwijn/RH/issues/324) | low | `rhanalyze.rhout` cannot read the output of a run with an active molecule. `molecule.read` uses `rhdir`, which is not defined there (NameError). Past that, the reader does not follow the files the C code writes: `molecules.out` holds the level energies `EvJ` (Nv×NJ) after `n` (`writemolec_xdr.c`), which `molecule.read` leaves unread, so every later molecule and `Hminus` are read out of place; `read_populations` reads `EvJ` from `pops_mol.<ID>.out`, which starts instead with the atmosphere ID, Nv and Nspace (`xdr_populations`), and runs out of data (EOFError). It also uses `self.NV` (no such attribute) in the spherical case and stores the populations in `self.Nv`, overwriting the level count | `python/rhanalyze/rhatmos.py:24-77` (`molecule`) | `1d_crd` with `CO_NLTE.molecule` ACTIVE (2026-11-12): NameError; with `rhdir` passed, EOFError in `read_populations`. Read by hand in the C layout, both files end exactly (`done()` passes) and the CO populations are sound (zero above Tmax, ≈ LTE deep down). Passive molecules have Nv = 0 and are not affected. Found in the review of #318. Same code upstream (c95bdf6) | open |
| [B-52](https://github.com/dwijn/RH/issues/335) | high | 3-D long characteristics of the upward sweep sample S and χ at the horizontal mirror image of their crossings. `setStencil` tabulates the crossings of a long characteristic with the vertical planes at +(μx, μy)·s from the local point for both winds; the `DOWNWIND` stencils, which `ShortChar` uses upwind of the upward (to-observer) sweep, need −(μx, μy)·s, as their `longfracx`, `longfracy` and z fractions already have. Every intermediate point of an upward long characteristic is then read up to tens of columns away from where the ray crosses (0.14–0.88 for dx = 0.016 on a 64 × 64 mesh with layers 0.25 apart, i.e. 9–56 columns). Invisible on a horizontally uniform atmosphere, so the slab tests of B-24 did not see it; it affects the emergent intensity (the upward sweep) of every 3-D run with long characteristics in a structured atmosphere | `rhsc3d/geometry.c:412-430` (`setStencil`) | Probe interpolating the coordinates x, y, z with every long-characteristic stencil of `fillMesh` (five rays, all octants, both winds), 2026-11-13: `UPWIND` exact, `DOWNWIND` mirrored. On a slab with S = s0 + s1 τ + A sin 2π(x + y) and an exact solution, the emergent error drops from 0.34 to 0.076 (N = 8) and from 2.0×10⁻² to 1.1×10⁻² (N = 64) with the fix. Found while testing the N-39 fix. Same code upstream (c95bdf6) | fixed (#344) |
| [B-53](https://github.com/dwijn/RH/issues/336) | medium | 3-D long characteristics cross the vertical planes one column or row off for rays with a negative μx or μy. `setStencil` takes the plane index of a crossing as `round(y / dy)` (`round(x / dx)`) with a `round` of its own, `modf(x + 0.5)`, which truncates toward zero and so returns j + 1 for a negative integer j; S and χ at every crossing of such a ray are read from the neighbouring plane. The local `round` also replaces libm's for the whole `rhsc3d` program (no other caller) | `rhsc3d/geometry.c:518-523, 570-576` | Same probe as B-52: `UPWIND` long characteristics of the three rays with a negative direction cosine off by exactly dx, the others exact. Same code upstream (c95bdf6) | fixed (#343) |

### 3.2 Numerics and convergence

| ID | Sev | Issue | Where | Evidence | Status |
|---|---|---|---|---|---|
| [N-1](https://github.com/dwijn/RH/issues/16) | high | No convergence report: at `N_MAX_ITER` the loop ends silently and output is written as if converged; log has no converged/not-converged line; exit status 0 | `iterate.c:132`, `rhf1d/rhf1d.c:91-131` | all runs | fixed |
| [N-2](https://github.com/dwijn/RH/issues/17) | high | Ng acceleration without safeguards: one global coefficient set for all levels/depths weighted by 1/n; no positivity, residual or conditioning check; history not reset after ne or HSE updates (the population history is reset during relaxation since B-3, the PRD history never), and with `SOLVE_NE = ITERATION` or HSE the populations are rescaled in place after Ng stored them; a zero entry (active molecules outside their temperature range, `NG_MOLECULES`) raises SIGFPE in the 1/x weights (SOLVER.md §2.2) | `accelerate.c:88-122` | He_large Ng jump 0.34 → 20.7; FALC_shock H_6 SIGFPE with Ng, converges in 115 without. Bifrost column 36 (Ca II, angle-averaged PRD, after B-29): the Ng step at iteration 42 sends delta to 38, and two iterations later `exp(-dt)` overflows in `Bezier3_coeffs` (#17) | open |
| [N-3](https://github.com/dwijn/RH/issues/18) | high | Convergence criterion is not an accuracy guarantee: max relative change over all levels/depths, dominated by trace levels, measures Ng jumps, ignores PRD ρ, background J and ne | `maxchange.c:43`, `iterate.c:132` | Baseline 2026-10-15 (upstream): He_large "converged" to 1e-3 at 1 thread (108 iterations) and 4 threads (48) differ by 1.4 %; MgI+II at 1e-6 vs 1e-9 differ by 2e-5 | open |
| [N-4](https://github.com/dwijn/RH/issues/19) | high | PRD ρ iteration is Λ-iteration with populations fixed; early exit at `PRD_ITER_LIMIT` makes the map discontinuous → limit cycles; ρ never clipped (ρ < 0 possible after Ng) | `redistribute.c:74-105`, `scatter.c:346,765` | H_6 FALC at 1e-8 stalls at 3–9e-4 for 600 iterations; converges in 142 with fixed 3 sub-iterations; FALXCO H_6 (PRD) no convergence, H_6_CRD converges in 17. Bifrost cuts, Ca II, angle-averaged PRD, after B-29 (#19): 32 of 126 1-D columns end in a periodic limit cycle at 3e-3 to 1e-2 (the PRD sub-iteration exits after one step, its delta 6e-3 being below its 1e-2 limit); on 8 of them `PRD_N_MAX_ITER = 10` with `PRD_ITER_LIMIT = 1e-4` converges 8/8 in 22-31 iterations (old rho: 6/8, 2 FPE), either setting alone 6-7/8, Ng off 0/8. 2-D slice 1.3e-2 after 200 iterations, 3-D 24×24×89 2.6e-2 after 150; with hybrid PRD both converge | open |
| [N-5](https://github.com/dwijn/RH/issues/20) | high | Parabolic short-characteristics Ψ* unbounded (≈ 1 + 1/Δτ_d when Δτ_u ≫ Δτ_d); fallback to linear only on I < 0. `S_PARABOLIC` is the shipped 2-D default | `rhf1d/piecewise_1D.c:234`, `rhsc2d/piecewise_2D.c:309`, `rhsc2d/run/keyword.input:161` | Ψ* up to 12 (> 4000 before blow-up), 11 % of points outside [0,1]; CaII/MgII on FALC_v10 SIGFPE in 3–8 iterations; Bezier/linear converge in 14–28 | fixed |
| [N-6](https://github.com/dwijn/RH/issues/21) | high | No ne preconditioning: `SOLVE_NE` couples ne by plain iteration (dfjk = 0 for NLTE elements); Newton tolerance 1e-2 | `solvene.c:40,165` | (with B-1 patched) H_6_CRD period-2 oscillation at 0.22 after 80 iterations; 2-cycle on FALC_shock | open |
| [N-7](https://github.com/dwijn/RH/issues/22) | medium | Rate-matrix conditioning: LU with row scaling only; iterative improvement in the same precision; LTE relative to the neutral ground level (overflow for high stages); n* clipped to DBL_MIN | `ludcmp.c:73,142`, `ltepops.c:98` | CI–VI Judge (115 levels): 230 negative populations after 4 iterations, b up to 1e307, divergence, SIGFPE | open |
| [N-8](https://github.com/dwijn/RH/issues/23) | medium | No positivity enforcement or damping anywhere | — | grep | open |
| [N-9](https://github.com/dwijn/RH/issues/24) | medium | Starting solutions: `ZERO_RADIATION` has no photoionization (over-recombined start); no J = B or escape-probability option | `initial_xdr.c:21,588-640` | LTE start 40 iterations vs 22 (H_6 FALC) | open |
| [N-10](https://github.com/dwijn/RH/issues/25) | high | Multi-D atmospheres with a horizontally corrugated transition region do not converge: 3-D diverges and crashes, 2-D with the same gradients locks into a period-3 cycle | — | 3-D 8×8×82, H_6+CaII, upstream and WIP. Unperturbed and velocity-only variants converge. With the height/temperature perturbation alone it diverges (SIGFPE) with bicubic or linear 3-D interpolation, Bezier or linear S, PRD or CRD atoms, Ng on or off, and even with 4× gentler gradients (crash at iteration 28). 2-D with equally steep gradients: no crash, deltas cycle 23 → 1 → 2 with period 3 = `NG_PERIOD` for 40 iterations; 4× gentler 2-D converges slowly (1e-2 after 40). Cause not yet isolated: candidates are the 3-D short-characteristics Ψ*, Ng (N-2) and the lack of positivity safeguards (N-8). On `modernize` after #187 (2026-10-20) it still diverges from iteration 1 and overflows at iteration 5 (delta 1.7e6). Unchanged on 2026-11-08 (56c9aec and 6d1eeff, before and after B-40): the largest change of H is 1.4e3 at iteration 1 and 3.3e5 at iteration 4, and an overflow in `SolveLong` (`rhsc3d/solvelong.c:77`) stops it at iteration 5 | open |
| [N-11](https://github.com/dwijn/RH/issues/26) | medium | FPE traps on with no handler: any overflow kills the run with no output or diagnostic | `fpehandler.c:37` | all SIGFPE cases | fixed |
| [N-12](https://github.com/dwijn/RH/issues/27) | low | Collisions with H/protons use H populations frozen at setup; spline interpolation of rates can overshoot negative and clamps outside the table | `collision.c:685,726-746` | code reading | open |
| [N-13](https://github.com/dwijn/RH/issues/28) | low | `w3` uses `1 − exp(−x)` instead of `expm1` (rel. error ~1e-9 near x = 5e-4); Bezier Ψ* ignores dS/dS_k in the derivative | `w3.c`, `bezier_*.c` | code reading | open |
| [N-14](https://github.com/dwijn/RH/issues/81) | medium | Results depend on depth resolution: on FALC_v10 (±10 km/s) the linear and Bezier formal solvers disagree by 5–16 % throughout the chromosphere, so the grid does not resolve the velocity gradients; there is no grid refinement and no resolution check | — | CaII on FALC_v10, linear vs Bezier vs parabolic | open |
| [N-15](https://github.com/dwijn/RH/issues/82) | medium | 3-D horizontal interpolation of S, χ and I (`BICUBIC_3D`, the default) is not limited and can overshoot across a corrugated transition region (the 2-D `Quadr` has a monotonic option); possible contributor to N-10 | `rhsc3d/geometry.c:344,495,633` | code reading | open |
| [N-16](https://github.com/dwijn/RH/issues/114) | low | When the iteration stops without converging, ne (`SOLVE_NE = ITERATION`) and hydrostatic equilibrium are updated once more after the last population solve, so the output populations do not belong to the output ne/atmosphere | `iterate.c:150-171` | code review of N-1 | open |
| [N-17](https://github.com/dwijn/RH/issues/131) | medium | Cubic Bezier integration weights (`Bezier3_coeffs`, used by the default solver in all geometries and for Stokes) switch to a Taylor series at dt < 0.05 that stops at dt³: the weights are wrong by up to 1.2e-5 relative just below the switch, and the closed form above it loses digits to cancellation (9e-10 at dt = 0.056); I is off by about 1e-6 for a linear S, which the method should reproduce exactly | `bezier_aux.c:179` | `test_formal_1d` (Q-14) | fixed |
| [N-18](https://github.com/dwijn/RH/issues/132) | low | `w2`/`w3` (linear and long-characteristics weights) use a two-term series for dtau < 5e-4 (relative error up to 4e-8) and closed forms with cancellation just above | `w3.c:23,42` | `test_formal_1d` (Q-14) | fixed |
| [N-19](https://github.com/dwijn/RH/issues/133) | low | The Bezier solvers' diagonal operator Psi = beta + theta leaves out that the control point moves with S(k) (Steffen-limited derivative): it differs from the true dI(k)/dS(k) by up to 28% on a smooth 1-D slab. Changes MALI's convergence rate, not the solution; relevant to the preconditioner of the S track | `rhf1d/bezier_1D.c:437`, 2-D/3-D Bezier | `test_formal_1d` (Q-14) | wontfix |
| [N-20](https://github.com/dwijn/RH/issues/138) | low | Gauss-Legendre weights accurate to only ~1e-14 (P'(x) evaluated at the previous Newton iterate, stop at 3e-14): weights sum short by 2.4e-14 at n = 64, polynomial exactness 5e-13 | `gaussleg.c:49-56` | `test_quadrature` (Q-14) | fixed |
| [N-21](https://github.com/dwijn/RH/issues/139) | low | `SolveLinearEq(..., improve)`: iterative refinement with the residual in plain double gives no gain (Hilbert n = 8: 8.9e-7 with, 1.4e-7 without); a singular matrix with nonzero rows is not detected (zero pivot replaced by 1e-20) | `ludcmp.c` | `test_linalg` (Q-14) | open |
| [N-22](https://github.com/dwijn/RH/issues/140) | low | Special functions: `Planck` uses `exp(x) − 1` (relative error eps/x, 6e-10 at x = 1e-7, radio wavelengths); `E1`/`E2` return 0 above x = 80 and `E2` loses accuracy toward it (1.2e-6) | `planck.c:46`, `expint.c:51,73` | `test_functions` (Q-14) | fixed |
| [N-23](https://github.com/dwijn/RH/issues/143) | medium | In 1-D, every angle-independent wavelength (static atmosphere, no Stokes, no angle-dependent PRD) is solved with second-order Feautrier, whatever `S_INTERPOLATION` says; Bezier3 (4th order for the emergent intensity) is used only with velocities or polarization. The keyword is silently ignored, and static and moving atmospheres are solved with different methods | `rhf1d/formal.c:356-372` | `1d/nostokes` unchanged by the Bezier fix (N-17); `test_formal_1d` | fixed |
| [N-24](https://github.com/dwijn/RH/issues/150) | medium | 2-D and 3-D cubic Bezier solvers are second order: the upwind control point uses the secant slope (S_local − S_uw)/Δτ instead of a derivative at the upwind point (the 1-D solver carries the centred derivative from the previous point). On a uniform slab where every upwind point lies on a layer: emergent error 4.8e-4 at 16 points per decade, where 1-D Bezier gives 9e-7. Several ways to fix it (interpolated derivatives, a derivative from the upwind stencil); to be chosen | `rhsc2d/bezier_2D.c:510`, `rhsc3d/bezier_3D.c` | `test_formal_2d`, `test_formal_3d` (Q-14) | fixed |
| [N-25](https://github.com/dwijn/RH/issues/174) | low | The polarized DELO-Bezier solvers (1-D, 2-D, 3-D) limit the derivatives of the emission vector and of K element by element (Steffen; `cent_deriv_vec`, `cent_deriv_mat`), which does not commute with turning the (Q, U) frame: for a vertical ray on FAL-C (82 points) and a depth-dependent field, turning the field azimuth by 35 deg gives Q and U off the rotated solution by up to 7e-3 of their maximum, V by 2e-4, I by 1e-5; exact (4e-14) with the limiter replaced by the unlimited centred derivative. A discretization error, third order under refinement. A frame-invariant limiter (one factor for the (Q, U) pair of each quantity, from its length) would remove it | `bezier_aux.c:108-133` | `tests/e2e/stokes_1d.py` | fixed |
| [N-26](https://github.com/dwijn/RH/issues/192) | low | The 2-D and 3-D polarized solvers interpolate the upwind and downwind I, S and K component by component with a monotonic (clipping) interpolation, which, like the limiter of N-25, does not commute with turning the reference frame of Q and U, so 2-D and 3-D Stokes results still depend on it where the interpolant clips. 1-D is exact since N-25 | `rhsc2d/bezier_2D.c` (`Quadr`), `rhsc3d/bezier_3D.c` (`Interpolate_3D`) | reading (review of #191) | open |
| [N-27](https://github.com/dwijn/RH/issues/193) | low | Hybrid PRD in a moving atmosphere: `Jgas` is zeroed at every wavelength in every `solveSpectrum`, also in the PRD sub-iterations, which re-solve only the wavelengths of PRD lines. Until the next full formal solution the edge rows of a PRD line miss the gas-frame contributions of neighbouring wavelengths of other lines | `iterate.c:256-262` | reading (review of #188) | open |
| [N-28](https://github.com/dwijn/RH/issues/209) | medium | `He_large` (53 levels, three stages) does not converge to a tight limit, and `CI-VI_Judge` (115 levels, six stages) does not converge at all (stress case: the old Gamma overflows at iteration 123, the B-6 fix reaches 200 iterations diverging). At `ITER_LIMIT` 1e-10, 600 iterations, the old Gamma stalls at 2e-3 with a residual of 1.6 in the plain rate equations (the stress case "converges" at 1e-3 while stalling), and the per-transition preconditioning of the B-6 fix stalls at 2e-3 with a residual of 2-3e-4, with or without Ng. The behaviour flips with small changes: before N-20/N-22 the stress case took 39 iterations on the old Gamma and the B-6 fix crashed at iteration 13 with Ng (oscillating iterates extrapolated to negative populations, N-2); after them, 75 and 81 iterations. Full preconditioning (RH92 eq. 2.19) diverges from iteration 4 with or without Ng (about 7500 negative off-diagonal Gamma elements per iteration). Since #230 (float background file) the stress case `1d_falc_helarge` dies again with Ng: SIGFPE at iteration 13 at 1 thread (Ng jumps to delta 6.6e6 at iteration 11, then `exp(-dt)` overflows in `Bezier3_coeffs` on a negative opacity), at iteration 24 at 4 threads, converged in 126 at 2 threads; 81 iterations before #230, 93 with the background file in double. Without Ng it stalls at 2.3e-3 after 200 iterations. Another flip of the same chaotic trajectory, not a precision problem: the float file changes each opacity by at most 6e-8 | iteration (`iterate.c`, `accelerate.c`), see N-2 | FAL-C He_large, `ITER_LIMIT` 1e-10; `tests/se_residual.py` | open |
| [N-29](https://github.com/dwijn/RH/issues/216) | low | 2-D long characteristics: on a horizontally uniform periodic slab, where all columns should give the same emergent Stokes vector, columns differ by up to 2e-3 of I once rays need long characteristics (tan(theta) dz/dx > 1); residual after B-24, which added to it upstream (1e-2). Just above the threshold a contiguous band of 45-60 of 128 columns deviates, further from it all columns; the unpolarized solver shows the same, so it is the long-characteristics path in general, probably sensitive to where each periodic sweep starts (`lstart`). Which value is accurate is not known. Candidate unit test: a `test_formal_2d` case at a long-characteristics angle requiring all columns to agree to ~1e-9 | `rhsc2d/solvelong.c`, `rhsc2d/geometry.c` (`fillLongChar`, `lstart`), `rhsc2d/bezier_2D.c`; probably `rhsc3d` | 128 x 108 slab from one repeated MURaM column (dx = 32 km, dz = 16 km), full Stokes Fe I 617.3 nm, `solveray`: max column deviation 1.8e-7 at mu = 0.60 (short characteristics only), 3.0e-4 at 0.44, 1.1e-3 at 0.40, 1.9e-3 at 0.20 (upstream 7.8e-7 to 1.4e-2) | open |
| [N-30](https://github.com/dwijn/RH/issues/226) | low | The direct form of the short-characteristics solution (I and Ψ, used everywhere) has a precision floor: J − S and the preconditioned I − Ψη are small differences of nearly equal numbers in optically thick, strongly scattering lines, so rounding u in I gives population errors ~u/ε. In double that is 1e-16/ε. It may be what stops H_6 near 1e-8 (N-4) | `rhf1d/bezier_1D.c`, `fillgamma.c:178-233` | Two-level atom, ALI, 10 points per decade, long-double reference: double direct form 2.4e-14 at ε = 1e-2, 1.3e-8 at ε = 1e-8 (∝ 1/ε); the same solver carrying D = I − S stays at 1e-13–1e-11 for all ε (`tests/mixprec/twolevel.c`, branch `exp/mixed-precision`). Not yet tested in RH: that needs Ieff built from D and 1 − Ψχ instead of from I | open |
| [N-31](https://github.com/dwijn/RH/issues/262) | medium | The rate equations of each atom and molecule are solved by LU decomposition (`SolveLinearEq`: partial pivoting with row scaling and one refinement step in the same precision), which is normwise stable but not componentwise accurate. Where levels are tightly coupled by collisions, as in the fine-structure multiplets of `He_large` (split by 1e-4 cm⁻¹), each diagonal element −Σ_j Γ_ji is dominated by the large rates within the multiplet, and its rounding erases the small rates out of it: the populations of these levels are wrong by up to 4.4e-3. The componentwise (Skeel) condition number is 1e13–1e14 in the deep layers and stays there after row and column equilibration, so no scaling helps; long-double LU lowers the error only by the mantissa ratio (2048). The rate matrix is a Markov generator, so the GTH variant of Gaussian elimination (Grassmann, Taksar & Heyman 1985), which forms every pivot as a sum of off-diagonal rates and never subtracts, is componentwise accurate at the same cost; it needs non-negative off-diagonal rates | `statequil.c:91,159`, `ludcmp.c` | `1d_falc_helarge` at the N-28 stall state, Γ dumped and solved in 60-digit arithmetic: LU relative error up to 4.4e-3 (levels 26–28 and 36–42: 4F ³F°, 5D ³D, 5F ³F°), GTH ≤ 2e-15. Noise of G (ECnoise δ, SOLVER.md S-2) 1.7e-3 with LU, 5e-14 with GTH; Ca II 1.3e-10 → 5.6e-13, Mg II 5e-11 → 2e-13. With GTH the reference recipe (no Ng, `ITER_LIMIT` 1e-10) converges in 667 iterations where LU stalls at 2.3e-3, and the Ng configuration converges in 63 where LU stops with SIGFPE at iteration 13 (N-28). Answers SOLVER.md §9 Q12. Same code upstream | fixed |
| [N-32](https://github.com/dwijn/RH/issues/263) | medium | Tabulated atomic data are interpolated in ways that can give non-physical values. Collision rates (`OMEGA`, `CE`, `CI`, `CP`, `CH0`, `CH+`, `CH`) are interpolated in temperature with a natural cubic spline in linear T (linear for two points, clamped outside the table), and explicit bound-free cross-sections with a natural cubic spline in linear wavelength. A spline is neither local nor positivity-preserving: a non-monotone table, a single bad value or a resonance makes it overshoot, ring through the table and go negative, giving negative collisional rates (B-42, B-43, B-44), negative photoionization rates and negative bound-free opacity contributions. Collision rates also vary roughly as power laws in T, which a spline in linear T represents poorly even on clean tables. Interpolating log C in log T is positive by construction, and monotone cubic Hermite interpolation (PCHIP) in log-log adds no new extrema (chosen 2026-11-08). The scope (widened 2026-11-09) is every interpolation of a tabulated quantity that must stay positive or within physical bounds; before the fix the code is reviewed for other such places | `collision.c:705-733`, `sortlambda.c:347` (active atoms), `metal.c:96,171` (passive atoms) | Collision tables: leave-one-out test on the 6272 distinct non-constant tables of the 88 atom files (each interior node predicted from the others, 40357 nodes): median / 90% / 99% relative error 2.2e-2 / 0.16 / 0.47 for the current spline, with 15 predictions ≤ 0; 1.5e-4 / 3e-3 / 4e-2 for PCHIP in log-log, none ≤ 0 (cubic spline in log-log 1.5e-4 / 2e-3 / 3e-2, linear in log-log 1.4e-4 / 1.4e-2 / 0.13). Compared with the current interpolant through all nodes, 10% of the tables change by more than 5–8% somewhere in their range. Bound-free cross-sections: the spline of 225 of the 936 explicit tables in 28 atom files (99 distinct) goes negative inside the table, including those of the passive Al, Fe, He, Mg, Na, O, S and Si of the default run, where 215 evaluations (wavelength and continuum) are negative, 109 of them Mg I 15 → 17; the summed metal bound-free opacity stays positive there. In `CI-VI_Judge` the photoionization rate 49 → 60 is negative at k = 37–42, which sends the rate equations there to LU (N-31) and keeps the iteration from converging; with the cross-sections clipped at zero (scratch build) it converges to 1e-3 in 71 iterations without Ng (N-28). Same code upstream | fixed |
| [N-33](https://github.com/dwijn/RH/issues/266) | medium | The 2-D stress case `2d_pert_h6_ca_mg` (perturbed FAL-C, 32 columns, H_6 + Ca II + Mg II, 8 threads) stops with SIGFPE at iteration 6 since B-40 (#257); before, it converged in 82 iterations. Its iteration 5 is violent either way: the largest change of H is 83 before #257, after which the iteration recovers, and 424 after it, after which Ca II and Mg II change by 69 and 171 and an `exp` overflow in `Bezier3_coeffs` (a negative opacity) stops the run. Iterations 1–4 are nearly the same in both (H 356 / 344, 1.26 / 1.29, 2.35 / 2.27, 1.03 / 1.01). B-40 lowered the spurious H⁻ density where hydrogen is ionized; the iteration was already close to the edge, so the problem is its lack of robustness on this atmosphere, not B-40 | iteration (`iterate.c`); crash in `rhsc2d/bezier_2D.c:538` | Bisect over the merges since S-1, 2026-11-08: 8a921e3, cc0bbbd, 0ba542d, a4e1365 and 56c9aec converge in 82 iterations, 6d1eeff (#257) stops at iteration 6; the same with GTH (N-31). Upstream does not have B-40 | fixed (#312) |
| [N-34](https://github.com/dwijn/RH/issues/268) | low | Collision rates given by formulas instead of tables (`SHULL82`, `BADNELL`, `AR85-CDI`, `AR85-CEA`, `AR85-CHP`, `AR85-CHH`, `BURGESS`) are not checked for sign or for the temperature range of their fits, unlike the tables, which N-32 keeps positive. Several can go negative for some coefficient sets: the dielectronic term of `SHULL82` has the factor 1 + B exp(−T1/T), negative for B < −1; `AR85-CHH` has 1 + c exp(d T/10⁴), negative for c < −exp(−d T/10⁴); `BADNELL` sums the given c_i exp(−E_i/T) without a check; `AR85-CDI` sets a negative sum to zero without a message. The fits are used at any temperature, except `AR85-CHP` and `AR85-CHH`, which drop to zero outside [T1, T2]. `AR85-CHP` takes the proton density from `atmos.H->n[5]`, the continuum only of a six-level hydrogen model: with `H_20` it uses a bound level, with fewer than six levels it reads past the array. Not yet checked numerically | `collision.c:563-680` (reading), `collision.c:800-980` (evaluation) | Code reading, 2026-11-09, after the N-32 audit of tabulated data. Of the shipped atom files only `C_I+II_9.atom` uses formula records (`SHULL82`, `BADNELL`, `AR85-CDI` and `AR85-CEA` four times each, and `SUMMERS`); none uses `AR85-CHP`, `AR85-CHH` or `BURGESS`. Seen in passing, to check against the sources: in that file the records between C I and the two C II ²P° levels (`AR85-CDI`, `AR85-CEA`) and between C III and those levels (`SHULL82`, `BADNELL`) split the total 0.6 : 0.4 between ²P°1/2 (g = 2) and ²P°3/2 (g = 4), the reverse of the ratio of the statistical weights. Same code upstream | open |
| [N-35](https://github.com/dwijn/RH/issues/277) | low | The H⁻ bound-free cross-section (`Hminus_bf`) is the table of Geltman (1962), 34 nodes from 0 to the threshold at 1641.9 nm, interpolated with a natural cubic spline in wavelength. The photodetachment calculation of Wishart (1979) is accurate to about 1%, and John (1988, A&A 193, 189, Table 2) fits it to 1% between 0.125 µm and the threshold with σ = 10⁻¹⁸ λ³ (1/λ − 1/λ₀)^{3/2} Σₙ Cₙ (1/λ − 1/λ₀)^{(n−1)/2} cm² (λ in µm, λ₀ = 1.6419 µm), which has the threshold law built in. Geltman / John: 0.81 at 125 nm, 0.88 at 150, 0.94 at 300, 0.97 at 500, 0.98–0.99 from 600 to 1400, 0.96 at 1500, then 1.31 at 1600, 1.78 at 1620 and 3.1 at 1635 nm, where the table does not follow the (Δ)^{3/2} threshold law. So H⁻ bound-free, the main visible and near-infrared continuum opacity of the photosphere, is 1–3% low in the visible, 6% low at 300 nm and up to 19% low at 125 nm. Near the threshold it is too large, but there free-free dominates. RH already uses John (1988, Table 3a) for H⁻ free-free beyond 9113 nm (`Hminus_ff_long`). Replacing the table with the fit needs no interpolation; below 125 nm the fit is not valid and the table can stay. Changes most photospheric spectra slightly; effect not measured. Interpolation of the table stays as it is (2026-11-09: positive, and as accurate as PCHIP in a leave-one-out test; log-log is impossible with σ = 0 at both ends) | `hydrogen.c:279-336` | Table against John's fit (coefficients checked against the paper's scan on ADS), 2026-11-09. Same code upstream | open |
| [N-36](https://github.com/dwijn/RH/issues/285) | medium | With `SOLVE_NE = ITERATION` the default solver (ne updated outside G by `Background()` and `Solve_ne`) and the solvers on the map G with ne in the state (S-2b: `BackgroundOfState`, `Solve_ne_OfState`; `FIXED_POINT`, `ANDERSON`) converge to different solutions. On `1d_falc_solvene` (FAL-C, H_6 + Ca II, PRD, background in double), both run to `ITER_LIMIT` 10⁻¹⁰ with no acceleration, the emergent spectra differ by 1.7×10⁻⁴ (relative to the peak) at 163.2 nm, and by 4.7×10⁻⁶ in the lines. The populations differ by 4.8×10⁻⁶. Two converged solutions of the same equations should agree to about 10⁻¹⁰, so one formulation differs from the other: in what the background is computed from, the chemical equilibrium, or the final background. Which one is right is not known. A continuum at 163.2 nm points to a metal bound-free edge whose opacity depends on ne; first compare ne and the background opacity there between the two runs. For S-5 the ne cases are judged against references made with ne in the state (SOLVER.md, S-5 calibration) | `solver.c` (G with ne in the state), `background.c` (`BackgroundOfState`), `solvene.c`, `iterate.c` | References `claude-tmp/s1/refs2/out-ne` (default solver) and `out-ne-fp` (`FIXED_POINT`), both at f5e9847, 466 iterations each, F_SE ≤ 1e-10; found in the S-5 `ITER_LIMIT` calibration, 2026-11-09. The S-2b formulation exists only here | open |
| [N-37](https://github.com/dwijn/RH/issues/320) | medium | On `2d_pert_h6_ca_mg` the 2-D iteration does not converge: every solver, plain G included, stalls on a floor in the maximum relative change of 2–5×10⁻⁵ (PRD) or 5×10⁻⁴–10⁻³ (all CRD), so `ITER_LIMIT` 10⁻⁶ is never reached and 2-D comparisons are judged on the first iteration below 10⁻³ and 10⁻⁴. The floor is a cycle, not rounding. G is deterministic: a rerun matches line by line, and `BACKGROUND_DOUBLE` changes nothing. In the all-CRD run H n = 5 cycles with period 5 (max change 1.04, 0.54, 0.63, 0.60, 1.02×10⁻³) at one depth (z = 31, about 9 kK) in columns 15–20; Ca II and Mg II settle to 10⁻⁶. A restart from the end state returns to the same cycle. With the populations fixed, the J iteration of the final scattering sweeps alternates with period 2 (2.33/2.34×10⁻³ CRD, 1.53×10⁻⁴ PRD). The cause is the monotonic option of `Quadr`, the quadratic interpolation of upwind and downwind values between columns (Auer & Paletou 1994): where the parabola overshoots the three nodes it is replaced by the linear interpolation of the nearest pair, which makes I_uw, S_uw and χ_uw jump as functions of the values, so G is discontinuous and need not have a fixed point. Clamping the parabola to [v_min, v_max] instead (continuous) removes it: the J sweeps converge geometrically (5.6×10⁻¹⁰ after 50) and G converges at a ratio of 0.93 per iteration through the old floor. Replacing the Steffen limiter of the Bézier derivative by the harmonic mean keeps the cycle. Removing either limiter overflows in the first sweep. The clamped interpolation moves the populations (median 3×10⁻⁴ for H) and by up to 24 % (H II) and 74 % (Mg III) in the top layers, so a fix changes results. Not tested: 3-D, and whether the clamp is accurate enough. Related: N-3 (convergence criterion dominated by trace levels) | `rhsc2d/piecewise_2D.c` (`Quadr`), used by `bezier_2D.c`, `piecewise_2D.c`, `piecestokes_2D.c`, `solvelong.c` | `claude-tmp/s1/floor2d` (README): rerun, bgdbl, linear, crd, cycle (restarts mapping the cycle), smooth/smooth2 (limiter switches in a scratch build), upstream (c95bdf6, same inputs, all CRD: H and Ca II cycle with period 2, J sweeps alternate at 2.6×10⁻³). Also seen on `2d_bf_s12_prd` (about 3×10⁻⁵), 2026-11-12 | open |
| [N-38](https://github.com/dwijn/RH/issues/326) | medium | With PRD, the 2-D iteration on `2d_pert_h6_ca_mg` ends in a second cycle, independent of N-37: with the interpolation between columns made continuous (monotonic cubic Bézier or linear, branch `fix/N-37-bezier-interp`), plain G (`FIXED_POINT`, `PRD_N_MAX_ITER` 3, `PRD_ITER_LIMIT` 0) converges through the N-37 floor until Mg II settles into a period-2 cycle: the maximum change of Mg II alternates between 4.2555 and 4.2595×10⁻⁶ from iteration ~355 on with Bézier (H still falling geometrically, 2.3×10⁻⁵ at 400), and between 1.3096 and 1.3106×10⁻⁴, slowly growing, from about iteration 200 with linear interpolation. The PRD sub-iteration deltas of the Mg II lines repeat exactly every two outer iterations (Bézier: 3.74, 5.52, 6.03 then 4.55, 4.39, 4.08×10⁻⁶ for h). Without PRD there is no cycle: from the all-CRD state both interpolations converge (Mg II 8.6×10⁻⁷ and 2.0×10⁻⁷ after 200 iterations, still falling). Not yet known: whether G has an unstable direction there (an eigenvalue at or below −1 of the linearized G with its PRD sub-iterations, a flip, which damping or Anderson would remove) or whether something in the PRD path is discontinuous; which Mg II level and depths cycle; whether upstream has it (upstream with PRD stalls at 1.2–2.8×10⁻⁵, where N-37 dominates). Related: N-37, N-3 | PRD path of the iteration (`iterate.c`, PRD sub-iterations); not localized | `claude-tmp/s1/n37` (README): full/{bez,lin} (400 iterations from scratch, 8 threads, 2026-11-12), restart/{bez,lin} (from the all-CRD state, 200 iterations) | open |
| [N-39](https://github.com/dwijn/RH/issues/333) | medium | The 3-D stress case `3d_pert_h6_ca` (FAL-C perturbed over 8 × 8 columns by `tests/make_perturbed_atmos.py`, H_6 and Ca II with PRD, 8 threads) stops with SIGFPE (floating-point overflow) in the first iteration, before any solver acts. The lambda iteration of the background scattering at the start (`initScatter`, `N_MAX_SCATTER` 50) does not converge: it ends in a period-2 cycle whose largest relative change of J grows to 8.3×10⁵ and 3.8×10⁴ (Ca II K and H cores, λ# 293 and 338) with `INTERPOLATE_3D = BICUBIC_3D`, and to 1.5 and 1.2×10³ with `LINEAR_3D`. The first iteration then changes the populations by 10⁴⁸ (H) and 10⁴⁰ (Ca II), and the run dies. 4951ebb (before #328) and 41738a6 give bit-identical J changes, so it is not the interpolation change of N-37. With 1 thread the J changes equal those with 8, in both trees, so it is not a thread effect. Upstream c95bdf6 on the same input also ends `initScatter` in a period-2 cycle, but a bounded one: relative changes of 50–490 at the Lyα and Hα cores (λ# 129, 144, 492, 493) with either interpolation. Its first iteration changes the populations by 5×10² (H) and 3.7×10³ (Ca II), and the run completes. With the populations fixed, background scattering hardly matters in a strong line core, so J there should barely change between sweeps. Changes of 10²–10⁵ point to the 3-D formal solution or to points where J is near zero or negative. Not yet known: which points cycle; whether J is negative there; whether the generated atmosphere is reasonable (dx 100 km, temperature ±35 %, height shifts up to 250 km above the temperature minimum); why modernize grows where upstream stays bounded. Related: N-37, N-38 **Cause** (2026-11-13): the parabolic integration of S along the long characteristics (`SolveLong`). The weight of the downwind S in a parabolic step is (w₂ − w₁Δτ_uw)/(Δτ_dw(Δτ_uw + Δτ_dw)): negative once the upwind section is optically thicker than 2, and about −1/Δτ_dw for a thin downwind one. Since B-24 (#152) the downwind point of the last section is the local point, so J depends negatively on the local S and the pure lambda iteration of `initScatter` alternates and grows. `git bisect` from upstream to modernize: the first bad commit is the merge of #152; its parent gives upstream's numbers to every digit. With the last section piecewise linear the cycle is bounded as upstream's; with all sections linear the iteration converges (largest change 0.20, 0.13, 0.089, … 0.033 after 9 sweeps, ratio 0.72), and the full case converges (Ng, 173 iterations to 10⁻⁶). So the bounded cycle upstream has the same cause. Not float storage: `BACKGROUND_DOUBLE` and a build with double line profiles give the same J changes to every printed digit; not N-24: reverting it still diverges | `rhsc3d/solvelong.c`, `rhsc2d/solvelong.c` (`SolveLong`) | `claude-tmp/s1/n38/p3chk` (README): runs/{base,sweep} (4951ebb, 41738a6), up-bicubic, up-linear, sweep-linear, sweep1t and up1t (1 thread); `N_MAX_ITER` 1, 2026-11-13; `claude-tmp/s1/n38/p3chk` run5–9, bisect/, lcA/, lcB/ | open |
| [N-40](https://github.com/dwijn/RH/issues/338) | medium | A negative total opacity makes the formal solvers overflow and the run die with SIGFPE (or, without traps, carry inf and NaN into J and the populations) instead of stopping with a message: the intensity grows as exp(−τ) across a run of negative optical depth, and `exp` overflows below −709. Total opacity does go negative (strong inverted lines, SOLVER.md trap 5; `1d_falc_h20` under 10 % perturbations, S-2), and no check catches it: S-13 showed the line test of the guard does not prevent it (S-14). Plan: a detector in `Formal()` (every geometry) that computes, per wavelength and ray, the most negative slant optical depth summed over consecutive negative intervals and flags anything below −50; the formal solution of that wavelength is skipped, Anderson rejects the trial point, and every other case stops with a message naming the wavelength, depth and optical depth. Decided by the user 2026-11-13 | `rhf1d/formal.c`, `rhsc2d/formal.c`, `rhsc3d/formal.c`, `rhsphere/formal.c`; `solver.c` | S-13 option C measurements (`claude-tmp/s1/s13/chi`, README): worst run −0.094 (`falc_mgi66`), −4.9×10⁻⁶ (`falc_helarge`) | open |

### 3.3 Memory

Estimated for a realistic 3-D run: 256×256×200 grid (one array V = 105 MB),
24 rays, 5 atoms with ΣNlevel = 190, ~1000 transitions, Nspect = 4000.

| ID | Sev | Item | Formula | 3-D size | Where | Status |
|---|---|---|---|---|---|---|
| [M-1](https://github.com/dwijn/RH/issues/29) | high | Line profiles φ stored per ray (moving atmosphere). Measured 2026-11-06: 14.0 of the 16.3 GB peak RSS (86 %) of the 2-D 128×82 perf case (H_6, Ca II, Mg II, Fe33). Stored as float they halve, peak −43 %; populations change by ≤ 3e-6 (mixed-precision study, §2.4). Float since #230: peak 17.5 → 10.7 GB (−39 %); still stored per ray | 2·Nrays·ΣNλ·V | ~50 TB (1 TB static) | `profile.c:340-357` | open |
| [M-2](https://github.com/dwijn/RH/issues/30) | high | Γ and collision matrix C stored dense per point | 2·ΣNlevel²·V | 1.7 TB | `initial_xdr.c:567`, `readatom.c:569` | open |
| [M-3](https://github.com/dwijn/RH/issues/31) | high | Per-thread copies of Γ and Rij/Rji (WIP lock fix) | Nthreads·(ΣNlevel² + 2Ntrans)·V | ~1 TB per thread | `opacity.c:571`, `readatom.c:323,423` | open |
| [M-4](https://github.com/dwijn/RH/issues/32) | high | PRD redistribution weight files per point (also a disk problem: they filled /tmp). With the interpolation folded into the weights (P-8) the weights have fewer entries (Ca II H, K 5×, Lyα 1.2×, Mg II h, k 1.05×) and, as float, can live in memory: 2–10× smaller than the double files. Since #232 the weights have the interpolation folded in and are float: about 10× smaller (Ca II angle-averaged 2.5 → 0.25 MB per line, angle-dependent 253 → 24.5 MB per 1-D column, 5.5 → 0.56 GB for a 4×82 2-D cut); still a scratch file | ~149 KB/point per line | ~2 TB per PRD line; 1.8 GB per run on a 63×64 grid | `scatter.c` | open |
| [M-5](https://github.com/dwijn/RH/issues/33) | medium | J (and Jgas for hybrid PRD) | Nspect·V each | 420 GB each | `initial_xdr.c:88,103` | open |
| [M-6](https://github.com/dwijn/RH/issues/34) | medium | Hybrid-PRD tables; O(Nspect²) setup | — | > 10 TB | `initial_xdr.c:120-300` | open |
| [M-7](https://github.com/dwijn/RH/issues/35) | medium | Background opacity file. Float since #230 (emissivities scaled by 2^100), half the size; still a file re-read every iteration | ≥ 3·Nspect·V | ≥ 1.3 TB | `background.c` | open |
| [M-8](https://github.com/dwijn/RH/issues/36) | low | Per-thread work arrays; 6 full-grid malloc/free per Formal call | ~(3·Nrt + 3·ΣNlevel + 11)·V per thread | ~65 GB/thread | `rhsc2d/formal.c:105-123` | open |
| [M-9](https://github.com/dwijn/RH/issues/37) | low | WIP tree uses 5–30 % more peak memory than upstream in 1-D | — | measured | — | open |

### 3.4 Speed and threading

| ID | Sev | Issue | Where | Status |
|---|---|---|---|---|
| [P-1](https://github.com/dwijn/RH/issues/38) | high | Angle-independent opacity (continua, angle-averaged PRD) and the full χ/η sum recomputed for every ray and direction in moving atmospheres | `rhsc2d/formal.c:162-166`, `opacity.c:318-325` | open |
| [P-2](https://github.com/dwijn/RH/issues/39) | high | Statistical equilibrium serial over grid points with 5+ mallocs per point; O(Nspace·Nlevel³) (~1e3 s/iteration for 13 M points and 50 levels) | `statequil.c:58-96`, `ludcmp.c:51-58,101` | open |
| [P-3](https://github.com/dwijn/RH/issues/40) | high | 3-D interpolation does 8 integer `MODULO`s per call, 5 calls per point | `rhsc3d/geometry.c:618-642` | open |
| [P-4](https://github.com/dwijn/RH/issues/41) | medium | Static strided wavelength split; thread slot hard-wired to `nspect % Nthreads`, so no dynamic scheduling; load imbalance | `iterate.c:237-263`, `opacity.c:94,483`, `fillgamma.c:199,348,547` | open |
| [P-5](https://github.com/dwijn/RH/issues/42) | medium | Serial Γ/Rij merge and zeroing of per-thread copies (WIP): O(Nthreads·(ΣNlevel² + 2Ntrans)·Nspace) per iteration | `fillgamma.c:101-176,419-518` | open |
| [P-6](https://github.com/dwijn/RH/issues/43) | medium | `addtoGamma` and `addtoRates` recompute the same Vij·wla·wmu | `fillgamma.c` | open |
| [P-7](https://github.com/dwijn/RH/issues/44) | medium | Build flags: `-O2` only; gcc 8 does not vectorize at -O2. `-O3 -march=native` gives −14 %. The per-transition loops (`addPreconditionedRates`, `addtoRates`, `Opacity`) do not vectorize even at -O3: their arrays may alias. With `restrict` and -O3 the `addPreconditionedRates` loop is 1.4× faster in double and 2.8× in float (microbenchmark, §2.4). Fixed in #238: `restrict`, `#pragma omp simd` and -O3 by default; whole runs 1.15–1.23× faster, byte-identical | `Makefile` | fixed |
| [P-8](https://github.com/dwijn/RH/issues/45) | high | PRD scattering integral: O(Nspect) linear scans in hybrid PRD; weights re-read from disk every PRD iteration. Measured 2026-11-06: in angle-averaged and hybrid PRD the integrals are a third of a 1-D run, mostly spent re-interpolating J onto the fine absorption grid on every call (`Linear` + `Hunt` 33 % of a 1-D H_6 + Ca II + Mg II run, 2e8 `Hunt` calls in 10 iterations). The integral is linear in J on the line's own grid, so the interpolation can be folded into the weights once; every later call is a short dot product. Prototype (`prdfold`, branch `exp/mixed-precision`): populations as before to ≤ 1.4e-9, same iteration counts, integrals 12× faster (H + Ca + Mg angle-averaged), 11× (Mg II hybrid), 35× (Bifrost column Ca II hybrid), 8× (2-D Ca II hybrid); whole 1-D runs 1.2–1.5× faster. Float integrals are safe in all three PRD modes (populations ≤ 2e-5) | `scatter.c:282-347`, `formal.c:275-336` | fixed |
| [P-9](https://github.com/dwijn/RH/issues/46) | low | `time.out` uses `clock()` (CPU time summed over threads), so it never shows a parallel speedup | `getcpu.c:61` | open |
| [P-10](https://github.com/dwijn/RH/issues/47) | low | Writing `opacity.out` is serial and large (26 s of 69 s in the Bifrost 2-D run) | `writeopac_xdr.c` | open |
| [P-11](https://github.com/dwijn/RH/issues/48) | — | No MPI / spatial domain decomposition: a 256²×200 run cannot fit on one node | — | open |
| [P-12](https://github.com/dwijn/RH/issues/194) | medium | Hybrid PRD: building the gas-frame interpolation table in `initSolution` scans the whole wavelength grid four times for each of the 2·Nrays·Nspect·Nspace entries, in two passes, on one thread: O(N_lamuk·Nspect). On a 3-D Bifrost cut of 42×42×89 (2.2e9 entries) it takes most of 1 h 38 min. The grid is sorted, so a bracketing search (or the previous wavelength's bracket) and a loop threaded over points would make it negligible. Found with B-32. In 1-D too: 22 % of a moving FAL-C Mg II hybrid-PRD run to convergence (gprof, 2026-11-06) | `initial_xdr.c:180-330` | open |
| [P-13](https://github.com/dwijn/RH/issues/218) | high | Background opacity is serial over wavelength and redoes wavelength-independent work for every Kurucz line at every wavelength (and, in moving or polarized atmospheres, for every ray and direction): the partition function is interpolated over all points (`Linear`), `T^0.3` (Unsöld) and the level populations (two `exp`) are evaluated per point, and `Metal_bf`/`Hydrogen_bf` rebuild spline coefficients for every continuum at every wavelength. With large line lists this stage dominates: 255–285 nm Fe/Mg UV synthesis (60k-point wavelength table), ~2 min with 98k lines and ~20 min with 810k lines on a laptop, 24+ min on sundog1, while the threaded iterations after it are fast. Stand-in (20k wavelengths at 390–400 nm, 5235 lines, FAL-C, 1 thread): background 8.5 s of 15.5 s; `Voigt` 49 %, `pow` 15 %, rest of `RLKProfile` 12 %, `exp` 8 %, partition-function and spline interpolation 8 %, `pwrite` 3 %. Plan: (1) remove the shared state of T-5 and compute per-point tables (ln pf per element and stage, `T^0.3`, thermal widths) once; (2) thread over wavelengths with dynamically scheduled contiguous chunks, each worker reserving its records in `background.dat` with an atomic add (record order then varies with N > 1; readers go through `brs.out`, values are unchanged; byte-identical at N = 1); (3) verify record by record, serial against threaded, static/moving/Stokes, plus ThreadSanitizer. *Step 1 done 2026-11-04 (branch `perf/background-threads`, not yet merged):* `backgroundPrepare`/`rlk_prepare` do all one-time setup before the loop; reentrant `SplineTable` (coefficients per passive continuum, built once); per-point tables of ln pf, `T^0.3` and Doppler widths; `passive_bb` cache passed in explicitly. Output byte-identical to `modernize` in every case (`run_tests.sh --exact` 50/50 with ASan; background records of static, moving, full-Stokes, UV 77k-line and 150–400 nm cases). Background stage at 1 thread: static 9.6 → 5.9 s, moving 88.6 → 57.2 s, UV 77k lines 25.1 → 17.0 s, full Stokes 500 → 463 s (Stokes is dominated by `Voigt` over the Zeeman components). *Step 2 done 2026-11-06 (same branch):* the wavelength loop runs through `parallelForDynamic` (atomic chunk counter, contiguous chunks) with per-worker scratch arrays and `passive_bb` caches; each wavelength reserves its records with an atomic add, and `backgroundRecords()` in `readj.c` is the single definition of a wavelength's record count. Verified on sundog3/sundog2 (Linux gcc): N = 1 byte-identical to step 1 including full Stokes; at N = 3, 8, 16 and 32 every background record matches the serial run value for value (order-independent comparison through `brs.out`) for static, moving, full-Stokes, UV 77k-line and 150–400 nm cases; ThreadSanitizer clean. Background loop on the 810k-line, 60k-wavelength UV case (sundog3): 2538 s at 1 thread, 681 s at 4, 354 s at 8, 187 s at 16, 113 s at 32 (22× over 1 thread, 28× over the old serial 3158 s). Whole runs at 16 threads against the old serial code: static 53 → 6 s, moving 304 → 20.5 s, full Stokes 1490 → 114 s, UV 77k lines 273 → 41 s (the rest is mostly the threaded formal solutions; reading the 77k lines takes 0.25 s, 810k lines 1.3 s, on an M2 Pro) | `background.c:327-625`, `kurucz.c` (`rlk_opacity`, `RLKProfile`), `metal.c:126`, `hydrogen.c:310` | fixed |
| [P-14](https://github.com/dwijn/RH/issues/227) | high | Angle-dependent PRD (`PRDAngleScatter`) is serial over locations, re-interpolates the intensity of every absorption direction on every call, and re-reads all its weights from disk on every call. Moving FAL-C, Ca II, 10 iterations: 270 s in the scattering integrals against 1 s for everything else; a third of it computing R_II once, the rest per-call interpolation and loop overhead; the weights are 253 MB per line for one 82-point column, re-read 80 times per line. Thread over locations (each uses only its own data, as in `PRDScatter`), fold the interpolation into the weights (P-8), store them in float. In 2-D and 3-D the weights grow as Nspace·Nrays², so they must be streamed or recomputed | `scatter.c:553-834` | fixed |
| [P-15](https://github.com/dwijn/RH/issues/239) | medium | Background opacity of a static atmosphere with a magnetic field (`FULL_STOKES`, and `FIELD_FREE`, whose final formal solution is full Stokes) recomputes every Kurucz line profile 2 × Nrays times per wavelength: `Background` calls `rlk_opacity` once per ray and direction, because the polarized records differ between them, and each call evaluates `Voigt` with the Faraday function for every Zeeman component at every depth. Without velocities H and F depend only on the damping and the frequency offset; only the geometric factors (cos γ, cos 2χ, sin 2χ and the sign for the direction) change from ray to ray. FAL-C, Ca II, `gf0400.10` (5235 lines), 2291 wavelengths over 390–400 nm, 5 rays, 1 thread, M2 Pro: background 46.1 s with the 500 G field of `FALC_82_500G_45.B`, 47.6 s `FIELD_FREE`, 0.58 s without a field; 27.4 s with 3 rays instead of 5 (linear in Nrays); `Voigt` takes 87 % of the samples and `RLKProfile` 7 %. Compute the σ−, π, σ+ sums of H and F once per line and depth at each wavelength and apply the geometry per ray: 2 × Nrays times fewer `Voigt` calls, about 4–5× on the Stokes background. Moving atmospheres do not gain (the frequency offset depends on the ray) | `background.c:419-534`, `kurucz.c` (`rlk_opacity`, `RLKProfile:798-839`) | open |

### 3.5 Thread safety

| ID | Sev | Issue | Where | Status |
|---|---|---|---|---|
| [T-1](https://github.com/dwijn/RH/issues/49) | medium | Shared global `messageStr` written from threaded code paths | `rhf1d/formal.c:209,228`, `bezier_1D.c`, `piecestokes_1D.c:93`, `opacity.c:312,408`, `readj.c` | open |
| [T-2](https://github.com/dwijn/RH/issues/50) | medium | Function-static lazy initialization (data race): Voigt tables, `BiLinear` hunt indices, wigner, hydrogen, metal | `voigt.c:261,447`, `linear.c:59`, `wigner.c:166`, `hydrogen.c:336,548,693`, `metal.c:160` | open |
| [T-3](https://github.com/dwijn/RH/issues/51) | low | File-static spline state, not reentrant (only called serially today) | `spline.c:22-36`, `expspline.c:25-36` | open |
| [T-4](https://github.com/dwijn/RH/issues/179) | high | Hybrid PRD (`PRD_ANGLE_APPROX`) with `N_THREADS > 1`: the formal solution runs threads over wavelengths, and each adds its intensity into the gas-frame mean intensity `Jgas` of *neighbouring* wavelengths (`Jgas[iprdh[idx]][k] +=`, and the `prdh_limit_mem` branch) without a lock, so concurrent updates are lost. After one iteration on a Bifrost column, `Jgas` on 4 threads differs from 1 thread by up to 34% (30% of entries), two 4-thread runs by 51%. The corrupted J feeds the PRD emission profile: the same run converges in 30 iterations on 1 thread (1-D; 36 in 2-D, byte-identical when repeated) and in 42 to 78 on 4, or stalls for 200 (2-D); "converged" populations differ by 2.5%. Present in all three geometries and upstream | `rhf1d/formal.c:295-330`, `rhsc2d/formal.c:283-339`, `rhsc3d/formal.c:256-312` | fixed |
| [T-5](https://github.com/dwijn/RH/issues/219) | medium | Shared state in the background opacity loop, which blocks threading it (P-13): `rlk_opacity` allocates `element->n` and calls `LTEpops_elem` lazily and calls `RLKZeeman` on first use; `MolecularOpacity` calls `MolZeeman` on first use; `passive_bb` keeps a static list of lines (damping cache) across wavelengths; `Metal_bf` and `Hydrogen_bf` use the file-static spline (T-3); the free-free routines in `hydrogen.c` keep static `Hunt` indices and lazily allocated tables (T-2); records in `background.dat` are numbered by a running counter. All removed on `perf/background-threads`: one-time setup in `backgroundPrepare`/`rlk_prepare`, reentrant `SplineTable`, explicit `passive_bb` cache and free-free state (step 1); atomic record reservation (step 2) | `kurucz.c`, `opacity.c` (`MolecularOpacity`), `metal.c:156-158`, `hydrogen.c:336,547,690`, `background.c:307,569,622` | fixed |

### 3.6 Code quality and build

| ID | Sev | Issue | Status |
|---|---|---|---|
| [Q-1](https://github.com/dwijn/RH/issues/52) | high | Hand-written makefiles with archive-member rules: not safe under `make -j`; hand-maintained, partly wrong dependencies; four ~90 % identical geometry makefiles | fixed |
| [Q-2](https://github.com/dwijn/RH/issues/53) | high | Fortran still required by rhsc2d/rhsc3d/rhsphere/tools makefiles (`librh_f90.a`) after the WIP dropped it at top level | fixed |
| [Q-3](https://github.com/dwijn/RH/issues/54) | medium | 354 warnings with `-Wall -Wextra -std=c11 -pedantic` (138 maybe-uninitialized, 140 unused variables, 32 set-but-unused, 30 unused parameters, 8 sign-compare, 5 implicit declarations — pread/pwrite/fileno need `_POSIX_C_SOURCE`) | fixed |
| [Q-4](https://github.com/dwijn/RH/issues/55) | medium | Pervasive global state: 315 `extern` declarations in .c files; `atmos` referenced 2145 times in 108 files | open |
| [Q-5](https://github.com/dwijn/RH/issues/56) | medium | ~25 % of code in duplicated blocks; rhsc2d 53 %, rhsc3d 52 %, rhf1d 42 %, rhsphere 38 % | open |
| [Q-6](https://github.com/dwijn/RH/issues/57) | medium | Very long functions: 77 over 100 lines, 31 over 200 (`initSolution` 647, `SortLambda` 572, `Background` 558, `readAtom` 543, `CollisionRate` 494) | open |
| [Q-7](https://github.com/dwijn/RH/issues/58) | high | Error handling: `Error()` calls `exit()` (262 sites), so the library cannot be unit-tested; 395 `sprintf`, 0 `snprintf`; unbounded `%s` in `sscanf` | open |
| [Q-8](https://github.com/dwijn/RH/issues/59) | medium | Memory ownership implicit; 493 allocations vs 310 frees; `matrix_double` has no NULL check | open |
| [Q-9](https://github.com/dwijn/RH/issues/60) | medium | Headers not self-contained; reserved include-guard names; `rh.h` has none and pulls in `<rpc/types.h>` and `<pthread.h>` everywhere | open |
| [Q-10](https://github.com/dwijn/RH/issues/61) | low | Dead or bit-rotted code: `collision_15d.c` (does not compile), `collision_Oslo.c`, `rhsphere/readinput.c`, `rhsphere/solveplane.c`, `tools/constructatom.c` (does not compile upstream either), `tools/adjustbf.c` | fixed |
| [Q-11](https://github.com/dwijn/RH/issues/62) | low | Five Voigt algorithms plus Fortran variants | fixed |
| [Q-12](https://github.com/dwijn/RH/issues/63) | low | `tests/make_tiled_2d_atmos.py` uses the wrong array order | fixed |
| [Q-13](https://github.com/dwijn/RH/issues/64) | low | Python readers depend on `xdrlib` (removed in Python 3.13) | fixed |
| [Q-14](https://github.com/dwijn/RH/issues/65) | — | No unit tests (Voigt, LU, splines, w3j, formal solvers on analytic slabs), no CI | fixed |
| [Q-15](https://github.com/dwijn/RH/issues/83) | low | Default `CFLAGS` lack `-g`, so N-11 backtraces resolve to function names only, without file and line | fixed |
| [Q-16](https://github.com/dwijn/RH/issues/84) | medium | `Error(ERROR_LEVEL_1, …)` prints "Trying to continue" but terminates (exit status 2): the continue branch in `error.c` is unreachable. 74 call sites are written to recover afterwards (return FALSE, fall back, skip an output file), so their recovery code is dead, and `ERROR_LEVEL_1` and `_2` differ only in exit status (2 vs 3). Decide per site between a warning with recovery and a fatal error; define the exit statuses together with N-1 | open |
| [Q-17](https://github.com/dwijn/RH/issues/96) | medium | ~30 unchecked I/O results: `fread`/`fwrite`, `xdr_*`, `fscanf`/`fgets` return values assigned and never tested (e.g. a truncated Barklem table goes unnoticed) | fixed |
| [Q-18](https://github.com/dwijn/RH/issues/105) | medium | Multi-D test cases use perturbed FAL-C atmospheres; realistic cases need cut-outs of freely available 3-D RMHD snapshots (e.g. Bifrost, MURaM), small enough for the stress suite and CI, with converter scripts in `tests/` | fixed |
| [Q-19](https://github.com/dwijn/RH/issues/157) | medium | Manual out of date and without source: only `docs/rhmanual.pdf` (26 pages, November 2010, pdfTeX) exists, its LaTeX source is not in the repository or upstream. It describes the makefile build and Fortran, keywords that have since been removed (`S_PARABOLIC`, the Voigt choices), and nothing of CMake, threads, the tests or the solver changes. To be revised once the code settles, with its source in the repository; `pdftotext` extracts the text cleanly as a starting point | open |
| [Q-20](https://github.com/dwijn/RH/issues/158) | medium | Unit tests documented only in source comments: write up each problem, its analytic solution and the expected accuracy and order of convergence, with derivations, in a LaTeX document (`docs/unit-tests/`) | fixed |
| [Q-21](https://github.com/dwijn/RH/issues/269) | low | `tests/analyze_run.py compare` (the S-1 spectrum errors) and the tree-against-tree spectrum comparison of `run_stress.sh` read `spectrum.out` as unpolarized. A run that writes Q, U and V (a magnetic field was read, or `BACKGROUND_POLARIZATION = TRUE`) still passes the size check, since the four arrays together are also a multiple of Nspect, so I, the vacuum-to-air record and Q, U, V are read as one long I: the reported spectrum errors are garbage (values up to ±1.4e308 and NaN in the C I+II runs of B-42, FAL-C with a 500 G field). `run_tests.sh` passes `--stokes` when `STOKES_INPUT` is set and misses `BACKGROUND_POLARIZATION` | fixed |
| [Q-22](https://github.com/dwijn/RH/issues/272) | medium | `tests/analyze_run.py compare` groups the emergent intensities by the wrong axis. RH stores I with the wavelength as the slowest index (`I[Nspect][rays…]`, `initial_xdr.c`), as `compare_spectra.py` reads it, but `compare` reshapes it to (rays, Nspect), so each of its "wavelengths" mixes Nrays intensities from wavelengths spread over the whole spectrum. Its `spec_wl` normalizes by the peak of such a group, and `spec_line` looks up a line's window in these groups, so both are wrong for every run, and `worst_line` names the wrong line. Affected: the spectrum columns of the S-1 table in SOLVER.md and the spectrum effects measured with it since S-1 | fixed |

Output formats are a contract: the Python (`python/rhanalyze`) and IDL readers
unpack the XDR files field by field, and `rhopacity.py` seeks into raw
native-endian files by record. Layout changes must update or version them.


### 3.7 Wishlist

New capabilities, not problems with the existing code.

| ID | Sev | Item | Status |
|---|---|---|---|
| [F-1](https://github.com/dwijn/RH/issues/109) | medium | Tool for building RH model atoms from standard data sources: choosing the right model atom is one of the hard parts of running RH well. Assemble levels, lines, photoionization cross-sections and collisional data from standard sources (NIST, Kurucz, TOPbase/Opacity Project, CHIANTI, collision databases), with options to reduce a large model (level/line selection, superlevels) and checks of completeness. Replaces the removed `constructatom` (Kurucz lines only, levels and lines, Si only) | open |
| [F-2](https://github.com/dwijn/RH/issues/164) | low | Forbidden (non-E1) lines with polarization: RH treats every line as electric dipole (E1). Magnetic-dipole (M1) lines have the same rank-1 Zeeman pattern and strengths, but their linear polarization is rotated by 90° (eta_Q, eta_U, rho_Q, rho_U change sign; eta_V, rho_V do not). Electric-quadrupole (E2) lines need rank-2 strengths and |dM| = 2 components; mixed M1+E2 lines (e.g. [O I] 630.03 nm, blended with Ni I 630.04 nm) need both plus their interference. Scope undecided: a per-line multipole flag with the M1 sign change is small; E2 and mixtures are a larger change to `zeeman.c` and `profile.c`, and RH may not be the right tool for such lines at all | open |
| [F-3](https://github.com/dwijn/RH/issues/206) | low | Incomplete Paschen-Back effect: RH has only the linear Zeeman effect (`zeeman.c`: components from 3j symbols, g from LS coupling or a given g_eff, no mixing between the J levels of a term), which fails where the Zeeman splitting approaches the fine-structure splitting (He I 1083 nm and D3, Na D and Li I in strong fields, many molecular lines). Needs term-level atom input (L, S and the fine-structure energies of each term), diagonalization of the magnetic Hamiltonian of each term per M (LL04 Sect. 3.4), and line components with positions and strengths from the eigenvectors, in `zeeman.c`/`profile.c` and for the PRD and Stokes paths. `paschen.c` is an unfinished start, never called (`Paschen_Back()` is an empty stub); its matrix elements have bugs found in the review of #199 (`eigenval[]` accumulated on uninitialized memory, rank 2 instead of 1 in the off-diagonal 3j symbol, `E0[j]` indexed by position instead of by J) and would need rewriting | open |
| [F-4](https://github.com/dwijn/RH/issues/228) | medium | Single-precision, GPU-capable core: formal solver, opacities, line profiles and rate terms in float (consumer GPUs run FP64 at 1/32–1/64 of FP32, Apple GPUs have none). Feasible per the 2026-11-06 study (§2.4), with the conditions found there: a hybrid I / I − S formal solver, rate terms formed from D = I − S and 1 − Ψχ, every source function and emissivity computed once and shared by the formal solution and the rate equations, per-wavelength intensity scaling, cell widths instead of heights, statistical equilibrium and rate accumulators in double | open |
| [F-5](https://github.com/dwijn/RH/issues/282) | low | Structured input format for the run configuration (YAML or TOML), replacing the flat `KEY = value` of `keyword.input` and the separate `atoms.input`, `molecules.input` and `ray.input` lists. The keyword set keeps growing and has structure the flat file cannot show: solver blocks (`NG_*`, `PRD_*`, now `ANDERSON_*`, valid only with one `SOLVER`), per-atom options spread over two files, and conflicts that `readInput` has to check by hand. A schema would give names, types, defaults, allowed values, conflicts and the keyword documentation from one place, and nested sections would make most conflicts impossible to write. Needs: a reader for the old files (or a converter) so existing runs keep working; a small C parser (vendored or libyaml); a decision on scope (atom, molecule and atmosphere files stay as they are). Raised 2026-11-09 with the S-3 keywords | open |
| [F-6](https://github.com/dwijn/RH/issues/329) | low | Mg I electron collisions from Barklem et al. (2017, A&A 606, A11; VizieR J/A+A/606/A11): effective collision strengths from two R-matrix calculations (BSR and CCC, which agree to a median 0.02 dex) for all transitions among the 25 lowest Mg I terms (to 6p ¹P°), 500 K to 20 kK. The Mg I part of `MgI+II.atom` uses older CE rates. The atom groups ¹F and ³F of 4f and 5f into single g = 28 terms, so Barklem's pairs would be summed; the levels above 6p keep their current rates | open |
| [F-7](https://github.com/dwijn/RH/issues/330) | low | Inelastic collisions with neutral hydrogen including charge transfer with H⁻ (ion-pair production X + H → X⁺ + H⁻ and mutual neutralisation X⁺ + H⁻ → X + H), as in Barklem et al. (2012, A&A 541, A80) for Mg I + H. The excitation part fits the existing `CH` keyword; the charge-transfer part does not: `CH0` and `CH+` are charge exchange with H and H⁺, and the H⁻ density is not an atom level. It needs a keyword that takes n(H⁻) from the chemical equilibrium (`nHmin`, B-40). In the Mg+H data mutual neutralisation into 3s4s ¹S (about 5e-8 cm³ s⁻¹) is the largest rate, and charge transfer was the H process that mattered in Li and Na | open |

### 3.8 Solver redesign

Work items of the solver redesign; [`SOLVER.md`](https://github.com/dwijn/RH/blob/modernize/docs/SOLVER.md) §5 has goal,
interface, traps and acceptance criteria for each. S-1 and S-2 come first
(prerequisites for a fair comparison), S-5 is a decision gate, S-6 to S-9
depend on it. The N issues listed under *Addresses* stay open until the S
item that resolves them lands.

| ID | Sev | Item | Addresses | Status |
|---|---|---|---|---|
| [S-1](https://github.com/dwijn/RH/issues/117) | high | Residual norm, reference solutions and metrics independent of the solver: true statistical-equilibrium residual ‖F_SE‖ (formal solution with fixed populations), a residual norm, tight reference solutions for the S-5 cases, and reports of formal solutions, wall time and true error in populations and spectra | N-3 | fixed |
| [S-2](https://github.com/dwijn/RH/issues/118) | high | One MALI step as a deterministic, smooth function x ↦ G(x) with an explicit state vector (populations, PRD ρ, background-scattering J, ne): no Ng inside, fixed PRD sub-iteration count, no hidden history, deterministic at a fixed thread count. Split (2026-11-07) into S-2 (ne fixed) and S-2b (ne in the state, `Background()` split into setup and recompute), delivered in one pull request. B-38, B-39 and the `Background()` part of B-11, which G depends on, are fixed in their own pull requests first | N-4 (early exit) | fixed |
| [S-3](https://github.com/dwijn/RH/issues/119) | high | Stabilized Anderson acceleration replacing Ng: QR least squares, scaled variables, one coefficient set for the whole state, safeguard with fallback to the plain step, positivity, history reset after ne/HSE/relaxation updates (`SOLVER = ANDERSON`) | N-2, N-8 | fixed |
| [S-4](https://github.com/dwijn/RH/issues/120) | high | Newton–GMRES on R(x) = G(x) − x: finite-difference Jacobian products (one MALI step each), Eisenstat–Walker forcing, backtracking line search, stop on ‖F_SE‖ (`SOLVER = NEWTON_KRYLOV`) | N-8 | fixed |
| [S-5](https://github.com/dwijn/RH/issues/121) | high | Experiment and decision gate: Ng vs Anderson vs Newton on a fixed 1-D case set (easy, slow, failing; then 2-D and RMHD), success criteria fixed beforehand, measured in formal solutions and true error. Candidate for the operator: with the quadratic upwind slope of N-24 the upwind control point carries 2/3 of S_c, and Ψ = β + θ + (2/3)γ took `2d_pert_h6_ca_mg` from 94 to 67 iterations but a 325-point perturbed FAL-C from 44 to 76 (#190) | N-10 (as test case) | fixed |
| [S-6](https://github.com/dwijn/RH/issues/122) | medium | Analytic residual F = Γ(n, J[n])·n (Arramy et al. 2024) with MALI's Γ + Ψ* as explicit preconditioner; only if S-5 favours Newton | — | open |
| [S-7](https://github.com/dwijn/RH/issues/123) | medium | Monolithic electron density (charge conservation in the residual) and PRD (ρ in the residual or CRD-preconditioned inner solve); only if S-5 favours Newton | N-4, N-6 | open |
| [S-8](https://github.com/dwijn/RH/issues/124) | medium | Pseudo-transient continuation and physics continuation (CRD → +ne → +PRD → finer grid), warm starts from neighbouring columns or snapshots | N-9 (partly) | open |
| [S-9](https://github.com/dwijn/RH/issues/125) | low | Mesh-independent convergence: multigrid preconditioner or nonlinear multigrid (FAS), for fine-grid RMHD columns | N-14 (iteration count only) | open |
| [S-10](https://github.com/dwijn/RH/issues/293) | medium | `ANDERSON` the default solver of `rhf1d` and `rhsphere`; `rhsc2d` and `rhsc3d` stay on `NG` until S-11. Every keyword that is set is checked against the solver that runs, independent of the atoms (an unset `PRD_ITER_LIMIT` is 0 under a G solver); keyword files with Ng keywords and no `SOLVER` stop with a generic message, no fallback to `NG`; the default is printed; shipped 1-D keyword files without Ng keywords; test tools that mean Ng pin `SOLVER = NG`. Changes results (SOLVER.md, Anderson as the default solver, decided 2026-11-10) | — | fixed |
| [S-11](https://github.com/dwijn/RH/issues/294) | medium | Anderson's history in 2-D and 3-D: the state is 60 % J and 39 % PRD ρ in the shipped examples, so 40 history vectors take about 24 GB on the 2-D perf case (10.7 GB run) and about 17 TB on the 3-D estimate of §3.3. Extrapolate only the populations (perhaps ρ) and carry J from the last G, as Ng does; measure on the S-5 set in 1-D, then on `2d_pert_h6_ca_mg` after N-33. Prerequisite of `ANDERSON` as the multi-D default (SOLVER.md §9 Q1). `ANDERSON` with `POPULATIONS_RHO` is the default in 2-D and 3-D; 1-D keeps `ALL` (decided 2026-11-13). Changes results | — | fixed |
| [S-12](https://github.com/dwijn/RH/issues/306) | low | Decreasing-envelope safeguard for Anderson: with c = 1 (`ANDERSON_SAFEGUARD`; a point is accepted only if its residual is at most that of the newest accepted pair, `anderson.c:391`) `shock_h6` rejects 78 of 282 points at m = 5 (31 of 146 at m = 10), and m = 5 then needs 12 % more G evaluations than plain G (580 against 517; the rejected points' evaluations count). Try the decreasing envelope of Zhang, O'Donoghue & Boyd (2020), which also accepts a point whose residual is under a slowly shrinking bound; keep it only if it beats plain G on `shock_h6` without costing the easy controls of S-5 | SOLVER.md, Safeguards (S-3) and S-5 Results | open |
| [S-13](https://github.com/dwijn/RH/issues/307) | low | A margin in Anderson's inversion guard: the check for new population inversions (line test on n/g, continuum test with n/n* and the Boltzmann factor at the edge) damps or replaces extrapolated points from the first extrapolations on: on `falc_helarge` 105 damped and 44 plain steps, the same at m = 3, 5 and 10; also often on `falc_civi`, `t5750_h6` and `falc_mgi66`, though Anderson still wins clearly there. Cause not measured: first find which test fires, for which transitions and how far past the threshold (levels with nearly equal n/g, or n/n* across an edge, would flip on tiny changes); then try a small relative margin and measure guard actions and work on the S-5 set | SOLVER.md, Safeguards (S-3) and S-5 Results | open |
| [S-14](https://github.com/dwijn/RH/issues/339) | low | Drop the line test of Anderson's trial-point guard (`guardTrial`: a trial point must not invert a line at a depth where G's point does not), keeping the finite and negative-value checks. The test was meant to keep a strong inverted line from making the total opacity negative at its centre (SOLVER.md trap 5), but it does not: the total opacity goes negative at G's own points (the ZERO_RADIATION starts of `falc_helarge`, `falc_civi`, `falc_mgi66`; down to −44 χ_c on `falc_helarge`) and at `falc_helarge`'s converged solution (504 points, down to −6.6 χ_c, IR He I 4296 nm), and the deepest negative slant optical depth summed over consecutive intervals is −0.094 on `falc_mgi66` with the test and −0.026 without, far from the −709 where `exp` overflows. Without it all 11 S-5 cases converge with the same results; work changes by 0 to −8 % (fewer formal solutions on `falc_helarge`, `falc_mgi66`, `falc_civi`). Insurance against the overflow is N-40. Decided by the user 2026-11-13 | SOLVER.md S-13; `claude-tmp/s1/s13/chi` (README) | open |

### 3.9 Bugs also present in upstream RH

Bugs found here that also exist in upstream RH (han-uitenbroek/RH,
c95bdf6), to be reported there whatever becomes of this version: a wrong
answer that looks believable is worse than a crash, and the first one
here (B-24) went unnoticed until a test with an analytic solution. Every
new bug that also exists upstream gets a row; *Reported* records the
upstream issue or pull request once filed. Most severe first.

| ID | Sev | What happens upstream | Fixed here | Reported |
|---|---|---|---|---|
| [B-24](https://github.com/dwijn/RH/issues/149) | high | 2-D and 3-D long characteristics stop one section short: wrong intensities wherever short characteristics cross column lines or vertical planes before the next layer (deep or optically thick layers, fine horizontal grids). 7–11% error for a linear source function on a uniform slab; the shipped 2-D/3-D examples change by up to 25% at short wavelengths | [#152](https://github.com/dwijn/RH/pull/152) | — |
| [B-25](https://github.com/dwijn/RH/issues/155) | high | The polarized DELO-Bezier solvers (all geometries) pair the derivative terms of the control points with the wrong weights and have the wrong sign of dK/dτ at the central point: first order where K varies with depth. Converged Ca II full Stokes on FALC: Stokes V 38% from the converged answer, the 854.2 nm V amplitude 21% too large | [#161](https://github.com/dwijn/RH/pull/161) | — |
| [B-1](https://github.com/dwijn/RH/issues/2) | high | `SOLVE_NE = ITERATION` aborts on the second `Background()` call with active atoms | [#69](https://github.com/dwijn/RH/pull/69) | — |
| [B-3](https://github.com/dwijn/RH/issues/4) | high | Collisional relaxation stops after one iteration and returns LTE populations reported as converged | [#70](https://github.com/dwijn/RH/pull/70) | — |
| [B-2](https://github.com/dwijn/RH/issues/3) | high | `PRD_NG_ORDER > 0` with the default `PRD_NG_PERIOD = 0` divides by zero | [#67](https://github.com/dwijn/RH/pull/67) | — |
| [B-4](https://github.com/dwijn/RH/issues/5) | high | Atom data with zero or non-finite values crash with SIGFPE at setup (e.g. `CaII_large`) | [#71](https://github.com/dwijn/RH/pull/71) | — |
| [B-29](https://github.com/dwijn/RH/issues/172) | high | Ca II with PRD in a moving realistic atmosphere (a Bifrost snapshot): the PRD sub-iterations diverge, populations and opacity go negative and the run dies with a floating-point exception (1-D: 35 of 126 columns of a slice; 2-D and 3-D cuts) | [#187](https://github.com/dwijn/RH/pull/187) | — |
| [T-4](https://github.com/dwijn/RH/issues/179) | high | Hybrid PRD (`PRD_ANGLE_APPROX`) on more than one thread loses updates of the gas-frame mean intensity (a data race across wavelengths): results and convergence depend on thread timing | [#182](https://github.com/dwijn/RH/pull/182) | — |
| [N-17](https://github.com/dwijn/RH/issues/131) | medium | Cubic Bezier weights (default solver, all geometries) wrong by up to 1.2e-5 below dt = 0.05 (truncated series) and cancelling above; I off by ~1e-6 | [#147](https://github.com/dwijn/RH/pull/147) | — |
| [B-5](https://github.com/dwijn/RH/issues/6) | medium | Linear short characteristics compute the last point of each 1-D ray with stale weights | [#147](https://github.com/dwijn/RH/pull/147) | — |
| [B-16](https://github.com/dwijn/RH/issues/93) | medium | Collision data before a `TEMP` line read an uninitialized item count (crash) or inherit one (silently wrong rates) | [#99](https://github.com/dwijn/RH/pull/99) | — |
| [B-26](https://github.com/dwijn/RH/issues/156) | medium | The 2-D and 3-D polarized solvers compute the last layer of each ray, which includes the emergent Stokes vector, with a stale (2-D) or uninitialized (3-D) source vector | [#161](https://github.com/dwijn/RH/pull/161) | — |
| [B-28](https://github.com/dwijn/RH/issues/167) | low | A magnetic field exactly along an oblique ray gives NaN polarized profiles (0/0 in the field azimuth) | [#169](https://github.com/dwijn/RH/pull/169) | — |
| [B-27](https://github.com/dwijn/RH/issues/166) | low | A J = 0 -> 0 line in a full-Stokes run gives NaN profiles (Zeeman strengths normalized by zero) | [#169](https://github.com/dwijn/RH/pull/169) | — |
| [N-23](https://github.com/dwijn/RH/issues/143) | medium | 1-D static atmospheres are always solved with second-order Feautrier; `S_INTERPOLATION` is silently ignored there. Converged on FALC (82 points), Feautrier is up to 20% from the converged spectrum, Bezier3 3% | [#154](https://github.com/dwijn/RH/pull/154) | — |
| [B-6](https://github.com/dwijn/RH/issues/7) | high | The preconditioned rate equations converge to populations that violate statistical equilibrium wherever active transitions share a lower level, form a ladder or overlap across ionization stages: `Ieff` subtracts more emissivity than Gamma puts back (Mg II 4s and 3d populations off by up to 34 %, H by 1e-3) | [#211](https://github.com/dwijn/RH/pull/211) | — |
| [B-7](https://github.com/dwijn/RH/issues/8) | medium | 3-D `Jgas` possibly double counted (not yet confirmed) | — | — |
| [B-8](https://github.com/dwijn/RH/issues/9) | medium | `calcDistance()` (3-D) can fall off the end without a return value | [#76](https://github.com/dwijn/RH/pull/76) | — |
| [B-9](https://github.com/dwijn/RH/issues/10) | medium | Scratch-file `pread`/`pwrite` errors mishandled (signed/unsigned comparison) | [#76](https://github.com/dwijn/RH/pull/76) | — |
| [B-10](https://github.com/dwijn/RH/issues/11) | medium | 32-bit overflow of `Nlambda*Nspace` sizes on large grids | — | — |
| [N-18](https://github.com/dwijn/RH/issues/132) | low | `w2`/`w3` weights wrong by up to 2e-6 near dtau = 5e-4 | [#147](https://github.com/dwijn/RH/pull/147) | — |
| [B-20](https://github.com/dwijn/RH/issues/128) | low | `VoigtHui()` returns NaN at a = v = 0 (angle-dependent PRD with zero damping) | [#130](https://github.com/dwijn/RH/pull/130) (Voigt replaced, [#146](https://github.com/dwijn/RH/pull/146)) | — |
| [B-17](https://github.com/dwijn/RH/issues/94) | low | `writeMULTIatmos()` (tools): uninitialized heights for a geometric scale; velocities 1000× too large | [#101](https://github.com/dwijn/RH/pull/101) | — |
| [B-11](https://github.com/dwijn/RH/issues/12) | low | Memory leaks in `Background()` and 3-D Stokes setup | [#253](https://github.com/dwijn/RH/pull/253) | — |
| [B-12](https://github.com/dwijn/RH/issues/13) | low | AR85-CHP collisions assume the 6th hydrogen level is the proton | — | — |
| [B-14](https://github.com/dwijn/RH/issues/15) | low | `error.c` can overflow its 512-byte message buffer | — | — |
| [B-15](https://github.com/dwijn/RH/issues/80) | low | `Atoms/CaII_large.atom` contains `Inf` collision strengths | [#203](https://github.com/dwijn/RH/pull/203) | — |
| [B-18](https://github.com/dwijn/RH/issues/95) | low | `solve1d` does not validate its boundary-condition codes | [#196](https://github.com/dwijn/RH/pull/196) | — |
| [B-19](https://github.com/dwijn/RH/issues/113) | low | `StopRequested()` never closes the file it opens | [#197](https://github.com/dwijn/RH/pull/197) | — |
| [B-21](https://github.com/dwijn/RH/issues/135) | low | `BiLinear()` reads past its table at the upper edge | [#198](https://github.com/dwijn/RH/pull/198) | — |
| [B-22](https://github.com/dwijn/RH/issues/136) | low | `w3js()` accepts invalid arguments and overflows from j = 50 | [#199](https://github.com/dwijn/RH/pull/199) | — |
| [B-23](https://github.com/dwijn/RH/issues/137) | low | `Hunt()` one index low on exact node values (harmless for interpolation) | [#200](https://github.com/dwijn/RH/pull/200) | — |
| [B-30](https://github.com/dwijn/RH/issues/173) | low | `OH_bf_opac()`/`CH_bf_opac()` start `Hunt()` from uninitialised indices (undefined behaviour; results unaffected) | [#176](https://github.com/dwijn/RH/pull/176) | — |
| [B-32](https://github.com/dwijn/RH/issues/184) | medium | Hybrid PRD on large 3-D grids: the gas-frame interpolation table is sized and indexed with `int` arithmetic, which overflows | [#189](https://github.com/dwijn/RH/pull/189) | — |
| [B-31](https://github.com/dwijn/RH/issues/183) | low | Hybrid PRD with a single ray accumulates the gas-frame J over iterations (never zeroed); in 2-D and 3-D it is 2 J at wavelengths without PRD lines | [#188](https://github.com/dwijn/RH/pull/188) | — |
| [B-34](https://github.com/dwijn/RH/issues/215) | medium | 2-D and 3-D `solveray`: a derived direction cosine is NaN when the given ones' squares sum to just over 1 (no clamp, no check); in 2-D every polarized result is NaN, silently | — | — |
| [N-29](https://github.com/dwijn/RH/issues/216) | low | 2-D long characteristics: columns of a horizontally uniform periodic slab differ by up to 1e-2 of I (2e-3 here after B-24) | — | — |
| [N-31](https://github.com/dwijn/RH/issues/262) | medium | The rate equations are solved by LU, which loses componentwise accuracy for levels tightly coupled by collisions: populations of the `He_large` fine-structure multiplets are wrong by up to 4e-3, which stalls its convergence | [#265](https://github.com/dwijn/RH/pull/265) | — |
| [N-32](https://github.com/dwijn/RH/issues/263) | medium | Tabulated atomic data (collision rates in T, bound-free cross-sections in wavelength) are interpolated with natural cubic splines in linear coordinates: inaccurate for rates that vary as power laws, and negative where a table is non-monotone or has resonances (negative collisional and photoionization rates) | [#279](https://github.com/dwijn/RH/pull/279) | — |
| [B-33](https://github.com/dwijn/RH/issues/204) | low | A run that uses less than one clock tick (10 ms) of CPU dies with SIGFPE at the very end: `printTotalCPU()` divides by the total CPU time, 0/0 with FP traps enabled. Hits tiny runs on fast machines | [#224](https://github.com/dwijn/RH/pull/224) | — |
| [B-35](https://github.com/dwijn/RH/issues/221) | low | `rhanalyze` returns the Stokes-V emissivity as the background scattering opacity at polarized wavelengths (skips 3 records instead of 4) | [#222](https://github.com/dwijn/RH/pull/222) | — |
| [B-36](https://github.com/dwijn/RH/issues/235) | low | Collision keyword `SUMMERS` (density-dependent dielectronic recombination): for an ion with 55 to 86 electrons the periodic-table row in `rowcol` is never set, so the density factor uses an uninitialized value | [#238](https://github.com/dwijn/RH/pull/238) | — |
| [B-37](https://github.com/dwijn/RH/issues/241) | low | `chartflow` adds two never-written arrays instead of the active opacity and emissivity at wavelengths without angle-dependent opacity, so its flow maps leave out the active lines whenever `malloc` returns zeroed memory | [#238](https://github.com/dwijn/RH/pull/238) | — |
| [B-38](https://github.com/dwijn/RH/issues/243) | low | Active molecules with overtone lines (Δv ≥ 2): their radiative rates accumulate in Gamma over the iterations, because `initGammaMolecule` resets only the Δv = 1 entries | [#251](https://github.com/dwijn/RH/pull/251) | — |
| [B-40](https://github.com/dwijn/RH/issues/247) | low | The stored H⁻ density (`nHmin`, H⁻ bound-free opacity and `molecules.out`) uses all atomic hydrogen instead of the neutral part that the chemical-equilibrium equations use: too large by 1/fn0 where hydrogen is ionized | [#257](https://github.com/dwijn/RH/pull/257) | — |
| [B-41](https://github.com/dwijn/RH/issues/255) | low | A terminating error message longer than about 460 characters (with the routine name) overflows the 512-byte stack buffer in `ErrorExit`; latent, no current message is that long | — | — |
| [B-42](https://github.com/dwijn/RH/issues/258) | medium | `C_I+II_9.atom`: two C II fine-structure collision strengths (including the [C II] 158 µm line) have the node temperature, 9000, in place of the value at 9000 K; the interpolated rates ring from 1 to 19 kK, up to 4000 times too large or negative | [#275](https://github.com/dwijn/RH/pull/275) | — |
| [B-43](https://github.com/dwijn/RH/issues/259) | medium | `CaII.atom` and copies: Ca II collision strengths for 4s–3d and the 3d and 4p fine structure are up to 9 times below the R-matrix values of Meléndez et al. (2007); the interpolated 3d–3d rate is negative between 50 and 100 kK | [#276](https://github.com/dwijn/RH/pull/276) | — |
| [B-44](https://github.com/dwijn/RH/issues/260) | low | `MgII.atom` and copies: 9 of 45 Mg II collision strengths disagree with Sigut & Pradhan (1995), the source the file names: 4f–3d ten to a hundred times too small (the interpolated rate negative near 24 kK), 4s–3p, 5s–3p and 5s–3d not split by g(J) (1.5 to 3 times too large) | — | — |
| [B-45](https://github.com/dwijn/RH/issues/298) | medium | A negative `NG_ORDER` is accepted and corrupts the Ng allocation; the `rhf1d` example with `NG_ORDER = -2` stops in `readBackground` with "Error reading file" | [#302](https://github.com/dwijn/RH/pull/302) | — |
| [B-48](https://github.com/dwijn/RH/issues/316) | medium | `rhanalyze` cannot be imported on Python 3.13 and later unless the undeclared outside package `mda_xdrlib` is installed (`xdrlib` was removed from the standard library) | — | — |
| [B-49](https://github.com/dwijn/RH/issues/322) | medium | `rhanalyze.rhout` cannot read any 3-D output (`THREEE_D_PLANE` typo in `rhatmos.py`; UnboundLocalError) | — | — |
| [B-50](https://github.com/dwijn/RH/issues/323) | medium | `input_atmos.write` ignores its `Bfile` argument and overwrites the B file the atmosphere was read from | — | — |
| [B-51](https://github.com/dwijn/RH/issues/324) | low | `rhanalyze.rhout` cannot read the output of a run with an active molecule (NameError on `rhdir`; `molecules.out` and `pops_mol.*.out` read in a layout the C code does not write) | — | — |
| [B-52](https://github.com/dwijn/RH/issues/335) | high | 3-D long characteristics of the upward sweep read S and χ at the horizontal mirror image of their crossings with the vertical planes (wrong emergent intensities in structured 3-D atmospheres) | [#344](https://github.com/dwijn/RH/pull/344) | — |
| [B-53](https://github.com/dwijn/RH/issues/336) | medium | 3-D long characteristics of rays with a negative μx or μy cross the vertical planes one column or row off (a local `round` that truncates toward zero) | [#343](https://github.com/dwijn/RH/pull/343) | — |
| [N-34](https://github.com/dwijn/RH/issues/268) | low | Collision rates from formulas (`SHULL82`, `BADNELL`, `AR85-*`, `BURGESS`) are not checked for sign or for the validity range of their fits; `AR85-CHP` takes the protons from hydrogen level 5, right only for a six-level model | — | — |
| [N-35](https://github.com/dwijn/RH/issues/277) | low | H⁻ bound-free cross-section from Geltman (1962): 1–3% below the Wishart (1979) values (John 1988 fit) in the visible, 6–19% below between 300 and 125 nm, and without the threshold law near 1.64 µm | — | — |
| [N-37](https://github.com/dwijn/RH/issues/320) | medium | 2-D: the monotonic quadratic interpolation between columns (`Quadr`) switches to linear where the parabola overshoots, a jump that makes the iteration a discontinuous map. On a perturbed FAL-C (32 columns, H + Ca II + Mg II, all CRD) the populations end in a cycle with a maximum change of 10⁻⁴ (period 2) and never converge; the scattering sweeps with fixed populations alternate at 2.6×10⁻³. Clamping the parabola instead (tested here) converges | — | — |
| [N-39](https://github.com/dwijn/RH/issues/333) | medium | 3-D: on a FAL-C perturbed over 8 × 8 columns (H_6 and Ca II with PRD), the lambda iteration of the background scattering at the start does not converge but ends in a period-2 cycle, with relative changes of J of 50–490 per sweep at the Lyα and Hα cores, and the first iteration changes the populations by up to 3.7×10³. Here the cycle grows and the run stops with SIGFPE. Cause: the parabolic integration along the long characteristics, whose weight of the downwind S is negative for a thick upwind section and about −1/Δτ_dw for a thin downwind one; with linear sections the iteration converges | — | — |
| [N-40](https://github.com/dwijn/RH/issues/338) | medium | A negative total opacity (a strong inverted line) makes the formal solvers overflow; the run dies with SIGFPE or carries inf and NaN, with no message saying where | — | — |

---

## 4. Plan

Ordered so that every step can be verified against upstream: byte-identical
(`run_tests.sh --exact`) for pure refactors, tolerance mode for numerics, and
`run_stress.sh` for convergence and performance.

1. **Cheap, high-value fixes** — crash and silent-wrong-answer bugs (B-1 to
   B-4, B-8, B-9), convergence status reporting (N-1), clamp/fallback for
   parabolic Ψ* (N-5) and switch the 2-D example to Bezier.
2. **Build and hygiene** — drop Fortran (Q-2, Q-11); new build system with
   automatic dependencies and safe parallel builds (Q-1); warnings on and
   cleaned (Q-3); debug info in default builds (Q-15); dead code removed (Q-10); CI with the regression suite.
   *Done 2026-10-16:* CMake, CI on Linux/macOS with gcc and clang plus
   sanitizers, zero warnings with `-Werror` in CI. Q-11 (Voigt variants)
   moved to step 3, where the Voigt unit tests decide which to keep.
3. **Unit tests** for Voigt (against a Faddeeva reference; then reduce
   the Voigt variants, Q-11), LU, splines, formal solvers on analytic
   slabs — before touching numerics. Unit tests use problems with known
   answers only. Realistic multi-D atmospheres from RMHD snapshots (Q-18)
   belong to the stress tier, not the unit tests; check the data's
   redistribution terms first, and commit a fetch script with checksums
   instead of the data if needed.
   *Two test tiers from here on:* regression (`run_tests.sh`, must pass on
   every pull request; a non-default solver is checked on converged
   solutions against stored references, not on iteration 4) and stress
   (`run_stress.sh`, convergence and performance, including the RMHD
   cases). See [`SOLVER.md`](https://github.com/dwijn/RH/blob/modernize/docs/SOLVER.md) §6.
   *Progress 2026-10-17:* unit tests run by `ctest` (also under
   ASan/UBSan in CI) for Voigt, the numerical kernels and the 1-D formal
   solvers; one Voigt function replaces five (Q-11); Bezier and linear
   weights accurate to 1e-12 (B-5, N-17, N-18). CI now compares each
   change with its base commit (label `changes-results` for intended
   changes) and reports drift from upstream. Remaining: 2-D/3-D solvers on
   horizontally uniform slabs against the 1-D solutions, the spherical
   Feautrier solver, RMHD stress cases (Q-18).
   *Progress 2026-10-18:* multi-D and spherical solvers tested (B-24
   found and fixed); polarized solvers tested against the matrix form
   of the analytic solution (B-25, B-26 found); 1-D static atmospheres
   follow `S_INTERPOLATION` (N-23). The test problems and their
   derivations are written up in `docs/unit-tests/` (Q-20). Remaining:
   RMHD stress cases (Q-18).
   *Done 2026-10-21:* RMHD stress cases on Bifrost cuts (Q-18), unit
   tests and CI (Q-14) closed; the small findings of the unit tests fixed
   (B-18, B-19, B-21 to B-23, N-20, N-22); N-21 (`SolveLinearEq`) moved
   to the S track.
4. **Solver redesign** (S track, §3.8, [`SOLVER.md`](https://github.com/dwijn/RH/blob/modernize/docs/SOLVER.md)) replaces
   the earlier list of convergence patches (Ng safeguards, PRD
   sub-iteration rework, ne preconditioning, positivity): S-1 and S-2,
   then S-3 and S-4, then the S-5 decision gate; S-6 to S-9 depending on
   its outcome. Independent of the solver and still part of this step:
   rate-matrix scaling (N-7; the preconditioner needs it), starting
   solutions (N-9), monotonic 3-D interpolation (N-15), a depth-resolution
   test (N-14), N-12, N-13, N-16 and B-5 to B-7. N-10 becomes a test case
   for the new solvers. G must not keep hidden state, and with ne in the
   state it calls `Background()` once per evaluation, so B-38, B-39 and
   the `Background()` part of B-11 are fixed first, each in its own pull
   request (bug fixes never ride along with the solver work, so that
   they survive if it is reverted).
5. **Speed** — *Started early, 2026-11-04:* the background opacity stage
   (P-13, T-5) blocks large-line-list UV syntheses and does not depend on
   the solver, so it goes first: shared state removed and per-point tables
   (step 1), then threaded over wavelength (step 2).
   Hoist angle-independent work (P-1), 3-D interpolation (P-3),
   dynamic scheduling (P-4), flags (P-7). Parallel statistical equilibrium
   (P-2), the Γ merge (P-5) and fused accumulation (P-6) depend on how Γ is
   built and used by the new solver: decide after S-5.
6. **Memory** — profiles on the fly (M-1), PRD weights tabulated instead
   of per point (M-4), float storage where accuracy allows. Sparse Γ/C and
   no per-thread copies (M-2, M-3) after S-5, together with the storage of
   Krylov/Anderson history vectors (about 20 GB per vector in the 3-D
   estimate of §3.3; [`SOLVER.md`](https://github.com/dwijn/RH/blob/modernize/docs/SOLVER.md) §4.3).
7. **Structure** — shared geometry layer to remove duplication (Q-5),
   self-contained headers, reduced global state, error handling (Q-4, Q-7,
   Q-9, Q-16). Q-7 comes first: the code is to be ported into STIC, which
   calls RH as a library, so errors must return to the caller instead of
   calling `exit()`. MPI domain decomposition (P-11) is a separate, larger
   project.
   The manual is revised and rewritten with its source in the repository
   once the code has settled (Q-19).

---

## 5. Changelog

### 2026-11-10 — Negative NG_ORDER, NG_DELAY, PRD_NG_ORDER and PRD_NG_DELAY stop the run at input instead of crashing it (B-45)

readInput accepted negative Ng orders, and the run then crashed: `NG_ORDER = -2` with SIGFPE in the first iteration (`Accelerate` took the iteration count modulo `NG_ORDER + 2`), `PRD_NG_ORDER = -2` with heap corruption (`NgInit` allocated a history with `order + 2` = 0 rows; below -2 the count is negative); upstream c95bdf6 with `NG_ORDER = -2` stops in `readBackground` with "Error reading file". readInput now stops with an error when `NG_ORDER`, `NG_DELAY`, `PRD_NG_ORDER` or `PRD_NG_DELAY` is negative, under every solver; negative delays were harmless but meaningless. `NgInit` itself also stops on a negative order (or a period below 1 with a positive order), for the callers that do not go through readInput (`solve1d`, `solve2d`). The e2e test `solver_g` checks that all four keywords, negative, stop in readInput. Runs with valid values are unchanged.

### 2026-11-10 — The G solvers accept PRD_N_MAX_ITER < 1 (every line in CRD) (B-46)

`SOLVER = FIXED_POINT`, `ANDERSON` and `NEWTON_KRYLOV` refused `PRD_N_MAX_ITER` < 1 at input ("the map needs at least one PRD sub-iteration"). But `PRD_N_MAX_ITER` < 1 is RH's switch for treating every line in CRD: `readAtom` then marks no line PRD, and G has no PRD sub-iterations, so the input is valid and the check is removed. `PRD_ITER_LIMIT` other than 0 and `PRD_NG_ORDER` other than 0 are still refused under these solvers, also in a run without PRD lines. The e2e test `solver_g` now also runs `NG` and `FIXED_POINT` with `PRD_N_MAX_ITER = 0`: the two agree bit for bit after 5 iterations, and their populations differ from those of the runs with PRD. Runs with `SOLVER = NG` are unchanged.

### 2026-11-10 — ANDERSON is the default solver of rhf1d and rhsphere; every keyword that is set is checked against the solver (changes results, breaking) (S-10)

`rhf1d` and `rhsphere` now iterate with `SOLVER = ANDERSON` (stabilized Anderson acceleration of the map G, S-3, depth 10) unless the keyword file sets `SOLVER` (the other 1-D and spherical programs that read keyword.input, such as `solveray`, check the keywords against the same default); `rhsc2d` and `rhsc3d` keep `NG` until Anderson's history fits in memory there (S-11). Decided 2026-11-10 after S-5 (docs/SOLVER.md, Anderson as the default solver). Every run that iterates prints its solver when the iteration starts (" -- SOLVER = ANDERSON (the default)"). On the shipped examples Anderson converges in 10 iterations where Ng needed 11 (`rhf1d`, ITER_LIMIT 1e-2) and 11 where Ng needed 15 (`rhsphere`); 21 and 23 against 24 and 49 at 1e-5. Results at the shipped ITER_LIMIT of 1e-2 move by 1–3 % in the populations, about each solver's distance from the converged solution; converged (1e-5) the two agree to 8e-5.

**Breaking:** a keyword file that sets keywords of the Ng iteration and no `SOLVER` now stops at input, with one warning per keyword ("NG_DELAY is set, but SOLVER = ANDERSON (the default) does not use it. Set SOLVER = NG or remove the keyword.") and a closing error; there is no silent fallback to `NG`. Add `SOLVER = NG` to keep the old iteration (results unchanged bit for bit), or remove `NG_ORDER`, `NG_DELAY`, `NG_PERIOD`, `PRD_ITER_LIMIT` and `PRD_NG_*` to use Anderson. `COLL_RELAX_*`, `HYDROSTATIC = TRUE`, `LIMIT_MEMORY = TRUE`, and with `SOLVE_NE = ITERATION` also `OLD_BACKGROUND = TRUE`, `NON_ICE = TRUE` and `BACKGROUND_DOUBLE = FALSE`, need `SOLVER = NG`. The shipped `rhf1d` and `rhsphere` keyword files no longer set the Ng keywords and have a commented `SOLVER` line with the choices.

Keyword checks (step 1, item 2): every keyword that is set is checked against the solver that runs, and no check depends on the atoms; a default is never a conflict. Under `FIXED_POINT`, `ANDERSON` and `NEWTON_KRYLOV` an unset `PRD_ITER_LIMIT` is now 0 (`PRD_N_MAX_ITER` sub-iterations per G), where its default 1.0E-2, which belongs to `NG`, used to stop the run; set to anything but 0 it still stops the run. Keywords without an "off" value stop the run when they are set at all: `NG_DELAY` and `NG_PERIOD` under `ANDERSON` and `NEWTON_KRYLOV` (`FIXED_POINT` keeps the Ng keywords for its "Ng on G" control), and `PRD_NG_DELAY` and `PRD_NG_PERIOD` under every G solver; they used to be ignored silently. Keywords with an "off" value are accepted at that value, which states what the solver does: `NG_ORDER`, `PRD_NG_ORDER`, `PRD_ITER_LIMIT` and `COLL_RELAX_NSTEP` at 0, and `COLL_RELAX_FACTOR` at 1 (it used to be ignored at any value); other values stop the run. `-showkeywords` prints the `PRD_ITER_LIMIT` in use. Runs with `SOLVER = NG` are unchanged bit for bit. The test tools can remove a keyword from a shipped keyword file: `run_stress.sh --kw '-NG_DELAY'` (and in case keywords), `setkw(text, key, None)` in the e2e tests, whose Anderson and Newton runs now remove `NG_DELAY` and `NG_PERIOD`. New e2e checks: `FIXED_POINT` without `PRD_ITER_LIMIT` in the file repeats the history of `NG` with `PRD_ITER_LIMIT = 0` bit for bit; `NG_DELAY` with `ANDERSON`, `NG_PERIOD` with `NEWTON_KRYLOV`, `PRD_NG_DELAY`, `PRD_NG_PERIOD` and `COLL_RELAX_FACTOR = 0.01` with `FIXED_POINT` stop the run at input.

Test tools (step 2): `make_references.sh` and the stress cases that test the Ng iteration (`1d_falc_nong`, `1d_falc_h6_tight`, `1d_falc_prdng`, `1d_falxco_collrel`) pin `SOLVER = NG` with the former shipped Ng settings written out; the other stress cases, `run_tests.sh` and the e2e tests run each program's default. `run_tests.sh`'s comparison with upstream in 1-D and spherical geometry now includes the change of solver.

### 2026-11-10 — ANDERSON_EXTRAPOLATE = ALL | POPULATIONS | POPULATIONS_RHO: Anderson's history over part of the state, J carried from G (S-11 step 1, default ALL)

First step of S-11 (docs/SOLVER.md, S-11): the new keyword `ANDERSON_EXTRAPOLATE` chooses which blocks of the state Anderson extrapolates. `ALL` (the default) extrapolates every block, as before. `POPULATIONS` extrapolates the populations of the active atoms and molecules, and ne when it is in the state. `POPULATIONS_RHO` also extrapolates the ρ of the PRD lines. The other blocks (J, J20, and ρ under `POPULATIONS`) are carried: every point takes them from G at the last accepted point, as RH's Ng does with J. The history, the variables, the weights and the least squares cover only the extrapolated blocks. x, G(x) and G of the last accepted point stay full, so the memory is 3·N + (3m + 7)·N_x doubles plus (m + 4)·N_x bytes (N_x doubles more with `RELATIVE`), printed at the start (history and variables, and the three full states, separately) with the extrapolated and the carried blocks. On the shipped `rhf1d` example the history and variables take 12.4 MB with `ALL`, 0.15 MB with `POPULATIONS` and 5.0 MB with `POPULATIONS_RHO`; the full states 0.96 MB in each. Before every evaluation the solver checks bit for bit that the carried blocks equal those of G at the last accepted point, and reports the count at the end; the check covers the packing and the copy of the carried blocks, not G. Like every `ANDERSON_*` keyword, it stops the run under another solver. `POPULATIONS_RHO` on atoms without PRD lines runs as `POPULATIONS`, with a warning that says so; it does not stop the run, since no keyword check depends on the atoms. Results do not change: `ALL` is the default, and `run_tests.sh --exact` against modernize 4951ebb is byte-identical. Newton-Krylov, which shares the variable code with every block packed, gives bit-identical populations, J and spectra in its e2e runs. New checks in `e2e_solver_anderson`: with depth 0, `POPULATIONS` and `POPULATIONS_RHO` repeat `FIXED_POINT` bit for bit; a `POPULATIONS_RHO` run with `RELATIVE` and `PER_BLOCK` (three solves) is bit-reproducible; the banner names the extrapolated blocks (Ca II populations, 492 unknowns in 1 block; with ρ 16 072 in 3); with forced rejections (`ANDERSON_SAFEGUARD = 0.5`) the carried blocks agree in every evaluation after the first; `ANDERSON_EXTRAPOLATE` under `FIXED_POINT` and an invalid value stop the run at input; `POPULATIONS_RHO` without PRD lines (`PRD_N_MAX_ITER = 0`) warns and runs.

### 2026-11-09 — C_I+II_9.atom: two C II collision strengths no longer hold the temperature in place of the value (changes results) (B-42)

In `Atoms/C_I+II_9.atom` the 9000 K entries of `OMEGA 4 3` (2p ²P°1/2–²P°3/2, the [C II] 158 µm line) and `OMEGA 7 6` (2s2p² ²D) were 9.000E+03, the node temperature, and the 18000 K entry of `OMEGA 4 3` was 2.656 between 2.247 and 2.270. Through the spline in T they made both rates wrong from 1 to 19 kK, up to 4000 times too large or negative. They are now interpolated in log-log between their neighbours. The records give no source, and the published calculations disagree with each other (2D3/2–2D5/2 at 10⁴ K: 1.03, Wilson, Bell & Hudson 2005; 2.42, Tayal 2008), so the rest of the records is left as it was; the file lists the published values. Effect, FAL-C with `C_I+II_9` active at `ITER_LIMIT` 1e-4, the same binary with the old and the fixed file: upstream RH converges in 41 iterations instead of 176, and its emergent intensity changes by more than 1% at half of the 592 wavelengths, nearly all of them in the C II 133.5 nm multiplet (2s2p² ²D–2p ²P°, whose upper levels the `OMEGA 7 6` record couples), by up to 90% there; on this branch the run converges in 111 iterations where it did not converge in 500, and the same multiplet changes by up to 38%. The C II] 232.6 nm lines do not change. No test or example uses this atom.

### 2026-11-09 — Ca II collision strengths from Meléndez, Bautista & Badnell (2007) (changes results) (B-43)

The collision strengths among the five lowest Ca II levels in `CaII.atom` and its copies (`CaII_CRD`, `CaII_XRD`, `CaII_comp`, `CaII_fixed`, `CaI+II`, and the three-level `CaII_3D`) came from an unnamed source and disagreed with the R-matrix calculation of Meléndez, Bautista & Badnell (2007, A&A 469, 1203): up to 9 times lower for the 3d and 4p fine-structure transitions, half for 4s–3d, and with a 3d–3d row whose spline went negative between 50 and 100 kK. They are replaced by Meléndez et al.'s effective collision strengths (VizieR J/A+A/469/1203, 3000 to 38000 K; RH uses the 38000 K values above that). In `CaII_large`, `CaII_medium` and `CaII_med_oldbf` only the 3d–3d record changes, from a van Regemorter value of 2102, 140 times too large for this forbidden transition, to Meléndez et al.'s. Effect at `ITER_LIMIT` 1e-4, the same binary with the old and the new files: Ca II populations change by up to 30% (`1d_v10_ca_bez`, `1d_falc_h6_ca_mg`) and 22% (`1d_fala_h_ca`), the emergent intensity by up to 8%, 8% and 11% of its peak at that wavelength, and the Ca II lines by up to 2.4% (854.2 nm), 1.8% and 2.3% (866.2 nm) of their peak. Upstream RH and this branch change alike, the populations to 1e-3 and the spectra to 2e-3. Iteration counts change by at most 3.

### 2026-11-09 — Collision tables and bound-free cross-sections are interpolated monotonically in log-log; no more negative rates from tabulated data (changes results) (N-32)

Tabulated collision rates (`OMEGA`, `CE`, `CI`, `CP`, `CH0`, `CH+`, `CH`, in temperature) and explicit bound-free cross-sections (in wavelength, active atoms in `SortLambda`, passive ones in `Metal_bf`) were interpolated with natural cubic splines in linear coordinates. A spline is not local and not positivity-preserving: one steep or bad value rings through the table, and sharp turns make it negative. In the 88 atom files, 16 collision records and 225 of the 936 explicit cross-section tables went negative somewhere inside their range, giving negative collisional and photoionization rates (Ca II 3d–3d at the top of FAL-C; CI-VI_Judge's photoionization 49 → 60, which kept that case from converging) and negative bound-free opacity contributions (215 wavelength-continuum evaluations in the passive metals of the default run). Both are now interpolated by monotone piecewise cubic Hermite interpolation (`pchip.c`; Fritsch & Carlson 1980, slopes of Fritsch & Butland 1984, as in SciPy's PchipInterpolator) in log x and log y: every interval stays between its two node values, a node affects only the intervals within two nodes of it, a power law between nodes is reproduced exactly, and the result is positive. A table with a value ≤ 0 (five Mg cross-section tables contain zeros) is interpolated in y, where it stays non-negative. Outside a table the end values are used, as before. Leave-one-out tests on all atom files (each interior node predicted from the others) choose this interpolant: collision tables, median relative error 1.5e-4 against 2.2e-2 for the spline (99th percentile 4e-2 against 0.47); cross-sections, 99th percentile 5 against 35, and no negative predictions where the spline gave 407. The other interpolations of input data were checked and keep their values in range: `Linear` (fudge factors, partition functions) clamps, the molecular partition functions and equilibrium constants are fits in log space that are zero outside their range, the cubic convolution of the Barklem tables stays within 1% of the tables' range (cross-sections 126 to 3160 bohr², exponents 0.18 to 0.45), and the built-in H⁻ cross-section spline stays positive. Results change by the interpolation error of the old spline on sparse tables: the H_6 collision rates by up to 4% (14.5 kK) and Ca II's by up to 10% (28 kK) in FAL-C, converged populations in the stress set by 0.2–8% and emergent spectra by up to 4%, the regression suite by 2e-3 (Ca II K) to 2e-2 (Lyα). `CI-VI_Judge` without Ng now converges in 71 iterations (N-28; with N-31), where it never converged. New tests in `test_interp` (SciPy reference values, exact power laws, every interval bounded on random tables with spikes, locality, zeros, descending tables, ends, one- and two-point tables).

### 2026-11-09 — The stress tooling reads polarized spectra correctly (Q-21)

`tests/analyze_run.py compare` (the S-1 spectrum errors) and the tree-against-tree spectrum comparison of `tests/run_stress.sh` read every `spectrum.out` as unpolarized. A run that also writes Q, U and V (a magnetic field was read, or `BACKGROUND_POLARIZATION = TRUE`) passes the size check that way, since the four arrays together are again a multiple of Nspect, so they reported garbage (values up to ±1.4e308 and NaN for FAL-C with a 500 G field). `compare_spectra.load` now finds the layout from the vacuum-to-air record that RH writes right after I (a bool and the constant `VACUUM_TO_AIR_LIMIT`), which lies at a different place in the two layouts, and fails if neither or both fit; `--stokes` and the new `--no-stokes` state the layout and fail if the file has the other one. `run_tests.sh` keeps passing `--stokes` when `STOKES_INPUT` is set, and runs with `BACKGROUND_POLARIZATION` alone are now read correctly as well. Checked on 1-D unpolarized and polarized, 2-D and 3-D output: the right layout is found, and the wrong explicit one is refused.

### 2026-11-09 — analyze_run.py compare groups the intensities by wavelength (Q-22)

`tests/analyze_run.py compare` reshaped the emergent intensities of `spectrum.out` as (rays, Nspect), while RH stores them with the wavelength as the slowest index (`I[Nspect][rays, columns]`, `initial_xdr.c`), as `compare_spectra.py` reads them. Each of its "wavelengths" therefore mixed intensities from wavelengths spread over the whole spectrum, so `spec_wl`, `spec_line` and `worst_line` were wrong for every run. Fixed; `spec_wl` now agrees with `compare_spectra.py` (1.1e-2 for the S-1 Ng run of `1d_falc_h6_ca_mg`, at the Mg II k core). The S-1 table in SOLVER.md is corrected: Ng at `ITER_LIMIT` 1e-3 leaves spectrum errors of 5e-4 to 3e-2 in the runs it counts as converged (Lyα, Mg II k, 122.6 nm), not 1e-5 to 2e-3.

### 2026-11-09 — SOLVER = ANDERSON: stabilized Anderson acceleration of the map G, with keywords ANDERSON_DEPTH, ANDERSON_DELAY, ANDERSON_COND_MAX, ANDERSON_SAFEGUARD, ANDERSON_VARIABLES, ANDERSON_COEFFICIENTS (S-3)

Third item of the solver track (docs/SOLVER.md S-3), the replacement for Ng. `SOLVER = ANDERSON` iterates the map G of S-2 with Anderson acceleration (type II, mixing 1; Walker & Ni 2011) on the whole state: the populations of the active atoms and molecules, PRD rho, J, and with `SOLVE_NE = ITERATION` the electron density. The history slides by one every iteration over the last `ANDERSON_DEPTH` (default 5) differences, with one coefficient set for the whole state (`ANDERSON_COEFFICIENTS = GLOBAL`; `PER_BLOCK`, one set per block as Ng per atom, for comparison). The least squares is solved by QR (classical Gram–Schmidt with reorthogonalization, fixed summation order, bit-reproducible), never by normal equations; while the exact condition number of R exceeds `ANDERSON_COND_MAX` (default 1e8) the oldest columns are left out. Variables: `ANDERSON_VARIABLES = LOG` (default), where the residual is the relative change that `ITER_LIMIT` judges and populations stay positive, or `RELATIVE` (linear, scaled by fixed totals); every block of the state gets the same weight (1/√N per entry). Entries without a value in these variables take G's value. The first extrapolation comes after evaluation max(`ANDERSON_DELAY`, 2) (default 2, provisional until the delay scan of S-5). Safeguards: an extrapolated point that is not finite, negative (RELATIVE), or inverts a line or continuum where G's point does not is damped to 1/2 and 1/4 of the step, then replaced by the plain step with the history restarted; after the evaluation, a point whose weighted residual grew by more than `ANDERSON_SAFEGUARD` (default 1; 0 = off) times that of the last accepted point is rejected, the plain step from that point (already computed) is taken and the history restarts from it. One iteration is one G evaluation, rejected ones included. The run stops when the change of an evaluation (MaxChange's definition on x and G(x)) is below `ITER_LIMIT`; the output is that G(x), never an extrapolated point. The counts (extrapolated, plain, rejected, guard actions, columns dropped) are printed at the end and read by `tests/analyze_run.py`. readInput stops on `NG_ORDER > 0` with `ANDERSON`, on any `ANDERSON_*` keyword with another solver, and on the conflicts of `FIXED_POINT` (whose messages now name the solver in use). `ANDERSON_DEPTH = 0` is `FIXED_POINT` without Ng bit for bit. Default runs (`SOLVER = NG`) are unchanged bit for bit. New: `anderson.c` (generic, no RH globals), `iterateAnderson()` in `solver.c`, unit test `test_anderson` (Walker & Ni's equivalence with GMRES on a linear map to 5e-16, depth 0, collinear columns, the safeguard and restart, damping, masked entries, per-block coefficients, reproducibility) and e2e test `e2e_solver_anderson`. Measurements on the S-5 case set follow separately.

### 2026-11-09 — SOLVER = NEWTON_KRYLOV: inexact Newton-GMRES on R(x) = G(x) - x, with keywords NEWTON_NOISE, NEWTON_KRYLOV_MAX, NEWTON_ETA_MAX, NEWTON_VARIABLES (S-4)

Fourth item of the solver track (docs/SOLVER.md S-4, Form A): inexact Newton on the residual of the map G of S-2, R(u) = s(G(s⁻¹(u))) − u, in the variables of S-3 (`NEWTON_VARIABLES = LOG`, default, or `RELATIVE`) and with the same block-weighted norm. Each Newton step solves J d = −R with GMRES from d = 0 (Arnoldi by classical Gram–Schmidt with reorthogonalization in a fixed order, bit-reproducible; Givens rotations), at most `NEWTON_KRYLOV_MAX` (default 20) products and no restart; the Jacobian products are finite differences, one G evaluation each, with step h = √`NEWTON_NOISE` (default 1e-10, so h = 1e-5) in the scaled variables. The GMRES tolerance is η‖R‖ with the Eisenstat–Walker forcing term (choice 2, γ = 0.9, α = 2, their safeguard, η ≤ `NEWTON_ETA_MAX` = 0.9, and η ≥ ½ `ITER_LIMIT`/‖R‖ against oversolving). The step is globalized by backtracking on ‖R‖ (λ = 1 … 1/16, Armijo 1e-4), with Anderson's guard (finite, no new negative values, no new population inversions) before every trial; if no λ is accepted the plain step G(x) is taken. Entries without a value in the variables take G's value. The stopping test is that of the other solvers (MaxChange's definition on x and G(x) against `ITER_LIMIT`), `N_MAX_ITER` caps the G evaluations of every kind, and the output is G of the last base point. Every constant is printed at the start; every Newton step logs ‖R‖, η, the GMRES products and residual, λ and the evaluations so far; the counts (base, Krylov, line search, plain steps, guard halvings) are printed at the end and read by `tests/analyze_run.py`. readInput stops on `NG_ORDER > 0`, on `ANDERSON_*` with `NEWTON_KRYLOV`, on `NEWTON_*` with another solver, on values out of range, and on the conflicts of `FIXED_POINT`. The keyword setter for the variables is shared (`setSolverVariables`), and `iterateAnderson()` now shares the block and variable set-up and the per-atom change print with the new solver, with unchanged output. Default runs (`SOLVER = NG`) are unchanged bit for bit. New: `newton.c` (generic, no RH globals), `iterateNewton()` in `solver.c`, unit test `test_newton` (GMRES against a dense solve, the minimal-residual property and true residual for k < n, the tolerance stop, masks, an invariant space, reproducibility, the forcing term) and e2e test `e2e_solver_newton`. First run (Ca II with two PRD lines on FAL-C, H_6 passive, `ITER_LIMIT` 1e-6): converged in 124 G evaluations (18 Newton steps, 100 Krylov products, 23 line-search trials), against 89 for `ANDERSON`; h between 1e-6 and 1e-4 gives 119–124, larger h more (182 at 1e-3, 381 at 1e-2). Fixed during testing: the guard judges a Newton trial against the base point as well as G(x_k) (otherwise an inversion of x_k blocked every λ), and a trial whose log variables would overflow exp() is a guard failure instead of a SIGFPE; every Newton step logs |d|_max and the guard failures by reason. First results on the S-5 set: 1.4–3.2× the G evaluations of `ANDERSON` on the ten cases that converge, `1d_shock_h6` not converged after 3000 (docs/SOLVER.md S-4, Status). After review: the Armijo test does not accept a trial whose image loses one of the base point's entries (`1d_shock_h6crd_ne` 77 → 88 evaluations, the other cases unchanged), and when `N_MAX_ITER` ends a failed line search one more evaluation at the base point keeps G's side effects consistent with the output.

### 2026-11-09 — SOLVER = ANDERSON: ANDERSON_DEPTH defaults to 10 (S-5)

The S-5 solver comparison (docs/SOLVER.md S-5, MODERNIZATION.md §2.5) found that depth 10 needs less work than depth 5 in 9 of its 11 one-dimensional cases, most on the hard ones (`1d_t5750_h6` 300 against 432 formal-solution equivalents, `1d_shock_h6` 320 against 580), and more on one (`1d_falc_civi`, 151 against 133). Both met every S-5 criterion. The default of `ANDERSON_DEPTH` is therefore 10 instead of 5. The history then takes 40·N doubles instead of 25·N for a state of N unknowns; set `ANDERSON_DEPTH = 5` where memory is tight. Runs with `SOLVER = ANDERSON` that do not set `ANDERSON_DEPTH` change; the default solver (`NG`) is unaffected.

### 2026-11-08 — The H⁻ density comes from neutral hydrogen, as in the chemical equilibrium (changes results) (B-40)

ChemicalEquilibrium stored the H⁻ density (atmos.nHmin, used for the H⁻ bound-free opacity and written to molecules.out) as ne * n_H * Phi_H- with n_H all atomic hydrogen, while the hydrogen conservation equation it had just solved uses the neutral part, ne * fn0 * n_H * Phi_H-. The stored value is now the one in the equation, fHmin * n_H. Where hydrogen is ionized the old value was too large by 1/fn0, up to 10^11 at the top of FAL-C, so grazing rays, which cross the chromosphere, saw spurious H⁻ opacity and emission. Results change in every regression case, at the inclined rays only. In the 5-ray 1-D cases the rays near disk centre change by at most 7e-4. At the most inclined ray (mu ≈ 0.05) the continua below 200 nm change by up to 9% (FALC_82_s, 30 nm; the same with ITER_LIMIT = 1e-8, so this is not convergence noise) and 15% (1d_falc_h6_ca_mg, 22.8 nm). In 1d_t5750_h6 the two most inclined rays change by a median of 1%. In the spherical example the limb rays change by up to 17%. Every case still matches its own 1-thread result at 2, 3 and 8 threads. Correction to the log: B-40 said the emergent spectra change by at most 5e-4 (median 3e-5); that holds for the rays near disk centre only.

### 2026-11-08 — Rate equations solved by GTH elimination; He_large converges (changes results) (N-31)

The rate equations of every atom and molecule were solved by LU decomposition, which is normwise stable but loses the populations of levels tightly coupled by collisions: in He_large the fine-structure multiplets (4F, 5D, 5F, split by 1e-4 cm⁻¹) were wrong by up to 4.4e-3, and its iteration stalled at that level (N-28). They are now solved by GTH elimination (`gth.c`; Grassmann, Taksar & Heyman 1985), which works on the off-diagonal rates only, never subtracts, and gives every population to rounding at the same cost. The reference level (`I_SUM`, or the most populated level) is eliminated last. Where a rate is negative, not finite, or a level has no rates to the others, that location is solved with LU as before; the log names the first such rate of each call (with whether it is collisional or radiative), and the work summary counts the rate equations solved and how many needed LU. He_large: the noise of the iteration map drops from 1.7e-3 to 3.9e-14, and the reference recipe (no Ng, `ITER_LIMIT` 1e-10) converges in 667 iterations where it stalled at 2.3e-3; with Ng the outcome still depends on Ng's extrapolation (2 of 6 variants converge, LU 3 of 6), but no longer on the thread count. The noise of the other cases drops by a factor 100 to 3000 (Ca II 1.3e-10 to 5.6e-13). The regression suite changes by at most 9e-12. In the convergence stress set (`ITER_LIMIT` 1e-3) the populations of the converged cases change by 2e-11 to 2e-8, 1.6e-6 for Fe33 and 4.5e-7 for PRD Ca II, with the same iteration counts except He_large. Negative rates met so far: from atomic data, Ca II 3d–3d collisions at 60–71 kK (B-43, N-32) and a CI–VI photoionization cross-section that its spline makes negative; and radiative rates out of the He I ground level after an Ng extrapolation of He_large. None in a converging iteration without Ng. New unit test `test_gth` (detailed-balance rate matrices with known populations, including the multiplet case where LU is off by 6e-5). Answers SOLVER.md §9 Q12. Same code upstream.

### 2026-11-08 — The map G and SOLVER = FIXED_POINT, with ne in the state; gcheck programs; keyword BACKGROUND_DOUBLE (S-2)

Second item of the solver track (docs/SOLVER.md S-2 and S-2b). One MALI iteration is now available as a function G of an explicit state vector (populations of active atoms and molecules, PRD rho, J, and with SOLVE_NE = ITERATION the electron density), with no acceleration and a fixed number of PRD sub-iterations inside. The new keyword SOLVER selects NG (the default, unchanged bit for bit) or FIXED_POINT, which iterates x <- G(x), with RH's Ng on G's output when NG_ORDER > 0. Keyword settings that G cannot honour (PRD early exit or PRD Ng, collisional relaxation, hydrostatic equilibrium, LIMIT_MEMORY, and with ne in the state OLD_BACKGROUND and NON_ICE) stop the run when the keywords are read, with a message for each. With ne in the state the background is recomputed from the state at the start of every G: ntotal from its initial value, LTE quantities and collision rates, a chemical equilibrium with the neutral fractions of active atoms taken from the state and a fixed number of Newton steps, the opacities; after the formal solution and statistical equilibrium, ne follows from charge conservation in a fixed number of Newton steps from the state's ne. The new keyword BACKGROUND_DOUBLE stores the background file in double; with ne in the state the background is always stored in double, and BACKGROUND_DOUBLE = FALSE stops the run. New programs gcheck_1d, gcheck_2d and gcheck_3d check G on states a run reaches: G twice on the same state, G(x_a), G(x_b), G(x_a) for hidden state, the noise of G (ECnoise) and finite-difference Jacobian products, per block of the state; new ctest e2e_solver_g. Results: FIXED_POINT without Ng reproduces the default solver without Ng (and with PRD_ITER_LIMIT = 0) bit for bit; G is bit-reproducible and free of hidden state on all S-5 cases, in 2-D and 3-D at 4 threads, and with ne in the state (where repeated calls differed by up to 55 % before); 2000 calls keep the descriptors and memory flat; plain G converges 1d_falc_solvene in 859 iterations to a fixed point within 4e-8 of the current code's. The noise of G is 1e-13 to 2e-10 on converged states, except He_large (2e-3). With ne in the state, the float background limits finite-difference Jacobian products to ~3e-3 (double: ~1e-6), which is why double is required there.

### 2026-11-07 — No memory leak per Background() call; 3-D Stokes runs allocate the intensity once (B-11)

Background() (and readBRS() with OLD_BACKGROUND) allocated backgrflags and backgrrecno again on every call without freeing the previous arrays. With SOLVE_NE = ITERATION that is one pair per iteration: plain iteration of 1d_falc_solvene grew from 19.2 MB peak at 300 iterations to 28.4 MB at 1500 (7.6 KB per call); the previous arrays are now freed first, and the same runs peak at 16.5 MB both. In 3-D with a magnetic field (or polarized background) spectrum.I was allocated twice; it is now allocated once, as in 2-D (the first allocation was never written, so it cost address space rather than resident memory). Results are byte-identical: regression suite in exact mode, the 1500-iteration populations, and a 3-D FULL_STOKES run of FALC_04x04x82 with a uniform field.

### 2026-11-07 — Active molecules: Gamma starts from zero every iteration, so the rates of overtone lines no longer accumulate (B-38)

initGammaMolecule set only the Delta v = 1 entries of Gamma (the collisional rates), while the radiative rates of every vibration-rotation line are added with +=. The entries of lines with Delta v >= 2 (overtone bands) were never reset, so their radiative rates grew by one iteration's worth every iteration. Gamma now starts from zero, as initGammaAtom does for atoms. The shipped CO line lists have Delta v = 1 only and give byte-identical results (active CO on FALXCO_80 with vmax=3_Jmax=49_dv=1_26: spectrum, molecules.out and pops_mol.CO.out). With five synthetic 2-0 lines added, the old code never converges: the population change stays at 3.7e-5 per iteration for 40 iterations, drifting with the iteration count; with the fix it falls to 9e-15 in 6 iterations. No test runs an active molecule; the reproducer is in the pull request.

### 2026-11-07 — Background() can be called many times: the background file is closed before it is reopened, and the Kurucz lines are read once (B-39)

With SOLVE_NE = ITERATION (and in hydrostatic equilibrium) Background() runs once per iteration. Every call opened background.dat again without closing the previous descriptor, so a run ended when it reached the descriptor limit: with the common limit of 1024, the S-1 reference runs of 1d_falc_solvene and 1d_shock_h6crd_ne stopped at iteration ~1010 with "Unable to open output file brs.out"; under ulimit -n 64, 1d_falc_solvene stopped at iteration 52. The previous descriptor is now closed first, and the same run completes 100 iterations under that limit. Every call also read the Kurucz line lists and the three Barklem tables again. The lists do not depend on the atmosphere, so they are now read once and kept when Background() will be called again, like the passive atoms and the partition functions: 1d_falc_solvene with gf0400.10 (5235 lines, 27 calls) reads them once instead of 27 times, 8.0 s → 7.2 s. Results are byte-identical (regression suite in exact mode; the Kurucz run's spectrum and populations). Correction to the log: the line list was not leaked, Background() freed it at the end of every call; the memory that grows with the number of calls is backgrflags and backgrrecno (B-11, ~9.5 KB per call in 1d_falc_solvene), which this change leaves alone.

### 2026-11-07 — Work counters, SE residual norms and reference solutions for the solver comparison (S-1)

First item of the solver track (docs/SOLVER.md). Every run now ends with a "Work" block in the log: the wavelength solves per phase (initial background-scattering iteration, main iteration, final solution) and kind of sweep (with approximate operator, plain, PRD sub-iteration over the PRD wavelengths only), and the calls and wall time of Background() and Hydrostatic(). Results are byte-identical to before (run_tests.sh --exact). tests/se_residual.py reports the candidate norms of the SE residual (max relative imbalance over all levels and over levels above 1e-10, 1e-6 and 1e-3 of the total, RMS in units of the total, conservation); analyze_run.py adds the work counts and the residual to its summary and compares a run with a reference (populations with the same thresholds, spectra per wavelength and per line). run_stress.sh gains --kw and --ref, the CRD + ne case 1d_shock_h6crd_ne and the Bifrost PRD limit-cycle case 1d_bf_c100_prd; tests/make_references.sh runs the reference recipe (no Ng, fixed PRD sub-iterations, ITER_LIMIT = 1e-10) and keeps the outputs with a manifest. References exist for seven 1-D cases plus 1d_falc_solvene, which reaches 1e-10 only with the background file in double (2.3e-9 floor in float). Against them, Ng runs that converge to 1e-3 are 2e-4 to 4e-2 off in the populations of levels above 1e-6 of the total and 1e-5 to 2e-3 in the spectrum (N-3); the residual alone does not bound the error. Plain iteration with fixed PRD sub-iterations converges 1d_t5750_h6 (2527 iterations) and the Bifrost limit-cycle column (645), which Ng does not.

### 2026-11-06 — Per-transition loops vectorized, -O3 by default; rowcol covers all periodic-table rows; chartflow uses the active opacity (P-7, B-36, B-37)

Step 3 of the float plan (§2.4). The loops over depth of the active opacity and the radiative rates (Opacity, addtoRates, addPreconditionedRates; 30–36 % of NLTE runs) did not vectorize: their arrays may alias as far as the compiler can tell, and GCC before 12 does not vectorize at -O2. They now have restrict pointers and are marked RH_SIMD (simd.h, #pragma omp simd with -fopenmp-simd where the compiler takes it), and RelWithDebInfo, the default build type, uses -O3 (option RH_O3). Whole runs are 1.15–1.23× faster with bit-identical results (FAL-C H_6 + Ca II + Mg II PRD 1.16×, moving Mg II hybrid PRD 1.23×, 2-D 32×82 1.15×; one pinned CPU, Xeon E5-2697 v2, gcc 8.5). The loops alone give 1.08–1.11×, -O3 alone 1.07–1.14×. -O3 exposed B-36: rowcol (collision keyword SUMMERS) left the periodic-table row unset for 55 to 86 electrons; it now covers all six rows and stops with an error beyond them. The CI build with GCC 15 at -O3 also exposed B-37: chartflow, at wavelengths without angle-dependent opacity, added two never-written arrays instead of the active opacity and emissivity.

### 2026-11-06 — PRD scattering integrals with the interpolation folded into the weights; angle-dependent PRD threaded (P-8, P-14)

Step 2 of the float plan (§2.4). The scattering integrals of PRD lines interpolated J (Jgas, or the intensity of every absorption direction) from the line's grid onto a fine absorption grid on every call; interpolation and Hunt were a third of a 1-D PRD run. The integral is linear in the values on the line's grid, so the interpolation is now folded into the weights when they are computed, and every later call is a short dot product per emission wavelength (prdfold.c). The weights stay in the scratch file per line, in single precision (the dot product is in double), and are about 10× smaller than before (Ca II: 2.5 MB → 0.25 MB per line for FAL-C angle-averaged, 253 MB → 24.5 MB angle-dependent). Angle-dependent PRD (PRDAngleScatter) now reads the intensities of all absorption directions once per call instead of once per emission wavelength and direction, and is threaded over locations. On 1-D runs to convergence the integrals are 4× faster (H_6 + Ca II + Mg II angle-averaged 14.6 → 3.5 s, Mg II hybrid 7.5 → 1.9 s, Bifrost column 1.3 → 0.3 s), whole runs 1.3–1.4× (cross-redistribution 2.2–2.7×), and angle-dependent Ca II 316 → 75 s on one thread and 10.6 s on eight. Iteration counts are unchanged and populations agree with the previous code to ≤ 3e-8 (the 2-D hybrid case converges in 38 instead of 44 iterations, to within 4.6e-6). The spline representations, never used, are gone; RII (1-D) sets up its quadrature with pthread_once. New unit test test_prdfold: folded integral = interpolated-then-summed to 8e-16.

### 2026-11-06 — Line profiles and the background opacity file in single precision (M-1, M-7)

Step 1 of the float plan (§2.4). The line profiles of active atoms and molecules (phi and, in full Stokes, phi_Q/U/V and psi_Q/U/V) are stored as float, in memory and in the LIMIT_MEMORY files. Each value is still summed in double and rounded once, and wphi normalizes the stored values. The background opacity file is written in float. Its emissivities are scaled by 2^100: SI emissivities reach 1e-50 in FAL-C below 100 nm, and without the scale a C I-VI run lost them and trapped on J = 0. The format is recorded in two trailing ints of brs.out (4 bytes, exponent 100); without them, as in files from earlier versions, the file is read as unscaled double. rhanalyze and the IDL readers handle both formats. The 2-D 128×82 perf case (H_6, Ca II, Mg II, Fe33, 5 iterations) peaks at 10.7 GB instead of 17.5 GB (−39 %), and its background.dat is 1.03 GB instead of 2.05 GB. Results change at float-rounding level: the regression suite differs from modernize by at most 7e-8, and converged FAL-C NLTE cases (H_20, H_6 + Ca II + Mg II, Fe33, moving Ca II) take the same number of iterations, with populations within 8e-7 and spectra within 7e-7. Sunspot pixels with a polarized magneto-optical background (PHOMI Fe I 617.3) agree to 1.2e-7 of I_max in I, Q, U and V. He_large stays chaotic under Ng (N-2, N-3): float overflows at more thread counts than double, which also overflows at 3 threads. New unit test test_background. M-1 and M-7 stay open: the profiles are still stored per ray, and the background is still a file.

### 2026-11-06 — rhanalyze reads the background scattering opacity correctly at polarized wavelengths (B-35)

At a polarized wavelength a block of `background.dat` holds the opacity (4 records), the magneto-optical terms (3, if enabled), the emissivity (4) and the scattering opacity (1). `rhopacity.py` skipped only 3 records past the start of the emissivity to find the scattering opacity, so it returned the Stokes-V emissivity instead; it now skips 4. Checked on the 1-D Ca II example in full Stokes with the Kurucz list `gf0400.10` (94 polarized wavelengths, magneto-optical on): the returned `scatt` now equals the last record of each block exactly and agrees with the neighbouring unpolarized wavelengths to 5 %; before, it was negative at 91 of the 94. Readers in C and IDL were already correct.

### 2026-11-06 — No SIGFPE when a run uses less than one clock tick (B-33)

`printTotalCPU()` reported user and system time as percentages of their sum, from `times()`, which counts in clock ticks of 10 ms. A run that finishes in less than one tick of CPU time has a sum of 0, and the 0/0 raised FE_INVALID, which `SetFPEtraps()` turns into SIGFPE at the very end of the run. This was the intermittent failure of `e2e_stokes_1d` (run `rot_40_turned`, LTE with `N_MAX_ITER = 0`): it only happened on the fastest CI runners (AMD EPYC 9V45, Zen 5; 3 of 4 jobs on ubuntu-26.04, none of 14 on other CPUs), whenever the run happened to stay under 10 ms. The CPU and backtrace diagnostics added for B-33 on 2026-10-22 caught it on #222 and on `modernize` after #223. The percentages are now printed only when the total is positive, otherwise `(less than one clock tick)`. Checked on sundog2 (gcc 8.5) by preloading a `times()` that reports zero ticks: before, SIGFPE with `addr2line` pointing at `getcpu.c:101`; after, the run completes (3 of 3) and normal runs still print the percentages. The same code is in upstream RH.

### 2026-10-22 — Diagnostics for the intermittent SIGFPE in `e2e_stokes_1d` (B-33)

B-33 hit again on #211 (`build (ubuntu-24.04, gcc)`, same run `rot_40_turned`, status -8); the re-run passed. It does not reproduce locally on the #211 branch: gcc 8.5 30/30 and gcc 13.4 20/20 pass, `MALLOC_PERTURB_` with seven fill bytes passes, valgrind `--track-origins` reports 0 errors. Two compilers failing intermittently on the same deterministic single-threaded run points at the runner (CPU-dependent libm code paths in glibc 2.39) rather than the code, but that is unconfirmed. The e2e scripts now print the last 40 lines of `out.log` when a run fails, which includes the trapped exception and backtrace from `fpehandler.c`, and the CI jobs print the runner CPU model, its FMA/AVX2/AVX-512 flags and the libc version. Diagnostics only; B-33 stays open.

### 2026-10-21 — Converged populations satisfy statistical equilibrium when transitions overlap (changes results) (B-6)

The preconditioned rate equations subtracted Psi times the emissivity of all active transitions from the intensity (`Ieff`), but Gamma put only part of it back with the new populations. Pairs of levels without a transition were missing, one of the two opacity terms was missing on connected pairs, and the term for a lower level that is also an upper level was added once for every transition sharing it (the multi-ion double counting first suspected in B-6). Only lines sharing an upper level were exact. The iteration then converged to populations that do not satisfy the plain rate equations with their own radiation field. The residual of those equations after convergence to 1e-10 (FAL-C, CRD) was 1e-7 for Ca II (H and K wings), 7e-4 for H_6 (1-D, 2-D, 3-D, spherical) and 3e-2 for MgI+II; now it is 1e-11 to 5e-10 in all geometries. MgI+II: the Mg II 4s and 3d populations move by 2-4 % (median over depth), up to 22 % (34 % with PRD h and k); Mg III ground up to 15 %; h and k upper levels 1e-3. Lightweaver's form has no double count but the same fixed-point error. Full preconditioning (Rybicki & Hummer 1992, eq. 2.19) is exact but puts negative rates between unconnected levels and diverges on He_large and CI-VI. Each transition now gets its own `Ieff`, without the emission of its two levels, and Gamma puts exactly that emission back (`addPreconditionedRates` in `fillgamma.c`, also used by `rhsphere`; the unused `chi_up`/`chi_down` arrays are gone). The fixed point is exact by construction and Gamma has rates only between connected levels. The cost per iteration is unchanged (MgI+II, 8 iterations: 40.2 vs 39.5 s). Iterations: Ca II 42, H_6 67, MgI+II 108 vs 102; stress set unchanged except 1d_falxco_h6_ca 109 → 136, 1d_falc_h6_ca_mg 35 → 31, He_large 75 → 81, 2-D perturbed 76 → 74. He_large and CI-VI do not converge tightly with any of the forms (N-28). New: `tests/unit/test_coupling` (exact fixed point and linear form for random overlapping active sets; the old form is off by 4-33 %) and `tests/se_residual.py` (the residual of statistical equilibrium of a finished run, any geometry). Derivation in `docs/unit-tests`. Also present upstream.

### 2026-10-21 — CaII_large.atom without the Inf collision records (B-15)

The 13 `OMEGA` records of `Atoms/CaII_large.atom` that were `Inf` are removed, as was done by hand in `CaII_medium.atom`. All collision data in the file are van Regemorter estimates from `tools/impact.c`, which divides by the energy difference of the two levels; the model gives the two fine-structure levels of the nF terms and of the nD terms with n >= 11 the same energy (NIST lists 4f 7/2 and 5/2 at the same energy too), so these records were infinite. The collision rates between these levels are small compared with the radiative rates in and out of them. The model's other collision strengths between fine-structure levels are far too large as well (3d 3/2-5/2: 2100, against 15 from Meléndez et al. 2007), which is left as it is. `1d_falc_cail` now converges (23 iterations).

### 2026-10-21 — solve1d validates its boundary-condition codes (B-18)

`solve1d` read the top and bottom boundary-condition codes from standard input and stored them without checking them. It always runs the Bezier solver, where a bad top code already stopped with an error, but a bad bottom code fell through the switch and left the incoming boundary intensity uninitialized. Both codes are now checked against `enum boundcond` when they are read, and a bad one stops the run with an error naming it.

### 2026-10-21 — StopRequested no longer leaks a file handle (B-19)

StopRequested() tested for the STOP_RH file with fopen() and never closed the returned FILE, leaking a handle on every call while the file was present. It now uses access(F_OK), a pure existence test; the warning message and return value are unchanged.

### 2026-10-21 — BiLinear removed (no callers since the Voigt consolidation) (B-21)

`BiLinear()` in `linear.c` read one element past its tables at the upper edge, extrapolated asymmetrically and kept its hunt indices in function statics. Its only production caller was the Voigt LOOKUP table, removed in PR #146, so the function and its prototype are deleted rather than fixed, together with its part of `test_interp`.

### 2026-10-21 — w3js checks its arguments; 3j/6j overflow stops with an error (B-22)

`w3js()` accepted arguments that violate the selection rules for j - m and j1 + j2 + j3 being whole (`w3js(0.5, 0.5, 0.5, 0.5, -0.5, 0)` returned 1); it now returns 0 for them. `fact()` returned 1 for a non-integer argument and read past its table above 300, and every entry above 170! in it was infinite; it now holds 0! to 170! and stops with an error outside that range. The products of factorials in the Racah formulae overflow from j = 49, which gave +-inf (j = 50 to 56) or 0 (from j = 58); `w3js()` and `w6js()` now stop with an error instead. The only caller, `paschen.c`, uses small J, where nothing changes.

### 2026-10-21 — Hunt agrees with Locate on exact node values (B-23)

Hunt returned one index too low when the value equalled a table node and the search started on the wrong side of it (hunting downward in ascending tables, upward in descending ones), so it disagreed with Locate. The hunt loops now use the same "last node reached" condition as the bisection. Linear and spline interpolation are unchanged at an exact node, and regression output is byte-identical. The unit test no longer tolerates the deviation.

### 2026-10-21 — Gauss-Legendre weights accurate to rounding (N-20)

`GaussLeg()` computed each weight from P_n' at the previous Newton iterate instead of at the converged node, so the weights were off by up to 2e-12 (n = 64) and summed short of the interval length by up to 2.4e-14. The derivative is now evaluated at the node itself (the recurrence moved into a helper, `legendre()`). Against mpmath the weights are now accurate to 2e-15 at n = 8, 1e-14 at n = 24 to 32 and 1.2e-13 at n = 100 (the outermost node), the limit set by rounding the node to double; they sum to the length to 1e-15, and polynomial integrals are exact to 2.6e-14 instead of 5e-13. Angle quadratures (1-D rays, 2-D/3-D inclinations, angle-dependent PRD) change accordingly; regression spectra change by up to 4e-12.

### 2026-10-21 — Planck with expm1; E1 and E2 beyond x = 80, E2 without cancellation (N-22)

`Planck()` evaluated exp(x) - 1, whose relative error grows as eps/x at small x = h nu/kT (6e-10 at x = 1e-7, radio wavelengths); it now uses `expm1`, and the Rayleigh-Jeans limit is exact to rounding. `E1()` and `E2()` returned 0 above x = 80; `E1` now uses A&S 5.1.56 for all x >= 1, and both are zero only where exp(-x)/x underflows (x around 740) and for x = inf (`fixedrate.c` with T = 0). `E2` came from the recurrence exp(-x) - x E1, which cancels as x grows (1.2e-6 at x = 80); for x >= 1 it now comes from the continued fraction 5.1.22 (modified Lentz), accurate to 6e-15. In RH itself only `E1` is used (collision rates, fixed rates); with dE/kT > 80 a collisional rate that was zero is now computed. `E2` is used by `tools/make_h.c` only.

### 2026-10-20 — PRD no longer diverges on moving realistic atmospheres: the emission profile ratio stays positive (changes results) (B-29)

The ratio of emission to absorption profile in a PRD line was computed as rho = 1 + gamma (<J>_R - Jbar), gamma = n_i B_ij / (n_j P_j). That form uses the statistical-equilibrium identity n_j P_j = n_i B_ij Jbar + Y_j (Y_j: all other inflow into the upper level, plus n_j Q_elast). During the PRD sub-iterations the populations lag the radiation field, the identity does not hold, and rho goes negative once gamma Jbar > 1. The negative emissivity then gives negative J, negative populations, a negative opacity and a floating-point exception. This is what killed Ca II PRD on the Bifrost cuts. rho is now computed without the identity, as (1 + sum_x g_x <J>_R(x)) / (1 + sum_x g_x Jbar(x)) with g_x = n_i(x) B_ij(x) / Y_j, summed over the subordinate lines x with cross redistribution. This form is positive for any non-negative populations and equal to the old one in statistical equilibrium, so converged results are unchanged: on FAL-C, Ca II reaches the same populations to 6e-8 with angle-averaged, hybrid and angle-dependent PRD and with XRD. Angle-averaged and hybrid PRD take 37 iterations instead of 46; the other modes take about as many as before (40 -> 39, XRD 42 -> 43). Unconverged runs change: the 4-iteration regression cases differ by up to 16% in the PRD line cores. In the regression cases twoactive (Lyalpha PRD) and spherical, the old form stalls at a population change of 1e-5; the new one converges to 1e-8, and the two agree to 2e-5. On the Bifrost slice, hybrid PRD now converges in all 126 1-D columns (120 before) and in 2-D (126 x 89, 38 iterations, 4 min on 16 threads; before: floating-point exception at iteration 11). With angle-averaged PRD, 32 of the 126 columns end in a limit cycle at a change of 3e-3 to 1e-2. That is the early exit of the PRD sub-iterations (N-4): with `PRD_N_MAX_ITER = 10`, `PRD_ITER_LIMIT = 1e-4`, 8 of 8 tested columns converge, against 6 of 8 with the old rho (the other two crash). One column dies from an Ng step that overshoots (N-2). Also present upstream.

### 2026-10-20 — Hybrid PRD with a single ray uses the right gas-frame mean intensity (changes results) (B-31)

With `PRD_ANGLE_APPROX` and one ray (`NRAYS = 1` in rhf1d, or an angle set with one ray), the gas-frame mean intensity `Jgas` was wrong in two ways. It was zeroed before each formal solution only for more than one ray, and the formal solvers added to it only for more than one ray, except in the branch for wavelengths without moving or polarized line opacity, which added regardless. So at those wavelengths `Jgas` grew by J with every formal solution (on FAL-C, 27 to 34 times J after 30 iterations), and at the wavelengths of the PRD lines it stayed zero. The populations then differed by 25% from angle-averaged PRD, which must give the same answer in a static atmosphere. The test on the number of rays, meant to leave `Jgas` alone in the final single-ray solutions (`solveray`), is removed; those solutions do not update J, so `Jgas` is not touched there anyway. With one ray, `Jgas` now equals J in a static atmosphere and hybrid PRD gives the same populations as angle-averaged PRD, byte for byte; in a Bifrost column `Jgas` is filled at every wavelength. In 2-D and 3-D the wavelengths without PRD line opacity added each ray direction with the full ray weight instead of half of it, so `Jgas` was 2 J there; those rows are not used by the PRD scattering integral, so only `Jgas.dat` was affected, but they are now J as well. The 2-D Eddington solver now also fills `Jgas`. The e2e test `tests/e2e/prd_hybrid.py` checks `Jgas` = J and hybrid J = angle-averaged J in rhf1d (one and three rays) and rhsc2d (one and six rays) on static atmospheres. Runs with more than one ray are otherwise unchanged (regression suite byte-identical). Also present upstream.

### 2026-10-20 — Hybrid PRD on large grids no longer overflows its interpolation table (B-32)

The table of gas-frame interpolation coefficients of hybrid PRD (`PRD_ANGLE_APPROX`, built in `initSolution` unless `PRD_LIMIT_MEM` is set) has one entry per wavelength, ray direction and grid point, 2 x Nrays x Nspect x Nspace. Its size and index were computed in `int` and its cumulative counts stored as `int`, so beyond 2^31 entries the index wrapped and the run crashed. On a 3-D cut of the Bifrost snapshot (42 x 42 x 89, every 12th point) with Ca II, 286 wavelengths and 24 rays, 2.155e9 entries, rhsc3d died with a segmentation fault in `initSolution` after 3 minutes. The index and counts are now `long`, and the allocations are checked. The same run now builds the table and iterates (115 GB, 1 h 38 min for one iteration, most of it spent building the table on one thread). `matrix_int` and `matrix_char` computed their size the same way and now compute it in `size_t` (`matrix_double` already did); all three now stop with an error message when an allocation fails. Angle-dependent PRD with more than 2^31 values of rho per line, which the Ng acceleration of the PRD iteration cannot index, now stops with a message instead of overflowing. Results are unchanged (regression suite byte-identical). Also present upstream.

### 2026-10-20 — 2-D and 3-D Bezier solvers use the quadratic Bezier slope at the upwind point (changes results) (N-24)

The cubic Bezier solvers of rhsc2d and rhsc3d (unpolarized, and the polarized DELO-Bezier solvers) took the secant (S_c - S_u)/dtau as the derivative at the upwind point. That curve reproduces only linear source functions, so the solvers were second order where the 1-D solver, which carries a centred derivative to the upwind point, is fourth. The upwind slope is now 2 x secant - S'_c, the cubic form of the quadratic Bezier curve through the two points with the local derivative (one control point S_c - dtau/2 S'_c). It is exact for quadratic functions, costs nothing, and keeps the curve monotonic: the Steffen limiter keeps S'_c between 0 and twice the secant, so the upwind slope lies there too. The same is done for the opacity in the optical depth and, polarized, for the source vector and the propagation matrix. On the analytic slab of the unit tests ("layers"), the emergent error at 161 points drops from 4.8e-4 to 2.5e-5 and the order rises from 1.99 to 2.82 (2-D and 3-D; polarized 2-D 5.8e-4 to 3.3e-5, order 1.96 to 2.80); on "columns" the vertical interpolation of I limits it to 2.6. Interpolating derivatives to the upwind point instead (option A, measured on #150) reaches fourth order on "layers" but not on "columns", at 15% more time per iteration. On a perturbed FAL-C 2-D atmosphere refined from 82 to 325 depth points (Ca II, converged to 1e-5), the emergent continuum is 3 times closer to the refined solution (median 3.0e-3 vs 7.2e-3 at 82 points, 4.0e-4 vs 1.4e-3 at 163); in the line cores the sampling of the transition region dominates both solvers. Results of 2-D and 3-D runs change by up to the solvers' own discretization error, which on the Bifrost cuts and the perturbed atmosphere is 10-40% in the upper-chromosphere populations and up to 20% in individual line-core intensities at the resolution used.

### 2026-10-20 — The derivative limiter of the polarized Bezier solvers no longer depends on the reference direction of Q and U (changes results) (N-25)

The DELO-Bezier solvers (1-D, 2-D, 3-D) limit the derivatives of the source vector and of the propagation matrix with Steffen's limiter, which is nonlinear: zero at an extremum, otherwise the centred slope clipped to twice the smaller one-sided slope. Applied to Q and U, or to eta_Q and eta_U, one at a time, it does not commute with turning the reference frame, which mixes them: the result depended on the arbitrary reference direction of linear polarization. The pairs that turn with the frame, (S_Q, S_U), (eta_Q, eta_U) and (rho_Q, rho_U), are now limited as vectors: the derivative is zero where the one-sided slopes point in opposite directions (negative dot product), and otherwise the centred slope shortened to at most twice the shorter one-sided slope. This uses only lengths and a dot product, so it turns with the frame, and it is Steffen's limiter for a single component. For a vertical ray on FAL-C with a depth-dependent field, turning the field azimuth by 35 degrees now turns (Q, U) to 3e-14 (before: off by 7e-3 of their maximum, V by 2e-4, I by 1e-5); on a grey atmosphere to 1e-14 at every depth resolution (before 4.8e-3 at 40 points to 5.0e-6 at 320). Results change only where the direction of (Q, U) in the frame of the ray varies along the ray: for a vertical ray where the field azimuth varies, for an inclined ray also where the inclination varies. Where it is fixed the vector limiter gives the same derivatives as before (the regression cases, with a fixed field direction, agree to 1e-17). In 2-D and 3-D the solvers remain frame dependent in another place: they interpolate the upwind and downwind intensity, source vector and propagation matrix component by component with a monotonic (clipping) interpolation. The e2e test now requires exact rotation on FAL-C and the grey atmosphere, and a unit test checks that the limiter commutes with rotations on random data and reduces to Steffen's for one component.

### 2026-10-20 — Hybrid PRD on several threads no longer loses updates of the gas-frame mean intensity (T-4)

With `PRD_ANGLE_APPROX` the formal solution of each wavelength adds its intensity, shifted to the gas frame, to the gas-frame mean intensity `Jgas` of neighbouring wavelengths. Neighbouring wavelengths run on different threads at the same time, and the additions had no lock: updates were lost (`Jgas` off by up to 34% after one iteration with 4 threads), so results and convergence depended on thread timing. Each row of `Jgas` now has a lock. A ray's contributions are summed in a small buffer of the thread first and added row by row, so there are a few lock operations per ray: the lock costs about 1% of time and memory at 16 threads (2-D 128 x 82: 100 s vs 99 s; 3-D 24 x 24 x 82: 937 s vs 939 s, 29.25 vs 29.00 GB). Locking every single addition cost 5-25%. On a perturbed FAL-C 2-D case with hybrid Ca II PRD, 8-thread runs took 22 and 140 iterations and their populations differed from the 1-thread run by up to 0.9%; now all take 18 iterations and agree to 3e-10. On one thread the results are unchanged. Also present upstream.

### 2026-10-19 — OH and CH bound-free opacities no longer start Hunt() from uninitialised indices (B-30)

`OH_bf_opac()` and `CH_bf_opac()` start `Hunt()` from index 0. Before, the start indices were uninitialised: valgrind reported conditional jumps on and uses of uninitialised values in every run with molecules, and uninitialised bytes written to the background scratch file. `Hunt()` brackets correctly from any start, so results do not change (byte-identical spectrum).

### 2026-10-19 — End-to-end Stokes tests of rhf1d with a fictional atom (Q-14)

`tests/e2e/stokes_1d.py` (ctest label `e2e`, Python standard library only, about 10 s) runs `rhf1d` with a fictional two-line atom (`tests/e2e/TI_fict.atom`, anomalous Zeeman patterns with g_eff 4/3 and 7/6) in LTE, one formal solution with J = 0 in the background scattering, so every run solves the same linear problem. Checks: no field gives Q = U = V = 0 exactly; parity (mirror in the plane of the rays: gamma -> pi - gamma, chi -> pi - chi) gives (I, Q, -U, -V) to 3e-13 with and without magneto-optical effects; without them B -> -B gives (I, Q, U, -V) to 2e-13 (with them it is no symmetry: Faraday rotation is odd in B); the weak-field V against -g_eff u_B (cos theta / mu) dI/dvz, the derivative with respect to a uniform velocity (no wavelength derivative), converges at fourth order on a grey atmosphere, 1.9e-6 of max|V| at 320 points. Turning the field azimuth for a vertical ray must turn (Q, U); the polarized Bezier solver does so only to discretization error (7e-3 of max|Q,U| on FAL-C, 5e-6 on the 320-point grey atmosphere, third order), because it limits each Stokes component and each element of K separately; with the limiter off the rotation is exact to 4e-14. Derivations in `docs/unit-tests`, new section "End to end".

### 2026-10-19 — Stress cases on cuts of a Bifrost snapshot (Q-18)

`tests/make_bifrost_atmos.py` cuts a Bifrost snapshot (en024048_hion, Carlsson et al. 2016) into a 1-D column, a 2-D row or a 3-D box for rhf1d, rhsc2d and rhsc3d, with the field file, in a few seconds and without astropy. RH's horizontal boundaries are periodic, like the snapshot's, so a full row or the whole box sampled every S-th point (S dividing 504) stays periodic; windows are for timing only. `run_stress.sh --bifrost DIR` (or `RH_BIFROST_DIR`) adds five convergence cases and two timing cases; the data are not distributed and the cases run locally only. Ca II in CRD converges (column 260: 19 iterations; the 2-D row every 4th column, 126 x 89: 25 iterations, 68 s on 16 threads; the 3-D box every 21st point, 24 x 24 x 89: 25 iterations, 12 min, 7.5 GB). With PRD the iteration diverges (B-29); the `_prd` cases reproduce it.

### 2026-10-19 — Zeeman splitting and polarized profiles tested from first principles; NaN profiles for J = 0 -> 0 lines and fields along a ray (Q-14, B-27, B-28)

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

### 2026-10-18 — Polarized DELO-Bezier solvers fixed and tested (changes results) (B-25, B-26, Q-14)

- `test_formal_stokes` (new) solves the 1-D polarized transfer equation
  with a constant propagation matrix K = 1 + K' and an unpolarized
  polynomial source function, where the Stokes vector is known in closed
  form: the scalar solution with mu replaced by mu K^-1
  (`stokes_exact.h`; I- through a matrix exponential). It also requires
  the polarized solver to equal the unpolarized Bezier solver when
  K' = 0, and measures the order with a depth-dependent K from
  successive grids. `test_formal_2d` and `test_formal_3d` do the same for
  the 2-D and 3-D polarized solvers on their uniform slabs, which also
  tests `SolveLongStokes` (the polarized variant of B-24).
- B-25, all geometries: the derivative terms of the control points were
  paired with the wrong weights (gamma S'(c) - theta S'(u) instead of
  gamma S'(u) - theta S'(c)): second order, and with K' = 0 2.4% from the
  unpolarized solver at 81 points. The sign of dK'/dtau at the central
  point was wrong: first order (1.02) when K varies with depth. Now
  fourth order for the emergent Stokes vector, constant or varying K
  (9.9e-7 at 161 points, as unpolarized: 9.3e-7).
- B-26, 2-D and 3-D: in the last layer of each ray the local source
  vector was used before it was set (2-D: the last point of the layer
  below; 3-D: an uninitialized local). With K' = 0 the bottom layer of
  I- was 38% off; for linear S the emergent vector 2e-9 (exact now).
- Converged Ca II full Stokes (CRD, 500 G at 45 degrees, ITER_LIMIT
  1e-6) on FALC refined 1x to 8x, against the fixed solver at 8x:

  | grid | 82 | 163 | 325 | 649 |
  |---|---|---|---|---|
  | V error, before | 3.8e-1 | 1.7e-1 | 7.2e-2 | 2.7e-2 |
  | V error, now | 2.9e-2 | 5.9e-3 | 1.1e-3 | — |
  | iterations, before / now | 54 / 28 | 79 / 34 | 64 / 46 | 100 / 70 |

  (V error: largest |V - V_ref| over the spectrum and rays, relative to
  max |V_ref|.) The old solver converged at order ~1.2 to the same
  limit; it overestimated the 854.2 nm V amplitude by 21% at 82 points.
  The fixed solver also needs fewer iterations at every resolution.

### 2026-10-18 — 1-D static atmospheres solved with S_INTERPOLATION instead of always Feautrier (changes results) (N-23)

- In 1-D, every wavelength without velocities, polarization or
  angle-dependent PRD was solved with second-order Feautrier, whatever
  `S_INTERPOLATION` said; only moving or polarized cases used Bezier3.
  Now both directions of every ray are solved with the chosen method
  (`S_BEZIER3` by default, or `S_LINEAR`), as in the angle-dependent
  case. `Feautrier()` stays in use by `rhsphere` and `solve1d`.
- The regression cases change by up to 30% (Ca II H and Lyα cores), but
  they stop after 4 iterations, with changes of 0.25 left: the
  difference is mostly a different iteration path. Converged
  (ITER_LIMIT 1e-6) on FALC refined 1x, 2x, 4x and 8x (82 to 649
  points; Ca II active), against Bezier3 at 8x:

  | method | 82 points | 163 | 325 | 649 | order |
  |---|---|---|---|---|---|
  | Feautrier (before), worst over spectrum | 2.0e-1 | 4.8e-2 | 1.2e-2 | 3.3e-3 | 2 |
  | Bezier3 (now), worst over spectrum | 3.1e-2 | 1.2e-2 | 1.9e-3 | — | ~3 |
  | Feautrier, Ca II 8498 core | 3.9e-2 | 9.8e-3 | 2.5e-3 | 6.0e-4 | 2 |
  | Bezier3, Ca II 8498 core | 1.5e-2 | 1.9e-3 | 2.0e-4 | — | ~3 |

  Both converge to the same limit; on the standard 82-point grid
  Bezier3 is 4 to 6 times closer to it. Cost: Bezier3 needs two sweeps
  per ray instead of one Feautrier solve (1.6x per iteration on this
  case) and 28 iterations instead of 22 (its approximate diagonal
  operator, N-19): 0.73 s instead of 0.36 s.

### 2026-10-18 — The unit-test problems and their derivations, written up (Q-20)

- `docs/unit-tests/unit-tests.tex` (`make -C docs/unit-tests` builds the
  PDF with pdflatex) describes every unit test: the problem, its
  analytic solution or reference, why the expected accuracy and order of
  convergence are what they are, and the defects each test found. Voigt
  and Faraday-Voigt functions (the hybrid algorithm and its first-order
  expansion in a); LU, interpolation, quadrature and special functions;
  1-D formal solvers on slabs with polynomial S (I+ and I- in closed
  form, without cancellation), the cubic Bezier weights and why Bezier3
  is fourth order for the emergent intensity, the exact diagonal of the
  Bezier operator (N-19); multi-D solvers on uniform slabs and the long
  characteristics (B-24); the spherical shell; the polarized solvers
  against the matrix form of the analytic solution, DELO-Bezier and the
  errors found there (B-25, B-26).
- `make_bezier_reference.py`: its docstring had the weights of the two
  control points swapped (gamma belongs to the upwind one).

### 2026-10-17 — Long characteristics in 2-D and 3-D stopped one section short (changes results) (B-24, Q-14)

- `test_formal_2d` and `test_formal_3d` solve a horizontally uniform,
  periodic slab, whose every column must reproduce the 1-D analytic
  solution. Where the short characteristics cross column lines (2-D) or
  vertical planes (3-D) before reaching the next layer, RH integrates
  back along a long characteristic. B-24: `SolveLong` (and
  `SolveLongStokes`, both geometries) treated the last crossing as the
  local point, returned the intensity one crossing upwind of the true
  upwind point, and so skipped the last section. On a 2-D slab with
  dz·|mux/muz| > dx in the deep layers this made the emergent intensity
  7% wrong for a linear source function (which every interpolation
  reproduces exactly), 15% different between identical columns; in 3-D
  11%. Inherited from upstream RH. Now the local point is the downwind
  end of the last section: I+ for linear S exact to 2e-15 in both.
- Measured with the tests (Nz = 161 over six decades of optical depth,
  27 points per decade): on a mesh where every
  upwind point lies on a layer, the 2-D and 3-D linear solvers equal the
  1-D one; the 2-D and 3-D Bezier solvers are second order (4.8e-4 where
  1-D Bezier gives 9e-7), because the upwind control point uses the
  secant slope (N-24). Across column crossings 3-D converges at order
  ~1.5 (vertical interpolation between two layers).
- The fix in `SolveLongStokes` (2-D and 3-D) is the same change; no test
  exercises the polarized long characteristics yet.
- The 1-D formal-solver entry of 2026-10-17 quoted its errors "at 16
  points per decade"; that grid (Ndep = 161 over six decades) has 27.
- `test_formal_sphere`: the spherical Feautrier solver along the rays of
  `getRays()`, in a shell of constant opacity and source function where
  every ray has a closed-form solution. It converges at order 1.4 to 1.6
  (formally 2, and 4 for Hermite), because the rays are sampled where they
  cross the radii and so coarsely near their tangent points; the grazing
  rays near the top are not refined by a finer radial grid at all (the
  outermost always has 3 points), with errors up to 7e-3.
- `rh_unit_test()` takes the geometry directory, so that tests can use the
  2-D and 3-D `geometry.h`; the analytic slab solution is shared in
  `tests/unit/slab_exact.h`.

### 2026-10-17 — One Voigt function, accurate to 2e-6, for profiles with and without polarization (changes results) (Q-11)

- `Voigt(a, v, &F)` returns H and the Faraday-Voigt function F together,
  from a hybrid of Weideman's rational approximation (N = 32) inside
  |v| + a < 8 and the Laplace continued fraction (6 levels; Humlicek's
  region I is its first level) outside, after Schreier (2018); H =
  exp(-v²) exactly for a = 0. Against mpmath: relative error of H below
  1.6e-6 wherever H ≥ 1e-12 H(a, 0), 2.7e-10 relative to line center, F
  within 1e-9; about 50 ns per call.
- It replaces all five algorithms of RH 2 (ARMSTRONG for profiles without
  polarization, 1.2e-4; HUMLICEK with polarization, 8e-5; RYBICKI and
  HUI_ETAL in the PRD redistribution code; the unused lookup table), so a
  field-free run and a B = 0 Stokes run use the same profile. `humlicek.c`
  and the `VoigtAlgorithm` argument are gone.
- Results change by up to 2.3e-5 in I and 7.6e-5 in Q, U, V in the
  regression cases (relative to the peak), the size of the old profile
  errors.

### 2026-10-17 — 1-D formal solvers tested on analytic slabs; accurate Bezier and linear weights (Q-14, B-5, N-17, N-18)

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

### 2026-10-17 — Unit tests for the numerical kernels (Q-14)

- Four test programs against analytic results and mpmath references
  (`make_kernels_reference.py`): `test_linalg` (LU, pivoting, Hilbert,
  inverse), `test_interp` (Hunt/Locate, Linear, BiLinear, natural and
  exponential splines, cubic convolution), `test_quadrature` (Gauss-Legendre
  to n = 100, w2/w3/U3), `test_functions` (E1, E2, Planck, 3j and 6j
  symbols). Written by a subagent, checked against the suite; no library
  changes.
- Findings, logged as B-21 to B-23 and N-20 to N-22: `BiLinear` reads past
  the table at its upper edge; `w3js` accepts invalid arguments and
  overflows from j = 50; `Hunt` is one index low on exact node values;
  Gauss-Legendre weights only to ~1e-14; LU refinement gives no gain and
  misses singular matrices; Planck without `expm1`, E1/E2 zero above x = 80.

### 2026-10-17 — Unit tests, starting with the Voigt functions; VoigtHui finite at a = v = 0 (Q-14, B-20)

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

### 2026-10-17 — CI compares each change with its base commit; drift from upstream is reported (Q-14)

- The regression job builds the base commit of the change (the target
  branch before it) and runs `run_tests.sh` against it in tolerance mode.
  Results must match, unless the pull request carries the label
  `changes-results` (an intended change of results): then the differences
  from the base are reported in the job summary instead of failing, while
  sanitizer, thread-count and memory checks still fail.
- The comparison with upstream RH (c95bdf6) remains, as a report that never
  fails. Upstream is the reference we trust but check: every difference
  from it needs an explanation in this log (a known numerical change or a
  fixed bug). The analytic unit tests guard against slow drift.
- Why: the first fixes that change results on purpose (Bezier weights,
  N-17: up to 7e-7) would otherwise fail every later comparison with
  upstream, or force a tolerance too loose to catch accidents.

### 2026-10-16 — Runs report whether they converged, also in the exit status (breaking) (N-1)

- `Iterate()` now ends with "Converged in N iterations: max change … <
  ITER_LIMIT = …", or with a WARNING "NOT CONVERGED after N iterations: max
  change … >= ITER_LIMIT = …" saying that the output is not a converged
  solution. Before, a run that hit `N_MAX_ITER` ended silently and wrote
  its output as if it had converged. All four geometries.
- The exit status says how the run went (documented in `error.h` and the
  README):

  | Status | Meaning |
  |---|---|
  | 0 | completed and converged |
  | 1–63 | completed, with problems, as a bitmask; output is written. 1 = populations not converged |
  | 70 | fatal error (EX_SOFTWARE in sysexits.h); see the TERMINATING_ERROR message |
  | 128+n | killed by signal n (136: trapped floating-point exception) |

  Problems of a completed run can occur together, hence the bitmask; six
  bits are available, one is used so far (PRD and electron-density
  convergence get theirs when they are checked). Fatal errors exit with one
  code instead of 2 (`ERROR_LEVEL_1`) and 3 (`ERROR_LEVEL_2`): they are rare,
  a user cannot act on a category, and the message says what went wrong.
- A run with a deliberately small `N_MAX_ITER` now exits with 1; scripts
  that test for success should test for 0, or accept bit 1.
  `run_tests.sh` and `analyze_run.py` treat 1–63 as completed.
- Verified: no input → 70; CaII CRD at `N_MAX_ITER = 4` → 1 with the NOT
  CONVERGED warning, at 200 → 0 after 10 iterations; the regression suite
  is byte-identical to the previous build.

### 2026-10-16 — Solver redesign planned (S track); plan steps 3 to 6 revised

Plan step 4 (convergence robustness) is replaced by a redesign of the
iteration scheme instead of patches to it. Design, rationale, cost model,
traps and work items are in `docs/SOLVER.md`; the work items are S-1 to S-9
in §3.8.

- Why: RH's Ng is Anderson acceleration (type II) solved through normal
  equations without safeguards, weighted so that it fits a different
  quantity than the convergence test measures; the PRD and ne loops are a
  block splitting of one coupled system, which can cycle or diverge
  although the coupled problem has a solution (N-2, N-4, N-6, N-10).
- What: a monolithic residual solved by globalized inexact Newton–Krylov
  with MALI as preconditioner, built in stages. First a fair comparison
  (S-1 measuring stick, S-2 deterministic MALI map), then stabilized
  Anderson (S-3) against Newton–GMRES on the MALI map (S-4), decided by
  the S-5 experiment with criteria fixed beforehand.
- Plan step 3 now separates two test tiers: regression (must pass; a
  non-default solver is checked on converged solutions) and stress
  (convergence and performance, including RMHD cut-outs, Q-18).
- Steps 5 and 6: the items that depend on how Γ is built and stored (P-2,
  P-5, P-6, M-2, M-3) wait for S-5; the rest is unaffected.
- `docs/make_issues.py` knows the S area (label `solver`, plan step 4).

No code changes; the regression suite is unaffected.

### 2026-10-16 — One code style, enforced by clang-format

The code followed Emacs's C mode (2-space indentation with leading runs
of 8 spaces stored as tabs), which looked broken in any editor not set
to 8-column tabs. All C sources are now formatted by clang-format
(`.clang-format`): K&R braces, one tab per level shown as 4 columns,
spaces for alignment, aligned declarations and assignments, no line
limit yet (existing line breaks kept). `.editorconfig` sets tab width 4
for editors and GitHub, a CI job rejects unformatted code, and
`.git-blame-ignore-revs` lets `git blame` skip the reformat commit. The
reformat changes no machine code (all 27 programs compile identically).
With 4-column tabs, 960 lines exceed 80 columns (56 before), mostly in
deeply nested loops; they are left for the structural work (step 7).
Function and variable naming remains inconsistent (camelCase,
Capitalized, snake_case); also for step 7.

### 2026-10-16 — MULTI atmospheres from readtrujillo have velocities in km/s (B-17)

`writeMULTIatmos()` (used by `readtrujillo` and `readkurucz`) wrote the
macroscopic velocity in m/s into a column that RH reads as km/s, so a
`readtrujillo` atmosphere with a velocity field came out with velocities
1000 times too large (`readkurucz` writes zero velocities). For a geometric
depth scale it wrote an uninitialized height. Velocities are now written
in km/s and heights in km, as `rhf1d` reads them. With this, gcc reports
no warnings at all under `-Wall -Wextra`.

### 2026-10-16 — Dead code removed (Q-10)

Seven source files that no build used are gone (about 2,800 lines):
`collision_15d.c` and `collision_Oslo.c` (older variants of
`collision.c`; the first did not compile), `rhsphere/readinput.c` and
`rhsphere/solveplane.c`, and the tools `constructatom.c` (built model atoms
from a Kurucz line list, levels and lines only; it made `Atoms/SiI.atom`,
did not compile since about 2000, and hard-coded silicon ion labels),
`adjustbf.c` and `lande.c`. Earlier, `listincl` and `listincl_f` went with
the old makefiles. A new model-atom tool is on the wishlist (F-1).

### 2026-10-16 — RH compiles without warnings; CI treats warnings as errors (Q-3)

`-Wall -Wextra` gave 381 warnings with gcc 8 (more with clang and gcc 13).
All are gone, in gcc 8, gcc 13 and clang (Linux and macOS), and the CI
builds now use `-Werror` (new option `RH_WERROR`), so they stay gone.

- Unused variables and parameters removed (about 230, mostly by Sonnet
  subagents). `Error()` now tells the compiler that the error levels do
  not return. 120 false-positive "maybe uninitialized" warnings were
  checked one by one and initialized.
- The triage found real bugs, fixed separately: collision data before a
  `TEMP` line (B-16, silently wrong rates), velocities 1000 times too large
  in `readtrujillo` output (B-17), about 30 unchecked I/O results (Q-17),
  an unbounded `sprintf` into a label buffer in `topbase`, a NULL `qsort`
  in `readMolecule`, and missing validation of the radial grid size in
  `rhsphere`.

The regression suite stayed byte-identical through all of it.

### 2026-10-16 — Collision data before a TEMP line is an error; no inherited item counts (B-16)

Tabulated collision keywords (OMEGA, CE, CI, CP, CH0, CH+, CH) give one
value per temperature of the last `TEMP` grid, but read as many values as
the last keyword of any kind had set. Two consequences for atom files:

- Collision data before any `TEMP` line read an uninitialized count and a
  NULL temperature grid: a segmentation fault. It now stops with
  "OMEGA collision data for levels 1 0 before any TEMP line".
- After a fixed-size keyword (AR85-*, BURGESS, SHULL82, BADNELL, SUMMERS)
  a tabulated line read only that keyword's number of values and
  interpolated in a truncated grid, without any message. With a BURGESS
  line placed before the first OMEGA line of `CaII.atom` the Ca II
  849.8 nm line changed by 6.6 %. Tabulated lines now always use the
  `TEMP` grid size.

No shipped atom file is affected: the regression suite is byte-identical.

### 2026-10-16 — Continuous integration on Linux and macOS; Voigt complex arithmetic

GitHub Actions now builds every push and PR to `modernize` on Ubuntu (gcc,
clang) and macOS 14 (Apple silicon, clang) with `-Wall -Wextra`, runs the
smoke tests, and runs the regression suite against upstream c95bdf6 built
on the same runner (`tests/build_upstream.sh`), on Ubuntu and macOS and with
ASan+UBSan.

The first runs found:

- `tools/readtrujillo.c` did not compile with clang or on macOS
  (`void main(int, void *[])`, `<malloc.h>`).
- UBSan: `readMolecule()` passed a NULL array to `qsort` for a molecule
  without lines (undefined behaviour); fixed. Otherwise the sanitizers are
  clean, the first time they ran on RH.
- The Stokes cases on macOS were not byte-identical to upstream.

Complex arithmetic in the Voigt code: upstream divided complex numbers with
the textbook formula (`cmplx_div`). With that formula (an explicit
`cdiv()`, or gcc's `-fcx-limited-range`) gcc on Linux is byte-identical to
upstream; on macOS the Stokes cases differ even then (most likely FMA
contraction, which clang enables by default on arm64). Since exact
agreement across compilers is not a goal, the Voigt code now uses C's
default complex arithmetic and the CI regression runs in tolerance mode.
The Stokes cases then differ from upstream by at most 1e-10 relative
(I: 2e-11), the size of thread-count noise. For pure refactors, compare
exactly against a build of the PR's base instead of upstream.

### 2026-10-16 — CMake replaces the hand-written makefiles (breaking) (Q-1, Q-15)

RH now builds with CMake (3.20 or newer): `cmake -B build && cmake --build
build -j`. The seven makefiles and `tests/build_asan.sh` are gone.

- Parallel builds are safe (the archive-member rules raced under `make -j`
  and could leave a corrupt `librh.a`), and header dependencies are tracked
  by the compiler instead of by hand. A full build takes 3 s on 16 cores.
- Out-of-source builds: programs land in `build/rhf1d/rhf1d` etc. Installed,
  the auxiliary programs that several geometries share get a suffix
  (`solveray_1d`, `backgrcontr_2d`).
- The default build type is RelWithDebInfo (`-O2 -g`), so FPE backtraces
  (N-11) resolve to file and line (Q-15). Options: `RH_NATIVE`,
  `RH_SANITIZE`, `RH_WARNINGS`, `RH_BUILD_TOOLS`.
- Portability: XDR from libtirpc (found with pkg-config, or in
  `/usr/include/tirpc`) or from the C library (macOS); `feenableexcept`
  and `execinfo.h` are detected instead of keyed on `-DLinux`/`-DDarwin`,
  so traps work on any glibc system and are skipped cleanly elsewhere.
  Intel compilers get `-fp-model=precise`: with icx's default fast
  floating point the wavelength grids differed and the Kurucz case died
  with SIGFPE.
- `run_tests.sh` and `run_stress.sh` accept a CMake build directory or a
  source tree with `build/`, as well as in-place makefile builds (upstream).

Verified: `run_tests.sh --exact` against upstream is byte-identical in all
34 cases (this needed `-fcx-limited-range` for `voigt.c` and `humlicek.c`,
which the old Makefile set per file). Intel icx 2022.1 matches upstream
within the suite's tolerances (1e-9; 1e-7 in 2-D) with no warnings. The Unix
Makefiles generator with `make -j32` produces the same machine code as
Ninja. Sanitizer builds could not be run locally (no sanitizer runtime is
installed for gcc 8 or icx); CI will cover them.

### 2026-10-16 — I/O results are checked (Q-17)

About 30 places stored the result of an `fread`/`fwrite`, `xdr_*`,
`fscanf`/`sscanf` or `fgets` call and never looked at it, so a truncated
or malformed file, or a failed write, went unnoticed (e.g. a short
Barklem table in `readBarklemTable`). Each now stops with a message that
names the file or quotes the offending line. Mostly in the auxiliary
programs (`solve1d`, `solve2d`, `formal2d`, `formal3d`, `solveray`,
`chartflow`, `backgrcontr`), the 2-D/3-D atmosphere readers and geometry
writers, and the atomic-data tools. XDR calls that belong together are
checked once per group. Behaviour on valid input is unchanged (regression
suite byte-identical).

### 2026-10-15 — calcDistance always returns a value (B-8)

`calcDistance()` (3-D stencils) had no return for an unexpected
intersection plane (-Wreturn-type). All three planes were handled, so this
could not happen in practice; an invalid plane now stops with an error
instead of returning garbage. Verified: the warning is gone; the 3-D
regression case is byte-identical.

### 2026-10-15 — Scratch-file I/O checks every transfer (B-9)

The issue as logged (an `int`/`size_t` comparison mishandling -1) turned out
to be benign: -1 converts to SIZE_MAX and never matched the record size. The
real problems in `readj.c` were elsewhere:

- `pread_rh()`/`pwrite_rh()` (the chunked I/O used for the background
  opacity file) summed the return values, so an error (-1) wrapped the
  unsigned total, and partial transfers were never retried. They now share
  one loop that retries partial and interrupted transfers and reports
  failure, and all J, J20, Imu and profile records go through it too.
- Error messages now include the system error text.
- `readJgas()`/`writeJgas()` (hybrid PRD) did not check that `Jgas.dat`
  opened or that the full array was transferred.

Verified: `readj.c` has no warnings under `-Wall -Wextra` (12 before); the
regression suite is byte-identical to upstream.

### 2026-10-15 — Floating-point exceptions are reported with a backtrace (N-11)

On Linux RH traps division by zero, invalid operations and overflow, but a
trapped exception killed the run with only "Floating point exception". A
SIGFPE handler now prints the kind of exception, the faulting address and a
backtrace (resolve with `addr2line -f -e <program> <address>`), then
re-raises the signal so the exit status is unchanged. The handler only uses
async-signal-safe output.

Verified: the PRD_NG_PERIOD = 0 crash (B-2) now reports "integer division
by zero" with a backtrace that addr2line resolves to
Accelerate ← Redistribute ← Iterate; the regression suite is
byte-identical. File and line numbers need `-g` in the default compiler
flags (to be added with the build-system work, Q-1).

### 2026-10-15 — Parabolic interpolation removed (S_PARABOLIC, DELO_PARABOLIC) (N-5)

**Breaking change.** `S_INTERPOLATION = S_PARABOLIC` and
`S_INTERPOLATION_STOKES = DELO_PARABOLIC` are removed; either value now stops
the run with a message naming `S_BEZIER3`/`S_LINEAR` and `DELO_BEZIER3`.

- Why: the quadratic interpolation is not monotonic. It overshoots near steep
  gradients (transition region, shocks, velocity-shifted line cores), and its
  approximate operator grows like 1/dtau_dw, which makes the iteration
  diverge: CaII on a ±10 km/s FAL-C crashed with SIGFPE in 3 iterations, and
  bounding Psi* alone (PR #73, closed) did not cure it because the intensity
  itself overshoots. Monotonic cubic Bezier (de la Cruz Rodríguez & Piskunov
  2013) is the higher-order option, linear the robust one.
- Removed `Piecewise_1D/2D/3D` and `Piece_Stokes_1D/2D/3D` (the file
  `rhf1d/piecestokes_1D.c`); kept the shared helpers `Quadr`, `StokesK_2D`
  and `StokesK_3D`. `formal2d` and `chartflow` now use the 2-D Bezier solver.
  The 2-D example input uses `S_BEZIER3` (the default).
- Verified: all programs of all geometries build; the regression suite is
  byte-identical to upstream (the 2-D example now uses Bezier in both);
  stress case `1d_v10_ca_par` replaced by `1d_v10_ca_lin`, which converges
  in 19 iterations (Bezier: 14).

### 2026-10-15 — Memory test uses a valid tiled atmosphere (Q-12)

- `tests/make_tiled_2d_atmos.py` read the 2-D arrays as [Nx][Nz], but RH
  stores them with x fastest ([Nz][Nx]), so every tiled column mixed heights
  from different depths. The resulting atmosphere gave negative opacities, and
  whether a build survived the `2d/mem` case depended on luck: upstream itself
  crashed with SIGFPE on it once the 2-D example used Bezier.
- It now tiles along x; all 192 columns are proper FAL-C columns.
- Verified: `run_tests.sh --only 2d/mem` passes for this tree against upstream
  (208 MB vs 206 MB) and for upstream against itself.

### 2026-10-15 — SOLVE_NE = ITERATION works with active atoms (B-1)

Two bugs, the second hidden by the first, made every run with
`SOLVE_NE = ITERATION` and an active atom fail when the background was
recomputed after the first iteration:

- `SetLTEQuantities()` re-read the collisional data from the atom file
  without rewinding it, so the second read hit end of file ("Reached end of
  datafile"). It now seeks back to `atom->offset_coll`, which `readAtom()`
  saved for this but nothing used. `CollisionRate()` also allocated a new
  `atom->C` (Nlevel² × Nspace) on every call although `readAtom()` had
  already allocated it; it now reuses it, which removes that leak.
- `Background()` freed the Kurucz partition-function tables at the end of
  every call, but `Solve_ne()` needs them on the next one: segfault in
  `getKuruczpf()`. They are now kept when `Background()` will run again
  (electron density iteration or hydrostatic equilibrium), as the passive
  atoms already were.

Verified: stress case `1d_falc_solvene` (FALC, H_6 + CaII,
`SOLVE_NE = ITERATION`) converges in 29 iterations; before, it crashed after
the first. Hydrostatic runs give the same result as upstream. 1-D regression
cases are byte-identical to upstream.

### 2026-10-15 — PRD Ng acceleration no longer divides by zero (B-2)

- `PRD_NG_PERIOD` now defaults to 1, like `NG_PERIOD`; it defaulted to 0, so
  setting only `PRD_NG_ORDER` crashed with SIGFPE in `Accelerate()`.
- `readInput()` rejects `NG_PERIOD` or `PRD_NG_PERIOD` below 1 when the
  matching order is positive, with a clear message instead of a crash.
- Verified: stress case `1d_falc_prdng` (CaII, `PRD_NG_ORDER = 2`) now
  converges in 12 iterations instead of crashing after 1; an explicit
  `PRD_NG_PERIOD = 0` stops with the new message. 1-D regression cases are
  byte-identical to upstream (default runs do not use PRD Ng).

### 2026-10-15 — Collisional relaxation runs its full schedule (B-3)

- With `COLL_RELAX_NSTEP`/`COLL_RELAX_FACTOR` the collision rates are scaled
  by up to `factor^(1 - Nstep)` (1e38 for 0.01 and 20 steps) in the first
  iterations, so populations stay near LTE and hardly change. The convergence
  test fired at iteration 1 and the run "converged" to LTE. The test now
  starts at iteration `COLL_RELAX_NSTEP`, the first unscaled one.
- Ng acceleration restarts after every scaled iteration, so it only
  extrapolates from iterations of the actual problem.
- Verified: stress case `1d_falxco_collrel` (FALXCO, H_6, 20 steps, factor
  0.01) converges in 124 iterations; before, it stopped after 1 with LTE
  populations. The same problem without relaxation does not converge in 200
  iterations, so relaxation now rescues this case. Restarting the unrelaxed
  problem from the relaxed solution reconverges in 30 iterations to
  populations within 2.7 % and intensities within 0.27 % (PRD ρ and J start
  over on a restart; see N-3 for what the 1e-3 criterion guarantees). Default
  runs are byte-identical to upstream.

### 2026-10-15 — Clear errors instead of SIGFPE for bad atom data (B-4)

`CaII_large.atom` crashed RH at setup with a bare SIGFPE. It has two kinds
of bad data, and each one crashed a spline:

- Four bound-free continua list a cross-section row twice (same wavelength,
  same value). `readAtom()` now drops exact repeats with a warning, stops
  with a clear message if a repeated wavelength has a different cross
  section, and checks for strictly increasing wavelengths (the old check
  let equal neighbours through to a zero-width spline interval).
- 13 `OMEGA` collision records hold `Inf` at every temperature.
  `CollisionRate()` now checks all tabulated rates and stops with a message
  naming the atom, record and levels.

Verified: `1d_falc_cail` now stops at once with "Atom CA: collision record
OMEGA 11 10 has a non-finite value (inf) at T = 1264.9 K" instead of
SIGFPE. No other shipped atom file has repeated continuum wavelengths or
non-finite rates (CI-VI_Judge's are commented out). 1-D regression cases are
byte-identical to upstream. The atom file itself still needs its 13 `OMEGA`
records fixed before CaII_large can be used.

### 2026-10-15 — Fortran removed (Q-2, #53)

- Deleted `humlicek_.f90`, `hui_.f90`, the `HAVE_F90` branches in
  `voigt.c`, `librh_f90.a` and the gfortran rules; no makefile builds or
  links Fortran any more.
- Verified: `run_tests.sh --exact` against upstream built without
  `HAVE_F90` passes 34/34, byte-identical including full-Stokes and
  field-free. Against the default (Fortran) upstream build, Stokes Q/U
  differ by up to 1.5e-6 of the peak, which is the Fortran/C rounding
  difference.
- Code map: the "Files in and out" and drivers sections no longer need to
  mention `librh_f90.a`.

### 2026-10-15 — Baseline and analysis

- Ran the existing regression suite and the new stress suite against
  upstream; results in §2.
- Added `tests/run_stress.sh`, `tests/analyze_run.py`,
  `tests/make_perturbed_atmos.py`, `tests/atmos/`.
- Verified that the Stokes difference vs upstream comes only from the
  Fortran vs C Voigt routines (upstream built without `HAVE_F90` is
  bit-identical). Decision: phase out the Fortran code.

### Pre-existing uncommitted work (from an earlier session)

Present in the working tree at the start of this log; not yet committed.

- Per-thread Γ and rate accumulators merged after each sweep instead of
  mutex-protected updates (`fillgamma.c`, `iterate.c`, `atom.h`).
- PRD scattering integral and line profiles threaded with `parallelFor`
  (`parallel.c`, `scatter.c`, `profile.c`); PRD weights read in large windows.
- C99 `complex` replaces `complex.c`/`complex.h`; Fortran Voigt no longer
  built at top level.
- Kurucz line lists with a leading blank line no longer crash.
- Verified: bit-identical to upstream at 1 thread on all non-Stokes cases; 1-D
  wall time 2.4–3.6× faster at 4 threads; peak memory +5–30 % (M-9).
