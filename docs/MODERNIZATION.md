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

- **State of the code (2026-10-15):** the uncommitted threading work is
  correct. It is bit-identical to upstream at 1 thread in every regression and
  convergence case, and only the Fortran-vs-C Voigt choice changes Stokes
  results.
- **Robustness is the weakest area.** Of 19 one-dimensional convergence cases,
  8 converge cleanly, 1 converges slowly, 3 fail to converge, and 7 crash or
  return a wrong answer. Neither tree reports non-convergence. A 3-D
  chromospheric test diverges in both.
- **Memory blocks realistic 3-D runs:** line profiles stored per ray, dense
  rate and collision matrices, per-thread copies and per-point PRD weight files
  put a 256²×200 run with five atoms in the tens of TB.
- **Speed:** the formal solution and per-transition work dominate; about
  40 % of a PRD run is serial; `-O3 -march=native` alone gives 14 %.
- **Code quality:** 354 compiler warnings, a quarter of the code duplicated
  across the four geometries, global state everywhere, no unit tests or CI.
- **Next step:** fix the crashes and silent failures (plan step 1), then the
  build and Fortran removal, before any numerical or structural change.

Issue IDs: **B** bug or crash, **N** numerics/convergence, **P** speed,
**M** memory, **T** threading, **Q** code quality / build. Severity: high,
medium, low. Status: open, fixed (with changelog entry), wontfix.

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

_Perf-set results (2-D 128 columns and 3-D 24×24 / 16×16 with several
many-level atoms) go here (pending)._

---

## 3. Issues

### 3.1 Bugs and crashes

| ID | Sev | Issue | Where | Evidence | Status |
|---|---|---|---|---|---|
| [B-1](https://github.com/dwijn/RH/issues/2) | high | `SOLVE_NE = ITERATION` aborts on the second `Background()` call ("Reached end of datafile"): collision data re-read from an open file without rewinding; `offset_coll` saved but never used; `atom->C` reallocated (leak) | `ltepops.c:239`, `readatom.c:568-569`, `collision.c:482` | crashes upstream and WIP | open |
| [B-2](https://github.com/dwijn/RH/issues/3) | high | `PRD_NG_ORDER > 0` with default `PRD_NG_PERIOD = 0` divides by zero | `readinput.c:99`, `accelerate.c:88` | SIGFPE | open |
| [B-3](https://github.com/dwijn/RH/issues/4) | high | Collisional relaxation (`COLL_RELAX_NSTEP/FACTOR`) stops after one iteration with LTE populations: convergence test fires during the relaxation phase while collisions are scaled by up to 1e38 | `iterate.c:132`, `fillgamma.c:59` | FALXCO H_6: 1 iteration, delta 1.4e-12, LTE result | open |
| [B-4](https://github.com/dwijn/RH/issues/5) | high | `CaII_large` crashes at setup (SIGFPE in `splineCoef` from `SortLambda`) | `sortlambda.c:387` | upstream and WIP | open |
| [B-5](https://github.com/dwijn/RH/issues/6) | medium | Linear short-characteristics end point uses stale weights (`w2` not recomputed) | `rhf1d/piecewise_1D.c:126` | flux changes up to 1e-3 when fixed | open |
| [B-6](https://github.com/dwijn/RH/issues/7) | medium | Possible double counting of the cross-coupling term for multi-ion atoms (added once per transition sharing an upper level) | `fillgamma.c:279` | MgI+II converged populations shift 3e-3 at level 70 when changed; correct form not yet established | open — verify against an independent code (e.g. Lightweaver) |
| [B-7](https://github.com/dwijn/RH/issues/8) | medium | 3-D: `Jgas` accumulated with `wmu` inside the to_obs loop, looks double counted | `rhsc3d/formal.c:370-374` | code reading | open |
| [B-8](https://github.com/dwijn/RH/issues/9) | medium | `calcDistance()` falls off the end without a return value | `rhsc3d/geometry.c:265` | -Wreturn-type | open |
| [B-9](https://github.com/dwijn/RH/issues/10) | medium | `pread`/`pwrite` results compared as `int` against `size_t`; -1 error mishandled | `readj.c:58,82,146,171,200,229,413,442` | -Wsign-compare | open |
| [B-10](https://github.com/dwijn/RH/issues/11) | medium | 32-bit overflow of `Nlambda*Nspace` sizes (`NgInit(int N)`, `Nlamu`) once > 2³¹ (e.g. >163 PRD wavelengths at 13 M points) | `accelerate.c:34`, `redistribute.c` | code reading | open |
| [B-11](https://github.com/dwijn/RH/issues/12) | low | Leaks: `backgrflags`/`backgrrecno` re-malloc'd every `Background()`; `spectrum.I` allocated twice in 3-D with Stokes | `background.c:300,312`, `initial_xdr.c:386,392` | code reading | open |
| [B-12](https://github.com/dwijn/RH/issues/13) | low | AR85-CHP collisions hard-code `atmos.H->n[5]` as protons (wrong unless H has 6 levels; no shipped atom uses it) | `collision.c:883` | code reading | open |
| [B-13](https://github.com/dwijn/RH/issues/14) | low | `ERROR_LEVEL_1` for an unknown solver leaves I and Psi undefined and continues | `rhf1d/formal.c:209-231` | code reading | open |
| [B-14](https://github.com/dwijn/RH/issues/15) | low | `error.c` `sprintf`s a 512-byte message plus routine name into a 512-byte buffer | `error.c:51` | code reading | open |

### 3.2 Numerics and convergence

| ID | Sev | Issue | Where | Evidence | Status |
|---|---|---|---|---|---|
| [N-1](https://github.com/dwijn/RH/issues/16) | high | No convergence report: at `N_MAX_ITER` the loop ends silently and output is written as if converged; log has no converged/not-converged line; exit status 0 | `iterate.c:132`, `rhf1d/rhf1d.c:91-131` | all runs | open |
| [N-2](https://github.com/dwijn/RH/issues/17) | high | Ng acceleration without safeguards: one global coefficient set for all levels/depths weighted by 1/n; no positivity, residual or conditioning check; history not reset after relaxation, ne or HSE updates | `accelerate.c:88-122` | He_large Ng jump 0.34 → 20.7; FALC_shock H_6 SIGFPE with Ng, converges in 115 without | open |
| [N-3](https://github.com/dwijn/RH/issues/18) | high | Convergence criterion is not an accuracy guarantee: max relative change over all levels/depths, dominated by trace levels, measures Ng jumps, ignores PRD ρ, background J and ne | `maxchange.c:43`, `iterate.c:132` | He_large "converged" to 1e-3 at 1 thread (108 iterations) and 4 threads (48) differ by 1.4 %; MgI+II at 1e-6 vs 1e-9 differ by 2e-5 | open |
| [N-4](https://github.com/dwijn/RH/issues/19) | high | PRD ρ iteration is Λ-iteration with populations fixed; early exit at `PRD_ITER_LIMIT` makes the map discontinuous → limit cycles; ρ never clipped (ρ < 0 possible after Ng) | `redistribute.c:74-105`, `scatter.c:346,765` | H_6 FALC at 1e-8 stalls at 3–9e-4 for 600 iterations; converges in 142 with fixed 3 sub-iterations; FALXCO H_6 (PRD) no convergence, H_6_CRD converges in 17 | open |
| [N-5](https://github.com/dwijn/RH/issues/20) | high | Parabolic short-characteristics Ψ* unbounded (≈ 1 + 1/Δτ_d when Δτ_u ≫ Δτ_d); fallback to linear only on I < 0. `S_PARABOLIC` is the shipped 2-D default | `rhf1d/piecewise_1D.c:234`, `rhsc2d/piecewise_2D.c:309`, `rhsc2d/run/keyword.input:161` | Ψ* up to 12 (> 4000 before blow-up), 11 % of points outside [0,1]; CaII/MgII on FALC_v10 SIGFPE in 3–8 iterations; Bezier/linear converge in 14–28 | open |
| [N-6](https://github.com/dwijn/RH/issues/21) | high | No ne preconditioning: `SOLVE_NE` couples ne by plain iteration (dfjk = 0 for NLTE elements); Newton tolerance 1e-2 | `solvene.c:40,165` | (with B-1 patched) H_6_CRD period-2 oscillation at 0.22 after 80 iterations; 2-cycle on FALC_shock | open |
| [N-7](https://github.com/dwijn/RH/issues/22) | medium | Rate-matrix conditioning: LU with row scaling only; iterative improvement in the same precision; LTE relative to the neutral ground level (overflow for high stages); n* clipped to DBL_MIN | `ludcmp.c:73,142`, `ltepops.c:98` | CI–VI Judge (115 levels): 230 negative populations after 4 iterations, b up to 1e307, divergence, SIGFPE | open |
| [N-8](https://github.com/dwijn/RH/issues/23) | medium | No positivity enforcement or damping anywhere | — | grep | open |
| [N-9](https://github.com/dwijn/RH/issues/24) | medium | Starting solutions: `ZERO_RADIATION` has no photoionization (over-recombined start); no J = B or escape-probability option | `initial_xdr.c:21,588-640` | LTE start 40 iterations vs 22 (H_6 FALC) | open |
| [N-10](https://github.com/dwijn/RH/issues/25) | high | Multi-D atmospheres with a horizontally corrugated transition region do not converge: 3-D diverges and crashes, 2-D with the same gradients locks into a period-3 cycle | — | 3-D 8×8×82, H_6+CaII, upstream and WIP. Unperturbed and velocity-only variants converge. With the height/temperature perturbation alone it diverges (SIGFPE) with bicubic or linear 3-D interpolation, Bezier or linear S, PRD or CRD atoms, Ng on or off, and even with 4× gentler gradients (crash at iteration 28). 2-D with equally steep gradients: no crash, deltas cycle 23 → 1 → 2 with period 3 = `NG_PERIOD` for 40 iterations; 4× gentler 2-D converges slowly (1e-2 after 40). Cause not yet isolated: candidates are the 3-D short-characteristics Ψ*, Ng (N-2) and the lack of positivity safeguards (N-8) | open |
| [N-11](https://github.com/dwijn/RH/issues/26) | medium | FPE traps on with no handler: any overflow kills the run with no output or diagnostic | `fpehandler.c:37` | all SIGFPE cases | open |
| [N-12](https://github.com/dwijn/RH/issues/27) | low | Collisions with H/protons use H populations frozen at setup; spline interpolation of rates can overshoot negative and clamps outside the table | `collision.c:685,726-746` | code reading | open |
| [N-13](https://github.com/dwijn/RH/issues/28) | low | `w3` uses `1 − exp(−x)` instead of `expm1` (rel. error ~1e-9 near x = 5e-4); Bezier Ψ* ignores dS/dS_k in the derivative | `w3.c`, `bezier_*.c` | code reading | open |

### 3.3 Memory

Estimated for a realistic 3-D run: 256×256×200 grid (one array V = 105 MB),
24 rays, 5 atoms with ΣNlevel = 190, ~1000 transitions, Nspect = 4000.

| ID | Sev | Item | Formula | 3-D size | Where | Status |
|---|---|---|---|---|---|---|
| [M-1](https://github.com/dwijn/RH/issues/29) | high | Line profiles φ stored per ray (moving atmosphere) | 2·Nrays·ΣNλ·V | ~50 TB (1 TB static) | `profile.c:340-357` | open |
| [M-2](https://github.com/dwijn/RH/issues/30) | high | Γ and collision matrix C stored dense per point | 2·ΣNlevel²·V | 1.7 TB | `initial_xdr.c:567`, `readatom.c:569` | open |
| [M-3](https://github.com/dwijn/RH/issues/31) | high | Per-thread copies of Γ and Rij/Rji (WIP lock fix) | Nthreads·(ΣNlevel² + 2Ntrans)·V | ~1 TB per thread | `opacity.c:571`, `readatom.c:323,423` | open |
| [M-4](https://github.com/dwijn/RH/issues/32) | high | PRD redistribution weight files per point (also a disk problem: they filled /tmp) | ~149 KB/point per line | ~2 TB per PRD line; 1.8 GB per run on a 63×64 grid | `scatter.c` | open |
| [M-5](https://github.com/dwijn/RH/issues/33) | medium | J (and Jgas for hybrid PRD) | Nspect·V each | 420 GB each | `initial_xdr.c:88,103` | open |
| [M-6](https://github.com/dwijn/RH/issues/34) | medium | Hybrid-PRD tables; O(Nspect²) setup | — | > 10 TB | `initial_xdr.c:120-300` | open |
| [M-7](https://github.com/dwijn/RH/issues/35) | medium | Background opacity file | ≥ 3·Nspect·V | ≥ 1.3 TB | `background.c` | open |
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
| [P-7](https://github.com/dwijn/RH/issues/44) | medium | Build flags: `-O2` only; gcc 8 does not vectorize at -O2. `-O3 -march=native` gives −14 % | `Makefile` | open |
| [P-8](https://github.com/dwijn/RH/issues/45) | medium | PRD scattering integral: O(Nspect) linear scans in hybrid PRD; weights re-read from disk every PRD iteration | `scatter.c:282-347`, `formal.c:275-336` | open |
| [P-9](https://github.com/dwijn/RH/issues/46) | low | `time.out` uses `clock()` (CPU time summed over threads), so it never shows a parallel speedup | `getcpu.c:61` | open |
| [P-10](https://github.com/dwijn/RH/issues/47) | low | Writing `opacity.out` is serial and large (26 s of 69 s in the Bifrost 2-D run) | `writeopac_xdr.c` | open |
| [P-11](https://github.com/dwijn/RH/issues/48) | — | No MPI / spatial domain decomposition: a 256²×200 run cannot fit on one node | — | open |

### 3.5 Thread safety

| ID | Sev | Issue | Where | Status |
|---|---|---|---|---|
| [T-1](https://github.com/dwijn/RH/issues/49) | medium | Shared global `messageStr` written from threaded code paths | `rhf1d/formal.c:209,228`, `bezier_1D.c`, `piecestokes_1D.c:93`, `opacity.c:312,408`, `readj.c` | open |
| [T-2](https://github.com/dwijn/RH/issues/50) | medium | Function-static lazy initialization (data race): Voigt tables, `BiLinear` hunt indices, wigner, hydrogen, metal | `voigt.c:261,447`, `linear.c:59`, `wigner.c:166`, `hydrogen.c:336,548,693`, `metal.c:160` | open |
| [T-3](https://github.com/dwijn/RH/issues/51) | low | File-static spline state, not reentrant (only called serially today) | `spline.c:22-36`, `expspline.c:25-36` | open |

### 3.6 Code quality and build

| ID | Sev | Issue | Status |
|---|---|---|---|
| [Q-1](https://github.com/dwijn/RH/issues/52) | high | Hand-written makefiles with archive-member rules: not safe under `make -j`; hand-maintained, partly wrong dependencies; four ~90 % identical geometry makefiles | open |
| [Q-2](https://github.com/dwijn/RH/issues/53) | high | Fortran still required by rhsc2d/rhsc3d/rhsphere/tools makefiles (`librh_f90.a`) after the WIP dropped it at top level | fixed |
| [Q-3](https://github.com/dwijn/RH/issues/54) | medium | 354 warnings with `-Wall -Wextra -std=c11 -pedantic` (138 maybe-uninitialized, 140 unused variables, 32 set-but-unused, 30 unused parameters, 8 sign-compare, 5 implicit declarations — pread/pwrite/fileno need `_POSIX_C_SOURCE`) | open |
| [Q-4](https://github.com/dwijn/RH/issues/55) | medium | Pervasive global state: 315 `extern` declarations in .c files; `atmos` referenced 2145 times in 108 files | open |
| [Q-5](https://github.com/dwijn/RH/issues/56) | medium | ~25 % of code in duplicated blocks; rhsc2d 53 %, rhsc3d 52 %, rhf1d 42 %, rhsphere 38 % | open |
| [Q-6](https://github.com/dwijn/RH/issues/57) | medium | Very long functions: 77 over 100 lines, 31 over 200 (`initSolution` 647, `SortLambda` 572, `Background` 558, `readAtom` 543, `CollisionRate` 494) | open |
| [Q-7](https://github.com/dwijn/RH/issues/58) | medium | Error handling: `Error()` calls `exit()` (262 sites), so the library cannot be unit-tested; 395 `sprintf`, 0 `snprintf`; unbounded `%s` in `sscanf` | open |
| [Q-8](https://github.com/dwijn/RH/issues/59) | medium | Memory ownership implicit; 493 allocations vs 310 frees; `matrix_double` has no NULL check | open |
| [Q-9](https://github.com/dwijn/RH/issues/60) | medium | Headers not self-contained; reserved include-guard names; `rh.h` has none and pulls in `<rpc/types.h>` and `<pthread.h>` everywhere | open |
| [Q-10](https://github.com/dwijn/RH/issues/61) | low | Dead or bit-rotted code: `collision_15d.c` (does not compile), `collision_Oslo.c`, `rhsphere/readinput.c`, `rhsphere/solveplane.c`, `tools/constructatom.c` (does not compile upstream either), `tools/adjustbf.c` | open |
| [Q-11](https://github.com/dwijn/RH/issues/62) | low | Five Voigt algorithms plus Fortran variants | open |
| [Q-12](https://github.com/dwijn/RH/issues/63) | low | `tests/make_tiled_2d_atmos.py` uses the wrong array order | open |
| [Q-13](https://github.com/dwijn/RH/issues/64) | low | Python readers depend on `xdrlib` (removed in Python 3.13) | open |
| [Q-14](https://github.com/dwijn/RH/issues/65) | — | No unit tests (Voigt, LU, splines, w3j, formal solvers on analytic slabs), no CI | open |

Output formats are a contract: the Python (`python/rhanalyze`) and IDL readers
unpack the XDR files field by field, and `rhopacity.py` seeks into raw
native-endian files by record. Layout changes must update or version them.

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
   cleaned (Q-3); dead code removed (Q-10); CI with the regression suite.
3. **Unit tests** for Voigt (against a Faddeeva reference), LU, splines,
   formal solvers on analytic slabs — before touching numerics.
4. **Convergence robustness** — Ng safeguards (N-2), better convergence
   measure (N-3), PRD sub-iteration rework (N-4), positivity (N-8), rate-matrix
   scaling (N-7), ne preconditioning (N-6), collisional-radiative switching.
5. **Speed** — hoist angle-independent work (P-1), parallel statistical
   equilibrium (P-2), 3-D interpolation (P-3), dynamic scheduling (P-4), fused
   accumulation (P-6), flags (P-7).
6. **Memory** — profiles on the fly (M-1), sparse Γ/C by transition and no
   per-thread copies (M-2, M-3), PRD weights tabulated instead of per point
   (M-4), float storage where accuracy allows.
7. **Structure** — shared geometry layer to remove duplication (Q-5),
   self-contained headers, reduced global state, error handling (Q-4, Q-7,
   Q-9). MPI domain decomposition (P-11) is a separate, larger project.

---

## 5. Changelog

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
