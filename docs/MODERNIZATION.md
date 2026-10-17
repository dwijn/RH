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

---

## 3. Issues

### 3.1 Bugs and crashes

| ID | Sev | Issue | Where | Evidence | Status |
|---|---|---|---|---|---|
| [B-1](https://github.com/dwijn/RH/issues/2) | high | `SOLVE_NE = ITERATION` aborts on the second `Background()` call ("Reached end of datafile"): collision data re-read from an open file without rewinding; `offset_coll` saved but never used; `atom->C` reallocated (leak) | `ltepops.c:239`, `readatom.c:568-569`, `collision.c:482` | crashes upstream and WIP | fixed |
| [B-2](https://github.com/dwijn/RH/issues/3) | high | `PRD_NG_ORDER > 0` with default `PRD_NG_PERIOD = 0` divides by zero | `readinput.c:99`, `accelerate.c:88` | SIGFPE | fixed |
| [B-3](https://github.com/dwijn/RH/issues/4) | high | Collisional relaxation (`COLL_RELAX_NSTEP/FACTOR`) stops after one iteration with LTE populations: convergence test fires during the relaxation phase while collisions are scaled by up to 1e38 | `iterate.c:132`, `fillgamma.c:59` | FALXCO H_6: 1 iteration, delta 1.4e-12, LTE result | fixed |
| [B-4](https://github.com/dwijn/RH/issues/5) | high | `CaII_large` crashes at setup (SIGFPE in `splineCoef` from `SortLambda`) | `sortlambda.c:387` | upstream and WIP | fixed |
| [B-5](https://github.com/dwijn/RH/issues/6) | medium | Linear short-characteristics end point uses stale weights (`w2` not recomputed) | `rhf1d/piecewise_1D.c:126` | flux changes up to 1e-3 when fixed | open |
| [B-6](https://github.com/dwijn/RH/issues/7) | medium | Possible double counting of the cross-coupling term for multi-ion atoms (added once per transition sharing an upper level) | `fillgamma.c:279` | MgI+II converged populations shift 3e-3 at level 70 when changed; correct form not yet established | open — verify against an independent code (e.g. Lightweaver) |
| [B-7](https://github.com/dwijn/RH/issues/8) | medium | 3-D: `Jgas` accumulated with `wmu` inside the to_obs loop, looks double counted | `rhsc3d/formal.c:370-374` | code reading | open |
| [B-8](https://github.com/dwijn/RH/issues/9) | medium | `calcDistance()` falls off the end without a return value | `rhsc3d/geometry.c:265` | -Wreturn-type | fixed |
| [B-9](https://github.com/dwijn/RH/issues/10) | medium | `pread`/`pwrite` results compared as `int` against `size_t`; -1 error mishandled | `readj.c:58,82,146,171,200,229,413,442` | -Wsign-compare | fixed |
| [B-10](https://github.com/dwijn/RH/issues/11) | medium | 32-bit overflow of `Nlambda*Nspace` sizes (`NgInit(int N)`, `Nlamu`) once > 2³¹ (e.g. >163 PRD wavelengths at 13 M points) | `accelerate.c:34`, `redistribute.c` | code reading | open |
| [B-11](https://github.com/dwijn/RH/issues/12) | low | Leaks: `backgrflags`/`backgrrecno` re-malloc'd every `Background()`; `spectrum.I` allocated twice in 3-D with Stokes | `background.c:300,312`, `initial_xdr.c:386,392` | code reading | open |
| [B-12](https://github.com/dwijn/RH/issues/13) | low | AR85-CHP collisions hard-code `atmos.H->n[5]` as protons (wrong unless H has 6 levels; no shipped atom uses it) | `collision.c:883` | code reading | open |
| [B-13](https://github.com/dwijn/RH/issues/14) | low | `ERROR_LEVEL_1` for an unknown solver leaves I and Psi undefined and continues. **Not a bug:** `ERROR_LEVEL_1` terminates the run (see Q-16), and the keyword parser already rejects unknown solvers | `rhf1d/formal.c:209-231` | code reading | wontfix |
| [B-14](https://github.com/dwijn/RH/issues/15) | low | `error.c` `sprintf`s a 512-byte message plus routine name into a 512-byte buffer | `error.c:51` | code reading | open |
| [B-15](https://github.com/dwijn/RH/issues/80) | low | `Atoms/CaII_large.atom` has 13 `OMEGA` collision records equal to `Inf`; since B-4 the run stops at setup with a clear error instead of a SIGFPE | `Atoms/CaII_large.atom` | `1d_falc_cail`: "collision record OMEGA 11 10 has a non-finite value (inf)" | open |
| [B-16](https://github.com/dwijn/RH/issues/93) | medium | Collision data must follow a `TEMP` line, but this is not enforced: a tabulated keyword (OMEGA, CE, CI, CP, CH0, CH+, CH) before any `TEMP` reads an uninitialized item count and a NULL temperature grid (crash); after a fixed-size keyword (AR85-*, BURGESS, SHULL82, BADNELL, SUMMERS) without a new `TEMP` it inherits that keyword's item count, so a line with that many numbers gives silently wrong rates | `collision.c:469` | -Wmaybe-uninitialized triage; no shipped atom file is affected | fixed |
| [B-17](https://github.com/dwijn/RH/issues/94) | low | `writeMULTIatmos()` writes an uninitialized `height` for a GEOMETRIC depth scale (no branch for it); latent, both callers use other scales | `tools/writemulti.c:89` | -Wmaybe-uninitialized triage | fixed |
| [B-18](https://github.com/dwijn/RH/issues/95) | low | `solve1d` reads the boundary-condition codes `btop`/`bbot` from stdin without validating them; an out-of-range code now gives a zero boundary intensity | `rhf1d/solve1d.c` | -Wmaybe-uninitialized triage | open |

### 3.2 Numerics and convergence

| ID | Sev | Issue | Where | Evidence | Status |
|---|---|---|---|---|---|
| [N-1](https://github.com/dwijn/RH/issues/16) | high | No convergence report: at `N_MAX_ITER` the loop ends silently and output is written as if converged; log has no converged/not-converged line; exit status 0 | `iterate.c:132`, `rhf1d/rhf1d.c:91-131` | all runs | open |
| [N-2](https://github.com/dwijn/RH/issues/17) | high | Ng acceleration without safeguards: one global coefficient set for all levels/depths weighted by 1/n; no positivity, residual or conditioning check; history not reset after relaxation, ne or HSE updates | `accelerate.c:88-122` | He_large Ng jump 0.34 → 20.7; FALC_shock H_6 SIGFPE with Ng, converges in 115 without | open |
| [N-3](https://github.com/dwijn/RH/issues/18) | high | Convergence criterion is not an accuracy guarantee: max relative change over all levels/depths, dominated by trace levels, measures Ng jumps, ignores PRD ρ, background J and ne | `maxchange.c:43`, `iterate.c:132` | He_large "converged" to 1e-3 at 1 thread (108 iterations) and 4 threads (48) differ by 1.4 %; MgI+II at 1e-6 vs 1e-9 differ by 2e-5 | open |
| [N-4](https://github.com/dwijn/RH/issues/19) | high | PRD ρ iteration is Λ-iteration with populations fixed; early exit at `PRD_ITER_LIMIT` makes the map discontinuous → limit cycles; ρ never clipped (ρ < 0 possible after Ng) | `redistribute.c:74-105`, `scatter.c:346,765` | H_6 FALC at 1e-8 stalls at 3–9e-4 for 600 iterations; converges in 142 with fixed 3 sub-iterations; FALXCO H_6 (PRD) no convergence, H_6_CRD converges in 17 | open |
| [N-5](https://github.com/dwijn/RH/issues/20) | high | Parabolic short-characteristics Ψ* unbounded (≈ 1 + 1/Δτ_d when Δτ_u ≫ Δτ_d); fallback to linear only on I < 0. `S_PARABOLIC` is the shipped 2-D default | `rhf1d/piecewise_1D.c:234`, `rhsc2d/piecewise_2D.c:309`, `rhsc2d/run/keyword.input:161` | Ψ* up to 12 (> 4000 before blow-up), 11 % of points outside [0,1]; CaII/MgII on FALC_v10 SIGFPE in 3–8 iterations; Bezier/linear converge in 14–28 | fixed |
| [N-6](https://github.com/dwijn/RH/issues/21) | high | No ne preconditioning: `SOLVE_NE` couples ne by plain iteration (dfjk = 0 for NLTE elements); Newton tolerance 1e-2 | `solvene.c:40,165` | (with B-1 patched) H_6_CRD period-2 oscillation at 0.22 after 80 iterations; 2-cycle on FALC_shock | open |
| [N-7](https://github.com/dwijn/RH/issues/22) | medium | Rate-matrix conditioning: LU with row scaling only; iterative improvement in the same precision; LTE relative to the neutral ground level (overflow for high stages); n* clipped to DBL_MIN | `ludcmp.c:73,142`, `ltepops.c:98` | CI–VI Judge (115 levels): 230 negative populations after 4 iterations, b up to 1e307, divergence, SIGFPE | open |
| [N-8](https://github.com/dwijn/RH/issues/23) | medium | No positivity enforcement or damping anywhere | — | grep | open |
| [N-9](https://github.com/dwijn/RH/issues/24) | medium | Starting solutions: `ZERO_RADIATION` has no photoionization (over-recombined start); no J = B or escape-probability option | `initial_xdr.c:21,588-640` | LTE start 40 iterations vs 22 (H_6 FALC) | open |
| [N-10](https://github.com/dwijn/RH/issues/25) | high | Multi-D atmospheres with a horizontally corrugated transition region do not converge: 3-D diverges and crashes, 2-D with the same gradients locks into a period-3 cycle | — | 3-D 8×8×82, H_6+CaII, upstream and WIP. Unperturbed and velocity-only variants converge. With the height/temperature perturbation alone it diverges (SIGFPE) with bicubic or linear 3-D interpolation, Bezier or linear S, PRD or CRD atoms, Ng on or off, and even with 4× gentler gradients (crash at iteration 28). 2-D with equally steep gradients: no crash, deltas cycle 23 → 1 → 2 with period 3 = `NG_PERIOD` for 40 iterations; 4× gentler 2-D converges slowly (1e-2 after 40). Cause not yet isolated: candidates are the 3-D short-characteristics Ψ*, Ng (N-2) and the lack of positivity safeguards (N-8) | open |
| [N-11](https://github.com/dwijn/RH/issues/26) | medium | FPE traps on with no handler: any overflow kills the run with no output or diagnostic | `fpehandler.c:37` | all SIGFPE cases | fixed |
| [N-12](https://github.com/dwijn/RH/issues/27) | low | Collisions with H/protons use H populations frozen at setup; spline interpolation of rates can overshoot negative and clamps outside the table | `collision.c:685,726-746` | code reading | open |
| [N-13](https://github.com/dwijn/RH/issues/28) | low | `w3` uses `1 − exp(−x)` instead of `expm1` (rel. error ~1e-9 near x = 5e-4); Bezier Ψ* ignores dS/dS_k in the derivative | `w3.c`, `bezier_*.c` | code reading | open |
| [N-14](https://github.com/dwijn/RH/issues/81) | medium | Results depend on depth resolution: on FALC_v10 (±10 km/s) the linear and Bezier formal solvers disagree by 5–16 % throughout the chromosphere, so the grid does not resolve the velocity gradients; there is no grid refinement and no resolution check | — | CaII on FALC_v10, linear vs Bezier vs parabolic | open |
| [N-15](https://github.com/dwijn/RH/issues/82) | medium | 3-D horizontal interpolation of S, χ and I (`BICUBIC_3D`, the default) is not limited and can overshoot across a corrugated transition region (the 2-D `Quadr` has a monotonic option); possible contributor to N-10 | `rhsc3d/geometry.c:344,495,633` | code reading | open |

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
| [Q-1](https://github.com/dwijn/RH/issues/52) | high | Hand-written makefiles with archive-member rules: not safe under `make -j`; hand-maintained, partly wrong dependencies; four ~90 % identical geometry makefiles | fixed |
| [Q-2](https://github.com/dwijn/RH/issues/53) | high | Fortran still required by rhsc2d/rhsc3d/rhsphere/tools makefiles (`librh_f90.a`) after the WIP dropped it at top level | fixed |
| [Q-3](https://github.com/dwijn/RH/issues/54) | medium | 354 warnings with `-Wall -Wextra -std=c11 -pedantic` (138 maybe-uninitialized, 140 unused variables, 32 set-but-unused, 30 unused parameters, 8 sign-compare, 5 implicit declarations — pread/pwrite/fileno need `_POSIX_C_SOURCE`) | fixed |
| [Q-4](https://github.com/dwijn/RH/issues/55) | medium | Pervasive global state: 315 `extern` declarations in .c files; `atmos` referenced 2145 times in 108 files | open |
| [Q-5](https://github.com/dwijn/RH/issues/56) | medium | ~25 % of code in duplicated blocks; rhsc2d 53 %, rhsc3d 52 %, rhf1d 42 %, rhsphere 38 % | open |
| [Q-6](https://github.com/dwijn/RH/issues/57) | medium | Very long functions: 77 over 100 lines, 31 over 200 (`initSolution` 647, `SortLambda` 572, `Background` 558, `readAtom` 543, `CollisionRate` 494) | open |
| [Q-7](https://github.com/dwijn/RH/issues/58) | medium | Error handling: `Error()` calls `exit()` (262 sites), so the library cannot be unit-tested; 395 `sprintf`, 0 `snprintf`; unbounded `%s` in `sscanf` | open |
| [Q-8](https://github.com/dwijn/RH/issues/59) | medium | Memory ownership implicit; 493 allocations vs 310 frees; `matrix_double` has no NULL check | open |
| [Q-9](https://github.com/dwijn/RH/issues/60) | medium | Headers not self-contained; reserved include-guard names; `rh.h` has none and pulls in `<rpc/types.h>` and `<pthread.h>` everywhere | open |
| [Q-10](https://github.com/dwijn/RH/issues/61) | low | Dead or bit-rotted code: `collision_15d.c` (does not compile), `collision_Oslo.c`, `rhsphere/readinput.c`, `rhsphere/solveplane.c`, `tools/constructatom.c` (does not compile upstream either), `tools/adjustbf.c` | fixed |
| [Q-11](https://github.com/dwijn/RH/issues/62) | low | Five Voigt algorithms plus Fortran variants | open |
| [Q-12](https://github.com/dwijn/RH/issues/63) | low | `tests/make_tiled_2d_atmos.py` uses the wrong array order | fixed |
| [Q-13](https://github.com/dwijn/RH/issues/64) | low | Python readers depend on `xdrlib` (removed in Python 3.13) | open |
| [Q-14](https://github.com/dwijn/RH/issues/65) | — | No unit tests (Voigt, LU, splines, w3j, formal solvers on analytic slabs), no CI | open |
| [Q-15](https://github.com/dwijn/RH/issues/83) | low | Default `CFLAGS` lack `-g`, so N-11 backtraces resolve to function names only, without file and line | fixed |
| [Q-16](https://github.com/dwijn/RH/issues/84) | medium | `Error(ERROR_LEVEL_1, …)` prints "Trying to continue" but terminates (exit status 2): the continue branch in `error.c` is unreachable. 74 call sites are written to recover afterwards (return FALSE, fall back, skip an output file), so their recovery code is dead, and `ERROR_LEVEL_1` and `_2` differ only in exit status (2 vs 3). Decide per site between a warning with recovery and a fatal error; define the exit statuses together with N-1 | open |
| [Q-17](https://github.com/dwijn/RH/issues/96) | medium | ~30 unchecked I/O results: `fread`/`fwrite`, `xdr_*`, `fscanf`/`fgets` return values assigned and never tested (e.g. a truncated Barklem table goes unnoticed) | fixed |
| [Q-18](https://github.com/dwijn/RH/issues/105) | medium | Multi-D test cases use perturbed FAL-C atmospheres; realistic cases need cut-outs of freely available 3-D RMHD snapshots (e.g. Bifrost, MURaM), small enough for the stress suite and CI, with converter scripts in `tests/` | open |

Output formats are a contract: the Python (`python/rhanalyze`) and IDL readers
unpack the XDR files field by field, and `rhopacity.py` seeks into raw
native-endian files by record. Layout changes must update or version them.


### 3.7 Wishlist

New capabilities, not problems with the existing code.

| ID | Sev | Item | Status |
|---|---|---|---|
| [F-1](https://github.com/dwijn/RH/issues/109) | medium | Tool for building RH model atoms from standard data sources: choosing the right model atom is one of the hard parts of running RH well. Assemble levels, lines, photoionization cross-sections and collisional data from standard sources (NIST, Kurucz, TOPbase/Opacity Project, CHIANTI, collision databases), with options to reduce a large model (level/line selection, superlevels) and checks of completeness. Replaces the removed `constructatom` (Kurucz lines only, levels and lines, Si only) | open |
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
   slabs — before touching numerics. Realistic multi-D test atmospheres
   from RMHD snapshots (Q-18).
4. **Convergence robustness** — Ng safeguards (N-2), better convergence
   measure (N-3), PRD sub-iteration rework (N-4), positivity (N-8), rate-matrix
   scaling (N-7), ne preconditioning (N-6), collisional-radiative switching,
   monotonic 3-D interpolation (N-15) and a depth-resolution test (N-14).
5. **Speed** — hoist angle-independent work (P-1), parallel statistical
   equilibrium (P-2), 3-D interpolation (P-3), dynamic scheduling (P-4), fused
   accumulation (P-6), flags (P-7).
6. **Memory** — profiles on the fly (M-1), sparse Γ/C by transition and no
   per-thread copies (M-2, M-3), PRD weights tabulated instead of per point
   (M-4), float storage where accuracy allows.
7. **Structure** — shared geometry layer to remove duplication (Q-5),
   self-contained headers, reduced global state, error handling (Q-4, Q-7,
   Q-9, Q-16). MPI domain decomposition (P-11) is a separate, larger project.

---

## 5. Changelog

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
