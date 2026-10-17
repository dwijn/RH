# Solver redesign (S track)

Design document for replacing RH's iteration scheme (MALI with Ng
acceleration and nested consistency loops) by a solver that treats the
coupled NLTE problem as one nonlinear system. It records why, what is
decided, what is still a hypothesis, how the candidates will be compared,
and the work items (S-1 to S-9) that `docs/MODERNIZATION.md` tracks.

- **Status:** draft, 2026-10-16. Nothing here is implemented yet.
- **Replaces:** most of plan step 4 in `MODERNIZATION.md` (convergence
  robustness). Steps 5 and 6 (speed, memory) stay; the items in them that
  depend on the solver choice are listed in §8.
- **Labels used below.** **Decided**: agreed, implement as written.
  **Hypothesis**: expected, to be tested in S-5; do not treat as settled.
  **Open**: undecided; the implementer reports options and measurements
  rather than choosing silently.

For implementers: §2 explains the current scheme and why it fails, §3 the
target architecture, §4 the cost model that decides whether any of this
pays, §5 the work items with acceptance criteria, §6 the test strategy, §7
the traps that make an experiment measure noise instead of the method.
Read §7 before writing code.

---

## 1. Summary

RH solves statistical equilibrium (SE) and radiative transfer by MALI
(Rybicki & Hummer 1991, 1992), accelerated by Ng (1974), with separate,
nested loops for PRD (ρ), the electron density (`SOLVE_NE`) and hydrostatic
equilibrium. The baseline in `MODERNIZATION.md` §2.2 shows where this
breaks: PRD cases that do not converge (N-4), period-2 oscillation with
`SOLVE_NE` (N-6), Ng jumps and crashes (N-2), a convergence test that does
not bound the error (N-3), and multi-D cases that cycle or diverge (N-10).

Three facts frame the redesign:

1. **Ng is already a Krylov-type method.** RH's Ng is Anderson acceleration
   (type II, depth = `NG_ORDER`) applied periodically, solved through
   normal equations without safeguards (§2.2). For linear problems Anderson
   with full history equals GMRES (Walker & Ni 2011). Little asymptotic
   speed is left to gain on easy problems; robustness is the weak point.
2. **The nested loops are a block splitting** of one coupled system. A
   splitting can diverge or cycle while the coupled problem has a perfectly
   good solution (§2.3). No accelerator applied to the outer loop fixes
   this reliably.
3. **Newton–Krylov on the coupled system has been demonstrated** for
   multilevel NLTE in 1-D: Arramy, de la Cruz Rodríguez & Leenaarts (2024)
   needed up to 2× fewer formal solutions than their Rybicki & Hummer (1992)
   implementation, and reached a lower SE residual (FAL-C, six-level Ca II
   and H I, plane-parallel). The paper names PRD and charge conservation
   as extensions, not as done.

**Decided:** build toward a monolithic nonlinear formulation solved by a
globalized inexact Newton–Krylov method with MALI as the preconditioner, in
stages, each gated by measurement (S-5). The existing MALI step is kept: it
becomes the preconditioner and the reference.

**Hypothesis:** the gain is mainly robustness (cases that now fail or crawl
converge) and accuracy (stopping on a true residual); on easy CRD cases
expect at most ~2× fewer formal solutions, possibly none.

---

## 2. The current scheme and why it fails

### 2.1 What one iteration does

From `iterate.c` (`Iterate()`), one main iteration is:

1. `initGammaAtom`/`initGammaMolecule`: collisional part of Γ.
2. `solveSpectrum(eval_operator = TRUE)`: formal solution over all
   wavelengths and rays; accumulates radiative rates and the Ψ*-weighted
   operator terms into Γ (`addtoGamma`, `addtoRates`, `addtoCoupling`).
3. `updatePopulations()` (`statequil.c`): per active atom, `statEquil()`
   (LU solve per spatial point), then `Accelerate()` (Ng on that atom's
   population vector), then `MaxChange()`. Molecules likewise
   (`NG_MOLECULES`).
4. `Redistribute()` if PRD lines are active: up to `PRD_N_MAX_ITER`
   sub-iterations of `PRDScatter`/`PRDAngleScatter`, each followed by
   `Accelerate()` on ρ (PRD Ng) and `solveSpectrum(eval_operator = FALSE)`,
   with populations fixed; early exit when the ρ change falls below
   `PRD_ITER_LIMIT`.
5. Convergence test on the population `MaxChange`.
6. `Background()` if `SOLVE_NE = ITERATION` (new ne, background opacities);
   `Hydrostatic()` if hydrostatic equilibrium is on.

**State carried from one iteration to the next** (everything the next
iteration reads that a previous one wrote). This list matters for §3 and §7:

- level populations `atom->n` (active atoms) and `molecule->nv`;
- `spectrum.J`: the background scattering source uses the previous
  iteration's J (`Jdag` in `rhf1d/formal.c`, `sca_c·Jdag` in S);
- PRD ρ per PRD line (`line->rho_prd`), angle-averaged or angle-dependent;
- ne and the background opacities derived from it (with `SOLVE_NE`);
- the atmosphere (with hydrostatic equilibrium);
- the Ng histories (`atom->Ng_n`, `line->Ng_prd`), and the collisional
  relaxation schedule (depends on the iteration number).

S-2 must confirm this list from the code; it was compiled by reading
`iterate.c`, `statequil.c`, `redistribute.c` and `rhf1d/formal.c` only.

### 2.2 RH's Ng is Anderson acceleration

Notation: x_k is the population vector of one atom after iteration k, G the
MALI map (one iteration without acceleration), r_k = G(x_k) − x_k the
fixed-point residual. Without acceleration in between, successive
differences are residuals: x_{k+1} − x_k = r_k.

`Accelerate()` (`accelerate.c`), with m = `NG_ORDER` and c the current count:

1. Δ_i = x_{c−i} − x_{c−1−i} for i = 0…m, i.e. Δ_i = r_{c−1−i}.
2. With D_j = Δ_{j+1} − Δ_0 and weights w_k = 1/|x_k|, it forms
   A_ij = Σ_k w_k D_i,k D_j,k and b_j = −Σ_k w_k D_j,k Δ_0,k and solves
   A β = b (`SolveLinearEq`, LU).
3. These are the normal equations of
   min_β ‖Δ_0 + Σ_j β_j D_j‖_w = min ‖(1 − Σβ) r_{c−1} + Σ_j β_j r_{c−2−j}‖_w,
   a combination of the last m+1 residuals with coefficients summing to 1.
4. The new iterate is x = (1 − Σβ) x_c + Σ_j β_j x_{c−1−j}
   = (1 − Σβ) G(x_{c−1}) + Σ_j β_j G(x_{c−2−j}).

That is Anderson acceleration type II with mixing parameter 1 and depth m
(also known as Pulay mixing or DIIS), applied every `NG_PERIOD` iterations
after `NG_DELAY`. Walker & Ni (2011): for a linear map G(x) = Mx + b,
Anderson without truncation is equivalent to GMRES on (I − M)x = b. Ng is
therefore a short-window, restarted cousin of GMRES.

Consequences for robustness (each is a candidate cause of N-2 and N-10):

- **Normal equations square the condition number:** κ(A) = κ(D)². When one
  slow mode v with eigenvalue λ ≈ 1 dominates, r_k ≈ c λ^k v, so
  Δ_i ∝ v for all i: the columns D_j become nearly collinear, A becomes
  nearly singular, β is amplified. This happens exactly when acceleration
  is needed.
- **The weighting does not match the convergence test.**
  Σ_k Δ_k²/|x_k| = Σ_k |x_k| (Δ_k/x_k)²: relative changes weighted by the
  population itself, so the fit is dominated by the most populated entries
  (ground levels, deep layers). `MaxChange()` is a maximum relative change,
  dominated by trace levels (N-3). Ng minimizes one quantity; convergence
  is judged on another. (From code reading; not yet tested.)
- **No safeguards:** no check of the conditioning of A, no check that the
  extrapolated residual is smaller, no positivity, no damping.
- **Per atom, independently:** each atom's populations are extrapolated
  with their own coefficients although the atoms are coupled through
  overlapping transitions and background opacity.
- **History not reset** after ne or hydrostatic updates (N-2), so the
  extrapolation mixes iterates of different problems. (It is reset during
  collisional relaxation since B-3.)
- **PRD Ng** (`PRD_NG_ORDER`) is the same algorithm on ρ inside the
  sub-iteration, with the same issues.

### 2.3 Why the nested loops fail: block splitting

Write the coupled problem as F₁(x, y) = 0, F₂(x, y) = 0, with x the
populations and y the extra unknowns (ne, or PRD ρ). Jacobian blocks:
A = ∂F₁/∂x, B = ∂F₁/∂y, C = ∂F₂/∂x, D = ∂F₂/∂y. RH alternates:

1. solve F₁(x_{k+1}, y_k) = 0 → linearized A δx_{k+1} + B δy_k = 0;
2. solve F₂(x_{k+1}, y_{k+1}) = 0 → linearized C δx_{k+1} + D δy_{k+1} = 0;
3. eliminating δx: δy_{k+1} = D⁻¹C A⁻¹B δy_k ≡ T δy_k.

The alternation converges only if the spectral radius ρ(T) < 1, **even
with exact sub-solves**. The coupled problem is solvable whenever the Schur
complement S = D − C A⁻¹B = D(I − T) is nonsingular, i.e. whenever 1 is
not an eigenvalue of T. Between those two conditions (ρ(T) ≥ 1, but no
eigenvalue at 1) a solution exists that the alternation cannot reach and
Newton on the full system can. Negative real eigenvalues of T give
alternating (period-2) behaviour; complex ones give cycles.

**Toy estimate, H ionization ↔ ne** (local, photoionization balancing
recombination; illustrative only):

- proton fraction x with x/(1 − x) = Q/ne gives ∂ln x/∂ln ne = −(1 − x);
- ne = n_H x + n_metals gives ∂ln ne/∂ln x = b, the fraction of electrons
  from hydrogen;
- in log variables T = −(1 − x) b.

The sign predicts the period-2 oscillation seen with `SOLVE_NE` (N-6:
H_6_CRD oscillating with period 2 at a change of 0.22 after 80 iterations;
a 2-cycle on `FALC_shock`). The magnitude, however, is below 1 in this local toy. If the
ne loop actually diverges, the cause must lie in nonlocal coupling through
the Lyman and Balmer radiation field or in the inexact inner solves.
**Open:** which one. Logging the alternating residuals in the failing case
answers it.

PRD has the same structure with a nonlocal coupling: ρ depends on J, J on
the populations through the formal solution, the populations on ρ. N-4
documents the result: limit cycles caused partly by the discontinuous early
exit, and no convergence on `1d_falxco_h6_ca` where the same atoms in CRD
converge in 19 iterations.

### 2.4 Why the convergence test is not an accuracy test

For a linearly converging iteration with rate ρ, the error after iteration
k is e_k = Σ_{j≥k} Δ_j ≈ Δ_k/(1 − ρ). A criterion "max relative change
< 10⁻³" at ρ = 0.99 leaves ~10⁻¹ error. Ng reduces the effective ρ but makes
Δ_k erratic (it measures Ng jumps). Evidence in N-3: `1d_falc_helarge`
"converged" to 10⁻³ at 1 and 4 threads differs by 1.4 %. Newton-type methods
stop on the residual of the equations, which is the quantity that matters.

---

## 3. Target architecture

Layer by layer. Each layer is marked decided, hypothesis or open.

### 3.1 Residual (decided in principle; form open)

One residual for everything that must be consistent:
F(x) = 0 with x = (populations of all active atoms and molecules, and, as
stages are added, ne, PRD ρ or J, ...).

Two forms, used in sequence:

- **Form A, the preconditioned fixed-point residual:** R(x) = G(x) − x,
  where G is one MALI iteration *without* Ng (§2.1, steps 1–4/6). R = 0
  exactly at the solution, and MALI acts as the preconditioner. It needs no
  new physics code; each evaluation costs one MALI iteration. Used for the
  first experiment (S-4).
- **Form B, the analytic SE residual** (Arramy et al. 2024): for given x,
  compute J by one formal solution, then F(x) = Γ(x, J) n with the
  conservation row, per spatial point. Each evaluation skips the Ψ*
  operator terms, so it is cheaper (§4), and the Jacobian is that of the
  physical equations. MALI's Γ with Ψ* becomes an explicit preconditioner.
  Built only if S-5 favours Newton (S-6).

**Open: full space or reduced space.** Hidden state (§2.1) must either be
part of x (full space: ρ, J, ne are unknowns) or be recomputed inside every
evaluation of R or F (reduced space: more work per evaluation, smaller
system, nonlinearity hidden inside F). For background scattering J the full
space is expensive: J has Nspect·Nspace entries against ΣNlevel·Nspace for
the populations (e.g. 4000 vs 190 per point in the `MODERNIZATION.md` §3.3
3-D estimate). Recommended for the 1-D experiments: full space for
everything (simplest and exact); revisit for multi-D.

**Decided: scaling.** Unknowns are scaled before any Krylov, Anderson or
finite-difference operation: log populations, or populations relative to
the total population of the element at that depth. Populations span more
than ten decades; unscaled vectors make finite differences and inner
products meaningless.

### 3.2 Outer solver (decided: inexact Newton–Krylov; details open)

- Newton step: solve J(x_k) δ = −R(x_k) approximately with GMRES.
- Jacobian-vector products by finite differences:
  J v ≈ (R(x + εv) − R(x))/ε. One product = one evaluation of R.
- Inexact Newton: GMRES tolerance η_k‖R_k‖ with the Eisenstat–Walker
  forcing term, e.g. η_k = 0.9 (‖R_k‖/‖R_{k−1}‖)², capped at 0.9 or lower
  (Eisenstat & Walker 1996). Loose far from the solution, tight near it.
- Restarted GMRES(k) to bound memory (§4.3).
- Positivity: in log variables it is automatic; otherwise steps are
  truncated so no population goes non-positive.

### 3.3 Globalization (decided: needed; method open between two)

Newton converges quadratically only near the solution. Far from it (bad
start, steep transition region, strong PRD) it needs globalization:

- **Backtracking line search** on ‖R‖: accept x + λδ for the largest
  λ ∈ {1, 1/2, 1/4, …} with ‖R(x + λδ)‖ ≤ (1 − 10⁻⁴λ)‖R(x)‖. Simple, first
  choice for S-4.
- **Pseudo-transient continuation (PTC)**, physically natural here
  because it follows the time-dependent rate equations:
  solve (I/Δt_k + J) δ = −R(x_k), with Δt grown by switched evolution
  relaxation, Δt_k = Δt_{k−1}‖R_{k−1}‖/‖R_k‖ (Mulder & van Leer 1985;
  theory: Kelley & Keyes 1998). Small Δt: δ ≈ −Δt R, a damped relaxation
  along the physical evolution, robust and nearly positivity preserving.
  Large Δt: Newton. S-8 if the line search is not enough.

### 3.4 Preconditioner (decided: MALI; improvements open)

The Rybicki & Hummer (1992) linearization with local Λ* is an approximate
Jacobian: per spatial point it couples the levels, and it treats the local
part of the radiative coupling implicitly. In Form A it is built in; in
Form B it is an explicit operator, built once per Newton step (Γ with Ψ*,
LU per point) and applied by back-substitution.

For PRD, the demonstrated choice is to precondition with the CRD problem:
Benedusi et al. (2022) solved angle-dependent polarized PRD (two-level,
1-D) with a nested matrix-free GMRES preconditioned by the CRD
description, converging in a few iterations independent of the
discretization. **Open:** extension to multilevel atoms.

### 3.5 Mesh independence (deferred)

With a local Λ* the number of iterations grows with grid refinement
(roughly ∝ N_z² for Jacobi-type ALI, ∝ N_z for Krylov or optimal SOR).
Nonlinear multigrid with Gauss–Seidel smoothing removes the dependence:
Fabiani Bendicho, Trujillo Bueno & Auer (1997) solve multilevel problems in
fewer than 10 formal solutions, O(N) work. This matters for RMHD columns
whose transition region needs fine grids (N-14). Either a multigrid
preconditioner for Newton–Krylov, or FAS with a Newton smoother. S-9,
after S-5.

### 3.6 Initial guesses and continuation (hypothesis: large practical gain)

Newton is only as good as its starting point:

- **Physics continuation:** CRD with fixed ne → add charge conservation
  → add PRD → refine the grid (full multigrid). Each stage starts inside
  the convergence basin of the next.
- **Warm starts:** for RMHD time series and multi-column runs, start from
  the converged neighbour column or the previous snapshot.
- **Learned initial guesses:** SunnyNet (Chappell & Pereira 2022) maps LTE
  to NLTE populations ~10⁵× faster than a 3-D NLTE code, within 20–40 % of
  the true values when trained on other snapshots of the same simulation.
  Not accurate enough as an answer; possibly a good Newton start.
- **Better starting J** (N-9): `ZERO_RADIATION` lacks photoionization; an
  escape-probability or J = B option.

---

## 4. Cost model

Misfortune is conserved: fewer iterations only help if each one does not
cost correspondingly more. All comparisons in S-5 are in work, not
iteration counts.

### 4.1 Unit of work

**One formal solution (FS):** one call of `solveSpectrum()`, all
wavelengths and rays, with or without operator evaluation. A MALI iteration
with PRD costs 1 + (number of PRD sub-iterations) FS. Report FS and wall
time at 1 thread for every run.

Shares within one MALI iteration, from the 2-D PRD profile in
`MODERNIZATION.md` §2.3: formal solution ~36 %, `Opacity` 17 %,
`addtoGamma` 8 %, `addtoRates` 8 %, `addtoCoupling` 5.5 %, PRD scattering
integral ~12 %, rest ~14 % (Voigt setup and other one-time work).

### 4.2 Cost per Krylov iteration

- **Form A** (R = G − x): one Jacobian product is one MALI iteration,
  100 % by definition. Exact bookkeeping; no savings per iteration.
- **Form B** (analytic F): opacities with perturbed populations, one FS,
  plain radiative rates, PRD integral if ρ is eliminated; no Ψ* terms
  (`addtoGamma` + `addtoCoupling` ≈ 13.5 %). Roughly 36 + 17 + 8 + 12 ≈
  73 % of a MALI iteration, plus one preconditioner build per Newton step.
  An estimate from one profile; measure it.

Break-even for Form B:
0.73 × (total Krylov iterations) + (Newton steps) × (preconditioner build)
< (Ng-MALI iterations), all in units of one MALI iteration.

Worked example (hypothesis): 4 Newton steps × 10 Krylov iterations = 40
MALI-equivalents in Form A. That beats `1d_falc_helarge` (108 iterations)
but loses to the well-behaved references (16–26). The gain is expected
where RH now fails or crawls, not on average.

### 4.3 Memory

GMRES(k) stores k vectors of length N = Nspace·ΣNlevel (more with ρ, J or
ne in the state). For the 3-D estimate in `MODERNIZATION.md` §3.3
(256²×200 points, ΣNlevel = 190): N ≈ 2.5×10⁹, about 20 GB per vector in
double precision, ~400 GB for GMRES(20). Anderson(m) stores about 2m
vectors; RH's Ng stores `NG_ORDER` + 2 per atom. Small next to M-1 to M-4
but not negligible, and much larger if J is part of the state.

### 4.4 Where real savings can come from

- **Inexact Newton:** loose Krylov tolerances in early Newton steps.
- **Cheaper Jacobian products (hypothesis, untested in NLTE):** products
  only need to be approximately right, so a reduced angle set or frequency
  grid may serve for the preconditioner or early Krylov iterations.
- **Fewer outer iterations from good starts** (§3.6).
- **Cheaper FS for every method:** P-1, P-3, P-7, M-1 etc. in plan steps 5
  and 6. These pay off whatever solver wins.

---

## 5. Work items

Each item: goal, scope, interface, traps (see also §7), acceptance. IDs
follow `MODERNIZATION.md` (area **S**, solver). Order: S-1 and S-2 first
(they are prerequisites for any fair comparison), then S-3 and S-4 in
either order, then the S-5 gate. S-6 to S-9 depend on S-5.

### S-1 — Residual norm, reference solutions, metrics

**Goal:** a measuring stick independent of the solver (supersedes N-3).

- Define the true SE residual F_SE(n): one formal solution with the
  populations fixed (no operator), radiative rates from that J, then
  Γ(n, J) n per point with the conservation row, scaled per row. Report
  ‖F_SE‖ at the end of every run and, in experiment mode, every iteration.
- Define the norm. Proposal: RMS over all levels and points of the
  residual scaled by the total element population at that point, plus the
  maximum relative change restricted to levels with n/n_total > 10⁻¹⁰.
  **Open:** final choice; report both while deciding.
- Reference solutions for the S-5 case set: 1 thread, no Ng, collisional
  relaxation where needed, tight limit (‖F_SE‖ at round-off level, or
  max change ≤ 10⁻¹⁰), stored with the build hash and inputs.
- `tests/analyze_run.py` reports: FS count, wall time at 1 thread, final
  ‖F_SE‖, max relative population error against the reference (with the
  same level threshold), and emergent spectrum error against the
  reference relative to the line peak.

**Acceptance:** for the current Ng solver on the S-5 case set the report
shows all metrics; the He_large thread discrepancy (N-3) is visible as a
residual that is not small.

### S-2 — A deterministic, smooth map G and an explicit state vector

**Goal:** make one MALI step a well-defined function x ↦ G(x), so that
Anderson and finite-difference Jacobians see a function, not noise or a
hidden history.

- Enumerate all carried state (§2.1) by reading the code, and document it.
  Either include each item in x or reset it deterministically at the start
  of each evaluation of G. Recommended for 1-D experiments: include
  populations, ρ (PRD lines), J (background scattering) and ne.
- Provide `G_eval(x_in, x_out)`: copies x_in into the global state, runs
  steps 1–4 (and 6 when ne is in the state) of §2.1 **without** any Ng,
  copies the result to x_out, and leaves no other state changed. Snapshot
  and restore anything else it touches.
- PRD sub-iterations inside G: fixed count (`PRD_N_MAX_ITER`), no early
  exit on `PRD_ITER_LIMIT`, no PRD Ng. The early exit makes G
  discontinuous (N-4).
- Determinism: experiments at 1 thread, or deterministic reductions
  (fixed summation order regardless of thread count).
- Collisional relaxation and hydrostatic equilibrium excluded from the
  experiment (G must not depend on the iteration number).

**Acceptance:** G_eval called twice on the same x gives bit-identical
output; with Ng off, iterating x ← G_eval(x) reproduces the current code's
iteration history bit for bit on `1d_falc_h6_ca_mg` with
`NG_ORDER = 0` and fixed PRD sub-iterations (`PRD_ITER_LIMIT` small enough
that the current code never exits early); a finite-difference check of
J v at three step sizes on a converged solution shows the expected
O(ε) / O(δ/ε) behaviour (§7).

### S-3 — Stabilized Anderson acceleration (replaces Ng)

**Goal:** a robust drop-in for Ng on G; also the control for S-4
(supersedes N-2, part of N-8).

- Anderson type II on R = G(x) − x, sliding window m (default 5; report
  3, 5, 10), applied every iteration after a delay.
- Least squares by QR (or SVD with truncation at a relative singular value
  of ~10⁻¹⁰), never normal equations. Drop the oldest column when the
  condition estimate exceeds ~10⁸ (Open: thresholds).
- Scaled variables (§3.1); one coefficient set for the whole state vector
  (all atoms together), with an option for per-atom sets to compare.
- Safeguard: accept the accelerated iterate only if the residual of the
  result does not exceed the residual of the plain step by more than a
  factor (e.g. 1); otherwise take the plain step and restart the history
  (cf. Zhang, O'Donoghue & Boyd 2020 for a globally convergent variant).
- Positivity: automatic in log variables; otherwise reject or damp.
- History reset after any change of the problem (ne, hydrostatic,
  relaxation step).
- Keyword: `SOLVER = ANDERSON` with `ANDERSON_DEPTH`, `ANDERSON_DELAY`.

**Acceptance:** default runs (`SOLVER = NG`) unchanged bit for bit; with
`SOLVER = ANDERSON` the S-5 case set runs and reports all S-1 metrics.

### S-4 — Newton–GMRES on the MALI map (Form A)

**Goal:** the minimal-intrusion test of the Newton hypothesis: same G as
S-3, Newton outer loop instead of Anderson.

- Solve R(x) = G(x) − x = 0 with inexact Newton (§3.2), finite-difference
  Jacobian products (one G_eval each), restarted GMRES(k) (default k = 20),
  Eisenstat–Walker forcing, backtracking line search on ‖R‖ (§3.3).
- Step size for finite differences: ε = sqrt(δ)(1 + ‖x‖)/‖v‖ in scaled
  variables, where δ is the relative noise of G_eval (§7). Report δ as
  measured in S-2.
- Stop on ‖F_SE‖ (S-1), not on the step size.
- Every numerical parameter (ε, forcing constants, restart length,
  line-search constants, tolerances) is printed in the run log and stated
  in the pull request; none is buried as a default.
- Log per Newton step: ‖R‖, ‖F_SE‖, number of GMRES iterations, line
  search factor, FS count so far.
- Keyword: `SOLVER = NEWTON_KRYLOV`.

**Acceptance:** default runs unchanged; on the S-5 case set the solver runs
to its stopping criterion or reports failure clearly (exit status bit),
with all S-1 metrics.

### S-5 — Experiment and decision gate

**Goal:** decide by measurement which path to take.

**Case set** (from `tests/run_stress.sh --set conv`, all 1-D, 1 thread):

| Role | Cases |
|---|---|
| Easy controls | `1d_falc_h6_ca_mg`, `1d_falc_h20`, `1d_falc_mgi66` |
| Slow | `1d_falc_helarge` |
| Failing | `1d_falxco_h6_ca` (PRD), `1d_t5750_h6`, `1d_shock_h6` (fails with Ng) |
| ne coupling | `1d_falc_solvene` (converges in 29 since B-1), and the oscillating H_6_CRD + `SOLVE_NE` configuration of N-6 (add it to the stress set if it is not there) |
| Multi-D (only if 1-D is promising) | `2d_pert_h6_ca_mg` |
| RMHD (when available, Q-18) | 1-D columns and a small 2-D slice from an RMHD snapshot (§6) |

**Solvers:** Ng (current default), Anderson (S-3), Newton–GMRES (S-4).

**Proposed success criteria** (to be confirmed before the runs, not after):

- correctness: converged populations within 10⁻⁴ relative of the
  reference (levels above the S-1 threshold), emergent spectra within
  10⁻⁴ of the line peak;
- easy controls: no more than 20 % more FS than Ng;
- `1d_falc_helarge`: at least 1.5× fewer FS than Ng;
- failing cases: converge to the reference within the correctness
  tolerance.

**Report:** one table per case (FS, wall time, ‖F_SE‖, population error,
spectrum error, outcome) in `MODERNIZATION.md` §2, plus the log excerpts
showing residual histories.

**Decision rule** (proposal): if Newton meets the criteria where Anderson
does not, proceed to S-6/S-7. If Anderson meets all of them, make it the
default and keep Newton as an option, revisiting it for PRD and ne. If
neither does, analyse why (§2.3: which coupling dominates) before building
further.

### S-6 — Analytic residual (Form B), conditional on S-5

Implement F(x) = Γ(x, J[x]) n (Arramy et al. 2024) with MALI's Γ + Ψ* as
an explicit preconditioner built once per Newton step. Acceptance: same
answers as S-4 to the S-5 correctness tolerance, and a measured cost per
Krylov iteration (§4.2) below one MALI iteration.

### S-7 — Monolithic ne and PRD, conditional on S-5

- Charge conservation as equations in the residual, ne as unknowns
  (supersedes N-6).
- PRD ρ in the residual (full space) or eliminated (reduced space), with a
  CRD-preconditioned inner solve (§3.4) (supersedes N-4). **Open:** which.

### S-8 — Pseudo-transient continuation and physics continuation

If the line search in S-4 does not rescue the failing cases: PTC (§3.3),
then staged continuation (§3.6).

### S-9 — Mesh-independent convergence (deferred)

Multigrid preconditioner or FAS (§3.5), for fine-grid RMHD columns.

---

## 6. Test strategy

Two tiers, as agreed:

**Tier 1, regression (gate on every pull request).**

- With the default solver (`SOLVER = NG`), the existing `run_tests.sh`
  (4 iterations, compared with upstream) is unchanged.
- A new solver changes the iteration history by design, so its regression
  test compares **converged** solutions against stored references
  (populations, spectra) within a tolerance, on a small, fast set
  (e.g. `1d_falc_h6_ca_mg` and one PRD case), not iteration 4.
- Unit tests (plan step 3) are on analytic problems with known answers:
  constant-property slabs, two-level atom with known solution, LU,
  splines, Voigt. RMHD data is not used for unit tests.

**Tier 2, stress (convergence and performance).**

- `run_stress.sh --set conv` and `--set perf`, extended with the S-1
  metrics and run for each solver.
- RMHD cut-outs (Q-18): 1-D columns and small 2-D slices from a public RMHD
  snapshot (Bifrost or MURaM), chosen to include a steep transition region
  and shocks. **Open:** redistribution terms of the data must be checked
  before anything is committed. If redistribution is not allowed or the
  files are too large, commit a fetch script with checksums and the
  converter instead of the data.
- Before comparing solvers on an RMHD column, check depth resolution
  (N-14): if linear and Bezier formal solvers disagree strongly, the
  column is under-resolved and differences between solvers may be
  discretization, not convergence.

---

## 7. Traps

Each of these makes an experiment measure noise or a bug instead of the
method. They are easy to miss because the code still runs.

1. **Discontinuous G.** Any tolerance-based early exit inside G (PRD
   sub-iterations, inner ne Newton with tolerance 10⁻², hydrostatic
   iteration) makes G piecewise; finite-difference Jacobians across a jump
   are meaningless and Anderson sees a different map each step. Use fixed
   inner iteration counts inside G.
2. **Hidden history.** Ng (populations or PRD), lagged J, relaxation
   schedules: anything G reads that is not in x and not reset makes G
   depend on the past (§2.1, S-2).
3. **Noise versus finite-difference step.** With relative noise δ in
   G_eval and step ε, the error of a Jacobian product is about
   δ/ε + ε·(curvature); the best ε ≈ sqrt(δ) gives an error ≈ 2 sqrt(δ).
   At 1 thread δ is near round-off (ε ~ 10⁻⁸, error ~ 10⁻⁸). With
   thread-order noise δ ≈ 5×10⁻⁹ (`MODERNIZATION.md` §2.1), ε ≈ 7×10⁻⁵ and
   the error ≈ 10⁻⁴: usable for loose Krylov tolerances, but it limits the
   final accuracy. Measure δ (S-2) before choosing ε.
4. **Unscaled variables.** Populations over ten decades: inner products,
   least squares and finite differences in raw units are dominated by a
   few entries. Scale first (§3.1).
5. **Positivity.** Newton and Anderson steps can produce negative
   populations, which then give negative opacities and SIGFPE. Log
   variables or step truncation, designed in from the start.
6. **Comparing iterations instead of work.** Report FS and wall time
   (§4.1); a Newton "iteration" contains many FS.
7. **Comparing against an unconverged reference.** References must be far
   tighter than the criteria (S-1); the current 10⁻³ criterion is not
   (N-3).
8. **Developing on the test set.** Fix the S-5 criteria before the runs;
   keep at least the RMHD cases as held-out checks of anything tuned on
   the 1-D set.

---

## 8. Effect on plan steps 5 and 6

Independent of the solver (do them in any case): P-1, P-3, P-4, P-7, P-8,
P-9, P-10, M-1, M-4, M-5, M-6, M-7, M-8.

Depend on the solver choice (defer until S-5 decides):

- **P-2, P-5, P-6, M-2, M-3:** how Γ is built, stored and solved. Γ with
  Ψ* survives as the preconditioner, but with Form B it is built once per
  Newton step rather than every iteration, and its storage and threading
  needs change accordingly.
- **New memory item:** Krylov or Anderson history vectors (§4.3); large in
  3-D, especially if J or angle-dependent ρ is in the state vector.

Still needed whatever the solver (keep in step 4 under their N IDs): N-7
(rate-matrix scaling; the preconditioner needs it), N-9 (starting
solutions), N-12, N-13, N-14 (depth resolution), N-15 (monotonic 3-D
interpolation), N-16, and B-5 to B-7. N-10 becomes a test case for the new
solvers; its other candidate causes (3-D Ψ*, interpolation) remain.

---

## 9. Open questions

1. Full space or reduced space for J, ρ and ne (§3.1), especially in
   multi-D where J in the state is expensive.
2. Residual norm and level threshold for S-1.
3. Does the ne divergence come from local or nonlocal coupling (§2.3)?
4. Whether Form A's Newton gains survive the per-iteration cost in
   realistic PRD cases (§4.2), and whether Form B's ~27 % saving holds.
5. Multilevel CRD preconditioning for PRD (§3.4): untested.
6. Whether reduced angle or frequency sets work for Jacobian products
   (§4.4): untested in NLTE.
7. RMHD data: which snapshot, redistribution terms, cut-out size.
8. Whether the 1/|x| weighting in RH's Ng (§2.2) is a significant cause of
   its failures; a cheap test is Ng with relative (1/x²) weights.

---

## 10. References

Verified while writing this document are marked ✓; the others are cited
from memory and should be checked before being quoted.

- Arramy, D., de la Cruz Rodríguez, J., Leenaarts, J. 2024, A&A 690, A12:
  Jacobian-free Newton–Krylov for multilevel NLTE ✓
- Benedusi, P., Janett, G., Belluzzi, L., Krause, R. 2021,
  arXiv:2110.11873: Krylov methods and matrix-free implementation for
  polarized transfer ✓
- Benedusi, P., Janett, G., Riva, S., Krause, R., Belluzzi, L. 2022,
  A&A 664, A197: CRD-preconditioned GMRES for angle-dependent PRD ✓
- Benedusi, P., Riva, S., Zulian, P., Štěpán, J., Belluzzi, L. et al.
  2023, J. Comput. Phys. 479: matrix-free 3-D polarized transfer ✓
- Chappell, B., Pereira, T. M. D. 2022, A&A 658, A182: SunnyNet ✓
- Dullemond, C. P., Radiative transfer lecture notes (2017), chapter 4:
  Λ iteration, ALI, Ng, false convergence ✓
- Eisenstat, S. C., Walker, H. F. 1996, SIAM J. Sci. Comput. 17, 16:
  forcing terms for inexact Newton
- Fabiani Bendicho, P., Trujillo Bueno, J., Auer, L. 1997, A&A 324, 161:
  nonlinear multigrid for multilevel transfer ✓
- Hubeny, I., Mihalas, D. 2014, Theory of Stellar Atmospheres, Princeton,
  chapter 13
- Janett, G. et al. 2021, A&A 655, A13: stationary iterative methods for
  polarized transfer (paper I)
- Kelley, C. T., Keyes, D. E. 1998, SIAM J. Numer. Anal. 35, 508:
  convergence of pseudo-transient continuation
- Knoll, D. A., Keyes, D. E. 2004, J. Comput. Phys. 193, 357:
  Jacobian-free Newton–Krylov review ✓
- Lambert, J., Josselin, E., Ryde, N., Faure, A. 2015, A&A 580, A50:
  GMRES-based NLTE for cool stars ✓
- Mulder, W. A., van Leer, B. 1985, J. Comput. Phys. 59, 232: switched
  evolution relaxation
- Ng, K.-C. 1974, J. Chem. Phys. 61, 2680 ✓
- Olson, G. L., Auer, L. H., Buchler, J. R. 1986, JQSRT 35, 431 ✓
- Osborne, C. M. J., Milić, I. 2021, ApJ 917, 14: Lightweaver (RH-style,
  with charge conservation) ✓
- Rybicki, G. B., Hummer, D. G. 1991, A&A 245, 171; 1992, A&A 262, 209 ✓
- Saad, Y., Schultz, M. H. 1986, SIAM J. Sci. Stat. Comput. 7, 856: GMRES ✓
- Socas-Navarro, H., Trujillo Bueno, J. 1997, ApJ 490, 383: equivalence
  of complete linearization and ALI
- Toth, A., Kelley, C. T. 2015, SIAM J. Numer. Anal. 53, 805: convergence
  of Anderson acceleration
- Trujillo Bueno, J., Fabiani Bendicho, P. 1995, ApJ 455, 646:
  Gauss–Seidel and SOR for radiative transfer ✓
- Uitenbroek, H. 2001, ApJ 557, 389: RH ✓
- Walker, H. F., Ni, P. 2011, SIAM J. Numer. Anal. 49, 1715: Anderson
  acceleration, equivalence with GMRES ✓ (result; volume and page from
  memory)
- Zhang, J., O'Donoghue, B., Boyd, S. 2020, SIAM J. Optim. 30, 3170:
  globally convergent type-I Anderson acceleration ✓
