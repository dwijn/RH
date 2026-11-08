# Solver redesign (S track)

Design document for replacing RH's iteration scheme (MALI with Ng
acceleration and nested consistency loops) by a solver that treats the
coupled NLTE problem as one nonlinear system. It records why, what is
decided, what is still a hypothesis, how the candidates will be compared,
and the work items (S-1 to S-9) that `docs/MODERNIZATION.md` tracks.

- **Status:** draft, 2026-10-16; revised 2026-11-07 after an audit of this
  document against the code (b5aa20c; every claim below about the code was
  checked by reading it, and the numbers by runs). Nothing here is
  implemented yet.
- **Replaces:** most of plan step 4 in `MODERNIZATION.md` (convergence
  robustness). Steps 5 and 6 (speed, memory) stay; the items in them that
  depend on the solver choice are listed in §8.
- **Labels used below.** **Decided**: agreed, implement as written.
  **Hypothesis**: expected, to be tested in S-5; do not treat as settled.
  **Open**: undecided; the implementer reports options and measurements
  rather than choosing silently.
- **Revision 2026-11-07, decided with the user:** S-2 is split into S-2
  (ne fixed) and S-2b (ne in the state, with `Background()` split into a
  one-time setup and a recompute), delivered in one pull request; the
  background opacity file stays in float by default, with a keyword for
  double, and finite-difference steps come from the measured noise of G
  (§3.7); mixed-precision Newton is recorded as the direction for float and
  GPU work (§3.7, F-4); `1d_falc_helarge` moves to the failing cases of
  S-5 and is judged on the true SE residual only. Also new: the carried
  state of §2.1 confirmed from the code, the work unit of §4.1, two control
  solvers in S-5, and B-38 and B-39 as part of S-2.

For implementers: §2 explains the current scheme and why it fails, §3 the
target architecture (§3.7: precision), §4 the cost model that decides
whether any of this pays, §5 the work items with acceptance criteria, §6 the
test strategy, §7 the traps that make an experiment measure noise instead of
the method. Read §7 before writing code.

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

1. **Ng is already a Krylov-type method, on one block.** For one atom's
   populations with everything else frozen, RH's Ng is Anderson
   acceleration (type II, depth = `NG_ORDER`) applied periodically, solved
   through normal equations without safeguards (§2.2). For linear problems
   Anderson with full history equals GMRES (Walker & Ni 2011). Little
   asymptotic speed is left to gain on easy problems; robustness is the
   weak point.
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

1. `initGammaAtom`/`initGammaMolecule`: collisional and fixed rates into
   Γ; `zeroGammaThreads`/`zeroGammaThreadsMolecule` clear the per-thread
   copies (`N_THREADS` > 1). `initGammaMolecule` resets only the Δv = 1
   entries (B-38).
2. `solveSpectrum(eval_operator = TRUE)`: formal solution over all
   wavelengths and rays; accumulates the radiative rates (`addtoRates`,
   merged inside `solveSpectrum`) and the Ψ*-weighted operator terms into
   Γ (`addtoGamma`; since B-6 `addtoCoupling` only fills the per-wavelength
   emissivities `Uji_down` used by `addPreconditionedRates`). Then
   `mergeGammaThreads`/`mergeGammaThreadsMolecule` add the per-thread Γ
   copies in thread order.
3. `updatePopulations()` (`statequil.c`): per active atom, `statEquil()`
   (LU solve per spatial point; with `I_SUM = −1` the eliminated row is the
   most populated level of the input, chosen per point), then
   `Accelerate()` (Ng on that atom's population vector), then
   `MaxChange()`. Molecules likewise (`NG_MOLECULES`).
4. `Redistribute()` if PRD lines are active: up to `PRD_N_MAX_ITER`
   sub-iterations of `PRDScatter`/`PRDAngleScatter`, each followed by
   `Accelerate()` on ρ (PRD Ng) and `solveSpectrum(eval_operator = FALSE,
   redistribute = TRUE)`, which solves only the wavelengths that contain a
   PRD line (J there, rates of PRD lines only), with populations fixed.
   Early exit when the ρ change falls below `PRD_ITER_LIMIT`; a negative
   value means MAX(dpopsmax, |`PRD_ITER_LIMIT`|), so only
   `PRD_ITER_LIMIT = 0` gives a fixed count. In hybrid PRD
   (`PRD_ANGLE_APPROX`) the sub-iterations re-accumulate `Jgas` from PRD
   wavelengths only, so after the first one the line-edge rows lack the
   contributions of neighbouring wavelengths (N-27).
5. Convergence test on the population `MaxChange` (relative change between
   consecutive `statEquil()` outputs: it excludes ρ, J, ne and the step-6
   rescale). Skipped during collisional relaxation, which also
   re-initializes the population and molecule Ng histories; a `STOP_RH`
   file ends the loop.
6. `Background()` if `SOLVE_NE = ITERATION` (new ne, chemical equilibrium,
   LTE quantities, collision rates, background opacities; skipped after the
   converging iteration); `Hydrostatic()` if hydrostatic equilibrium is on.

G in §3 and §5 is steps 1–4, plus step 6 when ne is in the state (S-2b);
step 5, the Ng calls and the relaxation schedule are the loop around it.
The state at entry to `Iterate()` (starting populations, J after
`initScatter()`) is x₀. After `Iterate()` the output spectrum and rates come
from a final loop of 1 to `N_MAX_SCATTER` formal solutions without operator,
stopped when ΔJ ≤ `ITER_LIMIT`; every solver keeps it unchanged.

**State carried from one iteration to the next** (everything the next
iteration reads that a previous one wrote; confirmed from the code
2026-11-07). This list matters for §3 and §7:

- level populations `atom->n` (active atoms) and `molecule->nv`;
- `spectrum.J` (and `J20` with background polarization): the background
  scattering source uses the previous J (`Jdag` in `rhf1d/formal.c`,
  `sca_c·Jdag` in S). Each wavelength holds J from its most recent formal
  solution: the last PRD sub-iteration at PRD wavelengths, the main formal
  solution elsewhere. In `J_FILE` with `LIMIT_MEMORY`;
- PRD ρ per PRD line (`line->rho_prd`), angle-averaged or angle-dependent;
- with `SOLVE_NE = ITERATION`: ne, and everything `Background()` derives
  from ne and the active H populations: `nstar` of all atoms, the
  collision and fixed rates `atom->C`, `atom->ntotal`, `nHmin` and the
  molecular densities from the chemical equilibrium, and the background
  opacities in `background.dat`. `ChemicalEquilibrium()` also rescales the
  populations of active atoms in place by n_atomic/ntotal
  (`chemequil.c:322-332`), so step 6 changes x itself, by a factor that
  depends on the previous `ntotal`;
- the atmosphere (with hydrostatic equilibrium; `Hydrostatic()` also
  rescales the H populations in place);
- the Ng histories: `atom->Ng_n` and `molecule->Ng_nv` (reset during
  collisional relaxation since B-3), and `line->Ng_prd`, created at the
  first `Redistribute()` and never reset;
- the collisional relaxation schedule (depends on the iteration number);
- the molecular Γ entries of lines with Δv ≥ 2, which are never reset
  (B-38).

Recomputed before use, not carried: Γ of atoms, the radiative rates,
`Jgas`, the intensity scratch files.

Computed once at setup and never updated, although they depend on ne, also
with `SOLVE_NE = ITERATION`: line profiles (Stark and van der Waals
damping), `Qelast`, the PRD weights. With `SOLVE_NE` off the whole
background is computed once. Two conventions of the current ne coupling
that any reference computed with it embodies:

- `Solve_ne` uses NLTE populations only for elements whose model has
  `NLTEpops` set, which is H when the atmosphere supplies H populations;
  other active atoms contribute electrons in LTE (`solvene.c:163-179`,
  `hydrogen.c:89-100`).
- `ChemicalEquilibrium()` takes an active atom's neutral fraction from the
  LTE `nstar` in a fresh run and from the NLTE populations after
  `OLD_POPULATIONS` (forced by `STARTING_J = OLD_J`), so the map, and
  slightly its fixed point, depend on how the run was started (FALC,
  H_6 + Ca II: 3.7×10⁻⁸ between a fresh and a restarted run, both
  converged to 10⁻⁹). The first `Background()` call also solves ne from
  scratch, ignoring the ne of the atmosphere file and of an earlier run.

### 2.2 RH's Ng is Anderson acceleration

Notation: x_k is the population vector of one atom after iteration k, G the
MALI map (one iteration without acceleration), r_k = G(x_k) − x_k the
fixed-point residual. Without acceleration in between, successive
differences are residuals: x_{k+1} − x_k = r_k.

`Accelerate()` (`accelerate.c`), with m = `NG_ORDER` and c the number of
iterates since the last `NgInit` (x₀ is the vector given to `NgInit`):

1. Δ_i = x_{c−i} − x_{c−1−i} for i = 0…m, i.e. Δ_i = r_{c−1−i}.
2. With D_j = Δ_{j+1} − Δ_0 and weights w_ℓ = 1/|x_{c,ℓ}| per element ℓ
   (x_c = G(x_{c−1}), the newest unaccelerated iterate), it forms
   A_ij = Σ_ℓ w_ℓ D_i,ℓ D_j,ℓ and b_j = −Σ_ℓ w_ℓ D_j,ℓ Δ_0,ℓ and solves
   A β = b (`SolveLinearEq`, LU).
3. These are the normal equations of
   min_β ‖Δ_0 + Σ_j β_j D_j‖_w = min ‖(1 − Σβ) r_{c−1} + Σ_j β_j r_{c−2−j}‖_w,
   a combination of the last m+1 residuals with coefficients summing to 1.
4. The new iterate is x = (1 − Σβ) x_c + Σ_j β_j x_{c−1−j}
   = (1 − Σβ) G(x_{c−1}) + Σ_j β_j G(x_{c−2−j}). It replaces the newest
   stored vector (`accelerate.c:132`); the history is not reset.

That is Anderson acceleration type II with mixing parameter 1 and depth m
(also known as Pulay mixing or DIIS), applied every `NG_PERIOD` iterations;
the first extrapolation comes at iteration max(`NG_DELAY`, `NG_ORDER` + 2)
− 1, counted from the start of `Iterate()` or from the last relaxation
reset. All differences in the window are residuals, and Ng is exactly
Anderson(m), only if `NG_PERIOD` ≥ `NG_ORDER` + 1. That holds for the
shipped `NG_ORDER = 2`, `NG_PERIOD = 3`: each extrapolation uses only the
plain steps since the previous one, a restarted Anderson. With
`NG_PERIOD` ≤ `NG_ORDER` the window contains the previous extrapolation
step and Ng is not Anderson; PRD Ng is in this case with its default
`PRD_NG_PERIOD = 1`. Walker & Ni (2011): for a linear map G(x) = Mx + b,
Anderson without truncation is equivalent to GMRES on (I − M)x = b. Ng is
therefore a short-window, restarted cousin of GMRES.

This identification treats G as a map of one atom's populations. In RH, G
also reads and updates the lagged J, PRD ρ, ne and the background
(`SOLVE_NE`), the atmosphere (HSE) and the other atoms' populations, none of
which Ng extrapolates. Ng is the Anderson formula applied to one block of
the full state, so the GMRES equivalence holds for that block with the rest
frozen. Anderson on the full state of S-2 (S-3) is a different algorithm,
even at m = 2.

Consequences for robustness (each is a candidate cause of N-2 and N-10):

- **Normal equations square the condition number:** κ(A) = κ(D)². When one
  slow mode v with eigenvalue λ ≈ 1 dominates, r_k ≈ c λ^k v, so
  Δ_i ∝ v for all i: the columns D_j become nearly collinear, A becomes
  nearly singular, β is amplified. This happens exactly when acceleration
  is needed.
- **The weighting does not match the convergence test.**
  Σ_ℓ Δ_ℓ²/|x_ℓ| = Σ_ℓ |x_ℓ| (Δ_ℓ/x_ℓ)²: relative changes weighted by the
  population itself, so the fit is dominated by the most populated entries
  (ground levels, deep layers). `MaxChange()` is a maximum relative change
  over all entries; on `1d_falc_h6_ca_mg` its maximum sits at a trace entry
  (n/n_total ~ 10⁻⁹), though restricting it to n/n_total > 10⁻⁶ barely
  changes the value (Mg 8.9×10⁻⁴ → 7.5×10⁻⁴). Ng minimizes one quantity;
  convergence is judged on another. (From code reading; not yet tested.)
- **No safeguards:** no check of the conditioning of A, no check that the
  extrapolated residual is smaller, no positivity, no damping, and no
  guard against zero entries. Floating-point traps are on in all RH
  programs, so a zero population makes the weight 1/|x| raise SIGFPE at
  the first extrapolation. That happens for active molecules with
  `NG_MOLECULES`: `statEquilMolecule` sets the populations to exactly 0
  where the temperature is outside the molecule's range (checked with a
  driver; no stress case has an active molecule).
- **Per atom, independently:** each atom's populations are extrapolated
  with their own coefficients although the atoms are coupled through
  overlapping transitions and background opacity.
- **History not reset** after ne or hydrostatic updates (N-2), so the
  extrapolation mixes iterates of different problems. Worse, with
  `SOLVE_NE = ITERATION` (`ChemicalEquilibrium()`) and with HSE
  (`Hydrostatic()`) the populations are rescaled in place after Ng stored
  them, so the next stored difference contains that rescale, not a
  residual. The population and molecule histories are reset during
  collisional relaxation since B-3; the PRD history is not.
- **PRD Ng** (`PRD_NG_DELAY`/`ORDER`/`PERIOD`, defaults 0/0/1) applies the
  same `Accelerate()` to `line->rho_prd`. Its history is created at the
  first `Redistribute()` and never reset: it spans main iterations, mixing
  ρ iterates from different population states, and its schedule counts
  sub-iterations cumulatively.

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

The sign predicts the period-2 oscillation seen with `SOLVE_NE` (N-6,
baseline 2026-10-15: H_6_CRD oscillating with period 2 at a change of 0.22
after 80 iterations; a 2-cycle on `FALC_shock`). At b5aa20c `FALC_shock`
with H_6_CRD and `SOLVE_NE = ITERATION` alternates with period 2 without
Ng and cycles with period `NG_PERIOD` = 3 with Ng (11.8 after 200
iterations); `FALXCO_80` with H_6_CRD + Ca II and ne converges in 31. The
magnitude, however, is below 1 in this local toy. If the ne loop actually
diverges, the cause must lie in nonlocal coupling through the Lyman and
Balmer radiation field or in the inexact inner solves. **Open:** which one.
Logging the alternating residuals in the failing case (`1d_shock_h6crd_ne`,
S-5) answers it.

PRD has the same structure with a nonlocal coupling: ρ depends on J, J on
the populations through the formal solution, the populations on ρ. N-4
documents the result: limit cycles caused partly by the discontinuous early
exit (32 of 126 Bifrost 1-D columns with angle-averaged PRD). On
`1d_falxco_h6_ca` (b5aa20c) the same atoms in CRD converge in 25
iterations; with PRD, Ng-MALI needs 136, and plain MALI with the early exit
has not converged after 300 (9×10⁻³). With a fixed three sub-iterations
(`PRD_ITER_LIMIT = 0`) plain MALI is below 10⁻³ at iteration 112 and
reaches 7×10⁻¹¹ at 416: the early exit, not the coupling alone, stalls it.

### 2.4 Why the convergence test is not an accuracy test

For a linearly converging iteration with rate ρ, the error after iteration
k is e_k = Σ_{j≥k} Δ_j ≈ Δ_k/(1 − ρ). A criterion "max relative change
< 10⁻³" at ρ = 0.99 leaves ~10⁻¹ error. Ng reduces the effective ρ but makes
Δ_k erratic (it measures Ng jumps). Evidence in N-3 (baseline 2026-10-15,
upstream): `1d_falc_helarge` "converged" to 10⁻³ at 1 and 4 threads differs
by 1.4 %. (Since #230 its Ng trajectory crashes at 1 thread; N-28.)
Newton-type methods stop on the residual of the equations, which is the
quantity that matters.

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
3-D estimate). For the 1-D experiments: full space (S-2: populations, ρ, J;
S-2b adds ne); revisit for multi-D.

**Decided: scaling.** Unknowns are scaled before any Krylov, Anderson or
finite-difference operation: log populations, or populations relative to
the total population of the element at that depth. Populations span more
than ten decades; unscaled vectors make finite differences and inner
products meaningless. The scaling must handle exact zeros (molecular
populations outside their temperature range, §2.2): a floor in log
variables, or zero entries left out.

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
  the converged neighbour column or the previous snapshot. With
  `SOLVE_NE = ITERATION` a restart does not reproduce the converged ne
  today (§2.1); S-2b fixes that for G.
- **Learned initial guesses:** SunnyNet (Chappell & Pereira 2022) maps LTE
  to NLTE populations ~10⁵× faster than a 3-D NLTE code, within 20–40 % of
  the true values when trained on other snapshots of the same simulation.
  Not accurate enough as an answer; possibly a good Newton start.
- **Better starting J** (N-9): `ZERO_RADIATION` lacks photoionization; an
  escape-probability or J = B option.

### 3.7 Precision of G

**What is float today:** data that does not depend on x and is computed
once: line profiles of atoms and molecules (and the `LIMIT_MEMORY` profile
files), the PRD weights, and the background opacity file (since #230;
`Background()` computes in double and rounds once on writing). Everything
that depends on x is double. The background becomes a function of x only
when `Background()` runs inside G: `SOLVE_NE = ITERATION` today, ne in x in
S-2b. Its float rounding is then a step function of ne with steps of
6×10⁻⁸ relative, which reaches G diluted.

**Measured** (2026-11-07, b5aa20c):

- Noise δ of G (ECnoise, §7 trap 3; 1 thread, near convergence,
  max-relative): `1d_falc_h6_ca_mg` 1.5×10⁻¹⁰ (Ca II), 5×10⁻¹¹ (Mg II),
  10⁻¹¹ (H), set by the scattering cancellation of N-30 rather than by
  round-off. With the background recomputed inside G (`1d_falc_solvene`):
  about 4×10⁻⁹ max and 5×10⁻¹⁰ RMS with the float file, about 10⁻¹⁰ with
  it in double.
- Cost: `Background()` is one call and 0.7 % of `1d_falc_h6_ca_mg`; with
  `SOLVE_NE = ITERATION`, 27 calls and 25 % of `1d_falc_solvene`. Float or
  double storage makes no difference to wall time (1-D, 2-D and 3-D cases
  within 9 %, in both directions; `3d_perf_many` 6906 s float, 6682 s
  double) or to the iteration histories (identical to 4 digits). It halves
  the file (`3d_perf_many`: 4.8 vs 9.7 GB).

**Decided (2026-11-07):**

- The background stays in float storage by default. A keyword
  `BACKGROUND_DOUBLE` (default `FALSE`, S-2b) stores it in double, as an
  experiment control.
- Finite-difference steps come from the measured δ of each case
  (ε ≈ sqrt(δ) in scaled variables, S-4), never from the round-off of
  double. That is what makes the float background harmless.
- Computing the background itself in float is not pursued: it is under
  1 % of a run without `SOLVE_NE = ITERATION`, and every float operation
  would add noise where storage adds one rounding.

**Toy model** (`docs/solver/precision_toy.py`): G(x) = Mx + b + c sin x,
N = 200, slowest mode of M 0.995, c = 0.002; Newton–GMRES(40) with
finite-difference products. Lowest true ‖R‖ reached, and evaluations of R:

| Variant | lowest ‖R‖ | evaluations, η = 10⁻³ | η = 10⁻¹ |
|---|---|---|---|
| all double, ε = 1.5×10⁻⁸ | 9×10⁻¹⁵ | 240 | 189 |
| x-dependent term rounded to float ("float background"), ε = 1.5×10⁻⁸ | 5×10⁻¹⁰ | 2072 | 283 |
| same, ε = 10⁻⁴ (from the noise) | 5×10⁻¹⁰ | 235 | 188 |
| R in double, Jacobian products in float, ε = 10⁻³ | 8×10⁻¹⁵ | 6107 | 195 |
| all float, ε = 10⁻³ | 2×10⁻⁶ | 3841 | 155 |

The 5×10⁻¹⁰ of rows 2 and 3 is the offset of the rounded problem's
solution, not noise; in RH the corresponding offset is far below any
criterion. The counts in rows 2, 4 and 5 depend on rounding details (they
change by up to 1.7× with the number of BLAS threads); the lowest ‖R‖ and
the conclusions do not.

**Direction for F-4 and GPUs (recorded 2026-11-07, not scheduled):
mixed-precision Newton.** Evaluate the residual R(x), for the Newton
right-hand side and the stopping test, in double once per Newton step, and
the G calls of the Jacobian products in float. The accuracy of Newton is
set by the residual, and an inexact Newton step needs the Jacobian products
only to the forcing tolerance (Kelley 2022). Row 4 of the toy: double
accuracy at the cost of all-double, provided the forcing stays loose
(η = 10⁻¹ works; η = 10⁻³ asks GMRES for more than float products can give
and costs 25× more than all-double). Anderson has no such split: its iterates are
combinations of G values, so its accuracy is that of the precision G runs
in. Conditions:

- the forcing term bounded below by the error of the float products
  (about 10⁻² in the toy; measure for RH);
- the float G needs the changes found in the 2026-11-06 study
  (`MODERNIZATION.md` §2.4: hybrid I / I − S formal solver, rate terms from
  D = I − S and 1 − Ψχ, shared S and η, scaled intensities, SE solve and
  rate accumulators in double);
- ε from the float G's own δ;
- on GPUs with FP64 at 1/32–1/64 of FP32, one double G per Newton step can
  cost more than its ~10 float products: there the residual must come from
  the float G with those changes (within 6×10⁻⁶ of double in the study,
  inside the S-5 tolerance) or from the CPU.

---

## 4. Cost model

Misfortune is conserved: fewer iterations only help if each one does not
cost correspondingly more. All comparisons in S-5 are in work, not
iteration counts.

### 4.1 Unit of work

**One wavelength solve:** the formal solution at one wavelength, all rays,
in `solveSpectrum()`. A full formal solution (FS) is Nspect of them.
Counted separately: sweeps with and without operator; PRD sub-iterations,
which solve only the wavelengths that contain a PRD line (each costs
N_PRD/Nspect of an FS, and S-2's G pays `PRD_N_MAX_ITER` of them on every
call, where Ng-MALI often exits after one); `initScatter()`; the final
sweep after `Iterate()`. `Background()` and `Hydrostatic()` calls and
their time are counted separately as well: they make no FS call and lie
outside the iteration timer, but `Background()` is 25 % of
`1d_falc_solvene` (§3.7), so FS-only comparisons are invalid for the ne
cases. Report wall time at the run's fixed thread count with every count.
S-1 adds the counters.

Shares within one MALI iteration, from the 2-D PRD profile in
`MODERNIZATION.md` §2.3 (2026-10, before P-7, P-8 and B-6 changed the
costs): formal solution ~36 %, `Opacity` 17 %, `addtoGamma` 8 %,
`addtoRates` 8 %, `addtoCoupling` 5.5 %, PRD scattering integral ~12 %,
rest ~14 % (Voigt setup and other one-time work). Re-profile at the S-5
commit (a 1-D CRD case, a 1-D PRD case, `2d_pert_h6_ca_mg`) before using
them for S-6.

### 4.2 Cost per Krylov iteration

- **Form A** (R = G − x): one Jacobian product is one MALI iteration,
  100 % by definition. Exact bookkeeping; no savings per iteration.
- **Form B** (analytic F): opacities with perturbed populations, one FS,
  plain radiative rates, PRD integral if ρ is eliminated; no Ψ* terms
  (`addtoGamma` + `addtoCoupling` ≈ 13.5 %). Roughly 36 + 17 + 8 + 12 ≈
  73 % of a MALI iteration, plus one preconditioner build per Newton step.
  An estimate from an outdated profile (§4.1); measure it.

Break-even for Form B:
0.73 × (total Krylov iterations) + (Newton steps) × (preconditioner build)
< (Ng-MALI iterations), all in units of one MALI iteration.

Worked example (hypothesis): 4 Newton steps × 10 Krylov iterations = 40
MALI-equivalents in Form A. That beats a case that crawls with Ng
(`1d_falxco_h6_ca`, 136 iterations at b5aa20c) but loses to the
well-behaved references (13–31). The gain is expected where RH now fails
or crawls, not on average.

### 4.3 Memory

GMRES(k) stores k vectors of length N = Nspace·ΣNlevel (more with ρ, J or
ne in the state). For the 3-D estimate in `MODERNIZATION.md` §3.3
(256²×200 points, ΣNlevel = 190): N ≈ 2.5×10⁹, about 20 GB per vector in
double precision, ~400 GB for GMRES(20). Anderson(m) stores about 2m
vectors. RH's Ng keeps `NG_ORDER` + 2 vectors per atom and molecule
(2·`NG_ORDER` + 4 during an extrapolation), plus `PRD_NG_ORDER` + 2 copies
of ρ per PRD line (at least 2, also without PRD Ng; ×2·Nrays with
angle-dependent PRD). Small next to M-1 to M-4 but not negligible, and
much larger if J is part of the state.

### 4.4 Where real savings can come from

- **Inexact Newton:** loose Krylov tolerances in early Newton steps.
- **Cheaper Jacobian products (hypothesis, untested in NLTE):** products
  only need to be approximately right, so a reduced angle set or frequency
  grid, or float arithmetic (§3.7), may serve for the preconditioner or
  early Krylov iterations.
- **Fewer outer iterations from good starts** (§3.6).
- **Cheaper FS for every method:** P-1, P-3 etc. in plan steps 5 and 6.
  These pay off whatever solver wins.

---

## 5. Work items

Each item: goal, scope, interface, traps (see also §7), acceptance. IDs
follow `MODERNIZATION.md` (area **S**, solver). Order: S-1 and S-2 first
(they are prerequisites for any fair comparison), then S-3 and S-4 in
either order, then the S-5 gate. S-6 to S-9 depend on S-5.

### S-1 — Residual norm, reference solutions, metrics

**Goal:** a measuring stick independent of the solver (supersedes N-3).

- Start from `tests/se_residual.py` (B-6), which already checks the fixed
  point: decide its norm and row scaling against the proposal below, add
  the conservation row and the level threshold, and wire it into
  `tests/analyze_run.py`.
- Define the true SE residual F_SE(n): radiative rates from a formal
  solution with the populations fixed (no operator), then Γ(n, J) n per
  point with the conservation row, scaled per row. Background-scattering J
  and PRD ρ: taken from the state when they are unknowns (S-2), with their
  own residual rows (J − Λ[S(J)], ρ − ρ[J]); otherwise iterated to a stated
  tolerance before the rates are formed. State which one `se_residual.py`
  implements (today: the final sweep after `Iterate()`, §2.1, with the
  carried ρ).
- Define the norm. Proposal: RMS over all levels and points of the
  residual scaled by the total element population at that point, plus the
  maximum relative change restricted to levels with n/n_total above a
  threshold, reported at 10⁻¹⁰, 10⁻⁶ and 10⁻³ (on `1d_falc_h6_ca_mg` the
  binding entries have n/n_total ~ 10⁻⁹, §2.2). **Open:** final choice;
  report all while deciding.
- Work counters of §4.1, in `solveSpectrum()` and around `Background()`
  and `Hydrostatic()`.
- Reference solutions for the S-5 case set: 1 thread, no Ng
  (`NG_ORDER = 0`, `PRD_NG_ORDER = 0`), fixed PRD sub-iterations
  (`PRD_ITER_LIMIT = 0`, as in S-2's G), collisional relaxation where
  needed, tight limit (max change ≤ 10⁻¹⁰, or ‖F_SE‖ at the precision
  floor of the direct-form formal solution, N-30), stored with the build
  hash, the build options (`RH_NATIVE`, `RH_O3`, compiler) and the inputs.
  With this recipe `1d_falc_h6_ca_mg`, `1d_falxco_h6_ca` and `1d_shock_h6`
  reach max change ≤ 10⁻¹⁰ (SE residual ≤ 3×10⁻¹⁰). `1d_t5750_h6`
  (2–5×10⁻³ after 500 iterations) and `1d_falc_helarge` (stalls at
  ~2×10⁻³, N-28) do not: they have no reference, and S-5 judges them on
  ‖F_SE‖ alone. A solution from a new solver can become their reference
  once its ‖F_SE‖ reaches the level of the other references.
- `tests/analyze_run.py` reports: the work counts, wall time, final
  ‖F_SE‖, max relative population error against the reference (with the
  same level threshold; the thresholded comparison exists already), and the
  emergent spectrum error against the reference, both per wavelength
  (normalized by the maximum |I| over rays and columns at that wavelength,
  as `tests/compare_spectra.py` does) and per line (normalized by the line's
  peak over its wavelength window). **Open:** which one the S-5 criterion
  uses; report both.
- Add the CRD + ne case of S-5 (`1d_shock_h6crd_ne`) to
  `tests/run_stress.sh`.

**Acceptance:** for the current Ng solver on the S-5 case set the report
shows all metrics, and it shows what the 10⁻³ criterion leaves (N-3):
‖F_SE‖ and the population errors of the Ng runs against the references,
including a case that Ng converges at two thread counts (picked in S-1).

### S-2 and S-2b — A deterministic, smooth map G and an explicit state vector

**Goal:** make one MALI step a well-defined function x ↦ G(x), so that
Anderson and finite-difference Jacobians see a function, not noise or a
hidden history.

**Decided (2026-11-07): two stages, one pull request.** S-2 has ne fixed
(`SOLVE_NE` off or `ONCE`): x = populations of active atoms and molecules,
ρ of PRD lines, background-scattering J. S-2b adds ne to x, which needs a
different `Background()`. Each stage has its own acceptance tests; S-2's
pass first.

**S-2 (ne fixed):**

- Everything else in §2.1 is either computed once (unchanged) or reset at
  the start of every evaluation; document which.
- `G_eval(x_in, x_out)`: copies x_in into the global state, runs steps 1–4
  of §2.1 **without** any Ng (bypass `Accelerate()`, or snapshot and
  restore `Ng_n`, `Ng_nv` and `Ng_prd`, count and buffers), copies the
  result to x_out, and leaves no other state changed. A change metric it
  reports is computed from x_in and x_out with `MaxChange`'s definition
  (|new − old|/|new|, zeros skipped).
- `initGammaMolecule` resets all of the molecular Γ (fixes B-38).
- PRD inside G: `PRD_ITER_LIMIT = 0` exactly (negative values mean
  MAX(dpopsmax, |value|)), `PRD_N_MAX_ITER` ≥ 1 sub-iterations on every
  call, no PRD Ng. The early exit makes G discontinuous (N-4). The cost:
  on `1d_falc_h6_ca_mg` with 3 and 10⁻² the current code used 3
  sub-iterations in 18 of 31 iterations, 2 in 2 and 1 in 11.
- Determinism: at a fixed thread count G is bit-reproducible in `rhf1d`,
  `rhsc2d` and `rhsc3d` (per-thread Γ and rate copies merged in thread
  order, static split of the wavelengths). Two exceptions sum under locks
  in arrival order: hybrid PRD (`PRD_ANGLE_APPROX`) in moving atmospheres
  (`Jgas`; 1.7×10⁻¹¹ after 6 iterations) and `rhsphere` (Γ and rates;
  2–6×10⁻¹¹ after 4). Use 1 thread there, or per-thread buffers merged in
  thread order. Another thread count gives a slightly different G (about
  10⁻¹⁰ per step near convergence), not noise: run each experiment at one
  fixed thread count. Bitwise tests are built with `RH_NATIVE = OFF`.
- Collisional relaxation and hydrostatic equilibrium excluded from the
  experiment (G must not depend on the iteration number).

**Acceptance (S-2):**

- `G_eval` called twice on the same x gives bit-identical output;
- history independence: G(x_a), then G(x_b), then G(x_a) again gives the
  first G(x_a) bit for bit (the twice-identical test cannot see state that
  is consistent along the plain iteration path);
- with Ng off, iterating x ← G_eval(x) from x₀ reproduces the current
  code's iteration history bit for bit on `1d_falc_h6_ca_mg` with
  `NG_ORDER = 0`, `PRD_NG_ORDER = 0` and `PRD_ITER_LIMIT = 0`;
- δ measured with ECnoise (§7 trap 3) on a converged state of each S-5
  case that S-2 covers, per atom, in max and RMS norms;
- a finite-difference check of J v at several step sizes on a converged
  state shows the expected O(ε) / O(δ/ε) behaviour, with a fixed `I_SUM`
  and with `S_LINEAR` as a smooth control (§7 trap 3).

**S-2b (ne in x).** Requires `SOLVE_NE = ITERATION`. Step 6 recomputes, for
the next iteration, state that is not in x, and rescales x in place
(§2.1). Measured with an instrumented copy (H_6 + Ca II + Mg II on FALC,
`NG_ORDER = 0`, `PRD_ITER_LIMIT = 0`): after restoring (n, ρ, J, ne), a
second G(x₄) differs from the first by up to 55 % in n (35 % in ne, 77 % in
J); without `SOLVE_NE` the same test gives no differing entry. Copying ne
into `atmos.ne` is therefore not enough:

- Split `Background()` into a one-time setup (line lists, Barklem tables,
  the file and the per-wavelength arrays, the analysis output) and a
  recompute from a given ne and H populations: LTE quantities (`nstar`,
  `atom->C` with the fixed rates), chemical equilibrium (`ntotal`,
  molecular densities, `nHmin`) and the background opacities. The split
  fixes B-39 and the `Background()` part of B-11: G_eval runs the
  recompute thousands of times.
- `G_eval(x)` recomputes these from x.ne and x's H populations **before**
  step 1, so steps 1–4 see quantities that are functions of x only; step
  6's `Solve_ne` then gives ne_out. No in-place rescale by the previous
  `ntotal`: the conservation row of `statEquil()` uses the `ntotal`
  computed from x.ne.
- Inner iterations with fixed counts and fixed starts (§7 trap 1):
  `Solve_ne` (today Newton to 10⁻²) starts from x.ne; the chemical
  equilibrium (today Newton with Ng to 10⁻³, warm-started from the previous
  call) starts from a function of x only, with a fixed count and no Ng.
  Both report their final residual so that the counts can be checked.
- **Open (report):** which populations give the neutral fraction of active
  atoms in the chemical equilibrium. Today LTE in a fresh run and NLTE
  after `OLD_POPULATIONS` (§2.1); G must use one convention whatever the
  start. Recommendation: x's populations.
- Unchanged, listed as part of the model: profiles, `Qelast` and PRD
  weights stay at the startup ne (S-7 decides whether they follow ne);
  charge conservation uses NLTE populations of H only (§2.1).
- Background precision (§3.7): float storage by default; the keyword
  `BACKGROUND_DOUBLE` (default `FALSE`) stores the file in double
  (`backgrbytes = sizeof(double)`, `backgretaexp = 0`, which the reader
  already accepts). The finite-difference step comes from the measured δ,
  which the float file raises for these cases.
- Cost: one recompute per G_eval, counted (§4.1).

**Acceptance (S-2b),** on `1d_falc_solvene` and `1d_shock_h6crd_ne`:

- twice-identical and history-independence tests as for S-2, bit for bit,
  run in-process on states the solver reached (a restart of a `SOLVE_NE`
  run does not restore its ne, §2.1);
- the number of open descriptors and the memory stay flat over at least
  2000 calls of G_eval (more than the common descriptor limit of 1024;
  B-39);
- plain iteration of G converges on `1d_falc_solvene`, and its fixed point
  agrees with a tightly converged run of the current code (`NG_ORDER = 0`)
  to the S-5 correctness tolerance, with any larger difference explained
  (the history itself is not reproduced: the order of operations differs);
- δ measured with ECnoise for population and for ne directions
  separately, with the float and with the double background file, and the
  finite-difference check as for S-2.

### S-3 — Stabilized Anderson acceleration (replaces Ng)

**Goal:** a robust drop-in for Ng on G; also the control for S-4
(supersedes N-2, part of N-8).

- Anderson type II on R = G(x) − x, sliding window m (default 5; report
  3, 5, 10), applied every iteration after a delay. `ANDERSON_DELAY`
  counts like `NG_DELAY` (§2.2).
- Least squares by QR (or SVD with truncation at a relative singular value
  of ~10⁻¹⁰), never normal equations. Drop the oldest column when the
  condition estimate exceeds ~10⁸ (Open: thresholds).
- Scaled variables (§3.1, including exact zeros); one coefficient set for
  the whole state vector (all atoms and ρ together), with an option for
  per-atom sets to compare.
- Safeguard: compare ‖R(x_AA)‖, which the next G_eval computes anyway,
  with ‖R(x_k)‖. If it is larger by more than a factor (Open: start with
  1), discard x_AA, take the plain step G(x_k) (already computed) and
  restart the history. A rejection costs one wasted G_eval, which counts in
  the work totals; report the number of rejections (cf. Zhang, O'Donoghue
  & Boyd 2020 for a globally convergent variant).
- Positivity: automatic in log variables; otherwise reject or damp.
- History reset only for changes outside G: relaxation steps, and ne or
  hydrostatic updates when ne or the atmosphere is not in x. With ne in x
  (S-2b), `Background()` is part of G and nothing is reset. Running
  `SOLVE_NE = ITERATION` with ne outside x and resetting every iteration
  disables the acceleration; do not report that as an Anderson result.
- Keyword: `SOLVER = ANDERSON` with `ANDERSON_DEPTH`, `ANDERSON_DELAY`.

**Acceptance:** default runs (`SOLVER = NG`) unchanged bit for bit
(`run_tests.sh --exact` against the base commit); with
`SOLVER = ANDERSON` the S-5 case set runs and reports all S-1 metrics.

### S-4 — Newton–GMRES on the MALI map (Form A)

**Goal:** the minimal-intrusion test of the Newton hypothesis: same G as
S-3, Newton outer loop instead of Anderson.

- Solve R(x) = G(x) − x = 0 with inexact Newton (§3.2), finite-difference
  Jacobian products (one G_eval each), restarted GMRES(k) (default k = 20),
  Eisenstat–Walker forcing, backtracking line search on ‖R‖ (§3.3).
- Step size for finite differences: ε = sqrt(δ)(1 + ‖x‖)/‖v‖ in scaled
  variables, with δ measured per case in S-2 (§7 trap 3). **Decided:**
  never the textbook sqrt(u) of double: on a G with a float background it
  cost up to 9× more evaluations in the toy model of §3.7, while the step
  from the measured noise cost the same as all-double. Measured so far:
  δ ≈ 10⁻¹⁰ (populations, `1d_falc_h6_ca_mg`), ε ≈ 10⁻⁵; δ ≈ 4×10⁻⁹ with
  the float background inside G, ε ≈ 6×10⁻⁵.
- Stop on ‖F_SE‖ (S-1), not on the step size.
- Every numerical parameter (ε and the δ it came from, forcing constants,
  restart length, line-search constants, tolerances) is printed in the run
  log and stated in the pull request; none is buried as a default.
- Log per Newton step: ‖R‖, ‖F_SE‖, number of GMRES iterations, line
  search factor, work counts so far.
- Keyword: `SOLVER = NEWTON_KRYLOV`.

**Acceptance:** default runs unchanged; on the S-5 case set the solver runs
to its stopping criterion or reports failure clearly (exit status bit),
with all S-1 metrics.

### S-5 — Experiment and decision gate

**Goal:** decide by measurement which path to take.

**Case set** (from `tests/run_stress.sh --set conv`; 1-D at 1 thread;
iterations of the current default Ng at b5aa20c, all to be re-measured at
the commit S-5 runs on):

| Role | Case | PRD | ne | Ng at b5aa20c |
|---|---|---|---|---|
| Easy controls | `1d_falc_h6_ca_mg` | yes | — | 31 |
| | `1d_falc_h20` | — | — | 25 |
| | `1d_falc_mgi66` | — | — | 13 |
| Slow | `1d_falxco_h6_ca` | yes | — | 136 (plain MALI fails) |
| Failing | `1d_t5750_h6` | yes | — | not converged (4.5×10⁻³ after 200) |
| | `1d_shock_h6` | yes | — | SIGFPE at iteration 11 (141 without Ng) |
| | `1d_falc_helarge` | — | — | SIGFPE at iteration 13 (note 1) |
| | a PRD limit cycle (note 2) | yes | — | limit cycle |
| ne coupling | `1d_falc_solvene` | yes | yes | 27 |
| | `1d_shock_h6crd_ne` (note 3) | — | yes | period-3 cycle |
| Multi-D (only if 1-D is promising) | `2d_pert_h6_ca_mg` | yes | — | 82 at 8 threads (note 4) |
| RMHD | note 5 | | | |

1. He_large: SIGFPE since #230 (81 iterations before, 93 with the
   background file in double; N-28 comment of 2026-11-07); stalls at
   ~2×10⁻³ without Ng (N-28). No reference (S-1): judged on ‖F_SE‖ only.
2. One of the N-4 Bifrost 1-D columns that cycle with angle-averaged PRD
   (needs `--bifrost`), or `1d_falc_h6_tight` (H_6 to 10⁻⁸, N-4); chosen in
   S-1, before the runs.
3. `FALC_shock`, H_6_CRD, `SOLVE_NE = ITERATION`; added to
   `run_stress.sh` in S-1. Without Ng it alternates with period 2 (N-6).
   The only ne case without PRD.
4. At a fixed 8 threads G is bit-repeatable (angle-averaged PRD); run
   every solver at 8 threads.
5. Development cases: the Bifrost cases in `run_stress.sh` (column 260,
   cuts `s4` and `s21`). Held out (§7 trap 8): columns and cuts not used
   for B-29, N-4, T-4 or P-8, e.g. another row of the snapshot or a MURaM
   cut; named before the runs.

**Solvers and controls:**

- Ng, current default: the user-facing baseline, with the shipped
  `NG_ORDER = 2`, `NG_DELAY = 10`, `NG_PERIOD = 3` in 1-D (`NG_DELAY = 6`
  in 2-D and 3-D) and the PRD early exit of the stress keyword files.
  `SOLVER = NG` keeps the `NG_*` keywords and their defaults.
- plain G: S-2's map iterated without acceleration.
- Ng on G: S-2's map (fixed PRD sub-iterations) with Ng 2/10/3.
- Anderson (S-3) and Newton–GMRES (S-4), on G.

The success criteria compare against Ng on G, and both baselines are
reported. Without this control a gain would mix the effect of the
accelerator with that of the fixed PRD sub-iterations (§2.3: on
`1d_falxco_h6_ca` plain G converges where plain MALI with the early exit
does not).

**Proposed success criteria** (to be confirmed before the runs, not after):

- correctness: converged populations within 10⁻⁴ relative of the
  reference (levels above the S-1 threshold), emergent spectra within
  10⁻⁴ (normalization from S-1); for cases without a reference, ‖F_SE‖
  down to the level the cases with references reach;
- easy controls: no more than 20 % more work than Ng on G;
- slow case: at least 1.5× less work than Ng on G;
- failing cases: converge (correctness as above);
- work includes rejected steps (S-3) and `Background()` calls (§4.1).

**Report:** one table per case (work counts, wall time, ‖F_SE‖, population
error, spectrum error, outcome) in `MODERNIZATION.md` §2, plus the log
excerpts showing residual histories.

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
  (supersedes N-6). Decide whether active atoms other than H enter with
  NLTE populations, and whether profiles, `Qelast` and PRD weights follow
  ne (§2.1).
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
  (4 iterations; the CI gate compares with the base commit, `--exact` for
  refactors, and reports drift from upstream) is unchanged.
- A new solver changes the iteration history by design, so its regression
  test compares **converged** solutions against stored references
  (populations, spectra) within a tolerance, on a small, fast set
  (e.g. `1d_falc_h6_ca_mg`, which has PRD lines, and a CRD case such as
  `1d_falc_h20`), not iteration 4.
- Unit tests (plan step 3) are on analytic problems with known answers:
  constant-property slabs, two-level atom with known solution, LU,
  splines, Voigt. RMHD data is not used for unit tests.

**Tier 2, stress (convergence and performance).**

- `run_stress.sh --set conv` and `--set perf`, extended with the S-1
  metrics and run for each solver.
- RMHD cut-outs (Q-18) exist: `tests/make_bifrost_atmos.py` converts a
  public Bifrost snapshot, which is fetched by hand and not committed; a
  fetch script with checksums is still missing. S-5 note 5 separates
  development and held-out cases.
- Before comparing solvers on an RMHD column, check depth resolution
  (N-14): if linear and Bezier formal solvers disagree strongly, the
  column is under-resolved and differences between solvers may be
  discretization, not convergence.

---

## 7. Traps

Each of these makes an experiment measure noise or a bug instead of the
method. They are easy to miss because the code still runs.

1. **Discontinuous G.** Any tolerance-based early exit inside G (PRD
   sub-iterations, the inner ne Newton with tolerance 10⁻², the chemical
   equilibrium in `Background()` with Ng and tolerance 10⁻³, the
   hydrostatic iteration) makes G piecewise; finite-difference Jacobians
   across a jump are meaningless and Anderson sees a different map each
   step. Use fixed inner iteration counts inside G, and fixed starting
   points: the ne Newton, the chemical equilibrium and the hydrostatic
   iteration start from the previous state today.
2. **Hidden history.** Ng (populations, molecules, and PRD, whose history
   is never reset), lagged J, relaxation schedules, everything
   `Background()` derives from ne, the in-place rescale of the populations,
   the molecular Γ entries of B-38: anything G reads that is not in x and
   not reset makes G depend on the past (§2.1, S-2).
3. **Noise versus finite-difference step.** With relative noise δ in
   G_eval and step ε, the error of a Jacobian product is about
   δ/ε + ε·(curvature); the best ε ≈ sqrt(δ) gives an error ≈ 2 sqrt(δ).
   Measure δ per case with ECnoise (Moré & Wild 2011): about 8 points
   G(x + t_i v) along a direction v in scaled variables, with t_i spaced at
   least as wide as the intended ε, on a converged state, per atom in max
   and RMS norms. Do not use repeated calls (bit-identical at a fixed
   thread count, so δ = 0 even where G rounds data to float) or the
   difference between thread counts (a different G, not noise; it also
   misses the cancellation rounding). Measured (§3.7): at 1 thread δ is
   not round-off but rounding amplified by the scattering cancellation of
   N-30, 1.5×10⁻¹⁰ for Ca II on `1d_falc_h6_ca_mg` (ε ≈ 10⁻⁵, error
   ≈ 3×10⁻⁵), and N-30 suggests more for strongly scattering lines. With
   the float background inside G, about 4×10⁻⁹ (ε ≈ 6×10⁻⁵). Rounding-level
   switches add to it: the eliminated row with `I_SUM = −1` (argmax of n
   per point; use a fixed `I_SUM` for the finite-difference check), the
   pivoting in the LU, and the Steffen limiter of Bezier3 (kinks with
   coefficient jumps ≤ 4×10⁻¹²; `S_LINEAR` is a smooth control).
4. **Unscaled variables.** Populations over ten decades: inner products,
   least squares and finite differences in raw units are dominated by a
   few entries. Scale first (§3.1).
5. **Positivity.** Newton and Anderson steps can produce negative
   populations, which then give negative opacities and SIGFPE (He_large
   dies this way under Ng: `exp(−dt)` overflows in `Bezier3_coeffs`). Log
   variables or step truncation, designed in from the start.
6. **Comparing iterations instead of work.** Report the work counts and
   wall time (§4.1); a Newton "iteration" contains many FS, and an
   iteration with `SOLVE_NE` contains a `Background()` call.
7. **Comparing against an unconverged reference.** References must be far
   tighter than the criteria (S-1); the current 10⁻³ criterion is not
   (N-3), and some cases have no tight reference yet (N-28).
8. **Developing on the test set.** Fix the S-5 criteria before the runs;
   keep held-out RMHD inputs (S-5 note 5) as checks of anything tuned on
   the 1-D set. The Bifrost cases in `run_stress.sh` were used for B-29,
   N-4, T-4 and P-8 and no longer qualify.

---

## 8. Effect on plan steps 5 and 6

Independent of the solver (do them in any case): P-1, P-3, P-4, P-9, P-10,
P-12, P-15, M-5, M-6, M-8. Done since the draft: P-7, P-8, P-13, P-14;
partly done (float storage, #230 and #232): M-1, M-4, M-7.

Depend on the solver choice (defer until S-5 decides):

- **P-2, P-5, P-6, M-2, M-3:** how Γ is built, stored and solved. Γ with
  Ψ* survives as the preconditioner, but with Form B it is built once per
  Newton step rather than every iteration, and its storage and threading
  needs change accordingly.
- **F-4** (float, GPU-capable core): changes the precision of G, which
  sets δ and the accuracy a solver can reach; mixed-precision Newton
  (§3.7) is the recorded direction.
- **New memory item:** Krylov or Anderson history vectors (§4.3); large in
  3-D, especially if J or angle-dependent ρ is in the state vector.

Constraints from S-2 and S-4 on steps 5 and 6, until F-4 changes them:

- x and everything computed from x inside G stay double: populations, J,
  `Jgas`, intensities, ρ, Γ and the rate accumulators; the background is
  computed in double and stored in float (§3.7);
- every accumulation across wavelengths or rays keeps a summation order
  that does not depend on scheduling (P-4, P-5, M-3);
- float is fine for data computed once that does not depend on x: profiles
  (while they are not recomputed from x), PRD weights, hybrid-PRD tables.

Still needed whatever the solver (keep in step 4 under their IDs): N-7
(rate-matrix scaling; the preconditioner needs it), N-9 (starting
solutions), N-12, N-13, N-14 (depth resolution), N-15 (monotonic 3-D
interpolation), N-16, and B-7; B-38 and B-39 come with S-2. N-10 becomes a
test case for the new solvers; its other candidate causes (3-D Ψ*,
interpolation) remain. Inputs to the S track: N-21 (`SolveLinearEq`), N-28
(He_large and CI-VI_Judge), N-30 (the precision floor of the direct-form
formal solution), and the operator candidate of #190 (S-5).

---

## 9. Open questions

1. Full space or reduced space for J, ρ and ne in multi-D, where J in the
   state is expensive (§3.1; the 1-D experiments use the full space).
2. Residual norm, level threshold and spectrum normalization for S-1.
3. Does the ne divergence come from local or nonlocal coupling (§2.3)?
4. Whether Form A's Newton gains survive the per-iteration cost in
   realistic PRD cases (§4.2), and what Form B saves at the current costs
   (re-profile, §4.1).
5. Multilevel CRD preconditioning for PRD (§3.4): untested.
6. Whether reduced angle or frequency sets, or float, work for Jacobian
   products (§4.4): untested in NLTE.
7. RMHD: which held-out inputs (S-5 note 5), and a fetch script.
8. Whether the 1/|x| weighting in RH's Ng (§2.2) is a significant cause of
   its failures; a cheap test is Ng with relative (1/x²) weights, which
   must handle zero entries too.
9. The ne model with ne in x: neutral fractions in the chemical
   equilibrium (S-2b), NLTE electrons from atoms other than H, and
   profiles, `Qelast` and PRD weights that follow ne (S-7).
10. Mixed-precision Newton for F-4 (§3.7): the lower bound of the forcing
    term in RH, and where the residual is computed on a GPU.

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
- Kelley, C. T. 2022, SIAM Review 64, 191: Newton's method in mixed
  precision
- Kelley, C. T., Keyes, D. E. 1998, SIAM J. Numer. Anal. 35, 508:
  convergence of pseudo-transient continuation
- Knoll, D. A., Keyes, D. E. 2004, J. Comput. Phys. 193, 357:
  Jacobian-free Newton–Krylov review ✓
- Lambert, J., Josselin, E., Ryde, N., Faure, A. 2015, A&A 580, A50:
  GMRES-based NLTE for cool stars ✓
- Moré, J. J., Wild, S. M. 2011, SIAM J. Sci. Comput. 33, 1292:
  estimating computational noise (ECnoise)
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
