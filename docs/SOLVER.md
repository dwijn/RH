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
  solvers in S-5, and B-38 and B-39 (each fixed in its own pull request
  before S-2, so that the fixes stay if the S track is reverted).

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
- The same noise is a floor for plain iteration (S-1 references,
  2026-11-07): `1d_falc_solvene` without acceleration stalls at a max
  change of 2.3×10⁻⁹ from iteration ~400 to 3000 with the float file, and
  reaches 10⁻¹⁰ in 1001 iterations with it in double. The two end points
  agree to 2×10⁻⁸ in populations and 4×10⁻⁸ in the spectrum: the float
  file limits how far the ne iteration can converge, not the accuracy of
  what it converges to. Reference solutions of ne cases need
  `BACKGROUND_DOUBLE`.
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

**Decided (2026-11-08), revising the above for ne in the state:** with
`SOLVE_NE = ITERATION` and a solver on G (`SOLVER` other than `NG`) the
background is stored in double. The float file makes G a step function
of ne, which limits finite-difference Jacobian products to ~3×10⁻³, and
the measured δ does not show it (S-2b status). The cost is the larger
file only. `BACKGROUND_DOUBLE` may be left out or set to `TRUE` there;
`FALSE` stops the run at input.

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
S-1 added the counters (`work.c`): every run ends with a "Work" block in
the log, which `tests/analyze_run.py summary` parses. The time it reports
for `Iterate()` includes the `Background()` and `Hydrostatic()` calls made
inside its loop. A run that dies (SIGFPE) prints no counts. On
`1d_falc_h6_ca_mg` the PRD sub-iterations solve 48,852 wavelengths against
51,057 of the main sweeps: PRD is about half of the formal-solution work.

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
  Results of the recipe: see the status below. `1d_falc_helarge` and
  `1d_shock_h6crd_ne` have no reference, and S-5 judges them on ‖F_SE‖
  alone. A solution from a new solver can become their reference once its
  ‖F_SE‖ reaches the level of the other references.
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

**Status (2026-11-07, branch `s1/metrics`):** implemented as follows.

- `tests/se_residual.py` implements F_SE as the final sweep after
  `Iterate()`: rates from a formal solution with the populations fixed,
  background J iterated up to `N_MAX_SCATTER` sweeps to `ITER_LIMIT`, PRD ρ
  carried over (not recomputed). It reports all candidate norms: the
  maximum relative imbalance r_c = |F_c| / (n_c Σ_b P_cb) (the relative
  change of one Jacobi step) over all levels and over levels with
  n_c/n_tot ≥ 10⁻¹⁰, 10⁻⁶, 10⁻³; the RMS of r_c n_c/n_tot; and the
  conservation row max |Σn − n_tot|/n_tot (n_tot = Σ n*).
- `tests/analyze_run.py compare` gives population errors with the same
  thresholds and spectrum errors per wavelength and per line (the line's
  window is the Nlambda spectrum wavelengths nearest λ0, since RH does not
  write the line grids). `run_stress.sh --ref DIR` reports them;
  `tests/make_references.sh` runs the recipe and writes a manifest.
- New cases in `run_stress.sh`: `1d_shock_h6crd_ne` and, as the PRD
  limit cycle of S-5 note 2, `1d_bf_c100_prd` (Bifrost column 100, Ca II
  angle-averaged PRD): Ng repeats the deltas 3.4034×10⁻³ and 5.1739×10⁻³
  of N-4 exactly.
- References (build 1453cbe, gcc 8.5, `RH_O3`; in
  `claude-tmp/s1/refs`): max change ≤ 10⁻¹⁰ and ‖F_SE‖ ≤ 3.4×10⁻¹⁰ for
  `1d_falc_h6_ca_mg` (347 iterations), `1d_falc_h20` (547),
  `1d_falc_mgi66` (106), `1d_falxco_h6_ca` (415), `1d_shock_h6` (393),
  `1d_bf_c100_prd` (645) and `1d_t5750_h6` (2527: fixed PRD
  sub-iterations make it converge, as for `1d_falxco_h6_ca`, §2.3).
  `1d_falc_solvene` stalls at 2.3×10⁻⁹ with the float background and
  reaches 10⁻¹⁰ in 1001 iterations with it in double (§3.7); that run is
  its reference. Every `Background()` call leaks a descriptor (B-39), so
  ne runs beyond ~1000 iterations need a raised `ulimit -n` until S-2b.
  `1d_falc_helarge` stalls at 2–3×10⁻³ (N-28) and `1d_shock_h6crd_ne`
  diverges (max change ~10² after 3000 iterations).
- Ng at the defaults (`ITER_LIMIT = 10⁻³`, 1 thread) against the
  references, levels above 10⁻⁶ of the total:

  | Case | iterations, FS | max change | ‖F_SE‖ max | population error | spectrum error (λ / line) |
  |---|---|---|---|---|---|
  | `1d_falc_h6_ca_mg` | 31, 64 | 8.9×10⁻⁴ | 4.3×10⁻⁴ | 6.6×10⁻³ | 1.1×10⁻² / 3.2×10⁻³ |
  | `1d_falc_h20` | 25, 28 | 4.8×10⁻⁴ | 1.1×10⁻⁴ | 3.3×10⁻³ | 3.9×10⁻³ / 3.4×10⁻³ |
  | `1d_falc_mgi66` | 13, 17 | 1.0×10⁻³ | 3.3×10⁻⁴ | 2.4×10⁻⁴ | 2.8×10⁻² / 2.9×10⁻⁴ |
  | `1d_falxco_h6_ca` | 136, 218 | 9.0×10⁻⁴ | 7.2×10⁻⁵ | 4.0×10⁻² | 5.2×10⁻⁴ / 8.9×10⁻⁶ |
  | `1d_falc_solvene` | 27, 54 | 8.3×10⁻⁴ | 1.3×10⁻³ | 7.4×10⁻³ | 1.7×10⁻² / 1.5×10⁻² |
  | `1d_t5750_h6` | 200 (not conv.), 326 | 4.5×10⁻³ | 4.4×10⁻⁴ | 0.84 | 0.83 / 0.11 |
  | `1d_bf_c100_prd` | 200 (cycle), 305 | 5.2×10⁻³ | 1.4×10⁻⁴ | 2.8×10⁻³ | 5.1×10⁻⁴ / 1.0×10⁻⁵ |
  | `1d_shock_h6crd_ne` | 200 (cycle), 204 | 12 | 4.2 | — (no reference) | — |

  `1d_shock_h6` and `1d_falc_helarge` die with SIGFPE (iterations 11 and
  13). FS counts the wavelength solves of all phases in formal solutions.
  At 8 threads the three FAL-C controls give the same iterations and
  metrics, with populations within 10⁻⁹ of the 1-thread runs. The
  spectrum errors were corrected on 2026-11-09 (Q-22): the first version
  of this table grouped the intensities by the wrong axis, which made
  all but one of them too small, the per-wavelength errors by factors
  of 6 to 50 (10⁶ for `1d_t5750_h6`).
- What the 10⁻³ criterion leaves (N-3): population errors 2×10⁻⁴ to
  4×10⁻² on levels above 10⁻⁶ of the total in the runs that count as
  converged, up to 45 times their final max change (`1d_falxco_h6_ca`),
  and 5×10⁻⁴ to 3×10⁻² in the spectrum (largest at Lyα, the Mg II k
  core and 122.6 nm).
  ‖F_SE‖ does not bound the error either: `1d_falxco_h6_ca` has the
  smallest residual and the largest error, as expected where slow modes
  amplify the residual. Hence the S-5 correctness criterion is the error
  against a reference wherever one exists, and ‖F_SE‖ alone only where none
  does (He_large, the CRD + ne shock).

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
- G relies on `initGammaMolecule` resetting all of the molecular Γ
  (B-38, fixed on its own in #251).
- PRD inside G: `PRD_ITER_LIMIT = 0` exactly (negative values mean
  MAX(dpopsmax, |value|)), `PRD_N_MAX_ITER` sub-iterations on every
  call, no PRD Ng. The early exit makes G discontinuous (N-4). The cost:
  on `1d_falc_h6_ca_mg` with 3 and 10⁻² the current code used 3
  sub-iterations in 18 of 31 iterations, 2 in 2 and 1 in 11.
  (`PRD_N_MAX_ITER` < 1 treats every line in CRD, and G then has no PRD
  sub-iterations; B-46.)
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
  molecular densities, `nHmin`) and the background opacities. G_eval runs
  the recompute thousands of times, so it needs B-39 (descriptor leak,
  line lists re-read; fixed on its own in #250) and the `Background()`
  part of B-11 (per-call arrays; its own pull request before S-2b).
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
- Background precision (§3.7): double whenever ne is in the state
  (decided 2026-11-08, below); otherwise float storage by default, and
  the keyword `BACKGROUND_DOUBLE` (default `FALSE`) stores the file in double
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

**Status (2026-11-08, branch `s2/g-map`):** implemented as below; all
acceptance items of S-2 and S-2b pass.

- `SOLVER = NG | FIXED_POINT` (default `NG`, bit-identical to before on
  the regression suite and on `1d_falc_solvene`). `FIXED_POINT` iterates
  x ← G(x) (`solver.c`); `NG_ORDER` on G's output gives plain G (0) or
  Ng on G, the two S-5 controls (decided 2026-11-07). The state vector is
  the populations of active atoms and molecules, ρ of every PRD line, J
  (and J20), and with `SOLVE_NE = ITERATION` ne, each a copy of RH's own
  contiguous array.
- Keywords that G cannot honour stop the run in `readInput()` with one
  warning per conflict saying what to set (decided 2026-11-07): PRD early
  exit or PRD Ng, collisional relaxation,
  hydrostatic equilibrium, `LIMIT_MEMORY`, and with ne in the state
  `OLD_BACKGROUND` or `NON_ICE`.
- Checks: programs `gcheck_1d`, `gcheck_2d`, `gcheck_3d` (`gchecks.c`):
  twice, history (G(x_a), G(x_b), G(x_a)), ECnoise per block of x, FD
  products D(ε) per block, directions restricted to blocks with
  `--perturb` (decided 2026-11-07). ctest `e2e_solver_g` runs the twice and
  history checks (1-D, and 2-D at 2 threads), the bit-for-bit history
  against `NG` and the keyword validation; a mutation (one PRD
  sub-iteration fewer on every other call) fails all three.
- S-2b: `Background()` is split into the equilibria and
  `backgroundOpacities()`; `BackgroundOfState()` resets `ntotal` to its
  initial value (abundance × nHtot) so that nothing from an earlier call
  is used, computes LTE quantities and collision rates, runs
  `ChemicalEquilibriumOfState()` (neutral fractions of active atoms from
  x's populations, decided 2026-11-07; 20 Newton steps from complete
  dissociation; populations of x not rescaled) and the opacities. After
  steps 1–4 `Solve_ne_OfState()` does 10 Newton steps from x.ne. Both end
  at rounding level (last steps 2×10⁻¹⁶). Keyword `BACKGROUND_DOUBLE`.
  The final state's background is recomputed, with analysis output, for
  the final formal solution.

Results (gcheck at 1 thread unless stated; logs in
`claude-tmp/s1/s2/runs`):

- History: `FIXED_POINT` with `NG_ORDER = 0` reproduces `NG` with
  `NG_ORDER = PRD_NG_ORDER = PRD_ITER_LIMIT = 0` bit for bit over 40
  iterations of `1d_falc_h6_ca_mg` (populations, J, rates, spectrum, every
  logged change).
- Twice and history bit-identical on the converged references of all
  eight S-5 cases without ne, on the ne cases (float and double
  background), on the 2-D and 3-D examples at 4 threads, and after 2000
  calls. Hybrid PRD on `FALC_v10` is bit-identical at 1 thread and differs
  by ~10⁻¹⁵ at 4 (`Jgas` summed in arrival order), as expected (S-2
  bullets).
- δ (ECnoise, relative step 10⁻⁹, converged states, max over the block):
  populations 2×10⁻¹³ (Mg I) to 1.8×10⁻¹⁰ (Ca II), ρ 10⁻¹³ to 10⁻¹¹, J
  ≤ 1.5×10⁻¹¹; the 2-D example 10⁻⁹ (J), the 3-D example 4×10⁻¹¹.
  `1d_falc_solvene`: ≤ 2.7×10⁻¹⁰ (populations) and ≤ 5×10⁻¹² (ne) with the
  double background, ≤ 3×10⁻⁹ and ≤ 3.5×10⁻¹⁰ with float (J 3×10⁻⁸).
  `1d_falc_helarge`: 2×10⁻³, independent of the step from 10⁻¹¹ to 10⁻⁷
  and of `I_SUM`, while G is bit-reproducible. **Hypothesis:** G amplifies
  rounding by ~10⁸ for He_large (an ill-conditioned rate matrix, N-7), and
  its stall at 2–3×10⁻³ (N-28) is this floor, not slow convergence.
  **Confirmed (2026-11-08)**: the floor comes from the LU solve of the rate
  equations, and GTH elimination removes it (N-31, next section).
- FD check on `1d_falc_h6_ca_mg`: with `S_LINEAR` every block shows
  D(ε) − D(10ε) falling as O(ε) to ε ≈ 10⁻⁶ (error ~10⁻⁶) and rising as
  δ/ε below; with `S_BEZIER3` ρ of H and J level off at ~10⁻³ between
  ε = 10⁻⁵ and 10⁻⁸ (limiter kinks, trap 3), also with a fixed `I_SUM`.
  `1d_falc_h20`: perturbations of 10⁻¹ invert Rydberg levels (trap 5); at
  ε ≤ 10⁻³ J behaves, the populations never settle (‖Jv‖ ≪ ‖G‖ for
  collision-dominated levels: scale, trap 4).
- S-2b, `1d_falc_solvene`: plain G converges to 10⁻¹⁰ in 859 iterations
  with the double background (the current code: 1001), and stalls at
  9×10⁻⁹ with float (current code: 2.3×10⁻⁹). Its fixed point agrees with
  the reference of the current code to 4.1×10⁻⁸ in the populations and
  6×10⁻⁸ in the spectrum: the neutral-fraction convention (3.7×10⁻⁸,
  §2.1). Float and double fixed points agree to 3×10⁻⁸. Descriptors and
  memory are flat over 2004 calls (`ulimit -n 64`, peak 19.3 MB after 200
  and after 2000). `1d_shock_h6crd_ne`: plain G diverges, as plain MALI
  does.
- FD along ne directions: double background, clean O(ε) to ε = 10⁻⁷
  (error ~10⁻⁶); float background, best ~3×10⁻³ at ε = 10⁻⁴, and along
  population directions ~4×10⁻³ (double 6×10⁻⁵). The float rounding is a
  step function of ne (§3.7), correlated rather than random, so ECnoise
  at a small step (3×10⁻¹⁰ for ne) understates what an FD product sees.
  The decision of §3.7 (float storage by default, FD step from the
  measured δ) assumed δ describes the FD error; for the ne cases it does
  not. **Decided (2026-11-08):** with ne in the state the background is
  stored in double (§3.7); `readInput` sets it, and an explicit
  `BACKGROUND_DOUBLE = FALSE` stops the run. Checked: `1d_falc_solvene`
  with the keyword left out writes the double file and repeats the
  explicit-double run line for line over 30 iterations.

### Rate equations by GTH elimination (N-31, answers §9 Q12)

Measured 2026-11-08 on scratch builds (logs in `claude-tmp/s1/hl`), then
with the change itself.

- The noise of He_large's G comes from the solution of the rate
  equations. RH's LU (`SolveLinearEq`: partial pivoting with row scaling,
  one refinement step in double) was compared with a 60-digit solution
  of the same Γ, dumped at every depth at the N-28 stall state. It is
  wrong by up to 4.4×10⁻³ at levels 26–28 and 36–42: the 4F ³F°, 5D ³D and
  5F ³F° multiplets, whose levels `He_large.atom` splits by ~10⁻⁴ cm⁻¹.
  Collisions couple the levels of a multiplet with rates many orders of
  magnitude larger than the rates out of it. The diagonal −Σ_j Γ_ji is
  dominated by the rates within the multiplet, and its rounding is larger
  than the rates out of the multiplet, which set the populations.
- Scaling cannot help. The componentwise (Skeel) condition number is
  10¹³–10¹⁴ in the deep layers (k = 70–80), and still ~10¹⁴ after row and
  column equilibration. Column scaling by x's populations leaves δ at
  2.2×10⁻³; long-double LU with refinement lowers it to 9×10⁻⁷, the ratio
  of the mantissas.
- Γ is the generator of a Markov chain. GTH elimination (Grassmann,
  Taksar & Heyman 1985) works on the off-diagonal rates only and forms
  every pivot as a sum of rates, so nothing cancels. It is componentwise
  accurate (O'Cinneide 1993), within 2×10⁻¹⁵ of the 60-digit solution at
  every level and depth, at the cost of LU.

| `1d_falc_helarge` | LU | GTH |
|---|---|---|
| δ of the populations (ECnoise, 10⁻⁹) | 1.7×10⁻³ | 3.9×10⁻¹⁴ |
| Reference recipe (no Ng, `ITER_LIMIT` 10⁻¹⁰) | stalls at 2.3×10⁻³ from iteration ~150 | converges in 667 |
| Stress configuration (Ng), `NG_DELAY` 8, 10, 12 at 1 and 2 threads | converges in 3 of 6, outcome depends on the thread count | converges in 2 of 6, the same at 1 and 2 threads |

The two histories of the reference recipe agree to five digits up to
iteration ~100; then LU flattens. The N-28 stall of He_large is this
floor. With Ng the outcome depends on the extrapolation (N-2): after the
first Ng step at `NG_DELAY` 10 the populations give negative radiative
rates out of the He I ground level at k = 29–30, those locations fall
back to LU, and the run stops with SIGFPE at iteration 19 (LU: converges
in 108 at 1 thread, SIGFPE at 11 at 2). GTH makes the outcome
reproducible; it does not make Ng safe (S-3). GTH also lowers δ of the other S-5 references: Ca II 1.3×10⁻¹⁰ →
5.6×10⁻¹³, Mg II 5×10⁻¹¹ → 2×10⁻¹³, PRD Ca II 1.4×10⁻¹⁰ → 4×10⁻¹⁴. H stays
at ~10⁻¹¹, whose δ has another source (N-30).

GTH needs non-negative rates. In RH a rate can be negative for two
reasons. The preconditioned radiative rates are coefficients of a
linearization, not transition probabilities, and can go negative where
Ieff < 0, under a population inversion, or where transitions from the
same upper level overlap. And atomic data can interpolate to negative
values. The preconditioned rates have gone negative only from the states
that Ng extrapolates (He_large, above); in an iteration without Ng the
negative rates found so far are all of the second kind: Ca II's 3d ²D₃/₂–²D₅/₂ collisional rate at 60–71 kK
(B-43, N-32), and in `CI-VI_Judge` the photoionization rate 49 → 60,
whose explicit cross-section, full of resonances, the spline in λ makes
negative. At such locations the rate equations are solved with LU as
before. The log names the first negative rate of each call, with whether
it is collisional or radiative, and the work summary counts the rate
equations solved and how many needed LU, so negative radiative rates
from the preconditioning, which would be of physical interest, show up.

`CI-VI_Judge` (115 levels, six stages; `1d_falc_civi`), without Ng: with
LU the population changes swing between 10¹ and 10³² for 300 iterations.
With GTH the populations converge steadily (largest change 1.2×10⁻² after
40 iterations) everywhere except at k = 37–42, where the negative
photoionization rate 49 → 60 sends the solve to LU. There the populations
below ~10⁻¹⁸ of the total go negative and RH's change never drops below
1. With that cross-section clipped at zero (scratch build), GTH needs no
LU and the run converges to 10⁻³ in 71 iterations, the first time this
case has converged; with Ng it still stops at iteration 12. The
cross-sections need the same treatment as the collision tables (N-32).
With N-32 (monotone interpolation in log-log, no clipping) it converges
in 71 iterations without Ng and needs no LU (2026-11-09).

### S-3 — Stabilized Anderson acceleration (replaces Ng)

**Goal:** a robust drop-in for Ng on G; also the control for S-4
(supersedes N-2, part of N-8).

**Decided (2026-11-09):** the design below. Implementation on `s3/anderson`,
one pull request; the S-5 runs follow in a separate one, so that the code
is reviewed first.

**Algorithm.** Anderson type II (Walker & Ni 2011), mixing 1, on scaled
variables u = s(x):

- Each G evaluation gives (u_k, g_k) with g_k = s(G(x_k)) and the residual
  f_k = g_k − u_k. The history holds the last m differences
  ΔF = [f_{i+1} − f_i], ΔG = [g_{i+1} − g_i], sliding by one every
  iteration (no restarts every few iterations, as Ng has).
- γ = argmin ‖W(f_k − ΔF γ)‖₂ and u_{k+1} = g_k − ΔG γ. One coefficient
  set for the whole state (`ANDERSON_COEFFICIENTS = GLOBAL`);
  `PER_BLOCK` (one least squares per block, as Ng per atom) only for the
  comparison.
- Least squares by QR, recomputed every iteration by classical
  Gram–Schmidt with reorthogonalization in a fixed serial summation order
  (about 4m²N flops, under 1 % of a MALI iteration in 1-D; no LAPACK).
  Never normal equations. The condition number of R is computed exactly
  (‖R‖₁‖R⁻¹‖₁); while it exceeds `ANDERSON_COND_MAX` (default 10⁸, Open)
  the oldest column is dropped. Columns used, condition number and drops
  are logged every iteration.
- Where no extrapolation is made (no history yet, masked entries),
  x_{k+1} = G(x_k) is copied, not computed as s⁻¹(s(G)). So
  `ANDERSON_DEPTH = 0` reproduces `FIXED_POINT` with `NG_ORDER = 0` bit
  for bit.

**Variables and weights.**

- `ANDERSON_VARIABLES = LOG` (default): u = log x for populations, ρ, J
  and ne. f is then the relative change to first order, so Anderson
  minimizes (in ℓ₂) the quantity `MaxChange` and `ITER_LIMIT` judge
  (§2.2: the mismatch of Ng's weights with the convergence test is a
  suspected cause of N-2), and positivity is automatic. Entries where x or
  G(x) is ≤ 0 are masked: zero weight in the least squares, G's value in
  x_{k+1}. J20 is kept linear as J20/J.
- `ANDERSON_VARIABLES = RELATIVE`: linear variables divided by fixed
  totals (populations by the element's n_total at that depth, molecules by
  their total in x₀, ρ by 1, J by J of x₀, ne by n_H,tot). Weights by
  abundance and keeps the conservation sums exactly; positivity is left
  to the guard. Both are built; LOG stays the default unless the
  comparison in S-5 says otherwise.
- Block weights: W gives each block (each atom's and molecule's
  populations, each PRD line's ρ, J, ne) the weight 1/√N_block, so every
  block counts as its mean square. Per-entry weights would let J
  (~150 times the H_6 block on FAL-C) decide alone.

**Schedule.** `ANDERSON_DEPTH` m (default 10 since S-5; 5 before; report
3, 5, 10). History
from the first evaluation; the first extrapolation after evaluation
max(`ANDERSON_DELAY`, 2). The default delay is chosen from a scan of
{2, 5, 10} on the three easy controls and `1d_falxco_h6_ca`, and frozen
before the S-5 runs (§7 trap 8). One iteration is one G evaluation,
rejected ones included, so `N_MAX_ITER` caps the work and the log format
of `FIXED_POINT` is kept. No history resets other than after a rejection:
relaxation, hydrostatic equilibrium and ne outside x are already refused
for every solver other than NG, and with ne in x (S-2b) `Background()` is
part of G.

**Safeguards.**

- Before the evaluation (no G cost): the trial point must be finite,
  positive (RELATIVE), and must not create line inversions that g_k does
  not have: for every line, at every depth where n_l/g_l > n_u/g_u in
  g_k, also in the trial point. Inversions that G itself makes are left
  alone. Continua are not checked (S-13, below): the test on
  n_i/n_i* > (n_j/n_j*)·exp(−hν₀/kT) that was here blocked steps towards
  converged solutions that have inverted continua. On failure the extrapolation is damped
  (u = g_k − θΔGγ, θ = ½, then ¼); if it still fails, the plain step is
  taken and the history restarted. All actions are counted. (Trap 5:
  H_20's Rydberg levels, negative opacity, overflow in Bezier3. The
  inversion test stands in for the real condition, a negative total
  opacity, and since S-13 covers lines only.)
- After the evaluation: accept if ‖f_{k+1}‖_W ≤ c‖f_k‖_W with
  c = `ANDERSON_SAFEGUARD` (default 1; 0 turns it off). Otherwise take the
  plain step G(x_k), already computed, and restart the history from the
  last accepted pair (x_k, g_k). The wasted evaluation counts in the work
  totals; rejections are logged and counted. If c = 1 rejects often on
  the easy controls, the decreasing envelope of Zhang, O'Donoghue & Boyd
  (2020) replaces it (decided with the numbers). At m = 5, S-5 found
  few rejections on the easy controls, 27 of 215 points on `t5750_h6` and
  78 of 282 on `shock_h6`, where Anderson then does worse than plain G
  (S-5 Results); the envelope is now S-12.
- Negative radiative rates (the GTH → LU fallback of N-31) are counted
  per run and reported as an S-3 metric, not used as a guard.

**Stopping and output.** The change of an evaluation is `MaxChange`'s
definition on (x_k, G(x_k)), per atom, as in `FIXED_POINT`; the run stops
below `ITER_LIMIT`. The state left for the final formal solution and the
output is G(x_k) of that evaluation, never an extrapolated point.

**Keywords:** `SOLVER = ANDERSON`, `ANDERSON_DEPTH` (10),
`ANDERSON_DELAY` (from the scan), `ANDERSON_COND_MAX` (10⁸),
`ANDERSON_SAFEGUARD` (1), `ANDERSON_VARIABLES` (`LOG` | `RELATIVE`),
`ANDERSON_COEFFICIENTS` (`GLOBAL` | `PER_BLOCK`), and since S-11
`ANDERSON_EXTRAPOLATE` (`ALL` | `POPULATIONS` | `POPULATIONS_RHO`); all
printed at the start of the run with the memory, (3m + 10)·N doubles and
(m + 4)·N bytes with `ALL` (S-11, item 3, for the others)
(one more N doubles with `RELATIVE`). readInput stops
on: every `FIXED_POINT` conflict (messages name the solver in use);
`NG_ORDER > 0` with `ANDERSON` (two extrapolations); any `ANDERSON_*`
keyword with another solver (it would be ignored silently).

**Code:** `anderson.c` / `anderson.h`, a generic module without RH
globals (weighted, masked vectors; init, step, reset, free, counters),
unit-testable alone; `iterateAnderson()` in `solver.c` (transform per
block, guard, safeguard, log, summary; `G_eval` unchanged); keywords in
`readinput.c`, dispatch in `iterate.c`; the counters are printed at the
end of the iteration and read by `tests/analyze_run.py`.

**Tests:** unit (`test_anderson`): for a linear nonsymmetric map with
full depth, x_{k+1} equals G applied to the GMRES iterate (Walker & Ni) to
10⁻¹⁰; m = 0 equals plain iteration bit for bit; collinear columns are
dropped without NaN and κ ≤ `ANDERSON_COND_MAX` after; a map with a bad
region triggers rejection and restart with the right counts; masked
entries take G's value exactly; two runs are bit-identical. e2e:
`ANDERSON_DEPTH = 0` repeats `FIXED_POINT` with `NG_ORDER = 0` bit for bit
on `1d_falc_h6_ca_mg`; the keyword conflicts stop the run; a short
`ANDERSON` run is bit-reproducible.

**Acceptance:** default runs (`SOLVER = NG`) unchanged bit for bit
(`run_tests.sh --exact` against the base commit); with
`SOLVER = ANDERSON` the S-5 case set runs and reports all S-1 metrics,
plus rejections, guard actions and GTH fallbacks.

**Status (2026-11-09, branch `s3/anderson`):** implemented as above, with
these details fixed in the code:

- Columns left out for the condition number are left out of that solve
  only; the history keeps them until they slide out (at most m
  iterations). With `PER_BLOCK` each block's solve leaves out its own.
- The safeguard compares ‖f‖_W of the new evaluation (over the entries
  valid in it) with that of the newest accepted pair. Convergence is
  tested first: an evaluation below `ITER_LIMIT` ends the run even if
  the safeguard would have rejected its point.
- A run that stops without converging writes G of the last accepted
  point, not of a rejected one.
- `ANDERSON_DELAY` defaults to 2 until the delay scan.
- The guard flags only values that G's point does not have: a negative
  entry where G's is not negative (G can return negative populations
  where statistical equilibrium falls back to LU, N-31), and with
  `RELATIVE` a zero where G's is positive (in `LOG` a zero is an
  underflow of exp).
- With `PER_BLOCK`, a block whose solve uses no column takes G's values,
  copied, while other blocks extrapolate.
- The safeguard norm leaves out entries that are invalid in the new
  evaluation (in `LOG`, values G returns ≤ 0), so a trial point after
  which G returns non-positive values is judged on the rest only.

First run (Ca II with two PRD lines on FAL-C, H_6 passive, the e2e set-up,
`ITER_LIMIT` 10⁻⁶): `ANDERSON` converges in 89 G evaluations (76
extrapolated, 11 rejected), `RELATIVE` in 97; `FIXED_POINT` has not
converged after 200. The rejected points had ‖f‖_W only 0.2–8 % above
the last accepted one, the case where a strict c = 1 can discard useful
steps; S-5 decides with numbers whether the envelope of Zhang,
O'Donoghue & Boyd (2020) replaces it (Safeguards, above; after S-5 this
is S-12).

### S-4 — Newton–GMRES on the MALI map (Form A)

**Goal:** the minimal-intrusion test of the Newton hypothesis: same G as
S-3, Newton outer loop instead of Anderson. Also the scaffold for Form B
(S-6).

**Expectation (2026-11-09, before the runs).** Near the solution the
problem is effectively linear, and there Anderson with unlimited depth is
GMRES on the linearized problem (Walker & Ni 2011; `test_anderson`), as is
one Newton step. Newton discards its Krylov space every Newton step,
where Anderson carries its history over, and pays for line-search trials.
So Form A is expected at best to tie with Anderson on the easy cases, and
warm starts help both alike. Newton can win where the nonlinearity
dominates (the slow tail of `1d_shock_h6` and `1d_t5750_h6`, with many
Anderson rejections), and Form B (cheaper products, preconditioners not
tied to MALI) is where the method could pay off.

**Decided (2026-11-09):** the design below. Implementation on
`s4/newton`, one pull request; the runs on the S-5 set follow in a
separate one.

1. **Unknowns and residual.** u = s(x) in the variables of S-3
   (`NEWTON_VARIABLES = LOG`, default, or `RELATIVE`; J20 as J20/J),
   R(u) = s(G(s⁻¹(u))) − u, the block-weighted norm ‖·‖_W of S-3 for
   GMRES inner products, the forcing test and the line search. Entries
   with no value in the variables (x or G(x) ≤ 0 in `LOG`) are left out of
   the step at that Newton iteration and take G's value, as in Anderson.
2. **Finite-difference products.** J v ≈ (R(u + h v/‖v‖_W) − R(u))·‖v‖_W/h,
   one G evaluation each, h = √δ with δ = `NEWTON_NOISE` (default 10⁻¹⁰,
   measured on `1d_falc_h6_ca_mg`, §7 trap 3), so h = 10⁻⁵. In `LOG`
   variables a step in u is already relative, so the textbook factor
   (1 + ‖x‖) does not apply; never √u of double.
3. **GMRES.** Zero initial guess, at most k = `NEWTON_KRYLOV_MAX` (20)
   Krylov iterations per Newton step, no restart: if the forcing test is
   not met, the best step found is taken. Arnoldi by classical
   Gram–Schmidt with reorthogonalization in a fixed summation order
   (bit-reproducible, as Anderson's QR), Givens rotations. Forcing:
   Eisenstat & Walker (1996) choice 2, η_k = γ(‖R_k‖/‖R_{k−1}‖)^α with
   γ = 0.9, α = 2, their safeguard η_k ≥ γη_{k−1}^α when that exceeds 0.1,
   η_k ≤ η_max = `NEWTON_ETA_MAX` (0.9), η_0 = η_max, and against
   oversolving η_k ≥ ½·`ITER_LIMIT`/‖R_k‖_W (Kelley 1995).
4. **Globalization.** Backtracking on ‖R‖_W: trial u + λδ for
   λ = 1, ½, ¼, ⅛, 1/16, accepted if
   ‖R(u + λδ)‖_W ≤ (1 − 10⁻⁴λ)‖R(u)‖_W. The S-3 guard (finite, no new
   negatives, no new line inversions against G(x_k)) runs before
   each trial and a failure halves λ without a G evaluation. If no λ is
   accepted, the plain step x ← G(x_k) is taken; it costs nothing, as
   G(x_k) is known. The G evaluation of an accepted trial is the base of
   the next Newton step, never paid twice.
5. **Stopping and output.** As every other solver: `MaxChange`'s
   definition on (x_k, G(x_k)) at each base point, per atom, against
   `ITER_LIMIT` (S-5: 10⁻⁶); ‖F_SE‖ is reported by the S-1 metrics, not
   used to stop (this replaces "stop on ‖F_SE‖" in the first version of
   this item, for like-for-like comparisons in S-5). `N_MAX_ITER` caps G
   evaluations (all of them: base, Krylov, line search), so it bounds the
   work as for Anderson. The output is G of the last base point, never a
   probe or a trial.
6. **Accounting and log.** G evaluations split into base, Krylov and line
   search; per Newton step ‖R‖_W, η, GMRES iterations and the residual
   reached, λ (or "plain step"), guard halvings, and the evaluations so
   far; at the start every constant above (h and δ, γ, α, η_max, k, the
   line-search constants) and the memory. The per-atom `delta` lines of
   the other solvers are printed at every base point.

**Keywords:** `SOLVER = NEWTON_KRYLOV`, `NEWTON_NOISE` (10⁻¹⁰),
`NEWTON_KRYLOV_MAX` (20), `NEWTON_ETA_MAX` (0.9), `NEWTON_VARIABLES`
(`LOG` | `RELATIVE`). readInput stops on the `FIXED_POINT` conflicts
(map keywords, `NG_ORDER > 0`), on `NEWTON_*` with another solver and on
`ANDERSON_*` with `NEWTON_KRYLOV`.

**Code:** `newton.c` / `newton.h`, generic without RH globals
(matrix-free weighted GMRES with a product callback, Eisenstat–Walker
forcing), unit-tested alone; `iterateNewton()` in `solver.c`, sharing the
variable transform, guard and block setup with `iterateAnderson()`
(`G_eval` unchanged); keywords in `readinput.c`, dispatch in `iterate.c`,
counters read by `tests/analyze_run.py`.

**Tests:** unit (`test_newton`): GMRES on a nonsymmetric system matches a
dense solve with k = n, its residual history is non-increasing, a
tolerance stops it early, two runs are bit-identical; the forcing
sequence on a given residual history. e2e: a short `NEWTON_KRYLOV` run is
bit-reproducible and takes Newton steps; the keyword conflicts stop the
run.

**Acceptance:** default runs unchanged bit for bit; on the S-5 case set
the solver runs to its stopping criterion or reports failure clearly
(exit status bit), with all S-1 metrics and the counts of item 6.

**Status (2026-11-09, branch `s4/newton`):** implemented as above, with
these details fixed in the code:

- The finite-difference probes keep the base point's values outside the
  mask; a probe component whose image has no value in the variables is
  set to zero and counted. A trial point takes G(x_k)'s values outside
  the mask. GMRES is skipped (plain step) when the evaluation budget
  leaves no product and one trial.
- **Guard reference (changes item 4):** a trial x_k + λd lies near x_k
  for small λ, not near G(x_k). Judged against G(x_k) alone, as
  Anderson's extrapolations are, an inversion that x_k has and G(x_k)
  lacks blocked every λ: 155 of 243 plain steps on `1d_shock_h6`, all 16
  on `1d_falc_helarge`, and divergence by plain steps on
  `1d_shock_h6crd_ne`. A negative value or inversion now counts as new
  only where neither G(x_k) nor x_k has it.
- A trial with an entry above log(DBL_MAX) in `LOG` variables is a guard
  failure ("not finite") before exp() is called: Newton steps reach
  |d|_max ~10⁴ in J (`1d_falc_mgi66`), which trapped as SIGFPE.
- Every Newton step logs |d|_max and its block, and the guard failures by
  reason.
- The Armijo test compares norms over the same entries: a trial whose
  image loses an entry the base point has (G(trial) ≤ 0 in `LOG`) is not
  accepted, since dropping that entry from the norm could pass the test
  without a decrease ("lost an entry" in the end line).
- J does not model how the entries outside the mask move in a trial
  (to G(x_k)): a probe that moved them would carry an O(1)/ε jump. The
  line search sees it; the direction does not.
- If `N_MAX_ITER` ends the run after a failed line search, the last
  evaluation was a probe or a rejected trial. One more at the base point
  (G is deterministic: G(x_k) does not change) makes what G sets besides
  the state agree with the output ("restoring" in the end line), one
  evaluation past the cap.

First run (the e2e set-up: Ca II with two PRD lines on FAL-C, H_6
passive, `ITER_LIMIT` 10⁻⁶): converged in 124 G evaluations (18 Newton
steps; 100 Krylov products, 23 line-search trials), against 89 for
`ANDERSON`. The linear model often promised more than the step delivered
(step 7: GMRES predicted ‖R‖ 3.8×10⁻⁵, the step gave 3.2×10⁻³). It is
not the finite-difference step: h from 10⁻⁶ to 10⁻⁴ (`NEWTON_NOISE`
10⁻¹² to 10⁻⁸) gives 119–124 evaluations, while 10⁻³ and 10⁻² give 182
and 381.

**First results on the S-5 set** (2026-11-09, S-5 settings, data in
`claude-tmp/s1/s4`): G evaluations, `ANDERSON` (depth 10) → `NEWTON_KRYLOV`:
`h6_ca_mg` 39 → 56, `h20` 30 → 49, `mgi66` 24 → 77, `solvene` 40 → 69,
`falxco_h6_ca` 52 → 82, `t5750_h6` 145 → 378, `helarge` 88 → 206,
`civi` 51 → 117, `bf_c100_prd` 31 → 58, `shock_h6crd_ne` 57 → 88
(77 before the Armijo test kept the base point's entries: one trial lost
an entry and had passed);
`shock_h6` not converged after 3000 (‖R‖_W stalls near 0.07). Newton
needs 1.4–3.2× Anderson's evaluations and more wall time on every case;
all converged runs meet the S-5 accuracy criteria. On `shock_h6` the
Newton steps are huge in log variables (|d|_max 5–45 in PRD ρ and J),
the guard and line search cut λ to 1/16 or reject the step, and a cap on
‖d‖_∞ (1 or 3; scratch experiment) makes every case worse (`h6_ca_mg`
69–92, `helarge` 290–489, `t5750_h6` not converged in 600): far from the
solution the Newton direction itself is poor, not its length. As
expected before the runs, Form A does not beat Anderson. Side check
(depth 20 and 30): Anderson's slow cases are limited by history
restarts after rejections (mean 4–6 columns used at depth 30), not by
depth.

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
| Failing | `1d_t5750_h6` | yes | — | not converged (4.5×10⁻³ after 200; plain G converges in 2527) |
| | `1d_shock_h6` | yes | — | SIGFPE at iteration 11 (141 without Ng) |
| | `1d_falc_helarge` | — | — | SIGFPE at iteration 13 (note 1) |
| | `1d_falc_civi` (note 6) | — | — | SIGFPE at iteration 12 |
| | `1d_bf_c100_prd` (note 2) | yes | — | limit cycle (plain G converges in 645) |
| ne coupling | `1d_falc_solvene` | yes | yes | 27 |
| | `1d_shock_h6crd_ne` (note 3) | — | yes | period-3 cycle |
| Multi-D (only if 1-D is promising) | `2d_pert_h6_ca_mg` | yes | — | 82 at 8 threads before B-40; SIGFPE at iteration 6 since (note 4) |
| RMHD | note 5 | | | |

1. He_large: SIGFPE since #230 (81 iterations before, 93 with the
   background file in double; N-28 comment of 2026-11-07); stalls at
   ~2×10⁻³ without Ng (N-28). No reference (S-1): judged on ‖F_SE‖ only.
   With GTH (N-31) it converges in 63 iterations, and in 667 with the
   reference recipe, which gives it a reference (to be made with the
   regenerated S-5 references).
2. Chosen in S-1: Bifrost column 100 (N-4), Ca II with angle-averaged
   PRD; needs `--bifrost`. Without the snapshot, `1d_falc_h6_tight`
   (H_6 to 10⁻⁸, N-4) is the in-repository stand-in.
3. `FALC_shock`, H_6_CRD, `SOLVE_NE = ITERATION`; added to
   `run_stress.sh` in S-1. Without Ng it alternates with period 2 (N-6).
   The only ne case without PRD.
4. At a fixed 8 threads G is bit-repeatable (angle-averaged PRD); run
   every solver at 8 threads. Since B-40 (#257) the case stops with SIGFPE
   at iteration 6 (N-33); unless N-33 is fixed first it is out of the
   comparison.
5. Development cases: the Bifrost cases in `run_stress.sh` (column 260,
   cuts `s4` and `s21`). Held out (§7 trap 8): columns and cuts not used
   for B-29, N-4, T-4 or P-8, e.g. another row of the snapshot or a MURaM
   cut; named before the runs.
6. C I–VI (Judge model), added 2026-11-09: Ng stops at iteration 12;
   with N-32 (and GTH) it converges without Ng in 71 iterations of the
   default MALI loop.

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

**Success criteria** (confirmed 2026-11-09, before any run):

- convergence level: every solver runs every case to one common
  `ITER_LIMIT`, calibrated before the runs with Ng on G as the largest
  power of ten at which the easy controls meet the correctness bar below,
  then frozen. All solvers stop on the same measure (`MaxChange` on
  (x, G(x))), and how far that leaves a run from the solution depends on
  how slowly G converges, not on the accelerator; at the shipped 10⁻³ Ng
  leaves spectrum errors of 5×10⁻⁴ to 3×10⁻² (S-1, corrected for Q-22),
  so the bar cannot be met there. The outcome at 10⁻³ is reported too,
  for users;
- correctness: converged populations within 10⁻⁴ relative of the
  reference (levels above the S-1 threshold), emergent spectra within
  10⁻⁴ (normalization from S-1); for cases without a reference, ‖F_SE‖
  down to the level the cases with references reach;
- no regression: every case Ng on G converges, the candidate converges
  too;
- easy controls (including `1d_falc_solvene`): no more than 20 % more
  work than Ng on G;
- slow case: at least 1.5× less work than Ng on G;
- failing cases (`1d_t5750_h6`, `1d_shock_h6`, `1d_falc_helarge`,
  `1d_falc_civi`, `1d_bf_c100_prd` or its stand-in): converge
  (correctness as above);
- `1d_shock_h6crd_ne` (plain G diverges): reported, not a criterion for
  Anderson; it is the case most likely to need Newton (S-4), and is a
  failing case for S-4;
- work includes rejected steps (S-3) and `Background()` calls (§4.1).

**Calibration (2026-11-09, frozen):** `ITER_LIMIT` = 10⁻⁶ and
`N_MAX_SCATTER` = 50 for every solver and case, against the references at
f5e9847 (after N-31, N-32, B-42, B-43); for the ne cases the references
are made with ne in the state (`FIXED_POINT`, plain G), see below. Ng on G
on the easy controls, error against the references, populations (levels
above the S-1 threshold) / spectrum:

| `ITER_LIMIT` | `1d_falc_h20` | `1d_falc_h6_ca_mg` | `1d_falc_mgi66` | `1d_falc_solvene` |
|---|---|---|---|---|
| 10⁻⁴ | 7.0e-4 / 8.5e-4 | 4.4e-4 / 8.8e-4 | 1.2e-4 / (1.9e-2) | 1.5e-3 / 3.3e-3 |
| 10⁻⁵ | 7.9e-5 / 9.5e-5 | 4.5e-5 / 8.9e-5 | 4.8e-6 / 3.6e-5 | 1.9e-4 / 4.2e-4 |
| 10⁻⁶ | 5.6e-6 / 6.8e-6 | 3.6e-6 / 7.1e-6 | 6.6e-7 / (2.0e-3) | 5.1e-6 / 1.1e-5 |
| 10⁻⁷ | 2.1e-7 / 2.5e-7 | 4.4e-7 / 8.8e-7 | 8.3e-8 / (8.7e-4) | 5.6e-7 / 1.3e-6 |

10⁻⁶ is the largest power of ten at which all four meet the 10⁻⁴ bar
(at 10⁻⁵ solvene's populations do not, and H_20 and H_6 + Ca II + Mg II
pass with no margin). Two settings came out of it:

- `N_MAX_SCATTER` = 50 (the references' value; the stress cases use 2).
  The convergence test watches the populations only, and Ng accelerates
  only them; the background J at strongly scattering wavelengths
  converges by plain Λ iteration and can lag far behind. On
  `1d_falc_mgi66` the spectrum error at 122.6 nm stays at 10⁻³ (in
  parentheses above, with 2 final sweeps) however tight `ITER_LIMIT`;
  with 50 the final phase converges J in 24 sweeps and the error at
  10⁻⁵ is 3.6×10⁻⁵. Anderson, which extrapolates J too, reached 1.6×10⁻⁶
  with 2 final sweeps. The final sweeps count in the work totals, so a
  solver that leaves J unconverged pays for it there.
- References for the ne cases with ne in the state. The `1d_falc_solvene`
  references made with the default solver (ne outside G) and with
  `FIXED_POINT` (ne in the state, S-2b), both converged to 10⁻¹⁰, differ by
  1.7×10⁻⁴ in the spectrum at 163.2 nm and 4.8×10⁻⁶ in the populations
  (N-36). Every solver compared iterates the state formulation, so they
  are judged against its reference; against the other, Ng on G has an
  error floor of 1.7×10⁻⁴ at any `ITER_LIMIT`. The default Ng control is
  judged against the default-solver reference.

**Report:** one table per case (work counts, wall time, ‖F_SE‖, population
error, spectrum error, outcome) in `MODERNIZATION.md` §2, plus the log
excerpts showing residual histories.

**Decision rule** (proposal): if Newton meets the criteria where Anderson
does not, proceed to S-6/S-7. If Anderson meets all of them, make it the
default and keep Newton as an option, revisiting it for PRD and ne. If
neither does, analyse why (§2.3: which coupling dominates) before building
further.

**Results (2026-11-09):** the tables (work, errors, ablations, counters,
residual histories) are in `MODERNIZATION.md` §2.5. Anderson at the
default depth of the time (m = 5) meets every criterion, and so does
m = 10:

- correctness: every converged run within 10⁻⁴ of its reference (largest
  `t5750_h6`, 4.9×10⁻⁵ / 3.0×10⁻⁵);
- easy controls, work against Ng on G: `falc_h6_ca_mg` 144 / 130 against
  213, `falc_h20` 49 / 40 against 54, `falc_mgi66` 83 / 82 against 112,
  `falc_solvene` 147 / 121 against 519 (m = 5 / m = 10);
- slow case `falxco_h6_ca`: 140 / 131 against 405, 2.9× / 3.1× less;
- failing cases, all converged: `t5750_h6` 432 / 300 (Ng on G 2330, the
  default Ng not after 3000), `shock_h6` 580 / 320, `falc_helarge`
  143 / 139, `falc_civi` 133 / 151 (these three crash with Ng; plain G
  517, 416, 256), `bf_c100_prd` 89 / 77 (Ng on G 131);
- no regression: every case Ng on G converges, Anderson converges.

m = 3 fails two of them (`falc_h20` 24 % more work than Ng on G;
`t5750_h6` spectrum 1.1×10⁻⁴).

Findings:

- Depth: m = 10 needs less work than m = 5 in 9 of the 11 cases, most on
  the failing ones (`t5750_h6` 300 against 432, `shock_h6` 320 against
  580); more on `falc_civi` (151 against 133). History memory 40·N against
  25·N doubles.
- LOG beats RELATIVE except on `shock_h6`; GLOBAL beats PER_BLOCK
  everywhere (`falxco_h6_ca` 140 against 617 and 388). One coefficient
  set for the whole state is part of why Anderson beats Ng, which has
  one per atom.
- Safeguard: no effect on most cases, decisive on `shock_h6` (1721
  without against 580 with), slightly costly on `falc_helarge` (134
  against 143). On `shock_h6` c = 1 rejects 78 of 282 points at m = 5
  (31 of 146 at m = 10), and m = 5 does no better than plain G (580
  against 517): the case for the decreasing envelope of Zhang,
  O'Donoghue & Boyd (2020).
- Guard: the inversion check acts often on `falc_helarge` (105 damped, 44
  plain steps), `falc_civi`, `t5750_h6` and `falc_mgi66`, from the first
  extrapolations on, while Anderson still wins clearly there. A small
  margin in the test (levels whose n/g are nearly equal) is to be looked
  at. Done in S-13 (below): not a margin; the continuum test went.
- ne: on `falc_solvene` Anderson at m = 10 needs 4.3× less work than Ng
  on G, and 6× fewer G evaluations and `Background()` calls (40 against
  244). On `shock_h6crd_ne`, where
  plain G and Ng on G produce NaN and the reference never converged,
  Anderson converges to a solution with ‖F_SE‖ 5–7×10⁻⁷ and positive
  populations, the same for m = 3, 5, 10 to 10⁻⁵; statistical equilibrium
  falls back to LU at 49 locations there (negative preconditioned
  radiative rates).
- The default Ng control with the PRD early exit needs less work than Ng
  on G on `falc_h6_ca_mg` (163 against 213) and more on `falxco_h6_ca`
  (619 against 405): the early exit saves PRD sub-iterations where they
  are not needed, at the price of a discontinuous map (N-4).

**Decisions (2026-11-09):**

- `ANDERSON_DEPTH` defaults to 10 (code PR).
- Anderson becomes the default solver. This needs a design step before
  code: the shipped keyword files and keyword defaults that `ANDERSON`
  refuses (PRD early exit and PRD Ng, `NG_ORDER`, collisional relaxation,
  hydrostatic equilibrium, `LIMIT_MEMORY`), what replaces them, how runs
  that set them behave, the multi-D cases (not yet measured; the 2-D
  stress case is blocked by N-33), and the change of results for every
  run. The design is proposed below, under "Anderson as the default
  solver".
- S-4 (Newton–GMRES) is built anyway (decided 2026-11-09): Newton has no
  failing 1-D case left to win, but it measures Newton's outer
  convergence on the same G and is the scaffold for Form B (S-6). It
  replaces Anderson as the recommendation only if it wins in wall time;
  the expectation is in S-4.

**Open:** the envelope safeguard on `shock_h6` (S-12); N-36 (the two ne
formulations); multi-D (after N-33). The inversion guard is S-13 (below).

### Anderson as the default solver (design for the S-5 decision)

**Status:** agreed 2026-11-10, with items 2 and 3 decided stricter than
proposed (below). Work items: S-10 (the switch, items 1–6) and S-11 (the
history in multi-D, item 7). S-10 is done: the keyword checks (#296),
the Ng pins in the test tools (#297), and the switch (step 3, this
section's items 1, 3 and 4). S-11 is open. It answers the design step
that the S-5 decision asked for: the keywords and defaults that `ANDERSON`
refuses, the shipped keyword files, runs that set Ng-only keywords,
multi-D, and the change of results. Measurements are in
`claude-tmp/s1/default` (README). The build was `s4/newton` 73e1deb,
whose Anderson is bit-identical to `modernize` 0569220.

**Measured on the shipped examples.** Each program's own `run/`
directory, run to convergence: once with the shipped Ng settings
(`NG_ORDER` 2, `PRD_ITER_LIMIT` 10⁻²) and once with `SOLVER = ANDERSON`
(depth 10, Ng and the PRD early exit off). Iterations equal G
evaluations for Anderson. "History" is Anderson's own estimate of the
memory its history takes; RSS is the peak of the Ng run.

| program | unknowns | iterations at 10⁻² (Ng / AA) | at 10⁻⁵ | Ng RSS | AA history |
|---|---|---|---|---|---|
| `rhf1d` (FAL-C, Ca II) | 39 934 | 11 / 10 | 24 / 21 | 16 MB | 13 MB |
| `rhsphere` | 131 856 | 15 / 11 | 49 / 23 | 22 MB | 44 MB |
| `rhsc2d` | 118 572 | 12 / 10 | 27 / 22 | 16 MB | 40 MB |
| `rhsc3d` (4×4×82) | 629 760 | 13 / 10 | 27 / 21 | 38 MB | 210 MB |

- Anderson converges on every shipped example, including 2-D, 3-D and
  spherical geometry, which the S-5 set did not cover. It never needs more
  iterations than Ng, and needs half as many on `rhsphere` at 10⁻⁵.
- At 10⁻⁵ the two solvers agree to 2–8×10⁻⁵ in the populations. At the
  shipped `ITER_LIMIT` 10⁻² they differ by 1.1–3.0 %. That is the same
  size as each one's distance from the converged solution, measured
  against Ng at 10⁻⁵: Ng 0.9–4.5 %, Anderson 1.0–7.3 %.
- Wall time follows the iterations: `rhsphere` at 10⁻⁵ takes 16 s with
  Ng and 10 s with Anderson. The other runs take 1–37 s, too short to
  measure the overhead per iteration. So the change of default moves results by about
  the convergence error the run already accepts.
- **Memory is the obstacle in multi-D, and the cause is J.**
  - The state holds J (Nspect·Nspace) and PRD ρ beside the populations.
    In the `rhf1d` example the split is 1.2 % populations, 39 % ρ (Ca II
    H and K) and 60 % J (291 wavelengths). The 2-D (3×82) and 3-D
    (4×4×82) examples (Ca II, 286 wavelengths) have about the same
    split per point.
  - Anderson keeps about 3m + 10 vectors (40 at depth 10). RH's Ng keeps
    `NG_ORDER` + 2 copies of the populations, and `PRD_NG_ORDER` + 2
    copies of ρ (§4.3).
  - Shipped 3-D example: the history (210 MB) is 5.5× the whole Ng run.
  - 2-D perf case (128×82 points, Nspect 7175): J is 0.6 GB per vector,
    so about 24 GB of history against a 10.7 GB run.
  - 3-D estimate of `MODERNIZATION.md` §3.3: J is 420 GB per vector,
    about 17 TB of history.
  - In 1-D the history is at most 2× a small run's whole memory.

**Decisions** (proposed, then agreed 2026-11-10 as marked):

1. **Scope (agreed).** `ANDERSON` becomes the default of the 1-D programs (`rhf1d`,
   `rhsphere`). `rhsc2d` and `rhsc3d` keep `NG` as their default until
   item 7 is done. The default is chosen per program and printed at
   the start ("SOLVER = ANDERSON (default of rhf1d)"). STIC runs 1-D
   columns, so it gets Anderson either way. *Alternative:* one default for
   all programs, with multi-D users setting `SOLVER = NG`. That fails
   badly at realistic multi-D sizes (memory above), so I don't recommend
   it.
2. **Every keyword that is set is checked (decided).** The user's rule:
   always check the input and reject any input that is invalid, so that a
   run never surprises its user. As proposed, the checks look at whether a
   keyword is set in keyword.input, not at its default; but no check
   depends on the atoms.
   - A default is never a conflict. `PRD_ITER_LIMIT` defaults to 10⁻²
     for `NG`; under a G solver an unset `PRD_ITER_LIMIT` is 0 (a fixed
     `PRD_N_MAX_ITER` sub-iterations).
   - A setting that the solver would not honour stops the run, whether
     or not the atoms have PRD lines. Keywords of another solver:
     `ANDERSON_*` and `NEWTON_*` (as before); the Ng keywords under
     `ANDERSON` and `NEWTON_KRYLOV` (`FIXED_POINT` keeps them for the
     "Ng on G" control); `PRD_NG_*`, `PRD_ITER_LIMIT` and `COLL_RELAX_*`
     under every G solver.
   - A value that means "off" states what the solver does, and is
     accepted: `NG_ORDER`, `PRD_NG_ORDER`, `PRD_ITER_LIMIT` and
     `COLL_RELAX_NSTEP` = 0, `COLL_RELAX_FACTOR` = 1. `NG_DELAY`,
     `NG_PERIOD`, `PRD_NG_DELAY` and `PRD_NG_PERIOD` have no such value
     and stop the run when they are set.
   - `PRD_N_MAX_ITER` < 1 is valid under every solver: it is RH's
     switch for "every line in CRD" (`readAtom`), and G then has no PRD
     sub-iterations. S-2's check that refused it is removed (B-46).
     `PRD_ITER_LIMIT` and `PRD_NG_*` are still checked in such a run.
     `run_tests.sh`'s `1d/crd` needs no pin once item 4 has removed
     `PRD_ITER_LIMIT` from the shipped file.
   - The test tools gain a way to remove a keyword from a shipped file
     (`run_stress.sh --kw '-NG_DELAY'`), since `rhsc2d` and `rhsc3d`
     keep their Ng lines.
3. **Keyword files that set Ng-only keywords and no `SOLVER` (decided:
   no fallback to `NG`).**
   - readInput stops, by the checks of item 2. The messages are generic:
     per keyword, that it is set but `SOLVER = ANDERSON` (the default)
     does not use it, and to set `SOLVER` or remove the keyword; then
     one error that refers to the manual.
   - The details (which keywords, why the default changed, how to keep
     Ng) go in the manual (Q-19) and the release notes, not in the
     messages.
   - The default is printed at the start of every run (item 1).
4. **Shipped keyword files (agreed).**
   - `rhf1d/run` and `rhsphere/run` drop `NG_*` and `PRD_ITER_LIMIT`. They
     gain a commented `SOLVER` line that names the choices and points to
     SOLVER.md.
   - `rhsc2d/run` and `rhsc3d/run` are unchanged under item 1.
   - CI's drift report against upstream then runs upstream without Ng on
     those two files, so its numbers change once. That gets explained in
     the PR.
5. **`N_MAX_ITER` and `ITER_LIMIT` keep their meaning (agreed).** The stop is
   MaxChange against `ITER_LIMIT`, and one Anderson iteration is one G
   evaluation, as one Ng iteration is one MALI step. No keyword changes
   meaning. The exit status bits (N-1) are unchanged.
6. **Results change: the switch PR carries `changes-results` (agreed).**
   - CI's base comparison reports the change for the 1-D cases. After the
     merge the base runs Anderson as well.
   - The test tools that mean Ng pin it with `SOLVER = NG`:
     - `make_references.sh` (its plain-MALI recipe);
     - `run_stress.sh`'s Ng controls and the cases that test Ng
       (`1d_falc_prdng`, `1d_falc_nong`);
     - the e2e tests `prd_hybrid` and `stokes_1d`, which test properties
       after three iterations of the shipped rhf1d file. Each of these
       either still holds under Anderson or pins `SOLVER = NG`; this gets
       checked case by case.
   - Done (S-10 step 2): `make_references.sh` pins `SOLVER = NG`; the
     stress cases `1d_falc_nong`, `1d_falc_h6_tight`, `1d_falc_prdng`
     and `1d_falxco_collrel` pin `SOLVER = NG` with the shipped Ng
     settings written out, since item 4 removes them from the file.
     `prd_hybrid` holds under Anderson once item 4 has removed the Ng
     lines (all four checks pass with differences of 0, run with
     `SOLVER = ANDERSON` and those lines removed), and `stokes_1d` does
     no iterations (`N_MAX_ITER = 0`); neither needs a pin, but both
     need item 4 in the same PR as item 1. `run_tests.sh`'s `1d/crd`
     needs B-46 (#301).
   - The other stress cases run the default solver: they test what
     users get. Their group comments describe Ng's behaviour and are
     updated in step 3.
   - The manual's solver section waits for Q-19. The change fragment and
     the release notes say what changed and how to keep Ng.
7. **The way to a multi-D default (agreed; S-11).**
   - N-33 must be fixed first. Then `2d_pert_h6_ca_mg` gets measured.
   - The history must shrink. Candidates:
     - (a) Extrapolate only the populations (perhaps ρ too), and carry J
       (and ρ) from the last G, as Ng does. The history then shrinks to
       a few percent of today's. This needs measuring: it was not tested
       in S-3, and the carried part lags by one step. 1-D first, cheaply
       on the S-5 set.
     - (b) Store the J history in float: half the size, and not enough
       on its own.
     - (c) A smaller depth in multi-D. Depth 5 met every S-5 criterion,
       but it saves only 40 %: 25 vectors instead of 40.
   - (a) is §9 Q1 asked for Anderson; it is S-11 (design in the next
     section).

**Order of work.**
- One PR per step:
  1. Item 2: the checks on what is set, and the tools' way to remove a
     keyword. No results change under the current default.
  2. The `SOLVER = NG` pins in the tools of item 6. No results change.
  3. The switch: items 1, 3 and 4 and the change fragment, with the
     `changes-results` label.
- Log items S-10 (steps 1–3) and S-11 (item 7(a)).

### S-11 — Anderson's history in multi-D (design)

**Status:** agreed 2026-11-10, with the criteria and defaults of items 5
and 6 decided by the user (revised the same day: viable against Ng,
the 1-D switch on the mean ratio). It implements item 7(a) of
the default-solver design above (§9 Q1 for Anderson). Step 1 (the code)
is #310. Step 3 (the defaults): the user decided 2026-11-13 from the 2-D
measurements that `ANDERSON` with `POPULATIONS_RHO` is the default of
`rhsc2d` and `rhsc3d`, and that 1-D keeps `ALL` (the numbers are in
`docs/changes/S-11.md`). The write-up of step 2, with the 1-D
measurements redone after S-13, is still to come.

- The plumbing of decision 2 is checked by the solver itself. Before
  every evaluation after the first, the carried blocks of x are
  compared bit for bit with those of G at the last accepted point, and
  the end of the run reports in how many evaluations they agreed. The
  check covers the packing and the copy of the carried blocks in
  `fromVars`; it cannot see an error in G. The e2e test forces
  rejections (`ANDERSON_SAFEGUARD = 0.5`) and requires every
  evaluation to agree. A mutation that skips the copy of the carried
  blocks fails it.
- First look on the shipped `rhf1d` example (Ca II, two PRD lines, 1
  thread, `ITER_LIMIT` 10⁻⁵): `ALL` takes 21 iterations, with 12.4 MB
  of history and variables and 20.9 MB peak RSS. `POPULATIONS` takes 20
  iterations, 0.15 MB and 16.2 MB. `POPULATIONS_RHO` takes 21
  iterations, 5.0 MB and 16.1 MB. The three full states (x, G(x) and G
  at the last accepted point) add 0.96 MB in every setting. The
  comparison itself is step 2.

**Problem.** Anderson keeps (3m + 10) vectors of the whole state, 40 at
depth 10. The state is mostly J. In the shipped `rhf1d` example, 60 % of
it is J and 39 % is PRD ρ. On the 2-D perf case the J part of the
history alone would take about 24 GB (about 26 GB with ρ and the masks,
item 3), against 10.7 GB for the whole run. On the 3-D estimate of
`MODERNIZATION.md` §3.3 it would be about 17 TB.

**Idea.** Extrapolate only part of the state and carry the rest from
the last G, as Ng does. RH's Ng extrapolates populations only; J is
whatever the last formal solution gave, and ρ has its own Ng in the PRD
sub-iterations. A carried entry of the next point takes the value of G
at the last accepted point. That is G(x_k) after an accepted evaluation,
the plain-iteration value. After a rejection the whole next point is
G of the last accepted point, as now (decision 2). The code already does this for masked entries (`fromVars`), so
the change is mostly which entries Anderson's vectors hold.

**What changes in the method.**

- The least squares sees only the residuals of the extrapolated blocks.
  With `GLOBAL` coefficients, J's residual no longer steers γ. S-5 found
  `GLOBAL` better than `PER_BLOCK`, and part of that may be J's
  coupling. Only the measurement can tell.
- The safeguard norm ‖f‖_W is taken over the extrapolated blocks only.
  The convergence test is unchanged: `MaxChange` already looks only at
  populations.
- The guard is unchanged. Carried blocks hold G's values, so they cannot
  create new inversions or negative values.
- J is out of step with the populations. In plain iteration J and the
  populations of x_{k+1} come from the same G. Here the populations jump
  to the extrapolated point, while J stays G's J of the unextrapolated
  x_k. RH's Ng has the same mismatch: the next formal solution pairs
  the extrapolated populations with the J of the last one. G uses the J
  of x only in the background scattering source term (`sca_c`·J in
  `formal.c`, and J20 in the background polarization); the PRD
  sub-iterations use the J of G's own formal solution. Whether the
  mismatch costs iterations or rejections in scattering-dominated
  continua (large `sca_c`) is for the measurement to show. ρ of x, by contrast, enters the line emissivity
  of the first formal solution directly. That is why `POPULATIONS_RHO`
  is a candidate.

**Decisions** (agreed 2026-11-10):

1. **Keyword.** `ANDERSON_EXTRAPOLATE = ALL | POPULATIONS |
   POPULATIONS_RHO`.
   - `POPULATIONS` covers atoms, molecules and ne (when ne is in the
     state). ne is small and is the slow mode of the ne cases.
   - `POPULATIONS_RHO` adds the ρ of PRD lines.
   - J and J20 are always carried under the two reduced settings.
   - The default is `ALL` until the measurements, so runs are unchanged
     bit for bit. Item 6 gives the defaults after them.
   - Like every `ANDERSON_*` keyword, it stops the run under another
     solver.
   - `POPULATIONS_RHO` on atoms without PRD lines runs as
     `POPULATIONS`. It does not stop the run, because no keyword check
     depends on the atoms (S-10, item 2), but the solver warns that it is
     the same as `POPULATIONS` (decided by the user 2026-11-10). The start
     banner names the extrapolated blocks.
2. **Carried blocks on rejection.** A rejected point is followed by the
   plain step from the last accepted pair. That step is G of that pair,
   carried blocks included (`gacc`, as now). The rule is the same as
   for every other entry.
3. **Memory, first stage.**
   - The Anderson object and u, g(u), u_next, w and the masks cover the
     extrapolated blocks only, packed by a block list.
   - x, G(x) and `gacc` stay full: 3 full vectors. These copies are the
     floor of this stage. On the 2-D perf case they take about 1.9 GB;
     on the 3-D estimate they are 3 × 420 GB of J.
   - The printed memory becomes 3·N + (3m + 7)·N_x doubles plus
     (m + 4)·N_x bytes, where N_x is the length of the extrapolated part.
     With `RELATIVE` the scales add N_x doubles.
   - Estimate for the 2-D perf case at m = 10: `ALL` about 26 GB,
     `POPULATIONS` about 2.1 GB, `POPULATIONS_RHO` about 3.5 GB. The
     estimate assumes ρ is about 6 % of J per point (4 angle-averaged
     PRD lines); this gets measured.
4. **Memory, second stage (open; only if the 3 copies matter).** The
   carried blocks could stay in RH's own arrays: G_eval would scatter
   and gather only the extrapolated part.
   - That removes x and G(x) for J.
   - A rejection then needs J of `gacc`, which is one copy. Or it uses
     J of the rejected evaluation, which changes the safeguard's
     meaning.
   - Neither variant is proposed now.
5. **Measurement, with criteria frozen before the runs (§7 trap 8).**
   1. The S-5 set, at m = 10 with the frozen S-5 calibration (1 thread,
      `ITER_LIMIT` 10⁻⁶), with all three settings.
      - Reported: G evaluations, wall time, rejections, guard actions,
        and the S-1 metrics against `refs2`.
      - A reduced setting is **viable** (decided by the user) if it
        converges on every case where Ng (`SOLVER = NG`, the old
        default) converges, meets the S-5 correctness tolerance there,
        and beats Ng as in the S-5 criteria. It need not match `ALL`:
        beating Ng is enough for it to be an option.
   2. The shipped examples (`rhf1d`, `rhsphere`, `rhsc2d`, `rhsc3d`),
      as in the measurements above, plus peak RSS.
   3. The 2-D perf case: peak RSS over a few iterations under `NG` and
      the three settings, to check the estimate of item 3.
   4. `2d_pert_h6_ca_mg` once N-33 is fixed.
   5. Held-out inputs (§7 trap 8, S-5 note 5), named now: Bifrost row
      0, column 400 (`bifrost:1d:400`, Ca II with PRD), and the 2-D cut
      of every 12th column (`bifrost:2d:12`, Ca II with PRD). Neither was
      used for development. A setting that is viable on the S-5 set but
      fails on these is not viable.
6. **Defaults (decided by the user).** Sensible defaults, documented. No
   choice at run time from the available memory: users who need another
   setting set the keyword.
   - 2-D and 3-D (`rhsc2d`, `rhsc3d`): `POPULATIONS_RHO`, if it is
     viable there.
   - 1-D (`rhf1d`, `rhsphere`): `ALL`. Memory is almost never a problem
     in 1-D, so a reduced setting must not cost much work there.
     `POPULATIONS_RHO` becomes the 1-D default as well, so that one
     setting holds everywhere, if both of these hold on the 1-D cases
     (the S-5 set, the shipped 1-D examples and `bifrost:1d:400`):
     - it converges on every case where `ALL` converges;
     - on average it needs at most 1.1× the G evaluations of `ALL`.
       That is the mean over the cases where both converge of the
       per-case ratio (`POPULATIONS_RHO` / `ALL`), so that every case
       counts the same whatever its length.
   - `POPULATIONS` stays an option. If `POPULATIONS_RHO` is not viable
     in multi-D, the multi-D default goes back to the user with the
     numbers.
   - Making `ANDERSON` the default of `rhsc2d` and `rhsc3d` is a
     separate PR with `changes-results`, after N-33 and the 2-D
     measurements.

**Order of work.** One PR per step:
1. The code, with `ALL` as the default. No results change
   (`run_tests.sh --exact`).
2. The measurements of item 5, written up in this section and
   `MODERNIZATION.md` §2.
3. The defaults of item 6 (`changes-results`).

**Code.** All in `solver.c`.
- A per-block flag set from the keyword.
- Two offset tables. `blockStart` keeps indexing the full x, which
  `guardTrial`, `printChange` and the J offset for J20 in
  `toVars`/`fromVars` use. A new packed table indexes u.
- `toVars` maps full x to the packed u. `fromVars` maps the packed u and
  full G(x) to full x: carried blocks from G(x), the rest as now.
- `setVarScale` and the weights work over the packed blocks.
- `PER_BLOCK` solves one system per extrapolated block. `opt.Nblock`
  and `opt.start` come from the packed table, and `logSolve` names
  blocks through it.
- The safeguard norm and the rejection message use the packed length
  (`andersonNorm(N_x, …)`).
- The start banner names the extrapolated blocks and gives the memory
  of item 3.
- `anderson.c` stays generic and is unchanged.
- Newton-Krylov shares `setupVariables`, `toVars` and `fromVars`. It
  uses them with every block in the packed table (the `ALL` packing),
  and its results must not change.

**Tests.**
- e2e:
  - `ALL` repeats the current Anderson bit for bit.
  - With `ANDERSON_DEPTH = 0`, every setting repeats `FIXED_POINT` with
    `NG_ORDER = 0` bit for bit.
  - Under `POPULATIONS`, the J of every evaluated point equals the J of
    G at the last accepted point, bit for bit, rejections included
    (decision 2). The exact form of the check gets settled in the PR.
  - The `NEWTON_KRYLOV` e2e results are bit-identical before and after
    the change.
  - A short `POPULATIONS` run is bit-reproducible.
  - The keyword under another solver stops the run.
- Unit: none new. `anderson.c` does not change.

### S-13 — The inversion guard (measured, decided 2026-11-11)

The S-5 note above asked for a margin in the guard's inversion test. A
diagnostic build (scratch, not committed) printed every violation of each
rejected trial point (Anderson, m = 10, the S-5 calibration):

- `falc_helarge`: 149 rejected points, median 15 violations each, all
  continua: high He I levels to the ion at 9–18 kK (hν₀/kT 0.3–1) and
  He II levels to He III at 7–8 kK. G's own point is 3–12 % (median)
  from the threshold, the trial 3–9 % past it: not a borderline flip. At
  236 of the 246 (continuum, depth) pairs that fired, the converged
  solution is itself inverted (median ratio 0.69, the same in the
  plain-G reference converged to 10⁻¹⁰). At one depth G's iterates go
  1.135, 1.096, 1.067, 1.044, 1.027, 1.012, 1.0005 and then cross: the
  test rejected steps towards the solution until G got there by itself.
- `falc_mgi66`: 23 points, lines with small hν₀/kT (10⁻²–10⁻¹, nearly
  degenerate levels) near the threshold; 24 of the 47 (line, depth) pairs
  are inverted at convergence, the rest only crossed on the way.
- `t5750_h6`: all 11 guard actions are the H n = 4 continuum (hν₀/kT
  2.7) in the cool upper layers (3.7–3.9 kK); 26 of the 28 (continuum,
  depth) pairs are inverted at convergence (ratio 0.946, the same in the
  reference). `shock_h6`: no guard actions; its rejections are the
  residual safeguard (S-12).

Inversions of individual transitions are therefore physical, and RH's
MALI iteration never guarded them; what breaks the formal solver is a
negative total opacity (trap 5). A strong inverted line can make it
negative at line centre (H_20's Rydberg lines under 10 % FD perturbations,
S-2). No run with the continuum test off failed (the S-5 set;
`2d_pert_h6_ca_mg` for 904 iterations). That is evidence, not a bound:
an active atom's continua are not in the background, so an inverted
continuum of the dominant absorber could in principle make the total
opacity negative. Work (formal solutions) with
each part of the test, S-5 calibration:

| case | both (S-3) | lines only | none |
|---|---|---|---|
| `falc_helarge` | 139 | 109 | 105 |
| `falc_mgi66` | 82 | 82 | 76 |
| `t5750_h6` | 300 | 310 | 310 |
| `falc_civi` | 151 | 151 | 146 |
| `shock_h6crd_ne` | 109 | 95 | 95 |

The other six cases are the same in all three columns. With no test the
whole S-5 set converged, with errors within 10⁻⁴. Only `t5750_h6` costs
more without the continuum test (3 %), although its continua are
inverted at convergence too. **Decided (user):**
keep the line test, drop the continuum test. If a case needs the
continuum test after all, it comes back. Results of the change on the
S-5 set and `2d_pert_h6_ca_mg` are in the PR (`docs/changes/S-13.md`).
In 2-D it is neutral on `2d_pert_h6_ca_mg` (the first iteration below
10⁻³ and 10⁻⁴: 165 → 171 and 364 → 349; that case has a floor, N-37) and
lets Anderson extrapolate on `2d_perf_many`, where the old test rejected
every trial point.

**Follow-up (S-14, N-40, decided by the user 2026-11-13).** Would a check
on the total opacity itself do better than the inversion proxy? Counted
on the S-5 cases where the guard acts (scratch build, `claude-tmp/s1/s13/chi`):
the total opacity χ is ≤ 0 at G's own points (the ZERO_RADIATION starts
of `falc_helarge`, 118 000 (wavelength, depth) points down to −44 χ_c;
`falc_civi`; `falc_mgi66`) and at `falc_helarge`'s converged solution (504
points, down to −6.6 χ_c, IR He I 4296 nm), so rejecting χ ≤ 0 would stop
every solver there. What breaks the formal solver is `exp` overflowing,
at a negative optical depth of about −709; the deepest run of negative
slant optical depth was −0.094 (`falc_mgi66`, with the line test on;
−0.026 without), −4.9×10⁻⁶ on `falc_helarge`. The line test does not
prevent negative opacity, and no case needs it: it went too (S-14). In
its place a detector in the formal solution (N-40) stops the run with a
message, or rejects the trial point, where a run of negative optical
depth goes below −50.

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

- With each program's default solver (`ANDERSON` in 1-D and spherical
  geometry since S-10, `NG` in 2-D and 3-D), the existing `run_tests.sh`
  (4 iterations; the CI gate compares with the base commit, `--exact` for
  refactors, and reports drift from upstream) is unchanged. Its drift
  from upstream in 1-D and spherical geometry includes the change of
  solver.
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
   coefficient jumps ≤ 4×10⁻¹²; `S_LINEAR` is a smooth control). Measured
   in S-2: Bezier3 limits FD products of ρ and J to ~10⁻³, and with ne in
   the state the float background limits them to ~3×10⁻³, much worse than
   its ECnoise δ suggests (correlated steps, not noise).
4. **Unscaled variables.** Populations over ten decades: inner products,
   least squares and finite differences in raw units are dominated by a
   few entries. Scale first (§3.1).
5. **Positivity.** Newton and Anderson steps can produce negative
   populations, which then give negative opacities and SIGFPE (He_large
   dies this way under Ng: `exp(−dt)` overflows in `Bezier3_coeffs`). Log
   variables or step truncation, designed in from the start. Positive
   populations are not enough: independent relative perturbations of
   closely spaced levels invert them (H_20 Rydberg levels at 10⁻¹, S-2),
   which gives negative line opacity and the same overflow. Trial points
   of Newton and FD products must respect the ordering n_u g_l < n_l g_u
   where the levels are close, or the formal solver must accept masers.
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
interpolation), N-16, and B-7; B-38 and B-39 are fixed in their own pull requests (#251,
#250), ahead of S-2. N-10 becomes a
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
11. ~~Background precision with ne in the state of a Newton solver~~
    Decided 2026-11-08: double whenever ne is in the state (§3.7).
12. ~~He_large: is the δ ≈ 2×10⁻³ of its G an ill-conditioned rate matrix
    (N-7), and is that what stalls it (N-28)?~~ Answered 2026-11-08: yes,
    and yes. LU loses the populations of tightly coupled levels; GTH
    elimination does not (N-31, section after S-2).

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
- Grassmann, W. K., Taksar, M. I., Heyman, D. P. 1985, Oper. Res. 33,
  1107: regenerative analysis and steady state distributions for Markov
  chains (GTH elimination) ✓
- Hubeny, I., Mihalas, D. 2014, Theory of Stellar Atmospheres, Princeton,
  chapter 13
- Janett, G. et al. 2021, A&A 655, A13: stationary iterative methods for
  polarized transfer (paper I)
- Kelley, C. T. 1995, Iterative Methods for Linear and Nonlinear
  Equations, SIAM (Frontiers in Applied Mathematics 16), chapter 6:
  inexact Newton, the safeguard against oversolving
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
- O'Cinneide, C. A. 1993, Numer. Math. 65, 109: entrywise perturbation
  theory and error analysis for Markov chains ✓
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
