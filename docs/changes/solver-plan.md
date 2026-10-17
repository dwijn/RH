ids:
date: 2026-10-16
title: Solver redesign planned (S track); plan steps 3 to 6 revised

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
