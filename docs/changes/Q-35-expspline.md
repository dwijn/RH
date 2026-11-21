ids: Q-35
date: 2026-11-20
title: expspline.c and tools/hcollisions.c removed (follow-up to Q-35)

`expspline.c` (`exp_splineCoef`, `exp_splineEval`) was left in place by Q-35 because `tools/hcollisions.c` still called it. That tool is its only production user, and it cannot run as shipped: it reads two data files (`Giovanardi_etal.dat`, `Whelan.dat`) from hard-coded paths in the upstream author's home directory, and neither is in the repository. Both files are therefore removed.

`tools/hcollisions.c` computed hydrogen collisional excitation rates (the header cites Giovanardi et al. and Whelan, 1s-2s and 1s-2p transitions, printed as the coefficient CE of `C_ji = Ne CE gi/gj sqrt(T)`). Nothing in the repository ties it to the tables in the shipped atom files, and this change makes no claim that it produced them. Git history keeps both files: the last commit that has them is 77e50f7 (`modernize` before this change).

**Removed:** `expspline.c`, `tools/hcollisions.c`, the `expspline.c` entry in the `CMakeLists.txt` source list, the `hcollisions` target in `tools/CMakeLists.txt`, and the `exp_splineCoef`/`exp_splineEval` prototypes in `rh.h`. `tools/laguerre.c` stays (the `impact` tool uses it). A grep of C, CMake, docs, Python and IDL found no other reference apart from the logged history (`docs/MODERNIZATION.md`, changelog and older change notes, unchanged).

**Tests (`tests/unit/test_interp.c`):** `test_spline_linear_data`, `test_spline_reference` and `test_spline_limits` tested only the exponential spline and are gone. The rest of `test_interp` (Hunt/Locate, Linear, cubic convolution, pchip) is unchanged, and its header comment is updated. The spline reference data (`spline_ref_*`, including `spline_ref_natural`, and `expspline_ref_*`) had no other user, so it is removed from `tests/unit/kernels_reference.h`, together with the matching generators in `tests/unit/make_kernels_reference.py`. `tests/README.md` and `docs/unit-tests/unit-tests.tex` no longer describe the exponential-spline tests.

**Results:** no production code calls the removed functions, so results are unchanged: `run_tests.py --exact` is byte-identical.
