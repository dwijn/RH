ids:
date: 2026-10-16
title: Continuous integration on Linux and macOS; Voigt complex arithmetic

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
