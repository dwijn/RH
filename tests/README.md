# RH tests

Three tiers:

- **Unit tests** (`tests/unit/`, run by `ctest`): single routines on
  problems with known answers. Fast, no input files.
- **Regression** (`run_tests.sh`, below): whole runs of a candidate build
  against a baseline build. Must pass on every pull request.
- **Stress** (`run_stress.sh`, further below): convergence and performance
  on hard and large problems.

## Unit tests

```
cmake -B build && cmake --build build -j
ctest --test-dir build -L unit --output-on-failure
```

Each test is a small program linked against `librh` (see
`tests/unit/CMakeLists.txt`); `check.h` has the checks. References come from
analytic results or, where there is no closed form, from high-precision
computations committed as data together with the script that made them
(never from RH itself).

| Test | What it covers |
|---|---|
| `test_voigt` | `Voigt()`: H and the Faraday-Voigt function F against the Faddeeva function from mpmath (`make_voigt_reference.py`), relative error where H matters and error relative to line center; exact symmetry in v; finite values; continuity where the two approximations meet. `--table` prints the errors per damping parameter, `--bench` the time per call |
| `test_linalg` | `SolveLinearEq`, `LUdecomp`, `LUbacksubst`: random diagonally dominant systems (N up to 60), a pivoting case with exact integer answer and determinant, badly scaled rows, the Hilbert matrix scaled to integers (cond 1.5e10, with and without `improve`), the inverse |
| `test_interp` | `Hunt`/`Locate` (every starting guess, ascending and descending tables, values at, between and beyond the nodes), `Linear`, `BiLinear`, natural and exponential splines (exact for linear data; mpmath reference solutions; convergence order; tension limits), `cc_kernel`/`cubeconvol` (kernel identities, exact for biquadratics including boundary cells, order 3) |
| `test_quadrature` | `GaussLeg` (exact for x^k, k <= 2n-1, n up to 100; remainder term for x^2n; closed forms; node and weight properties), `w2`/`w3`/`U3` against incomplete-gamma references from mpmath across the series/exact/large-dtau switches |
| `test_functions` | `E1`/`E2` vs mpmath, `Planck` vs mpmath, Wien and Rayleigh-Jeans limits and the Stefan-Boltzmann integral, `w3js`/`w6js` vs closed forms, orthogonality, symmetries and sum rules. Reference tables for `test_interp`, `test_quadrature` and `test_functions` come from `make_kernels_reference.py` |

The other tests also take `--table` to print the measured errors. Error bounds in the tests are about twice what is measured, so that other
compilers and FMA contraction (Apple silicon) pass, while a real change in
accuracy fails.

# Regression tests

Compare a *candidate* build against an unmodified *baseline* build on small,
fast runs (the whole suite takes a few minutes).

```
cmake -B build && cmake --build build -j                  # the candidate
cmake -B build-asan -DRH_SANITIZE=address,undefined && cmake --build build-asan -j
tests/run_tests.sh --baseline /path/to/upstream-tree --candidate . [--asan build-asan] [--exact]
```

The usual baseline is upstream RH at c95bdf6, built without Fortran by
`tests/build_upstream.sh` (this is also what CI does). Against upstream use
tolerance mode: the Stokes cases differ by ~1e-10 because the Voigt code
uses the compiler's complex arithmetic. `--exact` is for refactors, with a
build of the branch's base (e.g. `modernize`) as the baseline:

```
git clone https://github.com/han-uitenbroek/RH upstream && git -C upstream checkout c95bdf6
tests/build_upstream.sh upstream
```

A build is given as a source tree with a CMake build in `build/`, as a CMake
build directory, or as a source tree built in place with the old makefiles
(upstream RH, the usual baseline). `rhf1d`, `rhsc2d`, `rhsc3d` and
`rhsphere` must be built. Run `run_tests.sh --help` for the options.

| Case | What it covers |
|---|---|
| `1d/nostokes`, `fieldfree`, `fullstokes` | all Stokes modes (FIELD_FREE iterates in 1 record then solves in 4) |
| `1d/crd` | no PRD lines |
| `1d/twoactive` | two active atoms |
| `2d/small`, `3d/small`, `sph/small` | the other geometries, from the shipped examples |
| `2d/mem` | peak memory on a larger tiled grid (buffers must not scale with Nspect) |
| `1d/kurucz_blank` | Kurucz list with a leading blank line (currently crashes upstream) |

Each case checks the candidate at 1 thread against the baseline, the candidate
at other thread counts against its own 1-thread result, and (with `--asan`)
that AddressSanitizer is clean. Outputs are compared per wavelength, relative
to the peak. The defaults are based on the noise of unmodified code against
itself (about 1e-11 in 1-D/3-D/spherical, 2e-9 in the ill-conditioned 2-D
example), since floating-point sums depend on thread count. Use `--exact`
for changes that must not alter a single bit (memory management).

The log, `keyword.input` and timing of every failed run are kept in
`WORK/failed/` (CI uploads them as an artifact), since run directories are
reused by later cases.

Running the suite with the same unmodified tree as baseline and candidate is
a useful self-test: everything passes except `1d/kurucz_blank`, which exposes
the upstream segfault.

## Stress, convergence and performance cases

`run_tests.sh` stops every case after 4 iterations. `run_stress.sh` runs
strongly NLTE problems **to convergence** (`--set conv`, 24 cases) or for a
fixed number of iterations on large multi-D grids (`--set perf`), on one or
more builds, and reports iterations, final population change, trend
(monotone / oscillating / stalled), wall time, peak memory, and the population
and spectrum differences between builds.

```
tests/run_stress.sh --tree upstream=/path/to/upstream --tree candidate=. \
    --set conv --jobs 8 --work /large/scratch/rh-stress
tests/run_stress.sh --list          # the cases and what each one exposes
```

The cases include known failures (non-convergence, crashes, wrong answers);
the case table in the script groups them by what they expose. 1-D cases run on
1 thread so that builds compare without summation-order noise. Multi-D
atmospheres are generated by `make_perturbed_atmos.py` (FAL-C with a
corrugated transition region and velocity field); `atmos/` holds two dynamic
1-D models. `analyze_run.py` does the log parsing and population comparison.
Use a `--work` directory on a large filesystem: PRD runs write GB-sized
scratch files.
