# RH tests

Four tiers:

- **Unit tests** (`tests/unit/`, run by `ctest`): single routines on
  problems with known answers. Fast, no input files.
- **End-to-end tests** (`tests/e2e/`, run by `ctest`): whole runs of a
  program on problems with exact properties.
- **Regression** (`run_tests.py`, below): whole runs of a candidate build
  against a baseline build. Must pass on every pull request.
- **Stress** (`run_stress.py`, further below): convergence and performance
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
| `test_coupling` | `addPreconditionedRates`: the radiative rates with per-transition preconditioning (B-6) for random active sets, including a shared lower level, a ladder and continua into an ion line. At the old populations the net rate into every level equals the plain rates with I (the fixed point is exact); at other populations the rates are the derived linear ones; no rates between levels without a transition; per-point angle weights |
| `test_negpop` | the backstop of `statEquil` for the LU fallback (N-8): negative or non-finite LU populations stop the run outside a trial evaluation (child process, ERROR_LEVEL_2) and invalidate the trial inside one (`negTauInvalid`, message with atom, level, location), non-negative rates do neither |
| `test_prdrho` | `PRDRhoCheck` (N-4): rho of a PRD line after a sub-iteration, with and without angle dependence: all positive passes; a zero, a negative value, a NaN and an infinity each stop the run outside a trial evaluation (child process, exit status 70) and invalidate the trial inside one, with a message naming atom, line, wavelength index (ray and direction), location, and the negative line emissivity; rho is not clipped; two bad values are counted |
| `test_interp` | `Hunt`/`Locate` (every starting guess, ascending and descending tables, values at, between and beyond the nodes), `Linear`, exponential splines (exact for linear data; mpmath reference solutions; tension limits), `cc_kernel`/`cubeconvol` (kernel identities, exact for biquadratics including boundary cells, order 3) |
| `test_quadrature` | `GaussLeg` (exact for x^k, k <= 2n-1, n up to 100; remainder term for x^2n; closed forms; node and weight properties), `w2`/`w3`/`U3` against incomplete-gamma references from mpmath across the series/exact/large-dtau switches |
| `test_functions` | `E1`/`E2` vs mpmath, `Planck` vs mpmath, Wien and Rayleigh-Jeans limits and the Stefan-Boltzmann integral, `w3js`/`w6js` vs closed forms, orthogonality, symmetries and sum rules. Reference tables for `test_interp`, `test_quadrature` and `test_functions` come from `make_kernels_reference.py` |
| `test_formal_1d` | 1-D formal solvers (linear and cubic Bezier short characteristics, Feautrier and Feautrier–Hermite) on slabs with a polynomial source function, where I± is known analytically: exactness for linear S, order of convergence under grid refinement with constant and exponentially varying opacity, the diagonal operator Psi against dI/dS, and the Bezier integration weights against mpmath (`make_bezier_reference.py`). `-v` prints the measured errors |
| `test_formal_2d` | 2-D linear and Bezier short characteristics on a horizontally uniform periodic slab, against the 1-D analytic solution along each ray (`slab_exact.h`): a mesh where all upwind points lie on layers (must match 1-D) and one where deep layers are crossed at column lines (long characteristics); exactness of I+ for linear S, convergence order, agreement between columns, Psi |
| `test_formal_3d` | the same for the 3-D solvers (`ShortChar`, linear and Bezier, bilinear and bicubic horizontal interpolation) |
| `test_formal_stokes` | the polarized 1-D DELO-Bezier solver with a constant (and a depth-dependent) propagation matrix and an unpolarized polynomial source function, against the matrix form of the analytic solution (`stokes_exact.h`): equal to the unpolarized Bezier solver when K' = 0, exactness of I+ for linear S, fourth order for the emergent Stokes vector. `test_formal_2d` and `test_formal_3d` check the 2-D and 3-D polarized solvers the same way on their slabs |
| `test_zeeman` | Zeeman splitting and polarized profiles of fictional lines (normal triplet, explicit g_eff, the anomalous patterns of Na D1 and D2, 5P2-5D2, 5P1-5D0, a g = 0 null line, a rejected J = 0-0 line): components, shifts and strengths against 3j symbols, g_eff, and phi_I,Q,U,V and psi_Q,U,V for several fields, two rays and both directions against the propagation matrix derived from the normal modes (RH writes V with the opposite sign to IAU). A Milne-Eddington slab built from RH's profiles and solved with the polarized solver against the Unno-Rachkovsky solution, and the weak-field relation for V |
| `test_formal_sphere` | spherical Feautrier along the rays of `getRays()` in a shell of constant opacity and S = B(T), against the closed-form intensities on shell rays (reflection at the tangent point) and core rays (thermalized core): convergence of P and the emergent intensity; grazing rays near the top bounded separately |

The other tests also take `--table` to print the measured errors. Error bounds in the tests are about twice what is measured, so that other
compilers and FMA contraction (Apple silicon) pass, while a real change in
accuracy fails.

## End-to-end tests

```
ctest --test-dir build -L e2e --output-on-failure
tests/e2e/stokes_1d.py --rhf1d build/rhf1d/rhf1d --work DIR -v   # keep the runs, print the errors
```

The scripts share a standard-library-only package, `tests/rhtest/` (no numpy:
the CI build jobs have none), which a script finds with a two-line
`sys.path` insertion of `tests/`, so that it also runs by hand:

| Module | What it holds |
|---|---|
| `cli.py` | `parse()`: `--PROGRAM` paths, `--src` (default: this source tree), `--work`, `-v` |
| `workdir.py` | `WorkDir`: the work directory with the Atoms/Atmos/Molecules links; `--work` is kept, a temporary directory under `$TMPDIR` is removed when the test passes and kept and named on a failure, `sys.exit` or an uncaught exception (`relink=True` replaces existing links, `announce=False` stops the message: `run_tests.py`) |
| `keywords.py` | `setkw()` (value `None` comments the key out), `unsetkw()` (every setting commented out, the line kept), `apply()`, `edit_file()`, `atoms_input()`, the `H6_CAII` atom set, `has_stokes_input()` |
| `run.py` | `prepare()` a run directory from `<geom>/run`, `copy_dir()` (a failed removal of an old directory raises), `execute()` with a timeout and `out.log`, `execute_streaming()` (line by line, with a callback, e.g. to create `STOP_RH`), `make_field()`, `Result.require()` (exit status and error text; 1 = not converged is expected after a few iterations), `RAN`, `EXIT_*`; for the regression suite `execute_to_log()` (output straight to a log, run under `/usr/bin/time -v -o time.txt` as the shell did, status as a shell reports it: 128 + signal), `peak_rss_mb()` (from `time.txt`), `completed()` (0-63), `log_problem()` (the reason a log shows a failed run: fatal error, sanitizer report, NaN) |
| `output.py` | readers (arrays of doubles come back as `array('d')` in native order, which numpy wraps with `np.frombuffer` without a copy) of `pops.*.out` (`read_pops_array(path, lte=True)`: the LTE populations), `atmos.out`, the first molecule of `molecules.out`, `J.dat`, `spectrum.out` (polarized or not; `read_spectrum_file` finds the layout) and the solveray output; `xdr_counted_string` (the binary formats live here) |
| `compare.py` | `compare_spectra()`: two `spectrum.out` files compared per wavelength relative to the peak, or byte-identical, without numpy; `compare_spectra.py` is a command line over it |
| `tree.py` | `tree_bin()`, `tree_src()`, `tree_sanitized()`, `missing_programs()`: where the programs and data of a build are (a CMake build directory, a source tree with `build/`, or one built in place with the old makefiles); `tree.sh` is the same for the shell script that still sources it (`make_references.sh`) |
| `check.py` | `Checker` (record, summary, exit status), `maxabs`, `maxdiff`, `maxrel`, `same_file` |

`tests/unit/test_rhtest.py` (ctest `unit_rhtest`, label `unit`) checks the readers (pops with the LTE block, `spectrum.out` of both layouts and the Q-21 size trap, solveray output), `compare_spectra()`, the command line of `compare_spectra.py` and the helpers of `run_tests.py` on synthetic files written with `struct`.

A new e2e script: `cli.parse`, `WorkDir`, `Checker`, then `prepare` and
`execute(...).require()` for each run, and `ck.finish(name, work)` last.

| Test | What it covers |
|---|---|
| `stokes_1d.py` | `rhf1d` with a fictional atom (`TI_fict.atom`: two lines with anomalous Zeeman patterns, g_eff 4/3 and 7/6) in LTE, one formal solution with J = 0 in the background scattering, so every run solves the same linear problem: no field gives Q = U = V = 0; parity (mirror in the plane of the rays) gives (I, Q, -U, -V) to 1e-12, with and without magneto-optical effects; without them B -> -B gives (I, Q, U, -V); turning the field azimuth for a vertical ray turns (Q, U), up to a discretization error from the polarized solver's componentwise limiter that must vanish at third order; the weak-field V against the derivative of I with respect to a uniform velocity (no wavelength derivative), which must converge at fourth order on a grey atmosphere. Python 3 standard library only, about 10 s |
| `prd_hybrid.py` | Ca II PRD, three iterations: `rhf1d` on FAL-C with one and three rays, `rhsc2d` on three FAL-C columns with one and six rays. In a static atmosphere the gas-frame mean intensity of hybrid PRD (`PRD_ANGLE_APPROX`) must equal J, and hybrid J must equal angle-averaged J, to 1e-12 and 1e-9 (B-31: with one ray Jgas was never reset; in 2-D and 3-D it was 2 J at wavelengths without PRD lines). About 20 s |
| `solveray_cosines.py` | `solveray` of `rhsc2d` and `rhsc3d` after a short FAL-C run with a field: the third direction cosine comes from the other two (B-34). 2-D: a pair with mux^2 + muz^2 = 1 + 6e-11 (mux written to 10 decimals) gives finite I, Q, U, V, equal to the ray with mux rounded down to 2e-4 of the peak (muy = 0 against 5e-6); mux = muz = 0.9 is an error that names the sum. 3-D (no main run: the checks stop before its results are read): sums of squares above 1 and equal to 1 (muz = 0) are errors. Before the fix the 2-D polarized output was all NaN and 3-D crashed. Python 3 standard library only, about 3 s |

Derivations: `docs/unit-tests`, section "End to end".

# Regression tests

Compare a *candidate* build against an unmodified *baseline* build on small,
fast runs (the whole suite takes a few minutes).

```
cmake -B build && cmake --build build -j                  # the candidate
cmake -B build-asan -DRH_SANITIZE=address,undefined && cmake --build build-asan -j
python3 tests/run_tests.py --baseline /path/to/upstream-tree --candidate . [--asan build-asan] [--exact]
```

CI compares every change with a build of its base commit (the target
branch before the change), in tolerance mode; results must match. A pull
request that changes results on purpose gets the label `changes-results`: then the
differences from the base are reported in the job summary instead of
failing (sanitizer, thread-count and memory checks still fail), and the PR
states what changed and by how much. CI also reports, without failing, the
differences from upstream RH at c95bdf6 (built without Fortran by
`tests/build_upstream.sh`): upstream is the reference we trust but check,
so every difference there needs an explanation in the log, either a known
numerical change or a fixed bug. On a pull request this report runs
only with the label `changes-results` or when the comparison with the base
fails: otherwise the results match the base, and so does their drift.
Pushes to `modernize` always report it. CI runs `2d/mem` on Linux only
(`--skip 2d/mem` on macOS, whose peak RSS is noisier). `--exact` is for refactors, with a build of
the branch's base (e.g. `modernize`) as the baseline. An upstream baseline
is built with:

```
git clone https://github.com/han-uitenbroek/RH upstream && git -C upstream checkout c95bdf6
tests/build_upstream.sh upstream
```

A build is given as a source tree with a CMake build in `build/`, as a CMake
build directory, or as a source tree built in place with the old makefiles
(upstream RH, the usual baseline). `rhf1d`, `rhsc2d`, `rhsc3d` and
`rhsphere` must be built. Run `tests/run_tests.py --help` for the options.

`run_tests.py` needs Python 3.8 or newer (CI has 3.12) and uses only the standard
library, with the `rhtest` package like the end-to-end tests. Programs run with their output in
`log.log` under `/usr/bin/time` (`-v` on Linux, `-l` on macOS), whose `time.txt`
the memory case reads. Python's own `wait4` would not do: a child's
`ru_maxrss` starts at the resident size of its parent when it is forked, so
the suite's own memory would be a floor under every measurement.

| Case | What it covers |
|---|---|
| `1d/nostokes`, `fieldfree`, `fullstokes` | all Stokes modes (FIELD_FREE iterates in 1 record then solves in 4) |
| `1d/crd` | no PRD lines |
| `1d/twoactive` | two active atoms |
| `2d/small`, `3d/small`, `sph/small` | the other geometries, from the shipped examples |
| `2d/fieldfree`, `3d/fieldfree` | the polarized 2-D and 3-D formal solvers, with the long characteristics, on the `small` atmospheres with a field that varies across the columns (written by `make_field.py` into the work directory): field-free iterations, then polarized final solutions (`N_MAX_SCATTER`) |
| `2d/lcbezier` | `2d/fieldfree` with Bézier long characteristics (`S_INTERPOLATION_LONGCHAR = S_BEZIER3`), unpolarized in the iterations and polarized in the final solutions |
| `2d/fullstokes` | `2d/fieldfree` with the polarized formal solution inside the iteration (`STOKES_MODE = FULL_STOKES`); diverged before B-55 |
| `2d/mem` | peak memory on a tiled grid of 96 columns (buffers must not scale with Nspect) |
| `1d/kurucz_blank` | Kurucz list with a leading blank line (currently crashes upstream) |

Each case checks the candidate at 1 thread against the baseline, the candidate
at other thread counts (3 and 8 by default) against its own 1-thread result,
and (with `--asan`) that AddressSanitizer is clean at 1 thread and the
largest count (3-D at the largest count only: the code that runs only with
one thread is shared by all geometries). Outputs are compared per wavelength, relative
to the peak. The defaults are based on the noise of unmodified code against
itself (about 1e-11 in 1-D/3-D/spherical, 2e-9 in the ill-conditioned 2-D
example), since floating-point sums depend on thread count. Use `--exact`
for changes that must not alter a single bit (memory management).

The log, `keyword.input` and timing of the runs behind every failure (both
sides of a comparison, in subdirectories named like the run directories; for
a failure raised as an error, every run of its case) are kept in
`WORK/failed/` (CI uploads them as an artifact), since run directories are
reused by later cases.

Running the suite with the same unmodified tree as baseline and candidate is
a useful self-test: everything passes except `1d/kurucz_blank`, which exposes
the upstream segfault.

## Stress, convergence and performance cases

`run_tests.py` stops every case after 4 iterations. `run_stress.py` runs
strongly NLTE problems **to convergence** (`--set conv`, 22 cases, 6 more on Bifrost cuts) or for a
fixed number of iterations on large multi-D grids (`--set perf`), on one or
more builds, and reports iterations, final population change, trend
(monotone / oscillating / stalled), wall time, peak memory, and the population
and spectrum differences between builds.

```
python3 tests/run_stress.py --tree upstream=/path/to/upstream --tree candidate=. \
    --set conv --jobs 8 --work /large/scratch/rh-stress
python3 tests/run_stress.py --list          # the cases and what each one exposes
```

The cases include known failures (non-convergence, crashes, wrong answers);
the case table in the script groups them by what they expose. 1-D cases run on
1 thread so that builds compare without summation-order noise. Multi-D
atmospheres are generated by `make_perturbed_atmos.py` (FAL-C with a
corrugated transition region and velocity field); `atmos/` holds two dynamic
1-D models. `run_stress.py` needs Python 3.8 or newer and the standard library only (`analyze_run.py`, which it calls, needs numpy). It needs GNU
`/usr/bin/time`. `--only` takes a `grep -E` pattern and is matched by `grep`. `analyze_run.py` does the log parsing and population comparison; it and
`compare_spectra.py` read the binary outputs with `rhtest.output` and keep numpy for the analysis.
`se_residual.py RUNDIR` checks the end point of a run in any geometry: the
residual of the plain rate equations, with the radiative rates RH wrote for
the final formal solution, at every point and level. It should be of order
`ITER_LIMIT` after convergence (needs `RADRATE_OUTPUT`, `COLLRATE_OUTPUT` and
`STOKES_MODE = NO_STOKES`).
Use a `--work` directory on a large filesystem: PRD runs write GB-sized
scratch files.

### Realistic atmospheres: Bifrost

The `bf_` cases use cuts of a Bifrost snapshot (Carlsson et al. 2016, A&A
585, A4: en024048_hion, snapshot 385, the FITS cubes of lgtg, lgne, ux, uy,
uz, bx, by, bz and lgn1..lgn6). The data are not distributed with RH and these
cases do not run on GitHub: download the cubes and point the script at them,

```
python3 tests/run_stress.py --tree candidate=. --bifrost /path/to/en024048_hion/atmos --only _bf_
```

(or set `RH_BIFROST_DIR`); without it they are skipped. `make_bifrost_atmos.py`
cuts the snapshot to z = -0.5 ... 3 Mm (every second point, 89 heights) as a
1-D column, a 2-D row or a 3-D box. RH's horizontal boundaries are periodic
(3-D allows nothing else; fixed 2-D boundaries need the incoming intensity at
every wavelength), and so is the snapshot: a full row or the whole box sampled
every S-th point, S a divisor of 504, stays periodic; a window does not and is
used for timing only. Ca II runs in CRD: with PRD the iteration diverges in
these moving atmospheres (B-29), which the `_prd` cases reproduce.
