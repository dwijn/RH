# RH tests

Four tiers:

- **Unit tests** (`tests/unit/`, run by `ctest`): single routines on
  problems with known answers. Fast, no input files.
- **End-to-end tests** (`tests/e2e/`, run by `ctest`): whole runs of a
  program on problems with exact properties.
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
| `test_interp` | `Hunt`/`Locate` (every starting guess, ascending and descending tables, values at, between and beyond the nodes), `Linear`, natural and exponential splines (exact for linear data; mpmath reference solutions; convergence order; tension limits), `cc_kernel`/`cubeconvol` (kernel identities, exact for biquadratics including boundary cells, order 3) |
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

| Test | What it covers |
|---|---|
| `stokes_1d.py` | `rhf1d` with a fictional atom (`TI_fict.atom`: two lines with anomalous Zeeman patterns, g_eff 4/3 and 7/6) in LTE, one formal solution with J = 0 in the background scattering, so every run solves the same linear problem: no field gives Q = U = V = 0; parity (mirror in the plane of the rays) gives (I, Q, -U, -V) to 1e-12, with and without magneto-optical effects; without them B -> -B gives (I, Q, U, -V); turning the field azimuth for a vertical ray turns (Q, U), up to a discretization error from the polarized solver's componentwise limiter that must vanish at third order; the weak-field V against the derivative of I with respect to a uniform velocity (no wavelength derivative), which must converge at fourth order on a grey atmosphere. Python 3 standard library only, about 10 s |
| `prd_hybrid.py` | Ca II PRD, three iterations: `rhf1d` on FAL-C with one and three rays, `rhsc2d` on three FAL-C columns with one and six rays. In a static atmosphere the gas-frame mean intensity of hybrid PRD (`PRD_ANGLE_APPROX`) must equal J, and hybrid J must equal angle-averaged J, to 1e-12 and 1e-9 (B-31: with one ray Jgas was never reset; in 2-D and 3-D it was 2 J at wavelengths without PRD lines). About 20 s |

Derivations: `docs/unit-tests`, section "End to end".

# Regression tests

Compare a *candidate* build against an unmodified *baseline* build on small,
fast runs (the whole suite takes a few minutes).

```
cmake -B build && cmake --build build -j                  # the candidate
cmake -B build-asan -DRH_SANITIZE=address,undefined && cmake --build build-asan -j
tests/run_tests.sh --baseline /path/to/upstream-tree --candidate . [--asan build-asan] [--exact]
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
numerical change or a fixed bug. `--exact` is for refactors, with a build of
the branch's base (e.g. `modernize`) as the baseline. An upstream baseline
is built with:

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
strongly NLTE problems **to convergence** (`--set conv`, 21 cases, 5 more on Bifrost cuts) or for a
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

### Realistic atmospheres: Bifrost

The `bf_` cases use cuts of a Bifrost snapshot (Carlsson et al. 2016, A&A
585, A4: en024048_hion, snapshot 385, the FITS cubes of lgtg, lgne, ux, uy,
uz, bx, by, bz and lgn1..lgn6). The data are not distributed with RH and these
cases do not run on GitHub: download the cubes and point the script at them,

```
tests/run_stress.sh --tree candidate=. --bifrost /path/to/en024048_hion/atmos --only _bf_
```

(or set `RH_BIFROST_DIR`); without it they are skipped. `make_bifrost_atmos.py`
cuts the snapshot to z = -0.5 ... 3 Mm (every second point, 89 heights) as a
1-D column, a 2-D row or a 3-D box. RH's horizontal boundaries are periodic
(3-D allows nothing else; fixed 2-D boundaries need the incoming intensity at
every wavelength), and so is the snapshot: a full row or the whole box sampled
every S-th point, S a divisor of 504, stays periodic; a window does not and is
used for timing only. Ca II runs in CRD: with PRD the iteration diverges in
these moving atmospheres (B-29), which the `_prd` cases reproduce.
