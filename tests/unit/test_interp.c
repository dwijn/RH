/* ------- file: -------------------------- test_interp.c -----------

	   Table look-up and interpolation: hunt.c (Hunt, Locate), linear.c
	   (Linear), spline.c, expspline.c (natural cubic and
	   exponential splines) and cubeconvol.c (cubic convolution).

	   -- Hunt, Locate: the bracketing index i with table[i] <= value <
		  table[i+1] for all values at, between, below and above the
		  nodes, for ascending and descending tables, tables of 1, 2, 3
		  and many elements, and every possible starting guess of Hunt;
	   -- Linear: exact for linear data, the end
		  values outside the table, error of at most h^2 max|f''|/8 and
		  order 2 for a smooth function;
	   -- splineCoef/splineEval, exp_splineCoef/exp_splineEval: exact for
		  linear data, the natural cubic spline and the exponential
		  spline (two tensions) against the solution of their defining
		  conditions with mpmath (kernels_reference.h, from
		  make_kernels_reference.py), the same for a descending table,
		  order 4 for sin(x) on [0, pi] (f'' = 0 at both ends, as the
		  natural boundary conditions assume), order 4 in the interior
		  but only 2 overall for exp(x), the limits of small and large
		  tension;
	   -- cc_kernel, cubeconvol: the kernel's partition of unity and
		  reproduction of quadratics, exact for biquadratic data
		  everywhere (also in the boundary cells where the data are
		  extrapolated), the end values outside the table and order 3 for
		  a smooth function.

	   Run with --table to print the measured errors.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "check.h"
#include "kernels_reference.h"

static bool_t table;

/* --- Hunt and Locate -- ------------------------------------------- */

/* --- The index the bracketing must give: for an ascending table the
	   last i with a[i] <= v, for a descending table the last i with
	   a[i] >= v, and 0 when there is none -- ------------------- */

static int expected_index(int n, const double *a, double v, bool_t ascend)
{
	int i = 0;

	for (int k = 0; k < n; k++)
		if (ascend ? (a[k] <= v) : (a[k] >= v)) i = k;
	return i;
}

static void check_bracketing(int n, const double *a)
{
	bool_t ascend = (n > 1 && a[n - 1] > a[0]);
	double lo     = (a[0] < a[n - 1]) ? a[0] : a[n - 1];
	double hi     = (a[0] < a[n - 1]) ? a[n - 1] : a[0];
	double values[3 * 64 + 4];
	int    nv = 0;

	values[nv++] = lo - 1.0;
	values[nv++] = lo - 1e-9;
	values[nv++] = hi + 1e-9;
	values[nv++] = hi + 7.0;
	for (int k = 0; k < n; k++) {
		values[nv++] = a[k];
		if (k > 0) values[nv++] = 0.5 * (a[k - 1] + a[k]);
		if (k > 0) values[nv++] = a[k - 1] + 0.999 * (a[k] - a[k - 1]);
	}
	for (int iv = 0; iv < nv; iv++) {
		int want = expected_index(n, a, values[iv], ascend), ilow;

		Locate(n, (double *) a, values[iv], &ilow);
		CHECK(ilow == want, "Locate: n = %d, value %.17g: got %d, want %d",
			  n, values[iv], ilow, want);

		/* --- Every starting guess, in and out of range -- ------- */

		for (int guess = -2; guess <= n + 2; guess++) {
			ilow = guess;
			Hunt(n, (double *) a, values[iv], &ilow);

			/* --- Hunt must agree with Locate and with the bracketing
				   convention, also when value equals a node exactly -- */

			CHECK(ilow == want, "Hunt: n = %d, value %.17g, guess %d: got %d, want %d",
				  n, values[iv], guess, ilow, want);
		}
	}
}

static void test_hunt_locate(void)
{
	static const int sizes[] = {1, 2, 3, 4, 7, 16, 37, 64};

	for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
		int    n = sizes[s];
		double up[64], down[64];

		/* --- Irregular spacing, so that midpoint and bisection
			   differ -- ---------------------------------------- */

		for (int k = 0; k < n; k++) up[k] = -3.0 + 0.5 * k + 0.37 * k * k / n;
		for (int k = 0; k < n; k++) down[k] = up[n - 1 - k];
		check_bracketing(n, up);
		check_bracketing(n, down);
	}

	/* --- Repeated hunting with the previous result as the guess,
		   as the callers do (random access, then monotonic) -- --- */

	{
		double a[50];
		int    ilow = 0;

		for (int k = 0; k < 50; k++) a[k] = k * k * 0.1;
		for (int k = 0; k < 400; k++) {
			double v    = 0.1 * ((k * 7919) % 2600) - 5.0;
			int    want = expected_index(50, a, v, TRUE);

			Hunt(50, a, v, &ilow);
			CHECK(ilow == want, "hunt sequence, value %g: got %d, want %d", v, ilow, want);
		}
		for (int k = 0; k < 600; k++) {
			double v    = 0.0043 * k * k / 6.0;
			int    want = expected_index(50, a, v, TRUE);

			Hunt(50, a, v, &ilow);
			CHECK(ilow == want, "monotonic hunt, value %g: got %d, want %d", v, ilow, want);
		}
	}
}

/* --- Linear -- --------------------------------------- */

static void test_linear(void)
{
	enum { N = 9,
		   M = 25 };
	double xt[N], yt[N], xd[N], yd[N], x[M], y[M], f;

	for (int k = 0; k < N; k++) {
		xt[k]         = 1.0 + k * k * 0.25 + 0.3 * k;
		yt[k]         = 3.0 * xt[k] - 2.0;
		xd[N - 1 - k] = xt[k];
		yd[N - 1 - k] = yt[k];
	}
	for (int k = 0; k < M; k++) x[k] = xt[0] - 1.0 + (xt[N - 1] - xt[0] + 2.0) * k / (M - 1.0);

	for (int hunt = 0; hunt < 2; hunt++) {
		for (int desc = 0; desc < 2; desc++) {
			Linear(N, desc ? xd : xt, desc ? yd : yt, M, x, y, hunt);
			for (int k = 0; k < M; k++) {
				/* --- linear data in the table, constant outside -- */

				f = 3.0 * MIN(MAX(x[k], xt[0]), xt[N - 1]) - 2.0;
				CHECK_CLOSE(y[k], f, 1e-14, 1e-14,
							"linear data, hunt = %d, desc = %d, x = %g", hunt, desc, x[k]);
			}
		}
	}

	/* --- Order 2: |error| <= h^2 max|f''| / 8 for f = sin(x) on [0, 3],
		   and the error shrinks by 4 when h is halved (measured
		   ratios 3.99 ..4.00) -- --------------------------------- */

	{
		double prev = 0.0;

		for (int n = 9; n <= 129; n = 2 * n - 1) {
			double xs[129], ys[129], xe[400], ye[400], err = 0.0,
													   h = 3.0 / (n - 1);

			for (int k = 0; k < n; k++) {
				xs[k] = 3.0 * k / (n - 1.0);
				ys[k] = sin(xs[k]);
			}
			for (int k = 0; k < 400; k++) xe[k] = 3.0 * (k + 0.5) / 400.0;
			Linear(n, xs, ys, 400, xe, ye, TRUE);
			for (int k = 0; k < 400; k++) err = MAX(err, fabs(ye[k] - sin(xe[k])));
			CHECK(err <= 1.0001 * h * h / 8.0, "n = %d: error %.3g exceeds h^2/8 = %.3g", n, err, h * h / 8.0);
			if (prev > 0.0)
				CHECK(prev / err > 3.8 && prev / err < 4.2,
					  "n = %d: error ratio %.3g, want ~4", n, prev / err);
			if (table) printf("Linear n = %3d: max error %.3e\n", n, err);
			prev = err;
		}
	}
}

/* --- Cubic splines -- --------------------------------------------- */

/* --- Both splines keep static state (issue T-3) and pointers to the
	   table, so the coefficients are always set just before the
	   evaluation, with the table still alive -- ----------------- */

static double max_rel_err(int n, const double *got, const double *want, double scale)
{
	double m = 0.0;

	for (int k = 0; k < n; k++) m = MAX(m, fabs(got[k] - want[k]) / scale);
	return m;
}

static void test_spline_linear_data(void)
{
	enum { N = 12,
		   M = 40 };
	double x[N], y[N], xd[N], yd[N], xe[M], ye[M];

	for (int k = 0; k < N; k++) {
		x[k]          = 0.5 + 0.1 * k * k + 0.2 * k;
		y[k]          = 5.0 - 2.0 * x[k];
		xd[N - 1 - k] = x[k];
		yd[N - 1 - k] = y[k];
	}
	for (int k = 0; k < M; k++) xe[k] = x[0] - 0.3 + (x[N - 1] - x[0] + 0.6) * k / (M - 1.0);

	for (int kind = 0; kind < 3; kind++) {     /* cubic, exp 0.7, exp 4 */
		for (int desc = 0; desc < 2; desc++) { /* ascending, descending */
			for (int hunt = 0; hunt < 2; hunt++) {
				double *xt = desc ? xd : x, *yt = desc ? yd : y, maxerr = 0.0;

				if (kind == 0) {
					splineCoef(N, xt, yt);
					splineEval(M, xe, ye, hunt);
				} else {
					exp_splineCoef(N, xt, yt, (kind == 1) ? 0.7 : 4.0);
					exp_splineEval(M, xe, ye, hunt);
				}
				for (int k = 0; k < M; k++) {
					double xc = MIN(MAX(xe[k], x[0]), x[N - 1]);

					maxerr = MAX(maxerr, fabs(ye[k] - (5.0 - 2.0 * xc)));
				}
				CHECK(maxerr < 8e-15, "kind %d, desc %d, hunt %d: max error %.3g", kind, desc, hunt, maxerr);
				if (table && !hunt && !desc)
					printf("spline kind %d, linear data: max error %.2e\n", kind, maxerr);
			}
		}
	}
}

/* --- Measured (gcc, x86-64), relative to max |y| = e: 1.6e-16 for all
	   (the error of the evaluation, 1 ulp of y); bound 4e-16. Linear data
	   (above): 3.6e-15 for y up to 25, bound 8e-15 --              - */

static void test_spline_reference(void)
{
	double xr[SPLINE_REF_N], yr[SPLINE_REF_N], ye[SPLINE_REF_NE], scale = exp(1.0);

	for (int desc = 0; desc < 2; desc++) {
		for (int k = 0; k < SPLINE_REF_N; k++) {
			xr[k] = desc ? spline_ref_x[SPLINE_REF_N - 1 - k] : spline_ref_x[k];
			yr[k] = desc ? spline_ref_y[SPLINE_REF_N - 1 - k] : spline_ref_y[k];
		}
		for (int hunt = 0; hunt < 2; hunt++) {
			double err;

			splineCoef(SPLINE_REF_N, xr, yr);
			splineEval(SPLINE_REF_NE, (double *) spline_ref_xe, ye, hunt);
			err = max_rel_err(SPLINE_REF_NE, ye, spline_ref_natural, scale);
			CHECK(err < 4e-16, "natural spline, desc %d, hunt %d: error %.3g", desc, hunt, err);
			if (table) printf("natural spline, desc %d hunt %d: max error %.2e\n", desc, hunt, err);

			for (int it = 0; it < EXPSPLINE_REF_NT; it++) {
				exp_splineCoef(SPLINE_REF_N, xr, yr, expspline_ref_tension[it]);
				exp_splineEval(SPLINE_REF_NE, (double *) spline_ref_xe, ye, hunt);
				err = max_rel_err(SPLINE_REF_NE, ye, expspline_ref_y[it], scale);
				CHECK(err < 4e-16, "exp spline, tension %g, desc %d, hunt %d: error %.3g",
					  expspline_ref_tension[it], desc, hunt, err);
				if (table)
					printf("exp spline, tension %g, desc %d hunt %d: max error %.2e\n",
						   expspline_ref_tension[it], desc, hunt, err);
			}
		}
	}

	/* --- The spline passes through the nodes, and is constant (the
		   end value) outside the table -- ------------------------ */

	{
		double xo[4] = {spline_ref_x[0] - 1.0, spline_ref_x[0], spline_ref_x[SPLINE_REF_N - 1],
						spline_ref_x[SPLINE_REF_N - 1] + 5.0},
			   yo[4];

		splineCoef(SPLINE_REF_N, (double *) spline_ref_x, (double *) spline_ref_y);
		splineEval(SPLINE_REF_N, (double *) spline_ref_x, ye, FALSE);
		for (int k = 1; k < SPLINE_REF_N - 1; k++)
			CHECK_CLOSE(ye[k], spline_ref_y[k], 1e-15, 0.0, "node %d", k);
		splineEval(4, xo, yo, TRUE);
		CHECK(yo[0] == spline_ref_y[0] && yo[1] == spline_ref_y[0] &&
				  yo[2] == spline_ref_y[SPLINE_REF_N - 1] && yo[3] == spline_ref_y[SPLINE_REF_N - 1],
			  "end values: %g %g %g %g", yo[0], yo[1], yo[2], yo[3]);
	}
}

/* --- Order of convergence of the natural spline -- ---------------- */

static double spline_error(int n, double (*f)(double), double x0, double x1, bool_t one_point)
{
	double xs[130], ys[130], xe[1000], ye[1000], err = 0.0;
	int    ne = 1000;

	for (int k = 0; k < n; k++) {
		xs[k] = x0 + (x1 - x0) * k / (n - 1.0);
		ys[k] = f(xs[k]);
	}
	if (one_point) {
		/* --- the middle of the cell just below the center of the
			   interval: the same position in the cell for every n -- */

		xe[0] = x0 + (x1 - x0) * (0.5 - 0.5 / (n - 1.0));
		ne    = 1;
	} else {
		for (int k = 0; k < ne; k++) xe[k] = x0 + (x1 - x0) * (k + 0.5) / ne;
	}
	splineCoef(n, xs, ys);
	splineEval(ne, xe, ye, TRUE);
	for (int k = 0; k < ne; k++) err = MAX(err, fabs(ye[k] - f(xe[k])));
	return err;
}

/* --- sin(x) on [0, pi] has f'' = 0 at both ends, as the natural
	   boundary conditions assume, so the error is O(h^4) everywhere
	   (ratio 16). For exp(x) on [0, 1] the natural conditions are wrong
	   by O(1) at the ends: the maximum error is O(h^2) (ratio 4) and,
	   because the effect decays by a factor of 3.7 per interval,
	   O(h^4) away from the ends. Measured ratios: see --table --  - */

static void test_spline_order(void)
{
	double prev = 0.0, prev_point = 0.0, prev_exp = 0.0;

	for (int n = 9; n <= 129; n = 2 * n - 1) {
		double e_sin = spline_error(n, sin, 0.0, acos(-1.0), FALSE);
		double e_exp = spline_error(n, exp, 0.0, 1.0, FALSE);
		double e_pt  = spline_error(n, exp, 0.0, 1.0, TRUE);

		if (table)
			printf("spline n = %3d: sin %.3e, exp %.3e, exp at one point %.3e\n",
				   n, e_sin, e_exp, e_pt);
		if (prev > 0.0) {
			CHECK(prev / e_sin > 14.0 && prev / e_sin < 17.0, "sin, n = %d: error ratio %.3g, want ~16", n, prev / e_sin);
			CHECK(prev_exp / e_exp > 3.7 && prev_exp / e_exp < 4.3, "exp, n = %d: error ratio %.3g, want ~4", n, prev_exp / e_exp);
			CHECK(n < 33 || (prev_point / e_pt > 13.0 && prev_point / e_pt < 19.0),
				  "exp at one point, n = %d: error ratio %.3g, want ~16", n, prev_point / e_pt);
		}
		prev       = e_sin;
		prev_exp   = e_exp;
		prev_point = e_pt;
	}
}

/* --- Small tension: the cubic spline, the difference O(sigma^2).
	   Large tension: piecewise linear interpolation, the difference
	   decreasing with tension. Measured (relative to e): 3.4e-8 at
	   tension 0.01, and 1.4e-2, 6.0e-3, 1.2e-3 against Linear for
	   tensions 2, 10, 50; bounds twice that -- */

static void test_spline_limits(void)
{
	double ys[SPLINE_REF_NE], yl[SPLINE_REF_NE], yc[SPLINE_REF_NE];
	double e = exp(1.0);

	splineCoef(SPLINE_REF_N, (double *) spline_ref_x, (double *) spline_ref_y);
	splineEval(SPLINE_REF_NE, (double *) spline_ref_xe, yc, FALSE);
	Linear(SPLINE_REF_N, (double *) spline_ref_x, (double *) spline_ref_y,
		   SPLINE_REF_NE, (double *) spline_ref_xe, yl, FALSE);

	exp_splineCoef(SPLINE_REF_N, (double *) spline_ref_x, (double *) spline_ref_y, 1e-2);
	exp_splineEval(SPLINE_REF_NE, (double *) spline_ref_xe, ys, FALSE);
	{
		double d = max_rel_err(SPLINE_REF_NE, ys, yc, e);

		CHECK(d < 7e-8, "small tension: difference from the cubic spline %.3g", d);
		if (table) printf("tension 1e-2 vs cubic spline: %.2e\n", d);
	}
	{
		static const double tension[3] = {2.0, 10.0, 50.0};
		static const double bound[3]   = {3e-2, 1.2e-2, 2.5e-3};
		double              prev       = 1.0;

		for (int it = 0; it < 3; it++) {
			double d;

			exp_splineCoef(SPLINE_REF_N, (double *) spline_ref_x, (double *) spline_ref_y, tension[it]);
			exp_splineEval(SPLINE_REF_NE, (double *) spline_ref_xe, ys, FALSE);
			d = max_rel_err(SPLINE_REF_NE, ys, yl, e);
			CHECK(d < bound[it] && d < prev, "tension %g: difference from linear %.3g (bound %.3g, previous %.3g)",
				  tension[it], d, bound[it], prev);
			if (table) printf("tension %g vs linear: %.2e\n", tension[it], d);
			prev = d;
		}
	}
}

/* --- Cubic convolution -- ----------------------------------------- */

static void test_cc_kernel(void)
{
	double u[4];

	/* --- Partition of unity and reproduction of 1, s, s^2 at the
		   positions -1, 0, 1, 2 of the four points -- ----------- */

	for (int k = 0; k <= 20; k++) {
		double s = k / 20.0, s0 = 0.0, s1 = 0.0, s2 = 0.0;

		cc_kernel(s, u);
		for (int m = 0; m < 4; m++) {
			s0 += u[m];
			s1 += u[m] * (m - 1.0);
			s2 += u[m] * (m - 1.0) * (m - 1.0);
		}
		CHECK_CLOSE(s0, 1.0, 1e-15, 0.0, "sum, s = %g", s);
		CHECK_CLOSE(s1, s, 1e-15, 1e-15, "first moment, s = %g", s);
		CHECK_CLOSE(s2, s * s, 1e-15, 1e-15, "second moment, s = %g", s);
	}
	cc_kernel(0.0, u);
	CHECK(u[0] == 0.0 && u[1] == 1.0 && u[2] == 0.0 && u[3] == 0.0, "kernel at 0: %g %g %g %g", u[0], u[1], u[2], u[3]);
	cc_kernel(1.0, u);
	CHECK(u[0] == 0.0 && u[1] == 0.0 && u[2] == 1.0 && u[3] == 0.0, "kernel at 1: %g %g %g %g", u[0], u[1], u[2], u[3]);
	cc_kernel(0.5, u);
	CHECK_CLOSE(u[0], -1.0 / 16.0, 1e-15, 0.0, "kernel at 1/2, m = %d", 0);
	CHECK_CLOSE(u[1], 9.0 / 16.0, 1e-15, 0.0, "kernel at 1/2, m = %d", 1);
	CHECK_CLOSE(u[2], 9.0 / 16.0, 1e-15, 0.0, "kernel at 1/2, m = %d", 2);
	CHECK_CLOSE(u[3], -1.0 / 16.0, 1e-15, 0.0, "kernel at 1/2, m = %d", 3);
}

/* --- f[j * Nx + i] with x along i, y along j --------------------- */

static double biquad(double x, double y)
{
	return 1.5 + 0.8 * x - 0.3 * y + 0.25 * x * x + 0.4 * y * y - 0.7 * x * y +
		   0.15 * x * x * y - 0.2 * x * y * y + 0.05 * x * x * y * y;
}

static void test_cubeconvol_exact(void)
{
	enum { NX = 7,
		   NY = 5 };
	double f[NX * NY], maxerr = 0.0, scale = 0.0;

	for (int j = 0; j < NY; j++)
		for (int i = 0; i < NX; i++) f[j * NX + i] = biquad(i, j);

	/* --- the nodes: in the interior, at edges and corners -- ---- */

	for (int j = 0; j < NY; j++)
		for (int i = 0; i < NX; i++) {
			double v = cubeconvol(NX, NY, f, i, j);

			CHECK_CLOSE(v, f[j * NX + i], 1e-14, 1e-14, "node %d %d", i, j);
		}

	/* --- biquadratic data at 0.1 spacing, including the boundary
		   cells where the missing points are extrapolated -- ------ */

	for (int jj = 0; jj <= 10 * (NY - 1); jj++)
		for (int ii = 0; ii <= 10 * (NX - 1); ii++) {
			double x = 0.1 * ii, y = 0.1 * jj, want = biquad(x, y);
			double v = cubeconvol(NX, NY, f, x, y);

			maxerr = MAX(maxerr, fabs(v - want));
			scale  = MAX(scale, fabs(want));
		}
	CHECK(maxerr < 1e-15 * scale, "biquadratic data: max error %.3g (scale %.3g)", maxerr, scale);
	if (table) printf("cubeconvol biquadratic data: max error %.2e (scale %.1f)\n", maxerr, scale);

	/* --- outside the table the value at the nearest edge -- ----- */

	CHECK_CLOSE(cubeconvol(NX, NY, f, -2.0, 1.5), cubeconvol(NX, NY, f, 0.0, 1.5), 1e-14, 0.0, "x below, code %d", 0);
	CHECK_CLOSE(cubeconvol(NX, NY, f, NX + 3.0, 1.5), cubeconvol(NX, NY, f, NX - 1.0, 1.5), 1e-14, 0.0, "x above, code %d", 1);
	CHECK_CLOSE(cubeconvol(NX, NY, f, 2.5, -0.1), cubeconvol(NX, NY, f, 2.5, 0.0), 1e-14, 0.0, "y below, code %d", 2);
	CHECK_CLOSE(cubeconvol(NX, NY, f, 2.5, NY + 0.1), cubeconvol(NX, NY, f, 2.5, NY - 1.0), 1e-14, 0.0, "y above, code %d", 3);
	CHECK_CLOSE(cubeconvol(NX, NY, f, -1.0, 1e3), f[(NY - 1) * NX + 0], 1e-14, 1e-14, "corner, code %d", 4);
}

static double smooth(double x, double y)
{
	return sin(1.3 * x + 0.4) * cos(0.9 * y - 0.2) + 0.1 * x * y;
}

/* --- Grid points at x = i h on [0, 2] x [0, 2]; order 3 (ratio of
	   the maximum error 8; measured 7.7, 8.3, 9.5) --                                      - */

static void test_cubeconvol_order(void)
{
	double prev = 0.0;

	for (int n = 11; n <= 81; n = 2 * n - 1) {
		double h = 2.0 / (n - 1), err = 0.0, *f = malloc((size_t) n * n * sizeof(double));

		for (int j = 0; j < n; j++)
			for (int i = 0; i < n; i++) f[j * n + i] = smooth(i * h, j * h);
		/* --- points at all phases within the cells (golden ratio
			   sequences), so that the maximum does not depend on n -- */

		for (int jj = 1; jj <= 60; jj++)
			for (int ii = 1; ii <= 60; ii++) {
				double xi = (n - 1) * fmod(ii * 0.6180339887498949, 1.0),
					   yi = (n - 1) * fmod(jj * 0.7548776662466927, 1.0);

				err = MAX(err, fabs(cubeconvol(n, n, f, xi, yi) - smooth(xi * h, yi * h)));
			}
		if (table) printf("cubeconvol n = %3d: max error %.3e\n", n, err);
		if (prev > 0.0)
			CHECK(prev / err > 6.5 && prev / err < 10.5, "n = %d: error ratio %.3g, want ~8", n, prev / err);
		prev = err;
		free(f);
	}
}

int main(int argc, char *argv[])
{
	table = (argc > 1 && strcmp(argv[1], "--table") == 0);

	test_hunt_locate();
	test_linear();
	test_spline_linear_data();
	test_spline_reference();
	test_spline_order();
	test_spline_limits();
	test_cc_kernel();
	test_cubeconvol_exact();
	test_cubeconvol_order();

	return check_summary("test_interp");
}
