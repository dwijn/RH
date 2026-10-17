/* ------- file: -------------------------- test_quadrature.c ------

	   Quadrature: gaussleg.c (GaussLeg) and the integration weights
	   of piecewise quadratic integration in w3.c (w2, w3, U3).

	   -- GaussLeg, n = 1 .. 64 and 100: nodes inside the interval, in
		  increasing order and symmetric, weights positive and summing
		  to the length of the interval, exact integration of x^k for
		  all k <= 2n-1 (on [-1, 1] and on [0.5, 3.2]); not exact for
		  x^(2n), where the error is known from the remainder term
		  2^(2n+1) (n!)^4 / ((2n+1) ((2n)!)^2); closed forms of the
		  nodes and weights for n = 2 .. 4; a reversed interval
		  (x1 > x2) gives the integral with the opposite sign;
	   -- w2, w3, U3 against w_k = gamma(k+1, dtau) (the weights
		  int_0^dtau t^k exp(-t) dt) and U_k = exp(-dtau) int_0^dtau
		  t^k exp(t) dt computed with mpmath (kernels_reference.h, from
		  make_kernels_reference.py) on dtau from 0 to 1000, spanning
		  the switches between the series for small dtau, the exact
		  expressions and the limits for large dtau.

	   Run with --table to print the measured errors.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "check.h"
#include "kernels_reference.h"

static bool_t table;

/* --- GaussLeg -- -------------------------------------------------- */

/* --- Integral of x^k over [a, b] -- ------------------------------- */

static double monomial_integral(int k, double a, double b)
{
	return (pow(b, k + 1.0) - pow(a, k + 1.0)) / (k + 1.0);
}

static void check_rule(int n, double a, double b, double *worst_poly)
{
	/* --- see test_gaussleg_rules for the bounds -- ------------ */

	double bound_nodes = 1.2e-15, bound_poly = (n <= 48) ? 1.3e-13 : 1.1e-12;
	double x[128], w[128], sum = 0.0, worst = 0.0;

	GaussLeg(a, b, x, w, n);

	for (int i = 0; i < n; i++) {
		CHECK(x[i] > a && x[i] < b, "n = %d: node %d = %.17g outside (%g, %g)", n, i, x[i], a, b);
		CHECK(w[i] > 0.0, "n = %d: weight %d = %.17g", n, i, w[i]);
		if (i > 0) CHECK(x[i] > x[i - 1], "n = %d: nodes %d, %d not increasing", n, i - 1, i);
		/* --- symmetric about the center -- ---------------------- */
		CHECK_CLOSE(x[i] - 0.5 * (a + b), -(x[n - 1 - i] - 0.5 * (a + b)), 0.0, bound_nodes,
					"n = %d: node symmetry, i = %d", n, i);
		CHECK_CLOSE(w[i], w[n - 1 - i], 0.0, 1e-15 * (b - a), "n = %d: weight symmetry, i = %d", n, i);
		sum += w[i];
	}
	CHECK_CLOSE(sum, b - a, 5e-14, 0.0, "n = %d: sum of the weights", n);

	/* --- Exact for degree <= 2n-1. The error relative to the integral
		   of |x|^k (which is what rounding is relative to) -- ----- */

	for (int k = 0; k <= 2 * n - 1; k++) {
		double q = 0.0, exact = monomial_integral(k, a, b), scale;

		for (int i = 0; i < n; i++) q += w[i] * pow(x[i], k);
		scale = (a < 0.0) ? monomial_integral(k + (k % 2), -1.0, 1.0) : fabs(exact);
		worst = MAX(worst, fabs(q - exact) / scale);
	}
	CHECK(worst < bound_poly, "n = %d on [%g, %g]: polynomials, worst relative error %.3g", n, a, b, worst);
	if (worst_poly) *worst_poly = MAX(*worst_poly, worst);
}

static void test_gaussleg_rules(void)
{
	static const int sizes[] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 32, 48, 64, 100};
	double           w1 = 0.0, w2 = 0.0;

	/* --- The weights are only accurate to ~1e-14 (relative), not to the
		   1e-16 of the arithmetic: they come out systematically too
		   small (the sum of the weights is below the length by up to
		   2.4e-14 at n = 64) and polynomial integrals are off by up to
		   6e-14 for n <= 48 and 5e-13 for n = 64 and 100 (relative to
		   the integral of |x|^k). Cause: gaussleg.c:49-56 computes the
		   weight from the derivative pp at the previous Newton iterate,
		   and iterates only until the step is below 3e-14. Bounds: about
		   twice the measured values; the nodes are symmetric to 1.2e-15 -- */

	for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
		check_rule(sizes[s], -1.0, 1.0, &w1);
		check_rule(sizes[s], 0.5, 3.2, &w2);
	}
	if (table) printf("GaussLeg: worst polynomial error %.2e on [-1,1], %.2e on [0.5,3.2]\n", w1, w2);
}

/* --- x^(2n) is not integrated exactly: on [-1, 1] the remainder
	   2^(2n+1) (n!)^4 / ((2n+1) ((2n)!)^2) is what is missing --- */

static void test_gaussleg_remainder(void)
{
	for (int n = 1; n <= 12; n++) {
		double x[16], w[16], q = 0.0, fact_n = 1.0, fact_2n = 1.0, rem;

		for (int i = 2; i <= n; i++) fact_n *= i;
		for (int i = 2; i <= 2 * n; i++) fact_2n *= i;
		rem = pow(2.0, 2.0 * n + 1.0) * pow(fact_n, 4.0) / ((2.0 * n + 1.0) * fact_2n * fact_2n);

		GaussLeg(-1.0, 1.0, x, w, n);
		for (int i = 0; i < n; i++) q += w[i] * pow(x[i], 2 * n);
		/* --- the quadrature gives less than the integral by rem;
			   the difference of two O(1) numbers, so the precision is
			   absolute: measured 7e-15 at most (see test_gaussleg_rules) -- ------------------------------ */
		CHECK_CLOSE(2.0 / (2.0 * n + 1.0) - q, rem, 0.0, 1.5e-14, "n = %d: remainder for x^(2n)", n);
		if (table) printf("GaussLeg n = %2d, x^(2n): remainder %.6e, error in it %.1e\n",
						  n, rem, fabs(2.0 / (2.0 * n + 1.0) - q - rem));
	}
}

/* --- Closed forms for n = 2, 3, 4. The weights are off by 1e-14 or
	   less (relative; see test_gaussleg_rules), measured 5e-15 at n = 3
	   and 1.8e-14 at n = 4 -- ----------------------------------- */

static void test_gaussleg_closed(void)
{
	double x[4], w[4];

	GaussLeg(-1.0, 1.0, x, w, 2);
	CHECK_CLOSE(x[0], -1.0 / sqrt(3.0), 1e-15, 0.0, "n = 2, node %d", 0);
	CHECK_CLOSE(x[1], 1.0 / sqrt(3.0), 1e-15, 0.0, "n = 2, node %d", 1);
	CHECK_CLOSE(w[0], 1.0, 1e-15, 0.0, "n = 2, weight %d", 0);
	CHECK_CLOSE(w[1], 1.0, 1e-15, 0.0, "n = 2, weight %d", 1);

	GaussLeg(-1.0, 1.0, x, w, 3);
	CHECK_CLOSE(x[0], -sqrt(0.6), 1e-15, 0.0, "n = 3, node %d", 0);
	CHECK(fabs(x[1]) < 1e-15, "n = 3, middle node %.3g", x[1]);
	CHECK_CLOSE(x[2], sqrt(0.6), 1e-15, 0.0, "n = 3, node %d", 2);
	CHECK_CLOSE(w[0], 5.0 / 9.0, 4e-14, 0.0, "n = 3, weight %d", 0);
	CHECK_CLOSE(w[1], 8.0 / 9.0, 4e-14, 0.0, "n = 3, weight %d", 1);
	CHECK_CLOSE(w[2], 5.0 / 9.0, 4e-14, 0.0, "n = 3, weight %d", 2);

	GaussLeg(-1.0, 1.0, x, w, 4);
	{
		double xo = sqrt(3.0 / 7.0 + 2.0 / 7.0 * sqrt(1.2)), xi = sqrt(3.0 / 7.0 - 2.0 / 7.0 * sqrt(1.2));
		double wo = (18.0 - sqrt(30.0)) / 36.0, wi = (18.0 + sqrt(30.0)) / 36.0;

		CHECK_CLOSE(x[0], -xo, 1e-15, 0.0, "n = 4, node %d", 0);
		CHECK_CLOSE(x[1], -xi, 1e-15, 0.0, "n = 4, node %d", 1);
		CHECK_CLOSE(x[2], xi, 1e-15, 0.0, "n = 4, node %d", 2);
		CHECK_CLOSE(x[3], xo, 1e-15, 0.0, "n = 4, node %d", 3);
		CHECK_CLOSE(w[0], wo, 4e-14, 0.0, "n = 4, weight %d", 0);
		CHECK_CLOSE(w[1], wi, 4e-14, 0.0, "n = 4, weight %d", 1);
		CHECK_CLOSE(w[2], wi, 4e-14, 0.0, "n = 4, weight %d", 2);
		CHECK_CLOSE(w[3], wo, 4e-14, 0.0, "n = 4, weight %d", 3);
	}
}

/* --- Reversed interval and a smooth function -- ------------------- */

static void test_gaussleg_other(void)
{
	double x[16], w[16], q = 0.0, qr = 0.0;

	/* --- int_0^pi sin x dx = 2, with 8 points; reversed: -2 -- --- */

	GaussLeg(0.0, acos(-1.0), x, w, 8);
	for (int i = 0; i < 8; i++) q += w[i] * sin(x[i]);
	CHECK_CLOSE(q, 2.0, 6e-15, 0.0, "int sin, n = %d", 8);
	if (table) printf("GaussLeg int_0^pi sin x, n = 8: error %.2e\n", fabs(q - 2.0));

	GaussLeg(acos(-1.0), 0.0, x, w, 8);
	for (int i = 0; i < 8; i++) qr += w[i] * sin(x[i]);
	CHECK_CLOSE(qr, -2.0, 6e-15, 0.0, "reversed int sin, n = %d", 8);
}

/* --- w2, w3, U3 -- ------------------------------------------------ */

/* --- Maximum relative errors by function and regime; see below -- */

enum { SERIES,
	   EXACT,
	   LIMIT,
	   NREGIME };

static const char *const regime_name[NREGIME] = {"series", "exact expressions", "large dtau"};

/* --- Bounds: about twice the maximum measured (gcc, x86-64), by
	   regime and element. Series: truncated after the cubic term (w) and
	   the quartic one (U), so the relative error grows to ~dtau^2/4 at
	   the switch (dtau = 5e-4 for w, 2e-4 for U): 6e-8 in w[1], 7e-8 in
	   w[2]. Exact expressions: cancellation at small dtau, which is
	   worst just above the switch, where it is much larger than the
	   series' error: w[2] is off by 9.8e-7 at dtau = 5e-4 and U[2] by
	   8.5e-6 at dtau = 2e-4. Large dtau: exact -- ----------------- */

static const double w_bound[NREGIME][3] = {
	/* series */ {1e-7, 1.2e-7, 1.5e-7},
	/* exact */ {1e-13, 3.5e-10, 2e-6},
	/* limit */ {1e-15, 1e-15, 1e-15}};

static const double u_bound[NREGIME][3] = {
	/* series */ {1.3e-8, 6.5e-9, 4e-9},
	/* exact */ {2e-13, 1.2e-9, 1.8e-5},
	/* limit */ {1e-15, 1e-15, 1e-15}};

static void test_w3_u3(void)
{
	double w_err[NREGIME][3] = {{0}}, u_err[NREGIME][3] = {{0}};

	for (int n = 0; n < DTAU_REF_N; n++) {
		double dt = dtau_ref[n], w[3], w_two[3] = {0.0, 0.0, 0.0}, u[3];
		double wref[3] = {w_ref_0[n], w_ref_1[n], w_ref_2[n]};
		double uref[3] = {u_ref_0[n], u_ref_1[n], u_ref_2[n]};
		int    rw      = (dt < 5e-4) ? SERIES : (dt > 50.0) ? LIMIT
															: EXACT;
		int    ru      = (dt < 2e-4) ? SERIES : (dt > 200.0) ? LIMIT
															 : EXACT;

		w3(dt, w);
		w2(dt, w_two);
		U3(dt, u);

		for (int k = 0; k < 3; k++) {
			double ew = (wref[k] > 0.0) ? fabs(w[k] - wref[k]) / wref[k] : fabs(w[k]);
			double eu = (uref[k] > 0.0) ? fabs(u[k] - uref[k]) / uref[k] : fabs(u[k]);

			w_err[rw][k] = MAX(w_err[rw][k], ew);
			u_err[ru][k] = MAX(u_err[ru][k], eu);
			CHECK(ew < w_bound[rw][k], "w3[%d](%g) = %.17g, want %.17g (rel. error %.3g)", k, dt, w[k], wref[k], ew);
			CHECK(eu < u_bound[ru][k], "U3[%d](%g) = %.17g, want %.17g (rel. error %.3g)", k, dt, u[k], uref[k], eu);
		}
		/* --- w2 is the first two elements of w3 -- ------------- */

		CHECK(w_two[0] == w[0] && w_two[1] == w[1], "w2(%g) = %.17g, %.17g differs from w3: %.17g, %.17g",
			  dt, w_two[0], w_two[1], w[0], w[1]);
	}
	if (table)
		for (int r = 0; r < NREGIME; r++)
			printf("%-18s w: %.2e %.2e %.2e   U: %.2e %.2e %.2e\n", regime_name[r],
				   w_err[r][0], w_err[r][1], w_err[r][2], u_err[r][0], u_err[r][1], u_err[r][2]);
}

/* --- Small dtau: w[k] = dtau^(k+1) / (k+1) (1 + O(dtau)) and
	   U[k] = dtau^(k+1) / (k+1) (1 + O(dtau)), also far below the
	   range of the reference table -- ----------------------------- */

static void test_w3_u3_tiny(void)
{
	for (double dt = 1e-300; dt < 1e-6; dt *= 1e7) {
		double w[3], u[3];

		w3(dt, w);
		U3(dt, u);
		for (int k = 0; k < 3; k++) {
			double lead = pow(dt, k + 1.0) / (k + 1.0);

			if (lead == 0.0) continue; /* underflow */
			CHECK_CLOSE(w[k], lead, 2e-6, 0.0, "w3[%d](%g)", k, dt);
			CHECK_CLOSE(u[k], lead, 2e-6, 0.0, "U3[%d](%g)", k, dt);
		}
	}
}

int main(int argc, char *argv[])
{
	table = (argc > 1 && strcmp(argv[1], "--table") == 0);

	test_gaussleg_rules();
	test_gaussleg_remainder();
	test_gaussleg_closed();
	test_gaussleg_other();
	test_w3_u3();
	test_w3_u3_tiny();

	return check_summary("test_quadrature");
}
