/* ------- file: -------------------------- test_functions.c -------

	   Special functions: expint.c (E1, E2), planck.c (Planck) and
	   wigner.c (w3js, w6js).

	   -- E1, E2 against mpmath (kernels_reference.h, from
		  make_kernels_reference.py) for 1e-6 <= x <= 80; the routines
		  are the approximations of Abramowitz & Stegun 5.1.53 and
		  5.1.56 and their accuracy is that of the approximations
		  (abs. 2e-7 for E1 at x <= 1, rel. 2e-8 for x e^x E1 above),
		  E2 loses relative accuracy at large x because of the
		  recurrence E2 = exp(-x) - x E1; both return 0 for x > 80;
	   -- Planck against the formula in high precision with the constants
		  of constant.h, the Wien and Rayleigh-Jeans limits, the
		  integral over all frequencies (Stefan-Boltzmann law, from
		  the constants and Gauss-Legendre quadrature of Planck in the
		  variable h nu / kT), the cut-off for h nu / kT > 150;
	   -- w3js against closed forms: (j j 0; m -m 0), (j1 j2 j3; 0 0 0),
		  the stretched case j3 = j1 + j2, tabulated values, the
		  orthogonality of the 3j symbols, the symmetries under
		  permutation of the columns and under m -> -m, the selection
		  rules; w6js against {a b c; 0 c b} and the sum rules.

	   Run with --table to print the measured errors.
	   --                                              -------------- */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "spectrum.h"
#include "constant.h"
#include "check.h"
#include "kernels_reference.h"

double E1(double x);
double E2(double x);

static bool_t table;

/* --- E1 and E2 -- ------------------------------------------------- */

/* --- Measured (gcc, x86-64): E1 absolute error at most 2.0e-7 for
	   x <= 1 (the stated accuracy of 5.1.53: 2e-7), relative error of
	   E1 at most 1.5e-7 over the grid and 2.0e-8 for x > 1 (5.1.56).
	   E2 relative error 1.3e-7 for x < 1 (from E1, through the
	   recurrence), and from the continued fraction 5.8e-15 at x = 1
	   and 6.4e-16 for x >= 10. Before N-22 E2 used the recurrence
	   everywhere and lost accuracy to the cancellation in exp(-x) - x E1
	   (1.2e-6 at x = 80), and both returned 0 beyond x = 80.
	   Bounds: about twice -- ----------------------------------- */

static void test_expint(void)
{
	double max_abs_small = 0.0, max_rel_large = 0.0, max_rel_e1 = 0.0,
		   max_rel_e2_small = 0.0, max_rel_e2_mid = 0.0, max_rel_e2_large = 0.0;

	for (int n = 0; n < EXP_REF_N; n++) {
		double x = exp_ref_x[n], e1 = E1(x), e2 = E2(x);
		double r1 = fabs(e1 - exp_ref_E1[n]) / exp_ref_E1[n];
		double r2 = fabs(e2 - exp_ref_E2[n]) / exp_ref_E2[n];

		max_rel_e1 = MAX(max_rel_e1, r1);
		if (x <= 1.0)
			max_abs_small = MAX(max_abs_small, fabs(e1 - exp_ref_E1[n]));
		else
			max_rel_large = MAX(max_rel_large, r1);

		if (x < 1.0)
			max_rel_e2_small = MAX(max_rel_e2_small, r2);
		else if (x < 10.0)
			max_rel_e2_mid = MAX(max_rel_e2_mid, r2);
		else
			max_rel_e2_large = MAX(max_rel_e2_large, r2);
	}
	CHECK(max_abs_small < 4.1e-7, "E1, x <= 1: max absolute error %.3g", max_abs_small);
	CHECK(max_rel_e1 < 3e-7, "E1: max relative error %.3g", max_rel_e1);
	CHECK(max_rel_large < 4e-8, "E1, x > 1: max relative error %.3g", max_rel_large);
	CHECK(max_rel_e2_small < 2.6e-7, "E2, x < 1: max relative error %.3g", max_rel_e2_small);
	CHECK(max_rel_e2_mid < 1.2e-14, "E2, 1 <= x < 10: max relative error %.3g", max_rel_e2_mid);
	CHECK(max_rel_e2_large < 1.3e-15, "E2, x >= 10: max relative error %.3g", max_rel_e2_large);
	if (table)
		printf("E1: abs. %.2e (x <= 1), rel. %.2e (x > 1), rel. %.2e (all)\n"
			   "E2: rel. %.2e (x < 1), %.2e (1 .. 10), %.2e (>= 10)\n",
			   max_abs_small, max_rel_large, max_rel_e1, max_rel_e2_small,
			   max_rel_e2_mid, max_rel_e2_large);

	/* --- Zero only where exp(-x)/x underflows (x ~ 740; subnormal from
		   x ~ 708), and for x = inf (fixedrate.c with T = 0); before N-22
		   both returned 0 beyond x = 80 -- ----------------------- */

	CHECK(E1(700.0) > 0.0 && E2(700.0) > 0.0, "E1(700) = %g, E2(700) = %g", E1(700.0), E2(700.0));
	CHECK(E1(800.0) == 0.0 && E2(800.0) == 0.0, "underflow: E1(800) = %g, E2(800) = %g", E1(800.0), E2(800.0));
	CHECK(E1(1e200) == 0.0 && E2(1e200) == 0.0, "E1(1e200) = %g, E2(1e200) = %g", E1(1e200), E2(1e200));
	CHECK(E1(HUGE_VAL) == 0.0 && E2(HUGE_VAL) == 0.0, "E1(inf) = %g, E2(inf) = %g", E1(HUGE_VAL), E2(HUGE_VAL));

	/* --- Properties: E1 and E2 positive and decreasing, E2 < E1
		   (for x > 0), E2(x) -> 1 as x -> 0 (E2(0) = 1) and
		   E1 ~ -ln x - gamma -- ------------------------------ */

	{
		double prev1 = HUGE_VAL, prev2 = HUGE_VAL;

		for (int k = 1; k <= 3000; k++) {
			double x = 0.2 * k, e1 = E1(x), e2 = E2(x);

			CHECK(e1 > 0.0 && e2 > 0.0 && e1 < prev1 && e2 < prev2 && e2 < e1,
				  "monotonic and positive at x = %g: E1 %g, E2 %g", x, e1, e2);
			prev1 = e1;
			prev2 = e2;
		}
	}
	CHECK_CLOSE(E2(1e-9), 1.0, 5e-8, 0.0, "E2 near 0, code %d", 0);
	CHECK_CLOSE(E1(1e-9), -log(1e-9) - 0.5772156649015329, 5e-8, 0.0, "E1 near 0, code %d", 0);
}

/* --- Planck -- ---------------------------------------------------- */

/* --- The argument of the exponential at lambda [nm] and T -- ----- */

static double x_of(double lambda, double T)
{
	return (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M * lambda * T);
}

/* --- The relative error is not the rounding error of one operation:
	   B ~ 1 / (exp(x) - 1), where x = h nu / kT carries a relative
	   error of a few eps, which becomes x eps in B (large x). Before
	   N-22 exp(x) - 1 was evaluated without expm1, which lost eps / x
	   at small x; the bound c eps (1 + x + 1/x) still allows for that,
	   over the reference grid (x from 1.4e-8 to 149); the largest
	   error / (eps (1 + x + 1/x)) measured is 0.77 (gcc, x86-64),
	   set by large x; c = 1.6 -- */

static void test_planck_reference(void)
{
	double worst = 0.0;

	for (int n = 0; n < PLANCK_REF_N; n++) {
		double B, x = x_of(planck_ref_lambda[n], planck_ref_T[n]);
		double scale = DBL_EPSILON * (1.0 + x + 1.0 / x), err;

		Planck(1, (double *) &planck_ref_T[n], planck_ref_lambda[n], &B);
		err   = fabs(B - planck_ref_B[n]) / planck_ref_B[n];
		worst = MAX(worst, err / scale);
		CHECK(err < 1.6 * scale, "lambda = %g nm, T = %g K, x = %g: relative error %.3g, bound %.3g",
			  planck_ref_lambda[n], planck_ref_T[n], x, err, 1.6 * scale);
	}
	if (table)
		printf("Planck: max relative error / (eps (1 + x + 1/x)) = %.2f over %d points\n",
			   worst, PLANCK_REF_N);
}

static void test_planck_limits(void)
{
	double lambda = 500.0, T, B;

	/* --- Wien: B = (2 h nu^3 / c^2) exp(-x) (1 + exp(-x) + ...); the
		   ratio to the Wien form is 1 / (1 - exp(-x)) -- ------- */

	for (double x = 5.0; x <= 140.0; x += 5.0) {
		double nu_c = 1.0 / (NM_TO_M * lambda);
		double wien = 2.0 * HPLANCK * CLIGHT * nu_c * nu_c * nu_c * exp(-x);

		T = x_of(lambda, 1.0) / x;
		Planck(1, &T, lambda, &B);
		/* --- the relative error is ~ x eps (see test_planck_reference): measured
		   at most 0.6 eps (1 + x) -- ---------------------------------- */

		CHECK_CLOSE(B / wien, 1.0 / (1.0 - exp(-x)), 1.2 * DBL_EPSILON * (1.0 + x), 0.0, "Wien, x = %g", x);
	}

	/* --- Rayleigh-Jeans: B = 2 nu^2 k T / c^2 (1 - x/2 + x^2/12 -
		   x^4/720 + ...). With expm1 the error is rounding plus the
		   truncated x^4/720 term: measured 1.6e-15 at x = 1e-3 (the
		   x^4 term) and 1.1e-16 at 1e-7; before N-22 (exp(x) - 1)
		   4.5e-14 and 5.7e-10. Bound 4 eps + x^4/360 -- --------- */

	for (int k = 0; k < 2; k++) {
		double x    = (k == 0) ? 1e-3 : 1e-7;
		double nu_c = 1.0 / (NM_TO_M * lambda), series = 1.0 - x / 2.0 + x * x / 12.0;

		T = x_of(lambda, 1.0) / x;
		Planck(1, &T, lambda, &B);
		if (table) printf("Planck, Rayleigh-Jeans x = %g: relative error %.2e\n", x, fabs(B / (2.0 * KBOLTZMANN * T * nu_c * nu_c * series) - 1.0));
		CHECK_CLOSE(B, 2.0 * KBOLTZMANN * T * nu_c * nu_c * series, 4.0 * DBL_EPSILON + x * x * x * x / 360.0, 0.0,
					"Rayleigh-Jeans, x = %g", x);
	}

	/* --- The cut-off: Planck is exactly zero for x > 150 (the true value
		   is ~1e-77 J m^-2 s^-1 Hz^-1 sr^-1 at 500 nm) -- --------- */

	T = x_of(lambda, 1.0) / 149.9;
	Planck(1, &T, lambda, &B);
	CHECK(B > 0.0 && B < 1e-70, "x = 149.9: B = %g", B);
	T = x_of(lambda, 1.0) / 150.1;
	Planck(1, &T, lambda, &B);
	CHECK(B == 0.0, "x = 150.1: B = %g", B);

	/* --- Arrays: all depths done, in order -- ----------------- */

	{
		double Ta[5] = {3000.0, 4500.0, 6000.0, 8000.0, 1e4}, Ba[5];

		Planck(5, Ta, 656.3, Ba);
		for (int k = 0; k < 5; k++) {
			double single;

			Planck(1, &Ta[k], 656.3, &single);
			CHECK(Ba[k] == single, "array element %d: %.17g vs %.17g", k, Ba[k], single);
		}
		for (int k = 1; k < 5; k++) CHECK(Ba[k] > Ba[k - 1], "B increases with T: %g %g", Ba[k - 1], Ba[k]);
	}
}

/* --- Stefan-Boltzmann: the integral of B_nu over frequency is
	   2 pi^4 k^4 T^4 / (15 h^3 c^2). Integrate over u = ln x,
	   x = h nu / kT from 1e-3 to 100 with 24-point Gauss-Legendre
	   in 60 panels, and add the Rayleigh-Jeans tail below
	   (int_0^x0 x^2 (1 - x/2) dx). Measured relative error 1.1e-14;
	   bound 3e-14 --                                               - */

static void test_planck_stefan_boltzmann(void)
{
	enum { NP = 60,
		   NG = 24 };
	double T = 5777.0, xg[NG], wg[NG], sum = 0.0, x0 = 1e-3, x1 = 100.0;
	double u0 = log(x0), u1 = log(x1), pi = acos(-1.0), exact, nu_unit, lambda_unit;

	GaussLeg(-1.0, 1.0, xg, wg, NG);

	/* --- B_nu(x) in J m^-2 s^-1 Hz^-1 sr^-1; frequency nu = x kT / h --- */

	nu_unit     = KBOLTZMANN * T / HPLANCK;
	lambda_unit = CLIGHT / nu_unit / NM_TO_M; /* lambda [nm] at x = 1 */
	for (int p = 0; p < NP; p++) {
		double ua = u0 + (u1 - u0) * p / NP, ub = u0 + (u1 - u0) * (p + 1) / NP;

		for (int g = 0; g < NG; g++) {
			double u = 0.5 * (ua + ub) + 0.5 * (ub - ua) * xg[g], x = exp(u), B;

			Planck(1, &T, lambda_unit / x, &B);
			sum += 0.5 * (ub - ua) * wg[g] * B * nu_unit * x; /* dnu = nu_unit dx = nu_unit x du */
		}
	}
	sum += (2.0 * HPLANCK / (CLIGHT * CLIGHT)) * pow(nu_unit, 3.0) * nu_unit * (x0 * x0 * x0 / 3.0 - x0 * x0 * x0 * x0 / 8.0);

	exact = 2.0 * pow(pi, 4.0) * pow(KBOLTZMANN * T, 4.0) / (15.0 * pow(HPLANCK, 3.0) * CLIGHT * CLIGHT);
	CHECK_CLOSE(sum, exact, 3e-14, 0.0, "integral of B_nu over frequency, T = %g", T);
	if (table) printf("Planck Stefan-Boltzmann: relative error %.2e\n", fabs(sum / exact - 1.0));
}

/* --- Wigner 3-j and 6-j symbols -- ------------------------------- */

static double sign(double k)
{
	return (fmod(fabs(k), 2.0) < 0.5) ? 1.0 : -1.0; /* (-1)^k, k integer */
}

/* --- (j j 0; m -m 0) = (-1)^(j-m) / sqrt(2j+1) -- ---------------- */

static void test_w3js_j_j_0(void)
{
	for (double j = 0.0; j <= 10.0; j += 0.5)
		for (double m = -j; m <= j + 1e-9; m += 1.0) {
			CHECK_CLOSE(w3js(j, j, 0.0, m, -m, 0.0), sign(j - m) / sqrt(2.0 * j + 1.0), 1e-14, 0.0,
						"(j j 0; m -m 0), j = %g, m = %g", j, m);
			CHECK_CLOSE(w3js(j, 0.0, j, m, 0.0, -m), sign(2.0 * j) * sign(j - m) / sqrt(2.0 * j + 1.0), 1e-14, 0.0,
						"(j 0 j; m 0 -m) (an odd permutation), j = %g, m = %g", j, m);
		}
}

/* --- (j1 j2 j3; 0 0 0) for integer j: zero if J = j1+j2+j3 is odd,
	   else (-1)^g sqrt[(2g-2j1)!(2g-2j2)!(2g-2j3)!/(2g+1)!]
	   g! / ((g-j1)!(g-j2)!(g-j3)!), 2g = J. Evaluated with lgamma
	   (accurate to ~1e-14 relative) -- --------------------------- */

static double w3js_zero_closed(int j1, int j2, int j3)
{
	int J = j1 + j2 + j3, g = J / 2;

	if (J % 2) return 0.0;
	return sign(g) * exp(0.5 * (lgamma(2.0 * g - 2.0 * j1 + 1.0) + lgamma(2.0 * g - 2.0 * j2 + 1.0) +
								lgamma(2.0 * g - 2.0 * j3 + 1.0) - lgamma(2.0 * g + 2.0)) +
						 lgamma(g + 1.0) - lgamma(g - j1 + 1.0) - lgamma(g - j2 + 1.0) - lgamma(g - j3 + 1.0));
}

/* --- Measured (gcc, x86-64): relative error at most 4e-14 for
	   j <= 24 (the alternating sum and the factorials lose digits);
	   see also the large-j test -- -------------------------------- */

static void test_w3js_zero(void)
{
	double worst = 0.0, worst_zero = 0.0;

	for (int j1 = 0; j1 <= 24; j1++)
		for (int j2 = 0; j2 <= 24; j2++)
			for (int j3 = abs(j1 - j2); j3 <= j1 + j2; j3++) {
				double got = w3js(j1, j2, j3, 0, 0, 0), want = w3js_zero_closed(j1, j2, j3);

				if (want != 0.0) {
					worst = MAX(worst, fabs(got - want) / fabs(want));
					CHECK_CLOSE(got, want, 3.5e-12, 0.0, "(%d %d %d; 0 0 0)", j1, j2, j3);
				} else {
					/* --- odd j1 + j2 + j3: zero, but the alternating
						   sum leaves up to 5.7e-14 at j <= 24 -- ------------- */

					worst_zero = MAX(worst_zero, fabs(got));
					CHECK(fabs(got) < 1.2e-13, "(%d %d %d; 0 0 0) = %g, want 0", j1, j2, j3, got);
				}
			}
	if (table)
		printf("w3js (j1 j2 j3; 0 0 0), j <= 24: max relative error %.2e, max |value| when it must be 0: %.2e\n",
			   worst, worst_zero);
}

/* --- Stretched case j3 = j1 + j2 (Edmonds 3.7.? ), M = m1 + m2:
	   (-1)^(j1-j2+M) sqrt[(2j1)!(2j2)!(j1+j2+M)!(j1+j2-M)! /
	   ((2j1+2j2+1)!(j1+m1)!(j1-m1)!(j2+m2)!(j2-m2)!)] -- ---------- */

static void test_w3js_stretched(void)
{
	double worst = 0.0;

	for (double j1 = 0.0; j1 <= 6.0; j1 += 0.5)
		for (double j2 = 0.0; j2 <= 6.0; j2 += 0.5)
			for (double m1 = -j1; m1 <= j1 + 1e-9; m1 += 1.0)
				for (double m2 = -j2; m2 <= j2 + 1e-9; m2 += 1.0) {
					double M = m1 + m2, j3 = j1 + j2;
					double want = sign(j1 - j2 + M) *
								  exp(0.5 * (lgamma(2 * j1 + 1) + lgamma(2 * j2 + 1) + lgamma(j3 + M + 1) + lgamma(j3 - M + 1) -
											 lgamma(2 * j3 + 2) - lgamma(j1 + m1 + 1) - lgamma(j1 - m1 + 1) - lgamma(j2 + m2 + 1) - lgamma(j2 - m2 + 1)));
					double got = w3js(j1, j2, j3, m1, m2, -M);

					worst = MAX(worst, fabs(got - want) / fabs(want));
					CHECK_CLOSE(got, want, 3e-14, 0.0, "stretched (%g %g %g; %g %g %g)", j1, j2, j3, m1, m2, -M);
				}
	if (table) printf("w3js stretched case: max relative error %.2e\n", worst);
}

/* --- A few tabulated values -- ------------------------------------ */

static void test_w3js_values(void)
{
	/* --- the (1/2 1/2 1) symbols follow from the Clebsch-Gordan
		   coefficients of two spins 1/2; (1 1 2) and (2 2 2) are
		   from the closed form for m = 0 and the standard tables -- */

	CHECK_CLOSE(w3js(0.5, 0.5, 1.0, 0.5, -0.5, 0.0), 1.0 / sqrt(6.0), 1e-14, 0.0, "(1/2 1/2 1; 1/2 -1/2 0), code %d", 0);
	CHECK_CLOSE(w3js(0.5, 0.5, 1.0, -0.5, 0.5, 0.0), 1.0 / sqrt(6.0), 1e-14, 0.0, "(1/2 1/2 1; -1/2 1/2 0), code %d", 1);
	CHECK_CLOSE(w3js(0.5, 0.5, 1.0, 0.5, 0.5, -1.0), -1.0 / sqrt(3.0), 1e-14, 0.0, "(1/2 1/2 1; 1/2 1/2 -1), code %d", 2);
	CHECK_CLOSE(w3js(0.5, 0.5, 1.0, -0.5, -0.5, 1.0), -1.0 / sqrt(3.0), 1e-14, 0.0, "(1/2 1/2 1; -1/2 -1/2 1), code %d", 3);
	CHECK_CLOSE(w3js(1.0, 1.0, 2.0, 0.0, 0.0, 0.0), sqrt(2.0 / 15.0), 1e-14, 0.0, "(1 1 2; 0 0 0), code %d", 0);
	CHECK_CLOSE(w3js(1.0, 1.0, 2.0, 1.0, -1.0, 0.0), sqrt(1.0 / 30.0), 1e-14, 0.0, "(1 1 2; 1 -1 0), code %d", 1);
	CHECK_CLOSE(w3js(2.0, 2.0, 2.0, 0.0, 0.0, 0.0), -sqrt(2.0 / 35.0), 1e-14, 0.0, "(2 2 2; 0 0 0), code %d", 2);
	CHECK_CLOSE(w3js(1.0, 1.0, 1.0, 1.0, -1.0, 0.0), 1.0 / sqrt(6.0), 1e-14, 0.0, "(1 1 1; 1 -1 0), code %d", 3);
}

/* --- Bounds in the tests of the 3j and 6j symbols (except where
	   stated): 5e-15, the measured errors being below 4e-15 -- */

/* --- Orthogonality: sum over m1, m2 of (2 j3 + 1) (j1 j2 j3; m1 m2 m3)
	   (j1 j2 j3'; m1 m2 m3') = delta(j3, j3') delta(m3, m3') -- --- */

static void test_w3js_orthogonality(void)
{
	static const double triples[][2] = {{0.5, 0.5}, {1.0, 0.5}, {1.5, 1.0}, {2.0, 2.0}, {2.5, 1.5}, {3.0, 4.0}, {6.5, 5.0}};
	double              worst        = 0.0;

	for (size_t t = 0; t < sizeof(triples) / sizeof(triples[0]); t++) {
		double j1 = triples[t][0], j2 = triples[t][1];

		for (double j3 = fabs(j1 - j2); j3 <= j1 + j2 + 1e-9; j3 += 1.0)
			for (double j3p = fabs(j1 - j2); j3p <= j1 + j2 + 1e-9; j3p += 1.0)
				for (double m3 = -j3; m3 <= j3 + 1e-9; m3 += 1.0) {
					double m3p = m3, sum = 0.0;

					if (fabs(m3p) > j3p) continue;
					for (double m1 = -j1; m1 <= j1 + 1e-9; m1 += 1.0)
						for (double m2 = -j2; m2 <= j2 + 1e-9; m2 += 1.0)
							sum += (2.0 * j3 + 1.0) * w3js(j1, j2, j3, m1, m2, -m3) * w3js(j1, j2, j3p, m1, m2, -m3p);
					worst = MAX(worst, fabs(sum - ((j3 == j3p) ? 1.0 : 0.0)));
					CHECK_CLOSE(sum, (j3 == j3p) ? 1.0 : 0.0, 5e-15, 5e-15,
								"orthogonality j1 = %g, j2 = %g, j3 = %g, j3' = %g, m3 = %g", j1, j2, j3, j3p, m3);
				}
	}
	if (table) printf("w3js orthogonality: max deviation %.2e\n", worst);
}

/* --- Symmetries: odd permutations of columns and m -> -m give
	   (-1)^(j1+j2+j3), even permutations leave the value -- ------- */

static void test_w3js_symmetry(void)
{
	static const double j[][3] = {{1.0, 1.0, 1.0}, {2.0, 1.5, 1.5}, {3.0, 2.0, 2.5 + 2.5}, {2.5, 4.0, 3.5}, {5.0, 6.0, 7.0}};

	for (size_t t = 0; t < sizeof(j) / sizeof(j[0]); t++) {
		double j1 = j[t][0], j2 = j[t][1], j3 = j[t][2], ph = sign(j1 + j2 + j3);

		if (j3 > j1 + j2 || j3 < fabs(j1 - j2) || fmod(j1 + j2 + j3, 1.0) != 0.0) continue;
		for (double m1 = -j1; m1 <= j1 + 1e-9; m1 += 1.0)
			for (double m2 = -j2; m2 <= j2 + 1e-9; m2 += 1.0) {
				double m3 = -m1 - m2, v = w3js(j1, j2, j3, m1, m2, m3);

				if (fabs(m3) > j3) continue;
				CHECK_CLOSE(w3js(j2, j3, j1, m2, m3, m1), v, 5e-15, 1e-15, "cyclic, j = %g %g %g, m = %g %g", j1, j2, j3, m1, m2);
				CHECK_CLOSE(w3js(j3, j1, j2, m3, m1, m2), v, 5e-15, 1e-15, "cyclic 2, j = %g %g %g, m = %g %g", j1, j2, j3, m1, m2);
				CHECK_CLOSE(w3js(j2, j1, j3, m2, m1, m3), ph * v, 5e-15, 1e-15, "swap 12, j = %g %g %g, m = %g %g", j1, j2, j3, m1, m2);
				CHECK_CLOSE(w3js(j1, j3, j2, m1, m3, m2), ph * v, 5e-15, 1e-15, "swap 23, j = %g %g %g, m = %g %g", j1, j2, j3, m1, m2);
				CHECK_CLOSE(w3js(j3, j2, j1, m3, m2, m1), ph * v, 5e-15, 1e-15, "swap 13, j = %g %g %g, m = %g %g", j1, j2, j3, m1, m2);
				CHECK_CLOSE(w3js(j1, j2, j3, -m1, -m2, -m3), ph * v, 5e-15, 1e-15, "m -> -m, j = %g %g %g, m = %g %g", j1, j2, j3, m1, m2);
			}
	}
}

static void test_w3js_selection(void)
{
	/* --- All of these are zero by the selection rules -- ------- */

	CHECK(w3js(1.0, 1.0, 3.0, 0.0, 0.0, 0.0) == 0.0, "triangle: j3 > j1 + j2, value %g", w3js(1.0, 1.0, 3.0, 0.0, 0.0, 0.0));
	CHECK(w3js(3.0, 1.0, 1.0, 0.0, 0.0, 0.0) == 0.0, "triangle: j1 > j2 + j3, value %g", w3js(3.0, 1.0, 1.0, 0.0, 0.0, 0.0));
	CHECK(w3js(1.0, 1.0, 1.0, 1.0, 1.0, 0.0) == 0.0, "m1 + m2 + m3 != 0, value %g", w3js(1.0, 1.0, 1.0, 1.0, 1.0, 0.0));
	CHECK(w3js(1.0, 1.0, 1.0, 2.0, -1.0, -1.0) == 0.0, "|m1| > j1, value %g", w3js(1.0, 1.0, 1.0, 2.0, -1.0, -1.0));
	CHECK(w3js(1.0, 1.0, 1.0, 0.0, 0.0, 0.0) == 0.0, "(1 1 1; 0 0 0) odd J, value %g", w3js(1.0, 1.0, 1.0, 0.0, 0.0, 0.0));
	CHECK(w3js(0.3, 1.0, 1.0, 0.0, 0.0, 0.0) == 0.0, "j1 not a half-integer, value %g", w3js(0.3, 1.0, 1.0, 0.0, 0.0, 0.0));
	CHECK(w3js(1.0, 1.0, 1.0, 0.25, -0.25, 0.0) == 0.0, "m1 not a half-integer, value %g", w3js(1.0, 1.0, 1.0, 0.25, -0.25, 0.0));
}

/* --- Large j: the symbol is a sum of alternating terms with products
	   of factorials. Closed form for (j j j; 0 0 0), j even. The
	   relative error grows with j, and depends on the rounding:
	   2e-13 at j = 20, 1e-9 at j = 40, 1.2e-8 at j = 44, 1.5e-7 at
	   j = 48. Above that it breaks down: the product of six
	   factorials in the square root overflows, giving +-inf at
	   j = 50 .. 56 (j! > 1e51 from j = 42, (j!)^6 > 1e308 from
	   j = 49), and 0 from j = 58, where the (3j+1)! of the triangle
	   coefficient overflows (170! is the largest finite factorial).
	   Bounds: 5e-13 for j <= 20 and 1e-6 for j = 32 .. 48 -- ---- */

static void test_w3js_large(void)
{
	double worst20 = 0.0, worst48 = 0.0;

	for (int j = 2; j <= 48; j += 2) {
		double got = w3js(j, j, j, 0, 0, 0), want = w3js_zero_closed(j, j, j);
		double err = fabs(got - want) / fabs(want);

		if (j <= 20) {
			worst20 = MAX(worst20, err);
			CHECK(err < 5e-13, "(j j j; 0 0 0), j = %d: got %.15g, want %.15g (rel. error %.3g)", j, got, want, err);
		}
		if (j >= 32) {
			worst48 = MAX(worst48, err);
			CHECK(err < 1e-6, "(j j j; 0 0 0), j = %d: got %.15g, want %.15g (rel. error %.3g)", j, got, want, err);
		}
	}
	if (table)
		printf("w3js (j j j; 0 0 0): max relative error %.2e for j <= 20, %.2e for 32 <= j <= 48\n", worst20, worst48);
}

/* --- 6-j symbols -- ----------------------------------------------- */

static void test_w6js(void)
{
	double worst = 0.0;

	/* --- {1 1 1; 1 1 1} = 1/6, {1/2 1/2 1; 1/2 1/2 1} = 1/6,
		   {a b c; 0 c b} = (-1)^(a+b+c) / sqrt((2b+1)(2c+1)) -- --- */

	CHECK_CLOSE(w6js(1.0, 1.0, 1.0, 1.0, 1.0, 1.0), 1.0 / 6.0, 1e-14, 0.0, "{1 1 1; 1 1 1}, code %d", 0);
	CHECK_CLOSE(w6js(0.5, 0.5, 1.0, 0.5, 0.5, 1.0), 1.0 / 6.0, 1e-14, 0.0, "{1/2 1/2 1; 1/2 1/2 1}, code %d", 1);
	for (double a = 0.0; a <= 5.0; a += 0.5)
		for (double b = 0.0; b <= 5.0; b += 0.5)
			for (double c = fabs(a - b); c <= a + b + 1e-9; c += 1.0) {
				double want = sign(a + b + c) / sqrt((2 * b + 1) * (2 * c + 1)), got = w6js(a, b, c, 0.0, c, b);

				worst = MAX(worst, fabs(got - want) / fabs(want));
				CHECK_CLOSE(got, want, 5e-15, 0.0, "{%g %g %g; 0 %g %g}", a, b, c, c, b);
			}

	/* --- Sum rules: sum_x (2x+1) {a b x; a b f} = (-1)^(2(a+b)) for f
		   with the triangle (a, b, f), and sum_x (-1)^(a+b+x) (2x+1)
		   {a b x; b a f} = sqrt((2a+1)(2b+1)) delta(f, 0) for f with
		   the triangles (a, a, f) and (b, b, f) -- ----------------- */

	for (double a = 0.0; a <= 4.0; a += 0.5)
		for (double b = 0.0; b <= 4.0; b += 0.5) {
			for (double f = fabs(a - b); f <= a + b + 1e-9; f += 1.0) {
				double s1 = 0.0;

				for (double x = fabs(a - b); x <= a + b + 1e-9; x += 1.0)
					s1 += (2.0 * x + 1.0) * w6js(a, b, x, a, b, f);
				CHECK_CLOSE(s1, sign(2.0 * (a + b)), 5e-15, 5e-15, "sum rule 1, a = %g, b = %g, f = %g", a, b, f);
			}
			for (double f = 0.0; f <= 2.0 * MIN(a, b) + 1e-9; f += 1.0) {
				double s2 = 0.0;

				for (double x = fabs(a - b); x <= a + b + 1e-9; x += 1.0)
					s2 += sign(a + b + x) * (2.0 * x + 1.0) * w6js(a, b, x, b, a, f);
				CHECK_CLOSE(s2, (f == 0.0) ? sqrt((2 * a + 1) * (2 * b + 1)) : 0.0, 5e-15, 5e-15,
							"sum rule 2, a = %g, b = %g, f = %g", a, b, f);
			}
		}

	/* --- Orthogonality: sum_x (2x+1)(2f+1) {a b x; c d f} {a b x; c d f'}
		   = delta(f, f') -- ------------------------------------- */

	{
		static const double abcd[][4] = {{1.0, 1.0, 1.0, 1.0}, {1.5, 1.0, 0.5, 2.0}, {2.0, 2.5, 1.5, 3.0}, {3.0, 3.0, 3.0, 3.0}};

		for (size_t t = 0; t < sizeof(abcd) / sizeof(abcd[0]); t++) {
			double a = abcd[t][0], b = abcd[t][1], c = abcd[t][2], d = abcd[t][3];

			for (double f = 0.0; f <= 12.0; f += 0.5)
				for (double fp = 0.0; fp <= 12.0; fp += 0.5) {
					double sum = 0.0;

					for (double x = 0.0; x <= 12.0; x += 0.5)
						sum += (2.0 * x + 1.0) * (2.0 * f + 1.0) * w6js(a, b, x, c, d, f) * w6js(a, b, x, c, d, fp);
					if (f == fp) {
						/* --- delta(f, f') only when f is allowed; else 0 */
						double nonzero = 0.0;

						for (double x = 0.0; x <= 12.0; x += 0.5) nonzero += fabs(w6js(a, b, x, c, d, f));
						if (nonzero == 0.0) continue;
					}
					CHECK_CLOSE(sum, (f == fp) ? 1.0 : 0.0, 5e-15, 5e-15, "6j orthogonality a b c d = %g %g %g %g, f = %g, f' = %g", a, b, c, d, f, fp);
				}
		}
	}
	/* --- selection rules and symmetry under exchange of columns and of
		   upper and lower entries of two columns -- ------------------ */

	CHECK(w6js(1.0, 1.0, 3.0, 1.0, 1.0, 1.0) == 0.0, "triangle (j1 j2 j3), value %g", w6js(1.0, 1.0, 3.0, 1.0, 1.0, 1.0));
	CHECK(w6js(1.0, 1.0, 1.0, 1.0, 1.0, 3.0) == 0.0, "triangle (J1 J2 j3), value %g", w6js(1.0, 1.0, 1.0, 1.0, 1.0, 3.0));
	CHECK(w6js(1.0, 1.0, 1.0, 1.0, 1.0, 0.3) == 0.0, "not half-integer, value %g", w6js(1.0, 1.0, 1.0, 1.0, 1.0, 0.3));
	{
		double j1 = 2.0, j2 = 1.5, j3 = 2.5, J1 = 1.0, J2 = 2.0, J3 = 1.5, v = w6js(j1, j2, j3, J1, J2, J3);

		CHECK_CLOSE(w6js(j2, j1, j3, J2, J1, J3), v, 5e-15, 0.0, "column exchange, code %d", 0);
		CHECK_CLOSE(w6js(j1, j3, j2, J1, J3, J2), v, 5e-15, 0.0, "column exchange, code %d", 1);
		CHECK_CLOSE(w6js(J1, J2, j3, j1, j2, J3), v, 5e-15, 0.0, "exchange of upper and lower, code %d", 0);
		CHECK_CLOSE(w6js(j1, J2, J3, J1, j2, j3), v, 5e-15, 0.0, "exchange of upper and lower, code %d", 1);
	}
	if (table) printf("w6js {a b c; 0 c b}: max relative error %.2e\n", worst);
}

int main(int argc, char *argv[])
{
	table = (argc > 1 && strcmp(argv[1], "--table") == 0);

	test_expint();
	test_planck_reference();
	test_planck_limits();
	test_planck_stefan_boltzmann();
	test_w3js_j_j_0();
	test_w3js_zero();
	test_w3js_stretched();
	test_w3js_values();
	test_w3js_orthogonality();
	test_w3js_symmetry();
	test_w3js_selection();
	test_w3js_large();
	test_w6js();

	return check_summary("test_functions");
}
