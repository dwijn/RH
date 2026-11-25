/* ------- file: -------------------------- test_negtau.c -----------

	   Detector for negative total opacity (N-40): the most negative
	   run of consecutive negative interval optical depths, the
	   threshold, the column bound of the plane-parallel geometries,
	   and what is recorded. The interval optical depths of the Bezier
	   solvers (B-68): the 1-D cubic on the column where NEWTON_KRYLOV
	   overflowed on 1d_shock_h6 (the value Piecewise_Bezier3_1D had,
	   from gdb), and both Bezier forms between the smallest and
	   largest chi of the interval times its path length.
	   --                                              -------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "check.h"
#include "negtau.h"
#include "state.h"

static void check_runs(void)
{
	int    kf, kl;
	double t;

	/* --- No negative interval: 0 and no run --       -------------- */
	{
		const double d[4] = {1.0, 2.0, 0.0, 3.0};

		t = negTauRun(4, d, &kf, &kl);
		CHECK(t == 0.0 && kf == -1 && kl == -1, "all positive: %g [%d, %d]", t, kf, kl);
	}
	/* --- Two stretches separated by a large positive interval; the
		   second is deeper and spans the +1 --          -------------- */
	{
		const double d[8] = {-1.0, -2.0, 5.0, -0.5, -1.0, -2.0, 1.0, -3.0};

		t = negTauRun(8, d, &kf, &kl);
		CHECK(t == -5.5 && kf == 3 && kl == 7, "two runs: %g [%d, %d], want -5.5 [3, 7]", t,
			  kf, kl);
	}
	/* --- Zero and slightly positive intervals do not stop the growth
		   of the intensity: the stretch spans them (code review) -- */
	{
		const double d[5] = {-30.0, 0.0, -30.0, 0.5, -30.0};

		t = negTauRun(5, d, &kf, &kl);
		CHECK(t == -89.5 && kf == 0 && kl == 4, "gaps: %g [%d, %d], want -89.5 [0, 4]", t, kf,
			  kl);
	}
	/* --- A run up to the last interval --            -------------- */
	{
		const double d[3] = {1.0, -2.0, -2.5};

		t = negTauRun(3, d, &kf, &kl);
		CHECK(t == -4.5 && kf == 1 && kl == 2, "run at the end: %g [%d, %d]", t, kf, kl);
	}
}

static void check_threshold(void)
{
	double d[3] = {-20.0, -20.0, -9.0};

	/* --- -49 does not fire, -51 does, and the deepest detection is
		   kept --                                        -------------- */

	negTauReset();
	CHECK(!negTauCheck(3, 0, TRUE, 3, d, -1) && !negTauFound(), "-49 fired");
	d[2] = -11.0;
	CHECK(negTauCheck(3, 1, TRUE, 3, d, -1) && negTauFound(), "-51 did not fire");
	d[2] = -30.0;
	CHECK(negTauCheck(5, 2, FALSE, 3, d, -1), "-70 did not fire");
	d[2] = -10.5;
	CHECK(negTauCheck(7, 0, TRUE, 3, d, -1), "-50.5 did not fire");
	{
		char msg[256];

		negTauMessage(msg, sizeof(msg));
		CHECK(strstr(msg, "wavelength 5") && strstr(msg, "ray 2 downward") &&
				  strstr(msg, "-70") && strstr(msg, "depth points 0 and 3"),
			  "message does not name the deepest run: %s", msg);
	}
	negTauReset();
	CHECK(!negTauFound(), "reset");
}

static void check_columns(void)
{
	/* --- 2 columns, 4 heights 1 apart, muz = 0.5: column 1 is
		   negative in the middle, column 0 never. Interval optical
		   depths along column 1: (1 + -30)/2 * 2 = -29, (-30 + -30)/2
		   * 2 = -60, then (-30 + 1)/2 * 2 = -29 --       -------------- */

	const double z[4]   = {3.0, 2.0, 1.0, 0.0};
	double       chi[8] = {1.0, 1.0, 1.0, -30.0, 1.0, -30.0, 1.0, 1.0};
	char         msg[256];

	negTauReset();
	CHECK(negTauColumns(4, 1, TRUE, 4, 2, chi, z, 0.5, NEGTAU_LINEAR), "column bound did not fire");
	negTauMessage(msg, sizeof(msg));
	CHECK(strstr(msg, "-118") && strstr(msg, "column 1"), "column message: %s", msg);

	/* --- Vertical (muz = 1): -14.5 - 30 - 14.5 = -59 still fires;
		   positive opacity everywhere never does --      -------------- */

	negTauReset();
	CHECK(negTauColumns(4, 0, TRUE, 4, 2, chi, z, 1.0, NEGTAU_LINEAR), "vertical: -59 did not fire");
	for (int n = 0; n < 8; n++) chi[n] = fabs(chi[n]);
	negTauReset();
	CHECK(!negTauColumns(4, 0, TRUE, 4, 2, chi, z, 0.1, NEGTAU_LINEAR) && !negTauFound(),
		  "positive opacity fired");

	/* --- Each column checked on its own (code review): negative
		   patches in different columns at different heights, -14.5 deep
		   in each column (vertical), do not add up to a stretch no ray
		   follows; a bound over the columns would give -29 + ... --- */
	{
		double c2[8] = {-14.5, 1.0, -14.5, 1.0, 1.0, -14.5, 1.0, -14.5};

		negTauReset();
		CHECK(!negTauColumns(4, 0, TRUE, 4, 2, c2, z, 0.5, NEGTAU_LINEAR) && !negTauFound(),
			  "patches in different columns fired");
	}
	negTauReset();
}

static void check_bezier(void)
{
	/* --- 1d_shock_h6, Newton step 28 trial, wavelength 141, ray 0
		   incoming, depths 40-47 (B-68): chi changes sign at depth 43.
		   The trapezoid of interval 43-44 is positive, the solver's
		   cubic Bezier gives dtau_uw = -16179.700117111772 there --- */

	const double chi[8] = {0.44360712188715712, 0.68639103822200853, 1.0416570617988259,
						   -0.026253741537111054, 0.027579353482152538, 0.70205974842824337,
						   5.2578860339143683, 7.7587051857347085};
	const double z[8]   = {1452762.3336470015, 1348498.8796282359, 1255861.4158466596,
						   1159446.9315932011, 1068105.4240777378, 977318.19180535316,
						   908847.29709734255, 849773.28763613408};
	const double muz    = 1.0 / 21.317381323979454;
	double       dtau[7];

	negTauIntervals(NEGTAU_LINEAR, 8, 1, 0, chi, z, muz, dtau);
	CHECK(dtau[3] > 1.0E3, "trapezoid of the sign change: %g, expected +1.3e3", dtau[3]);
	negTauIntervals(NEGTAU_BEZIER3_1D, 8, 1, 0, chi, z, muz, dtau);
	CHECK(fabs(dtau[3] / -16179.700117111772 - 1.0) < 1.0E-12,
		  "cubic Bezier of the sign change: %.17g, the solver had -16179.700117111772",
		  dtau[3]);

	negTauReset();
	CHECK(!negTauColumns(141, 0, FALSE, 8, 1, chi, z, muz, NEGTAU_LINEAR) && !negTauFound(),
		  "trapezoid fired on the B-68 column");
	negTauReset();
	CHECK(negTauColumns(141, 0, FALSE, 8, 1, chi, z, muz, NEGTAU_BEZIER3_1D),
		  "cubic Bezier did not fire on the B-68 column");
	negTauReset();
	CHECK(negTauColumns(141, 0, FALSE, 8, 1, chi, z, muz, NEGTAU_BEZIER_SC),
		  "quadratic Bezier did not fire on the B-68 column");
	negTauReset();

	/* --- Random columns (2 of them, uneven spacing, values of both
		   signs over 6 decades): every Bezier interval lies between
		   the smallest and the largest chi of its ends times its path
		   length, the property the fast path (no negative chi, no
		   check) relies on --                           -------------- */
	{
		enum { NZ = 40,
			   NC = 2 };
		unsigned long long r = 88172645463325252ULL;
		double             c[NZ * NC], zz[NZ], d[NZ - 1], lo, hi, h, worst = 0.0;

		for (int trial = 0; trial < 200; trial++) {
			zz[0] = 0.0;
			for (int k = 0; k < NZ; k++) {
				for (int col = 0; col < NC; col++) {
					r ^= r << 13;
					r ^= r >> 7;
					r ^= r << 17;
					double u = (r >> 11) * (1.0 / 9007199254740992.0);

					c[k * NC + col] = ((u < 0.2) ? -1.0 : 1.0) * pow(10.0, 6.0 * u - 3.0);
				}
				if (k > 0) {
					r ^= r << 13;
					r ^= r >> 7;
					r ^= r << 17;
					zz[k] = zz[k - 1] - (0.1 + (r >> 11) * (1.0 / 9007199254740992.0));
				}
			}
			for (int q = NEGTAU_BEZIER3_1D; q <= NEGTAU_BEZIER_SC; q++)
				for (int col = 0; col < NC; col++) {
					negTauIntervals((enum negtau_quadrature) q, NZ, NC, col, c, zz, 0.3, d);
					for (int k = 0; k < NZ - 1; k++) {
						h     = fabs(zz[k] - zz[k + 1]) / 0.3;
						lo    = fmin(c[k * NC + col], c[(k + 1) * NC + col]) * h;
						hi    = fmax(c[k * NC + col], c[(k + 1) * NC + col]) * h;
						worst = fmax(worst, fmax(lo - d[k], d[k] - hi) / (hi - lo));
					}
				}
		}
		CHECK(worst < 1.0E-12, "a Bezier interval outside [min, max] chi ds (by %.2e)", worst);
	}
}

int main(void)
{
	double lambda[8] = {100, 200, 300, 400, 500, 600, 700, 800};

	rh->spectrum.lambda = lambda;
	rh->spectrum.Nspect = 8;

	check_runs();
	check_threshold();
	check_columns();
	check_bezier();

	return check_summary("test_negtau");
}
/* ------- end ---------------------------- test_negtau.c ----------- */
