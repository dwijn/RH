/* ------- file: -------------------------- test_negtau.c -----------

	   Detector for negative total opacity (N-40): the most negative
	   run of consecutive negative interval optical depths, the
	   threshold, the column bound of the plane-parallel geometries,
	   and what is recorded.
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

extern Spectrum spectrum;

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
	CHECK(negTauColumns(4, 1, TRUE, 4, 2, chi, z, 0.5), "column bound did not fire");
	negTauMessage(msg, sizeof(msg));
	CHECK(strstr(msg, "-118") && strstr(msg, "column 1"), "column message: %s", msg);

	/* --- Vertical (muz = 1): -14.5 - 30 - 14.5 = -59 still fires;
		   positive opacity everywhere never does --      -------------- */

	negTauReset();
	CHECK(negTauColumns(4, 0, TRUE, 4, 2, chi, z, 1.0), "vertical: -59 did not fire");
	for (int n = 0; n < 8; n++) chi[n] = fabs(chi[n]);
	negTauReset();
	CHECK(!negTauColumns(4, 0, TRUE, 4, 2, chi, z, 0.1) && !negTauFound(),
		  "positive opacity fired");

	/* --- Each column checked on its own (code review): negative
		   patches in different columns at different heights, -14.5 deep
		   in each column (vertical), do not add up to a stretch no ray
		   follows; a bound over the columns would give -29 + ... --- */
	{
		double c2[8] = {-14.5, 1.0, -14.5, 1.0, 1.0, -14.5, 1.0, -14.5};

		negTauReset();
		CHECK(!negTauColumns(4, 0, TRUE, 4, 2, c2, z, 0.5) && !negTauFound(),
			  "patches in different columns fired");
	}
	negTauReset();
}

int main(void)
{
	double lambda[8] = {100, 200, 300, 400, 500, 600, 700, 800};

	spectrum.lambda = lambda;
	spectrum.Nspect = 8;

	check_runs();
	check_threshold();
	check_columns();

	return check_summary("test_negtau");
}
/* ------- end ---------------------------- test_negtau.c ----------- */
