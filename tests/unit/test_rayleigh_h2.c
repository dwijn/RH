/* ------- file: -------------------------- test_rayleigh_h2.c ------

	   Rayleigh scattering by H2 (hydrogen.c: Rayleigh_H2), the table
	   of Victor & Dalgarno (1969) from 121.57 to 632.8 nm and the
	   formula of Tarafdar & Vardya (1973) above it:

	   -- the table nodes are reproduced;
	   -- between nodes log sigma is linear in log lambda (N-48): at
		  the geometric midpoint of each interval sigma is the
		  geometric mean of the nodes, where linear interpolation in
		  sigma was up to 14 % higher (1.138 at the midpoint of
		  435.96-546.23 nm, 1.143 at 497 nm);
	   -- sigma falls monotonically and stays positive;
	   -- the table and the formula meet at 632.8 nm to 0.3 %;
	   -- zero (FALSE) below 121.57 nm; scatt = sigma n(H2).
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "background.h"
#include "constant.h"
#include "check.h"
#include "state.h"

#define NT 21

static const double lam[NT] = {
	121.57, 130.00, 140.00, 150.00, 160.00, 170.00, 185.46,
	186.27, 193.58, 199.05, 230.29, 237.91, 253.56, 275.36,
	296.81, 334.24, 404.77, 407.90, 435.96, 546.23, 632.80};

static const double sig[NT] = {
	2.35E-06, 1.22E-06, 6.80E-07, 4.24E-07, 2.84E-07, 2.00E-07, 1.25E-07,
	1.22E-07, 1.00E-07, 8.70E-08, 4.29E-08, 3.68E-08, 2.75E-08, 1.89E-08,
	1.36E-08, 8.11E-09, 3.60E-09, 3.48E-09, 2.64E-09, 1.04E-09, 5.69E-10};

static double sigma_at(double lambda)
{
	double s = -1.0;

	CHECK(Rayleigh_H2(rh, lambda, &s), "Rayleigh_H2 FALSE at %g nm", lambda);
	return s / MEGABARN_TO_M2; /* n(H2) = 1 */
}

int main(void)
{
	Molecule H2  = {0};
	double   nH2 = 1.0, prev = HUGE_VAL, worst_lin = 0.0, s;
	int      ok_mono = 1;

	H2.n             = &nH2;
	rh->atmos.H2     = &H2;
	rh->atmos.Nspace = 1;

	for (int i = 0; i < NT; i++)
		CHECK_CLOSE(sigma_at(lam[i]), sig[i], 1.0e-12, 0.0, "node %g nm", lam[i]);

	for (int i = 0; i < NT - 1; i++) {
		double lm = sqrt(lam[i] * lam[i + 1]), want = sqrt(sig[i] * sig[i + 1]),
			   lin = sig[i] + (sig[i + 1] - sig[i]) * (lm - lam[i]) / (lam[i + 1] - lam[i]);

		CHECK_CLOSE(sigma_at(lm), want, 1.0e-12, 0.0, "geometric midpoint %g nm", lm);
		worst_lin = MAX(worst_lin, lin / want);
	}
	for (double l = 121.57; l <= 2000.0; l *= 1.001) {
		s = sigma_at(l);
		if (!(s > 0.0 && s <= prev)) ok_mono = 0;
		prev = s;
	}
	CHECK(ok_mono, "sigma not positive and falling from 121.57 nm to 2 micron");
	CHECK(worst_lin > 1.13 && worst_lin < 1.15, "linear / log-log at the midpoints, largest %.4f",
		  worst_lin);

	double above = sigma_at(632.80 * (1.0 + 1.0e-12));
	CHECK(fabs(above / sig[NT - 1] - 1.0) < 3.0e-3, "table and formula at 632.8 nm: %.5f",
		  above / sig[NT - 1]);

	s = -1.0;
	CHECK(!Rayleigh_H2(rh, 121.0, &s) && s == -1.0, "below 121.57 nm");

	return check_summary("test_rayleigh_h2");
}
/* ------- end ---------------------------- test_rayleigh_h2.c ------ */
