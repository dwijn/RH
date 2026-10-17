/* ------- file: -------------------------- w3.c --------------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Fri Jan  4 13:30:36 2013 --

	   --------------------------                      ----------RH-- */

/* --- Functions for piecewise quadratic integration -- ------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rh.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* --- Integration weights of the short-characteristics solvers over
	   an interval of optical depth dtau:

	   w[0] = 1 - exp(-dtau)                          = P(1, dtau),
	   w[1] = 1 - exp(-dtau) (1 + dtau)               = P(2, dtau),
	   w[2] = 2 - exp(-dtau) (2 + 2 dtau + dtau^2)    = 2 P(3, dtau),

	   with P(m, t) = 1 - exp(-t) sum_{k<m} t^k/k! the regularized
	   incomplete gamma function. For dtau < 1 the closed forms lose
	   digits to cancellation (w[1] ~ dtau^2/2, w[2] ~ dtau^3/3), so P
	   is summed from its series there (also for small negative dtau,
	   which interpolated opacities can give in 2-D and 3-D). Above
	   dtau = 50 the weights are 1, 1 and 2 to double precision, also
	   for infinite dtau --                            -------------- */

#define W_SERIES_MAX 1.0

/* --- P(m, t) = exp(-t) sum_{k>=m} t^k/k! for |t| < 1 and m <= 3,
	   given expt = exp(-t). For t > 0 all terms are positive, accurate
	   to rounding. They drop below 1e-17 of the sum within 19 terms,
	   so k stays below 32. P(m - 1) = P(m) + exp(-t) t^(m-1)/(m-1)!
	   then adds positive terms only --                -------------- */

static const double inv_k[32] = {
	0.0, 1.0, 1.0 / 2, 1.0 / 3, 1.0 / 4, 1.0 / 5, 1.0 / 6, 1.0 / 7,
	1.0 / 8, 1.0 / 9, 1.0 / 10, 1.0 / 11, 1.0 / 12, 1.0 / 13, 1.0 / 14,
	1.0 / 15, 1.0 / 16, 1.0 / 17, 1.0 / 18, 1.0 / 19, 1.0 / 20, 1.0 / 21,
	1.0 / 22, 1.0 / 23, 1.0 / 24, 1.0 / 25, 1.0 / 26, 1.0 / 27, 1.0 / 28,
	1.0 / 29, 1.0 / 30, 1.0 / 31};

static double gammaP_series(int m, double t, double expt)
{
	double term = 1.0, sum = 0.0;
	int    k;

	for (k = 1; k <= m; k++)
		term *= t * inv_k[k];
	for (; k < 32 && fabs(term) > 1.0e-17 * fabs(sum); k++) {
		sum += term;
		term *= t * inv_k[k];
	}
	return expt * sum;
}

/* ------- begin -------------------------- w2.c -------------------- */

void w2(double dtau, double *w)
{
	double expdt;

	if (fabs(dtau) < W_SERIES_MAX) {
		expdt = exp(-dtau);
		w[1]  = gammaP_series(2, dtau, expdt);
		w[0]  = w[1] + expdt * dtau;
	} else if (dtau > 50.0) {
		w[0] = w[1] = 1.0;
	} else {
		expdt = exp(-dtau);
		w[0]  = 1.0 - expdt;
		w[1]  = w[0] - dtau * expdt;
	}
}
/* ------- end ---------------------------- w2.c -------------------- */

/* ------- begin -------------------------- w3.c -------------------- */

void w3(double dtau, double *w)
{
	double expdt, P3;

	if (fabs(dtau) < W_SERIES_MAX) {
		expdt = exp(-dtau);
		P3    = gammaP_series(3, dtau, expdt);
		w[2]  = 2.0 * P3;
		w[1]  = P3 + expdt * 0.5 * SQ(dtau);
		w[0]  = w[1] + expdt * dtau;
	} else if (dtau > 50.0) {
		w[0] = w[1] = 1.0;
		w[2]        = 2.0;
	} else {
		expdt = exp(-dtau);
		w[0]  = 1.0 - expdt;
		w[1]  = w[0] - dtau * expdt;
		w[2]  = 2.0 * w[1] - SQ(dtau) * expdt;
	}
}
/* ------- end ---------------------------- w3.c -------------------- */

/* ------- begin -------------------------- U3.c -------------------- */

void U3(double dtau, double *U)
{
	double delta;

	if (dtau < 2.0E-4) {
		U[0]  = dtau * (1.0 - 0.5 * dtau);
		delta = SQ(dtau);
		U[1]  = 0.5 * delta * (1.0 - dtau / 3.0);
		delta *= dtau;
		U[2] = (delta / 3.0) * (1.0 - 0.25 * dtau);
	} else if (dtau > 200.0) {
		U[0] = 1.0;
		U[1] = dtau - U[0];
		U[2] = SQ(dtau) - 2.0 * U[1];
	} else {
		U[0] = 1.0 - exp(-dtau);
		U[1] = dtau - U[0];
		U[2] = SQ(dtau) - 2.0 * U[1];
	}
}
/* ------- end ---------------------------- U3.c -------------------- */
