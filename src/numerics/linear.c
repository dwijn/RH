/* ------- file: -------------------------- linear.c ----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Jun 14 09:53:15 2001 --

	   --------------------------                      ----------RH-- */

/* --- Linear interpolation routine. --                -------------- */

#include <math.h>

#include "rh.h"
#include "error.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern char messageStr[];

/* ------- begin -------------------------- Linear.c ---------------- */

void Linear(int Ntable, double *xtable, double *ytable,
			int N, double *x, double *y, bool_t hunt)
{
	register int n;

	bool_t ascend;
	int    j = 0;
	double xmin, xmax, fx;

	ascend = (xtable[1] > xtable[0]) ? TRUE : FALSE;
	xmin   = (ascend) ? xtable[0] : xtable[Ntable - 1];
	xmax   = (ascend) ? xtable[Ntable - 1] : xtable[0];

	for (n = 0; n < N; n++) {
		if (x[n] <= xmin)
			y[n] = (ascend) ? ytable[0] : ytable[Ntable - 1];
		else if (x[n] >= xmax)
			y[n] = (ascend) ? ytable[Ntable - 1] : ytable[0];
		else {
			if (hunt)
				Hunt(Ntable, xtable, x[n], &j);
			else
				Locate(Ntable, xtable, x[n], &j);

			fx   = (xtable[j + 1] - x[n]) / (xtable[j + 1] - xtable[j]);
			y[n] = fx * ytable[j] + (1 - fx) * ytable[j + 1];
		}
	}
}
/* ------- end ---------------------------- Linear.c ---------------- */

/* ------- begin -------------------------- heightAtUnityTau.c ------ */

/* --- Height z1 at which tau_ref = 1, from a table z(tau_ref) with
	   tau_ref strictly increasing: linear interpolation in ln tau_ref
	   between the two bracketing points. Points with tau_ref <= 0 (the
	   top point of a column-mass scale can have none) are not used.
	   Returns TAU_UNITY_OK, TAU_NOT_INCREASING or, when tau_ref = 1 is
	   not inside the range of the positive points, TAU_NOT_BRACKETED
	   (z1 is not set then). --                        -------------- */

int heightAtUnityTau(int N, const double *tau_ref, const double *z, double *z1)
{
	int k;

	for (k = 1; k < N; k++)
		if (!(tau_ref[k] > tau_ref[k - 1])) return TAU_NOT_INCREASING;

	for (k = 0; k < N - 1; k++) {
		if (tau_ref[k] > 0.0 && tau_ref[k] <= 1.0 && tau_ref[k + 1] >= 1.0) {
			double f = log(1.0 / tau_ref[k]) /
					   log(tau_ref[k + 1] / tau_ref[k]);
			*z1 = (1.0 - f) * z[k] + f * z[k + 1];
			return TAU_UNITY_OK;
		}
	}
	return TAU_NOT_BRACKETED;
}
/* ------- end ---------------------------- heightAtUnityTau.c ------ */

/* ------- begin -------------------------- zeroPointAtUnityTau.c --- */

/* --- As heightAtUnityTau, but stops with ERROR_LEVEL_2 and a message
	   when tau_ref is not strictly increasing or does not reach
	   tau_ref = 1. --                                 -------------- */

double zeroPointAtUnityTau(const char *routineName, const char *atmosID,
						   int N, const double *tau_ref, const double *z)
{
	double z1 = 0.0;

	switch (heightAtUnityTau(N, tau_ref, z, &z1)) {
	case TAU_NOT_INCREASING:
		sprintf(messageStr,
				"Atmosphere %s: reference optical depth tau_ref is not "
				"strictly increasing with depth index; cannot set the "
				"zero of the height scale at tau_ref = 1",
				atmosID);
		Error(ERROR_LEVEL_2, routineName, messageStr);
		break;
	case TAU_NOT_BRACKETED:
		sprintf(messageStr,
				"Atmosphere %s: tau_ref = 1 lies outside the atmosphere "
				"(tau_ref from %g to %g); cannot set the zero of the height "
				"scale. Use a GEOMETRIC scale or extend the atmosphere",
				atmosID, tau_ref[0], tau_ref[N - 1]);
		Error(ERROR_LEVEL_2, routineName, messageStr);
		break;
	default: break;
	}
	return z1;
}
/* ------- end ---------------------------- zeroPointAtUnityTau.c --- */
