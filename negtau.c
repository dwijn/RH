/* ------- file: -------------------------- negtau.c ----------------

	   Detector for negative total opacity in the formal solution (N-40),
	   see negtau.h.
	   --                                              -------------- */

#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "error.h"
#include "negtau.h"

/* --- Global variables --                             -------------- */

extern Spectrum spectrum;
extern char     messageStr[];

/* --- The most negative detection since the last reset, filled by the
	   threads of solveSpectrum under the mutex --    -------------- */

static struct {
	bool_t found;
	int    nspect, mu, kfirst, klast;
	int    to_obs;
	long   column;
	double tau;
} detected = {FALSE, -1, -1, -1, -1, -1, -1, 0.0};

static pthread_mutex_t lock  = PTHREAD_MUTEX_INITIALIZER;
static bool_t          trial = FALSE;

/* ------- begin -------------------------- negTauRun.c ------------- */

double negTauRun(int Ninterval, const double *dtau, int *kfirst, int *klast)
{
	register int k;

	int    start = 0;
	double run = 0.0, runmin = 0.0;

	/* --- The most negative sum of dtau over a contiguous stretch of
		   intervals (Kadane's minimum-subarray sum), and its first and
		   last interval (-1 if no interval is negative). Interval k lies
		   between points k and k + 1. The intensity along the ray grows
		   by exp(-dtau) per interval, so this sum, not a run of
		   negative intervals only, sets how far it can grow: a zero or
		   slightly positive interval between negative ones does not
		   stop the growth --                            -------------- */

	*kfirst = *klast = -1;
	for (k = 0; k < Ninterval; k++) {
		if (k == 0 || run > 0.0) { /* a positive sum only adds */
			run   = dtau[k];
			start = k;
		} else
			run += dtau[k];
		if (run < runmin) {
			runmin  = run;
			*kfirst = start;
			*klast  = k;
		}
	}
	return runmin;
}
/* ------- end ---------------------------- negTauRun.c ------------- */

/* ------- begin -------------------------- negTauCheck.c ----------- */

bool_t negTauCheck(int nspect, int mu, int to_obs, int Ninterval,
				   const double *dtau, long column)
{
	int    kfirst, klast;
	double tau = negTauRun(Ninterval, dtau, &kfirst, &klast);

	/* --- dtau[k] is the optical depth of interval k along ray mu;
		   to_obs is the direction (-1: none, for an angle-independent
		   check or a spherical ray, whose Feautrier solution covers
		   both) and column the column in 2-D and 3-D (-1 in 1-D and
		   spherical geometry). Records a stretch below NEG_TAU_LIMIT and
		   returns TRUE for it --                        -------------- */

	if (!(tau < NEG_TAU_LIMIT)) return FALSE;

	pthread_mutex_lock(&lock);
	if (!detected.found || tau < detected.tau) {
		detected.found  = TRUE;
		detected.nspect = nspect;
		detected.mu     = mu;
		detected.to_obs = to_obs;
		detected.kfirst = kfirst;
		detected.klast  = klast;
		detected.column = column;
		detected.tau    = tau;
	}
	pthread_mutex_unlock(&lock);

	return TRUE;
}
/* ------- end ---------------------------- negTauCheck.c ----------- */

/* ------- begin -------------------------- negTauColumns.c --------- */

bool_t negTauColumns(int nspect, int mu, int to_obs, int Nz, long Ncol,
					 const double *chi, const double *z, double muz)
{
	register long c, n;
	register int  k;

	int     kf, kl;
	long    cmin = -1;
	double *dtau = NULL, tau, taumin = 0.0;

	/* --- Plane-parallel geometries: chi[k * Ncol + c] at height z[k]
		   in column c (1-D: Ncol = 1); also a spherical ray, with its
		   path positions as z, Ncol = 1 and muz = 1. Each column is
		   checked along the vertical scaled by 1 / muz. A ray in 2-D or
		   3-D crosses columns, but the threshold is more than ten times
		   short of an overflow, and a bound over all columns (each
		   interval in the column where it is most negative) could stop
		   a run that would finish. Nothing to do where no total
		   opacity is negative, the common case: one pass over chi -- */

	for (n = 0; n < Nz * Ncol; n++)
		if (chi[n] < 0.0) break;
	if (n == Nz * Ncol) return FALSE;

	dtau = (double *) malloc((Nz - 1) * sizeof(double));
	for (c = 0; c < Ncol; c++) {
		for (k = 0; k < Nz; k++)
			if (chi[k * Ncol + c] < 0.0) break;
		if (k == Nz) continue;

		for (k = 0; k < Nz - 1; k++)
			dtau[k] = 0.5 * (chi[k * Ncol + c] + chi[(k + 1) * Ncol + c]) *
					  fabs(z[k] - z[k + 1]) / muz;
		tau = negTauRun(Nz - 1, dtau, &kf, &kl);
		if (tau < taumin) {
			taumin = tau;
			cmin   = c;
		}
	}
	if (!(taumin < NEG_TAU_LIMIT)) {
		free(dtau);
		return FALSE;
	}
	for (k = 0; k < Nz - 1; k++)
		dtau[k] = 0.5 * (chi[k * Ncol + cmin] + chi[(k + 1) * Ncol + cmin]) *
				  fabs(z[k] - z[k + 1]) / muz;
	negTauCheck(nspect, mu, to_obs, Nz - 1, dtau, (Ncol > 1) ? cmin : -1);
	free(dtau);

	return TRUE;
}
/* ------- end ---------------------------- negTauColumns.c --------- */

/* ------- begin -------------------------- negTauTrial.c ----------- */

void negTauTrial(bool_t is_trial)
{
	/* --- Called by the solvers around G_eval at a trial point, from
		   the main thread between sweeps --             -------------- */

	trial = is_trial;
}
/* ------- end ---------------------------- negTauTrial.c ----------- */

/* ------- begin -------------------------- negTauFound.c ----------- */

bool_t negTauFound(void)
{
	bool_t found;

	/* --- Also read by the threads of solveSpectrum, to skip the
		   wavelengths left once one has been found --   -------------- */

	pthread_mutex_lock(&lock);
	found = detected.found;
	pthread_mutex_unlock(&lock);

	return found;
}
/* ------- end ---------------------------- negTauFound.c ----------- */

/* ------- begin -------------------------- negTauReset.c ----------- */

void negTauReset(void)
{
	detected.found = FALSE;
	detected.tau   = 0.0;
}
/* ------- end ---------------------------- negTauReset.c ----------- */

/* ------- begin -------------------------- negTauMessage.c --------- */

void negTauMessage(char *message, size_t size)
{
	char where[64];

	if (detected.column >= 0)
		snprintf(where, sizeof(where), ", column %ld", detected.column);
	else
		where[0] = '\0';

	snprintf(message, size,
			 "negative total opacity at lambda = %.6f nm (wavelength %d),"
			 " ray %d%s: optical depth %.4g between depth points %d and %d%s"
			 " (limit %g)",
			 spectrum.lambda[detected.nspect], detected.nspect, detected.mu,
			 (detected.to_obs < 0) ? "" : (detected.to_obs) ? " upward"
															: " downward",
			 detected.tau,
			 detected.kfirst, detected.klast + 1, where, NEG_TAU_LIMIT);
}
/* ------- end ---------------------------- negTauMessage.c --------- */

/* ------- begin -------------------------- negTauStop.c ------------ */

void negTauStop(void)
{
	const char routineName[] = "negTauStop";
	char       what[256];

	/* --- Called by solveSpectrum after its threads have joined: stops
		   the run unless a solver is evaluating a trial point -- ---- */

	if (!detected.found || trial) return;

	negTauMessage(what, sizeof(what));
	sprintf(messageStr,
			"\n  %s.\n  The formal solution would overflow there. A negative"
			" total opacity\n  comes from a strong inverted line (population"
			" inversion) in\n  the current populations (SOLVER.md, trap 5)",
			what);
	Error(ERROR_LEVEL_2, routineName, messageStr);
}
/* ------- end ---------------------------- negTauStop.c ------------ */
