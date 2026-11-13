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

	/* --- The most negative sum of dtau over a run of consecutive
		   intervals with dtau < 0, and the first and last interval of
		   that run (-1 if there is none). Interval k lies between
		   points k and k + 1 --                         -------------- */

	*kfirst = *klast = -1;
	for (k = 0; k < Ninterval; k++) {
		if (dtau[k] < 0.0) {
			if (run == 0.0) start = k;
			run += dtau[k];
			if (run < runmin) {
				runmin  = run;
				*kfirst = start;
				*klast  = k;
			}
		} else
			run = 0.0;
	}
	return runmin;
}
/* ------- end ---------------------------- negTauRun.c ------------- */

/* ------- begin -------------------------- negTauCheck.c ----------- */

bool_t negTauCheck(int nspect, int mu, int to_obs, int Ninterval,
				   const double *dtau, const long *column)
{
	int    kfirst, klast;
	double tau = negTauRun(Ninterval, dtau, &kfirst, &klast);

	/* --- dtau[k] is the optical depth of interval k along the ray
		   (for 2-D and 3-D the most negative over the columns, which
		   column[k] names; NULL in 1-D and spherical geometry). to_obs
		   is the direction, -1 for none (a spherical ray, whose
		   Feautrier solution covers both). Records
		   a run below NEG_TAU_LIMIT and returns TRUE for it -- ----- */

	if (!(tau < NEG_TAU_LIMIT)) return FALSE;

	pthread_mutex_lock(&lock);
	if (!detected.found || tau < detected.tau) {
		detected.found  = TRUE;
		detected.nspect = nspect;
		detected.mu     = mu;
		detected.to_obs = to_obs;
		detected.kfirst = kfirst;
		detected.klast  = klast;
		detected.column = (column) ? column[kfirst] : -1;
		detected.tau    = tau;
	}
	pthread_mutex_unlock(&lock);

	return TRUE;
}
/* ------- end ---------------------------- negTauCheck.c ----------- */

/* ------- begin -------------------------- negTauColumns.c --------- */

bool_t negTauColumns(int nspect, int mu, bool_t to_obs, int Nz, long Ncol,
					 const double *chi, const double *z, double muz)
{
	register long c, n;
	register int  k;

	bool_t  hit;
	long   *column;
	double *dtau, d;

	/* --- Plane-parallel geometries: chi[k * Ncol + c] at height z[k]
		   in column c (1-D: Ncol = 1). The optical depth of interval k
		   along a ray with direction cosine muz is taken in the column
		   where it is most negative, a bound for a ray that crosses
		   columns in 2-D and 3-D. Nothing to do where no total opacity
		   is negative, the common case: one pass over chi -- ------- */

	for (n = 0; n < Nz * Ncol; n++)
		if (chi[n] < 0.0) break;
	if (n == Nz * Ncol) return FALSE;

	dtau   = (double *) malloc((Nz - 1) * sizeof(double));
	column = (long *) malloc((Nz - 1) * sizeof(long));
	for (k = 0; k < Nz - 1; k++) {
		dtau[k]   = HUGE_VAL;
		column[k] = 0;
		for (c = 0; c < Ncol; c++) {
			d = chi[k * Ncol + c] + chi[(k + 1) * Ncol + c];
			if (d < dtau[k]) {
				dtau[k]   = d;
				column[k] = c;
			}
		}
		dtau[k] *= 0.5 * fabs(z[k] - z[k + 1]) / muz;
	}
	hit = negTauCheck(nspect, mu, to_obs, Nz - 1, dtau, (Ncol > 1) ? column : NULL);
	free(dtau);
	free(column);

	return hit;
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
	return detected.found;
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
			 " ray %d%s: optical depth %.4g over depth intervals %d-%d%s"
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
