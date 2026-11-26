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
#include "bezier.h"
#include "state.h"

/* --- Global variables --                             -------------- */

extern char messageStr[];

/* --- The most negative detection since the last reset, filled by the
	   threads of solveSpectrum under the mutex --    -------------- */

static struct {
	bool_t found;
	int    nspect, mu, kfirst, klast;
	int    to_obs;
	long   column;
	double tau;
} detected = {FALSE, -1, -1, -1, -1, -1, -1, 0.0};

/* --- A trial point invalid for another reason (negTauInvalid), with
	   the message of the first one since the last reset -- -------- */

static bool_t invalid = FALSE;
static char   invalidMessage[512];

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
		   check) and column the column in 2-D and 3-D (-1 in 1-D).
		   Records a stretch below NEG_TAU_LIMIT and
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

/* ------- begin -------------------------- negTauIntervals.c ------- */

/* --- chi, the path length and the derivative of chi along increasing
	   k in column c: the Steffen derivative at the inner points, the
	   secant at the two ends (the 1-D Bezier solvers) -- ---------- */

static double colChi(const double *chi, long Ncol, long c, int k)
{
	return chi[(long) k * Ncol + c];
}

static double colH(const double *z, int k)
{
	return fabs(z[k] - z[k + 1]);
}

static double colDeriv(const double *chi, const double *z, int Nz, long Ncol, long c,
					   int k)
{
	if (k == 0) return (colChi(chi, Ncol, c, 1) - colChi(chi, Ncol, c, 0)) / colH(z, 0);
	if (k == Nz - 1)
		return (colChi(chi, Ncol, c, Nz - 1) - colChi(chi, Ncol, c, Nz - 2)) /
			   colH(z, Nz - 2);
	return cent_deriv(colH(z, k - 1), colH(z, k), colChi(chi, Ncol, c, k - 1),
					  colChi(chi, Ncol, c, k), colChi(chi, Ncol, c, k + 1));
}

void negTauIntervals(enum negtau_quadrature quad, int Nz, long Ncol, long c,
					 const double *chi, const double *z, double muz, double *dtau)
{
	register int k;

	double h, x0, x1, d0, d1, up, down;

	/* --- The optical depths dtau[0 .. Nz - 2] of the intervals of
		   column c, chi[k * Ncol + c] at height z[k], along a ray with
		   cosine muz, integrated as quad says (negtau.h). The cubic
		   does not depend on the direction of the ray; the quadratic
		   does, and is taken both ways. A Steffen derivative is
		   homogeneous in the path lengths, so they are taken along the
		   vertical and the result is scaled by 1 / muz -- ---------- */

	for (k = 0; k < Nz - 1; k++) {
		h  = colH(z, k);
		x0 = colChi(chi, Ncol, c, k);
		x1 = colChi(chi, Ncol, c, k + 1);

		switch (quad) {
		case NEGTAU_BEZIER3_1D:
			d0      = colDeriv(chi, z, Nz, Ncol, c, k);
			d1      = colDeriv(chi, z, Nz, Ncol, c, k + 1);
			dtau[k] = 0.25 * h * (2.0 * (x0 + x1) + h / 3.0 * (d0 - d1));
			break;

		case NEGTAU_BEZIER_SC:
			/* --- Along increasing k the local (downwind) point is
				   k + 1, along decreasing k it is k; at the end of the
				   ray the trapezoid --              -------------- */

			up = down = 0.5 * h * (x0 + x1);
			if (k + 1 < Nz - 1) {
				d1 = colDeriv(chi, z, Nz, Ncol, c, k + 1);
				up = h * (0.5 * (x0 + x1) + (x1 - x0) / 6.0 - h * d1 / 6.0);
			}
			if (k > 0) {
				d0   = colDeriv(chi, z, Nz, Ncol, c, k);
				down = h * (0.5 * (x0 + x1) + (x0 - x1) / 6.0 + h * d0 / 6.0);
			}
			dtau[k] = fmin(up, down);
			break;

		default:
			dtau[k] = 0.5 * h * (x0 + x1);
		}
		dtau[k] /= muz;
	}
}
/* ------- end ---------------------------- negTauIntervals.c ------- */

/* ------- begin -------------------------- negTauColumns.c --------- */

bool_t negTauColumns(int nspect, int mu, int to_obs, int Nz, long Ncol,
					 const double *chi, const double *z, double muz,
					 enum negtau_quadrature quad)
{
	register long c, n;
	register int  k;

	int     kf, kl;
	long    cmin = -1;
	double *dtau = NULL, tau, taumin = 0.0;

	/* --- Plane-parallel geometries: chi[k * Ncol + c] at height z[k]
		   in column c (1-D: Ncol = 1). Each column is
		   checked along the vertical scaled by 1 / muz. A ray in 2-D or
		   3-D crosses columns, but the threshold is more than ten times
		   short of an overflow, and a bound over all columns (each
		   interval in the column where it is most negative) could stop
		   a run that would finish. The intervals are integrated as
		   the caller's formal solver does (quad, B-68). Nothing to do
		   where no total opacity is negative, the common case: one
		   pass over chi (the Bezier curves keep their control points
		   between the values at the ends, so they are not negative
		   either, Q-34) --                              -------------- */

	for (n = 0; n < Nz * Ncol; n++)
		if (chi[n] < 0.0) break;
	if (n == Nz * Ncol) return FALSE;

	dtau = (double *) malloc((Nz - 1) * sizeof(double));
	for (c = 0; c < Ncol; c++) {
		for (k = 0; k < Nz; k++)
			if (chi[k * Ncol + c] < 0.0) break;
		if (k == Nz) continue;

		negTauIntervals(quad, Nz, Ncol, c, chi, z, muz, dtau);
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
	negTauIntervals(quad, Nz, Ncol, cmin, chi, z, muz, dtau);
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

/* ------- begin -------------------------- negTauInTrial.c --------- */

bool_t negTauInTrial(void)
{
	return trial;
}
/* ------- end ---------------------------- negTauInTrial.c --------- */

/* ------- begin -------------------------- negTauInvalid.c --------- */

void negTauInvalid(const char *message)
{
	/* --- Marks the trial point invalid, like a detection: negTauFound
		   is TRUE until the next negTauReset and negTauMessage gives
		   this message (the first one, unless an opacity detection
		   came first). Only during a trial, from the main thread;
		   anywhere else it is a programming error and stops -- ----                                      -------------- */

	if (!trial) {
		Error(ERROR_LEVEL_2, "negTauInvalid",
			  "called outside a trial evaluation of G (programming error)");
		return;
	}
	pthread_mutex_lock(&lock);
	if (!invalid) {
		invalid = TRUE;
		snprintf(invalidMessage, sizeof(invalidMessage), "%s", message);
	}
	pthread_mutex_unlock(&lock);
}
/* ------- end ---------------------------- negTauInvalid.c --------- */

/* ------- begin -------------------------- negTauFound.c ----------- */

bool_t negTauFound(void)
{
	bool_t found;

	/* --- Also read by the threads of solveSpectrum, to skip the
		   wavelengths left once one has been found --   -------------- */

	pthread_mutex_lock(&lock);
	found = detected.found || invalid;
	pthread_mutex_unlock(&lock);

	return found;
}
/* ------- end ---------------------------- negTauFound.c ----------- */

/* ------- begin -------------------------- negTauReset.c ----------- */

void negTauReset(void)
{
	detected.found = FALSE;
	detected.tau   = 0.0;
	invalid        = FALSE;
}
/* ------- end ---------------------------- negTauReset.c ----------- */

/* ------- begin -------------------------- negTauMessage.c --------- */

void negTauMessage(char *message, size_t size)
{
	char where[64];

	if (!detected.found && invalid) {
		snprintf(message, size, "%s", invalidMessage);
		return;
	}
	if (detected.column >= 0)
		snprintf(where, sizeof(where), ", column %ld", detected.column);
	else
		where[0] = '\0';

	snprintf(message, size,
			 "negative total opacity at lambda = %.6f nm (wavelength %d),"
			 " ray %d%s: optical depth %.4g between depth points %d and %d%s"
			 " (limit %g)",
			 rh->spectrum.lambda[detected.nspect], detected.nspect, detected.mu,
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

	if (!detected.found || trial) return; /* (invalid only in a trial) */

	negTauMessage(what, sizeof(what));
	sprintf(messageStr,
			"\n  %s.\n  The formal solution would overflow there. A negative"
			" total opacity\n  comes from a strong inverted line (population"
			" inversion) in\n  the current populations (SOLVER.md, trap 5)",
			what);
	Error(ERROR_LEVEL_2, routineName, messageStr);
}
/* ------- end ---------------------------- negTauStop.c ------------ */
