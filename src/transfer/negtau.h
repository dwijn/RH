/* ------- file: -------------------------- negtau.h ----------------

	   Detector for negative total opacity in the formal solution (N-40).

	   A negative total opacity (a strong inverted line, SOLVER.md trap
	   5) makes the intensity grow as exp(-tau) across a run of negative
	   optical depth, and the formal solvers overflow when the run is
	   deeper than about -709. Each Formal() checks the opacity of every
	   wavelength and ray before solving: the most negative optical depth
	   summed over a contiguous stretch of intervals (Kadane), and
	   skips the rest of that wavelength when it is below NEG_TAU_LIMIT.

	   Outside a trial point, solveSpectrum then stops the run with a
	   message naming the wavelength, ray and depth (negTauStop). While a
	   solver evaluates G at a trial point (negTauTrial(TRUE)), the
	   detection is only recorded: G_eval returns FALSE and the solver
	   rejects the point.

	   The same mechanism invalidates a trial point for other reasons
	   (negTauInvalid, only inside a trial): negative or non-finite
	   populations out of the LU fallback of statEquil (N-8). The caller
	   outside a trial stops the run itself. negTauFound therefore means
	   "the trial point is invalid, for any reason", and negTauMessage
	   says which.
	   --                                              -------------- */

#ifndef NEGTAU_H
#define NEGTAU_H

#include <stddef.h> /* size_t; bool_t comes from rh.h, included first */

/* --- Ten times short of where exp overflows: the threshold only has
	   to separate the harmless negative opacities of converged and
	   starting solutions (runs no deeper than -0.1 on the S-5 set) from
	   an overflow --                                  -------------- */

#define NEG_TAU_LIMIT -50.0

/* --- How the formal solver of the caller integrates chi over an
	   interval (B-68): the detector takes the optical depth of each
	   interval as that solver does, since a Bezier curve through a
	   sign change of chi can make an interval negative where the
	   trapezoid is positive.

	   NEGTAU_LINEAR      the trapezoid: piecewise linear solvers and
						  Feautrier;
	   NEGTAU_BEZIER3_1D  the cubic Bezier curve of the 1-D solvers
						  (bezier_1D.c, also Stokes): Steffen
						  derivatives at the inner points, the secant
						  at the two ends;
	   NEGTAU_BEZIER_SC   the quadratic Bezier curve of the 2-D and 3-D
						  short (and long) characteristics: the Steffen
						  derivative at the downwind end of the
						  interval, the trapezoid where that end is
						  the end of the ray. A column stands for the
						  ray (see negTauColumns), so both directions
						  are taken and the smaller kept -- --------- */

enum negtau_quadrature { NEGTAU_LINEAR,
						 NEGTAU_BEZIER3_1D,
						 NEGTAU_BEZIER_SC };

double negTauRun(int Ninterval, const double *dtau, int *kfirst, int *klast);
bool_t negTauCheck(RHState *rh, int nspect, int mu, int to_obs, int Ninterval,
				   const double *dtau, long column);
void   negTauIntervals(enum negtau_quadrature quad, int Nz, long Ncol, long c,
					   const double *chi, const double *z, double muz, double *dtau);
bool_t negTauColumns(RHState *rh, int nspect, int mu, int to_obs, int Nz, long Ncol,
					 const double *chi, const double *z, double muz,
					 enum negtau_quadrature quad);

void   negTauTrial(RHState *rh, bool_t trial);
bool_t negTauInTrial(RHState *rh);
void   negTauInvalid(RHState *rh, const char *message);
bool_t negTauFound(RHState *rh);
void   negTauReset(RHState *rh);
void   negTauMessage(RHState *rh, char *message, size_t size);
void   negTauStop(RHState *rh);

/* --- The state of the detector between calls, its part of RHState
	   (rh->negtau, Q-46). The threads of solveSpectrum write detected
	   and read invalid under the mutex; reset and the messages run in
	   the main thread between sweeps -- --------------- */

typedef struct {
	bool_t found;
	int    nspect, mu, kfirst, klast;
	int    to_obs;
	long   column;
	double tau;
} NegTauDetection;

typedef struct {
	NegTauDetection detected;            /* the most negative since the last reset */
	bool_t          invalid;             /* a trial point invalid for another reason */
	char            invalidMessage[512]; /* the message of the first one */
	pthread_mutex_t lock;
	bool_t          trial; /* G is evaluated at a trial point */
} NegTauState;

#define NEGTAU_STATE_INIT {.detected = {FALSE, -1, -1, -1, -1, -1, -1, 0.0}, \
						   .invalid  = FALSE,                                \
						   .lock     = PTHREAD_MUTEX_INITIALIZER,            \
						   .trial    = FALSE}

#endif /* !NEGTAU_H */

/* ------- end ---------------------------- negtau.h ---------------- */
