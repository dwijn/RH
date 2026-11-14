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
	   --                                              -------------- */

#ifndef NEGTAU_H
#define NEGTAU_H

#include <stddef.h> /* size_t; bool_t comes from rh.h, included first */

/* --- Ten times short of where exp overflows: the threshold only has
	   to separate the harmless negative opacities of converged and
	   starting solutions (runs no deeper than -0.1 on the S-5 set) from
	   an overflow --                                  -------------- */

#define NEG_TAU_LIMIT -50.0

double negTauRun(int Ninterval, const double *dtau, int *kfirst, int *klast);
bool_t negTauCheck(int nspect, int mu, int to_obs, int Ninterval,
				   const double *dtau, long column);
bool_t negTauColumns(int nspect, int mu, int to_obs, int Nz, long Ncol,
					 const double *chi, const double *z, double muz);

void   negTauTrial(bool_t trial);
bool_t negTauFound(void);
void   negTauReset(void);
void   negTauMessage(char *message, size_t size);
void   negTauStop(void);

#endif /* !NEGTAU_H */

/* ------- end ---------------------------- negtau.h ---------------- */
