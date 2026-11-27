/* ------- file: -------------------------- redistribute.c ----------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Mon Feb 26 16:40:54 2024 --

	   --------------------------                      ----------RH-- */

/* --- Administers iterations to redistribute intensity in PRD line while
	   keeping the population number fixed. --         -------------- */

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "error.h"
#include "inputs.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- PRDDamp.c --------------- */

/* --- Damped update of rho after PRDScatter or PRDAngleScatter
	   (N-38): rho <- old + PRD_DAMPING (rho - old), with old the N
	   values of rho before the scatter integral. Damped and undamped
	   sub-iterations have the same fixed points. Undamped, a
	   sub-iteration can flip the sign of the change of rho, and an
	   odd number of them carries the flip into the main iteration:
	   plain G on 2d_pert_h6_ca_mg with the Mg II atom of before B-44
	   ends in a period-2 cycle at 4e-6.

	   Returns the largest relative change of the undamped update,
	   max |rho - old| / |rho| over rho != 0 (the residual of the
	   sub-iteration, which the PRD_ITER_LIMIT test and the log use;
	   the damped step is PRD_DAMPING times it), or 0 without old --- */

double PRDDamp(RHState *rh, AtomicLine *line, const double *old, long N)
{
	const double w   = rh->input.PRD_damping;
	double      *rho = line->rho_prd[0], dmax = 0.0;

	if (!rh->redistribute.dampingNoted) {
		Message(rh, " PRD sub-iterations: PRD_DAMPING = %g\n", w);
		rh->redistribute.dampingNoted = TRUE;
	}
	if (old == NULL) return 0.0;
	for (long k = 0; k < N; k++) {
		if (rho[k]) dmax = MAX(dmax, fabs((rho[k] - old[k]) / rho[k]));
	}
	if (w != 1.0) {
		for (long k = 0; k < N; k++) rho[k] = old[k] + w * (rho[k] - old[k]);
	}
	return dmax;
}
/* ------- end ---------------------------- PRDDamp.c --------------- */

/* ------- begin -------------------------- Redistribute.c ---------- */

void Redistribute(RHState *rh, int NmaxIter, double iterLimit)
{
	const char   routineName[] = "Redistribute";
	register int kr, nact;

	char        label[MAX_MESSAGE_LENGTH], status[MAX_MESSAGE_LENGTH];
	bool_t      quiet, accel, eval_operator, redistribute;
	int         niter, Nlamu;
	long        NlamuL;
	double      drho, dres = 0.0, drhomax, drhomaxa, *old = NULL;
	Atom       *atom;
	AtomicLine *line;

	for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
		atom = rh->atmos.activeatoms[nact];

		/* --- Initialize structures for Ng acceleration PRD iteration -- */

		for (kr = 0; kr < atom->Nline; kr++) {
			line = &atom->line[kr];
			if (line->PRD && line->Ng_prd == NULL) {
				NlamuL = PRDRhoSize(rh, line);
				if (NlamuL > INT_MAX) {
					Error(rh, ERROR_LEVEL_2, routineName, "PRD line %d -> %d: %ld values of rho, too many for"
														  " the Ng acceleration of the PRD iteration",
						  line->j, line->i, NlamuL);
				}
				Nlamu = (int) NlamuL;

				line->Ng_prd = NgInit(rh, Nlamu, rh->input.PRD_Ngdelay,
									  rh->input.PRD_Ngorder, rh->input.PRD_Ngperiod,
									  line->rho_prd[0]);
			}
		}
	}
	/* --- Iterate over scattering integral while keeping populations
		   fixed --                                      -------------- */

	niter = 1;
	while (niter <= NmaxIter) {

		drhomaxa = 0.0;
		for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
			atom = rh->atmos.activeatoms[nact];

			drhomax = 0.0;
			for (kr = 0; kr < atom->Nline; kr++) {
				line = &atom->line[kr];
				if (line->PRD) {
					if (rh->input.PRD_damping != 1.0) {
						old = (double *) realloc(old, (size_t) line->Ng_prd->N * sizeof(double));
						memcpy(old, line->rho_prd[0], (size_t) line->Ng_prd->N * sizeof(double));
					}
					if (rh->input.PRD_angle_dep == PRD_ANGLE_DEP)
						PRDAngleScatter(rh, line);
					else
						PRDScatter(rh, line);
					dres = PRDDamp(rh, line, old, line->Ng_prd->N);

					/* --- Damped, the change of rho is PRD_DAMPING times
						   the residual; the log and the PRD_ITER_LIMIT
						   test use the residual (N-38) --  -------------- */

					ExitOnError(rh, Accelerate(rh, line->Ng_prd, line->rho_prd[0], &accel));

					/* --- Stops the run if rho <= 0 or not finite (N-4).
						   SOLVER = NG has no trial evaluations, so a
						   FALSE return would be a bug -- ------------- */

					if (!PRDRhoCheck(rh, atom, kr))
						Error(rh, ERROR_LEVEL_2, routineName,
							  "rho of a PRD line is invalid in a trial evaluation");
					snprintf(label, sizeof(label),
							 "  PRD: iter #%d, atom %s, line %d,",
							 line->Ng_prd->count - 1, atom->ID, kr);
					if (rh->input.PRD_damping == 1.0)
						drho = MaxChange(rh, line->Ng_prd, label, quiet = FALSE);
					else {
						/* --- After an accepted Ng step, its jump counts
							   too (code review) --          -------------- */

						drho = accel ? MAX(dres, MaxChange(rh, line->Ng_prd, label, quiet = TRUE))
									 : dres;
						fprintf(rh->commandline.logfile, "%s delta = %6.4E", label, drho);
					}
					Message(rh, "%s", NgStatus(line->Ng_prd, accel, status, sizeof(status)));

					drhomax = MAX(drho, drhomax);
				}
				drhomaxa = MAX(drhomax, drhomaxa);
			}
		}
		/* --- Solve transfer equation with fixed populations -- -------- */

		solveSpectrum(rh, eval_operator = FALSE, redistribute = TRUE);

		if (drhomaxa < iterLimit) break;
		niter++;
	}
	free(old);
}
/* ------- end ---------------------------- Redistribute.c ---------- */
