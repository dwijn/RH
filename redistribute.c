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

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere  atmos;
extern Spectrum    spectrum;
extern InputData   input;
extern CommandLine commandline;
extern char        messageStr[];

/* ------- begin -------------------------- PRDDamp.c --------------- */

/* --- Damped update of rho after PRDScatter or PRDAngleScatter
	   (N-38): rho <- old + PRD_DAMPING (rho - old), with old the N
	   values of rho before the scatter integral. Damped and undamped
	   sub-iterations have the same fixed points. Undamped, a
	   sub-iteration can flip the sign of the change of rho, and an
	   odd number of them carries the flip into the main iteration:
	   plain G on 2d_pert_h6_ca_mg with the Mg II atom of before B-44
	   ends in a period-2 cycle at 4e-6 --            -------------- */

void PRDDamp(AtomicLine *line, const double *old, long N)
{
	const char    routineName[] = "PRDDamp";
	const double  w             = input.PRD_damping;
	double       *rho           = line->rho_prd[0];
	static bool_t told          = FALSE;

	if (!told) {
		sprintf(messageStr, " PRD sub-iterations: PRD_DAMPING = %g\n", w);
		Error(MESSAGE, routineName, messageStr);
		told = TRUE;
	}
	if (w == 1.0) return;
	for (long k = 0; k < N; k++) rho[k] = old[k] + w * (rho[k] - old[k]);
}
/* ------- end ---------------------------- PRDDamp.c --------------- */

/* ------- begin -------------------------- Redistribute.c ---------- */

void Redistribute(int NmaxIter, double iterLimit)
{
	const char   routineName[] = "Redistribute";
	register int kr, nact;

	bool_t      quiet, accel, eval_operator, redistribute;
	int         niter, Nlamu;
	long        NlamuL;
	double      drho, drhomax, drhomaxa, *old = NULL;
	Atom       *atom;
	AtomicLine *line;

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];

		/* --- Initialize structures for Ng acceleration PRD iteration -- */

		for (kr = 0; kr < atom->Nline; kr++) {
			line = &atom->line[kr];
			if (line->PRD && line->Ng_prd == NULL) {
				NlamuL = (input.PRD_angle_dep == PRD_ANGLE_DEP)
							 ? 2L * atmos.Nrays * line->Nlambda * atmos.Nspace
							 : (long) line->Nlambda * atmos.Nspace;
				if (NlamuL > INT_MAX) {
					sprintf(messageStr,
							"PRD line %d -> %d: %ld values of rho, too many for"
							" the Ng acceleration of the PRD iteration",
							line->j, line->i, NlamuL);
					Error(ERROR_LEVEL_2, routineName, messageStr);
				}
				Nlamu = (int) NlamuL;

				line->Ng_prd = NgInit(Nlamu, input.PRD_Ngdelay,
									  input.PRD_Ngorder, input.PRD_Ngperiod,
									  line->rho_prd[0]);
			}
		}
	}
	/* --- Iterate over scattering integral while keeping populations
		   fixed --                                      -------------- */

	niter = 1;
	while (niter <= NmaxIter) {

		drhomaxa = 0.0;
		for (nact = 0; nact < atmos.Nactiveatom; nact++) {
			atom = atmos.activeatoms[nact];

			drhomax = 0.0;
			for (kr = 0; kr < atom->Nline; kr++) {
				line = &atom->line[kr];
				if (line->PRD) {
					if (input.PRD_damping != 1.0) {
						old = (double *) realloc(old, (size_t) line->Ng_prd->N * sizeof(double));
						memcpy(old, line->rho_prd[0], (size_t) line->Ng_prd->N * sizeof(double));
					}
					if (input.PRD_angle_dep == PRD_ANGLE_DEP)
						PRDAngleScatter(line);
					else
						PRDScatter(line);
					PRDDamp(line, old, line->Ng_prd->N);

					accel = Accelerate(line->Ng_prd, line->rho_prd[0]);
					sprintf(messageStr, "  PRD: iter #%d, atom %s, line %d,",
							line->Ng_prd->count - 1, atom->ID, kr);
					drho = MaxChange(line->Ng_prd, messageStr, quiet = FALSE);
					sprintf(messageStr, "%s", NgStatus(line->Ng_prd, accel));
					Error(MESSAGE, routineName, messageStr);

					drhomax = MAX(drho, drhomax);
				}
				drhomaxa = MAX(drhomax, drhomaxa);
			}
		}
		/* --- Solve transfer equation with fixed populations -- -------- */

		solveSpectrum(eval_operator = FALSE, redistribute = TRUE);

		if (drhomaxa < iterLimit) break;
		niter++;
	}
	free(old);
}
/* ------- end ---------------------------- Redistribute.c ---------- */
