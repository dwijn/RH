/* ------- file: -------------------------- initscatter.c -----------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Feb 20 13:35:29 2024 --

	   --------------------------                      ----------RH-- */

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"
#include "bgscatter.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern char messageStr[];

/* ------- begin -------------------------- initScatter.c ----------- */

void initScatter()
{
	const char   routineName[] = "initScatter";
	register int n, kr, niter;

	double      dJmax;
	Atom       *atom;
	AtomicLine *line;

	/* --- Fill the radiative rates from previous solution if
		   PRD lines are present --                      -------------- */

	if (rh->atmos.NPRDactive > 0 && rh->input.startJ == OLD_J) {
		for (n = 0; n < rh->atmos.Natom; n++) {
			atom = &rh->atmos.atoms[n];
			if (atom->active) {
				readRadRate(atom);
				for (kr = 0; kr < atom->Nline; kr++) {
					line = atom->line + kr;
					if (line->PRD) {
						if (rh->input.PRD_angle_dep == PRD_ANGLE_DEP)
							PRDAngleScatter(line);
						else
							PRDScatter(line);
					}
				}
			}
		}
	}
	/* --- Whether the background scattering is iterated at each
		   wavelength (N-51) --                         -------------- */

	bgScatterSetup();

	/* --- Iterate the scattering in the background. Only iterate
		   when we are also planning to do main iterations -- --------- */

	if (rh->spectrum.updateJ) {
		sprintf(messageStr,
				" %s: Lambda iterating background scattering...\n\n",
				routineName);
		Error(MESSAGE, routineName, messageStr);

		niter = 0;
		while (rh->input.NmaxIter && niter < rh->input.NmaxScatter) {
			dJmax = solveSpectrum(FALSE, FALSE);
			if (dJmax < rh->input.iterLimit) break;
			niter++;
		}
	}
	getCPU(2, TIME_POLL, "Initial solution");
}
/* ------- end ---------------------------- initScatter.c ----------- */
