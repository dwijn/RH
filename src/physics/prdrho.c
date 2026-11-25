/* ------- file: -------------------------- prdrho.c ----------------

	   Check of rho of the PRD lines after a sub-iteration (N-4).
	   --                                              -------------- */

#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "error.h"
#include "negtau.h"
#include "state.h"

/* --- Global variables --                             -------------- */

extern char messageStr[];

/* ------- begin -------------------------- PRDRhoSize.c ----------- */

long PRDRhoSize(const AtomicLine *line)
{
	/* --- The number of values of rho_prd[0][.] of a PRD line: a row
		   per wavelength, or per wavelength, ray and direction with
		   PRD_ANGLE_DEP, times the locations -- ---------------------- */

	long Nlamu = (rh->input.PRD_angle_dep == PRD_ANGLE_DEP)
					 ? 2L * rh->atmos.Nrays * line->Nlambda
					 : (long) line->Nlambda;

	return Nlamu * rh->atmos.Nspace;
}
/* ------- end ---------------------------- PRDRhoSize.c ----------- */

/* ------- begin -------------------------- PRDRhoCheck.c ----------- */

bool_t PRDRhoCheck(const Atom *atom, int kr)
{
	const char        routineName[] = "PRDRhoCheck";
	const AtomicLine *line          = atom->line + kr;
	const long        Ns            = rh->atmos.Nspace;
	long              n, N, Nbad = 0, nbad = -1;
	int               la, mu, to_obs;
	char              where[96], what[440];

	/* --- After a PRD sub-iteration (N-4): rho, the ratio of the PRD to
		   the CRD profile, must be positive and finite in every value:
		   the line emissivity is proportional to it. Returns TRUE if it
		   is. Otherwise a trial evaluation of G is invalidated
		   (negTauInvalid, as for a negative total opacity) and FALSE is
		   returned; outside a trial the run stops. Nothing is
		   clipped --                                   -------------- */

	N = PRDRhoSize(line);
	for (n = 0; n < N; n++) {
		double rho = line->rho_prd[0][n];

		if (!(rho > 0.0) || !isfinite(rho)) {
			if (Nbad++ == 0) nbad = n;
		}
	}
	if (Nbad == 0) return TRUE;

	la = (int) (nbad / Ns);
	if (rh->input.PRD_angle_dep == PRD_ANGLE_DEP) {
		to_obs = la % 2;
		mu     = (la / 2) % rh->atmos.Nrays;
		la     = la / (2 * rh->atmos.Nrays);
		snprintf(where, sizeof(where), "wavelength index %d (%.6f nm), ray %d%s", la,
				 line->lambda[la], mu, (to_obs) ? " upward" : " downward");
	} else
		snprintf(where, sizeof(where), "wavelength index %d (%.6f nm)", la,
				 line->lambda[la]);
	snprintf(what, sizeof(what),
			 "PRD line %d of atom %s (%d -> %d): rho = %.4E at %s, location %ld"
			 " (%ld of %ld values of rho are not positive and finite); a"
			 " non-positive rho gives a negative line emissivity",
			 kr, atom->ID, line->j, line->i, line->rho_prd[0][nbad], where, nbad % Ns,
			 Nbad, N);

	if (negTauInTrial()) {
		negTauInvalid(what);
		return FALSE;
	}
	snprintf(messageStr, MAX_MESSAGE_LENGTH, "\n  %s", what);
	Error(ERROR_LEVEL_2, routineName, messageStr);
	return FALSE;
}
/* ------- end ---------------------------- PRDRhoCheck.c ----------- */
