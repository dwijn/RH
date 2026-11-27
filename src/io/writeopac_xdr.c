/* ------- file: -------------------------- writeopac_xdr.c ---------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Fri May 15 03:09:38 2009 --

	   --------------------------                      ----------RH-- */

/* --- Writes opacity and emissivity of active set of transitions for
	   each wavelength in the spectrum to file. If atmos.moving == TRUE
	   only write the quantities for the upward ray to file.

	   Write the record structure of the output file to a separate file
	   (ASRS_DOT_OUT). The record structure is an integer array of size
	   spectrum.Nspect in a static atmosphere and spectrum.Nspect *
	   atmos.Nrays otherwise. Whenever a wavelength has no overlap with
	   an active transition the record number for that wavelength is set
	   to -1;

	   XDR (external data representation) version. --  -------------- */

#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "state.h"

#define ASRS_DOT_OUT "asrs.out"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- writeOpacity.c ---------- */

void writeOpacity(RHState *rh)
{
	const char   routineName[] = "writeOpacity";
	register int nspect, mu;

	bool_t to_obs, initialize, crosscoupling, boundbound, polarized,
		PRD_angle_dep, result = TRUE;
	int        Nspace = rh->atmos.Nspace, record, *as_rn, Nrecord;
	FILE      *fp_out;
	XDR        xdrs;
	ActiveSet *as;

	if (!strcmp(rh->input.opac_output, "none")) return;

	if ((fp_out = fopen(rh->input.opac_output, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open output file %s",
			  rh->input.opac_output);
		return;
	}
	xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

	if (rh->atmos.moving || rh->atmos.Stokes ||
		(rh->atmos.NPRDactive > 0 && rh->input.PRD_angle_dep))
		Nrecord = rh->atmos.Nrays * rh->spectrum.Nspect;
	else
		Nrecord = rh->spectrum.Nspect;
	as_rn = (int *) malloc(Nrecord * sizeof(int));

	record = -1;
	for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++) {
		as = &rh->spectrum.as[nspect];

		if (containsActive(rh, as)) {
			alloc_as(rh, nspect, crosscoupling = FALSE);

			/* --- Check whether current active set includes a bound-bound
				   and/or polarized transition and/or angledependent PRD
				   transition. Otherwise, only angle-independent opacity and
				   source functions are needed --            -------------- */

			boundbound    = containsBoundBound(rh, as);
			PRD_angle_dep = (containsPRDline(rh, as) && rh->input.PRD_angle_dep);
			polarized     = containsPolarized(rh, as);

			/* --- Case of angle-dependent opacity and source function -- - */

			if (polarized || PRD_angle_dep || (rh->atmos.moving && boundbound)) {
				for (mu = 0; mu < rh->atmos.Nrays; mu++) {
					initialize = (mu == 0);
					Opacity(rh, nspect, mu, to_obs = TRUE, initialize);
					result &= xdr_vector(&xdrs, (char *) as->chi, Nspace,
										 sizeof(double), (xdrproc_t) rh_xdr_double);
					result &= xdr_vector(&xdrs, (char *) as->eta, Nspace,
										 sizeof(double), (xdrproc_t) rh_xdr_double);
					as_rn[nspect * rh->atmos.Nrays + mu] = ++record;
				}
			} else {
				Opacity(rh, nspect, 0, to_obs = TRUE, initialize = TRUE);
				result &= xdr_vector(&xdrs, (char *) as->chi, Nspace,
									 sizeof(double), (xdrproc_t) rh_xdr_double);
				result &= xdr_vector(&xdrs, (char *) as->eta, Nspace,
									 sizeof(double), (xdrproc_t) rh_xdr_double);
				if (rh->atmos.moving || rh->atmos.Stokes ||
					(rh->atmos.NPRDactive > 0 && rh->input.PRD_angle_dep)) {
					record++;
					for (mu = 0; mu < rh->atmos.Nrays; mu++)
						as_rn[nspect * rh->atmos.Nrays + mu] = record;
				} else
					as_rn[nspect] = ++record;
			}
			free_as(rh, nspect, crosscoupling = FALSE);
		} else {
			if (rh->atmos.moving || rh->atmos.Stokes ||
				(rh->atmos.NPRDactive > 0 && rh->input.PRD_angle_dep))
				for (mu = 0; mu < rh->atmos.Nrays; mu++)
					as_rn[nspect * rh->atmos.Nrays + mu] = -1;
			else
				as_rn[nspect] = -1;
		}
	}
	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to write opacities to file %s",
			  rh->input.opac_output);
	}
	xdr_destroy(&xdrs);
	fclose(fp_out);

	if ((fp_out = fopen(ASRS_DOT_OUT, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open output file %s", ASRS_DOT_OUT);
		return;
	}
	xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

	result &= xdr_vector(&xdrs, (char *) as_rn, Nrecord,
						 sizeof(int), (xdrproc_t) rh_xdr_int);
	free(as_rn);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to write record numbers to file %s",
			  ASRS_DOT_OUT);
	}
	xdr_destroy(&xdrs);
	fclose(fp_out);
}
/* ------- end ---------------------------- writeOpacity.c ---------- */
