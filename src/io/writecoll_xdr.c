/* ------- file: -------------------------- writecoll_xdr.c ---------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Fri Jul 26 13:38:58 2002 --

	   --------------------------                      ----------RH-- */

/* --- Routines for writing collisional rates to file.
	   XDR (external data representation) version. --  -------------- */

#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

bool_t xdr_collrate(RHState *rh, XDR *xdrs, Atom *atom);

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- writeCollisionRate.c ---- */

bool_t writeCollisionRate(RHState *rh, Atom *atom)
{
	const char routineName[] = "writeCollisionRate";

	char  ratesfile[16];
	FILE *fp_rate;
	XDR   xdrs;

	/* --- Write collision rates --                      -------------- */

	if (!strcmp(rh->input.collrateFile, "none")) return FALSE;

	sprintf(ratesfile, (atom->ID[1] == ' ') ? "collrate.%.1s.out" : "collrate.%.2s.out", atom->ID);

	if ((fp_rate = fopen(ratesfile, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open output file %s", ratesfile);
		return FALSE;
	}
	xdrstdio_create(&xdrs, fp_rate, XDR_ENCODE);

	if (!xdr_collrate(rh, &xdrs, atom)) {
		Warning(rh, routineName, "Unable to write to output file %s", ratesfile);
	}

	xdr_destroy(&xdrs);
	fclose(fp_rate);
	return TRUE;
}
/* ------- end ---------------------------- writeCollisionRate.c ---- */

/* ------- begin -------------------------- readCollisionRate.c ----- */

bool_t readCollisionRate(RHState *rh, Atom *atom)
{
	const char routineName[] = "readCollisionRate";

	char  ratesfile[16];
	FILE *fp_rate;
	XDR   xdrs;

	if (!strcmp(rh->input.collrateFile, "none")) return FALSE;

	sprintf(ratesfile, (atom->ID[1] == ' ') ? "collrate.%.1s.out" : "collrate.%.2s.out", atom->ID);

	/* --- Read collision rates --                       -------------- */

	if ((fp_rate = fopen(ratesfile, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open input file %s", ratesfile);
		return FALSE;
	}
	xdrstdio_create(&xdrs, fp_rate, XDR_DECODE);

	if (!xdr_collrate(rh, &xdrs, atom)) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read from input file %s", ratesfile);
	}
	xdr_destroy(&xdrs);
	fclose(fp_rate);
	return TRUE;
}
/* ------- end ---------------------------- readCollisionRate.c ----- */

/* ------- begin -------------------------- xdr_collrate.c ---------- */

bool_t xdr_collrate(RHState *rh, XDR *xdrs, Atom *atom)
{
	register int ij;

	bool_t result = TRUE;

	for (ij = 0; ij < SQ(atom->Nlevel); ij++) {
		result &= xdr_vector(xdrs, (char *) atom->C[ij], rh->atmos.Nspace,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
	}
	return result;
}
/* ------- end ---------------------------- xdr_collrate.c ---------- */
