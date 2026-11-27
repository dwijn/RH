/* ------- file: -------------------------- radrate_xdr.c -----------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Feb 28 22:21:24 2024 --

	   --------------------------                      ----------RH-- */

/* --- Routines for writing radiative rates to file.
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

bool_t xdr_radrate(XDR *xdrs, Atom *atom);

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- writeRadRate.c - -------- */

bool_t writeRadRate(Atom *atom)
{
	const char routineName[] = "writeRadRate";

	char  ratesfile[15];
	FILE *fp;
	XDR   xdrs;

	/* --- Write radiatve rates for transitions treated in detail.-- -- */

	if (!strcmp(rh->input.radrateFile, "none") &&
		atom->Nprd == 0) return FALSE;

	sprintf(ratesfile, (atom->ID[1] == ' ') ? "radrate.%.1s.out" : "radrate.%.2s.out", atom->ID);

	if ((fp = fopen(ratesfile, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open output file %s", ratesfile);
		return FALSE;
	}
	xdrstdio_create(&xdrs, fp, XDR_ENCODE);

	if (!xdr_radrate(&xdrs, atom)) {
		Warning(rh, routineName, "Unable to write to output file %s", ratesfile);
	}

	xdr_destroy(&xdrs);
	fclose(fp);
	return TRUE;
}
/* ------- end ---------------------------- writeRadRate.c ---------- */

/* ------- begin -------------------------- readRadRate.c ----------- */

bool_t readRadRate(Atom *atom)
{
	const char routineName[] = "readRadRate";

	char  ratesfile[15];
	FILE *fp;
	XDR   xdrs;

	if (!strcmp(rh->input.radrateFile, "none") &&
		atom->Nprd == 0) return FALSE;

	sprintf(ratesfile, (atom->ID[1] == ' ') ? "radrate.%.1s.out" : "radrate.%.2s.out", atom->ID);

	/* --- Read radiative rates for transitions treated in detail. -- - */

	if ((fp = fopen(ratesfile, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to open input file %s", ratesfile);
		return FALSE;
	}
	xdrstdio_create(&xdrs, fp, XDR_DECODE);

	if (!xdr_radrate(&xdrs, atom)) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read from input file %s", ratesfile);
	}
	xdr_destroy(&xdrs);
	fclose(fp);
	return TRUE;
}
/* ------- end ---------------------------- readRadRate.c ----------- */

/* ------- begin -------------------------- xdr_radrate.c ----------- */

bool_t xdr_radrate(XDR *xdrs, Atom *atom)
{
	register int kr;

	bool_t           result = TRUE;
	AtomicLine      *line;
	AtomicContinuum *continuum;

	for (kr = 0; kr < atom->Nline; kr++) {
		line = &atom->line[kr];
		result &= xdr_vector(xdrs, (char *) line->Rij, rh->atmos.Nspace,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		result &= xdr_vector(xdrs, (char *) line->Rji, rh->atmos.Nspace,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
	}
	for (kr = 0; kr < atom->Ncont; kr++) {
		continuum = &atom->continuum[kr];
		result &= xdr_vector(xdrs, (char *) continuum->Rij, rh->atmos.Nspace,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		result &= xdr_vector(xdrs, (char *) continuum->Rji, rh->atmos.Nspace,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
	}
	return result;
}
/* ------- end ---------------------------- xdr_radrate.c ----------- */
