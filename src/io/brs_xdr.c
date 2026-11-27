/* ------- file: -------------------------- brs_xdr.c ---------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Jan  3 14:35:46 2008 --

	   --------------------------                      ----------RH-- */

/* --- Routines to write and read background record structure to file.

	   XDR (external data representation) version. --  -------------- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "error.h"
#include "xdr.h"
#include "state.h"

#define BRS_DOT_OUT "brs.out"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- writeBRS.c -------------- */

void writeBRS(void)
{
	const char routineName[] = "writeBRS";

	FILE *fp_out;
	XDR   xdrs;

	if ((fp_out = fopen(BRS_DOT_OUT, "w")) == NULL) {
		Error(ERROR_LEVEL_1, routineName, "Unable to open output file %s", BRS_DOT_OUT);
		return;
	}
	xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

	if (!xdr_BRS(&xdrs)) {
		Error(ERROR_LEVEL_1, routineName, "Unable to write to output file %s", BRS_DOT_OUT);
	}

	xdr_destroy(&xdrs);
	fclose(fp_out);
}
/* ------- end ---------------------------- writeBRS.c -------------- */

/* ------- begin -------------------------- readBRS.c --------------- */

void readBRS(void)
{
	const char routineName[] = "readBRS";

	FILE *fp_in;
	XDR   xdrs;

	if ((fp_in = fopen(BRS_DOT_OUT, "r")) == NULL) {
		Error(ERROR_LEVEL_2, routineName, "Unable to open input file %s", BRS_DOT_OUT);
	}
	xdrstdio_create(&xdrs, fp_in, XDR_DECODE);

	if (!xdr_BRS(&xdrs)) {
		Error(ERROR_LEVEL_2, routineName, "Unable to read from input file %s", BRS_DOT_OUT);
	}

	xdr_destroy(&xdrs);
	fclose(fp_in);
}
/* ------- end ---------------------------- readBRS.c --------------- */

/* ------- begin -------------------------- xdr_BRS.c --------------- */

bool_t xdr_BRS(XDR *xdrs)
{
	const char   routineName[] = "xdr_BRS";
	register int nspect;

	char  *atmosID;
	bool_t result = TRUE, *hasline, *ispolarized;
	int    Nspace, Nsp, Nrecno;

	hasline =
		(bool_t *) malloc(rh->spectrum.Nspect * sizeof(bool_t));
	ispolarized =
		(bool_t *) malloc(rh->spectrum.Nspect * sizeof(bool_t));

	if (rh->atmos.moving || rh->atmos.Stokes)
		Nrecno = 2 * rh->spectrum.Nspect * rh->atmos.Nrays;
	else
		Nrecno = rh->spectrum.Nspect;

	if (xdrs->x_op == XDR_ENCODE) {
		atmosID = rh->atmos.ID;
		result &= xdr_counted_string(xdrs, &atmosID);
		result &= xdr_int(xdrs, &rh->atmos.Nspace);
		result &= xdr_int(xdrs, &rh->spectrum.Nspect);

		/* --- Flags for presence of background line --    -------------- */

		for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++)
			hasline[nspect] = rh->atmos.backgrflags[nspect].hasline;
		result &= xdr_vector(xdrs, (char *) hasline, rh->spectrum.Nspect,
							 sizeof(bool_t), (xdrproc_t) rh_xdr_bool);

		/* --- Flags for presence of polarized background -- ------------ */

		for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++)
			ispolarized[nspect] = rh->atmos.backgrflags[nspect].ispolarized;
		result &= xdr_vector(xdrs, (char *) ispolarized, rh->spectrum.Nspect,
							 sizeof(bool_t), (xdrproc_t) rh_xdr_bool);
	} else {
		free(rh->atmos.backgrflags);
		free(rh->atmos.backgrrecno);
		rh->atmos.backgrflags = (flags *) malloc(rh->spectrum.Nspect * sizeof(flags));
		rh->atmos.backgrrecno = (int *) malloc(Nrecno * sizeof(int));

		result &= xdr_counted_string(xdrs, &atmosID);
		if (!strstr(atmosID, rh->atmos.ID)) {
			Error(WARNING, routineName, "Input was derived from different atmosphere (%s) than current",
				  atmosID);
		}
		free(atmosID);

		result &= xdr_int(xdrs, &Nspace);
		result &= xdr_int(xdrs, &Nsp);
		if (Nspace != rh->atmos.Nspace || Nsp != rh->spectrum.Nspect) {
			free(hasline);
			return FALSE;
		}
		/* --- Flags for presence of background line --    -------------- */

		result &= xdr_vector(xdrs, (char *) hasline, rh->spectrum.Nspect,
							 sizeof(bool_t), (xdrproc_t) rh_xdr_bool);
		for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++)
			rh->atmos.backgrflags[nspect].hasline = hasline[nspect];

		/* --- Flags for presence of polarized background -- ------------ */

		result &= xdr_vector(xdrs, (char *) ispolarized, rh->spectrum.Nspect,
							 sizeof(bool_t), (xdrproc_t) rh_xdr_bool);
		for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++)
			rh->atmos.backgrflags[nspect].ispolarized = ispolarized[nspect];
	}

	result &= xdr_vector(xdrs, (char *) rh->atmos.backgrrecno, Nrecno,
						 sizeof(int), (xdrproc_t) rh_xdr_int);

	/* --- Bytes per value in the background file (4 since M-1) and
		   the power of two by which the emissivities are scaled in it.
		   Files written by older versions end before these and hold
		   unscaled doubles --                            ------------ */

	if (xdrs->x_op == XDR_ENCODE) {
		result &= xdr_int(xdrs, &rh->atmos.backgrbytes);
		result &= xdr_int(xdrs, &rh->atmos.backgretaexp);
	} else if (result) {
		if (!xdr_int(xdrs, &rh->atmos.backgrbytes)) {
			rh->atmos.backgrbytes  = sizeof(double);
			rh->atmos.backgretaexp = 0;
		} else
			result &= xdr_int(xdrs, &rh->atmos.backgretaexp);

		if (result && !((rh->atmos.backgrbytes == sizeof(float)) ||
						(rh->atmos.backgrbytes == sizeof(double) &&
						 rh->atmos.backgretaexp == 0))) {
			Error(ERROR_LEVEL_2, routineName, "Unknown format of the background file: %d bytes"
											  " per value, emissivity scale 2^%d",
				  rh->atmos.backgrbytes, rh->atmos.backgretaexp);
		}
	}

	free(hasline);
	free(ispolarized);

	return result;
}
/* ------- end ---------------------------- xdr_BRS.c --------------- */
