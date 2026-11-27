/* ------- file: -------------------------- writemetal_xdr.c --------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Apr 21 15:39:29 2009 --

	   --------------------------                      ----------RH-- */

/* --- Writes atomic data of background opacity package to output file.
	   XDR (external data representation) version. --  -------------- */

#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "xdr.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- writeMetals.c ----------- */

bool_t writeMetals(RHState *rh, char *fileName)
{
	const char   routineName[] = "writeMetals";
	register int n;

	int   result = TRUE, Nmetal;
	Atom *atom;
	FILE *fp_out;
	XDR   xdrs;

	if (!strcmp(fileName, "none")) return result;

	if ((fp_out = fopen(fileName, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open output file %s", fileName);
		return FALSE;
	}
	xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

	xdr_atom(rh, &xdrs, rh->atmos.H);

	/* --- Only write active atoms, write Hydrogen seperately -- ------ */

	Nmetal = rh->atmos.Natom - 1 - rh->atmos.Nactiveatom;
	result &= xdr_int(&xdrs, &Nmetal);
	for (n = 1; n < rh->atmos.Natom; n++) {
		atom = &rh->atmos.atoms[n];
		if (!atom->active)
			result &= xdr_atom(rh, &xdrs, atom);
	}

	xdr_destroy(&xdrs);
	fclose(fp_out);

	return result;
}
/* ------- end ---------------------------- writeMetals.c ----------- */
