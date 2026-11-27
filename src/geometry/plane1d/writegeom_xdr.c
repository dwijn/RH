/* ------- file: -------------------------- writegeom_xdr.c ---------

	   Version:       rh2.0, 1-D plane-parallel
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Jan 26 17:34:29 2000 --

	   --------------------------                      ----------RH-- */

/* --- Writes geometry data for 1-D plane-parallel version to output.
	   XDR (external data representation) version. --  -------------- */

#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- writeGeometry.c --------- */

void writeGeometry(RHState *rh, Geometry *geometry)
{
	const char routineName[] = "writeGeometry";

	bool_t result      = TRUE;
	int    theTopology = (int) rh->topology;
	FILE  *fp_out;
	XDR    xdrs;

	if (!strcmp(rh->input.geometry_output, "none")) return;

	if ((fp_out = fopen(rh->input.geometry_output, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open output file %s",
			  rh->input.geometry_output);
		return;
	}
	xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

	result &= xdr_int(&xdrs, &theTopology);
	result &= xdr_int(&xdrs, &geometry->Nrays);
	result &= xdr_int(&xdrs, &geometry->Ndep);

	result &= xdr_vector(&xdrs, (char *) geometry->muz, geometry->Nrays,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->wmu, geometry->Nrays,
						 sizeof(double), (xdrproc_t) rh_xdr_double);

	result &= xdr_vector(&xdrs, (char *) geometry->height, geometry->Ndep,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->cmass, geometry->Ndep,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->tau_ref, geometry->Ndep,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->vel, geometry->Ndep,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to write geometry to file %s",
			  rh->input.geometry_output);
	}

	xdr_destroy(&xdrs);
	fclose(fp_out);
}
/* ------- end ---------------------------- writeGeometry.c --------- */
