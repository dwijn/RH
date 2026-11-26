/* ------- file: -------------------------- writegeom_xdr.c ---------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Mar  2 15:19:00 2004 --

	   --------------------------                      ----------RH-- */

/* --- Writes geometry data for 2-D Cartesian version to output file.
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

void writeGeometry(Geometry *geometry)
{
	const char routineName[] = "writeGeometry";

	int result, Nrays = geometry->Nrays,
				Nspace = (geometry->Nx * geometry->Nz);
	FILE *fp_out;
	XDR   xdrs;

	if (!strcmp(rh->input.geometry_output, "none")) return;

	if ((fp_out = fopen(rh->input.geometry_output, "w")) == NULL) {
		Error(ERROR_LEVEL_1, routineName, "Unable to open output file %s",
			  rh->input.geometry_output);
		return;
	}
	xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

	result = xdr_enum(&xdrs, (enum_t *) &rh->topology);
	result &= xdr_int(&xdrs, &geometry->Nrays);
	result &= xdr_int(&xdrs, &geometry->Nx);
	result &= xdr_int(&xdrs, &geometry->Nz);
	result &= xdr_enum(&xdrs, (enum_t *) &rh->atmos.angleSet.set);

	result &= xdr_vector(&xdrs, (char *) geometry->mux, Nrays,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->muy, Nrays,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->wmu, Nrays,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->x, geometry->Nx,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->z, geometry->Nz,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->vx, Nspace,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->vz, Nspace,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	if (!result) {
		Error(ERROR_LEVEL_2, routineName, "Unable to write geometry to file %s", rh->input.geometry_output);
	}

	xdr_destroy(&xdrs);
	fclose(fp_out);
}
/* ------- end ---------------------------- writeGeometry.c --------- */
