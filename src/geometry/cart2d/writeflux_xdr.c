/* ------- file: -------------------------- writeflux_xdr.c ---------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Nov  7 11:39:47 2000 --

	   --------------------------                      ----------RH-- */

#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "constant.h"
#include "error.h"
#include "xdr.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- writeFlux.c ------------- */

bool_t writeFlux(RHState *rh, char *fileName)
{
	const char   routineName[] = "writeFlux";
	register int k, mu, nspect;

	bool_t  result = TRUE;
	int     Nrays = rh->geometry->Nrays, Nx = rh->geometry->Nx;
	double *flux, *wmuz;
	FILE   *fp_flux;
	XDR     xdrs;

	/* --- Write the radiative flux in the z-direction -- ------------ */

	if ((fp_flux = fopen(fileName, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open output file %s", fileName);
		return FALSE;
	}
	xdrstdio_create(&xdrs, fp_flux, XDR_ENCODE);

	flux = (double *) malloc(Nx * sizeof(double));
	wmuz = (double *) malloc(Nrays * sizeof(double));
	for (mu = 0; mu < Nrays; mu++)
		wmuz[mu] = rh->geometry->muz[mu] * rh->geometry->wmu[mu];

	for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++) {
		for (k = 0; k < Nx; k++) flux[k] = 0.0;
		for (mu = 0; mu < rh->geometry->Nrays; mu++) {
			for (k = 0; k < Nx; k++)
				flux[k] += rh->spectrum.I[nspect * Nrays + mu][k] * wmuz[mu];
		}
		for (k = 0; k < Nx; k++) flux[k] *= 2.0 * PI;
		result &= xdr_vector(&xdrs, (char *) flux, Nx,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
	}

	xdr_destroy(&xdrs);
	fclose(fp_flux);
	free(wmuz);
	free(flux);

	return result;
}
/* ------- end ---------------------------- writeFlux.c ------------- */
