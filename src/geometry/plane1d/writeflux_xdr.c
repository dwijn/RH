/* ------- file: -------------------------- writeflux_xdr.c ---------

	   Version:       rh2.0, 1-D plane-parallel
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Jan 20 09:29:21 2004 --

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

bool_t writeFlux(char *flux_output)
{
	const char   routineName[] = "writeFlux";
	register int mu, nspect;

	bool_t  result = TRUE;
	double *flux, *wmuz;
	FILE   *fp_flux;
	XDR     xdrs;

	/* --- Write the radiative flux in the z-direction -- ------------ */

	if ((fp_flux = fopen(flux_output, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open outputfile %s", flux_output);
	}
	xdrstdio_create(&xdrs, fp_flux, XDR_ENCODE);

	wmuz = (double *) malloc(rh->atmos.Nrays * sizeof(double));
	for (mu = 0; mu < rh->atmos.Nrays; mu++)
		wmuz[mu] = rh->geometry->muz[mu] * rh->geometry->wmu[mu];

	flux = (double *) calloc(rh->spectrum.Nspect, sizeof(double));
	for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++) {
		for (mu = 0; mu < rh->atmos.Nrays; mu++)
			flux[nspect] += rh->spectrum.I[nspect][mu] * wmuz[mu];
		flux[nspect] *= 2.0 * PI;
	}
	result &= xdr_vector(&xdrs, (char *) flux, rh->spectrum.Nspect,
						 sizeof(double), (xdrproc_t) rh_xdr_double);

	xdr_destroy(&xdrs);
	fclose(fp_flux);
	free(wmuz);
	free(flux);

	return result;
}
/* ------- end ---------------------------- writeFlux.c ------------- */
