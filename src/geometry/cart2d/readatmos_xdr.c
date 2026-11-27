/* ------- file: -------------------------- readatmos_xdr.c ---------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Jan 28 1997

	   --------------------------                      ----------RH-- */

/* --- Reads the atmospheric data and fills Atmos structure.
	   XDR (external data representation) version. --  -------------- */

/* --- Required input:

		 Nx, Nz      -- dimensions of grid            -- 2 4-byte integers
		 NHydr       -- number of Hydrogen levels     -- 4-byte integer
						(this includes the proton number)
		 hboundary   -- horizontal boundary condition -- 4-byte integer
						FIXED=0, PERIODIC=1
		 bvalue      -- vertical boundary values      -- 2 4-byte integers
						IRRADIATED=0, ZERO=1, THERMALIZED=2
						[TOPvalue, BOTTOMvalue]
		 dx          -- horizontal grid distances     -- double[Nx]
						(if hboundary == FIXED set the last value to 0.0)
		 z           -- vertical grid [km]            -- double[Nz]
		 T           -- temperature [K]               -- double[Nx*Nz]
		 ne          -- electron density [m^-3]       -- double[Nx*Nz]
		 vturb       -- microturbulence [km/s]        -- double[Nx*Nz]
		 vx          -- velocity field in x [km/s]    -- double[Nx*Nz]
		 vz          -- velocity field in z [km/s]    -- double[Nx*Nz]
		 nH          -- Hydrogen populations [m^-3]   -- double[Nx*Nz*NHydr]
						(see background.c in case of LTE hydrogen populations)


	   Conventions:

		 -- Indexing starts at the TOP of the atmosphere (in the vertical
			direction) and the left (horizontal direction).
		 -- z increases from the BOTTOM to the TOP.
		 -- Positive vz implies upward motion, positive vx implies
			rightward motion.


  See: geometry.c for the geometry of rays used by the code.
	   --                                              -------------- */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "constant.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

static XDR xdrs;

/* ------- begin -------------------------- readAtmos.c ------------- */

void readAtmos(RHState *rh, Atmosphere *atmos, Geometry *geometry)
{
	const char   routineName[] = "readAtmos";
	register int n, k, l;

	char       *filename;
	int         result, Nx, Nz, Nspace, bc[3];
	struct stat statBuffer;

	getCPU(rh, 2, TIME_START, NULL);

	/* --- Get abundances of background elements --      -------------- */

	readAbundance(rh, atmos);

	/* --- Open input file for model atmosphere --       -------------- */

	if ((atmos->fp_atmos = fopen(rh->input.atmos_input, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to open inputfile %s", rh->input.atmos_input);
	}
	/* --- Construct atmosID from filename and last modification date - */

	stat(rh->input.atmos_input, &statBuffer);
	if ((filename = strrchr(rh->input.atmos_input, '/')) != NULL)
		filename++;
	else
		filename = rh->input.atmos_input;
	sprintf(atmos->ID, "%s (%.24s)", filename,
			asctime(localtime(&statBuffer.st_mtime)));

	xdrstdio_create(&xdrs, atmos->fp_atmos, XDR_DECODE);

	/* --- Read size of horizontal and vertical grids, --  ------------ */

	result = xdr_int(&xdrs, &Nx);
	result &= xdr_int(&xdrs, &Nz);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read grid size from file %s", rh->input.atmos_input);
	}
	Nspace = Nx * Nz;

	/* --- Number of hydrogen levels to be used in background -- ------ */

	result = xdr_int(&xdrs, &atmos->NHydr);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read NHydr from file %s", rh->input.atmos_input);
	}
	if (atmos->NHydr < 2 && !atmos->H_LTE) {
		Error(rh, ERROR_LEVEL_2, routineName, "NHydr has to be at least 2, not %d to run with\n"
											  " NLTE hydrogen populations in background\n",
			  atmos->NHydr);
	}

	/* --- Keep duplicates of some of the geometrical dimensions in
		   Atmos structure --                            -------------- */

	atmos->Ndim = 2;
	atmos->N    = (int *) malloc(atmos->Ndim * sizeof(int));
	atmos->N[0] = geometry->Nx = Nx;
	atmos->N[1] = geometry->Nz = Nz;
	atmos->Nspace              = Nspace;

	/* --- Get the boundary conditions --                  ------------ */

	result = xdr_vector(&xdrs, (char *) bc, 3,
						sizeof(int), (xdrproc_t) rh_xdr_int);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read boundary conditions from file %s", rh->input.atmos_input);
	}
	geometry->hboundary = (enum boundary) bc[0];
	for (n = 0; n < 2; n++)
		geometry->bvalue[n] = (enum boundval) bc[1 + n];

	/* --- Check the validity of boundary conditions and values -- ---- */

	switch (geometry->hboundary) {
	case FIXED: break;
	case PERIODIC: break;
	default:
		Error(rh, ERROR_LEVEL_2, routineName, "Invalid horizontal boundary condition: %d",
			  geometry->hboundary);
		break;
	}

	switch (geometry->bvalue[TOP]) {
	case IRRADIATED: break;
	case ZERO: break;
	case THERMALIZED: break;
	default:
		Error(rh, ERROR_LEVEL_2, routineName, "Invalid boundary value at TOP: %d",
			  geometry->bvalue[TOP]);
		break;
	}
	switch (geometry->bvalue[BOTTOM]) {
	case IRRADIATED: break;
	case ZERO: break;
	case THERMALIZED: break;
	default:
		Error(rh, ERROR_LEVEL_2, routineName, "Invalid boundary value at BOTTOM: %d",
			  geometry->bvalue[BOTTOM]);
		break;
	}

	/* --- Get increments in x, store and check for monotonicity -- --- */

	geometry->dx = (double *) malloc(Nx * sizeof(double));
	result       = xdr_vector(&xdrs, (char *) geometry->dx, Nx,
							  sizeof(double), (xdrproc_t) rh_xdr_double);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read x grid from file %s", rh->input.atmos_input);
	}

	for (l = 0; l < ((geometry->hboundary == PERIODIC) ? Nx : Nx - 1); l++) {
		geometry->dx[l] *= KM_TO_M;
		if (geometry->dx[l] <= 0.0) {
			Error(rh, ERROR_LEVEL_2, routineName, "At l = %d: \n  x does not increase "
												  "strictly monotonically towards the right",
				  l);
		}
	}
	geometry->x    = (double *) malloc(Nx * sizeof(double));
	geometry->x[0] = 0.0;
	for (l = 0; l < Nx - 1; l++)
		geometry->x[l + 1] = geometry->x[l] + geometry->dx[l];

	/* --- Get vertical grid --                          -------------- */

	geometry->z = (double *) malloc(Nz * sizeof(double));
	result      = xdr_vector(&xdrs, (char *) geometry->z, Nz,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read z grid from file %s", rh->input.atmos_input);
	}
	for (k = 0; k < Nz; k++) geometry->z[k] *= KM_TO_M;

	geometry->dz = (double *) malloc(Nz * sizeof(double));
	for (k = 0; k < Nz - 1; k++) {
		geometry->dz[k] = geometry->z[k] - geometry->z[k + 1];
		if (geometry->dz[k] <= 0.0) {
			Error(rh, ERROR_LEVEL_2, routineName, "At k = %d:\n z does not decrease "
												  "strictly monotonically towards the bottom",
				  k);
		}
	}
	geometry->dz[Nz - 1] = 0.0;

	/* --- The atmospheric physics --                      ------------ */

	atmos->T     = (double *) malloc(Nspace * sizeof(double));
	atmos->ne    = (double *) malloc(Nspace * sizeof(double));
	atmos->vturb = (double *) malloc(Nspace * sizeof(double));
	geometry->vx = (double *) malloc(Nspace * sizeof(double));
	geometry->vz = (double *) malloc(Nspace * sizeof(double));

	/* --- Read quantities from atmosFile --               ------------ */

	result = xdr_vector(&xdrs, (char *) atmos->T, Nspace,
						sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) atmos->ne, Nspace,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) atmos->vturb, Nspace,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->vx, Nspace,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	result &= xdr_vector(&xdrs, (char *) geometry->vz, Nspace,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read atmospheric quantities from file %s", rh->input.atmos_input);
	}

	/* --- Convert input velocities from km/s to m/s --    ------------ */

	for (k = 0; k < Nspace; k++) {
		atmos->vturb[k] *= KM_TO_M;
		geometry->vx[k] *= KM_TO_M;
		geometry->vz[k] *= KM_TO_M;
	}
	/* --- Determine whether line profiles should account for shift due
		   to macroscopic velocities --                    ------------ */

	atmos->moving = FALSE;
	for (k = 0; k < Nspace; k++) {
		if ((fabs(geometry->vx[k]) > atmos->vmacro_tresh) ||
			(fabs(geometry->vz[k]) > atmos->vmacro_tresh)) {
			atmos->moving = TRUE;
			break;
		}
	}
	/* --- Finally, the hydrogen densities --            -------------- */

	atmos->nH = matrix_double(rh, atmos->NHydr, Nspace);
	result    = xdr_vector(&xdrs, (char *) atmos->nH[0], atmos->NHydr * Nspace,
						   sizeof(double), (xdrproc_t) rh_xdr_double);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read hydrogen density from file %s", rh->input.atmos_input);
	}
	atmos->nHtot = (double *) calloc(Nspace, sizeof(double));
	for (n = 0; n < atmos->NHydr; n++) {
		for (k = 0; k < Nspace; k++) {
			atmos->nHtot[k] += atmos->nH[n][k];
		}
	}
	/* --- Get angle-quadrature and copy geometry independent quantities
		   Nrays and wmu to atmos structure. --          -------------- */

	getAngleQuadr(rh, geometry);
	atmos->Nrays = geometry->Nrays;
	atmos->wmu   = geometry->wmu;

	/* --- Magnetic field is read here. --               -------------- */

	atmos->Stokes = readB(rh, atmos);

	getCPU(rh, 2, TIME_POLL, "Read Atmosphere");
}
/* ------- end ---------------------------- readAtmos.c ------------- */

/* ------- begin -------------------------- getBoundary.c ----------- */

void getBoundary(RHState *rh, Atmosphere *atmos, Geometry *geometry)
{
	const char routineName[] = "getBoundary";

	bool_t result = TRUE;
	int    Nx = geometry->Nx, Nz = geometry->Nz, Nspect = rh->spectrum.Nspect;

	/* --- Read boundary conditions and irradiation. First VERTICAL
		   boundary conditions --                        -------------- */

	switch (geometry->bvalue[TOP]) {
	case IRRADIATED:
		geometry->Itop = matrix_double(rh, Nspect, Nx);
		result &= xdr_vector(&xdrs, (char *) geometry->Itop[0], Nspect * Nx,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		break;
	default:
		geometry->Itop = NULL;
	}

	switch (geometry->bvalue[BOTTOM]) {
	case IRRADIATED:
		geometry->Ibottom = matrix_double(rh, Nspect, Nx);
		result &= xdr_vector(&xdrs, (char *) geometry->Ibottom[0], Nspect * Nx,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		break;
	default:
		geometry->Ibottom = NULL;
	}
	/* --- Get the HORIZONTAL boundary conditions --        ----------- */

	switch (geometry->hboundary) {
	case FIXED:
		geometry->Ileft  = matrix_double(rh, Nspect, Nz);
		geometry->Iright = matrix_double(rh, Nspect, Nz);

		result &= xdr_vector(&xdrs, (char *) geometry->Ileft[0], Nspect * Nz,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		result &= xdr_vector(&xdrs, (char *) geometry->Iright[0], Nspect * Nz,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		break;
	case PERIODIC:
		geometry->Ileft = NULL;
		break;
	}
	if (!result) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to read proper amount from input atmosphere");
	}

	xdr_destroy(&xdrs);
	fclose(atmos->fp_atmos);
}
/* ------- end ---------------------------- getBoundary.c ----------- */
