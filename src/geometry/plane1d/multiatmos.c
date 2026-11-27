/* ------- file: -------------------------- multiatmos.c ------------

	   Version:       rh2.0, 1-D plane-parallel
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Mon Sep 22 16:45:38 2025 --

	   --------------------------                      ----------RH-- */

/* --- Reads atmospheric model in MULTI format. --     -------------- */

#include <ctype.h>
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
#include "error.h"
#include "inputs.h"
#include "statistics.h"
#include "xdr.h"
#include "state.h"

#define MULTI_COMMENT_CHAR "*"
#define N_HYDROGEN_MULTI   6

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- MULTIatmos.c ------------ */

void MULTIatmos(RHState *rh, Atmosphere *atmos, Geometry *geometry)
{
	const char   routineName[] = "MULTIatmos";
	register int k, n;

	char    scaleStr[21], inputLine[MAX_LINE_SIZE], *filename;
	bool_t  exit_on_EOF, enhanced_atmos_ID = FALSE;
	int     Nread, Ndep, Nrequired, checkPoint;
	double *dscale, turbpress, turbelecpress, nbaryon, meanweight;

	/* --- P/Pe total ionized gas. Edited: BRC --          ------------ */

	double threshold_ion = 1.9082806;
	double gas_pressure, elec_pressure;

	struct stat statBuffer;

	getCPU(rh, 2, TIME_START, NULL);

	/* --- Get abundances of background elements --        ------------ */

	readAbundance(rh, atmos);

	/* --- Open the input file for model atmosphere in MULTI format - - */

	if ((atmos->fp_atmos = fopen(rh->input.atmos_input, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to open inputfile %s", rh->input.atmos_input);
	} else {
		Message(rh, "\n -- reading input file: %s\n\n",
				rh->input.atmos_input);
	}

	atmos->NHydr = N_HYDROGEN_MULTI;

	/* --- Boundary condition at TOP of atmosphere --      ------------ */

	if (strcmp(rh->input.Itop, "none"))
		geometry->vboundary[TOP] = IRRADIATED;
	else
		geometry->vboundary[TOP] = ZERO;

	/* --- Boundary condition at BOTTOM of atmosphere --   ------------ */

	geometry->vboundary[BOTTOM] = THERMALIZED;

	/* --- Read atmos ID, scale type, gravity, and number of depth
		   points --                                       ------------ */

	getLine(rh, atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
	if (enhanced_atmos_ID) {

		/* --- Construct atmosID from filename and last modification date */

		stat(rh->input.atmos_input, &statBuffer);
		if ((filename = strrchr(rh->input.atmos_input, '/')) != NULL)
			filename++;
		else
			filename = rh->input.atmos_input;
		sprintf(atmos->ID, "%s (%.24s)", filename,
				asctime(localtime(&statBuffer.st_mtime)));
		Nread = 1;
	} else
		Nread = sscanf(inputLine, "%s", atmos->ID);

	getLine(rh, atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
	Nread += sscanf(inputLine, "%20s", scaleStr);
	getLine(rh, atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
	Nread += sscanf(inputLine, "%lf", &atmos->gravity);
	getLine(rh, atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
	Nread += sscanf(inputLine, "%d", &geometry->Ndep);
	checkNread(rh, Nread, Nrequired = 4, routineName, checkPoint = 1);

	/* --- Keep duplicates of some of the geometrical quantities in
		   Atmos structure --                            -------------- */

	atmos->Ndim   = 1;
	atmos->N      = (int *) malloc(atmos->Ndim * sizeof(int));
	atmos->Nspace = Ndep = geometry->Ndep;
	atmos->N[0]          = Ndep;

	atmos->gravity = POW10(atmos->gravity) * CM_TO_M;

	/* --- Allocate space for arrays that define structure -- --------- */

	geometry->tau_ref = (double *) malloc(Ndep * sizeof(double));
	geometry->cmass   = (double *) malloc(Ndep * sizeof(double));
	geometry->height  = (double *) malloc(Ndep * sizeof(double));
	atmos->T          = (double *) malloc(Ndep * sizeof(double));
	atmos->ne         = (double *) malloc(Ndep * sizeof(double));
	atmos->vturb      = (double *) malloc(Ndep * sizeof(double));
	geometry->vel     = (double *) malloc(Ndep * sizeof(double));

	dscale = (double *) malloc(Ndep * sizeof(double));
	for (k = 0; k < Ndep; k++) {
		getLine(rh, atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
		Nread = sscanf(inputLine, "%lf %lf %lf %lf %lf",
					   &dscale[k], &atmos->T[k], &atmos->ne[k],
					   &geometry->vel[k], &atmos->vturb[k]);
		checkNread(rh, Nread, Nrequired = 5, routineName, checkPoint = 2);
	}

	switch (toupper(scaleStr[0])) {
	case 'M':
		geometry->scale = COLUMN_MASS;
		for (k = 0; k < Ndep; k++)
			geometry->cmass[k] = POW10(dscale[k]) * (G_TO_KG / SQ(CM_TO_M));
		break;
	case 'T':
		geometry->scale = TAU500;
		for (k = 0; k < Ndep; k++) geometry->tau_ref[k] = POW10(dscale[k]);
		break;
	case 'H':
		geometry->scale = GEOMETRIC;
		for (k = 0; k < Ndep; k++) geometry->height[k] = dscale[k] * KM_TO_M;
		break;
	default:
		Error(rh, ERROR_LEVEL_2, routineName, "Unknown depth scale string in file %s: %s",
			  rh->input.atmos_input, scaleStr);
	}
	free(dscale);

	for (k = 0; k < Ndep; k++) {
		geometry->vel[k] *= KM_TO_M;
		atmos->vturb[k] *= KM_TO_M;
		atmos->ne[k] /= CUBE(CM_TO_M);
	}
	atmos->moving = FALSE;
	for (k = 0; k < Ndep; k++) {
		if (fabs(geometry->vel[k]) > atmos->vmacro_tresh) {
			atmos->moving = TRUE;
			break;
		}
	}
	/* --- Get angle-quadrature and copy geometry independent quantity
		   wmu to atmos structure. --                    -------------- */

	getAngleQuad(rh, geometry);
	atmos->wmu = geometry->wmu;

	/* --- Magnetic field is read here. --               -------------- */

	atmos->Stokes = readB(rh, atmos);

	/* --- Read Hydrogen populations if present --       -------------- */

	atmos->nH = matrix_double(rh, atmos->NHydr, Ndep);
	for (k = 0; k < Ndep; k++) {
		if (getLine(rh, atmos->fp_atmos, MULTI_COMMENT_CHAR, inputLine,
					exit_on_EOF = FALSE) == EOF) break;
		Nread = sscanf(inputLine, "%lf %lf %lf %lf %lf %lf",
					   &atmos->nH[0][k], &atmos->nH[1][k], &atmos->nH[2][k],
					   &atmos->nH[3][k], &atmos->nH[4][k], &atmos->nH[5][k]);
		checkNread(rh, Nread, Nrequired = 6, routineName, checkPoint = 3);
	}
	if (k > 0 && k < Ndep) {
		Error(rh, ERROR_LEVEL_2, routineName, "Reached end of input file %s before all data was read",
			  rh->input.atmos_input);
	} else if (k == 0) {

		/* --- No hydrogen populations supplied: use LTE populations
			   like MULTI does --                          -------------- */

		if (geometry->scale != COLUMN_MASS) {
			Error(rh, ERROR_LEVEL_2, routineName, "Height scale should be COLUMNMASS when nH not supplied: "
												  "File %s",
				  rh->input.atmos_input);
		}
		atmos->nHtot = (double *) calloc(Ndep, sizeof(double));

		atmos->H_LTE = TRUE;
		meanweight   = atmos->avgMolWght * AMU;
		for (k = 0; k < Ndep; k++) {
			turbpress     = 0.5 * meanweight * SQ(atmos->vturb[k]);
			turbelecpress = 0.5 * M_ELECTRON * SQ(atmos->vturb[k]);

			/* --- Edits to prevent negative pressures --    -------------- */

			gas_pressure  = atmos->gravity * geometry->cmass[k];     /* BRC */
			elec_pressure = atmos->ne[k] * KBOLTZMANN * atmos->T[k]; /* BRC */
			if (gas_pressure < threshold_ion * elec_pressure) {
				gas_pressure = threshold_ion * elec_pressure; /* BRC */
			}
			nbaryon = (gas_pressure -
					   atmos->ne[k] * (KBOLTZMANN * atmos->T[k] + turbelecpress)); /* BRC */

			atmos->nHtot[k] = nbaryon /
							  (atmos->totalAbund * (KBOLTZMANN * atmos->T[k] + turbpress));
		}
	} else if (k == Ndep) {
		atmos->nHtot = (double *) calloc(Ndep, sizeof(double));
		for (n = 0; n < atmos->NHydr; n++) {
			for (k = 0; k < Ndep; k++) {
				atmos->nH[n][k] /= CUBE(CM_TO_M);
				atmos->nHtot[k] += atmos->nH[n][k];
			}
		}
	}

	getCPU(rh, 2, TIME_POLL, "Read Atmosphere");
}
/* ------- end ---------------------------- MULTIatmos.c ------------ */

/* ------- begin -------------------------- convertScales.c --------- */

void convertScales(RHState *rh, Atmosphere *atmos, Geometry *geometry)
{
	register int k;

	int        ref_index, Ndep = geometry->Ndep;
	double    *rho, *height, *cmass, *tau_ref, h_zero;
	ActiveSet *as;

	/* --- Convert between different depth scales --       ------------ */

	height  = geometry->height;
	cmass   = geometry->cmass;
	tau_ref = geometry->tau_ref;

	rho = (double *) malloc(Ndep * sizeof(double));
	for (k = 0; k < Ndep; k++)
		rho[k] = (AMU * atmos->wght_per_H) * atmos->nHtot[k];

	/* --- Get opacity of reference wavelength --          ------------ */

	Locate(rh->spectrum.Nspect, rh->spectrum.lambda, atmos->lambda_ref, &ref_index);

	as = &rh->spectrum.as[ref_index];
	alloc_as(rh, ref_index, FALSE);
	readBackground(rh, ref_index, 0, 0);

	/* --- Convert to missing depth scales --              ------------ */

	switch (geometry->scale) {
	case COLUMN_MASS:
		height[0]  = 0.0;
		tau_ref[0] = as->chi_c[0] / rho[0] * cmass[0];
		for (k = 1; k < Ndep; k++) {
			height[k] = height[k - 1] - 2.0 * (cmass[k] - cmass[k - 1]) /
											(rho[k - 1] + rho[k]);
			tau_ref[k] = tau_ref[k - 1] + 0.5 * (as->chi_c[k - 1] + as->chi_c[k]) *
											  (height[k - 1] - height[k]);
		}
		break;
	case TAU500:
		height[0] = 0.0;
		cmass[0]  = (tau_ref[0] / as->chi_c[0]) * rho[0];
		for (k = 1; k < Ndep; k++) {
			height[k] = height[k - 1] - 2.0 * (tau_ref[k] - tau_ref[k - 1]) /
											(as->chi_c[k - 1] + as->chi_c[k]);
			cmass[k] = cmass[k - 1] + 0.5 * (rho[k - 1] + rho[k]) *
										  (height[k - 1] - height[k]);
		}
		break;
	case GEOMETRIC:
		cmass[0] = (atmos->nHtot[0] * atmos->totalAbund + atmos->ne[0]) *
				   (KBOLTZMANN * atmos->T[0] / atmos->gravity);
		tau_ref[0] = 0.5 * as->chi_c[0] * (height[0] - height[1]);
		if (tau_ref[0] > 1.0) tau_ref[0] = 0.0;
		for (k = 1; k < Ndep; k++) {
			cmass[k] = cmass[k - 1] + 0.5 * (rho[k - 1] + rho[k]) *
										  (height[k - 1] - height[k]);
			tau_ref[k] = tau_ref[k - 1] + 0.5 * (as->chi_c[k - 1] + as->chi_c[k]) *
											  (height[k - 1] - height[k]);
		}
		break;
	}
	free_as(rh, ref_index, FALSE);

	if ((geometry->scale == COLUMN_MASS) || (geometry->scale == TAU500)) {
		ExitOnError(rh, zeroPointAtUnityTau(rh, "convertScales", atmos->ID, Ndep, tau_ref,
											height, &h_zero));
		for (k = 0; k < Ndep; k++) height[k] = height[k] - h_zero;
	}

	free(rho);
}
/* ------- end ---------------------------- convertScales.c --------- */

/* ------- begin -------------------------- getBoundary.c ----------- */

void getBoundary(RHState *rh, Geometry *geometry)
{
	const char routineName[] = "getBoundary";

	bool_t result = TRUE;
	FILE  *fp_Itop;
	XDR    xdrs;

	switch (geometry->vboundary[TOP]) {
	case ZERO: break;
	case THERMALIZED: break;
	case IRRADIATED:

		Message(rh, "\n -- reading irradiance input file: %s\n\n",
				rh->input.Itop);

		geometry->Itop = matrix_double(rh, rh->spectrum.Nspect, geometry->Nrays);

		/* --- Open input file for irradiation at TOP --     -------------- */

		if ((fp_Itop = fopen(rh->input.Itop, "r")) == NULL) {
			Error(rh, ERROR_LEVEL_2, routineName, "Unable to open inputfile %s", rh->input.Itop);
		}
		xdrstdio_create(&xdrs, fp_Itop, XDR_DECODE);

		result &= xdr_vector(&xdrs, (char *) geometry->Itop[0],
							 rh->spectrum.Nspect * geometry->Nrays,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		if (!result) {
			Error(rh, ERROR_LEVEL_2, routineName, "Unable to read irradiation data at TOP of atmosphere");
		}
		xdr_destroy(&xdrs);
		fclose(fp_Itop);
		break;
	case REFLECTIVE:
		break;
	default:
		Error(rh, ERROR_LEVEL_2, routineName,
			  "Invalid boundary condition at the TOP of atmosphere");
	}

	switch (geometry->vboundary[BOTTOM]) {
	case ZERO: break;
	case THERMALIZED: break;
	case IRRADIATED:
		geometry->Ibottom = matrix_double(rh, rh->spectrum.Nspect, geometry->Nrays);

		/* --- Infalling intensities at BOTTOM should be read here -- --- */

		Error(rh, ERROR_LEVEL_1, routineName,
			  "Boundary condition IRRADIATED at BOTTOM not yet implemented");
		break;
	case REFLECTIVE:
		break;
	default:
		Error(rh, ERROR_LEVEL_2, routineName,
			  "Invalid boundary condition at the BOTTOM of atmosphere");
	}
}
/* ------- end ---------------------------- getBoundary.c ----------- */
