/* ------- file: -------------------------- solveray.c --------------

	   Version:       rh2.0, 3-D Cartesian, short characteristics
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Fri May  1 10:19:36 2020 --

	   --------------------------                      ----------RH-- */

/* --- Solves radiative transfer for given atmosphere and model atom
	   along a ray with arbitrary \mu_x and \mu_y, assuming the atom's
	   population numbers and angle-averaged radiation field are given.
	   The equation of transfer is solved by solveSpectrum, with the
	   short characteristics of the main program (Q-34; this said
	   Feautrier with Fourier interpolation in the horizontal planes).

	   Expects input file ``ray.input'' containing two lines of the form

		 mux  muy
		 Nspect  wave_index1  ....   wave_indexNspect
	   --                                              -------------- */

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "inputs.h"
#include "error.h"
#include "xdr.h"
#include "state.h"

#define COMMENT_CHAR  "#"
#define RAY_DOT_INPUT "ray.input"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

static Geometry geometry_data;
static RHState  rh_state = RH_STATE_INIT(THREE_D_PLANE, &geometry_data);
RHState *const  rh       = &rh_state;

char messageStr[MAX_LINE_SIZE];

/* ------- begin -------------------------- solveray.c -------------- */

int main(int argc, char *argv[])
{
	register int n, k;

	char   rayFileName[22], format[21], inputLine[MAX_LINE_SIZE];
	bool_t result, exit_on_EOF, crosscoupling, initialize, analyze_output,
		to_obs, equilibria_only;
	int Nspect, Nread, Nrequired, checkPoint,
		*wave_index = NULL;
	double     mux, muy, muz, sumsq, *S = NULL, *chi = NULL, *J = NULL;
	FILE      *fp_out, *fp_ray;
	XDR        xdrs;
	ActiveSet *as;

	setOptions(argc, argv);
	getCPU(0, TIME_START, NULL);
	SetFPEtraps();

	/* --- Read input data and initialize --             -------------- */

	readInput();
	rh->spectrum.updateJ = FALSE;

	/* --- Read input data for atmosphere --             -------------- */

	getCPU(1, TIME_START, NULL);
	readAtmos(&rh->atmos, rh->geometry);

	/* --- Read direction cosines for ray --             -------------- */

	if ((fp_ray = fopen(RAY_DOT_INPUT, "r")) == NULL) {
		sprintf(messageStr, "Unable to open inputfile %s", RAY_DOT_INPUT);
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}

	getLine(fp_ray, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
	Nread = sscanf(inputLine, "%lf %lf", &mux, &muy);
	checkNread(Nread, Nrequired = 2, argv[0], checkPoint = 1);

	if (mux <= -1.0 || mux >= 1.0) {
		sprintf(messageStr,
				"Value of mux = %f does not lie in interval <-1.0, 1.0>\n", mux);
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	if (muy <= -1.0 || muy >= 1.0) {
		sprintf(messageStr,
				"Value of muy = %f does not lie in interval <-1.0, 1.0>\n", muy);
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}

	/* --- The third cosine follows from the first two. Rays point
		   upward, and muz divides the grid spacing in the stencil
		   distances (rz = muz/dz in geometry.c), so muz must be positive:
		   no clamping to zero, a sum of squares at or above one (also by
		   rounding of the input) is an error (B-34). --  -------------- */

	sumsq = SQ(mux) + SQ(muy);
	if (sumsq >= 1.0) {
		sprintf(messageStr,
				"mux = %.12g and muy = %.12g leave no upward component: "
				"mux^2 + muy^2 = %.12g is not below 1\n",
				mux, muy, sumsq);
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	muz = sqrt(1.0 - sumsq);

	if ((rh->atmos.Stokes && rh->input.StokesMode == FIELD_FREE) ||
		rh->input.backgr_pol) {

		/* --- Want formal solution to be polarized in these cases -- --- */

		rh->input.StokesMode = FULL_STOKES;
	}
	/* --- Redefine geometry for just this one ray --    -------------- */

	rh->atmos.Nrays = rh->geometry->Nrays = 1;
	rh->geometry->mux[0]                  = mux;
	rh->geometry->muy[0]                  = muy;
	rh->geometry->muz[0]                  = muz;
	rh->geometry->wmu[0]                  = 1.0;

	if (rh->atmos.Stokes) Bproject();
	fillMesh(rh->geometry);

	rh->input.startJ = OLD_J;

	readAtomicModels();
	readMolecularModels();
	SortLambda();

	/* --- Open file with background opacities --        -------------- */

	if (rh->atmos.moving || rh->atmos.Stokes) {
		strcpy(rh->input.background_File, rh->input.background_ray_File);
		Background(analyze_output = FALSE, equilibria_only = FALSE);
	} else {
		Background(analyze_output = FALSE, equilibria_only = TRUE);

		if ((rh->atmos.fd_background =
				 open(rh->input.background_File, O_RDONLY, 0)) == -1) {
			sprintf(messageStr, "Unable to open inputfile %s",
					rh->input.background_File);
			Error(ERROR_LEVEL_2, argv[0], messageStr);
		}
		readBRS();
	}

	getProfiles();
	initSolution();
	initScatter();

	getCPU(1, TIME_POLL, "Total initialize");

	/* --- Solve radiative transfer equations --     ------------------ */

	solveSpectrum(FALSE, FALSE);

	/* --- Write emergent spectrum to output file --  ----------------- */

	strcpy(format, "spectrum_%4.2f_%4.2f");
	if (mux < 0.0) format[10] = '5';
	if (muy < 0.0) format[16] = '5';
	sprintf(rayFileName, format, mux, muy);
	if ((fp_out = fopen(rayFileName, "w")) == NULL) {
		sprintf(messageStr, "Unable to open output file %s", rayFileName);
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

	result = xdr_double(&xdrs, &mux);
	result &= xdr_double(&xdrs, &muy);
	result &= xdr_vector(&xdrs, (char *) rh->spectrum.I[0],
						 rh->spectrum.Nspect * rh->geometry->Nplane,
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	if (!result) {
		sprintf(messageStr, "Unable to write spectrum to file %s", rayFileName);
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}

	/* --- Read wavelength indices for which chi and S are to be
		   written out for the specified direction --    -------------- */

	Nread = fscanf(fp_ray, "%d", &Nspect);
	checkNread(Nread, 1, argv[0], checkPoint = 2);

	if (Nspect > 0) {
		wave_index = (int *) malloc(Nspect * sizeof(int));
		Nread      = 0;
		while (fscanf(fp_ray, "%d", &wave_index[Nread]) != EOF) Nread++;
		checkNread(Nread, Nspect, argv[0], checkPoint = 3);
		fclose(fp_ray);

		chi = (double *) malloc(rh->atmos.Nspace * sizeof(double));
		if (rh->atmos.Stokes)
			S = (double *) malloc(4 * rh->atmos.Nspace * sizeof(double));
		else
			S = (double *) malloc(rh->atmos.Nspace * sizeof(double));
	}
	result = xdr_int(&xdrs, &Nspect);
	if (!result) {
		sprintf(messageStr, "Unable to write number of wavelengths to file %s", rayFileName);
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}

	/* --- Go through the list of wavelengths --         -------------- */

	if (Nspect > 0 && rh->input.limit_memory)
		J = (double *) malloc(rh->atmos.Nspace * sizeof(double));

	for (n = 0; n < Nspect; n++) {
		if (wave_index[n] < 0 || wave_index[n] >= rh->spectrum.Nspect) {
			sprintf(messageStr, "Illegal value of wave_index[n]: %4d\n"
								"Value has to be between 0 and %4d\n",
					wave_index[n], rh->spectrum.Nspect);
			Error(ERROR_LEVEL_2, argv[0], messageStr);
			continue;
		}
		sprintf(messageStr, "Processing n = %4d, lambda = %9.3f [nm]\n",
				wave_index[n], rh->spectrum.lambda[wave_index[n]]);
		Error(MESSAGE, NULL, messageStr);

		as = &rh->spectrum.as[wave_index[n]];
		alloc_as(wave_index[n], crosscoupling = FALSE);
		Opacity(wave_index[n], 0, to_obs = TRUE, initialize = TRUE);
		readBackground(wave_index[n], 0, to_obs = TRUE);

		if (rh->input.limit_memory) {
			readJlambda(wave_index[n], J);
		} else
			J = rh->spectrum.J[wave_index[n]];

		/* --- Add the continuum opacity and emissivity -- -------------- */

		for (k = 0; k < rh->atmos.Nspace; k++) {
			chi[k] = as->chi[k] + as->chi_c[k];
			S[k]   = (as->eta[k] + as->eta_c[k] + as->sca_c[k] * J[k]) / chi[k];
		}
		result = xdr_int(&xdrs, &wave_index[n]);
		result &= xdr_vector(&xdrs, (char *) chi, rh->atmos.Nspace,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		result &= xdr_vector(&xdrs, (char *) S, rh->atmos.Nspace,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		if (!result) {
			sprintf(messageStr, "Unable to write chi and S to file %s", rayFileName);
			Error(ERROR_LEVEL_2, argv[0], messageStr);
		}

		free_as(wave_index[n], crosscoupling = FALSE);
	}
	/* --- If magnetic fields are present --             -------------- */

	if (rh->atmos.Stokes || rh->input.backgr_pol) {
		result = xdr_vector(&xdrs, (char *) rh->spectrum.Stokes_Q[0],
							rh->spectrum.Nspect * rh->geometry->Nplane, sizeof(double),
							(xdrproc_t) rh_xdr_double);
		result &= xdr_vector(&xdrs, (char *) rh->spectrum.Stokes_U[0],
							 rh->spectrum.Nspect * rh->geometry->Nplane, sizeof(double),
							 (xdrproc_t) rh_xdr_double);
		result &= xdr_vector(&xdrs, (char *) rh->spectrum.Stokes_V[0],
							 rh->spectrum.Nspect * rh->geometry->Nplane, sizeof(double),
							 (xdrproc_t) rh_xdr_double);
		if (!result) {
			sprintf(messageStr, "Unable to write Stokes spectra to file %s", rayFileName);
			Error(ERROR_LEVEL_2, argv[0], messageStr);
		}
	}

	if (Nspect > 0 && rh->input.limit_memory)
		free(J);

	xdr_destroy(&xdrs);
	fclose(fp_out);
	printTotalCPU();
}
/* ------- end ---------------------------- solveray.c -------------- */
