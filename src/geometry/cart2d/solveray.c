/* ------- file: -------------------------- solveray.c --------------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Fri May  1 10:20:24 2020 --

	   --------------------------                      ----------RH-- */

/* --- Solves radiative transfer for given atmosphere and model atom
	   along a ray with arbitrary \mu_x and \mu_z, assuming the atom's
	   population numbers and angle-averaged radiation field is given.
	   The equation of transfer is solved with the shortcharacteristics
	   method.

  See: P. B. Kunasz and L. H. Auer 1988, JQSRT 39, p. 67


	   Expects input file ``ray.input'' containing two lines of the form

		 mux  muz
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
#include "atom.h"
#include "background.h"
#include "statistics.h"
#include "inputs.h"
#include "error.h"
#include "xdr.h"
#include "state.h"

#define COMMENT_CHAR   "#"
#define RAY_INPUT_FILE "ray.input"
#define MU_SUMSQ_TOL   1.0E-9

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

static Geometry       geometry_data;
static RHState        rh_state = RH_STATE_INIT(TWO_D_PLANE, &geometry_data);
static RHState *const rh       = &rh_state;

/* ------- begin -------------------------- solveray.c -------------- */

int main(int argc, char *argv[])
{
	register int n, k;

	char   rayFileName[22], inputLine[MAX_LINE_SIZE];
	bool_t result, exit_on_EOF, to_obs, crosscoupling, initialize,
		analyze_output, equilibria_only;
	int        Nspect, Nread, Nrequired, checkPoint, *wave_index = NULL;
	double     mux, muy, muz, sumsq, *S = NULL, *chi = NULL, *J = NULL;
	ActiveSet *as;
	FILE      *fp_out, *fp_ray;
	XDR        xdrs;

	setOptions(rh, argc, argv);
	getCPU(rh, 0, TIME_START, NULL);
	SetFPEtraps();

	/* --- Read input data and initialize --             -------------- */

	readInput(rh);
	rh->spectrum.updateJ = FALSE;

	/* --- Read input data for atmosphere --             -------------- */

	getCPU(rh, 1, TIME_START, NULL);
	readAtmos(rh, &rh->atmos, rh->geometry);

	/* --- Read direction cosines for ray --             -------------- */

	if ((fp_ray = fopen(RAY_INPUT_FILE, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Unable to open inputfile %s", RAY_INPUT_FILE);
	}

	getLine(rh, fp_ray, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
	Nread = sscanf(inputLine, "%lf %lf", &mux, &muz);
	checkNread(rh, Nread, Nrequired = 2, argv[0], checkPoint = 1);

	if (mux <= -1.0 || mux >= 1.0) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Value of mux = %f does not lie in interval <-1.0, 1.0>\n", mux);
	}
	if (muz <= 0.0 || muz > 1.0) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Value of muz = %f does not lie in interval <0.0, 1.0]\n", muz);
	}

	/* --- The third cosine follows from the first two. A sum of squares
		   above one by rounding of the input (e.g. mux written to 10
		   decimals, which exceeds 1 by up to 3e-10) is clamped, muy = 0
		   (ray in the x-z plane) being a valid direction; a larger
		   excess is an inconsistent input.
		   A NaN here made all Stokes parameters at line wavelengths NaN
		   (B-34). --                                    -------------- */

	sumsq = SQ(mux) + SQ(muz);
	if (sumsq > 1.0 + MU_SUMSQ_TOL) {
		Error(rh, ERROR_LEVEL_2, argv[0], "mux = %.12g and muz = %.12g are inconsistent: "
										  "mux^2 + muz^2 = %.12g exceeds 1\n",
			  mux, muz, sumsq);
	}
	muy = sqrt(fmax(0.0, 1.0 - sumsq));

	if ((rh->atmos.Stokes && rh->input.StokesMode == FIELD_FREE) ||
		rh->input.backgr_pol) {

		/* --- Want formal solution to be polarized in these cases -- --- */

		rh->input.StokesMode = FULL_STOKES;
	}

	/* --- Redefine geometry for just this one ray and fill mesh -- --- */

	rh->atmos.Nrays = rh->geometry->Nrays = 1;
	rh->geometry->mux[0]                  = mux;
	rh->geometry->muy[0]                  = muy;
	rh->geometry->muz[0]                  = muz;
	rh->geometry->wmu[0]                  = 1.0;

	if (rh->atmos.Stokes) Bproject(rh);
	fillMesh(rh, rh->geometry);

	rh->input.startJ = OLD_J;

	readAtomicModels(rh);
	readMolecularModels(rh);
	SortLambda(rh);

	getBoundary(rh, &rh->atmos, rh->geometry);

	/* --- Open file with background opacities --        -------------- */

	if (rh->atmos.moving || rh->atmos.Stokes) {
		strcpy(rh->input.background_File, rh->input.background_ray_File);
		Background(rh, analyze_output = FALSE, equilibria_only = FALSE);
	} else {
		Background(rh, analyze_output = FALSE, equilibria_only = TRUE);

		if ((rh->atmos.fd_background =
				 open(rh->input.background_File, O_RDONLY, 0)) == -1) {
			Error(rh, ERROR_LEVEL_2, argv[0], "Unable to open inputfile %s",
				  rh->input.background_File);
		}
		readBRS(rh);
	}

	getProfiles(rh);
	initSolution(rh);
	initScatter(rh);
	getCPU(rh, 1, TIME_POLL, "Total initialize");

	/* --- Solve radiative transfer equations --     ------------------ */

	solveSpectrum(rh, FALSE, FALSE);

	/* --- Write emergent spectrum to output file --  ----------------- */

	sprintf(rayFileName,
			(mux >= 0.0) ? "spectrum_%4.2f_%4.2f" : "spectrum_%5.2f_%4.2f",
			mux, muz);
	if ((fp_out = fopen(rayFileName, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Unable to open output file %s", rayFileName);
	}
	xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

	result = xdr_double(&xdrs, &mux);
	result &= xdr_double(&xdrs, &muz);
	result &= xdr_vector(&xdrs, (char *) rh->spectrum.I[0],
						 rh->spectrum.Nspect * rh->atmos.N[0],
						 sizeof(double), (xdrproc_t) rh_xdr_double);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Unable to write spectrum to file %s", rayFileName);
	}

	/* --- Read wavelength indices for which chi and S are to be
		   written out for the specified direction --    -------------- */

	Nread = fscanf(fp_ray, "%d", &Nspect);
	checkNread(rh, Nread, 1, argv[0], checkPoint = 2);

	if (Nspect > 0) {
		wave_index = (int *) malloc(Nspect * sizeof(int));
		Nread      = 0;
		while (fscanf(fp_ray, "%d", &wave_index[Nread]) != EOF) Nread++;
		checkNread(rh, Nread, Nspect, argv[0], checkPoint = 3);
		fclose(fp_ray);

		chi = (double *) malloc(rh->atmos.Nspace * sizeof(double));
		if (rh->atmos.Stokes)
			S = (double *) malloc(4 * rh->atmos.Nspace * sizeof(double));
		else
			S = (double *) malloc(rh->atmos.Nspace * sizeof(double));
	}
	result = xdr_int(&xdrs, &Nspect);
	if (!result) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Unable to write number of wavelengths to file %s", rayFileName);
	}

	/* --- Go through the list of wavelengths --         -------------- */

	if (Nspect > 0 && rh->input.limit_memory)
		J = (double *) malloc(rh->atmos.Nspace * sizeof(double));

	for (n = 0; n < Nspect; n++) {
		if (wave_index[n] < 0 || wave_index[n] >= rh->spectrum.Nspect) {
			Error(rh, ERROR_LEVEL_2, argv[0], "Illegal value of wave_index[n]: %4d\n"
											  "Value has to be between 0 and %4d\n",
				  wave_index[n], rh->spectrum.Nspect);
			continue;
		}
		Message(rh, "Processing n = %4d, lambda = %9.3f [nm]\n",
				wave_index[n], rh->spectrum.lambda[wave_index[n]]);

		as = &rh->spectrum.as[wave_index[n]];
		alloc_as(rh, wave_index[n], crosscoupling = FALSE);
		Opacity(rh, wave_index[n], 0, to_obs = TRUE, initialize = TRUE);
		readBackground(rh, wave_index[n], 0, to_obs = TRUE);

		if (rh->input.limit_memory) {
			readJlambda(rh, wave_index[n], J);
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
			Error(rh, ERROR_LEVEL_2, argv[0], "Unable to write chi and S to file %s", rayFileName);
		}

		free_as(rh, wave_index[n], crosscoupling = FALSE);
	}
	/* --- If magnetic fields are present --             -------------- */

	if (rh->atmos.Stokes || rh->input.backgr_pol) {
		result = xdr_vector(&xdrs, (char *) rh->spectrum.Stokes_Q[0],
							rh->spectrum.Nspect * rh->atmos.N[0], sizeof(double),
							(xdrproc_t) rh_xdr_double);
		result &= xdr_vector(&xdrs, (char *) rh->spectrum.Stokes_U[0],
							 rh->spectrum.Nspect * rh->atmos.N[0], sizeof(double),
							 (xdrproc_t) rh_xdr_double);
		result &= xdr_vector(&xdrs, (char *) rh->spectrum.Stokes_V[0],
							 rh->spectrum.Nspect * rh->atmos.N[0], sizeof(double),
							 (xdrproc_t) rh_xdr_double);
		if (!result) {
			Error(rh, ERROR_LEVEL_2, argv[0], "Unable to write Stokes spectra to file %s", rayFileName);
		}
	}

	if (Nspect > 0 && rh->input.limit_memory)
		free(J);

	xdr_destroy(&xdrs);
	fclose(fp_out);
	printTotalCPU(rh);
}
/* ------- end ---------------------------- solveray.c -------------- */
