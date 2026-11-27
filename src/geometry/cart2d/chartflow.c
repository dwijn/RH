/* ------- file: -------------------------- chartflow.c -------------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Fri May 25 14:00:54 2018 --

	   --------------------------                      ----------RH-- */

/* --- Compute the x- and z components of the flux vector on the
	   whole grid for given wavelength indices nspect1  ..... nspectN.

	   The equation of transfer is solved with the shortcharacteristics
	   method.

  See: P. B. Kunasz and L. H. Auer 1988, JQSRT 39, p. 67


	   Expects input file ``flow.input'' containing one line of the form

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

#define COMMENT_CHAR   "#"
#define FLOW_DOT_INPUT "flow.input"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

static Geometry       geometry_data;
static RHState        rh_state = RH_STATE_INIT(TWO_D_PLANE, &geometry_data);
static RHState *const rh       = &rh_state;

/* ------- begin -------------------------- chartflow.c ------------- */

int main(int argc, char *argv[])
{
	register int k, mu, n;

	char   chartFileName[10], inputLine[MAX_LINE_SIZE];
	bool_t result, exit_on_EOF, crosscoupling, to_obs, initialize,
		boundbound, PRD_angle_dep, polarized_as,
		polarized_c, angle_dep, analyze_output, equilibria_only;
	int        Nspect, Nread, *wave_index, checkPoint;
	double     wmux, wmuz, *I, *J, *Psi, *S, *chi, *Fx, *Fz;
	FILE      *fp_out = NULL, *fp_flow = NULL;
	XDR        xdrs;
	ActiveSet *as;

	setOptions(rh, argc, argv);
	getCPU(rh, 0, TIME_START, NULL);
	SetFPEtraps(rh);

	/* --- Read input data and initialize --             -------------- */

	readInput(rh);
	rh->spectrum.updateJ = FALSE;

	getCPU(rh, 1, TIME_START, NULL);
	readAtmos(rh, &rh->atmos, rh->geometry);
	if (rh->atmos.Stokes) Bproject(rh);
	fillMesh(rh, rh->geometry);

	rh->input.startJ = OLD_J;

	readAtomicModels(rh);
	readMolecularModels(rh);
	SortLambda(rh);

	getBoundary(rh, &rh->atmos, rh->geometry);

	/* --- Solve chemical equilibrium and open background opacity file- */

	Background(rh, analyze_output = FALSE, equilibria_only = TRUE);

	if ((rh->atmos.fd_background =
			 open(rh->input.background_File, O_RDONLY, 0)) == -1) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Unable to open inputfile %s",
			  rh->input.background_File);
	}
	readBRS(rh);

	getProfiles(rh);
	initSolution(rh);
	initScatter(rh);

	/* --- Get wavelength indices for flow field --      -------------- */

	if ((fp_flow = fopen(FLOW_DOT_INPUT, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Unable to open inputfile %s", FLOW_DOT_INPUT);
	}

	getLine(rh, fp_flow, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
	Nread      = sscanf(strtok(inputLine, " "), "%d", &Nspect);
	wave_index = (int *) malloc(Nspect * sizeof(int));
	for (n = 0; n < Nspect; n++)
		Nread += sscanf(strtok(NULL, " "), "%d", wave_index + n);
	checkNread(rh, Nread, Nspect + 1, argv[0], checkPoint = 2);
	fclose(fp_flow);

	I   = (double *) malloc(rh->atmos.Nspace * sizeof(double));
	Psi = (double *) malloc(rh->atmos.Nspace * sizeof(double));
	chi = (double *) malloc(rh->atmos.Nspace * sizeof(double));
	S   = (double *) malloc(rh->atmos.Nspace * sizeof(double));

	Fx = (double *) malloc(rh->atmos.Nspace * sizeof(double));
	Fz = (double *) malloc(rh->atmos.Nspace * sizeof(double));

	/* --- Solve radiative transfer equation at specified wavelengths - */

	getCPU(rh, 1, TIME_POLL, "Total initialize");

	for (n = 0; n < Nspect; n++) {
		if (wave_index[n] < 0 || wave_index[n] >= rh->spectrum.Nspect) {
			Error(rh, ERROR_LEVEL_2, argv[0], "Illegal value of wave_index[n]: %4d\n"
											  "Value has to be between 0 and %4d\n",
				  wave_index[n], rh->spectrum.Nspect);
			continue;
		}
		Message(rh, "Processing n = %4d, lambda = %9.3f [nm]\n",
				wave_index[n], rh->spectrum.lambda[wave_index[n]]);

		/* --- Retrieve active set of transitions as and read current
			   mean intensities at wavelength nspect. --   -------------- */

		as = &rh->spectrum.as[wave_index[n]];
		alloc_as(rh, wave_index[n], crosscoupling = FALSE);

		if (rh->input.limit_memory) {
			J = (double *) malloc(rh->atmos.Nspace * sizeof(double));
			readJlambda(rh, wave_index[n], J);
		} else {
			J = rh->spectrum.J[wave_index[n]];
		}

		/* --- Check whether current active set includes a bound-bound
			   and/or polarized transition and/or angledependent PRD
			   transition. Otherwise, only angle-independent opacity and
			   source functions are needed --              -------------- */

		boundbound    = containsBoundBound(rh, as);
		PRD_angle_dep = (containsPRDline(rh, as) && rh->input.PRD_angle_dep);
		polarized_as  = containsPolarized(rh, as);
		polarized_c   = rh->atmos.backgrflags[wave_index[n]].ispolarized;
		angle_dep     = (polarized_as || polarized_c || PRD_angle_dep ||
                     (rh->atmos.moving &&
                      (boundbound || rh->atmos.backgrflags[wave_index[n]].hasline)));

		if (angle_dep) {
			for (mu = 0; mu < rh->geometry->Nrays; mu++) {
				wmux = rh->geometry->mux[mu] * rh->geometry->wmu[mu];
				wmuz = rh->geometry->muz[mu] * rh->geometry->wmu[mu];
				for (to_obs = 0; to_obs <= 1; to_obs++) {
					initialize = (mu == 0 && to_obs == 0);

					if (initialize || rh->atmos.backgrflags[wave_index[n]].hasline)
						readBackground(rh, wave_index[n], mu, to_obs);

					if (initialize || boundbound)
						Opacity(rh, wave_index[n], mu, to_obs, initialize);

					for (k = 0; k < rh->atmos.Nspace; k++) {
						chi[k] = as->chi[k] + as->chi_c[k];
						S[k]   = (as->eta[k] + as->eta_c[k] +
                                as->sca_c[k] * J[k]) /
							   chi[k];
					}
					Piecewise_Bezier3_2D(rh, rh->geometry, wave_index[n], mu, to_obs,
										 chi, S, I, Psi);

					if (to_obs) {
						for (k = 0; k < rh->atmos.Nspace; k++) {
							Fx[k] += I[k] * wmux;
							Fz[k] += I[k] * wmuz;
						}
					} else {
						for (k = 0; k < rh->atmos.Nspace; k++) {
							Fx[k] -= I[k] * wmux;
							Fz[k] -= I[k] * wmuz;
						}
					}
				}
			}
		} else {
			Opacity(rh, wave_index[n], 0, 0, initialize = TRUE);
			readBackground(rh, wave_index[n], 0, 0);

			/* --- The active opacity and emissivity are in as->chi and
				   as->eta, as in the angle-dependent case above. This
				   used chi_as and eta_as, which were allocated and never
				   filled (B-37) --                       -------------- */

			for (k = 0; k < rh->atmos.Nspace; k++) {
				chi[k] = as->chi[k] + as->chi_c[k];
				S[k]   = (as->eta[k] + as->eta_c[k] +
                        as->sca_c[k] * J[k]) /
					   chi[k];
			}
			for (mu = 0; mu < rh->geometry->Nrays; mu++) {
				wmux = rh->geometry->mux[mu] * rh->geometry->wmu[mu];
				wmuz = rh->geometry->muz[mu] * rh->geometry->wmu[mu];
				for (to_obs = 0; to_obs <= 1; to_obs++) {
					Piecewise_Bezier3_2D(rh, rh->geometry, wave_index[n], mu,
										 to_obs, chi, S, I, Psi);

					if (to_obs) {
						for (k = 0; k < rh->atmos.Nspace; k++) {
							Fx[k] += I[k] * wmux;
							Fz[k] += I[k] * wmuz;
						}
					} else {
						for (k = 0; k < rh->atmos.Nspace; k++) {
							Fx[k] -= I[k] * wmux;
							Fz[k] -= I[k] * wmuz;
						}
					}
				}
			}
		}
		free_as(rh, wave_index[n], crosscoupling = FALSE);

		/* --- Write the radiation flow to output file --  -------------- */

		sprintf(chartFileName, "flow_%4.4d", wave_index[n]);
		if ((fp_out = fopen(chartFileName, "w")) == NULL) {
			Error(rh, ERROR_LEVEL_2, argv[0], "Unable to open output file %s", chartFileName);
		}
		xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

		result = xdr_int(&xdrs, &wave_index[n]);
		result &= xdr_vector(&xdrs, (char *) Fx, rh->atmos.Nspace,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		result &= xdr_vector(&xdrs, (char *) Fz, rh->atmos.Nspace,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		if (!result) {
			Error(rh, ERROR_LEVEL_2, argv[0], "Unable to write flow data to file %s", chartFileName);
		}

		if (rh->input.limit_memory) free(J);
	}
	xdr_destroy(&xdrs);
	fclose(fp_out);
	printTotalCPU(rh);
}
/* ------- end ---------------------------- chartflow.c ------------- */
