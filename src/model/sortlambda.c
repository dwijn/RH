/* ------- file: -------------------------- sortlambda.c ------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Fri Dec  6 09:35:11 2019 --

	   Modifications by Chris Osborne (CMO), Glasgow University,
		to prevent interpolation problems at ACTIVE Bound-free edges.

	   --------------------------                      ----------RH-- */

/* --- Sorts wavelengths and determines what transitions are active at
	   which wavelength. --                            -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"
#include "xdr.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- SortLambda.c ------------ */

void SortLambda(RHState *rh)
{
	const char   routineName[] = "SortLambda";
	register int kr, n, m, nspect, la, nact;

	bool_t           unique, result;
	int              Nred, Nspectrum, Nlambda_original, Z, i = 0, j = 0, Nwave;
	double          *alpha_original, gbf_0, n_eff, *wavetable = NULL;
	PchipTable       bftable;
	ActiveSet       *as;
	Atom            *atom;
	Molecule        *molecule;
	AtomicLine      *line;
	AtomicContinuum *continuum;
	MolecularLine   *mrt;
	FILE            *fp_wavetable;
	XDR              xdrs;

	getCPU(rh, 2, TIME_START, NULL);

	/* --- First read the wavelength table if specified -- ------------ */

	result = TRUE;

	if (strcmp(rh->input.wavetable_input, "none")) {
		if ((fp_wavetable = fopen(rh->input.wavetable_input, "r")) == NULL) {
			Error(rh, ERROR_LEVEL_2, routineName, "Unable to open input file %s",
				  rh->input.wavetable_input);
		}
		xdrstdio_create(&xdrs, fp_wavetable, XDR_DECODE);

		result &= xdr_int(&xdrs, &Nwave);
		wavetable = (double *) malloc(Nwave * sizeof(double));
		result &= xdr_vector(&xdrs, (char *) wavetable, Nwave,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		if (!result) {
			Error(rh, ERROR_LEVEL_2, routineName, "Unable to read from input file %s",
				  rh->input.wavetable_input);
		}
		xdr_destroy(&xdrs);
		fclose(fp_wavetable);
	} else
		Nwave = 0;

	/* --- Add reference wavelength if necessary --      -------------- */

	Nspectrum = (rh->atmos.lambda_ref > 0.0) ? Nwave + 1 : Nwave;

	/* --- Symmetric lines with wavelengths of other transitions, of
		   the table or the reference wavelength in their blue half get
		   full profiles (N-55) --                       -------------- */

	if (rh->atmos.Natom > 0) {
		double *extra = (double *) malloc((Nwave + 1) * sizeof(double));

		for (kr = 0; kr < Nwave; kr++) extra[kr] = wavetable[kr];
		if (rh->atmos.lambda_ref > 0.0) extra[Nwave] = rh->atmos.lambda_ref;
		asymmetricBlends(rh, extra, Nspectrum);
		free(extra);
	}

	/* --- Go through the active atoms to collect number of wavelengths
		   for the lines and continua treated in detail -- ------------ */

	rh->atmos.Nactiveatom = 0;
	for (n = 0; n < rh->atmos.Natom; n++) {
		atom = &rh->atmos.atoms[n];
		if (atom->active) {
			for (kr = 0; kr < atom->Ncont; kr++) {
				Nspectrum += atom->continuum[kr].Nlambda;

				/* --- Add one more for edge wavelength (CMO) -- ------------ */

				Nspectrum++;
			}
			for (kr = 0; kr < atom->Nline; kr++)
				Nspectrum += atom->line[kr].Nlambda;

			atom->activeindex = rh->atmos.Nactiveatom;
			rh->atmos.Nactiveatom++;
		}
	}
	/* --- Store the pointers to the active atoms, so that they can be
		   enumerated --                                 -------------- */

	if (rh->atmos.Nactiveatom > 0) {
		rh->atmos.activeatoms = (Atom **) malloc(rh->atmos.Nactiveatom *
												 sizeof(Atom *));
		for (n = 0; n < rh->atmos.Natom; n++) {
			atom = &rh->atmos.atoms[n];
			if (atom->active)
				rh->atmos.activeatoms[atom->activeindex] = atom;
		}
	} else
		rh->atmos.activeatoms = NULL;

	/* --- Add number of wavelengths for lines in active molecules -- - */

	rh->atmos.Nactivemol = 0;
	for (n = 0; n < rh->atmos.Nmolecule; n++) {
		molecule = &rh->atmos.molecules[n];
		if (molecule->active) {
			for (kr = 0; kr < molecule->Nrt; kr++)
				Nspectrum += molecule->mrt[kr].Nlambda;

			molecule->activeindex = rh->atmos.Nactivemol;
			rh->atmos.Nactivemol++;
		}
	}
	/* --- Store the pointers to the active molecules, so that they can be
		  enumerated --                                  -------------- */

	if (rh->atmos.Nactivemol > 0) {
		rh->atmos.activemols = (Molecule **) malloc(rh->atmos.Nactivemol *
													sizeof(Molecule *));
		for (n = 0; n < rh->atmos.Nmolecule; n++) {
			molecule = &rh->atmos.molecules[n];
			if (molecule->active)
				rh->atmos.activemols[molecule->activeindex] = molecule;
		}
	} else
		rh->atmos.activemols = NULL;

	/* --- Fill the wavelength array --                  -------------- */

	nspect              = 0;
	rh->spectrum.lambda = (double *) malloc(Nspectrum * sizeof(double));

	/* --- First the referenece wavelength if specified -- ------------ */

	if (rh->atmos.lambda_ref > 0.0)
		rh->spectrum.lambda[nspect++] = rh->atmos.lambda_ref;

	/* --- Then the wavelength table --                  -------------- */

	for (kr = 0; kr < Nwave; kr++)
		rh->spectrum.lambda[nspect++] = wavetable[kr];

	/* --- Finally, all the detailed radiative transitions -- --------- */

	rh->atmos.NPRDactive = 0;

	for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
		atom = rh->atmos.activeatoms[nact];
		for (kr = 0; kr < atom->Ncont; kr++) {
			continuum = &atom->continuum[kr];

			/* --- NOTE(CMO): Modified the following to insert multiple copies
			   of lambdaEdge if we go over it --         -------------- */

			rh->spectrum.lambda[nspect++] = continuum->lambda0;
			for (la = 0; la < continuum->Nlambda; la++)
				if (continuum->lambda[la] <= continuum->lambda0)
					rh->spectrum.lambda[nspect++] = continuum->lambda[la];
				else
					rh->spectrum.lambda[nspect++] = continuum->lambda0;
		}
		for (kr = 0; kr < atom->Nline; kr++) {
			line = &atom->line[kr];
			for (la = 0; la < line->Nlambda; la++)
				rh->spectrum.lambda[nspect++] = line->lambda[la];

			if (line->PRD)
				rh->atmos.NPRDactive++;
		}
	}
	/* --- Active molecular lines --                     -------------- */

	for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
		molecule = rh->atmos.activemols[nact];
		for (kr = 0; kr < molecule->Nrt; kr++) {
			mrt = &molecule->mrt[kr];
			for (la = 0; la < mrt->Nlambda; la++)
				rh->spectrum.lambda[nspect++] = mrt->lambda[la];
		}
	}
	/* --- Sort the wavelengths in ascending order --    -------------- */

	qsort(rh->spectrum.lambda, Nspectrum, sizeof(double), qsascend);

	/* --- Check for duplicate wavelengths --            -------------- */

	rh->spectrum.Nspect = 1;
	for (nspect = 1; nspect < Nspectrum; nspect++) {
		if (rh->spectrum.lambda[nspect] > rh->spectrum.lambda[nspect - 1] + 1e-10) {
			rh->spectrum.lambda[rh->spectrum.Nspect] = rh->spectrum.lambda[nspect];
			rh->spectrum.Nspect++;
		}
	}
	Message(rh, "\n %s: Found %d unique wavelengths\n",
			routineName, rh->spectrum.Nspect);
	if (rh->spectrum.Nspect < Nspectrum) {
		Message(rh, " %s: Eliminated %d duplicate wavelengths\n\n",
				routineName, Nspectrum - rh->spectrum.Nspect);
	}
	/* --- Allocate space for wavelength array and active sets -- ----- */

	rh->spectrum.lambda = (double *) realloc(rh->spectrum.lambda,
											 rh->spectrum.Nspect * sizeof(double));
	rh->spectrum.as     = (ActiveSet *) calloc(rh->spectrum.Nspect, sizeof(ActiveSet));

	/* --- Go through each established wavelength and gather active
		   transitions --                                -------------- */

	for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++) {
		as = &rh->spectrum.as[nspect];

		/* --- as->art and as->mrt store the arrays of active
			   transitions for each active atom
			   and molecule at this wavelength seperately -- ------------ */

		if (rh->atmos.Nactiveatom > 0) {
			as->Nactiveatomrt =
				(int *) malloc(rh->atmos.Nactiveatom * sizeof(int));
			as->art = (AtomicTransition **)
				malloc(rh->atmos.Nactiveatom * sizeof(AtomicTransition *));

			for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
				as->Nactiveatomrt[nact] = 0;
				as->art[nact]           = (AtomicTransition *)
					malloc(N_MAX_OVERLAP * sizeof(AtomicTransition));
			}
		} else {
			as->Nactiveatomrt = NULL;
			as->art           = NULL;
		}

		if (rh->atmos.Nactivemol > 0) {
			as->Nactivemolrt = (int *)
				malloc(rh->atmos.Nactivemol * sizeof(int));
			as->mrt = (MolTransition **)
				malloc(rh->atmos.Nactivemol * sizeof(MolTransition *));

			for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
				as->Nactivemolrt[nact] = 0;
				as->mrt[nact]          = (MolTransition *)
					malloc(N_MAX_OVERLAP * sizeof(MolTransition));
			}
		} else {
			as->Nactivemolrt = NULL;
			as->mrt          = NULL;
		}
	}
	/* --- Determine what transitions are active at which
		   wavelengths and store the pointers to those transitions. - - */

	for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
		atom = rh->atmos.activeatoms[nact];

		/*--- Go through the continua first --             -------------- */

		Nred = 0;
		for (kr = 0; kr < atom->Ncont; kr++) {
			continuum = &atom->continuum[kr];

			/* --- Store the original wavelength array size -- ------------ */

			Nlambda_original = continuum->Nlambda;

			/* --- Find the indices of the lowest (Nblue) and highest (Nred)
				   wavelength of the current transition --     ------------ */

			Hunt(rh->spectrum.Nspect, rh->spectrum.lambda,
				 continuum->lambda[0], &continuum->Nblue);

			/* --- NOTE(CMO): Changed from lambda[Nlambda-1] -- ----------- */

			Hunt(rh->spectrum.Nspect, rh->spectrum.lambda,
				 continuum->lambda0, &Nred);
			continuum->Nlambda = Nred - continuum->Nblue + 1;

			/* --- Store the pointer to the current transition in the
			   active set (as) at each wavelength covered by the current
				   transition. Calculate wavelength integration weights - - */

			for (nspect = continuum->Nblue; nspect <= Nred; nspect++) {
				as   = &rh->spectrum.as[nspect];
				nact = atom->activeindex;

				as->art[nact][as->Nactiveatomrt[nact]].type            = ATOMIC_CONTINUUM;
				as->art[nact][as->Nactiveatomrt[nact]].ptype.continuum = continuum;
				as->Nactiveatomrt[nact]++;

				if (as->Nactiveatomrt[nact] == N_MAX_OVERLAP) {
					Error(rh, ERROR_LEVEL_2, routineName, "\n Too many overlapping transitions (> %d) "
														  "for atom %s and nspect = %d\n",
						  as->Nactiveatomrt[nact], atom->ID, nspect);
				}
			}
			/* --- In case of Bound-Free transition compute absorption
				   cross-section if wavelength dependence is hydrogenic,
			   interpolate if wavelength dependence is given
				   explicitly --                             -------------- */

			if (continuum->hydrogenic) {
				free(continuum->lambda);
				continuum->lambda = rh->spectrum.lambda + continuum->Nblue;
				continuum->alpha =
					(double *) realloc(continuum->alpha,
									   continuum->Nlambda * sizeof(double));

				Z     = atom->stage[continuum->j];
				n_eff = Z * sqrt(E_RYDBERG /
								 (atom->E[continuum->j] - atom->E[continuum->i]));
				gbf_0 = Gaunt_bf(continuum->lambda0, n_eff, Z);

				for (la = 0; la < continuum->Nlambda; la++) {
					continuum->alpha[la] = continuum->alpha0 *
										   Gaunt_bf(continuum->lambda[la], n_eff, Z) / gbf_0 *
										   CUBE(continuum->lambda[la] / continuum->lambda0);
				}
			} else {
				/* --- Monotone cubic Hermite in log lambda and log alpha
					   (pchip.c, N-32): positive, no new extrema. A spline
					   went negative around resonances --  -------------- */

				alpha_original = continuum->alpha;
				ExitOnError(rh, pchipTableInit(rh, &bftable, Nlambda_original, continuum->lambda,
											   alpha_original, PCHIP_LOGLOG));

				continuum->alpha =
					(double *) malloc(continuum->Nlambda * sizeof(double));
				pchipTableEval(&bftable, continuum->Nlambda,
							   rh->spectrum.lambda + continuum->Nblue, continuum->alpha);
				pchipTableFree(&bftable);

				free(continuum->lambda);
				continuum->lambda = rh->spectrum.lambda + continuum->Nblue;
				free(alpha_original);
			}
		}
		/* --- Then go through the bound-bound transitions -- ----------- */

		for (kr = 0; kr < atom->Nline; kr++) {
			line = &atom->line[kr];

			/* --- Store the original wavelength array size --  ----------- */

			Nlambda_original = line->Nlambda;

			/* --- Find the indices of the lowest (Nblue) and highest (Nred)
				   wavelength of the current transition --      ----------- */

			Hunt(rh->spectrum.Nspect, rh->spectrum.lambda, line->lambda[0],
				 &line->Nblue);
			Hunt(rh->spectrum.Nspect, rh->spectrum.lambda,
				 line->lambda[line->Nlambda - 1], &Nred);
			line->Nlambda = Nred - line->Nblue + 1;

			/* --- Store the pointer to the current transition in the
				   active set (as) at each wavelength covered by the current
				   transition. Calculate wavelength integration weights - - */

			for (nspect = line->Nblue; nspect <= Nred; nspect++) {
				as   = &rh->spectrum.as[nspect];
				nact = atom->activeindex;

				as->art[nact][as->Nactiveatomrt[nact]].type       = ATOMIC_LINE;
				as->art[nact][as->Nactiveatomrt[nact]].ptype.line = line;
				as->Nactiveatomrt[nact]++;

				if (as->Nactiveatomrt[nact] == N_MAX_OVERLAP) {
					Error(rh, ERROR_LEVEL_2, routineName, "\n Too many overlapping transitions (> %d) "
														  "for atom %s and nspect = %d\n",
						  as->Nactiveatomrt[nact], atom->ID, nspect);
				}
			}
			free(line->lambda);
			line->lambda = rh->spectrum.lambda + line->Nblue;
		}
	}

	for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
		molecule = rh->atmos.activemols[nact];

		Nred = 0;
		for (kr = 0; kr < molecule->Nrt; kr++) {
			mrt = &molecule->mrt[kr];

			/* --- Find the indices of the lowest (Nblue) and highest
				   wavelength of the current transition and allocate memory
				   for the wavelength integration weights -- -------------- */

			Hunt(rh->spectrum.Nspect, rh->spectrum.lambda,
				 mrt->lambda[0], &mrt->Nblue);
			Hunt(rh->spectrum.Nspect, rh->spectrum.lambda,
				 mrt->lambda[mrt->Nlambda - 1], &Nred);
			mrt->Nlambda = Nred - mrt->Nblue + 1;

			/* --- Repoint to proper position in wavelength array -- ------ */

			free(mrt->lambda);
			mrt->lambda = rh->spectrum.lambda + mrt->Nblue;

			/* --- Store the transition number of the current transition
				   in the active molecular set (as) at each wavelength
				   covered by the current transition. Calculate wavelength
				   integration weights --                      ------------ */

			for (nspect = mrt->Nblue;
				 nspect < mrt->Nblue + mrt->Nlambda; nspect++) {
				as   = &rh->spectrum.as[nspect];
				nact = molecule->activeindex;

				as->mrt[nact][as->Nactivemolrt[nact]].type         = mrt->type;
				as->mrt[nact][as->Nactivemolrt[nact]].ptype.vrline = mrt;
				as->Nactivemolrt[nact]++;

				if (as->Nactivemolrt[nact] == N_MAX_OVERLAP) {
					Error(rh, ERROR_LEVEL_2, routineName, "\n Too many overlapping transitions (> %d) "
														  "for molecule %s and nspect = %d\n",
						  as->Nactivemolrt[nact], molecule->ID, nspect);
				}
			}
		}
	}

	for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++) {
		as = &rh->spectrum.as[nspect];

		/* --- For each wavelength in the spectrum gather the unique set of
			   lower and upper levels involved in active atomic transitions at that
			   wavelength. These are needed in the calculation of the cross
			   coupling coefficients in the approximate lambda iteration -- */

		as->Nlower = (int *) malloc(rh->atmos.Nactiveatom * sizeof(int));
		as->Nupper = (int *) malloc(rh->atmos.Nactiveatom * sizeof(int));

		as->upper_levels =
			(int **) malloc(rh->atmos.Nactiveatom * sizeof(int *));
		as->lower_levels =
			(int **) malloc(rh->atmos.Nactiveatom * sizeof(int *));

		for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
			atom = rh->atmos.activeatoms[nact];

			/* --- First, for each wavelength in the spectrum gather
		   the unique set of lower levels of active transitions
		   within each atomic model --               -------------- */

			as->Nlower[nact] = 0;
			as->Nupper[nact] = 0;
			if (as->Nactiveatomrt[nact] > 0) {
				as->lower_levels[nact] =
					(int *) malloc(as->Nactiveatomrt[nact] * sizeof(int));
				as->upper_levels[nact] =
					(int *) malloc(as->Nactiveatomrt[nact] * sizeof(int));
			} else {
				as->lower_levels[nact] = NULL;
				as->upper_levels[nact] = NULL;
			}

			for (n = 0; n < as->Nactiveatomrt[nact]; n++) {
				unique = TRUE;
				switch (as->art[nact][n].type) {
				case ATOMIC_LINE:
					i = as->art[nact][n].ptype.line->i;
					j = as->art[nact][n].ptype.line->j;
					break;
				case ATOMIC_CONTINUUM:
					i = as->art[nact][n].ptype.continuum->i;
					j = as->art[nact][n].ptype.continuum->j;
					break;
				default:;
				}
				for (m = 0; m < as->Nlower[nact]; m++) {
					if (i == as->lower_levels[nact][m]) {
						unique = FALSE;
						break;
					}
				}
				if (unique) {
					as->lower_levels[nact][as->Nlower[nact]] = i;
					as->Nlower[nact]++;
				}

				/* --- Then add the upper level if unique --   -------------- */

				unique = TRUE;
				for (m = 0; m < as->Nupper[nact]; m++) {
					if (j == as->upper_levels[nact][m]) {
						unique = FALSE;
						break;
					}
				}
				if (unique) {
					as->upper_levels[nact][as->Nupper[nact]] = j;
					as->Nupper[nact]++;
				}
			}
		}
		/* --- Reallocate space for the atomic transition arrays -- ------- */

		for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
			if (as->Nactiveatomrt[nact] > 0) {
				as->art[nact] = (AtomicTransition *)
					realloc(as->art[nact],
							as->Nactiveatomrt[nact] * sizeof(AtomicTransition));

				as->lower_levels[nact] =
					(int *) realloc(as->lower_levels[nact],
									as->Nlower[nact] * sizeof(int));
				as->upper_levels[nact] =
					(int *) realloc(as->upper_levels[nact],
									as->Nupper[nact] * sizeof(int));
			}
		}
		/* --- Reallocate space for the molecular transition arrays -- -- */

		for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
			if (as->Nactivemolrt[nact] > 0) {
				as->mrt[nact] = (MolTransition *)
					realloc(as->mrt[nact],
							as->Nactivemolrt[nact] * sizeof(MolTransition));
			}
		}
	}

	/* --- Transitions that span the blue half of a symmetric line
		   integrate across it with the red half mirrored (N-62) -- */

	mirrorQuadrature(rh);

	getCPU(rh, 2, TIME_POLL, "SortLambda");
}
/* ------- end ---------------------------- SortLambda.c ------------ */
