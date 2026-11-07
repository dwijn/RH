/* ------- file: -------------------------- background.c ------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu May 20 13:55:29 2021 --

	   --------------------------                      ----------RH-- */

/* Driving subroutine for background opacity sources.

 * Included at the moment:

  ++ Thomson scattering by free electrons

  ++ Hydrogen:
	-- Bound-free absorption and emission by neutral Hydrogen
	-- Free-free absorption and emission by neutral Hydrogen
	-- Rayleigh scattering by neutral Hydrogen and Helium
	-- Rayleigh scattering by molecular Hydrogen (H2)
	-- Bound-free absorption and emission by Hminus (H^-)
	-- Free-free absorption and emission by Hminus (H + e)
	-- Free-free absorption and emission by H2minus (H2 + e)
	-- Free-free absorption and emission by H2plus (H + H^+)

  ++ Metals:
	-- Bound-free absorption and emission from metals specified in
	   file background.input.
	-- Bound-bound absorption and emission from metals specified in
	   file background.input.
	-- LTE Bound-bound absorption and emission all elements from
	   a Kurucz line list.

  ++ Molecules:
	-- Chemical equilibrium is calculated for molecules specified
	   in file background.input and populations of constituent atoms
	   are reduced accordingly.
	-- molecular opacities (LTE) may be taken into account by specifying
	   data files with transition lists in the molecular input files.
	-- Bound-free absorption and emission by OH and CH molecules.

 * Atomic models are specified in atoms.input, molecules in
   molecules.input

 * Entries for the atoms.input and molecules.input files should have
   the form, respectively:

	-------------------------------------------------------------------
   |                                                                  |
   |   Nmetal                                                         |
   |                                                                  |
   |   model file ACTIVE/PASSIVE  INITIAL_SOLUTION   population file  |
   |                          .                                       |
   |                          .                                       |
   |                                                                  |
   |   Nmolecule                                                      |
   |                                                                  |
   |   model file ACTIVE/PASSIVE  INITIAL_SOLUTION                    |
   |                 .                                                |
   |                 .                                                |
	-------------------------------------------------------------------


   Nmetal and Nmolecule are the number of metal and molecule entries.
   metalID is the two-character atomID, the next entry is either
   set to LTE or NLTE, model_file is the input file containing atomic
   data for this metal (generic atomic input data format), and
   population_file is the input file containing the NLTE population
   numbers from a previous calculation. This last entry is only read when
   the second entry is set to NLTE.

   -- Units:
	  Wavelengths are given in nm, densities in m^-3, opacities in m^2,
	  and emissivities in J s^-1 Hz^-1 sr^-1.

 Note: The model atom file for hydrogen is specified in keyword.input.
	   If H_LTE = TRUE is specified there LTE hydrogen populations are
	   used. See: distribute_nH in the file hydrogen.c

 Note: Scattering opacity is added to total opacity after all
	   contributions have been computed.

 Note: The quantities chi_ai and eta_ai store the angle-inpendent
	   opacities and emissivities in case atmos.moving == TRUE.
	   In static atmospheres these quantities are just mapped to
	   atmos.chi_c and atmos.eta_c to save memory space.

 Note: If write_analyze_output == FALSE the auxiliary output files for
	   the Background Record Structure (BRS), metals, and molecules
	   are NOT written. This option is used when Background is called
	   from programs like solveray (formal solution along one specific
	   ray) in cases with moving atmospheres (angle-dependent opacity).

 Note: If equilibria_only is set to TRUE only the electron density,
	   LTE populations and collisions, and chemical equilibria are
	   evaluated.

 Note: Record numbers stored in atmos.backgrrecno refer to records
	   of the size atmos.Nspace. If a polarized line is present 9
	   (4 + 4 + 1, no magneto-optical effects), or 12 (7 + 4 +1, with
	   magneto-optical effects) records are used, otherwise 3 (1 + 1 + 1)
	   records.
	   If the atmosphere is moving, or if a polarized line is present
	   data is stored for each angle and wavelength, otherwise data
	   is stored once for each wavelength only.
	   --                                              -------------- */

#include <fcntl.h>
#include <limits.h>
#include <math.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "background.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"
#include "xdr.h"
#include "parallel.h"

#define COMMENT_CHAR "#"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum   spectrum;
extern InputData  input;
extern char       messageStr[];

/* ------- begin -------------------------- backgroundRange.c ------- */

/* --- Per-worker scratch space of the main loop in Background(). The
	   quantities have the roles described in Background(). In the static
	   unpolarized case chi_ai, eta_ai, and sca_ai alias chi_c, eta_c, and
	   sca_c. pending holds the indices into atmos.backgrrecno of the
	   (mu, to_obs) entries of the current wavelength that have not been
	   written yet. Memory per worker, in units of Nspace doubles:
	   7 for a static unpolarized atmosphere, 10 if it is moving, 22 with
	   Stokes (28 with magneto-optical effects), plus the PassiveBBCache
	   and 2 * Nrays ints. --                          -------------- */

typedef struct {
	double *chi, *eta, *scatt, *chip, *chi_c, *eta_c, *sca_c, *chip_c,
		*chi_ai, *eta_ai, *sca_ai, *Bnu;
	PassiveBBCache *passivecache;
	int            *pending, Npending;
} BackgroundScratch;

/* --- Data shared by the workers: read-only inputs, the scratch space
	   per worker, and the counter of records in background.dat. -- --- */

typedef struct {
	BackgroundScratch **scratch;
	bool_t              do_fudge;
	int                 Nfudge;
	double             *lambda_fudge, **fudge, *thomson;
	Atom               *He;
	atomic_llong        backgrrecno;
} BackgroundJob;

static BackgroundScratch *backgroundScratchCreate(void)
{
	int                NrecStokes;
	BackgroundScratch *sc;

	/* --- Allocate the temporary storage space of one worker. See
		   Background() for what the quantities are used for -- ------- */

	NrecStokes = (atmos.Stokes) ? 4 : 1;

	sc = (BackgroundScratch *) malloc(sizeof(BackgroundScratch));

	sc->chi_c = (double *) malloc(NrecStokes * atmos.Nspace * sizeof(double));
	sc->eta_c = (double *) malloc(NrecStokes * atmos.Nspace * sizeof(double));
	sc->sca_c = (double *) malloc(atmos.Nspace * sizeof(double));

	sc->chi   = (double *) malloc(NrecStokes * atmos.Nspace * sizeof(double));
	sc->eta   = (double *) malloc(NrecStokes * atmos.Nspace * sizeof(double));
	sc->scatt = (double *) malloc(atmos.Nspace * sizeof(double));

	if (atmos.Stokes && input.magneto_optical) {
		sc->chip   = (double *) malloc(3 * atmos.Nspace * sizeof(double));
		sc->chip_c = (double *) malloc(3 * atmos.Nspace * sizeof(double));
	} else {
		sc->chip   = NULL;
		sc->chip_c = NULL;
	}

	if (atmos.moving || atmos.Stokes) {
		sc->chi_ai = (double *) malloc(atmos.Nspace * sizeof(double));
		sc->eta_ai = (double *) malloc(atmos.Nspace * sizeof(double));
		sc->sca_ai = (double *) malloc(atmos.Nspace * sizeof(double));
	} else {
		sc->chi_ai = sc->chi_c;
		sc->eta_ai = sc->eta_c;
		sc->sca_ai = sc->sca_c;
	}
	sc->Bnu = (double *) malloc(atmos.Nspace * sizeof(double));

	sc->passivecache = passiveBBCacheCreate();
	sc->pending      = (int *) malloc(2 * atmos.Nrays * sizeof(int));
	sc->Npending     = 0;

	return sc;
}

static void backgroundScratchFree(BackgroundScratch *sc)
{
	free(sc->chi);
	free(sc->eta);
	free(sc->scatt);
	free(sc->Bnu);
	free(sc->chi_c);
	free(sc->eta_c);
	free(sc->sca_c);

	if (atmos.moving || atmos.Stokes) {
		free(sc->chi_ai);
		free(sc->eta_ai);
		free(sc->sca_ai);
	}
	if (atmos.Stokes && input.magneto_optical) {
		free(sc->chip);
		free(sc->chip_c);
	}
	passiveBBCacheFree(sc->passivecache);
	free(sc->pending);
	free(sc);
}

static int reserveRecords(BackgroundJob *job, int nspect)
{
	const char routineName[] = "reserveRecords";

	long long recno;
	int       Nrec;

	/* --- Reserve the records for a block of background data of
		   wavelength nspect and return the number of the first. The
		   counter is 64-bit, but atmos.backgrrecno is int. With one
		   worker the blocks follow each other in wavelength order. -- */

	Nrec  = backgroundRecords(nspect);
	recno = atomic_fetch_add(&job->backgrrecno, Nrec);
	if (recno + Nrec > INT_MAX)
		Error(ERROR_LEVEL_2, routineName,
			  "Record number of background data exceeds INT_MAX");

	return (int) recno;
}

static void backgroundRange(int begin, int end, int nt, void *arg)
{
	const char   routineName[] = "backgroundRange";
	register int k, nspect, n, mu, to_obs;

	int                index, recno, Nwrite, NrecStokes;
	double             wavelength, Hmin_fudge, scatt_fudge, metal_fudge;
	flags              backgrflags;
	BackgroundJob     *job = (BackgroundJob *) arg;
	BackgroundScratch *sc  = job->scratch[nt];

	bool_t  do_fudge     = job->do_fudge;
	int     Nfudge       = job->Nfudge;
	double *lambda_fudge = job->lambda_fudge, **fudge = job->fudge,
		   *thomson = job->thomson;
	Atom   *He      = job->He;
	double *chi = sc->chi, *eta = sc->eta, *scatt = sc->scatt, *chip = sc->chip,
		   *chi_c = sc->chi_c, *eta_c = sc->eta_c, *sca_c = sc->sca_c,
		   *chip_c = sc->chip_c, *chi_ai = sc->chi_ai, *eta_ai = sc->eta_ai,
		   *sca_ai = sc->sca_ai, *Bnu = sc->Bnu;

	for (nspect = begin; nspect < end; nspect++) {
		sc->Npending = 0;
		wavelength   = spectrum.lambda[nspect];

		/* --- The Planck function at this wavelength --   -------------- */

		Planck(atmos.Nspace, atmos.T, wavelength, Bnu);

		/* --- Initialize the flags for this wavelength -- -------------- */

		atmos.backgrflags[nspect].hasline     = FALSE;
		atmos.backgrflags[nspect].ispolarized = FALSE;

		/* --- Initialize angle-independent quantities --  -------------- */

		for (k = 0; k < atmos.Nspace; k++) {
			chi_ai[k] = 0.0;
			eta_ai[k] = 0.0;
			sca_ai[k] = thomson[k];
		}
		/* --- Negative hydrogen ion, bound-free and free-free -- ------- */

		if (Hminus_bf(wavelength, chi, eta)) {
			for (k = 0; k < atmos.Nspace; k++) {
				chi_ai[k] += chi[k];
				eta_ai[k] += eta[k];
			}
		}
		if (Hminus_ff(wavelength, chi)) {
			for (k = 0; k < atmos.Nspace; k++) {
				chi_ai[k] += chi[k];
				eta_ai[k] += chi[k] * Bnu[k];
			}
		}
		/* --- Opacity fudge factors, applied to Hminus opacity -- ------ */

		if (do_fudge) {
			Linear(Nfudge, lambda_fudge, fudge[0],
				   1, &wavelength, &Hmin_fudge, FALSE);
			for (k = 0; k < atmos.Nspace; k++) {
				chi_ai[k] *= Hmin_fudge;
				eta_ai[k] *= Hmin_fudge;
			}
		}
		/* --- Opacities from bound-free transitions in OH and CH -- ---- */

		if (OH_bf_opac(wavelength, chi, eta)) {
			for (k = 0; k < atmos.Nspace; k++) {
				chi_ai[k] += chi[k];
				eta_ai[k] += eta[k];
			}
		}
		if (CH_bf_opac(wavelength, chi, eta)) {
			for (k = 0; k < atmos.Nspace; k++) {
				chi_ai[k] += chi[k];
				eta_ai[k] += eta[k];
			}
		}
		/* --- Neutral Hydrogen Bound-Free and Free-Free --  ------------ */

		if (Hydrogen_bf(wavelength, chi, eta)) {
			for (k = 0; k < atmos.Nspace; k++) {
				chi_ai[k] += chi[k];
				eta_ai[k] += eta[k];
			}
		}
		Hydrogen_ff(wavelength, chi);
		for (k = 0; k < atmos.Nspace; k++) {
			chi_ai[k] += chi[k];
			eta_ai[k] += chi[k] * Bnu[k];
		}
		/* --- Rayleigh scattering by neutral hydrogen --  -------------- */

		if (Rayleigh(wavelength, atmos.H, scatt)) {
			for (k = 0; k < atmos.Nspace; k++) {
				sca_ai[k] += scatt[k];
			}
		}
		/* --- Rayleigh scattering by neutral helium --    -------------- */

		if (He && Rayleigh(wavelength, He, scatt)) {
			for (k = 0; k < atmos.Nspace; k++) {
				sca_ai[k] += scatt[k];
			}
		}
		/* --- Absorption by H + H^+ (referred to as H2plus free-free) -- */

		if (H2plus_ff(wavelength, chi)) {
			for (k = 0; k < atmos.Nspace; k++) {
				chi_ai[k] += chi[k];
				eta_ai[k] += chi[k] * Bnu[k];
			}
		}
		/* --- Rayleigh scattering and free-free absorption by
			   molecular hydrogen --                       -------------- */

		if (Rayleigh_H2(wavelength, scatt)) {
			for (k = 0; k < atmos.Nspace; k++) {
				sca_ai[k] += scatt[k];
			}
		}
		if (H2minus_ff(wavelength, chi)) {
			for (k = 0; k < atmos.Nspace; k++) {
				chi_ai[k] += chi[k];
				eta_ai[k] += chi[k] * Bnu[k];
			}
		}
		/* --- Bound-Free opacities due to ``metals'' --   -------------- */

		if (do_fudge) {
			Linear(Nfudge, lambda_fudge, fudge[2],
				   1, &wavelength, &metal_fudge, FALSE);
		} else {
			metal_fudge = 1.0;
		}
		/* --- Note: Hydrogen bound-free opacities are calculated in
			   routine Hydrogen_bf --                      -------------- */

		Metal_bf(wavelength, atmos.Natom - 1, atmos.atoms + 1, chi, eta);
		for (k = 0; k < atmos.Nspace; k++) {
			chi_ai[k] += chi[k] * metal_fudge;
			eta_ai[k] += eta[k] * metal_fudge;
		}
		/* --- Add the scattering opacity to the absorption part to store
			   the total opacity --                        -------------- */

		if (do_fudge) {
			Linear(Nfudge, lambda_fudge, fudge[1],
				   1, &wavelength, &scatt_fudge, FALSE);
		} else {
			scatt_fudge = 1.0;
		}
		for (k = 0; k < atmos.Nspace; k++) {
			sca_ai[k] *= scatt_fudge;
			chi_ai[k] += sca_ai[k];
		}
		/* --- Now the contributions that may be angle-dependent due to the
			   presence of atomic or molecular lines --    -------------- */

		if (atmos.moving || atmos.Stokes) {
			for (mu = 0; mu < atmos.Nrays; mu++) {
				for (to_obs = 0; to_obs <= 1; to_obs++) {
					index = 2 * (nspect * atmos.Nrays + mu) + to_obs;

					/* --- First, copy the angle-independent parts -- --------- */

					for (k = 0; k < atmos.Nspace; k++) {
						chi_c[k] = chi_ai[k];
						eta_c[k] = eta_ai[k];
						sca_c[k] = sca_ai[k];
					}

					/* --- Zero the polarized quantities, if necessary -- ----- */

					if (atmos.Stokes) {
						for (k = atmos.Nspace; k < 4 * atmos.Nspace; k++) {
							chi_c[k] = 0.0;
							eta_c[k] = 0.0;
						}
						if (input.magneto_optical)
							for (k = 0; k < 3 * atmos.Nspace; k++) chip_c[k] = 0.0;
					}
					/* --- Add opacity from passive atomic lines (including
						   hydrogen) --                          -------------- */

					if (input.allow_passive_bb) {
						backgrflags = passive_bb(sc->passivecache, wavelength, nspect, mu,
												 to_obs, chi, eta, chip);
						if (backgrflags.hasline) {
							atmos.backgrflags[nspect].hasline = TRUE;
							if (backgrflags.ispolarized) {
								NrecStokes                            = 4;
								atmos.backgrflags[nspect].ispolarized = TRUE;
								if (input.magneto_optical) {
									for (k = 0; k < 3 * atmos.Nspace; k++)
										chip_c[k] += chip[k];
								}
							} else
								NrecStokes = 1;

							for (k = 0; k < NrecStokes * atmos.Nspace; k++) {
								chi_c[k] += chi[k];
								eta_c[k] += eta[k];
							}
						}
					}
					/* --- Add opacity from Kurucz line list --  -------------- */

					if (atmos.Nrlk > 0) {
						backgrflags = rlk_opacity(wavelength, nspect, mu, to_obs,
												  chi, eta, scatt, chip);
						if (backgrflags.hasline) {
							atmos.backgrflags[nspect].hasline = TRUE;
							if (backgrflags.ispolarized) {
								NrecStokes                            = 4;
								atmos.backgrflags[nspect].ispolarized = TRUE;
								if (input.magneto_optical) {
									for (k = 0; k < 3 * atmos.Nspace; k++)
										chip_c[k] += chip[k];
								}
							} else
								NrecStokes = 1;

							for (k = 0; k < NrecStokes * atmos.Nspace; k++) {
								chi_c[k] += chi[k];
								eta_c[k] += eta[k];
							}
							if (input.rlkscatter) {
								for (k = 0; k < atmos.Nspace; k++) {
									sca_c[k] += scatt[k];
									chi_c[k] += scatt[k];
								}
							}
						}
					}
					/* --- Add opacity from molecular lines --   -------------- */

					backgrflags = MolecularOpacity(wavelength, nspect, mu, to_obs,
												   chi, eta, chip);
					if (backgrflags.hasline) {
						atmos.backgrflags[nspect].hasline = TRUE;
						if (backgrflags.ispolarized) {
							NrecStokes                            = 4;
							atmos.backgrflags[nspect].ispolarized = TRUE;
							if (input.magneto_optical) {
								for (k = 0; k < 3 * atmos.Nspace; k++)
									chip_c[k] += chip[k];
							}
						} else
							NrecStokes = 1;

						for (k = 0; k < NrecStokes * atmos.Nspace; k++) {
							chi_c[k] += chi[k];
							eta_c[k] += eta[k];
						}
					}
					/* --- Store angle-dependent results only if at least one
							   line was found at this wavelength --  -------------- */

					if ((mu == atmos.Nrays - 1 && to_obs) ||
						(atmos.backgrflags[nspect].hasline &&
						 (atmos.moving || atmos.backgrflags[nspect].ispolarized))) {

						/* --- Reserve the records of this block. Entries
							   that were not written since the previous
							   block of this wavelength point to this one,
							   as they did when the counter ran in
							   wavelength order -- --------------------- */

						recno = reserveRecords(job, nspect);
						for (n = 0; n < sc->Npending; n++)
							atmos.backgrrecno[sc->pending[n]] = recno;
						sc->Npending             = 0;
						atmos.backgrrecno[index] = recno;
						Nwrite                   = writeBackground(nspect, mu, to_obs,
																   chi_c, eta_c, sca_c, chip_c);
						if (Nwrite != backgroundRecords(nspect))
							Error(ERROR_LEVEL_2, routineName,
								  "Reserved and written number of records differ");
					} else
						sc->pending[sc->Npending++] = index;
				}
			}
		} else {
			/* --- Angle-independent case. First, add opacity from passive
			   atomic lines (including hydrogen) --      -------------- */

			if (input.allow_passive_bb) {
				backgrflags = passive_bb(sc->passivecache, wavelength, nspect, 0,
										 TRUE, chi, eta, NULL);
				if (backgrflags.hasline) {
					atmos.backgrflags[nspect].hasline = TRUE;
					for (k = 0; k < atmos.Nspace; k++) {
						chi_c[k] += chi[k];
						eta_c[k] += eta[k];
					}
				}
			}
			/* --- Add opacity from Kurucz line list --      -------------- */

			if (atmos.Nrlk > 0) {
				backgrflags = rlk_opacity(wavelength, nspect, 0, TRUE,
										  chi, eta, scatt, NULL);
				if (backgrflags.hasline) {
					atmos.backgrflags[nspect].hasline = TRUE;
					for (k = 0; k < atmos.Nspace; k++) {
						chi_c[k] += chi[k];
						eta_c[k] += eta[k];
					}
					if (input.rlkscatter) {
						for (k = 0; k < atmos.Nspace; k++) {
							sca_c[k] += scatt[k];
							chi_c[k] += scatt[k];
						}
					}
				}
			}
			/* --- Add opacity from molecular lines --       -------------- */

			backgrflags = MolecularOpacity(wavelength, nspect, 0, TRUE,
										   chi, eta, NULL);
			if (backgrflags.hasline) {
				atmos.backgrflags[nspect].hasline = TRUE;
				for (k = 0; k < atmos.Nspace; k++) {
					chi_c[k] += chi[k];
					eta_c[k] += eta[k];
				}
			}
			/* --- Store results --                          -------------- */

			atmos.backgrrecno[nspect] = reserveRecords(job, nspect);
			Nwrite                    = writeBackground(nspect, 0, 0,
														chi_c, eta_c, sca_c, NULL);
			if (Nwrite != backgroundRecords(nspect))
				Error(ERROR_LEVEL_2, routineName,
					  "Reserved and written number of records differ");
		}
	}
}
/* ------- end ---------------------------- backgroundRange.c ------- */

/* ------- begin -------------------------- Background.c ------------ */

void Background(bool_t write_analyze_output, bool_t equilibria_only)
{
	const char   routineName[] = "Background";
	register int nspect, n;

	static int    ne_iter = 0;
	char          inputLine[MAX_LINE_SIZE];
	bool_t        exit_on_EOF, do_fudge, fromscratch;
	int           Nfudge, Nworker, chunk;
	double       *thomson, *lambda_fudge = NULL, **fudge = NULL;
	Atom         *He;
	FILE         *fp_fudge;
	BackgroundJob job;

	getCPU(2, TIME_START, NULL);

	if (input.solve_ne == ONCE || input.solve_ne == ITERATION) {
		fromscratch = (input.solve_ne == ONCE ||
					   (input.solve_ne == ITERATION && ne_iter == 0))
						  ? TRUE
						  : FALSE;
		Solve_ne(atmos.ne, fromscratch);
		ne_iter++;
	}
	SetLTEQuantities();

	if (input.NonICE)
		readMolecules(MOLECULAR_CONCENTRATION_FILE);
	else
		ChemicalEquilibrium(N_MAX_CHEM_ITER, CHEM_ITER_LIMIT);

	if (equilibria_only) {

		/* --- If we only need ne, LTE populations and collisions, and
			   chemical equilibrium leave here --          -------------- */

		getCPU(2, TIME_POLL, "Total Background");
		return;
	}

	if (input.old_background) {
		if ((atmos.fd_background =
				 open(input.background_File, O_RDONLY, 0)) == -1) {
			sprintf(messageStr, "Unable to open input file %s",
					input.background_File);
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
		readBRS();
		return;
	}

	getCPU(3, TIME_START, NULL);

	do_fudge = FALSE;

	if (strcmp(input.fudgeData, "none")) {
		do_fudge = TRUE;

		/* --- Read wavelength-dependent fudge factors to compensate for
			   missing UV backround line haze --           -------------- */

		if ((fp_fudge = fopen(input.fudgeData, "r")) == NULL) {
			sprintf(messageStr, "Unable to open input file %s", input.fudgeData);
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
		sprintf(messageStr,
				"\n-Fudging background opacities with file\n  %s\n\n",
				input.fudgeData);
		Error(MESSAGE, routineName, messageStr);

		getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
		sscanf(inputLine, "%d", &Nfudge);
		lambda_fudge = (double *) malloc(Nfudge * sizeof(double));
		fudge        = matrix_double(3, Nfudge);
		for (n = 0; n < Nfudge; n++) {
			getLine(fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
			sscanf(inputLine, "%lf %lf %lf %lf", &lambda_fudge[n],
				   &fudge[0][n], &fudge[1][n], &fudge[2][n]);
		}
		for (n = 0; n < 3 * Nfudge; n++) fudge[0][n] += 1.0;
		fclose(fp_fudge);
	}

	/* --- The temporary storage space is allocated per worker, further
		   down, in backgroundScratchCreate. The quantities are used
		   for the following purposes:

		 - chi, eta, scatt: Get contributions to opacity, emissivity,
		   and scattering opacity, respectively, from a specific process
		   for a given wavelength and possibly angle.

		 - chi_c, eta_c, sca_c: Collect the total opacity, emissivity
		   and scattering opacity for a given wavelength and possibly
		   angle.

		 - chi_ai, eta_ai: Collect the angle-independent part of
		   opacity and emissivity for each wavelength so that these
		   need not be recalculated in an angle-dependent case.
		   When the atmosphere is not moving and has no magnetic fields
		   these just point to the total quantities chi_c and eta_c.

	 Note: In case of magnetic fields in the atmosphere chi, eta and
		   chip, and chi_c, eta_c and chip_c contain all four Stokes
		   parameters, and should be allocated a 4 and 3 times larger
		   storage space, respectively.
		   --                                            -------------- */

	/* --- Thomson scattering by free electrons is wavelength independent
		   in non-relativistic limit so we compute it only once -- ---- */

	thomson = (double *) malloc(atmos.Nspace * sizeof(double));
	Thomson(thomson);

	/* --- Check whether an atomic model is present for He -- --------- */

	He = (atmos.elements[1].model) ? atmos.elements[1].model : NULL;

	/* --- Read background files from Kurucz data file -- ------------- */

	atmos.Nrlk = 0;
	readKuruczLines(input.KuruczData);
	if (atmos.Nrlk > 0) {
		qsort(atmos.rlk_lines, atmos.Nrlk, sizeof(RLK_Line), rlk_ascend);
	}
	/* --- Allocate memory for the boolean array that stores whether
		   a wavelength overlaps with a Bound-Bound transition in the
		   background, or whether it is polarized --     -------------- */

	atmos.backgrflags = (flags *) malloc(spectrum.Nspect * sizeof(flags));
	for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
		atmos.backgrflags[nspect].hasline     = FALSE;
		atmos.backgrflags[nspect].ispolarized = FALSE;
	}
	/* --- Allocate memory for the list of record numbers that specifies
		   for each wavelength where to find the background opacity,
		   scattering opacity, and emissivity --         -------------- */

	if (atmos.moving || atmos.Stokes) {
		atmos.backgrrecno =
			(int *) malloc(2 * spectrum.Nspect * atmos.Nrays * sizeof(int));
	} else
		atmos.backgrrecno = (int *) malloc(spectrum.Nspect * sizeof(int));

	/* --- Open output file for background opacity, emissivity,
		   scattering. The values are stored as float (M-1); truncate,
		   so that a file left by an earlier run (possibly in double)
		   does not keep its size --                     -------------- */

	atmos.backgrbytes  = sizeof(float);
	atmos.backgretaexp = BACKGR_ETA_EXP;
	if ((atmos.fd_background =
			 open(input.background_File, O_RDWR | O_CREAT | O_TRUNC,
				  PERMISSIONS)) == -1) {
		sprintf(messageStr, "Unable to open output file %s",
				input.background_File);
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	/* --- One-time setup for the opacity routines called in the main
		   loop, and the cache for passive_bb --         -------------- */

	backgroundPrepare();
	if (atmos.Nrlk > 0) rlk_prepare();

	/* --- Go through the spectrum and add the different opacity and
		   emissivity contributions. This is the main loop, over
		   wavelengths. Each worker has its own scratch space, and the
		   wavelengths are handed out in chunks with dynamic scheduling,
		   because the cost per wavelength varies a lot (dense line
		   forests versus continuum). The chunk is about 1/8 of a
		   thread's share of the spectrum, at most 32 wavelengths: small
		   enough to balance the load, large enough to keep the counter
		   uncontended. With input.Nthreads <= 1 everything runs in this
		   thread, in wavelength order, and the layout of background.dat
		   is the sequential one. With more threads the order of the
		   blocks in the file varies from run to run; readers find them
		   through atmos.backgrrecno. --                 -------------- */

	Nworker = (input.Nthreads > 1) ? input.Nthreads : 1;
	chunk   = spectrum.Nspect / (8 * Nworker);
	if (chunk > 32) chunk = 32;
	if (chunk < 1) chunk = 1;

	job.scratch      = (BackgroundScratch **) malloc(Nworker * sizeof(BackgroundScratch *));
	job.do_fudge     = do_fudge;
	job.Nfudge       = Nfudge;
	job.lambda_fudge = lambda_fudge;
	job.fudge        = fudge;
	job.thomson      = thomson;
	job.He           = He;
	atomic_init(&job.backgrrecno, 0);
	for (n = 0; n < Nworker; n++) job.scratch[n] = backgroundScratchCreate();

	parallelForDynamic(spectrum.Nspect, chunk, Nworker, backgroundRange, &job);

	if (atmos.Nrlk > 0) rlk_cleanup();

	if (write_analyze_output) {
		/* --- Write background record structure --          ------------ */

		writeBRS();

		/* --- Write out the metals and molecules --         ------------ */

		writeMetals("metals.out");
		writeMolecules(MOLECULAR_CONCENTRATION_FILE);
	}
	/* --- Clean up but keep H, H2, and active atom and/or molecule
		   if appropriate --                               ------------ */

	if (atmos.Natom > 1) {
		for (n = 1; n < atmos.Natom; n++)
			if (!atmos.atoms[n].active &&
				!atmos.hydrostatic &&
				input.solve_ne < ITERATION)
				freeAtom(&atmos.atoms[n]);
	}
	if (atmos.Nmolecule > 1) {
		for (n = 1; n < atmos.Nmolecule; n++)
			if (!atmos.molecules[n].active &&
				!atmos.hydrostatic &&
				input.solve_ne < ITERATION)
				freeMolecule(&atmos.molecules[n]);
	}

	/* --- The partition functions are needed again when Background is
		   called again (electron density iteration, hydrostatic
		   equilibrium), like the passive atoms above --  ------------- */

	if (strcmp(input.KuruczData, "none") &&
		!atmos.hydrostatic && input.solve_ne < ITERATION) {
		free(atmos.Tpf);
		atmos.Tpf = NULL;
		for (n = 0; n < atmos.Nelem; n++) {
			free(atmos.elements[n].ionpot);
			freeMatrix((void **) atmos.elements[n].pf);
			if (atmos.elements[n].n)
				freeMatrix((void **) atmos.elements[n].n);
		}
	}
	getCPU(3, TIME_POLL, "Background Opacity");

	/* --- Free the temporary space set up by backgroundPrepare -- --- */

	backgroundCleanup();

	for (n = 0; n < Nworker; n++) backgroundScratchFree(job.scratch[n]);
	free(job.scratch);
	free(thomson);

	if (do_fudge) {
		free(lambda_fudge);
		freeMatrix((void **) fudge);
	}

	for (n = 0; n < atmos.Nrlk; n++) {
		if (atmos.rlk_lines[n].zm != NULL) freeZeeman(atmos.rlk_lines[n].zm);
	}
	free(atmos.rlk_lines);

	getCPU(2, TIME_POLL, "Total Background");
}
/* ------- end ---------------------------- Background.c ------------ */

/* ------- begin -------------------------- backgroundPrepare.c ----- */

void backgroundPrepare(void)
{
	/* --- One-time setup of the state that the continuum and passive
		   line opacity routines would otherwise initialize lazily at
		   the first wavelength. After this they only read shared data,
		   and can be called concurrently for different wavelengths.
		   Call after the atmosphere (T), the passive atoms, and the
		   molecules are in place. --                  -------------- */

	hydrogenPrepare();
	metalPrepare();
	molecularPrepare();
}
/* ------- end ---------------------------- backgroundPrepare.c ----- */

/* ------- begin -------------------------- backgroundCleanup.c ----- */

void backgroundCleanup(void)
{
	/* --- Free the temporary space allocated by backgroundPrepare. The
		   spline coefficients of the passive atoms are freed with the
		   atoms. --                                   -------------- */

	hydrogenCleanup();
}
/* ------- end ---------------------------- backgroundCleanup.c ----- */
