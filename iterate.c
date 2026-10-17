/* ------- file: -------------------------- iterate.c ---------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek  (huitenbroek@nso.edu)
	   Last modified: Tue Feb 20 13:13:59 2024 --

	   --------------------------                      ----------RH-- */

/* --- Main iteration routine --                       -------------- */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "accelerate.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"

typedef struct {
	bool_t eval_operator, redistribute;
	int    nt, Nthreads, lambda_max;
	double dJmax;
} threadinfo;

/* --- Function prototypes --                          -------------- */

void *Formal_pthread(void *argument);

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum   spectrum;
extern InputData  input;
extern char       messageStr[];

/* ------- begin -------------------------- Iterate.c --------------- */

void Iterate(int NmaxIter, double iterLimit)
{
	const char   routineName[] = "Iterate";
	register int niter, nact;

	bool_t eval_operator, write_analyze_output, equilibria_only,
		converged = FALSE, stopped = FALSE;
	int       Ngorder;
	double    dpopsmax = 0.0, PRDiterlimit;
	Atom     *atom;
	Molecule *molecule;

	if (NmaxIter <= 0) return;
	getCPU(1, TIME_START, NULL);

	/* --- Initialize structures for Ng acceleration of population
		   convergence --                                  ------------ */

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom       = atmos.activeatoms[nact];
		atom->Ng_n = NgInit(atom->Nlevel * atmos.Nspace, input.Ngdelay,
							input.Ngorder, input.Ngperiod, atom->n[0]);
	}
	for (nact = 0; nact < atmos.Nactivemol; nact++) {
		molecule = atmos.activemols[nact];
		Ngorder  = (input.accelerate_mols) ? input.Ngorder : 0;

		molecule->Ng_nv = NgInit(molecule->Nv * atmos.Nspace, input.Ngdelay,
								 Ngorder, input.Ngperiod, molecule->nv[0]);
	}

	/* --- Start of the main iteration loop --             ------------ */

	niter = 1;
	while (niter <= NmaxIter ||
		   (input.CR_Nstep > 0 && niter <= input.CR_Nstep)) {

		if (StopRequested()) {
			stopped = TRUE;
			break;
		}

		getCPU(2, TIME_START, NULL);

		for (nact = 0; nact < atmos.Nactiveatom; nact++) {
			initGammaAtom(atmos.activeatoms[nact], niter);
			zeroGammaThreads(atmos.activeatoms[nact]);
		}
		for (nact = 0; nact < atmos.Nactivemol; nact++) {
			initGammaMolecule(atmos.activemols[nact]);
			zeroGammaThreadsMolecule(atmos.activemols[nact]);
		}

		/* --- Formal solution for all wavelengths --      -------------- */

		solveSpectrum(eval_operator = TRUE, FALSE);

		/* --- Merge each thread's Gamma contributions (accumulated with
			   no locking during the sweep above) into the shared Gamma,
			   which already carries this iteration's collisional seed
			   from initGammaAtom/initGammaMolecule above --  ---------- */

		for (nact = 0; nact < atmos.Nactiveatom; nact++)
			mergeGammaThreads(atmos.activeatoms[nact]);
		for (nact = 0; nact < atmos.Nactivemol; nact++)
			mergeGammaThreadsMolecule(atmos.activemols[nact]);

		/* --- Solve statistical equilibrium equations --  -------------- */

		sprintf(messageStr, "\n -- Iteration %3d\n", niter);
		Error(MESSAGE, routineName, messageStr);
		dpopsmax = updatePopulations(niter);

		if (atmos.NPRDactive > 0) {

			/* --- Redistribute intensity in PRD lines if necessary -- ---- */

			if (input.PRDiterLimit < 0.0)
				PRDiterlimit = MAX(dpopsmax, -input.PRDiterLimit);
			else
				PRDiterlimit = input.PRDiterLimit;

			Redistribute(input.PRD_NmaxIter, PRDiterlimit);
		}

		sprintf(messageStr, "Total Iteration %3d", niter);
		getCPU(2, TIME_POLL, messageStr);

		/* --- With collisional relaxation the collision rates are scaled
			   up until iteration CR_Nstep, so populations close to LTE
			   change little and are not converged. Test for convergence
			   only from iteration CR_Nstep on, where the scaling factor is
			   1, and restart Ng acceleration after every scaled iteration
			   so that it only extrapolates from iterations of the actual
			   problem --                                   ------------ */

		if (niter < input.CR_Nstep) {
			for (nact = 0; nact < atmos.Nactiveatom; nact++) {
				atom = atmos.activeatoms[nact];
				NgFree(atom->Ng_n);
				atom->Ng_n = NgInit(atom->Nlevel * atmos.Nspace, input.Ngdelay,
									input.Ngorder, input.Ngperiod, atom->n[0]);
			}
			for (nact = 0; nact < atmos.Nactivemol; nact++) {
				molecule = atmos.activemols[nact];
				Ngorder  = (input.accelerate_mols) ? input.Ngorder : 0;
				NgFree(molecule->Ng_nv);
				molecule->Ng_nv = NgInit(molecule->Nv * atmos.Nspace, input.Ngdelay,
										 Ngorder, input.Ngperiod, molecule->nv[0]);
			}
		} else if (dpopsmax < iterLimit) {
			converged = TRUE;
			break;
		}
		niter++;

		if (input.solve_ne == ITERATION)
			Background(write_analyze_output = TRUE, equilibria_only = FALSE);

		if (atmos.hydrostatic) {
			if (!atmos.atoms[0].active) {
				sprintf(messageStr, "Can only perform hydrostatic equilibrium"
									" for hydrogen active");
				Error(ERROR_LEVEL_2, routineName, messageStr);
			}
			Hydrostatic(N_MAX_HSE_ITER, HSE_ITER_LIMIT);
		}
	}

	/* --- Say whether the populations converged. A run that stops at
		   N_MAX_ITER still writes its output, so make it clear that
		   this output is not a converged solution --   -------------- */

	if (converged) {
		sprintf(messageStr, "\n -- Converged in %d iterations: max change"
							" %.4E < ITER_LIMIT = %.4E\n",
				niter, dpopsmax, iterLimit);
		Error(MESSAGE, routineName, messageStr);
	} else if (stopped) {
		sprintf(messageStr, "Stopped on request (STOP_RH) after %d iterations:"
							" NOT CONVERGED.\n The output is not a converged solution",
				niter - 1);
		Error(WARNING, routineName, messageStr);
		SetExitFlag(RH_EXIT_NOT_CONVERGED);
	} else {
		sprintf(messageStr, "NOT CONVERGED after %d iterations: max change"
							" %.4E >= ITER_LIMIT = %.4E.\n The output is not a converged"
							" solution; increase N_MAX_ITER or see what limits convergence",
				niter - 1, dpopsmax, iterLimit);
		Error(WARNING, routineName, messageStr);
		SetExitFlag(RH_EXIT_NOT_CONVERGED);
	}

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];
		freeMatrix((void **) atom->Gamma);
		NgFree(atom->Ng_n);
	}
	for (nact = 0; nact < atmos.Nactivemol; nact++) {
		molecule = atmos.activemols[nact];
		freeMatrix((void **) molecule->Gamma);
		NgFree(molecule->Ng_nv);
	}

	getCPU(1, TIME_POLL, "Iteration Total");
}
/* ------- end ---------------------------- Iterate.c --------------- */

/* ------- begin -------------------------- solveSpectrum.c --------- */

double solveSpectrum(bool_t eval_operator, bool_t redistribute)
{
	register int nspect, nt, k;

	int         lambda_max;
	double      dJ, dJmax;
	pthread_t  *thread_id;
	threadinfo *ti;

	/* --- Administers the formal solution for each wavelength. When
		   input.Nthreads > 1 the solutions are performed concurrently
		   in Nthreads threads. These are POSIX style threads.

	  See: - David R. Butenhof, Programming with POSIX threads,
			 Addison & Wesley.

		   - Multithreaded Programming Guide, http://sun.docs.com
			 (search for POSIX threads).

		   With eval_operator == TRUE and more than one thread, the
		   threads leave their Gamma contributions in their own copies
		   (rhth[nt].Gamma); the caller has to run mergeGammaThreads and
		   mergeGammaThreadsMolecule before using Gamma (Iterate does).
		   The radiative rates Rij and Rji are merged here. With one
		   thread everything is accumulated in place, as in a serial
		   code, and the merges do nothing.

		   When solveSpectrum is called with redistribute == TRUE only
		   wavelengths that have an active PRD line are solved. The
		   redistribute key is passed to the addtoRates routine via
		   Formal so that only the radiative rates of PRD lines are
		   updated. These are needed for the emission profile ratio \rho.
		   --                                            -------------- */

	getCPU(3, TIME_START, NULL);

	/* --- First zero the radiative rates --             -------------- */

	zeroRates(redistribute);
	lambda_max = 0;
	dJmax      = 0.0;

	/* --- Zero out J in gas parcel's frame --           -------------- */

	if (spectrum.updateJ &&
		input.PRD_angle_dep == PRD_ANGLE_APPROX &&
		atmos.Nrays > 1 && atmos.NPRDactive > 0) {
		for (k = 0; k < atmos.Nspace; k++) {
			for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
				spectrum.Jgas[nspect][k] = 0.0;
			}
		}
	}

	if (input.Nthreads > 1) {

		/* --- Spawn input.Nthreads worker threads ONCE for this call,
			   rather than once per batch of input.Nthreads wavelengths.
			   Each thread nt processes its own strided slice of the
			   wavelength grid (nt, nt+Nthreads, nt+2*Nthreads, ...) and
			   is joined only after it has worked through all of it.

			   The earlier per-batch create/join meant spectrum.Nspect
			   pthread_create/pthread_join cycles per call to
			   solveSpectrum, regardless of input.Nthreads (batches *
			   threads-per-batch = total wavelengths, always) -- measured
			   to completely swamp any benefit from the parallelism, with
			   wall-clock time flat-to-worse as thread count increased and
			   system time (thread creation/teardown) growing by two
			   orders of magnitude from 1 to 12 threads. This reduces
			   thread creation to input.Nthreads per call. --  ---------- */

		ti        = (threadinfo *) malloc(input.Nthreads * sizeof(threadinfo));
		thread_id = (pthread_t *) malloc(input.Nthreads * sizeof(pthread_t));

		for (nt = 0; nt < input.Nthreads; nt++) {
			ti[nt].eval_operator = eval_operator;
			ti[nt].redistribute  = redistribute;
			ti[nt].nt            = nt;
			ti[nt].Nthreads      = input.Nthreads;

			if (pthread_create(&thread_id[nt], &input.thread_attr,
							   Formal_pthread, &ti[nt])) {

				/* --- Could not start the thread: do its slice here --- ---- */

				Error(WARNING, "solveSpectrum",
					  "Could not create a thread, running its share in the caller");
				Formal_pthread(&ti[nt]);
				thread_id[nt] = 0;
			}
		}
		for (nt = 0; nt < input.Nthreads; nt++) {
			if (thread_id[nt]) pthread_join(thread_id[nt], NULL);
			if (ti[nt].dJmax > dJmax) {
				dJmax      = ti[nt].dJmax;
				lambda_max = ti[nt].lambda_max;
			}
		}
		free(thread_id);
		free(ti);
	} else {

		/* --- Else call the solution for wavelengths sequentially -- --- */

		for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
			if (!redistribute ||
				(redistribute && containsPRDline(&spectrum.as[nspect]))) {
				dJ = Formal(nspect, eval_operator, redistribute);
				if (dJ > dJmax) {
					dJmax      = dJ;
					lambda_max = nspect;
				}
			}
		}
	}

	/* --- Merge each thread's Rij/Rji contributions (accumulated with
		   no locking during the sweep above) into the shared Rij/Rji,
		   before anything downstream (PRDScatter, RADRATE_OUTPUT, ...)
		   reads them --                                  -------------- */

	mergeRates(redistribute);

	sprintf(messageStr, " Spectrum max delta J = %6.4E (lambda#: %d)\n",
			dJmax, lambda_max);
	Error(MESSAGE, NULL, messageStr);

	getCPU(3, TIME_POLL,
		   (eval_operator) ? "Spectrum & Operator" : "Solve Spectrum");

	return dJmax;
}
/* ------- end ---------------------------- solveSpectrum.c --------- */

/* ------- begin -------------------------- Formal_pthread.c -------- */

void *Formal_pthread(void *argument)
{
	threadinfo *ti = (threadinfo *) argument;
	int         nspect;
	double      dJ;

	/* --- Thread wrapper around Formal: works through this thread's
		   strided slice of the wavelength grid --        -------------- */

	ti->dJmax      = 0.0;
	ti->lambda_max = -1;

	for (nspect = ti->nt; nspect < spectrum.Nspect; nspect += ti->Nthreads) {
		if (!ti->redistribute ||
			(ti->redistribute && containsPRDline(&spectrum.as[nspect]))) {
			dJ = Formal(nspect, ti->eval_operator, ti->redistribute);
			if (dJ > ti->dJmax) {
				ti->dJmax      = dJ;
				ti->lambda_max = nspect;
			}
		}
	}

	return (NULL);
}
/* ------- end ---------------------------- Formal_pthread.c -------- */
