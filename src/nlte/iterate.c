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
#include "solver.h"
#include "inputs.h"
#include "negtau.h"
#include "bgscatter.h"
#include "state.h"

typedef struct {
	RHState *rh;
	bool_t   eval_operator, redistribute;
	int      nt, Nthreads, lambda_max;
	double   dJmax;
} threadinfo;

/* --- Function prototypes --                          -------------- */

void *Formal_pthread(void *argument);

/* --- Global variables --                             -------------- */

static void   iterateNg(RHState *rh, int NmaxIter, double iterLimit);
static bool_t moreIterations(RHState *rh, int niter, int NmaxIter);

/* ------- begin -------------------------- Iterate.c --------------- */

/* --- Whether the loop of iterateNg runs iteration niter: up to
	   N_MAX_ITER, and up to CR_Nstep if that is larger -- ------------ */

static bool_t moreIterations(RHState *rh, int niter, int NmaxIter)
{
	return (niter <= NmaxIter || (rh->input.CR_Nstep > 0 && niter <= rh->input.CR_Nstep));
}

void Iterate(RHState *rh, int NmaxIter, double iterLimit)
{
	double t_start;

	/* --- The main iteration with the solver selected by keyword
		   SOLVER: RH's Ng-accelerated MALI iteration (NG, the default in
		   2-D and 3-D), the fixed-point iteration of the map G (FIXED_POINT), or
		   Anderson acceleration of it (ANDERSON, the default in 1-D,
		   S-10) or inexact Newton on
		   G(x) - x (NEWTON_KRYLOV; see solver.c) --                                   -------------- */

	if (NmaxIter <= 0) {
		workPhase(rh, WORK_FINAL);
		return;
	}
	Message(rh, " -- SOLVER = %s%s\n", solverName(rh->input.solver),
			rh->input.solver_default ? " (the default)" : "");

	getCPU(rh, 1, TIME_START, NULL);
	workPhase(rh, WORK_ITERATE);
	t_start = wallClock();

	switch (rh->input.solver) {
	case SOLVER_FIXED_POINT:
		iterateFixedPoint(rh, NmaxIter, iterLimit);
		break;
	case SOLVER_ANDERSON:
		iterateAnderson(rh, NmaxIter, iterLimit);
		break;
	case SOLVER_NEWTON_KRYLOV:
		iterateNewton(rh, NmaxIter, iterLimit);
		break;
	default:
		iterateNg(rh, NmaxIter, iterLimit);
		break;
	}

	workAdd(rh, WORK_ITERATION, wallClock() - t_start);
	workPhase(rh, WORK_FINAL);
	getCPU(rh, 1, TIME_POLL, "Iteration Total");
}
/* ------- end ---------------------------- Iterate.c --------------- */

/* ------- begin -------------------------- iterateNg.c ------------- */

static void iterateNg(RHState *rh, int NmaxIter, double iterLimit)
{
	const char   routineName[] = "Iterate";
	register int niter, nact;

	char   label[MAX_MESSAGE_LENGTH];
	bool_t eval_operator, write_analyze_output, equilibria_only,
		converged = FALSE, stopped = FALSE, more;
	int       Ngorder;
	double    dpopsmax = 0.0, PRDiterlimit;
	Atom     *atom;
	Molecule *molecule;

	/* --- Initialize structures for Ng acceleration of population
		   convergence --                                  ------------ */

	for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
		atom       = rh->atmos.activeatoms[nact];
		atom->Ng_n = NgInit(rh, atom->Nlevel * rh->atmos.Nspace, rh->input.Ngdelay,
							rh->input.Ngorder, rh->input.Ngperiod, atom->n[0]);
	}
	for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
		molecule = rh->atmos.activemols[nact];
		Ngorder  = (rh->input.accelerate_mols) ? rh->input.Ngorder : 0;

		molecule->Ng_nv = NgInit(rh, molecule->Nv * rh->atmos.Nspace, rh->input.Ngdelay,
								 Ngorder, rh->input.Ngperiod, molecule->nv[0]);
	}

	/* --- Start of the main iteration loop --             ------------ */

	niter = 1;
	while (moreIterations(rh, niter, NmaxIter)) {

		/* --- A stop request is read once per iteration: here before the
			   first one, afterwards at the end of the previous iteration
			   (below), before its hydrostatic step -- ------------------- */

		if (niter == 1 && StopRequested(rh)) {
			stopped = TRUE;
			break;
		}

		getCPU(rh, 2, TIME_START, NULL);

		for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
			initGammaAtom(rh, rh->atmos.activeatoms[nact], niter);
			zeroGammaThreads(rh, rh->atmos.activeatoms[nact]);
		}
		for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
			initGammaMolecule(rh, rh->atmos.activemols[nact]);
			zeroGammaThreadsMolecule(rh, rh->atmos.activemols[nact]);
		}

		/* --- Formal solution for all wavelengths --      -------------- */

		solveSpectrum(rh, eval_operator = TRUE, FALSE);

		/* --- Merge each thread's Gamma contributions (accumulated with
			   no locking during the sweep above) into the shared Gamma,
			   which already carries this iteration's collisional seed
			   from initGammaAtom/initGammaMolecule above --  ---------- */

		for (nact = 0; nact < rh->atmos.Nactiveatom; nact++)
			mergeGammaThreads(rh, rh->atmos.activeatoms[nact]);
		for (nact = 0; nact < rh->atmos.Nactivemol; nact++)
			mergeGammaThreadsMolecule(rh, rh->atmos.activemols[nact]);

		/* --- Solve statistical equilibrium equations --  -------------- */

		Message(rh, "\n -- Iteration %3d\n", niter);
		dpopsmax = updatePopulations(rh, niter);

		if (rh->atmos.NPRDactive > 0) {

			/* --- Redistribute intensity in PRD lines if necessary -- ---- */

			if (rh->input.PRDiterLimit < 0.0)
				PRDiterlimit = MAX(dpopsmax, -rh->input.PRDiterLimit);
			else
				PRDiterlimit = rh->input.PRDiterLimit;

			Redistribute(rh, rh->input.PRD_NmaxIter, PRDiterlimit);
		}

		snprintf(label, sizeof(label), "Total Iteration %3d", niter);
		getCPU(rh, 2, TIME_POLL, label);

		/* --- With collisional relaxation the collision rates are scaled
			   up until iteration CR_Nstep, so populations close to LTE
			   change little and are not converged. Test for convergence
			   only from iteration CR_Nstep on, where the scaling factor is
			   1, and restart Ng acceleration after every scaled iteration
			   so that it only extrapolates from iterations of the actual
			   problem --                                   ------------ */

		if (niter < rh->input.CR_Nstep) {
			for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
				atom = rh->atmos.activeatoms[nact];
				NgFree(rh, atom->Ng_n);
				atom->Ng_n = NgInit(rh, atom->Nlevel * rh->atmos.Nspace, rh->input.Ngdelay,
									rh->input.Ngorder, rh->input.Ngperiod, atom->n[0]);
			}
			for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
				molecule = rh->atmos.activemols[nact];
				Ngorder  = (rh->input.accelerate_mols) ? rh->input.Ngorder : 0;
				NgFree(rh, molecule->Ng_nv);
				molecule->Ng_nv = NgInit(rh, molecule->Nv * rh->atmos.Nspace, rh->input.Ngdelay,
										 Ngorder, rh->input.Ngperiod, molecule->nv[0]);
			}
		} else if (dpopsmax < iterLimit) {
			converged = TRUE;
			break;
		}
		niter++;

		/* --- Whether another iteration follows; if a stop is requested
			   none does --                               -------------- */

		more = moreIterations(rh, niter, NmaxIter);
		if (more && StopRequested(rh)) {
			stopped = TRUE;
			more    = FALSE;
		}

		if (rh->input.solve_ne == ITERATION)
			Background(rh, write_analyze_output = TRUE, equilibria_only = FALSE);

		if (rh->atmos.hydrostatic) {
			if (!rh->atmos.atoms[0].active) {
				Error(rh, ERROR_LEVEL_2, routineName, "Can only perform hydrostatic equilibrium"
													  " for hydrogen active");
			}

			/* --- After the last iteration of a run that did not converge
				   (N_MAX_ITER reached, or a stop requested) skip the
				   hydrostatic step. It rescales only hydrogen to the new
				   nHtot, so the other active atoms would be off by the last
				   HSE change in the output (N-16). A converged run never
				   gets here --                                ------------ */

			if (more) Hydrostatic(rh, N_MAX_HSE_ITER, HSE_ITER_LIMIT);
		}
		if (stopped) break;
	}

	/* --- Say whether the populations converged. A run that stops at
		   N_MAX_ITER still writes its output, so make it clear that
		   this output is not a converged solution --   -------------- */

	if (converged) {
		Message(rh, "\n -- Converged in %d iterations: max change"
					" %.4E < ITER_LIMIT = %.4E\n",
				niter, dpopsmax, iterLimit);
	} else if (stopped) {
		Warning(rh, routineName, "Stopped on request (STOP_RH) after %d iterations:"
								 " NOT CONVERGED.\n The output is not a converged solution",
				niter - 1);
		SetExitFlag(rh, RH_EXIT_NOT_CONVERGED);
	} else {
		Warning(rh, routineName, "NOT CONVERGED after %d iterations: max change"
								 " %.4E >= ITER_LIMIT = %.4E.\n The output is not a converged"
								 " solution; increase N_MAX_ITER or see what limits convergence",
				niter - 1, dpopsmax, iterLimit);
		SetExitFlag(rh, RH_EXIT_NOT_CONVERGED);
	}

	for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
		atom = rh->atmos.activeatoms[nact];
		freeMatrix(rh, (void **) atom->Gamma);
		NgFree(rh, atom->Ng_n);
	}
	for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
		molecule = rh->atmos.activemols[nact];
		freeMatrix(rh, (void **) molecule->Gamma);
		NgFree(rh, molecule->Ng_nv);
	}
}
/* ------- end ---------------------------- Iterate.c --------------- */

/* ------- begin -------------------------- solveSpectrum.c --------- */

/* --- The note on the wavelengths solved linear (N-52), for the
	   background scattering line of the sweep (main thread) -- ---- */

static const char *linearNote(const BgScatterStats *bs, char *note, size_t size)
{
	note[0] = '\0';
	if (bs->Nlinear > 0)
		snprintf(note, size, " (%ld solved linear; +%ld -%ld)",
				 bs->Nlinear, bs->Nto, bs->Nback);
	return note;
}

double solveSpectrum(RHState *rh, bool_t eval_operator, bool_t redistribute)
{
	register int nspect, nt;

	int         lambda_max, Nsolved;
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

	getCPU(rh, 3, TIME_START, NULL);

	/* --- First zero the radiative rates --             -------------- */

	zeroRates(rh, redistribute);
	lambda_max = 0;
	dJmax      = 0.0;

	/* --- Count the wavelength solves of this sweep (S-1) --  -------- */

	if (redistribute) {
		Nsolved = 0;
		for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++)
			if (containsPRDline(rh, &rh->spectrum.as[nspect])) Nsolved++;
	} else
		Nsolved = rh->spectrum.Nspect;
	workSweep(rh, eval_operator, redistribute, Nsolved);

	/* --- Zero out J in gas parcel's frame. The PRD sub-iterations
		   re-solve only the wavelengths with a PRD line and keep the
		   contributions of the other wavelengths (N-27) -- ---------- */

	if (rh->spectrum.updateJ &&
		rh->input.PRD_angle_dep == PRD_ANGLE_APPROX &&
		rh->atmos.NPRDactive > 0) {
		beginJgas(rh, redistribute);
	}

	if (rh->input.Nthreads > 1) {

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

		ti        = (threadinfo *) mallocOrAbort(rh->input.Nthreads * sizeof(threadinfo));
		thread_id = (pthread_t *) mallocOrAbort(rh->input.Nthreads * sizeof(pthread_t));

		for (nt = 0; nt < rh->input.Nthreads; nt++) {
			ti[nt].rh            = rh;
			ti[nt].eval_operator = eval_operator;
			ti[nt].redistribute  = redistribute;
			ti[nt].nt            = nt;
			ti[nt].Nthreads      = rh->input.Nthreads;

			if (pthread_create(&thread_id[nt], &rh->input.thread_attr,
							   Formal_pthread, &ti[nt])) {

				/* --- Could not start the thread: do its slice here --- ---- */

				Warning(rh, "solveSpectrum",
						"Could not create a thread, running its share in the caller");
				Formal_pthread(&ti[nt]);
				thread_id[nt] = 0;
			}
		}
		for (nt = 0; nt < rh->input.Nthreads; nt++) {
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

		for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++) {
			if (negTauFound(rh)) break; /* the sweep is void (N-40) */
			if (!redistribute ||
				(redistribute && containsPRDline(rh, &rh->spectrum.as[nspect]))) {
				dJ = Formal(rh, nspect, eval_operator, redistribute);
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

	mergeRates(rh, redistribute);

	if (rh->spectrum.updateJ &&
		rh->input.PRD_angle_dep == PRD_ANGLE_APPROX &&
		rh->atmos.NPRDactive > 0)
		endJgas(rh, redistribute);

	/* --- A wavelength whose total opacity was negative enough to
		   overflow the formal solution was skipped (N-40): stop with a
		   message, unless a solver is evaluating a trial point, which it
		   then rejects --                                -------------- */

	negTauStop(rh);

	Message(rh, " Spectrum max delta J = %6.4E (lambda#: %d)\n",
			dJmax, lambda_max);

	/* --- The background scattering iteration of the sweep (N-51) -- */

	if (bgScatterActive(rh)) {
		BgScatterStats bs;
		char           note[64];

		bgScatterStats(rh, &bs, TRUE);
		if (bs.Niterated > 0) {
			Message(rh, "  Background scattering: %ld wavelengths iterated,"
						" %ld solves%s%s",
					bs.Niterated, bs.Nsweeps, linearNote(&bs, note, sizeof(note)),
					(bs.Nunconverged > 0) ? "" : "\n");
			if (bs.Nunconverged > 0) {
				Message(rh, "; %ld not converged at BACKGROUND_N_MAX_ITER"
							" (estimated error up to %.2E)\n",
						bs.Nunconverged, bs.rmax);
			}
			workBgScatter(rh, bs.Nsweeps);
		}
	}

	getCPU(rh, 3, TIME_POLL,
		   (eval_operator) ? "Spectrum & Operator" : "Solve Spectrum");

	return dJmax;
}
/* ------- end ---------------------------- solveSpectrum.c --------- */

/* ------- begin -------------------------- Formal_pthread.c -------- */

void *Formal_pthread(void *argument)
{
	threadinfo *ti = (threadinfo *) argument;
	RHState    *rh = ti->rh;
	int         nspect;
	double      dJ;

	/* --- Thread wrapper around Formal: works through this thread's
		   strided slice of the wavelength grid --        -------------- */

	ti->dJmax      = 0.0;
	ti->lambda_max = -1;

	for (nspect = ti->nt; nspect < rh->spectrum.Nspect; nspect += ti->Nthreads) {
		if (negTauFound(rh)) break; /* the sweep is void (N-40) */
		if (!ti->redistribute ||
			(ti->redistribute && containsPRDline(rh, &rh->spectrum.as[nspect]))) {
			dJ = Formal(rh, nspect, ti->eval_operator, ti->redistribute);
			if (dJ > ti->dJmax) {
				ti->dJmax      = dJ;
				ti->lambda_max = nspect;
			}
		}
	}
	freeJgasBuffer();

	return (NULL);
}
/* ------- end ---------------------------- Formal_pthread.c -------- */
