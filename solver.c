/* ------- file: -------------------------- solver.c ----------------

	   The map G of docs/SOLVER.md (S-2) and the fixed-point driver
	   (SOLVER = FIXED_POINT).

	   --------------------------                      ----------RH-- */

/* --- G is one MALI iteration as a function of an explicit state
	   vector x (S-2):

	   x = (populations of the active atoms, populations of the active
				molecules, rho of every PRD line, the mean intensity J at
				all wavelengths [and J20 with background polarization]),

	   in that order, each block a copy of the contiguous array RH
	   keeps it in. G_eval(x_in, x_out) copies x_in into those arrays,
	   then does steps 1-4 of SOLVER.md section 2.1:

		 1. collisional and fixed rates into Gamma,
		 2. formal solution with the approximate operator,
		 3. statistical equilibrium for every active atom and molecule,
			without Ng,
		 4. PRD_N_MAX_ITER sub-iterations of the PRD rho, each followed
			by a formal solution of the PRD wavelengths, without Ng and
			without early exit,

	   and copies the arrays into x_out. Everything else the iteration
	   reads is either computed once at setup and never changed (the
	   atmosphere, ne, LTE populations, collision rates, background
	   opacities, line profiles, PRD weights), or recomputed from x in
	   every call before it is used (Gamma, the radiative rates, Jgas,
	   the intensity scratch files). readInput() rejects the keywords
	   that would break this (collisional relaxation, hydrostatic
	   equilibrium, SOLVE_NE = ITERATION, LIMIT_MEMORY, PRD early exit
	   and PRD Ng).

	   The fixed-point driver iterates x <- G(x). With NG_ORDER > 0 it
	   applies RH's Ng acceleration to the populations of each atom
	   (and with NG_MOLECULES of each molecule) after G, as the default
	   solver does: NG_ORDER = 0 gives plain G, the defaults give "Ng
	   on G", the two controls of S-5.
	   --                                              -------------- */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"
#include "background.h"
#include "solver.h"

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum   spectrum;
extern InputData  input;
extern char       messageStr[];

/* --- One block of the state vector: a contiguous array of RH -- --- */

typedef struct {
	char    name[MAX_MESSAGE_LENGTH];
	double *data;
	long    N;
} StateBlock;

static StateBlock *blocks  = NULL;
static int         Nblocks = 0;
static long        Nstate  = 0;

static long   Gcount   = 0;
static double neChange = 0.0, chemChange = 0.0;

/* ------- begin -------------------------- rhoSize.c --------------- */

static long rhoSize(AtomicLine *line)
{
	long Nlamu = (input.PRD_angle_dep == PRD_ANGLE_DEP)
					 ? 2L * atmos.Nrays * line->Nlambda
					 : (long) line->Nlambda;

	return Nlamu * atmos.Nspace;
}
/* ------- end ---------------------------- rhoSize.c --------------- */

/* ------- begin -------------------------- addBlock.c -------------- */

static void addBlock(const char *name, double *data, long N)
{
	blocks = (StateBlock *) realloc(blocks, (Nblocks + 1) * sizeof(StateBlock));
	snprintf(blocks[Nblocks].name, MAX_MESSAGE_LENGTH, "%s", name);
	blocks[Nblocks].data = data;
	blocks[Nblocks].N    = N;
	Nblocks++;
	Nstate += N;
}
/* ------- end ---------------------------- addBlock.c -------------- */

/* ------- begin -------------------------- stateInit.c ------------- */

long stateInit(void)
{
	register int nact, kr;

	char        name[MAX_MESSAGE_LENGTH];
	Atom       *atom;
	Molecule   *molecule;
	AtomicLine *line;

	/* --- Set up the list of blocks of the state vector. Call after
		   initSolution() (populations, rho) and the allocation of J,
		   before the first G_eval. Returns the length of x -- ------- */

	free(blocks);
	blocks  = NULL;
	Nblocks = 0;
	Nstate  = 0;

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];
		sprintf(name, "populations %s", atom->ID);
		addBlock(name, atom->n[0], (long) atom->Nlevel * atmos.Nspace);
	}
	for (nact = 0; nact < atmos.Nactivemol; nact++) {
		molecule = atmos.activemols[nact];
		sprintf(name, "populations %s", molecule->ID);
		addBlock(name, molecule->nv[0], (long) molecule->Nv * atmos.Nspace);
	}
	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];
		for (kr = 0; kr < atom->Nline; kr++) {
			line = atom->line + kr;
			if (line->PRD) {
				sprintf(name, "rho %s %d-%d", atom->ID, line->j, line->i);
				addBlock(name, line->rho_prd[0], rhoSize(line));
			}
		}
	}
	addBlock("J", spectrum.J[0], (long) spectrum.Nspect * atmos.Nspace);
	if (spectrum.J20 != NULL)
		addBlock("J20", spectrum.J20[0], (long) spectrum.Nspect * atmos.Nspace);

	if (input.solve_ne == ITERATION) addBlock("ne", atmos.ne, atmos.Nspace);

	return Nstate;
}
/* ------- end ---------------------------- stateInit.c ------------- */

/* ------- begin -------------------------- stateBlocks.c ----------- */

/* --- The blocks of x, for reports per block (gchecks.c) -- ------- */

int stateNblocks(void)
{
	return Nblocks;
}

const char *stateBlockName(int n)
{
	return blocks[n].name;
}

long stateBlockSize(int n)
{
	return blocks[n].N;
}
/* ------- end ---------------------------- stateBlocks.c ----------- */

/* ------- begin -------------------------- stateGather.c ----------- */

void stateGather(double *x)
{
	register int n;

	for (n = 0; n < Nblocks; n++) {
		memcpy(x, blocks[n].data, blocks[n].N * sizeof(double));
		x += blocks[n].N;
	}
}
/* ------- end ---------------------------- stateGather.c ----------- */

/* ------- begin -------------------------- stateScatter.c ---------- */

void stateScatter(const double *x)
{
	register int n;

	for (n = 0; n < Nblocks; n++) {
		memcpy(blocks[n].data, x, blocks[n].N * sizeof(double));
		x += blocks[n].N;
	}
}
/* ------- end ---------------------------- stateScatter.c ---------- */

/* ------- begin -------------------------- relChange.c ------------- */

static double relChange(const double *old, const double *new, long N)
{
	register long k;

	double dmax = 0.0;

	/* --- MaxChange's definition: max |new - old| / |new|, zeros of
		   new skipped --                                -------------- */

	for (k = 0; k < N; k++) {
		if (new[k]) dmax = MAX(dmax, fabs((new[k] - old[k]) / new[k]));
	}
	return dmax;
}
/* ------- end ---------------------------- relChange.c ------------- */

/* ------- begin -------------------------- prdSubIterations.c ------ */

static double prdSubIterations(int Nsub)
{
	register int nact, kr;

	bool_t      eval_operator, redistribute;
	int         nsub;
	long        N;
	double     *old = NULL, drho, drhomax = 0.0;
	Atom       *atom;
	AtomicLine *line;

	/* --- Step 4: a fixed number of sub-iterations of rho with the
		   populations fixed, no Ng, no early exit (N-4). Returns the
		   largest relative change of rho in the last sub-iteration,
		   for the log --                                -------------- */

	for (nsub = 1; nsub <= Nsub; nsub++) {
		drhomax = 0.0;
		for (nact = 0; nact < atmos.Nactiveatom; nact++) {
			atom = atmos.activeatoms[nact];
			for (kr = 0; kr < atom->Nline; kr++) {
				line = atom->line + kr;
				if (!line->PRD) continue;

				N   = rhoSize(line);
				old = (double *) realloc(old, N * sizeof(double));
				memcpy(old, line->rho_prd[0], N * sizeof(double));

				if (input.PRD_angle_dep == PRD_ANGLE_DEP)
					PRDAngleScatter(line);
				else
					PRDScatter(line);

				drho = relChange(old, line->rho_prd[0], N);
				sprintf(messageStr, "  PRD: sub-iteration %d, atom %s, line %d,"
									" delta = %6.4E\n",
						nsub, atom->ID, kr, drho);
				Error(MESSAGE, NULL, messageStr);
				drhomax = MAX(drho, drhomax);
			}
		}
		solveSpectrum(eval_operator = FALSE, redistribute = TRUE);
	}
	free(old);

	return drhomax;
}
/* ------- end ---------------------------- prdSubIterations.c ------ */

/* ------- begin -------------------------- G_eval.c ---------------- */

void G_eval(const double *x_in, double *x_out)
{
	register int nact;

	bool_t eval_operator, redistribute;

	/* --- x_in and x_out may be the same array --     -------------- */

	stateScatter(x_in);

	/* --- 0. With ne in the state (S-2b): everything the background
		   derives from ne and the populations, from those of x_in:
		   LTE populations, collision and fixed rates, chemical
		   equilibrium, background opacities --        -------------- */

	if (input.solve_ne == ITERATION) BackgroundOfState(FALSE);

	/* --- 1. Collisional and fixed rates into Gamma -- -------------- */

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		initGammaAtom(atmos.activeatoms[nact], 0);
		zeroGammaThreads(atmos.activeatoms[nact]);
	}
	for (nact = 0; nact < atmos.Nactivemol; nact++) {
		initGammaMolecule(atmos.activemols[nact]);
		zeroGammaThreadsMolecule(atmos.activemols[nact]);
	}
	/* --- 2. Formal solution with the approximate operator -- ------- */

	solveSpectrum(eval_operator = TRUE, redistribute = FALSE);
	for (nact = 0; nact < atmos.Nactiveatom; nact++)
		mergeGammaThreads(atmos.activeatoms[nact]);
	for (nact = 0; nact < atmos.Nactivemol; nact++)
		mergeGammaThreadsMolecule(atmos.activemols[nact]);

	/* --- 3. Statistical equilibrium, no acceleration -- ------------ */

	for (nact = 0; nact < atmos.Nactiveatom; nact++)
		statEquil(atmos.activeatoms[nact], input.isum);
	for (nact = 0; nact < atmos.Nactivemol; nact++)
		statEquilMolecule(atmos.activemols[nact], 0);

	/* --- 4. PRD sub-iterations --                     -------------- */

	if (atmos.NPRDactive > 0) prdSubIterations(input.PRD_NmaxIter);

	/* --- 5. With ne in the state: ne from charge conservation with
		   the new populations, Newton from the ne of x_in -- ------- */

	if (input.solve_ne == ITERATION) {
		neChange   = Solve_ne_OfState(atmos.ne, N_NE_ITER_OF_STATE);
		chemChange = chemLastChangeReset();
	}
	stateGather(x_out);
	Gcount++;
}
/* ------- end ---------------------------- G_eval.c ---------------- */

/* ------- begin -------------------------- G_count.c --------------- */

long G_count(void)
{
	return Gcount;
}
/* ------- end ---------------------------- G_count.c --------------- */

/* ------- begin -------------------------- iterateFixedPoint.c ----- */

void iterateFixedPoint(int NmaxIter, double iterLimit)
{
	const char   routineName[] = "iterateFixedPoint";
	register int niter, nact;

	bool_t    accel, quiet, converged = FALSE, stopped = FALSE;
	int       Ngorder;
	double   *x, dpops, dpopsmax = 0.0;
	Atom     *atom;
	Molecule *molecule;

	stateInit();
	x = (double *) malloc(Nstate * sizeof(double));

	sprintf(messageStr,
			"\n -- SOLVER = FIXED_POINT: x <- G(x), %ld unknowns;"
			" PRD sub-iterations per G: %d;\n    Ng on G: order %d,"
			" delay %d, period %d%s\n",
			Nstate, (atmos.NPRDactive > 0) ? input.PRD_NmaxIter : 0,
			input.Ngorder, input.Ngdelay, input.Ngperiod,
			(input.Ngorder > 0) ? "" : " (none: plain G)");
	Error(MESSAGE, routineName, messageStr);

	/* --- Ng histories on the populations, as in the default solver.
		   With NG_ORDER = 0 they only keep the previous iterate for
		   the change measure --                         -------------- */

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
	stateGather(x);

	niter = 1;
	while (niter <= NmaxIter) {
		if (StopRequested()) {
			stopped = TRUE;
			break;
		}
		getCPU(2, TIME_START, NULL);

		G_eval(x, x);

		/* --- Change measure (and optional Ng) per atom, printed in
			   the format of the default solver --        -------------- */

		sprintf(messageStr, "\n -- Iteration %3d\n", niter);
		Error(MESSAGE, routineName, messageStr);
		dpopsmax = 0.0;
		for (nact = 0; nact < atmos.Nactiveatom; nact++) {
			atom  = atmos.activeatoms[nact];
			accel = Accelerate(atom->Ng_n, atom->n[0]);
			sprintf(messageStr, " %s,", atom->ID);
			dpops = MaxChange(atom->Ng_n, messageStr, quiet = FALSE);
			Error(MESSAGE, NULL, (accel) ? " (accelerated)\n" : "\n");
			dpopsmax = MAX(dpops, dpopsmax);
		}
		for (nact = 0; nact < atmos.Nactivemol; nact++) {
			molecule = atmos.activemols[nact];
			accel    = Accelerate(molecule->Ng_nv, molecule->nv[0]);
			sprintf(messageStr, " %s ,", molecule->ID);
			dpops = MaxChange(molecule->Ng_nv, messageStr, quiet = FALSE);
			Error(MESSAGE, NULL, (accel) ? " (accelerated)\n" : "\n");
			dpopsmax = MAX(dpops, dpopsmax);
		}
		if (input.solve_ne == ITERATION) {
			sprintf(messageStr, " ne: last Newton step %.2E, chemical"
								" equilibrium: last Newton step %.2E\n",
					neChange, chemChange);
			Error(MESSAGE, NULL, messageStr);
		}
		/* --- Ng changes the populations in place: take them back
			   into x --                                     -------------- */

		stateGather(x);

		sprintf(messageStr, "Total Iteration %3d", niter);
		getCPU(2, TIME_POLL, messageStr);

		if (dpopsmax < iterLimit) {
			converged = TRUE;
			break;
		}
		niter++;
	}

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
	sprintf(messageStr, " -- G evaluations: %ld\n", Gcount);
	Error(MESSAGE, routineName, messageStr);

	/* --- With ne in the state: the background of the final state, for
		   the final formal solution and the output --   -------------- */

	if (input.solve_ne == ITERATION) BackgroundOfState(TRUE);

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
	free(x);
}
/* ------- end ---------------------------- iterateFixedPoint.c ----- */
