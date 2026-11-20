/* ------- file: -------------------------- solver.c ----------------

	   The map G of docs/SOLVER.md (S-2) and the solvers that iterate
	   it: the fixed-point driver (SOLVER = FIXED_POINT) and Anderson
	   acceleration (SOLVER = ANDERSON, S-3).

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

#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "inputs.h"
#include "negtau.h"
#include "error.h"
#include "statistics.h"
#include "background.h"
#include "anderson.h"
#include "newton.h"
#include "mpe.h"
#include "solver.h"

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum   spectrum;
extern InputData  input;
extern char       messageStr[];

/* --- One block of the state vector: a contiguous array of RH -- --- */

enum BlockKind { BLOCK_ATOM,
				 BLOCK_MOLECULE,
				 BLOCK_RHO,
				 BLOCK_J,
				 BLOCK_J20,
				 BLOCK_NE };

typedef struct {
	char           name[MAX_MESSAGE_LENGTH];
	double        *data;
	long           N;
	enum BlockKind kind;
	void          *owner; /* the Atom or Molecule of a population or rho block */
} StateBlock;

static StateBlock *blocks  = NULL;
static int         Nblocks = 0;
static long        Nstate  = 0;

/* --- Nblocks + 1 offsets of the blocks into x, set by addBlock -- -- */

static long *blockStart = NULL;

/* --- Largest u with exp(u) finite: log(DBL_MAX) rounded down -- */

#define LOG_DBL_MAX 709.0

static long   Gcount   = 0;
static double neChange = 0.0, chemChange = 0.0;

/* ------- begin -------------------------- addBlock.c -------------- */

static void addBlock(const char *name, double *data, long N,
					 enum BlockKind kind, void *owner)
{
	blocks = (StateBlock *) realloc(blocks, (Nblocks + 1) * sizeof(StateBlock));
	snprintf(blocks[Nblocks].name, MAX_MESSAGE_LENGTH, "%s", name);
	blocks[Nblocks].data    = data;
	blocks[Nblocks].N       = N;
	blocks[Nblocks].kind    = kind;
	blocks[Nblocks].owner   = owner;
	blockStart              = (long *) realloc(blockStart, (Nblocks + 2) * sizeof(long));
	blockStart[0]           = 0;
	blockStart[Nblocks + 1] = blockStart[Nblocks] + N;
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
		addBlock(name, atom->n[0], (long) atom->Nlevel * atmos.Nspace,
				 BLOCK_ATOM, atom);
	}
	for (nact = 0; nact < atmos.Nactivemol; nact++) {
		molecule = atmos.activemols[nact];
		sprintf(name, "populations %s", molecule->ID);
		addBlock(name, molecule->nv[0], (long) molecule->Nv * atmos.Nspace,
				 BLOCK_MOLECULE, molecule);
	}
	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];
		for (kr = 0; kr < atom->Nline; kr++) {
			line = atom->line + kr;
			if (line->PRD) {
				sprintf(name, "rho %s %d-%d", atom->ID, line->j, line->i);
				addBlock(name, line->rho_prd[0], PRDRhoSize(line), BLOCK_RHO, atom);
			}
		}
	}
	addBlock("J", spectrum.J[0], (long) spectrum.Nspect * atmos.Nspace, BLOCK_J,
			 NULL);
	if (spectrum.J20 != NULL)
		addBlock("J20", spectrum.J20[0], (long) spectrum.Nspect * atmos.Nspace,
				 BLOCK_J20, NULL);

	if (input.solve_ne == ITERATION)
		addBlock("ne", atmos.ne, atmos.Nspace, BLOCK_NE, NULL);

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

/* ------- begin -------------------------- residualHistory.c ------- */

/* --- The residuals of the iterations (largest relChange of the
	   populations) of one run of SOLVER = FIXED_POINT, the last
	   RES_HIST_MAX of them, for the estimate of the contraction at the
	   end of the run (N-3). Only there the residuals follow the
	   contraction of G: for ANDERSON and NEWTON_KRYLOV the same estimate
	   was measured to understate the error by up to 20 times, and is
	   not made -- ----------------------------------------------- */

#define RES_HIST_MAX 11

static double resHist[RES_HIST_MAX];
static int    NresHist = 0;

static void resHistReset(void)
{
	NresHist = 0;
}

static void resHistAdd(double r)
{
	if (NresHist == RES_HIST_MAX) {
		memmove(resHist, resHist + 1, (RES_HIST_MAX - 1) * sizeof(double));
		NresHist--;
	}
	resHist[NresHist++] = r;
}

static void printContraction(const char *routineName)
{
	int    m;
	double lambda, rlast;

	/* --- The geometric mean of the ratio of successive residuals
		   over the last m (at most 10) iterations. Without the
		   evidence (fewer than four iterations, no decrease) say so.
		   The residual bounds the error through the contraction:
		   error ~ r / (1 - lambda) --                   -------------- */

	m = NresHist - 1;
	if (m < 2 || !(resHist[0] > 0.0) || !(resHist[NresHist - 1] >= 0.0) ||
		!isfinite(resHist[0]) || !isfinite(resHist[NresHist - 1])) {
		sprintf(messageStr, " -- Estimated contraction: no estimate (fewer than 4"
							" iterations, or no residual to compare)\n");
		Error(MESSAGE, routineName, messageStr);
		return;
	}
	rlast  = resHist[NresHist - 1];
	lambda = pow(rlast / resHist[0], 1.0 / m);
	if (!(lambda < 1.0)) {
		sprintf(messageStr, " -- Estimated contraction: no estimate (the residual"
							" did not decrease over the last %d iterations)\n",
				m);
		Error(MESSAGE, routineName, messageStr);
		return;
	}
	sprintf(messageStr, " -- Estimated contraction %.3f per iteration (last %d"
						" iterations); the error of the populations may be up to"
						" about %.1E (residual / (1 - contraction))\n",
			lambda, m, rlast / (1.0 - lambda));
	Error(MESSAGE, routineName, messageStr);
}
/* ------- end ---------------------------- residualHistory.c ------- */

/* ------- begin -------------------------- printRhoNe.c ------------ */

static long rhoNeSize(void)
{
	register int n;

	long N = 0;

	for (n = 0; n < Nblocks; n++)
		if (blocks[n].kind == BLOCK_RHO || blocks[n].kind == BLOCK_NE) N += blocks[n].N;
	return N;
}

static void rhoNeSave(const double *x, double *buf)
{
	register int n;

	/* --- The rho and ne blocks of x, one after the other, for a
		   solver that overwrites x in place --         -------------- */

	for (n = 0; n < Nblocks; n++) {
		if (blocks[n].kind != BLOCK_RHO && blocks[n].kind != BLOCK_NE) continue;
		memcpy(buf, x + blockStart[n], blocks[n].N * sizeof(double));
		buf += blocks[n].N;
	}
}

static void printRhoNe(const double *x, const double *saved, const double *gx)
{
	register int n, nact;

	const double *old;
	double        dpops;

	/* --- Reporting only (N-3): the residual of the PRD rho blocks, one
		   per atom with PRD lines, and of ne when it is in the state,
		   relChange's measure on (x, G(x)). They are not part of the
		   ITER_LIMIT test. With saved (rhoNeSave) the old values are
		   taken from there instead of x -- ------------- */

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		Atom  *atom = atmos.activeatoms[nact];
		bool_t any  = FALSE;

		dpops = 0.0;
		old   = saved;
		for (n = 0; n < Nblocks; n++) {
			if (blocks[n].kind != BLOCK_RHO && blocks[n].kind != BLOCK_NE) continue;
			if (blocks[n].kind == BLOCK_RHO && blocks[n].owner == (void *) atom) {
				dpops = MAX(dpops, relChange((saved) ? old : x + blockStart[n],
											 gx + blockStart[n], blocks[n].N));
				any   = TRUE;
			}
			if (saved) old += blocks[n].N;
		}
		if (any) {
			sprintf(messageStr, " rho %s: delta = %6.4E\n", atom->ID, dpops);
			Error(MESSAGE, NULL, messageStr);
		}
	}
	old = saved;
	for (n = 0; n < Nblocks; n++) {
		if (blocks[n].kind != BLOCK_RHO && blocks[n].kind != BLOCK_NE) continue;
		if (blocks[n].kind == BLOCK_NE) {
			dpops = relChange((saved) ? old : x + blockStart[n], gx + blockStart[n],
							  blocks[n].N);
			sprintf(messageStr, " ne: delta = %6.4E\n", dpops);
			Error(MESSAGE, NULL, messageStr);
		}
		if (saved) old += blocks[n].N;
	}
}
/* ------- end ---------------------------- printRhoNe.c ------------ */

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
		   populations fixed, no Ng, no early exit (N-4), each update
		   damped by PRD_DAMPING (N-38). Returns the largest residual of
		   rho in the last sub-iteration, the undamped change (the step
		   taken is PRD_DAMPING times it), for the log -- ---------- */

	for (nsub = 1; nsub <= Nsub; nsub++) {
		drhomax = 0.0;
		for (nact = 0; nact < atmos.Nactiveatom; nact++) {
			atom = atmos.activeatoms[nact];
			for (kr = 0; kr < atom->Nline; kr++) {
				line = atom->line + kr;
				if (!line->PRD) continue;

				N   = PRDRhoSize(line);
				old = (double *) realloc(old, N * sizeof(double));
				memcpy(old, line->rho_prd[0], N * sizeof(double));

				if (input.PRD_angle_dep == PRD_ANGLE_DEP)
					PRDAngleScatter(line);
				else
					PRDScatter(line);
				drho = PRDDamp(line, old, N);
				sprintf(messageStr, "  PRD: sub-iteration %d, atom %s, line %d,"
									" delta = %6.4E\n",
						nsub, atom->ID, kr, drho);
				Error(MESSAGE, NULL, messageStr);
				drhomax = MAX(drho, drhomax);

				/* --- rho <= 0 or not finite: stops the run, or at a
					   trial point invalidates it (N-4) --  -------------- */

				if (!PRDRhoCheck(atom, kr)) break; /* G_eval gives up (N-4) */
			}
			if (negTauFound()) break;
		}
		if (negTauFound()) break;
		solveSpectrum(eval_operator = FALSE, redistribute = TRUE);
		if (negTauFound()) break; /* G_eval gives up (N-40) */
	}
	free(old);

	return drhomax;
}
/* ------- end ---------------------------- prdSubIterations.c ------ */

/* ------- begin -------------------------- G_eval.c ---------------- */

bool_t G_eval(const double *x_in, double *x_out)
{
	register int nact;

	bool_t eval_operator, redistribute;

	/* --- x_in and x_out may be the same array. Returns FALSE, with
		   x_out left as it was, when a formal solution met a total
		   opacity negative enough to overflow (N-40). That can only
		   happen at a trial point (negTauTrial(TRUE)): anywhere else
		   solveSpectrum stops the run. Also FALSE, in the same way,
		   when statEquil found negative or non-finite populations out
		   of its LU fallback (N-8; the state is then that of x_in
		   again), or a PRD sub-iteration gave rho <= 0 or not finite
		   (N-4) --                                     -------------- */

	stateScatter(x_in);

	/* --- A trial point (an extrapolated Anderson point, a Newton
		   probe) brings its own rho: not positive and finite is invalid,
		   before the formal solution makes a negative emissivity of it
		   (N-4). The guard of the Anderson trial leaves zero and the
		   values where G's point is negative alone, and the probes of
		   Newton have none --                           -------------- */

	if (negTauInTrial() && atmos.NPRDactive > 0) {
		for (nact = 0; nact < atmos.Nactiveatom; nact++) {
			Atom *atom = atmos.activeatoms[nact];

			for (int kr = 0; kr < atom->Nline; kr++) {
				if (atom->line[kr].PRD && !PRDRhoCheck(atom, kr)) {
					Gcount++;
					return FALSE;
				}
			}
		}
	}
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
	if (negTauFound()) {
		Gcount++;
		return FALSE;
	}
	for (nact = 0; nact < atmos.Nactiveatom; nact++)
		mergeGammaThreads(atmos.activeatoms[nact]);
	for (nact = 0; nact < atmos.Nactivemol; nact++)
		mergeGammaThreadsMolecule(atmos.activemols[nact]);

	/* --- 3. Statistical equilibrium, no acceleration -- ------------ */

	for (nact = 0; nact < atmos.Nactiveatom && !negTauFound(); nact++)
		statEquil(atmos.activeatoms[nact], input.isum);
	for (nact = 0; nact < atmos.Nactivemol && !negTauFound(); nact++)
		statEquilMolecule(atmos.activemols[nact], 0);
	if (negTauFound()) { /* bad LU populations (N-8) */

		/* --- statEquil has overwritten the populations: put the state
			   of x_in back, as it is when the formal solution gives up -- */

		stateScatter(x_in);
		Gcount++;
		return FALSE;
	}

	/* --- 4. PRD sub-iterations --                     -------------- */

	if (atmos.NPRDactive > 0) prdSubIterations(input.PRD_NmaxIter);
	if (negTauFound()) {
		stateScatter(x_in); /* rho of the sub-iterations goes too */
		Gcount++;
		return FALSE;
	}

	/* --- 5. With ne in the state: ne from charge conservation with
		   the new populations, Newton from the ne of x_in -- ------- */

	if (input.solve_ne == ITERATION) {
		neChange   = Solve_ne_OfState(atmos.ne, N_NE_ITER_OF_STATE);
		chemChange = chemLastChangeReset();
	}
	stateGather(x_out);
	Gcount++;

	return TRUE;
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
	double   *x, *saved, dpops, dpopsmax = 0.0;
	Atom     *atom;
	Molecule *molecule;

	stateInit();
	x     = (double *) malloc(Nstate * sizeof(double));
	saved = (double *) malloc(rhoNeSize() * sizeof(double) + 1);
	resHistReset();

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

		rhoNeSave(x, saved);
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
			Error(MESSAGE, NULL, NgStatus(atom->Ng_n, accel));
			dpopsmax = MAX(dpops, dpopsmax);
		}
		for (nact = 0; nact < atmos.Nactivemol; nact++) {
			molecule = atmos.activemols[nact];
			accel    = Accelerate(molecule->Ng_nv, molecule->nv[0]);
			sprintf(messageStr, " %s ,", molecule->ID);
			dpops = MaxChange(molecule->Ng_nv, messageStr, quiet = FALSE);
			Error(MESSAGE, NULL, NgStatus(molecule->Ng_nv, accel));
			dpopsmax = MAX(dpops, dpopsmax);
		}
		if (input.solve_ne == ITERATION) {
			sprintf(messageStr, " ne: last Newton step %.2E, chemical"
								" equilibrium: last Newton step %.2E\n",
					neChange, chemChange);
			Error(MESSAGE, NULL, messageStr);
		}
		printRhoNe(NULL, saved, x);
		if (niter > 1) resHistAdd(dpopsmax); /* not the jump from the start */
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
	/* --- With Ng the residuals are those of accelerated iterates:
		   no estimate, as for ANDERSON and NEWTON_KRYLOV --  ------- */

	if (input.Ngorder == 0) printContraction(routineName);
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
	free(saved);
}
/* ------- end ---------------------------- iterateFixedPoint.c ----- */

/* ------- begin -------------------------- toVars.c ---------------- */

/* --- SOLVER = ANDERSON (S-3 in docs/SOLVER.md) and NEWTON_KRYLOV
	   (S-4) work on scaled variables u = s(x), one block of x at a
	   time:

		 LOG       u = log x (populations, rho, J, ne);
		 RELATIVE  u = x / scale, scale fixed at the start (the
				   element's total population at that depth, the
				   molecule's total in x_0, 1 for rho, J of x_0 for J,
				   n_H,tot for ne);

	   J20 (background polarization, signed) is J20 / J in both. An
	   entry that has no value in these variables (x <= 0 in LOG) is
	   invalid: it takes no part in the least squares (the Newton
	   step) and takes G's value.

	   The variables need not cover all of x. Anderson can extrapolate
	   only some blocks (ANDERSON_EXTRAPOLATE, S-11); the others are
	   carried: they take G's value, like an invalid entry. The packed
	   blocks are listed in packBlock[], in the order of x, and u holds
	   them one after the other (packStart[]). blockStart[] indexes the
	   full x. With every block packed (ALL, and Newton-Krylov) the two
	   tables are the same. --                         -------------- */

static int                  *packBlock  = NULL; /* block of each packed block */
static int                  *carryBlock = NULL; /* the other blocks: carried */
static int                   Ncarry     = 0;
static long                 *packStart  = NULL; /* Npack + 1 offsets into u */
static int                   Npack      = 0;
static long                  Npacked    = 0;    /* length of u */
static long                  offJ       = -1;   /* offset of J in x */
static double               *varScale   = NULL; /* RELATIVE: fixed scale per entry of u */
static enum solver_variables varKind    = VARIABLES_LOG;

static void toVars(const double *x, double *u, unsigned char *valid)
{
	register long i, j;
	register int  p, n;

	long   off;
	double J;

	/* --- x is the full state, u and valid the packed variables -- */

	for (p = 0; p < Npack; p++) {
		n   = packBlock[p];
		off = blockStart[n];
		for (i = off, j = packStart[p]; i < off + blocks[n].N; i++, j++) {
			if (blocks[n].kind == BLOCK_J20) {
				J        = x[offJ + (i - off)];
				valid[j] = (J > 0.0 && isfinite(x[i]) && isfinite(J));
				u[j]     = (valid[j]) ? x[i] / J : 0.0;
			} else if (varKind == VARIABLES_LOG) {
				valid[j] = (x[i] > 0.0 && isfinite(x[i]));
				u[j]     = (valid[j]) ? log(x[i]) : 0.0;
			} else {
				valid[j] = (varScale[j] > 0.0 && isfinite(x[i]));
				u[j]     = (valid[j]) ? x[i] / varScale[j] : 0.0;
			}
		}
	}
}
/* ------- end ---------------------------- toVars.c ---------------- */

/* ------- begin -------------------------- fromVars.c -------------- */

static void fromVars(const double *u, const double *gx,
					 const unsigned char *mask, double *x)
{
	register long i, j;
	register int  p, n;

	long off;

	/* --- x = s^-1(u) where the entry took part, G's value (copied,
		   not transformed) where it did not and in the carried blocks.
		   u and mask are packed, gx and x full. J comes before J20,
		   and J20 is packed only when J is (ALL) -- ---------------- */

	for (p = 0; p < Ncarry; p++) {
		n = carryBlock[p];
		memcpy(x + blockStart[n], gx + blockStart[n], blocks[n].N * sizeof(double));
	}
	for (p = 0; p < Npack; p++) {
		n   = packBlock[p];
		off = blockStart[n];
		for (i = off, j = packStart[p]; i < off + blocks[n].N; i++, j++) {
			if (!mask[j])
				x[i] = gx[i];
			else if (blocks[n].kind == BLOCK_J20)
				x[i] = u[j] * x[offJ + (i - off)];
			else if (varKind == VARIABLES_LOG)
				x[i] = exp(u[j]);
			else
				x[i] = u[j] * varScale[j];
		}
	}
}
/* ------- end ---------------------------- fromVars.c -------------- */

/* ------- begin -------------------------- setVarScale.c ----------- */

static void setVarScale(const double *x0)
{
	register long i;
	register int  p, n, l;

	long      off, poff, k, Nspace = atmos.Nspace;
	double    total;
	Atom     *atom;
	Molecule *molecule;

	/* --- x0 is the full state, the scales are packed -- -------- */

	for (p = 0; p < Npack; p++) {
		n    = packBlock[p];
		off  = blockStart[n];
		poff = packStart[p];
		switch (blocks[n].kind) {
		case BLOCK_ATOM:
			atom = (Atom *) blocks[n].owner;
			for (i = 0; i < blocks[n].N; i++)
				varScale[poff + i] = atom->ntotal[i % Nspace];
			break;
		case BLOCK_MOLECULE:
			molecule = (Molecule *) blocks[n].owner;
			for (k = 0; k < Nspace; k++) {
				total = 0.0;
				for (l = 0; l < molecule->Nv; l++) total += x0[off + l * Nspace + k];
				for (l = 0; l < molecule->Nv; l++) varScale[poff + l * Nspace + k] = total;
			}
			break;
		case BLOCK_J:
			for (i = 0; i < blocks[n].N; i++) varScale[poff + i] = x0[off + i];
			break;
		case BLOCK_NE:
			for (i = 0; i < blocks[n].N; i++) varScale[poff + i] = atmos.nHtot[i];
			break;
		default: /* rho, J20 (J20 / J is used instead) */
			for (i = 0; i < blocks[n].N; i++) varScale[poff + i] = 1.0;
			break;
		}
	}
}
/* ------- end ---------------------------- setVarScale.c ----------- */

/* ------- begin -------------------------- setupVariables.c -------- */

static bool_t extrapolated(enum BlockKind kind, enum anderson_extrapolate which)
{
	switch (which) {
	case EXTRAPOLATE_POPULATIONS:
		return (kind == BLOCK_ATOM || kind == BLOCK_MOLECULE || kind == BLOCK_NE);
	case EXTRAPOLATE_POPULATIONS_RHO:
		return (kind == BLOCK_ATOM || kind == BLOCK_MOLECULE || kind == BLOCK_NE ||
				kind == BLOCK_RHO);
	default:
		return TRUE;
	}
}

static double *setupVariables(enum solver_variables     kind,
							  enum anderson_extrapolate which, double *x)
{
	register int n;

	long    i;
	double *w;

	/* --- After stateInit(): the block offsets in x, the packed
		   blocks (which: those the variables cover), the weights of
		   the packed blocks (every block counts as its mean square:
		   1/sqrt(N)), x_0 in x, and the variables. Returns the
		   weights, packed --                          -------------- */

	offJ = -1;
	for (n = 0; n < Nblocks; n++)
		if (blocks[n].kind == BLOCK_J) offJ = blockStart[n];
	packBlock    = (int *) realloc(packBlock, Nblocks * sizeof(int));
	carryBlock   = (int *) realloc(carryBlock, Nblocks * sizeof(int));
	packStart    = (long *) realloc(packStart, (Nblocks + 1) * sizeof(long));
	Npack        = 0;
	Ncarry       = 0;
	packStart[0] = 0;
	for (n = 0; n < Nblocks; n++) {
		if (!extrapolated(blocks[n].kind, which)) {
			carryBlock[Ncarry++] = n;
			continue;
		}
		packBlock[Npack]     = n;
		packStart[Npack + 1] = packStart[Npack] + blocks[n].N;
		Npack++;
	}
	Npacked = packStart[Npack];

	w = (double *) malloc(MAX(Npacked, 1) * sizeof(double));
	for (n = 0; n < Npack; n++)
		for (i = packStart[n]; i < packStart[n + 1]; i++)
			w[i] = 1.0 / sqrt((double) blocks[packBlock[n]].N);

	varKind = kind;
	stateGather(x);
	if (kind == VARIABLES_RELATIVE) {
		varScale = (double *) realloc(varScale, MAX(Npacked, 1) * sizeof(double));
		setVarScale(x);
	}
	return w;
}
/* ------- end ---------------------------- setupVariables.c -------- */

/* ------- begin -------------------------- guardTrial.c ------------ */

enum GuardResult { GUARD_OK,
				   GUARD_NONFINITE,
				   GUARD_NEGATIVE,
				   GUARD_NRESULT };

static const char *guardName[GUARD_NRESULT] = {"ok", "not finite", "negative"};

static enum GuardResult guardTrial(const double *xt, const double *gx,
								   const double *xb)
{
	long i;
	int  n;

	/* --- The extrapolated point must be finite and must not be
		   negative where G's point is not (and with RELATIVE not zero
		   where it is positive; in LOG a zero is an underflow of exp,
		   not a sign change). Negative values G itself makes (the LU
		   fallback of statEquil, N-31) are left alone.

		   Population inversions are not checked. Converged solutions
		   have inverted continua and lines (S-13), and a test on line
		   inversions did not keep the total opacity from going
		   negative, which is what breaks the formal solver (SOLVER.md
		   trap 5); that is detected where it happens, in the formal
		   solution (N-40, merged with this change). The continuum
		   test went in S-13, the line test in S-14.

		   With a base point xb (Newton: the trial is xb + lambda d), a
		   value only counts as new where neither G's point nor xb has
		   it, so that a small enough lambda always passes; Anderson
		   passes NULL. --                               -------------- */

	for (i = 0; i < Nstate; i++)
		if (!isfinite(xt[i])) return GUARD_NONFINITE;

	for (n = 0; n < Nblocks; n++) {
		long off = blockStart[n];

		if (blocks[n].kind != BLOCK_J20) {
			for (i = off; i < off + blocks[n].N; i++)
				if ((xt[i] < 0.0 && !(gx[i] < 0.0) && !(xb && xb[i] < 0.0)) ||
					(varKind == VARIABLES_RELATIVE && xt[i] == 0.0 &&
					 gx[i] > 0.0 && !(xb && !(xb[i] > 0.0))))
					return GUARD_NEGATIVE;
		}
	}
	return GUARD_OK;
}
/* ------- end ---------------------------- guardTrial.c ------------ */

/* ------- begin -------------------------- logSolve.c -------------- */

static void logSolve(int Nsolve, bool_t per_block, const AndersonSolve *info,
					 double fnorm)
{
	int b;

	/* --- Columns used, condition number and columns left out, every
		   iteration that solved (also when all were left out) -- --- */

	for (b = 0; b < Nsolve; b++) {
		sprintf(messageStr,
				"  Anderson%s%s: %d column%s, cond %.2E, %d dropped;"
				" |f|_W %.3E\n",
				(per_block) ? " " : "", (per_block) ? blocks[packBlock[b]].name : "",
				info[b].Ncols, (info[b].Ncols == 1) ? "" : "s", info[b].cond,
				info[b].Ndrop, fnorm);
		Error(MESSAGE, NULL, messageStr);
	}
}
/* ------- end ---------------------------- logSolve.c -------------- */

/* ------- begin -------------------------- printChange.c ----------- */

static double printChange(const char *routineName, int niter, const double *x,
						  const double *gx)
{
	register int n;

	double dpops, dpopsmax = 0.0;

	/* --- Change of an evaluation per atom and molecule, MaxChange's
		   definition on (x, G(x)), in the format of the other
		   solvers. Returns the largest --            -------------- */

	sprintf(messageStr, "\n -- Iteration %3d\n", niter);
	Error(MESSAGE, routineName, messageStr);
	for (n = 0; n < Nblocks; n++) {
		if (blocks[n].kind != BLOCK_ATOM && blocks[n].kind != BLOCK_MOLECULE)
			continue;
		dpops = relChange(x + blockStart[n], gx + blockStart[n], blocks[n].N);
		sprintf(messageStr,
				(blocks[n].kind == BLOCK_ATOM) ? " %s, delta = %6.4E\n"
											   : " %s , delta = %6.4E\n",
				(blocks[n].kind == BLOCK_ATOM) ? ((Atom *) blocks[n].owner)->ID
											   : ((Molecule *) blocks[n].owner)->ID,
				dpops);
		Error(MESSAGE, NULL, messageStr);
		dpopsmax = MAX(dpops, dpopsmax);
	}
	if (input.solve_ne == ITERATION) {
		sprintf(messageStr, " ne: last Newton step %.2E, chemical"
							" equilibrium: last Newton step %.2E\n",
				neChange, chemChange);
		Error(MESSAGE, NULL, messageStr);
	}
	printRhoNe(x, NULL, gx);
	if (input.errorLimit > 0.0) {
		for (n = 0; n < Nblocks; n++) {
			if (blocks[n].kind != BLOCK_J) continue;
			sprintf(messageStr, " J: delta = %6.4E\n",
					relChange(x + blockStart[n], gx + blockStart[n], blocks[n].N));
			Error(MESSAGE, NULL, messageStr);
		}
	}
	return dpopsmax;
}
/* ------- end ---------------------------- printChange.c ----------- */

/* ------- begin -------------------------- errorTest.c ------------- */

/* --- ERROR_LIMIT (N-3): the residual r of G at x does not bound the
	   error of G(x) when G contracts slowly (error ~ r / (1 - lambda)),
	   and the residual history of ANDERSON and NEWTON_KRYLOV does not
	   give lambda. Once the residual is small enough, ERROR_STEPS plain
	   G steps from G(x) are taken, and minimal polynomial extrapolation
	   (mpe.c) of them, in the logarithms of the populations, rho, ne
	   and J, estimates the fixed point s; the estimated error is
	   max |s / G(x) - 1| over those entries. The run stops when it is
	   below ERROR_LIMIT (and the residual of the populations below
	   ITER_LIMIT). Otherwise the residual target of the next estimate
	   is lowered in proportion, 0.8 r ERROR_LIMIT / e, at most 0.8 r.

	   Measured against converged runs (docs/changes/N-3.md): the
	   estimate is 0.6 to 1.4 times the actual error with 6 steps, in
	   1-D, 2-D and 3-D, also where plain G is not a contraction at the
	   solution (1d_shock_h6crd_ne). G is a function of x alone (S-2),
	   so the steps change nothing the solver carries on with; before
	   a stop G(x) is evaluated again, so that what G sets besides the
	   state (Gamma, the radiative rates, the intensities) is that of
	   the output --                                -------------- */

#define ERROR_STEPS 6

typedef struct {
	double limit;    /* ERROR_LIMIT, 0: off */
	double target;   /* residual at which the next estimate is made */
	double residual; /* r and the estimate of the last one, -1: none */
	double estimate;
	int    Nestimate;
	long   Geval;
} ErrorTest;

static bool_t errorBlock(enum BlockKind kind)
{
	return (kind != BLOCK_J20);
}

static long errorSize(void)
{
	register int n;

	long N = 0;

	for (n = 0; n < Nblocks; n++)
		if (errorBlock(blocks[n].kind)) N += blocks[n].N;
	return N;
}

static void errorTestInit(ErrorTest *t, const char *routineName)
{
	t->limit     = input.errorLimit;
	t->target    = input.errorLimit;
	t->residual  = -1.0;
	t->estimate  = -1.0;
	t->Nestimate = 0;
	t->Geval     = 0;
	if (t->limit <= 0.0) return;

	sprintf(messageStr,
			" -- ERROR_LIMIT = %.2E: the run stops when the estimated error"
			" of the\n    populations, rho, ne and J is below it%s; each"
			" estimate takes %d plain G\n    steps (memory %.1f MB)\n",
			t->limit,
			(input.iterLimit_set) ? " and the residual of the populations below"
									" ITER_LIMIT"
								  : "",
			ERROR_STEPS,
			((double) ERROR_STEPS * errorSize() + 2.0 * Nstate) * sizeof(double) /
				1.0E6);
	Error(MESSAGE, routineName, messageStr);
}

/* --- The residual of rho, ne and J (relChange) -- ----------------- */

static double errorResidual(const double *x, const double *gx)
{
	register int n;

	double r = 0.0;

	for (n = 0; n < Nblocks; n++)
		if (blocks[n].kind == BLOCK_RHO || blocks[n].kind == BLOCK_NE ||
			blocks[n].kind == BLOCK_J)
			r = MAX(r, relChange(x + blockStart[n], gx + blockStart[n], blocks[n].N));
	return r;
}

/* --- log y - log y0 for the entries in the estimate, 0 where either
	   is not positive. Both logarithms are of positive numbers in
	   every lane: a vectorizing compiler may evaluate both sides of
	   the select (FPE traps, as in N-37) --         -------------- */

static void errorLogDiff(const double *y, const double *y0, double *out)
{
	register int  n;
	register long i;

	for (n = 0; n < Nblocks; n++) {
		if (!errorBlock(blocks[n].kind)) continue;
		for (i = blockStart[n]; i < blockStart[n + 1]; i++)
			*out++ = (y[i] > 0.0 && y0[i] > 0.0)
						 ? log(fmax(y[i], DBL_MIN)) - log(fmax(y0[i], DBL_MIN))
						 : 0.0;
	}
}

static double errorEstimate(ErrorTest *t, const double *gx)
{
	register int  j;
	register long i;

	bool_t  Gok = TRUE;
	long    N   = errorSize();
	double *y, *gy, *D[ERROR_STEPS], e = -1.0;

	/* --- ERROR_STEPS plain steps from G(x) as trial points (a step
		   that fails, N-40, N-8, N-4, gives no estimate instead of
		   stopping the run). The differences of the logarithms, in
		   terms of G(x): D_j = log(y_(j+1) / G(x)) - log(y_j / G(x)).
		   Returns the estimated error of G(x), or -1 -- ------------- */

	y  = (double *) malloc(Nstate * sizeof(double));
	gy = (double *) malloc(Nstate * sizeof(double));
	for (j = 0; j < ERROR_STEPS; j++) D[j] = (double *) malloc(N * sizeof(double));

	memcpy(y, gx, Nstate * sizeof(double));
	for (j = 0; j < ERROR_STEPS && Gok; j++) {
		negTauReset();
		negTauTrial(TRUE);
		Gok = G_eval(y, gy);
		negTauTrial(FALSE);
		t->Geval++;
		if (!Gok) break;
		errorLogDiff(gy, gx, D[j]);
		memcpy(y, gy, Nstate * sizeof(double));
	}
	negTauReset();
	if (Gok) {
		/* --- D[j] holds log(y_(j+1) / G(x)): differences in place -- */

		for (j = ERROR_STEPS - 1; j > 0; j--)
			for (i = 0; i < N; i++) D[j][i] -= D[j - 1][i];
		if (mpeExtrapolate(ERROR_STEPS, N, D)) {
			e = 0.0;
			for (i = 0; i < N; i++) e = MAX(e, fabs(expm1(D[ERROR_STEPS - 1][i])));
			if (!isfinite(e)) e = -1.0;
		}
	}
	for (j = 0; j < ERROR_STEPS; j++) free(D[j]);
	free(y);
	free(gy);
	t->Nestimate++;

	return e;
}

/* --- Whether the run stops at x, with dpopsmax the residual of the
	   populations: below ITER_LIMIT, and with ERROR_LIMIT an estimated
	   error below it -- ------------------------------------------ */

static bool_t errorTestStop(ErrorTest *t, const char *routineName, int niter,
							const double *x, const double *gx, double dpopsmax,
							double iterLimit)
{
	double r, e, *gy;

	if (!(dpopsmax < iterLimit)) return FALSE;
	if (t->limit <= 0.0) return TRUE;
	r = MAX(dpopsmax, errorResidual(x, gx));
	if (!(r < t->target)) return FALSE;

	e           = errorEstimate(t, gx);
	t->residual = r;
	t->estimate = e;
	if (e >= 0.0 && e < t->limit) {
		sprintf(messageStr,
				" -- Iteration %d: residual %.2E, estimated error %.2E <"
				" ERROR_LIMIT\n",
				niter, r, e);
		Error(MESSAGE, routineName, messageStr);

		/* --- What G sets besides the state, that of G(x) -- -------- */

		gy = (double *) malloc(Nstate * sizeof(double));
		G_eval(x, gy);
		t->Geval++;
		free(gy);
		return TRUE;
	}
	t->target = (e > 0.0) ? MIN(0.8 * r * t->limit / e, 0.8 * r) : 0.5 * r;
	if (e >= 0.0)
		sprintf(messageStr,
				" -- Iteration %d: residual %.2E, estimated error %.2E >="
				" ERROR_LIMIT;\n    next estimate at residual %.2E\n",
				niter, r, e, t->target);
	else
		sprintf(messageStr,
				" -- Iteration %d: residual %.2E, no error estimate (a plain"
				" step failed,\n    or no finite limit); next estimate at"
				" residual %.2E\n",
				niter, r, t->target);
	Error(MESSAGE, routineName, messageStr);
	return FALSE;
}

/* --- The estimated error of the output gx of a run that did not
	   converge (not after a stop on request) -- ------------------- */

static void errorTestFinal(ErrorTest *t, const double *gx)
{
	if (t->limit <= 0.0) return;
	t->estimate = errorEstimate(t, gx);
}

/* --- The limits part of the convergence messages -- -------------- */

static void errorTestLimits(const ErrorTest *t, char *s, size_t n, double dpopsmax,
							double iterLimit, bool_t converged)
{
	char est[32];

	if (t->limit <= 0.0) {
		snprintf(s, n, "max change %.4E %s ITER_LIMIT = %.4E", dpopsmax,
				 (converged) ? "<" : ">=", iterLimit);
		return;
	}
	if (t->estimate >= 0.0)
		snprintf(est, sizeof(est), "%.2E %s", t->estimate,
				 (t->estimate < t->limit) ? "<" : ">=");
	else
		snprintf(est, sizeof(est), "unknown,");
	if (input.iterLimit_set)
		snprintf(s, n,
				 "max change %.4E %s ITER_LIMIT = %.4E,\n    estimated error %s"
				 " ERROR_LIMIT = %.2E",
				 dpopsmax, (dpopsmax < iterLimit) ? "<" : ">=", iterLimit, est,
				 t->limit);
	else
		snprintf(s, n, "max change %.4E, estimated error %s ERROR_LIMIT = %.2E",
				 dpopsmax, est, t->limit);
}

static void errorTestSummary(const ErrorTest *t, const char *routineName)
{
	if (t->limit <= 0.0) return;
	sprintf(messageStr, " -- ERROR_LIMIT: %d error estimate%s, %ld G evaluations\n",
			t->Nestimate, (t->Nestimate == 1) ? "" : "s", t->Geval);
	Error(MESSAGE, routineName, messageStr);
}
/* ------- end ---------------------------- errorTest.c ------------- */

/* ------- begin -------------------------- iterateAnderson.c ------- */

static const char *extrapolateName[] = {"ALL", "POPULATIONS", "POPULATIONS_RHO"};

static bool_t carriedEqual(const double *x, const double *y)
{
	register int n, p;

	/* --- The carried blocks agree bit for bit -- --------------- */

	for (p = 0; p < Ncarry; p++) {
		n = carryBlock[p];
		if (memcmp(x + blockStart[n], y + blockStart[n], blocks[n].N * sizeof(double)))
			return FALSE;
	}
	return TRUE;
}

void iterateAnderson(int NmaxIter, double iterLimit)
{
	const char   routineName[] = "iterateAnderson";
	register int niter, n;

	bool_t    converged = FALSE, stopped = FALSE, Gok;
	char      what[256], limits[256];
	int       m = input.Anderson_depth, Nsolve, delay, t;
	ErrorTest et;
	long      i, Ndamped = 0, Nguardplain = 0, Ncarried = 0, Ncarried_ok = 0,
			Nbadtrial = 0;
	double             *x, *gx, *gacc, *u, *gu, *un, *w, *gamma, dpopsmax = 0.0, theta;
	unsigned char      *valid, *validg;
	enum AndersonAction action;
	enum GuardResult    guard;
	AndersonOptions     opt;
	AndersonSolve      *info;
	Anderson           *a;

	stateInit();
	/* --- x, G(x) and G of the last accepted point are full states;
		   the variables and the history cover the extrapolated blocks
		   only (ANDERSON_EXTRAPOLATE, S-11) --         -------------- */

	x      = (double *) malloc(Nstate * sizeof(double));
	gx     = (double *) malloc(Nstate * sizeof(double));
	gacc   = (double *) malloc(Nstate * sizeof(double));
	w      = setupVariables(input.Anderson_variables, input.Anderson_extrapolate, x);
	u      = (double *) malloc(Npacked * sizeof(double));
	gu     = (double *) malloc(Npacked * sizeof(double));
	un     = (double *) malloc(Npacked * sizeof(double));
	valid  = (unsigned char *) malloc(Npacked);
	validg = (unsigned char *) malloc(Npacked);
	Nsolve = (input.Anderson_coefficients == ANDERSON_PER_BLOCK) ? Npack : 1;
	gamma  = (double *) calloc((size_t) Nsolve * (m > 0 ? m : 1), sizeof(double));
	info   = (AndersonSolve *) calloc(Nsolve, sizeof(AndersonSolve));
	a      = andersonInit(Npacked, m);

	opt.Nblock    = Npack;
	opt.start     = packStart;
	opt.per_block = (input.Anderson_coefficients == ANDERSON_PER_BLOCK);
	opt.cond_max  = input.Anderson_cond_max;
	opt.safeguard = input.Anderson_safeguard;
	delay         = MAX(input.Anderson_delay, 2);

	sprintf(messageStr,
			"\n -- SOLVER = ANDERSON: %ld unknowns in %d blocks; PRD"
			" sub-iterations per G: %d;\n    ANDERSON_DEPTH = %d,"
			" ANDERSON_DELAY = %d, ANDERSON_COND_MAX = %.1E,\n"
			"    ANDERSON_SAFEGUARD = %g, ANDERSON_VARIABLES = %s,"
			" ANDERSON_COEFFICIENTS = %s,\n    ANDERSON_EXTRAPOLATE = %s%s:"
			" %ld unknowns in %d blocks extrapolated;\n    memory: history"
			" and variables %.1f MB, x, G(x) and G at the last accepted"
			" point %.1f MB\n",
			Nstate, Nblocks, (atmos.NPRDactive > 0) ? input.PRD_NmaxIter : 0, m,
			input.Anderson_delay, input.Anderson_cond_max,
			input.Anderson_safeguard,
			(varKind == VARIABLES_LOG) ? "LOG" : "RELATIVE",
			(opt.per_block) ? "PER_BLOCK" : "GLOBAL",
			extrapolateName[input.Anderson_extrapolate],
			input.Anderson_extrapolate_default ? " (the default)" : "", Npacked, Npack,
			((3.0 * m + 7 + (varKind == VARIABLES_RELATIVE)) * sizeof(double) + m + 4) *
				Npacked / 1.0E6,
			3.0 * sizeof(double) * Nstate / 1.0E6);
	Error(MESSAGE, routineName, messageStr);
	/* --- Without PRD lines there is no rho to extrapolate: say so
		   when the keyword file asked for POPULATIONS_RHO, the run goes
		   on as POPULATIONS (S-11). The 2-D and 3-D default says
		   nothing: the banner lists the extrapolated blocks -- ----- */

	if (input.Anderson_extrapolate == EXTRAPOLATE_POPULATIONS_RHO &&
		!input.Anderson_extrapolate_default) {
		for (n = 0; n < Nblocks && blocks[n].kind != BLOCK_RHO; n++);
		if (n == Nblocks) {
			sprintf(messageStr,
					"ANDERSON_EXTRAPOLATE = POPULATIONS_RHO, but there are no PRD"
					" lines:\n the same as POPULATIONS");
			Error(WARNING, routineName, messageStr);
		}
	}
	/* --- The carried blocks, several lines if needed -- -------- */

	for (t = 0; t < Ncarry;) {
		sprintf(messageStr, "    %s", (t == 0) ? "carried from G:" : "               ");
		do {
			sprintf(messageStr + strlen(messageStr), " %s%s", blocks[carryBlock[t]].name,
					(t < Ncarry - 1) ? ";" : "");
			t++;
		} while (t < Ncarry &&
				 strlen(messageStr) + strlen(blocks[carryBlock[t]].name) + 4 < 80);
		strcat(messageStr, "\n");
		Error(MESSAGE, NULL, messageStr);
	}
	errorTestInit(&et, routineName);

	niter = 1;
	while (niter <= NmaxIter) {
		if (StopRequested()) {
			stopped = TRUE;
			break;
		}
		getCPU(2, TIME_START, NULL);

		/* --- The carried blocks of every point are those of G at the
			   last accepted point (S-11), checked bit for bit -- ---- */

		if (Ncarry > 0 && niter > 1) {
			Ncarried++;
			if (carriedEqual(x, gacc)) Ncarried_ok++;
		}

		/* --- An extrapolated point whose total opacity is negative
			   enough to overflow the formal solution is rejected like
			   one that fails the safeguard (N-40) --      -------------- */

		negTauReset();
		negTauTrial(a->trial);
		Gok = G_eval(x, gx);
		negTauTrial(FALSE);
		if (!Gok) {
			negTauMessage(what, sizeof(what));
			sprintf(messageStr,
					"\n -- Iteration %3d\n  Anderson: %s at the extrapolated point:"
					" point rejected; plain step from the last accepted point\n",
					niter, what);
			Error(MESSAGE, NULL, messageStr);
			Nbadtrial++;
			negTauReset();
			andersonRejectTrial(a);
			memcpy(x, gacc, Nstate * sizeof(double));
			sprintf(messageStr, "Total Iteration %3d", niter);
			getCPU(2, TIME_POLL, messageStr);
			niter++;
			continue;
		}

		dpopsmax = printChange(routineName, niter, x, gx);
		if (errorTestStop(&et, routineName, niter, x, gx, dpopsmax, iterLimit))
			converged = TRUE;
		sprintf(messageStr, "Total Iteration %3d", niter);
		getCPU(2, TIME_POLL, messageStr);
		if (converged) break;
		/* --- The next point --                          -------------- */

		toVars(x, u, valid);
		toVars(gx, gu, validg);
		for (i = 0; i < Npacked; i++) valid[i] &= validg[i];

		action = andersonStep(a, &opt, u, gu, valid, w, (m > 0 && niter >= delay),
							  un, gamma, info);
		switch (action) {
		case ANDERSON_REJECTED:
			sprintf(messageStr,
					"  Anderson: point rejected, |f|_W %.3E > %g x %.3E;"
					" plain step from the last accepted point\n",
					andersonNorm(Npacked, u, gu, valid, w), opt.safeguard,
					a->fnorm_acc);
			Error(MESSAGE, NULL, messageStr);
			memcpy(x, gacc, Nstate * sizeof(double));
			break;

		case ANDERSON_EXTRAPOLATED:
			logSolve(Nsolve, opt.per_block, info, a->fnorm_acc);
			memcpy(gacc, gx, Nstate * sizeof(double));
			fromVars(un, gx, andersonMask(a), x);
			guard = guardTrial(x, gx, NULL);
			for (t = 1, theta = 0.5; guard != GUARD_OK && t <= 2; t++, theta *= 0.5) {
				sprintf(messageStr, "  Anderson: guard: %s; step damped to %g\n",
						guardName[guard], theta);
				Error(MESSAGE, NULL, messageStr);
				Ndamped++;
				andersonDamp(a, &opt, gamma, info, theta, un);
				fromVars(un, gx, andersonMask(a), x);
				guard = guardTrial(x, gx, NULL);
			}
			if (guard != GUARD_OK) {
				sprintf(messageStr, "  Anderson: guard: %s; plain step,"
									" history restarted\n",
						guardName[guard]);
				Error(MESSAGE, NULL, messageStr);
				Nguardplain++;
				andersonAbandon(a, un);
				memcpy(x, gx, Nstate * sizeof(double));
				break;
			}
			break;

		default: /* plain step, copied --                 -------------- */
			if (m > 0 && niter >= delay && andersonNcols(a) > 0)
				logSolve(Nsolve, opt.per_block, info, a->fnorm_acc);
			memcpy(gacc, gx, Nstate * sizeof(double));
			memcpy(x, gx, Nstate * sizeof(double));
			break;
		}
		niter++;
	}

	if (converged) {
		errorTestLimits(&et, limits, sizeof(limits), dpopsmax, iterLimit, TRUE);
		sprintf(messageStr, "\n -- Converged in %d iterations: %s\n", niter, limits);
		Error(MESSAGE, routineName, messageStr);
	} else {
		/* --- The output is G of the last accepted point: the state
			   arrays may hold the image of a rejected one -- ------- */

		if (!stopped && niter > 1) errorTestFinal(&et, gacc);
		if (niter > 1) stateScatter(gacc);
		if (stopped) {
			sprintf(messageStr, "Stopped on request (STOP_RH) after %d iterations:"
								" NOT CONVERGED.\n The output is not a converged"
								" solution",
					niter - 1);
		} else {
			errorTestLimits(&et, limits, sizeof(limits), dpopsmax, iterLimit, FALSE);
			sprintf(messageStr, "NOT CONVERGED after %d iterations: %s.\n The output"
								" is not a converged solution; increase N_MAX_ITER or"
								" see what limits convergence",
					niter - 1, limits);
		}
		Error(WARNING, routineName, messageStr);
		SetExitFlag(RH_EXIT_NOT_CONVERGED);
	}
	sprintf(messageStr,
			" -- G evaluations: %ld\n -- Anderson: extrapolated %ld, plain %ld,"
			" rejected %ld; guard: damped %ld, plain %ld; columns dropped %ld\n",
			Gcount, a->Nextrapolated, a->Nplain, a->Nrejected, Ndamped, Nguardplain,
			a->Ndropped);
	Error(MESSAGE, routineName, messageStr);
	errorTestSummary(&et, routineName);
	if (Nbadtrial > 0) {
		sprintf(messageStr,
				" -- Anderson: %ld of the rejected points were invalid trial points"
				" (negative total opacity, negative populations or non-positive PRD rho;"
				" N-40, N-8, N-4)\n",
				Nbadtrial);
		Error(MESSAGE, routineName, messageStr);
	}
	if (Ncarry > 0) {
		sprintf(messageStr,
				" -- Anderson: carried blocks equal to G at the last accepted"
				" point in %ld of %ld evaluations\n",
				Ncarried_ok, Ncarried);
		Error(MESSAGE, routineName, messageStr);
	}

	if (input.solve_ne == ITERATION) BackgroundOfState(TRUE);

	for (n = 0; n < atmos.Nactiveatom; n++)
		freeMatrix((void **) atmos.activeatoms[n]->Gamma);
	for (n = 0; n < atmos.Nactivemol; n++)
		freeMatrix((void **) atmos.activemols[n]->Gamma);

	andersonFree(a);
	free(x);
	free(gx);
	free(gacc);
	free(u);
	free(gu);
	free(un);
	free(w);
	free(valid);
	free(validg);
	free(gamma);
	free(info);
}
/* ------- end ---------------------------- iterateAnderson.c ------- */

/* ------- begin -------------------------- fdProduct.c ------------- */

/* --- SOLVER = NEWTON_KRYLOV (S-4 in docs/SOLVER.md): inexact Newton
	   on R(u) = s(G(s^-1(u))) - u, with GMRES (newton.c) on
	   finite-difference products of the Jacobian, Eisenstat-Walker
	   forcing and a backtracking line search guarded as Anderson's
	   steps. --                                       -------------- */

typedef struct {
	const double        *u, *r, *x, *w; /* base point: u, R(u), x */
	const unsigned char *mask;          /* entries in the step */
	double              *up, *xp, *gp, *gup, h;
	unsigned char       *validp;
	long                 Nprobe, Ninvalid, Nbadtrial;
} FDContext;

static void fdProduct(void *vctx, const double *v, double *y)
{
	register long i;

	FDContext *c = (FDContext *) vctx;
	double     vnorm, eps;

	/* --- y = (R(u + eps v) - R(u)) / eps, eps = h / ||v||_W: one G
		   evaluation. Outside the mask the probe keeps the base x: with
		   G(x) there, as a trial has, the difference would carry an
		   O(1) / eps jump. So J does not see how those entries move in a
		   trial; the line search does. A component whose probe image
		   has no value in the variables (G <= 0 in LOG) is set to 0 and
		   counted -- --                                -------------- */

	vnorm = weightedNorm(Nstate, v, c->mask, c->w);
	eps   = c->h / ((vnorm > 0.0) ? vnorm : 1.0);
	for (i = 0; i < Nstate; i++) c->up[i] = (c->mask[i]) ? c->u[i] + eps * v[i] : 0.0;
	fromVars(c->up, c->x, c->mask, c->xp);

	/* --- A probe whose total opacity is negative enough to overflow
		   (N-40) has no image: its product is zero, as for invalid
		   entries, and GMRES works with the directions it has -- --- */

	negTauReset();
	negTauTrial(TRUE);
	if (!G_eval(c->xp, c->gp)) {
		negTauTrial(FALSE);
		negTauReset();
		c->Nprobe++;
		c->Nbadtrial++;
		for (i = 0; i < Nstate; i++) y[i] = 0.0;
		return;
	}
	negTauTrial(FALSE);
	c->Nprobe++;
	toVars(c->gp, c->gup, c->validp);
	for (i = 0; i < Nstate; i++) {
		if (!c->mask[i]) {
			y[i] = 0.0;
		} else if (!c->validp[i]) {
			y[i] = 0.0;
			c->Ninvalid++;
		} else {
			y[i] = ((c->gup[i] - c->up[i]) - c->r[i]) / eps;
		}
	}
}
/* ------- end ---------------------------- fdProduct.c ------------- */

/* ------- begin -------------------------- residualOf.c ------------ */

static double residualOf(const double *x, const double *gx, double *u,
						 double *gu, unsigned char *valid,
						 unsigned char *validg, double *r, const double *w)
{
	register long i;

	/* --- u = s(x), r = s(G(x)) - u on the entries valid in both
		   (valid on return), 0 elsewhere; returns ||r||_W -- ------ */

	toVars(x, u, valid);
	toVars(gx, gu, validg);
	for (i = 0; i < Nstate; i++) {
		valid[i] &= validg[i];
		r[i] = (valid[i]) ? gu[i] - u[i] : 0.0;
	}
	return weightedNorm(Nstate, r, valid, w);
}
/* ------- end ---------------------------- residualOf.c ------------ */

/* ------- begin -------------------------- iterateNewton.c --------- */

void iterateNewton(int NmaxIter, double iterLimit)
{
	const char   routineName[] = "iterateNewton";
	register int niter;

	bool_t converged = FALSE, stopped = FALSE, accepted, Gok;
	char   limits[256];
	int    k = input.Newton_krylov_max, kmax, nback, Nhalve, n, ninf,
		Nreason[GUARD_NRESULT], Nbadtrial = 0;
	ErrorTest et;
	long      i, Neval = 0, Nbase = 0, Nkrylov = 0, Nline = 0, Nplain = 0,
			Nguard = 0, Nsteps = 0, Nlost = 0, Nrestore = 0;
	double *x, *gx, *u, *gu, *r, *d, *xt, *gt, *ut, *gut, *rt, *w;
	double  dpopsmax = 0.0, rnorm, rnorm_prev = 0.0, eta = 0.0, rtnorm, dinf,
		   lambda, h = sqrt(input.Newton_noise);
	unsigned char   *valid, *validg, *validt, *validgt;
	enum GuardResult guard;
	GmresResult      gres;
	FDContext        ctx;
	Gmres           *gm;

	stateInit();
	x       = (double *) malloc(Nstate * sizeof(double));
	gx      = (double *) malloc(Nstate * sizeof(double));
	u       = (double *) malloc(Nstate * sizeof(double));
	gu      = (double *) malloc(Nstate * sizeof(double));
	r       = (double *) malloc(Nstate * sizeof(double));
	d       = (double *) malloc(Nstate * sizeof(double));
	xt      = (double *) malloc(Nstate * sizeof(double));
	gt      = (double *) malloc(Nstate * sizeof(double));
	ut      = (double *) malloc(Nstate * sizeof(double));
	gut     = (double *) malloc(Nstate * sizeof(double));
	rt      = (double *) malloc(Nstate * sizeof(double));
	w       = setupVariables(input.Newton_variables, EXTRAPOLATE_ALL, x);
	valid   = (unsigned char *) malloc(Nstate);
	validg  = (unsigned char *) malloc(Nstate);
	validt  = (unsigned char *) malloc(Nstate);
	validgt = (unsigned char *) malloc(Nstate);
	gm      = gmresInit(Nstate, k);

	/* --- The probes of fdProduct use the trial buffers: they are
		   never alive at the same time --              -------------- */

	ctx.u         = u;
	ctx.r         = r;
	ctx.x         = x;
	ctx.w         = w;
	ctx.mask      = valid;
	ctx.up        = ut;
	ctx.xp        = xt;
	ctx.gp        = gt;
	ctx.gup       = gut;
	ctx.validp    = validgt;
	ctx.h         = h;
	ctx.Nprobe    = 0;
	ctx.Ninvalid  = 0;
	ctx.Nbadtrial = 0;

	sprintf(messageStr,
			"\n -- SOLVER = NEWTON_KRYLOV: %ld unknowns in %d blocks; PRD"
			" sub-iterations per G: %d;\n    NEWTON_VARIABLES = %s,"
			" NEWTON_NOISE = %.2E (finite-difference step h = %.2E),\n"
			"    NEWTON_KRYLOV_MAX = %d (no restart), NEWTON_ETA_MAX = %g;\n"
			"    forcing: Eisenstat-Walker choice 2, gamma = %g, alpha = %g,"
			" safeguard above %g,\n    floor 0.5 x ITER_LIMIT / |R|_W;"
			" line search: lambda = 1 ... 1/%d, Armijo %.0E;\n"
			"    N_MAX_ITER caps G evaluations; memory %.1f MB\n",
			Nstate, Nblocks, (atmos.NPRDactive > 0) ? input.PRD_NmaxIter : 0,
			(varKind == VARIABLES_LOG) ? "LOG" : "RELATIVE", input.Newton_noise,
			h, k, input.Newton_eta_max, NEWTON_EW_GAMMA, NEWTON_EW_ALPHA,
			NEWTON_EW_SAFE, 1 << NEWTON_N_BACKTRACK, NEWTON_ARMIJO,
			((k + 1.0 + 12 + (varKind == VARIABLES_RELATIVE)) * sizeof(double) + 4) *
				Nstate / 1.0E6);
	Error(MESSAGE, routineName, messageStr);

	errorTestInit(&et, routineName);

	/* --- The first base point --                     -------------- */

	G_eval(x, gx);
	Neval++;
	Nbase++;

	niter = 1;
	for (;;) {
		dpopsmax = printChange(routineName, niter, x, gx);
		if (errorTestStop(&et, routineName, niter, x, gx, dpopsmax, iterLimit)) {
			converged = TRUE;
			break;
		}
		if (StopRequested()) {
			stopped = TRUE;
			break;
		}
		if (Neval >= NmaxIter) break;

		/* --- Newton step: J d = -R by GMRES, at most k products and
			   one evaluation left for the line search -- ---------- */

		getCPU(2, TIME_START, NULL);
		rnorm = residualOf(x, gx, u, gu, valid, validg, r, w);
		eta   = newtonForcing(Nsteps == 0, eta, rnorm, rnorm_prev,
							  input.Newton_eta_max,
                            (et.limit > 0.0) ? MIN(iterLimit, et.target) : iterLimit);
		kmax  = (int) MIN((long) k, NmaxIter - Neval - 1);
		for (i = 0; i < Nstate; i++) rt[i] = -r[i];
		memset(&gres, 0, sizeof(gres));
		gres.rnorm0 = gres.rnorm = rnorm;
		memset(d, 0, Nstate * sizeof(double));
		if (kmax >= 1) {
			ctx.Nprobe = 0;
			gres       = gmresSolve(gm, fdProduct, &ctx, rt, valid, w, kmax, eta * rnorm, d);
			Neval += ctx.Nprobe;
			Nkrylov += ctx.Nprobe;
		}
		Nsteps++;

		/* --- Backtracking on ||R||_W, the guard before every trial,
			   the plain step if no lambda is accepted -- ---------- */

		accepted  = FALSE;
		Nhalve    = 0;
		Nbadtrial = 0;
		memset(Nreason, 0, sizeof(Nreason));
		dinf = 0.0;
		ninf = 0;
		for (n = 0; n < Nblocks; n++)
			for (i = blockStart[n]; i < blockStart[n + 1]; i++)
				if (valid[i] && fabs(d[i]) > dinf) {
					dinf = fabs(d[i]);
					ninf = n;
				}
		lambda = 1.0;
		rtnorm = 0.0;
		if (kmax >= 1 && gres.Niter > 0) {
			for (nback = 0; nback <= NEWTON_N_BACKTRACK && Neval < NmaxIter;
				 nback++, lambda *= 0.5) {
				guard = GUARD_OK;
				for (n = 0; n < Nblocks; n++)
					for (i = blockStart[n]; i < blockStart[n + 1]; i++) {
						ut[i] = (valid[i]) ? u[i] + lambda * d[i] : 0.0;
						if (varKind == VARIABLES_LOG && blocks[n].kind != BLOCK_J20 &&
							ut[i] > LOG_DBL_MAX)
							guard = GUARD_NONFINITE; /* exp would overflow */
					}
				if (guard == GUARD_OK) {
					fromVars(ut, gx, valid, xt);
					guard = guardTrial(xt, gx, x);
				}
				if (guard != GUARD_OK) {
					Nhalve++;
					Nguard++;
					Nreason[guard]++;
					continue;
				}
				negTauReset();
				negTauTrial(TRUE);
				Gok = G_eval(xt, gt);
				negTauTrial(FALSE);
				Neval++;
				if (!Gok) { /* negative opacity (N-40): halve */
					Nhalve++;
					Nguard++;
					Nbadtrial++;
					negTauReset();
					continue;
				}
				Nline++;
				rtnorm = residualOf(xt, gt, ut, gut, validt, validgt, rt, w);

				/* --- Both norms over the same entries: a trial that
					   loses one of the base's (G(trial) <= 0 in LOG)
					   is not accepted, or dropping it could pass --  */

				for (i = 0; i < Nstate; i++)
					if (valid[i] && !validt[i]) break;
				if (i < Nstate) {
					Nlost++;
					continue;
				}
				if (rtnorm <= (1.0 - NEWTON_ARMIJO * lambda) * rnorm) {
					accepted = TRUE;
					break;
				}
			}
		}
		if (accepted) {
			sprintf(messageStr,
					"  Newton step %ld: |R|_W %.3E, eta %.2E; GMRES %d products,"
					" |R + J d|_W %.3E%s;\n  |d|_max %.2E (%s); lambda %g"
					" (guard halvings %d: %d not finite, %d negative, %d invalid trial),"
					" |R(trial)|_W %.3E; G evaluations %ld\n",
					Nsteps, rnorm, eta, gres.Niter, gres.rnorm,
					(gres.breakdown) ? " (invariant space)" : "", dinf,
					blocks[ninf].name, lambda, Nhalve, Nreason[GUARD_NONFINITE],
					Nreason[GUARD_NEGATIVE], Nbadtrial, rtnorm, Neval);
			Error(MESSAGE, NULL, messageStr);
			memcpy(x, xt, Nstate * sizeof(double));
			memcpy(gx, gt, Nstate * sizeof(double));
		} else {
			sprintf(messageStr,
					"  Newton step %ld: |R|_W %.3E, eta %.2E; GMRES %d products,"
					" |R + J d|_W %.3E%s;\n  |d|_max %.2E (%s); no lambda accepted"
					" (guard halvings %d: %d not finite, %d negative, %d invalid trial):"
					" plain step; G evaluations %ld\n",
					Nsteps, rnorm, eta, gres.Niter, gres.rnorm,
					(gres.breakdown) ? " (invariant space)" : "", dinf,
					blocks[ninf].name, Nhalve, Nreason[GUARD_NONFINITE],
					Nreason[GUARD_NEGATIVE], Nbadtrial, Neval);
			Error(MESSAGE, NULL, messageStr);
			Nplain++;
			if (Neval >= NmaxIter) {

				/* --- The last evaluation was a probe or a rejected
					   trial: one more at the base point (G is
					   deterministic, so gx does not change) makes what
					   G sets besides the state agree with the output,
					   one evaluation past the cap -- ------------- */

				G_eval(x, gx);
				Neval++;
				Nrestore++;
				getCPU(2, TIME_POLL, "Newton step");
				break;
			}
			memcpy(x, gx, Nstate * sizeof(double));
			G_eval(x, gx);
			Neval++;
			Nbase++;
		}
		rnorm_prev = rnorm;
		getCPU(2, TIME_POLL, "Newton step");
		niter++;
	}

	/* --- The output is G of the last base point. After the error
		   estimate of a run that did not converge, G there again, for
		   what G sets besides the state -- -------------------------- */

	if (!converged && !stopped && et.limit > 0.0) {
		errorTestFinal(&et, gx);
		G_eval(x, gx);
		et.Geval++;
	}
	stateScatter(gx);
	if (converged) {
		errorTestLimits(&et, limits, sizeof(limits), dpopsmax, iterLimit, TRUE);
		sprintf(messageStr, "\n -- Converged at base point %d: %s\n", niter, limits);
		Error(MESSAGE, routineName, messageStr);
	} else {
		if (stopped) {
			sprintf(messageStr, "Stopped on request (STOP_RH) after %ld G"
								" evaluations: NOT CONVERGED.\n The output is not a"
								" converged solution",
					Neval);
		} else {
			errorTestLimits(&et, limits, sizeof(limits), dpopsmax, iterLimit, FALSE);
			sprintf(messageStr, "NOT CONVERGED after %ld G evaluations"
								" (N_MAX_ITER = %d): %s.\n The output is not a converged"
								" solution; increase N_MAX_ITER or see what limits"
								" convergence",
					Neval, NmaxIter, limits);
		}
		Error(WARNING, routineName, messageStr);
		SetExitFlag(RH_EXIT_NOT_CONVERGED);
	}
	sprintf(messageStr,
			" -- G evaluations: %ld\n -- Newton: steps %ld, G evaluations base %ld,"
			" Krylov %ld, line search %ld (lost an entry %ld); plain steps %ld;"
			" guard halvings %ld; invalid probe entries %ld; restoring %ld\n",
			Neval, Nsteps, Nbase, Nkrylov, Nline, Nlost, Nplain, Nguard,
			ctx.Ninvalid, Nrestore);
	Error(MESSAGE, routineName, messageStr);
	errorTestSummary(&et, routineName);
	if (ctx.Nbadtrial > 0) {
		sprintf(messageStr,
				" -- Newton: %ld Jacobian probes were invalid trial points (negative"
				" total opacity, populations or PRD rho, N-40, N-8, N-4) and gave a zero product\n",
				ctx.Nbadtrial);
		Error(MESSAGE, routineName, messageStr);
	}

	if (input.solve_ne == ITERATION) BackgroundOfState(TRUE);

	for (i = 0; i < atmos.Nactiveatom; i++)
		freeMatrix((void **) atmos.activeatoms[i]->Gamma);
	for (i = 0; i < atmos.Nactivemol; i++)
		freeMatrix((void **) atmos.activemols[i]->Gamma);

	gmresFree(gm);
	free(x);
	free(gx);
	free(u);
	free(gu);
	free(r);
	free(d);
	free(xt);
	free(gt);
	free(ut);
	free(gut);
	free(rt);
	free(w);
	free(valid);
	free(validg);
	free(validt);
	free(validgt);
}
/* ------- end ---------------------------- iterateNewton.c --------- */
