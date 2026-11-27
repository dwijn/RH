/* ------- file: -------------------------- solver.h ----------------

	   The map G of docs/SOLVER.md (S-2) and the solvers built on it.

	   --------------------------                      ----------RH-- */

#ifndef __SOLVER_H__
#define __SOLVER_H__

#include "inputs.h" /* enum solver_variables; needs rh.h first */

/* --- State vector x (populations, molecular populations, PRD rho, J)
	   and the map G, see solver.c --                 -------------- */

long        stateInit(RHState *rh);
void        stateGather(RHState *rh, double *x);
void        stateScatter(RHState *rh, const double *x);
int         stateNblocks(RHState *rh);
const char *stateBlockName(RHState *rh, int n);
long        stateBlockSize(RHState *rh, int n);
bool_t      G_eval(RHState *rh, const double *x_in, double *x_out);
long        G_count(RHState *rh);

/* --- The state of the solver between calls, its part of RHState
	   (rh->solver, Q-46) --                          -------------- */

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

#define RES_HIST_MAX 11 /* residuals kept for the contraction estimate (N-3) */

typedef struct {
	StateBlock *blocks; /* the blocks of x, rebuilt by stateInit */
	int         Nblocks;
	long        Nstate;     /* length of x */
	long       *blockStart; /* Nblocks + 1 offsets of the blocks into x */
	long        Gcount;     /* G evaluations */
	double      neChange, chemChange;
	double      resHist[RES_HIST_MAX]; /* residuals of FIXED_POINT (N-3) */
	int         NresHist;

	/* --- Variables u of the solvers (setupVariables) -- ---------- */

	int                  *packBlock;  /* block of each packed block */
	int                  *carryBlock; /* the other blocks: carried */
	int                   Ncarry;
	long                 *packStart; /* Npack + 1 offsets into u */
	int                   Npack;
	long                  Npacked;  /* length of u */
	long                  offJ;     /* offset of J in x */
	double               *varScale; /* RELATIVE: fixed scale per entry of u */
	enum solver_variables varKind;
} SolverState;

#define SOLVER_STATE_INIT {.offJ = -1, .varKind = VARIABLES_LOG}

void iterateFixedPoint(RHState *rh, int NmaxIter, double iterLimit);
void iterateAnderson(RHState *rh, int NmaxIter, double iterLimit);
void iterateNewton(RHState *rh, int NmaxIter, double iterLimit);

/* --- Acceptance checks of G (gchecks.c), run by the gcheck programs */

typedef struct {
	int    Nsteps, lag, Necn;
	double hEcnoise, fdmax;
	bool_t fd;
	char  *perturb;
} GCheckOptions;

void gcheckParse(int *argc, char *argv[], GCheckOptions *opt);
void gcheckRun(RHState *rh, const GCheckOptions *opt);

#endif /* !__SOLVER_H__ */

/* ------- end ---------------------------- solver.h ---------------- */
