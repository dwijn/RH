/* ------- file: -------------------------- solver.h ----------------

	   The map G of docs/SOLVER.md (S-2) and the solvers built on it.

	   --------------------------                      ----------RH-- */

#ifndef __SOLVER_H__
#define __SOLVER_H__

/* --- State vector x (populations, molecular populations, PRD rho, J)
	   and the map G, see solver.c --                 -------------- */

long        stateInit(void);
void        stateGather(double *x);
void        stateScatter(const double *x);
int         stateNblocks(void);
const char *stateBlockName(int n);
long        stateBlockSize(int n);
void        G_eval(const double *x_in, double *x_out);
long        G_count(void);

void iterateFixedPoint(int NmaxIter, double iterLimit);

/* --- Acceptance checks of G (gchecks.c), run by the gcheck programs */

typedef struct {
	int    Nsteps, lag, Necn;
	double hEcnoise, fdmax;
	bool_t fd;
	char  *perturb;
} GCheckOptions;

void gcheckParse(int *argc, char *argv[], GCheckOptions *opt);
void gcheckRun(const GCheckOptions *opt);

#endif /* !__SOLVER_H__ */

/* ------- end ---------------------------- solver.h ---------------- */
