/* ------- file: -------------------------- solver.h ----------------

	   The map G of docs/SOLVER.md (S-2) and the solvers built on it.

	   --------------------------                      ----------RH-- */

#ifndef __SOLVER_H__
#define __SOLVER_H__

/* --- State vector x (populations, molecular populations, PRD rho, J)
	   and the map G, see solver.c --                 -------------- */

long stateInit(void);
void stateGather(double *x);
void stateScatter(const double *x);
void G_eval(const double *x_in, double *x_out);
long G_count(void);

void iterateFixedPoint(int NmaxIter, double iterLimit);

#endif /* !__SOLVER_H__ */

/* ------- end ---------------------------- solver.h ---------------- */
