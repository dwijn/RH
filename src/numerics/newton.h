/* ------- file: -------------------------- newton.h ----------------

	   Matrix-free GMRES in a weighted, masked inner product and the
	   Eisenstat-Walker forcing term, for the inexact Newton solver
	   of S-4 in docs/SOLVER.md. See newton.c.

	   --------------------------                      ----------RH-- */

#ifndef __NEWTON_H__
#define __NEWTON_H__

/* --- Constants of the forcing term (Eisenstat & Walker 1996,
	   choice 2) and of the backtracking line search --  ------------ */

#define NEWTON_EW_GAMMA    0.9
#define NEWTON_EW_ALPHA    2.0
#define NEWTON_EW_SAFE     0.1 /* their safeguard applies above this */
#define NEWTON_ARMIJO      1.0E-4
#define NEWTON_N_BACKTRACK 4 /* lambda = 1, 1/2, ..., 1/16 */

/* --- The product y = A v of the linear system. v has unit weighted
	   norm and is zero outside the mask; y is read on the mask only -- */

typedef void (*GmresProduct)(void *ctx, const double *v, double *y);

typedef struct {
	long     N; /* length of the vectors */
	int      k; /* at most k iterations */
	double **V; /* [k + 1][N]: Krylov basis */
	double  *H; /* [(k + 1) * k]: Hessenberg matrix, column major */
	double  *cs, *sn, *g, *y, *h;
	double  *resid; /* [k + 1]: residual norm after each iteration */
} Gmres;

typedef struct {
	int    Niter;     /* products made */
	double rnorm0;    /* ||b||_W */
	double rnorm;     /* ||b - A x||_W, from the rotations */
	int    breakdown; /* the Krylov space became invariant */
} GmresResult;

Gmres      *gmresInit(long N, int k);
void        gmresFree(Gmres *gm);
double      weightedNorm(long N, const double *v, const unsigned char *mask,
						 const double *w);
GmresResult gmresSolve(Gmres *gm, GmresProduct A, void *ctx, const double *b,
					   const unsigned char *mask, const double *w, int kmax,
					   double tol, double *x);
double      newtonForcing(int first, double eta_prev, double rnorm,
						  double rnorm_prev, double eta_max, double floor_abs);

#endif /* !__NEWTON_H__ */

/* ------- end ---------------------------- newton.h ---------------- */
