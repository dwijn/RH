/* ------- file: -------------------------- bgscatter.h -------------

	   Background scattering at one wavelength with the populations
	   fixed: Jacobi (diagonal operator) with Anderson acceleration
	   (N-51), see bgscatter.c.

	   --------------------------                      ----------RH-- */

#ifndef __BGSCATTER_H__
#define __BGSCATTER_H__

/* --- One formal solution of all rays at the wavelength, with the
	   background scattering source sca_c u: returns the new mean
	   intensity J and the tridiagonal of dJ/du, A[3k + i] =
	   dJ_k / du_k-1+i (the angle sum of w (dI_k/dS_j) sca_c,j /
	   chi_j), without side effects on rates or output -- ---- */

typedef void (*BgSweep)(void *ctx, const double *u, double *J, double *A,
						double *albedo);

/* --- Per-sweep totals over all wavelengths (main thread, between
	   sweeps) --                                      -------------- */

typedef struct {
	long   Niterated;    /* wavelengths that needed the iteration */
	long   Nsweeps;      /* formal solutions of a wavelength it took */
	long   Nunconverged; /* stopped at BACKGROUND_N_MAX_ITER */
	double rmax;         /* largest final residual of those */
	long   Nto, Nback;   /* switches to the linear solver and back */
	long   Nlinear;      /* wavelengths solved linear now (N-52;
							not reset) */
} BgScatterStats;

int  bgScatterActive(void);
void bgScatterSetup(void);

/* --- Defined by each geometry (its formal.c): NULL where the
	   iteration is implemented for the current input, otherwise why
	   not --                                          -------------- */

const char *bgScatterSupported(void);
int         bgScatterIterate(int Nspace, const double *Jdag, const double *Jfs,
							 const double *A, const double *albedo, double *J,
							 BgSweep sweep, void *ctx);
void        bgScatterStats(BgScatterStats *s, int reset);

/* --- Iterate before the main sweep (N-51): bgScatterPending says
	   whether wavelength nspect is to be iterated first; then
	   bgScatterBefore gives the source J* of the main sweep (from Jdag,
	   with one more formal solution, sweep). Otherwise, after the main
	   sweep, bgScatterAfter decides from its J (Jfs, from u), A and
	   albedo whether the next sweep iterates first -- ----------- */

bool_t bgScatterPending(int nspect);
void   bgScatterBefore(int nspect, int Nspace, const double *Jdag,
					   double *Jstar, BgSweep sweep, void *ctx);
void   bgScatterAfter(int nspect, int Nspace, const double *u,
					  const double *Jfs, const double *A, const double *albedo);
/* --- With S_BEZIER3, whether wavelength nspect is solved with the
	   linear solver (N-52): where albedo exceeds
	   BACKGROUND_LINEAR_ALBEDO; albedo < 0 keeps the last choice - */

bool_t bgScatterLinear(int nspect, double albedo);

#endif /* !__BGSCATTER_H__ */

/* ------- end ---------------------------- bgscatter.h ------------- */
