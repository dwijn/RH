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
	double rmax;         /* largest estimated error of those (N-66) */
	long   Nto, Nback;   /* switches to the linear solver and back */
	long   Nlinear;      /* wavelengths solved linear now (N-52;
							not reset) */
} BgScatterStats;

/* --- Whether the geometry has the iteration, and whether it falls back
	   to the linear solver (N-52), it says in geometryOps
	   (bg_supported, bg_linear_fallback; geometryops.h) -- ------- */

int  bgScatterActive(void);
void bgScatterSetup(void);
int  bgScatterIterate(int Nspace, const double *Jdag, const double *Jfs,
					  const double *A, const double *albedo, double *J,
					  double *kappa, BgSweep sweep, void *ctx);
void bgScatterStats(BgScatterStats *s, int reset);

/* --- The state of the background scattering between calls, its part
	   of RHState (rh->bgscatter, Q-46). The per-wavelength arrays are
	   made by bgScatterSetup; one thread works on a wavelength at a
	   time, so only the counters need the mutex -- -------------- */

typedef struct {
	pthread_mutex_t lock;    /* guards stats and Nlinear */
	BgScatterStats  stats;   /* since the last bgScatterStats(.., TRUE) */
	unsigned char  *pending; /* per wavelength: iterate before the sweep */
	unsigned char  *linear;  /* per wavelength, with S_BEZIER3 (N-52) */
	float          *kappa;   /* per wavelength: the last kappa, 0: none */
	long            Nlinear;
} BgScatterState;

#define BGSCATTER_STATE_INIT {.lock = PTHREAD_MUTEX_INITIALIZER}

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
