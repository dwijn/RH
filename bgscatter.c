/* ------- file: -------------------------- bgscatter.c -------------

	   Background scattering at one wavelength, with the populations
	   fixed (N-51).

	   --------------------------                      ----------RH-- */

/* --- The background source function holds sca_c J, scattering by
	   electrons, Rayleigh scattering and, with RLK_SCATTER, the
	   scattering part (1 - eps) of every Kurucz line. Formal() used
	   the J of the previous formal solution for it: a plain lambda
	   iteration, which needs of order 1/eps sweeps where a line with
	   eps = 1e-4 .. 1e-6 is optically thick. Nothing else converges
	   it: the operator of the active atoms leaves it out, and the
	   outer solver sees only one lambda step per iteration.

	   With the populations fixed, J at one wavelength solves the
	   linear problem J = Lambda[(eta + sca_c J) / chi]. From a formal
	   solution J_fs from u, the Jacobi step with the diagonal operator

		 L_k = sum over rays of w Psi_k sca_c,k / chi_k   (dJ_k / du_k)

	   is g(u) = u + (J_fs(u) - u) / (1 - L) (or the step with the
	   tridiagonal of dJ/du). Its fixed point is that of the lambda
	   iteration. bgScatterIterate repeats it, with the caller's sweep
	   function for the formal solutions, accelerated by Anderson
	   (anderson.c, depth BACKGROUND_ANDERSON_DEPTH, 0 for plain
	   Jacobi), until the largest |g - u| / g is below
	   BACKGROUND_ITER_LIMIT or after BACKGROUND_N_MAX_ITER sweeps.

	   The iteration runs before the main sweep of a wavelength
	   (bgScatterBefore), from the J of the previous sweep, and the
	   main sweep starts from the iterated J*. J, the rates of the
	   active atoms, Jgas of hybrid PRD and the emergent intensity then
	   all come from that one formal solution, J = Lambda[S(J*)], which
	   differs from J* by about BACKGROUND_ITER_LIMIT. Every wavelength
	   is iterated first in the first sweep, and after that while the
	   iteration still does something: after a main sweep from the
	   plain J (bgScatterAfter), where the Jacobi step would change J by
	   more than BACKGROUND_ITER_LIMIT times the scattering albedo, the
	   next sweep iterates first; where the iteration from J needed no
	   step, the next sweep does not. The cost is one more formal
	   solution per wavelength iterated. BACKGROUND_N_MAX_ITER = 0
	   gives the plain lambda iteration of before.

	   L is clamped to [0, L_MAX], and a point where the Jacobi step or
	   the extrapolation would make J negative takes the plain lambda
	   value J_fs, or g, instead.

	   Memory: per wavelength being iterated (and per thread), about
	   3 m + 9 vectors of Nspace doubles, freed at the end --
	   --                                              -------------- */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "spectrum.h"
#include "anderson.h"
#include "error.h"
#include "bgscatter.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern InputData input;
extern Spectrum  spectrum;
extern char      messageStr[];

#define L_MAX    0.9999
#define COND_MAX 1.0E+08 /* as ANDERSON_COND_MAX */

static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static BgScatterStats  stats;   /* since the last bgScatterStats(.., TRUE) */
static unsigned char  *pending; /* per wavelength: iterate before the sweep */

/* ------- begin -------------------------- bgScatterActive.c ------- */

int bgScatterActive(void)
{
	return input.bg_NmaxIter > 0;
}
/* ------- end ---------------------------- bgScatterActive.c ------- */

/* ------- begin -------------------------- bgScatterSetup.c -------- */

void bgScatterSetup(void)
{
	const char  routineName[] = "bgScatterSetup";
	const char *why;

	/* --- Once, before the first sweep (initScatter): the iteration
		   runs only where the geometry implements it; elsewhere the
		   background scattering keeps the lambda iteration, with a
		   warning if BACKGROUND_N_MAX_ITER was set -- -------------- */

	if (input.bg_NmaxIter == 0) return;
	if ((why = bgScatterSupported()) != NULL) {
		if (input.bg_NmaxIter_set) {
			sprintf(messageStr,
					"BACKGROUND_N_MAX_ITER = %d is ignored: %s.\n"
					" The background scattering is lambda iterated.",
					input.bg_NmaxIter, why);
			Error(WARNING, routineName, messageStr);
		}
		input.bg_NmaxIter = 0;
		return;
	}
	sprintf(messageStr,
			" Background scattering iterated at each wavelength (N-51):"
			" BACKGROUND_OPERATOR = %s,\n  BACKGROUND_ANDERSON_DEPTH = %d,"
			" BACKGROUND_ITER_LIMIT = %.2E, BACKGROUND_N_MAX_ITER = %d\n",
			(input.bg_operator == BG_OPERATOR_TRIDIAGONAL) ? "TRIDIAGONAL"
														   : "DIAGONAL",
			input.bg_depth, input.bg_iterLimit, input.bg_NmaxIter);
	Error(MESSAGE, routineName, messageStr);

	/* --- The first sweep iterates every wavelength first -- ------ */

	if (pending == NULL) {
		pending = (unsigned char *) malloc(spectrum.Nspect);
		memset(pending, 1, spectrum.Nspect);
	}
}
/* ------- end ---------------------------- bgScatterSetup.c -------- */

/* ------- begin -------------------------- tridiagSolve.c ---------- */

static int tridiagSolve(int N, double *lo, double *d, double *up, double *x)
{
	register int k;

	double *up2, f, t;

	/* --- Solve the tridiagonal system (lo[k] x_k-1 + d[k] x_k + up[k]
		   x_k+1 = b_k, b in x on entry) by Gaussian elimination with
		   partial pivoting (as LAPACK dgtsv). Overwrites lo, d, up.
		   Returns FALSE when singular or not finite -- ------------ */

	up2 = (double *) calloc(N, sizeof(double));
	for (k = 0; k < N - 1; k++) {
		if (fabs(d[k]) >= fabs(lo[k + 1])) {
			if (d[k] == 0.0) {
				free(up2);
				return FALSE;
			}
			f = lo[k + 1] / d[k];
			d[k + 1] -= f * up[k];
			x[k + 1] -= f * x[k];
			up2[k] = 0.0;
		} else {
			/* --- Swap rows k and k + 1 --           -------------- */

			f        = d[k] / lo[k + 1];
			d[k]     = lo[k + 1];
			t        = d[k + 1];
			d[k + 1] = up[k] - f * t;
			up[k]    = t;
			up2[k]   = (k < N - 2) ? up[k + 1] : 0.0;
			if (k < N - 2) up[k + 1] = -f * up2[k];
			t        = x[k];
			x[k]     = x[k + 1];
			x[k + 1] = t - f * x[k + 1];
		}
	}
	if (d[N - 1] == 0.0) {
		free(up2);
		return FALSE;
	}
	x[N - 1] /= d[N - 1];
	if (N > 1) x[N - 2] = (x[N - 2] - up[N - 2] * x[N - 1]) / d[N - 2];
	for (k = N - 3; k >= 0; k--)
		x[k] = (x[k] - up[k] * x[k + 1] - up2[k] * x[k + 2]) / d[k];
	free(up2);
	for (k = 0; k < N; k++)
		if (!isfinite(x[k])) return FALSE;
	return TRUE;
}
/* ------- end ---------------------------- tridiagSolve.c ---------- */

/* ------- begin -------------------------- precondStep.c ----------- */

static void precondStep(int N, const double *u, const double *Jfs,
						const double *A, double *g, double *work)
{
	register int k;

	double l, *lo = work, *d = work + N, *up = work + 2 * N, *x = work + 3 * N;

	/* --- g = u + M^-1 (Jfs - u): M = I - A with A the tridiagonal of
		   dJ/du (BACKGROUND_OPERATOR = TRIDIAGONAL), or its diagonal
		   (DIAGONAL: Jacobi, A clamped to [0, L_MAX]). Where the
		   tridiagonal system is singular, the Jacobi step. A point that
		   would get J < 0 takes the lambda value Jfs -- ---------- */

	if (input.bg_operator == BG_OPERATOR_TRIDIAGONAL) {
		for (k = 0; k < N; k++) {
			lo[k] = (k > 0) ? -A[3 * k] : 0.0;
			d[k]  = 1.0 - A[3 * k + 1];
			up[k] = (k < N - 1) ? -A[3 * k + 2] : 0.0;
			x[k]  = Jfs[k] - u[k];
		}
		if (tridiagSolve(N, lo, d, up, x)) {
			for (k = 0; k < N; k++) {
				g[k] = u[k] + x[k];
				if (!(g[k] >= 0.0)) g[k] = Jfs[k];
			}
			return;
		}
	}
	for (k = 0; k < N; k++) {
		l = A[3 * k + 1];
		if (!(l > 0.0)) l = 0.0;
		if (l > L_MAX) l = L_MAX;
		g[k] = u[k] + (Jfs[k] - u[k]) / (1.0 - l);
		if (!(g[k] >= 0.0)) g[k] = Jfs[k];
	}
}
/* ------- end ---------------------------- precondStep.c ----------- */

/* ------- begin -------------------------- relChange.c ------------- */

static double relChange(int N, const double *u, const double *g,
						const double *albedo)
{
	register int k;

	double r = 0.0, d;

	/* --- Largest |g - u| / g, weighted by albedo if given. g = 0 where
		   every emissivity underflowed (far EUV) counts as no change,
		   as in Formal -- ---                          -------------- */

	for (k = 0; k < N; k++) {
		if (g[k] != 0.0)
			d = fabs(1.0 - u[k] / g[k]);
		else
			d = (u[k] != 0.0) ? 1.0 : 0.0;
		if (albedo) d *= albedo[k];
		if (d > r) r = d;
	}
	return r;
}
/* ------- end ---------------------------- relChange.c ------------- */

/* ------- begin -------------------------- bgScatterIterate.c ------ */

int bgScatterIterate(int Nspace, const double *Jdag, const double *Jfs,
					 const double *A, const double *albedo, double *J,
					 BgSweep sweep, void *ctx)
{
	register int k;

	int                 m = input.bg_depth, niter = 0, converged = FALSE;
	long                start[2] = {0, Nspace};
	double             *u, *g, *un, *Jn, *An, *w, *gamma, *work, r;
	enum AndersonAction action;
	AndersonOptions     opt;
	AndersonSolve       info;
	Anderson           *a = NULL;

	/* --- Jdag: the J the main sweep started from; Jfs: its formal
		   solution; A: the tridiagonal of dJ/du, A[3k + i] = dJ_k/
		   du_k-1+i; albedo: the angle average of sca_c / chi. J (may
		   be Jfs) gets the result -- --                  -------------- */

	u    = (double *) malloc(6 * Nspace * sizeof(double));
	g    = u + Nspace;
	work = g + Nspace;
	memcpy(u, Jdag, Nspace * sizeof(double));
	precondStep(Nspace, u, Jfs, A, g, work);

	/* --- Nothing to iterate where the scattering cannot change J by
		   more than the limit --                       -------------- */

	if (relChange(Nspace, u, g, albedo) <= input.bg_iterLimit) {
		memcpy(J, g, Nspace * sizeof(double));
		free(u);
		return 0;
	}
	un    = (double *) malloc(6 * Nspace * sizeof(double));
	Jn    = un + Nspace;
	An    = Jn + Nspace;
	w     = An + 3 * Nspace;
	gamma = (double *) calloc((m > 0) ? m : 1, sizeof(double));
	for (k = 0; k < Nspace; k++) w[k] = (g[k] > 0.0) ? 1.0 / g[k] : 0.0;
	if (m > 0) {
		a             = andersonInit(Nspace, m);
		opt.Nblock    = 1;
		opt.start     = start;
		opt.per_block = FALSE;
		opt.cond_max  = COND_MAX;
		opt.safeguard = 1.0;
		opt.window    = 1; /* the newest accepted pair, as before S-12 */
	}
	r = relChange(Nspace, u, g, NULL);
	for (;;) {
		if (r <= input.bg_iterLimit) {
			converged = TRUE;
			break;
		}
		if (niter >= input.bg_NmaxIter) break;

		/* --- The next point: Anderson, or the Jacobi step. An
			   extrapolated J that is not positive takes g -- ----- */

		if (a) {
			action = andersonStep(a, &opt, u, g, NULL, w, TRUE, un, gamma, &info);
			if (action == ANDERSON_EXTRAPOLATED) {
				for (k = 0; k < Nspace; k++)
					if (!(un[k] > 0.0)) un[k] = g[k];
			}
		} else
			memcpy(un, g, Nspace * sizeof(double));
		memcpy(u, un, Nspace * sizeof(double));

		sweep(ctx, u, Jn, An, NULL);
		precondStep(Nspace, u, Jn, An, g, work);
		niter++;
		r = relChange(Nspace, u, g, NULL);
	}
	memcpy(J, g, Nspace * sizeof(double));

	pthread_mutex_lock(&lock);
	stats.Niterated++;
	stats.Nsweeps += niter;
	if (!converged) {
		stats.Nunconverged++;
		if (r > stats.rmax) stats.rmax = r;
	}
	pthread_mutex_unlock(&lock);

	andersonFree(a);
	free(gamma);
	free(un);
	free(u);
	return niter;
}
/* ------- end ---------------------------- bgScatterIterate.c ------ */

/* ------- begin -------------------------- bgScatterPending.c ------ */

bool_t bgScatterPending(int nspect)
{
	return pending != NULL && pending[nspect];
}
/* ------- end ---------------------------- bgScatterPending.c ------ */

/* ------- begin -------------------------- bgScatterBefore.c ------- */

void bgScatterBefore(int nspect, int Nspace, const double *Jdag,
					 double *Jstar, BgSweep sweep, void *ctx)
{
	register int k;

	int     niter;
	double *Jfs, *A, *albedo;

	/* --- One formal solution from Jdag, the iteration from there to
		   J*, the source of the main sweep. A wavelength the scattering
		   no longer changes by more than the limit is not iterated
		   first in the next sweep (J* is then the Jacobi step g). A
		   wavelength is solved by one thread at a time, so its flag
		   needs no lock --                              -------------- */

	Jfs    = (double *) malloc(5 * Nspace * sizeof(double));
	A      = Jfs + Nspace;
	albedo = A + 3 * Nspace;
	sweep(ctx, Jdag, Jfs, A, albedo);
	niter = bgScatterIterate(Nspace, Jdag, Jfs, A, albedo, Jstar, sweep, ctx);
	if (niter == 0) pending[nspect] = 0;

	/* --- Opacities negative enough to overflow are caught by the main
		   sweep (N-40); keep the source finite until then -- ------- */

	for (k = 0; k < Nspace; k++) {
		if (!isfinite(Jstar[k])) {
			memcpy(Jstar, Jdag, Nspace * sizeof(double));
			break;
		}
	}
	pthread_mutex_lock(&lock);
	stats.Nsweeps++;
	pthread_mutex_unlock(&lock);
	free(Jfs);
}
/* ------- end ---------------------------- bgScatterBefore.c ------- */

/* ------- begin -------------------------- bgScatterAfter.c -------- */

void bgScatterAfter(int nspect, int Nspace, const double *u,
					const double *Jfs, const double *A, const double *albedo)
{
	double *g;

	if (pending == NULL) return;
	g = (double *) malloc(5 * Nspace * sizeof(double)); /* g, work */
	precondStep(Nspace, u, Jfs, A, g, g + Nspace);
	if (relChange(Nspace, u, g, albedo) > input.bg_iterLimit)
		pending[nspect] = 1;
	free(g);
}
/* ------- end ---------------------------- bgScatterAfter.c -------- */

/* ------- begin -------------------------- bgScatterStats.c -------- */

void bgScatterStats(BgScatterStats *s, int reset)
{
	pthread_mutex_lock(&lock);
	*s = stats;
	if (reset) memset(&stats, 0, sizeof(stats));
	pthread_mutex_unlock(&lock);
}
/* ------- end ---------------------------- bgScatterStats.c -------- */

/* ------- end ---------------------------- bgscatter.c ------------- */
