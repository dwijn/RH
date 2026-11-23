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
	   Jacobi), until the error of J estimated from the step is below
	   BACKGROUND_ITER_LIMIT (N-66, below) or after
	   BACKGROUND_N_MAX_ITER sweeps.

	   The iteration runs before the main sweep of a wavelength
	   (bgScatterBefore), from the J of the previous sweep, and the
	   main sweep starts from the iterated J*. J, the rates of the
	   active atoms, Jgas of hybrid PRD and the emergent intensity then
	   all come from that one formal solution, J = Lambda[S(J*)], which
	   is within about BACKGROUND_ITER_LIMIT of the solution. Every
	   wavelength is iterated first in the first sweep, and after that
	   while the iteration still does something: after a main sweep from
	   the plain J (bgScatterAfter), where the error estimated from the
	   Jacobi step exceeds BACKGROUND_ITER_LIMIT, the next sweep
	   iterates first; where the iteration from J needed no step, the
	   next sweep does not. The cost is one more formal
	   solution per wavelength iterated. BACKGROUND_N_MAX_ITER = 0
	   gives the plain lambda iteration of before.

	   L is clamped to [0, L_MAX], and a point where the Jacobi step or
	   the extrapolation would make J negative takes the plain lambda
	   value J_fs, or g, instead.

	   Memory: per wavelength being iterated (and per thread), about
	   3 m + 9 + 4 KAPPA_PAIRS vectors of Nspace doubles, freed at the
	   end --                                          -------------- */

/* --- The stopping test (N-66). A small step is not a small error:
	   with M the Jacobian of the step map g(u), the error of g is
	   M (I - M)^-1 (g - u), of order rho / (1 - rho) times the step,
	   rho the spectral radius of M. In the core of a strong scattering
	   line rho is up to 0.999, and a test on the step alone left
	   errors of 170-400 times the limit. The iteration therefore stops
	   where F max |g - u| / g <= BACKGROUND_ITER_LIMIT, with

		 F = min(max(kappa, 1), a / (1 - a)):

	   kappa, the secant estimate of ||(I - M)^-1||, from the last
	   KAPPA_PAIRS evaluated pairs (u_i, g_i). At fixed populations the
	   map is affine, so the differences du_i = u_i - u_new and df_i =
	   (g_i - u_i) - (g_new - u_new) obey df_i = (M - I) du_i, and
	   1 / kappa is the smallest singular value of M - I on the span of
	   the du_i, in the weights of the iteration (1 / J at its start).
	   It underestimates the norm where that span misses the slowest
	   mode; Anderson's own history (m + 1 pairs) is too short for it.
	   a, the largest albedo: the iteration solves with the linear
	   short characteristics (also with S_BEZIER3, N-52), where
	   Lambda sigma >= 0 and the step is a regular splitting of
	   I - Lambda sigma, so that rho <= ||Lambda sigma||_inf <= a.
	   kappa is computed only once the step is below the limit, and is
	   kept per wavelength as the largest so far: a warm start takes a
	   few steps only, whose span misses the slowest mode. The first
	   check of bgScatterIterate and bgScatterAfter use it, or a / (1 -
	   a) where there is none yet. MURaM column of N-65 with the populations frozen: the
	   largest error of J is then 0.4-0.6 times the limit for limits
	   from 1e-3 to 1e-8 (the step alone: 170-400 times), for 12 %
	   more formal solutions in the first sweep (N-66 evidence) -- */

/* --- The operator of S_BEZIER3 is not positive (N-52): where the
	   Steffen limiter takes the downwind slope, I_k has a negative
	   weight on S_k+1, and with the velocity field the up- and
	   downgoing rays no longer cancel it in J. Its spectral radius at
	   strongly scattering wavelengths exceeds 1 and the iteration
	   diverges there. With S_BEZIER3, a wavelength is therefore solved
	   with the linear solver, for the whole formal solution, where its
	   largest albedo below vertical optical depth BACKGROUND_LINEAR_TAU
	   exceeds BACKGROUND_LINEAR_ALBEDO (bgScatterLinear; the geometry
	   decides). The choice is made afresh in every sweep that updates
	   J, as the opacities of the active atoms change; a sweep that does
	   not keeps the last one. Only those wavelengths iterate the
	   background scattering; the others keep the cubic Bezier solver
	   and the lambda iteration, which converges where the albedo is
	   small. In a MURaM column the switches back to Bezier die out
	   within the first 20 or so iterations, and a restart picks the
	   same wavelengths (N-52 evidence) -- --       -------------- */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
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

#define L_MAX       0.9999
#define COND_MAX    1.0E+08 /* as ANDERSON_COND_MAX */
#define KAPPA_PAIRS 10      /* evaluated pairs for kappa (N-66) */

static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
static BgScatterStats  stats;   /* since the last bgScatterStats(.., TRUE) */
static unsigned char  *pending; /* per wavelength: iterate before the sweep */
static unsigned char  *linear;  /* per wavelength, with S_BEZIER3 (N-52) */
static float          *kappa;   /* per wavelength: the last kappa, 0: none */
static long            Nlinear;

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
		kappa = (float *) calloc(spectrum.Nspect, sizeof(float));
	}

	if (input.S_interpolation != S_LINEAR && linear == NULL) {
		linear = (unsigned char *) calloc(spectrum.Nspect, 1);
		sprintf(messageStr,
				"  with S_INTERPOLATION = S_BEZIER3: wavelengths whose albedo"
				" exceeds BACKGROUND_LINEAR_ALBEDO = %.2f\n  below vertical"
				" optical depth BACKGROUND_LINEAR_TAU = %.2f use S_LINEAR"
				" (N-52)\n",
				input.bg_linearAlbedo, input.bg_linearTau);
		Error(MESSAGE, routineName, messageStr);
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

static double relChange(int N, const double *u, const double *g)
{
	register int k;

	double r = 0.0, d;

	/* --- Largest |g - u| / g. g = 0 where
		   every emissivity underflowed (far EUV) counts as no change,
		   as in Formal -- ---                          -------------- */

	for (k = 0; k < N; k++) {
		if (g[k] != 0.0)
			d = fabs(1.0 - u[k] / g[k]);
		else
			d = (u[k] != 0.0) ? 1.0 : 0.0;
		if (d > r) r = d;
	}
	return r;
}
/* ------- end ---------------------------- relChange.c ------------- */

/* ------- begin -------------------------- maxAlbedo.c ------------- */

static double maxAlbedo(int N, const double *albedo)
{
	register int k;

	double a = 0.0;

	for (k = 0; k < N; k++)
		if (albedo[k] > a) a = albedo[k];
	return a;
}
/* ------- end ---------------------------- maxAlbedo.c ------------- */

/* ------- begin -------------------------- errFactor.c ------------- */

static double errFactor(double kap, double amax)
{
	double bound = amax / fmax(1.0 - amax, DBL_MIN); /* no select (SIGFPE, N-37) */

	/* --- F of the stopping test: the error of J is about F times the
		   relative step. kap <= 0: no estimate yet -- -------------- */

	if (!(kap > 0.0)) return bound;
	return fmin(fmax(kap, 1.0), bound);
}
/* ------- end ---------------------------- errFactor.c ------------- */

/* ------- begin -------------------------- secantKappa.c ----------- */

static double secantKappa(int N, int Npairs, int inew, const double *hu,
						  const double *hf, const double *w, double *work)
{
	register int k;

	int    i, j, l, q = 0, it;
	double C[KAPPA_PAIRS * KAPPA_PAIRS], L[KAPPA_PAIRS * KAPPA_PAIRS],
		x[KAPPA_PAIRS], y[KAPPA_PAIRS], c, n0, nrm, lam = 0.0, *Q = work,
													*B = work + KAPPA_PAIRS * N, *d, *e;

	/* --- The estimate kappa of ||(I - M)^-1|| from the pairs held,
		   hu[i N ...] the iterate u_i and hf[i N ...] its step g_i -
		   u_i, i = inew the newest (see the stopping test above).
		   Gram-Schmidt on the weighted du_i (Q), with the same column
		   operations on the df_i (B = (M - I) Q), dropping a du_i
		   that is (nearly) in the span of the earlier ones; then the
		   smallest eigenvalue of B^T B by inverse iteration. Returns
		   0 without differences, and a huge kappa where B^T B is not
		   positive definite in floating point. work: 2 KAPPA_PAIRS N
		   doubles --                                   -------------- */

	for (i = 0; i < Npairs; i++) {
		if (i == inew) continue;
		d  = Q + (long) q * N;
		e  = B + (long) q * N;
		n0 = 0.0;
		for (k = 0; k < N; k++) {
			d[k] = (hu[(long) i * N + k] - hu[(long) inew * N + k]) * w[k];
			e[k] = (hf[(long) i * N + k] - hf[(long) inew * N + k]) * w[k];
			n0 += d[k] * d[k];
		}
		if (!(n0 > 0.0)) continue;
		for (l = 0; l < q; l++) {
			c = 0.0;
			for (k = 0; k < N; k++) c += Q[(long) l * N + k] * d[k];
			for (k = 0; k < N; k++) {
				d[k] -= c * Q[(long) l * N + k];
				e[k] -= c * B[(long) l * N + k];
			}
		}
		nrm = 0.0;
		for (k = 0; k < N; k++) nrm += d[k] * d[k];
		if (!(nrm > 1.0E-16 * n0)) continue;
		nrm = sqrt(nrm);
		for (k = 0; k < N; k++) {
			d[k] /= nrm;
			e[k] /= nrm;
		}
		q++;
	}
	if (q == 0) return 0.0;

	for (i = 0; i < q; i++) {
		for (j = 0; j <= i; j++) {
			c = 0.0;
			for (k = 0; k < N; k++) c += B[(long) i * N + k] * B[(long) j * N + k];
			C[i * q + j] = C[j * q + i] = c;
		}
	}
	for (j = 0; j < q; j++) {
		c = C[j * q + j];
		for (l = 0; l < j; l++) c -= L[j * q + l] * L[j * q + l];
		if (!(c > 0.0)) return 1.0E+30;
		L[j * q + j] = sqrt(c);
		for (i = j + 1; i < q; i++) {
			c = C[i * q + j];
			for (l = 0; l < j; l++) c -= L[i * q + l] * L[j * q + l];
			L[i * q + j] = c / L[j * q + j];
		}
	}
	for (i = 0; i < q; i++) x[i] = 1.0;
	for (it = 0; it < 50; it++) {
		for (i = 0; i < q; i++) {
			c = x[i];
			for (l = 0; l < i; l++) c -= L[i * q + l] * y[l];
			y[i] = c / L[i * q + i];
		}
		for (i = q - 1; i >= 0; i--) {
			c = y[i];
			for (l = i + 1; l < q; l++) c -= L[l * q + i] * x[l];
			x[i] = c / L[i * q + i];
		}
		nrm = 0.0;
		for (i = 0; i < q; i++) nrm += x[i] * x[i];
		nrm = sqrt(nrm);
		for (i = 0; i < q; i++) x[i] /= nrm;
		lam = 1.0 / nrm;
	}
	return 1.0 / sqrt(fmax(lam, 1.0E-60));
}
/* ------- end ---------------------------- secantKappa.c ----------- */

/* ------- begin -------------------------- pushPair.c -------------- */

static void pushPair(int N, int *Npairs, int *inew, const double *u,
					 const double *g, double *hu, double *hf)
{
	register int k;

	/* --- A ring of the last KAPPA_PAIRS evaluated pairs -- ------- */

	*inew = (*Npairs == 0) ? 0 : (*inew + 1) % KAPPA_PAIRS;
	if (*Npairs < KAPPA_PAIRS) (*Npairs)++;
	for (k = 0; k < N; k++) {
		hu[(long) *inew * N + k] = u[k];
		hf[(long) *inew * N + k] = g[k] - u[k];
	}
}
/* ------- end ---------------------------- pushPair.c -------------- */

/* ------- begin -------------------------- bgScatterIterate.c ------ */

int bgScatterIterate(int Nspace, const double *Jdag, const double *Jfs,
					 const double *A, const double *albedo, double *J,
					 double *kap, BgSweep sweep, void *ctx)
{
	register int k;

	int m = input.bg_depth, niter = 0, converged = FALSE,
		Npairs = 0, inew = 0, fresh = FALSE;
	long    start[2] = {0, Nspace};
	double *u, *g, *un, *Jn, *An, *w, *gamma, *work, *hu, *hf, r,
		amax, kn = 0.0;
	enum AndersonAction action;
	AndersonOptions     opt;
	AndersonSolve       info;
	Anderson           *a = NULL;

	/* --- Jdag: the J the main sweep started from; Jfs: its formal
		   solution; A: the tridiagonal of dJ/du, A[3k + i] = dJ_k/
		   du_k-1+i; albedo: the angle average of sca_c / chi. J (may
		   be Jfs) gets the result. kap: in, the last kappa of this
		   wavelength (<= 0: none); out, the new one (N-66) -- ---- */

	u    = (double *) malloc(6 * Nspace * sizeof(double));
	g    = u + Nspace;
	work = g + Nspace;
	memcpy(u, Jdag, Nspace * sizeof(double));
	precondStep(Nspace, u, Jfs, A, g, work);
	amax = maxAlbedo(Nspace, albedo);
	r    = relChange(Nspace, u, g);

	/* --- Nothing to iterate where the error estimated from the step,
		   with the last kappa, is below the limit -- -------------- */

	if (errFactor(*kap, amax) * r <= input.bg_iterLimit) {
		memcpy(J, g, Nspace * sizeof(double));
		free(u);
		return 0;
	}
	un    = (double *) malloc((6 + 4 * KAPPA_PAIRS) * Nspace * sizeof(double));
	Jn    = un + Nspace;
	An    = Jn + Nspace;
	w     = An + 3 * Nspace;
	hu    = w + Nspace;
	hf    = hu + KAPPA_PAIRS * Nspace;
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
	pushPair(Nspace, &Npairs, &inew, u, g, hu, hf);
	for (;;) {

		/* --- Converged where the step is below the limit and the
			   error estimated from it too (kappa only then: work in
			   hf + KAPPA_PAIRS Nspace onwards). kappa never falls
			   below the last one of this wavelength: from a J close
			   to the solution the iteration takes a few steps only,
			   and the span of so few differences misses the slowest
			   mode -- --                                -------------- */

		if (r <= input.bg_iterLimit) {
			kn    = fmax(secantKappa(Nspace, Npairs, inew, hu, hf, w,
									 hf + KAPPA_PAIRS * Nspace),
						 *kap);
			fresh = TRUE;
			if (errFactor(kn, amax) * r <= input.bg_iterLimit) {
				converged = TRUE;
				break;
			}
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
		r = relChange(Nspace, u, g);
		pushPair(Nspace, &Npairs, &inew, u, g, hu, hf);
		fresh = FALSE;
	}
	if (!fresh)
		kn = fmax(secantKappa(Nspace, Npairs, inew, hu, hf, w,
							  hf + KAPPA_PAIRS * Nspace),
				  *kap);
	if (kn > 0.0) *kap = kn;
	memcpy(J, g, Nspace * sizeof(double));

	pthread_mutex_lock(&lock);
	stats.Niterated++;
	stats.Nsweeps += niter;
	if (!converged) {
		stats.Nunconverged++;
		r *= errFactor(kn, amax);
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
	double *Jfs, *A, *albedo, kap = kappa[nspect];

	/* --- One formal solution from Jdag, the iteration from there to
		   J*, the source of the main sweep. A wavelength the scattering
		   no longer changes by more than the limit is not iterated
		   first in the next sweep (J* is then the Jacobi step g). A
		   wavelength is solved by one thread at a time, so its flag
		   and kappa need no lock --                     -------------- */

	Jfs    = (double *) malloc(5 * Nspace * sizeof(double));
	A      = Jfs + Nspace;
	albedo = A + 3 * Nspace;
	sweep(ctx, Jdag, Jfs, A, albedo);
	niter         = bgScatterIterate(Nspace, Jdag, Jfs, A, albedo, Jstar, &kap,
									 sweep, ctx);
	kappa[nspect] = (float) kap;
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
	if (errFactor(kappa[nspect], maxAlbedo(Nspace, albedo)) *
			relChange(Nspace, u, g) >
		input.bg_iterLimit)
		pending[nspect] = 1;
	free(g);
}
/* ------- end ---------------------------- bgScatterAfter.c -------- */

/* ------- begin -------------------------- bgScatterStats.c -------- */

void bgScatterStats(BgScatterStats *s, int reset)
{
	pthread_mutex_lock(&lock);
	*s         = stats;
	s->Nlinear = Nlinear;
	if (reset) memset(&stats, 0, sizeof(stats));
	pthread_mutex_unlock(&lock);
}
/* ------- end ---------------------------- bgScatterStats.c -------- */

/* ------- begin -------------------------- bgScatterLinear.c ------- */

bool_t bgScatterLinear(int nspect, double albedo)
{
	bool_t lin;

	/* --- A wavelength is solved by one thread at a time, so its flag
		   needs no lock; the counts do -- --            -------------- */

	if (linear == NULL) return FALSE;
	if (albedo < 0.0) return linear[nspect];
	lin = (albedo > input.bg_linearAlbedo);
	if (lin != linear[nspect]) {
		pthread_mutex_lock(&lock);
		if (lin) {
			Nlinear++;
			stats.Nto++;
		} else {
			Nlinear--;
			stats.Nback++;
		}
		pthread_mutex_unlock(&lock);
		linear[nspect] = lin;
	}
	return lin;
}
/* ------- end ---------------------------- bgScatterLinear.c ------- */

/* ------- end ---------------------------- bgscatter.c ------------- */
