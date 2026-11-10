/* ------- file: -------------------------- newton.c ----------------

	   Matrix-free GMRES and the forcing term of inexact Newton, for
	   SOLVER = NEWTON_KRYLOV (S-4 in docs/SOLVER.md).

	   --------------------------                      ----------RH-- */

/* --- Generic: no RH globals. The caller supplies the product
	   y = A v (in S-4 a finite-difference Jacobian product, one G
	   evaluation each), the weights and the mask.

	   gmresSolve() solves A x = b approximately from x_0 = 0 (Saad &
	   Schultz 1986) in the inner product

		 <a, b>_W = sum over the mask of (w_i a_i) (w_i b_i),

	   the norm of S-3 (each block counts as its mean square). Entries
	   outside the mask are zero in every basis vector and in x. The
	   Arnoldi basis is orthogonalized by classical Gram-Schmidt with
	   reorthogonalization (CGS2), in a fixed serial order, so reruns are
	   bit-identical; the least-squares problem is reduced by Givens
	   rotations, which give ||b - A x_j||_W after every product without
	   forming x_j. It stops after kmax products, when that norm is
	   <= tol, or when the Krylov space becomes invariant (then x is the
	   solution in it). There is no restart: the step is the last
	   iterate, the best GMRES found, as the norm never increases.

	   newtonForcing() is the relative tolerance eta_k of the Newton
	   step k: Eisenstat & Walker (1996) choice 2 with their safeguard,
	   capped at eta_max, and kept above 0.5 floor_abs / ||R_k|| so that
	   the last steps do not solve the linear problem further than the
	   stopping tolerance floor_abs needs (Kelley 1995, ch. 6). --  ---- */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "newton.h"

/* ------- begin -------------------------- gmresInit.c ------------- */

Gmres *gmresInit(long N, int k)
{
	int    j;
	Gmres *gm = (Gmres *) calloc(1, sizeof(Gmres));

	gm->N = N;
	gm->k = (k > 0) ? k : 1;
	gm->V = (double **) malloc((gm->k + 1) * sizeof(double *));
	for (j = 0; j <= gm->k; j++) gm->V[j] = (double *) calloc(N, sizeof(double));
	gm->H     = (double *) calloc((size_t) (gm->k + 1) * gm->k, sizeof(double));
	gm->cs    = (double *) calloc(gm->k, sizeof(double));
	gm->sn    = (double *) calloc(gm->k, sizeof(double));
	gm->g     = (double *) calloc(gm->k + 1, sizeof(double));
	gm->y     = (double *) calloc(gm->k, sizeof(double));
	gm->h     = (double *) calloc(gm->k + 1, sizeof(double));
	gm->resid = (double *) calloc(gm->k + 1, sizeof(double));

	return gm;
}
/* ------- end ---------------------------- gmresInit.c ------------- */

/* ------- begin -------------------------- gmresFree.c ------------- */

void gmresFree(Gmres *gm)
{
	int j;

	if (gm == NULL) return;
	for (j = 0; j <= gm->k; j++) free(gm->V[j]);
	free(gm->V);
	free(gm->H);
	free(gm->cs);
	free(gm->sn);
	free(gm->g);
	free(gm->y);
	free(gm->h);
	free(gm->resid);
	free(gm);
}
/* ------- end ---------------------------- gmresFree.c ------------- */

/* ------- begin -------------------------- wdot.c ------------------ */

static double wdot(long N, const double *a, const double *b,
				   const unsigned char *mask, const double *w)
{
	register long i;

	double s = 0.0;

	for (i = 0; i < N; i++)
		if (mask[i]) s += (w[i] * a[i]) * (w[i] * b[i]);
	return s;
}
/* ------- end ---------------------------- wdot.c ------------------ */

/* ------- begin -------------------------- weightedNorm.c ---------- */

double weightedNorm(long N, const double *v, const unsigned char *mask,
					const double *w)
{
	return sqrt(wdot(N, v, v, mask, w));
}
/* ------- end ---------------------------- weightedNorm.c ---------- */

/* ------- begin -------------------------- gmresSolve.c ------------ */

#define H_(i, j) gm->H[(size_t) (j) * (gm->k + 1) + (i)]

GmresResult gmresSolve(Gmres *gm, GmresProduct A, void *ctx, const double *b,
					   const unsigned char *mask, const double *w, int kmax,
					   double tol, double *x)
{
	register long i;
	register int  j, l, pass;

	long        N = gm->N;
	int         k = (kmax < gm->k) ? kmax : gm->k, Ncol = 0;
	double      beta, anorm, hn, t, denom;
	GmresResult res = {0, 0.0, 0.0, 0};

	memset(x, 0, N * sizeof(double));
	beta       = weightedNorm(N, b, mask, w);
	res.rnorm0 = res.rnorm = beta;
	gm->resid[0]           = beta;
	if (!(beta > 0.0) || k <= 0) return res;

	for (i = 0; i < N; i++) gm->V[0][i] = (mask[i]) ? b[i] / beta : 0.0;
	gm->g[0] = beta;

	for (j = 0; j < k; j++) {
		A(ctx, gm->V[j], gm->V[j + 1]);
		res.Niter++;
		for (i = 0; i < N; i++)
			if (!mask[i]) gm->V[j + 1][i] = 0.0;
		anorm = weightedNorm(N, gm->V[j + 1], mask, w);

		/* --- CGS2: all projections of a pass from the same vector,
			   then subtracted; twice --               -------------- */

		for (l = 0; l <= j; l++) H_(l, j) = 0.0;
		for (pass = 0; pass < 2; pass++) {
			for (l = 0; l <= j; l++) gm->h[l] = wdot(N, gm->V[l], gm->V[j + 1], mask, w);
			for (l = 0; l <= j; l++) {
				H_(l, j) += gm->h[l];
				for (i = 0; i < N; i++) gm->V[j + 1][i] -= gm->h[l] * gm->V[l][i];
			}
		}
		hn           = weightedNorm(N, gm->V[j + 1], mask, w);
		H_(j + 1, j) = hn;

		/* --- Previous rotations on the new column, then a new one
			   that zeroes H(j + 1, j) --              -------------- */

		for (l = 0; l < j; l++) {
			t            = gm->cs[l] * H_(l, j) + gm->sn[l] * H_(l + 1, j);
			H_(l + 1, j) = -gm->sn[l] * H_(l, j) + gm->cs[l] * H_(l + 1, j);
			H_(l, j)     = t;
		}
		denom = hypot(H_(j, j), H_(j + 1, j));
		if (!(denom > 0.0)) { /* A v_j = 0 on the basis: singular */
			res.breakdown = 1;
			break;
		}
		gm->cs[j]    = H_(j, j) / denom;
		gm->sn[j]    = H_(j + 1, j) / denom;
		H_(j, j)     = denom;
		H_(j + 1, j) = 0.0;
		gm->g[j + 1] = -gm->sn[j] * gm->g[j];
		gm->g[j]     = gm->cs[j] * gm->g[j];

		Ncol             = j + 1;
		res.rnorm        = fabs(gm->g[j + 1]);
		gm->resid[j + 1] = res.rnorm;

		if (hn <= 1.0E-14 * anorm) { /* invariant: exact in the space */
			res.breakdown = 1;
			break;
		}
		for (i = 0; i < N; i++) gm->V[j + 1][i] /= hn;
		if (res.rnorm <= tol) break;
	}
	/* --- y = R^-1 g, x = V y --                       -------------- */

	for (l = Ncol - 1; l >= 0; l--) {
		t = gm->g[l];
		for (j = l + 1; j < Ncol; j++) t -= H_(l, j) * gm->y[j];
		gm->y[l] = t / H_(l, l);
	}
	for (l = 0; l < Ncol; l++)
		for (i = 0; i < N; i++) x[i] += gm->y[l] * gm->V[l][i];

	return res;
}
#undef H_
/* ------- end ---------------------------- gmresSolve.c ------------ */

/* ------- begin -------------------------- newtonForcing.c --------- */

double newtonForcing(int first, double eta_prev, double rnorm,
					 double rnorm_prev, double eta_max, double floor_abs)
{
	double eta, safe;

	if (first || !(rnorm_prev > 0.0)) {
		eta = eta_max;
	} else {
		eta  = NEWTON_EW_GAMMA * pow(rnorm / rnorm_prev, NEWTON_EW_ALPHA);
		safe = NEWTON_EW_GAMMA * pow(eta_prev, NEWTON_EW_ALPHA);
		if (safe > NEWTON_EW_SAFE && safe > eta) eta = safe;
	}
	if (eta > eta_max) eta = eta_max;
	if (rnorm > 0.0 && eta < 0.5 * floor_abs / rnorm) eta = 0.5 * floor_abs / rnorm;
	if (eta > eta_max) eta = eta_max;

	return eta;
}
/* ------- end ---------------------------- newtonForcing.c --------- */
