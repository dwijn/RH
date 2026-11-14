/* ------- file: -------------------------- anderson.c --------------

	   Anderson acceleration, type II with mixing 1 (Walker & Ni 2011),
	   on a history of iterate pairs (u_i, g_i = G(u_i)). S-3 in
	   docs/SOLVER.md.

	   --------------------------                      ----------RH-- */

/* --- Generic: no RH globals. The caller chooses the variables (log,
	   scaled), the weights and what to do with the result.

	   With the residuals f_i = g_i - u_i of the m_k + 1 newest pairs,
	   the columns are differences of consecutive pairs, newest first:

		 dF_j = f_(k-j) - f_(k-j-1),   dG_j = g_(k-j) - g_(k-j-1),

	   j = 0 .. m_k - 1. andersonSolve() finds

		 gamma = argmin || W (f_k - dF gamma) ||_2

	   over the entries [lo, hi) by a QR factorization of W dF,
	   classical Gram-Schmidt with reorthogonalization (CGS2), recomputed
	   from scratch every call, in a fixed serial order (bit-identical
	   reruns). The condition number of R is computed exactly as
	   ||R||_1 ||R^-1||_1. While it exceeds cond_max the oldest column is
	   left out; with the columns newest first that is the last one, and
	   the leading block of R and of R^-1 stays valid. andersonApply()
	   then gives

		 u_(k+1) = g_k - theta dG gamma

	   (theta = 1 is Anderson's step, smaller values damp it).

	   Entries can be invalid in a pair (log of a value <= 0, say). An
	   entry takes part only if it is valid in every pair held, and only
	   if its range used at least one column: andersonMask is 1 exactly
	   where andersonApply extrapolated. Elsewhere it returns g_k
	   unchanged, and the caller should copy G's value itself. --      -------------- */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "anderson.h"

/* ------- begin -------------------------- pairIndex.c ------------- */

static int pairIndex(const Anderson *a, int j)
{
	/* --- Ring index of the j-th newest pair (j = 0: newest) -- ---- */

	return (a->newest - j + 2 * (a->m + 1)) % (a->m + 1);
}
/* ------- end ---------------------------- pairIndex.c ------------- */

/* ------- begin -------------------------- andersonInit.c ---------- */

Anderson *andersonInit(long N, int m)
{
	int       j;
	Anderson *a = (Anderson *) calloc(1, sizeof(Anderson));

	a->N      = N;
	a->m      = (m > 0) ? m : 0;
	a->Npairs = 0;
	a->newest = 0;

	a->u     = (double **) malloc((a->m + 1) * sizeof(double *));
	a->g     = (double **) malloc((a->m + 1) * sizeof(double *));
	a->valid = (unsigned char **) malloc((a->m + 1) * sizeof(unsigned char *));
	for (j = 0; j <= a->m; j++) {
		a->u[j]     = (double *) malloc(N * sizeof(double));
		a->g[j]     = (double *) malloc(N * sizeof(double));
		a->valid[j] = (unsigned char *) malloc(N);
	}
	a->mask = (unsigned char *) calloc(N, 1);
	if (a->m > 0) {
		a->A = (double **) malloc(a->m * sizeof(double *));
		for (j = 0; j < a->m; j++)
			a->A[j] = (double *) malloc(N * sizeof(double));
		a->R    = (double *) calloc(a->m * a->m, sizeof(double));
		a->Rinv = (double *) calloc(a->m * a->m, sizeof(double));
		a->r    = (double *) calloc(a->m, sizeof(double));
	}
	a->b = (double *) malloc(N * sizeof(double));

	return a;
}
/* ------- end ---------------------------- andersonInit.c ---------- */

/* ------- begin -------------------------- andersonFree.c ---------- */

void andersonFree(Anderson *a)
{
	int j;

	if (a == NULL) return;
	for (j = 0; j <= a->m; j++) {
		free(a->u[j]);
		free(a->g[j]);
		free(a->valid[j]);
	}
	free(a->u);
	free(a->g);
	free(a->valid);
	free(a->mask);
	for (j = 0; j < a->m; j++) free(a->A[j]);
	free(a->A);
	free(a->R);
	free(a->Rinv);
	free(a->r);
	free(a->b);
	free(a);
}
/* ------- end ---------------------------- andersonFree.c ---------- */

/* ------- begin -------------------------- andersonReset.c --------- */

void andersonReset(Anderson *a, int keep_newest)
{
	/* --- Forget the history, or all of it but the newest pair -- -- */

	a->Npairs = (keep_newest && a->Npairs > 0) ? 1 : 0;
}
/* ------- end ---------------------------- andersonReset.c --------- */

/* ------- begin -------------------------- andersonPush.c ---------- */

void andersonPush(Anderson *a, const double *u, const double *g,
				  const unsigned char *valid)
{
	/* --- Add the pair (u, g = G(u)); valid == NULL: all entries
		   valid. The oldest pair drops out beyond m + 1 -- ------- */

	a->newest = (a->Npairs == 0) ? 0 : (a->newest + 1) % (a->m + 1);
	memcpy(a->u[a->newest], u, a->N * sizeof(double));
	memcpy(a->g[a->newest], g, a->N * sizeof(double));
	if (valid)
		memcpy(a->valid[a->newest], valid, a->N);
	else
		memset(a->valid[a->newest], 1, a->N);
	if (a->Npairs < a->m + 1) a->Npairs++;
}
/* ------- end ---------------------------- andersonPush.c ---------- */

/* ------- begin -------------------------- andersonNcols.c --------- */

int andersonNcols(const Anderson *a)
{
	return (a->Npairs > 1) ? a->Npairs - 1 : 0;
}
/* ------- end ---------------------------- andersonNcols.c --------- */

/* ------- begin -------------------------- andersonMask.c ---------- */

const unsigned char *andersonMask(const Anderson *a)
{
	return a->mask;
}
/* ------- end ---------------------------- andersonMask.c ---------- */

/* ------- begin -------------------------- dot.c ------------------- */

static double dot(const double *x, const double *y, long n)
{
	register long i;

	double s = 0.0;

	for (i = 0; i < n; i++) s += x[i] * y[i];
	return s;
}
/* ------- end ---------------------------- dot.c ------------------- */

/* ------- begin -------------------------- andersonSolve.c --------- */

AndersonSolve andersonSolve(Anderson *a, long lo, long hi, const double *w,
							double cond_max, double *gamma)
{
	register long i;
	register int  j, l, pass;

	int           m = a->m, Nmax = andersonNcols(a), Ncols, p, q;
	long          len = hi - lo;
	double       *r   = a->r, *v, norm, Rnorm, Rinvnorm, s, cond;
	AndersonSolve res = {0, 0, 0.0};

	if (len <= 0) return res;
	if (Nmax == 0) {
		memset(a->mask + lo, 0, len);
		return res;
	}

	/* --- Entries valid in every pair held --         -------------- */

	for (i = lo; i < hi; i++) {
		a->mask[i] = 1;
		for (j = 0; j < a->Npairs; j++) {
			if (!a->valid[pairIndex(a, j)][i]) {
				a->mask[i] = 0;
				break;
			}
		}
	}
	/* --- Weighted residual and difference columns --  -------------- */

	p = pairIndex(a, 0);
	for (i = lo; i < hi; i++)
		a->b[i - lo] = (a->mask[i]) ? w[i] * (a->g[p][i] - a->u[p][i]) : 0.0;
	for (j = 0; j < Nmax; j++) {
		p = pairIndex(a, j);
		q = pairIndex(a, j + 1);
		for (i = lo; i < hi; i++) {
			a->A[j][i - lo] = (a->mask[i])
								  ? w[i] * ((a->g[p][i] - a->u[p][i]) -
											(a->g[q][i] - a->u[q][i]))
								  : 0.0;
		}
	}
	/* --- QR by CGS2; a column that is zero after orthogonalization
		   (or not finite) ends the usable set --       -------------- */

	memset(a->R, 0, m * m * sizeof(double));
	Ncols = Nmax;
	for (j = 0; j < Nmax; j++) {
		v = a->A[j];
		for (pass = 0; pass < 2; pass++) {
			for (l = 0; l < j; l++) r[l] = dot(a->A[l], v, len);
			for (l = 0; l < j; l++) {
				a->R[l * m + j] += r[l];
				for (i = 0; i < len; i++) v[i] -= r[l] * a->A[l][i];
			}
		}
		norm = sqrt(dot(v, v, len));
		if (!(norm > 0.0) || !isfinite(norm)) {
			Ncols = j;
			break;
		}
		a->R[j * m + j] = norm;
		for (i = 0; i < len; i++) v[i] /= norm;
	}
	/* --- R^-1 by back substitution, column by column. Its leading
		   blocks are the inverses of the leading blocks of R -- ---- */

	memset(a->Rinv, 0, m * m * sizeof(double));
	for (j = 0; j < Ncols; j++) {
		a->Rinv[j * m + j] = 1.0 / a->R[j * m + j];
		for (l = j - 1; l >= 0; l--) {
			s = 0.0;
			for (p = l + 1; p <= j; p++) s += a->R[l * m + p] * a->Rinv[p * m + j];
			a->Rinv[l * m + j] = -s / a->R[l * m + l];
		}
	}
	/* --- Leave out the oldest columns while kappa_1(R) > cond_max -- */

	cond = 0.0;
	while (Ncols > 0) {
		Rnorm = Rinvnorm = 0.0;
		for (j = 0; j < Ncols; j++) {
			s = 0.0;
			for (l = 0; l <= j; l++) s += fabs(a->R[l * m + j]);
			Rnorm = fmax(Rnorm, s);
			s     = 0.0;
			for (l = 0; l <= j; l++) s += fabs(a->Rinv[l * m + j]);
			Rinvnorm = fmax(Rinvnorm, s);
		}
		cond = Rnorm * Rinvnorm;
		if (isfinite(cond) && cond <= cond_max) break;
		Ncols--;
	}
	res.Ncols = Ncols;
	res.Ndrop = Nmax - Ncols;
	res.cond  = (Ncols > 0) ? cond : 0.0;
	if (Ncols == 0) { /* nothing extrapolated in [lo, hi) */
		memset(a->mask + lo, 0, len);
		return res;
	}

	/* --- gamma = R^-1 Q^T b --                        -------------- */

	for (l = 0; l < Ncols; l++) r[l] = dot(a->A[l], a->b, len);
	for (l = 0; l < Ncols; l++) {
		s = 0.0;
		for (j = l; j < Ncols; j++) s += a->Rinv[l * m + j] * r[j];
		gamma[l] = s;
	}
	return res;
}
/* ------- end ---------------------------- andersonSolve.c --------- */

/* ------- begin -------------------------- andersonApply.c --------- */

void andersonApply(const Anderson *a, long lo, long hi, const double *gamma,
				   int Ncols, double theta, double *u_out)
{
	register long i;
	register int  j;

	int    p0 = pairIndex(a, 0), p, q;
	double s;

	/* --- u_out = g_k - theta dG gamma on [lo, hi), g_k where an
		   entry is masked. Call after andersonSolve on the same
		   range --                                      -------------- */

	for (i = lo; i < hi; i++) {
		if (Ncols == 0 || !a->mask[i]) {
			u_out[i] = a->g[p0][i];
			continue;
		}
		s = 0.0;
		for (j = 0; j < Ncols; j++) {
			p = pairIndex(a, j);
			q = pairIndex(a, j + 1);
			s += gamma[j] * (a->g[p][i] - a->g[q][i]);
		}
		u_out[i] = a->g[p0][i] - theta * s;
	}
}
/* ------- end ---------------------------- andersonApply.c --------- */

/* ------- begin -------------------------- andersonNorm.c ---------- */

double andersonNorm(long N, const double *u, const double *g,
					const unsigned char *valid, const double *w)
{
	register long i;

	double s = 0.0, d;

	/* --- ||W (g - u)||_2 over the valid entries --   -------------- */

	for (i = 0; i < N; i++) {
		if (valid && !valid[i]) continue;
		d = w[i] * (g[i] - u[i]);
		s += d * d;
	}
	return sqrt(s);
}
/* ------- end ---------------------------- andersonNorm.c ---------- */

/* ------- begin -------------------------- solveAll.c -------------- */

static int solveAll(Anderson *a, const AndersonOptions *opt, const double *w,
					double *gamma, AndersonSolve *info, double theta,
					double *u_next)
{
	int b, Nsolve = (opt->per_block) ? opt->Nblock : 1, used = 0;

	/* --- One least squares over the whole vector, or one per block;
		   gamma[b * m ...] and info[b] per solve. Returns the largest
		   number of columns used --                     -------------- */

	for (b = 0; b < Nsolve; b++) {
		long lo = (opt->per_block) ? opt->start[b] : 0,
			 hi = (opt->per_block) ? opt->start[b + 1] : a->N;

		info[b] = andersonSolve(a, lo, hi, w, opt->cond_max, gamma + b * a->m);
		andersonApply(a, lo, hi, gamma + b * a->m, info[b].Ncols, theta, u_next);
		a->Ndropped += info[b].Ndrop;
		if (info[b].Ncols > used) used = info[b].Ncols;
	}
	return used;
}
/* ------- end ---------------------------- solveAll.c -------------- */

/* ------- begin -------------------------- andersonStep.c ---------- */

enum AndersonAction andersonStep(Anderson *a, const AndersonOptions *opt,
								 const double *u, const double *g,
								 const unsigned char *valid, const double *w,
								 int extrapolate, double *u_next,
								 double *gamma, AndersonSolve *info)
{
	double fnorm = andersonNorm(a->N, u, g, valid, w);

	/* --- One iteration, given the newest evaluation g = G(u):

		   1. If u was an extrapolated point and ||W f|| grew by more
			  than the factor c over the newest accepted pair: reject
			  it. The history restarts from that accepted pair, and
			  u_next is its image (the plain step, already computed;
			  the caller should copy that in its own variables).
		   2. Otherwise the pair joins the history. With extrapolate
			  and at least one column, u_next = g - dG gamma;
			  otherwise u_next = g (the plain step).

		   gamma needs m * Nblock and info Nblock entries (one with
		   GLOBAL coefficients) --                       -------------- */

	if (a->trial && opt->safeguard > 0.0 &&
		!(fnorm <= opt->safeguard * a->fnorm_acc)) {
		andersonRejectTrial(a);
		memcpy(u_next, a->g[a->newest], a->N * sizeof(double));
		return ANDERSON_REJECTED;
	}
	andersonPush(a, u, g, valid);
	a->fnorm_acc = fnorm;

	if (extrapolate && andersonNcols(a) > 0 &&
		solveAll(a, opt, w, gamma, info, 1.0, u_next) > 0) {
		a->trial = 1;
		a->Nextrapolated++;
		return ANDERSON_EXTRAPOLATED;
	}
	a->trial = 0;
	a->Nplain++;
	memcpy(u_next, g, a->N * sizeof(double));
	return ANDERSON_PLAIN;
}
/* ------- end ---------------------------- andersonStep.c ---------- */

/* ------- begin -------------------------- andersonRejectTrial.c --- */

void andersonRejectTrial(Anderson *a)
{
	/* --- Reject the extrapolated point last returned: the history
		   restarts from the newest accepted pair, whose image is the
		   next point (the caller has it; andersonStep copies it into
		   u_next). Also for a point whose G could not be evaluated
		   (negative opacity, N-40) --                   -------------- */

	a->trial = 0;
	a->Nrejected++;
	andersonReset(a, 1);
}
/* ------- end ---------------------------- andersonRejectTrial.c --- */

/* ------- begin -------------------------- andersonDamp.c ---------- */

void andersonDamp(const Anderson *a, const AndersonOptions *opt,
				  const double *gamma, const AndersonSolve *info,
				  double theta, double *u_next)
{
	int b, Nsolve = (opt->per_block) ? opt->Nblock : 1;

	/* --- The last extrapolation again, as g - theta dG gamma -- --- */

	for (b = 0; b < Nsolve; b++) {
		long lo = (opt->per_block) ? opt->start[b] : 0,
			 hi = (opt->per_block) ? opt->start[b + 1] : a->N;

		andersonApply(a, lo, hi, gamma + b * a->m, info[b].Ncols, theta, u_next);
	}
}
/* ------- end ---------------------------- andersonDamp.c ---------- */

/* ------- begin -------------------------- andersonAbandon.c ------- */

void andersonAbandon(Anderson *a, double *u_next)
{
	/* --- Take the plain step instead of the last extrapolation and
		   restart the history from the newest pair -- -------------- */

	a->trial = 0;
	a->Nextrapolated--;
	a->Nplain++;
	andersonReset(a, 1);
	memcpy(u_next, a->g[a->newest], a->N * sizeof(double));
}
/* ------- end ---------------------------- andersonAbandon.c ------- */

/* ------- end ---------------------------- anderson.c -------------- */
