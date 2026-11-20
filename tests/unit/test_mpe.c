/* ------- file: -------------------------- test_mpe.c --------------

	   Minimal polynomial extrapolation (mpe.c, the error estimate of
	   ERROR_LIMIT, N-3).

	   -- A linear iteration x <- M x + b with M symmetric, 5 distinct
		  eigenvalues (0.98, 0.5, -0.6, 0.1 and 1.85, so the iteration
		  diverges): from 6 points (K = 6, 5 columns) the limit is
		  exact, to 1e-9 of |s - V_0|; from 5 (K = 5) it is not.
		  A single mode with lambda = 0.99 and K = 2 is exact.
	   -- All differences zero: TRUE with s - V_0 = 0, no NaN.
	   -- A mode with lambda = 1 (constant differences): no finite
		  limit, FALSE, D[K-1] unchanged.
	   -- Parallel columns (a single mode, K = 6): the dependent ones
		  are dropped and the limit is still exact.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "mpe.h"
#include "check.h"

#define N 30

/* --- Deterministic pseudo-random numbers in [0, 1) -- ------------- */

static unsigned long long rng_state = 88172645463325252ULL;

static double urand(void)
{
	rng_state ^= rng_state << 13;
	rng_state ^= rng_state >> 7;
	rng_state ^= rng_state << 17;
	return (rng_state >> 11) * (1.0 / 9007199254740992.0);
}

static double M[N][N], bvec[N], sfix[N];

/* --- M = H diag(lam) H with H a Householder reflection; the fixed
	   point s of x <- M x + b is set first, b = s - M s -- -------- */

static void makeMap(const double *lam)
{
	double v[N], vv = 0.0;

	for (int i = 0; i < N; i++) {
		v[i] = urand() - 0.5;
		vv += v[i] * v[i];
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			double s = 0.0;

			for (int k = 0; k < N; k++) {
				double hik = ((i == k) ? 1.0 : 0.0) - 2.0 * v[i] * v[k] / vv;
				double hkj = ((k == j) ? 1.0 : 0.0) - 2.0 * v[k] * v[j] / vv;

				s += hik * lam[k] * hkj;
			}
			M[i][j] = s;
		}
	}
	for (int i = 0; i < N; i++) sfix[i] = urand() - 0.5;
	for (int i = 0; i < N; i++) {
		bvec[i] = sfix[i];
		for (int j = 0; j < N; j++) bvec[i] -= M[i][j] * sfix[j];
	}
}

/* --- K + 1 points from x0, their K differences in D; returns
	   max |d - (s - x0)| / max |s - x0| with d from mpeExtrapolate,
	   or -1 when it returned FALSE -- ---------------------------- */

static double runMpe(int K, const double *x0, double **D)
{
	double x[N], y[N], err = 0.0, ref = 0.0;

	memcpy(x, x0, sizeof(x));
	for (int k = 0; k < K; k++) {
		for (int i = 0; i < N; i++) {
			y[i] = bvec[i];
			for (int j = 0; j < N; j++) y[i] += M[i][j] * x[j];
			D[k][i] = y[i] - x[i];
		}
		memcpy(x, y, sizeof(x));
	}
	if (!mpeExtrapolate(K, N, D)) return -1.0;
	for (int i = 0; i < N; i++) {
		err = fmax(err, fabs(D[K - 1][i] - (sfix[i] - x0[i])));
		ref = fmax(ref, fabs(sfix[i] - x0[i]));
	}
	return err / ref;
}

int main(void)
{
	double       lam[N], x0[N], *D[8], e;
	const double five[5] = {0.98, 0.5, -0.6, 0.1, 1.85};
	bool_t       ok, nan;

	for (int k = 0; k < 8; k++) D[k] = (double *) malloc(N * sizeof(double));

	/* --- Five distinct eigenvalues, one of them > 1 -- ------------ */

	for (int i = 0; i < N; i++) lam[i] = five[i % 5];
	makeMap(lam);
	for (int i = 0; i < N; i++) x0[i] = sfix[i] + urand() - 0.5;
	e = runMpe(6, x0, D);
	CHECK(e >= 0.0 && e < 1.0E-9,
		  "5 eigenvalues (one 1.85), K = 6: limit exact (error %.2e of |s - V_0|)", e);
	e = runMpe(5, x0, D);
	CHECK(e < 0.0 || e > 1.0E-6,
		  "5 eigenvalues, K = 5: not exact (error %.2e of |s - V_0|)", e);

	/* --- One mode, lambda = 0.99 -- ------------------------------- */

	for (int i = 0; i < N; i++) lam[i] = 0.99;
	makeMap(lam);
	e = runMpe(2, x0, D);
	CHECK(e >= 0.0 && e < 1.0E-10, "one mode 0.99, K = 2: exact (%.2e)", e);

	/* --- The same with K = 6: four dependent columns dropped -- --- */

	e = runMpe(6, x0, D);
	CHECK(e >= 0.0 && e < 1.0E-8, "one mode 0.99, K = 6 (dependent columns): exact (%.2e)",
		  e);

	/* --- Zero differences: converged, s - V_0 = 0 -- ------------- */

	for (int k = 0; k < 6; k++) memset(D[k], 0, N * sizeof(double));
	ok  = mpeExtrapolate(6, N, D);
	nan = FALSE;
	e   = 0.0;
	for (int i = 0; i < N; i++) {
		if (!isfinite(D[5][i])) nan = TRUE;
		e = fmax(e, fabs(D[5][i]));
	}
	CHECK(ok && !nan && e == 0.0, "zero differences: TRUE, s - V_0 = 0 (max %.1e)", e);

	/* --- lambda = 1: constant differences, no finite limit -- ---- */

	for (int k = 0; k < 4; k++)
		for (int i = 0; i < N; i++) D[k][i] = 1.0 + 0.1 * i;
	ok  = mpeExtrapolate(4, N, D);
	nan = FALSE;
	for (int i = 0; i < N; i++)
		if (D[3][i] != 1.0 + 0.1 * i) nan = TRUE;
	CHECK(!ok && !nan, "constant differences (lambda = 1): FALSE, D[K-1] unchanged");

	CHECK(!mpeExtrapolate(1, N, D), "K = 1: FALSE");

	for (int k = 0; k < 8; k++) free(D[k]);
	return check_summary("test_mpe");
}
