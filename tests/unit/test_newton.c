/* ------- file: -------------------------- test_newton.c -----------

	   Matrix-free GMRES and the forcing term (newton.c, S-4 in
	   docs/SOLVER.md).

	   -- with k = n and tolerance 0, GMRES solves a nonsymmetric system
		  to rounding (against Gaussian elimination), whatever the
		  weights;
	   -- after k < n products, x_k minimizes the weighted residual over
		  the Krylov space: <A v_l, b - A x_k>_W = 0 for every basis
		  vector, and the norm from the rotations is the true one;
	   -- the residual history never increases, and a tolerance stops
		  the iteration at the first product that meets it;
	   -- with a mask, the result is that of the system restricted to
		  the masked entries, and zero elsewhere;
	   -- A = I: one product, an invariant space, x = b to rounding;
	   -- two runs are bit-identical;
	   -- the forcing term: the first step, Eisenstat-Walker choice 2,
		  their safeguard, the cap and the floor.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "newton.h"
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

static double A[N][N], b[N], w[N];
static int    Nproducts;

static void product(void *ctx, const double *v, double *y)
{
	(void) ctx;
	Nproducts++;
	for (int i = 0; i < N; i++) {
		double s = 0.0;
		for (int j = 0; j < N; j++) s += A[i][j] * v[j];
		y[i] = s;
	}
}

static void identity(void *ctx, const double *v, double *y)
{
	(void) ctx;
	Nproducts++;
	memcpy(y, v, N * sizeof(double));
}

/* --- Dense solve of the rows and columns in mask, by Gaussian
	   elimination with partial pivoting; zero elsewhere -- -------- */

static void denseSolve(const unsigned char *mask, double *x)
{
	int    idx[N], n = 0;
	double a[N][N + 1];

	for (int i = 0; i < N; i++)
		if (mask[i]) idx[n++] = i;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) a[i][j] = A[idx[i]][idx[j]];
		a[i][n] = b[idx[i]];
	}
	for (int c = 0; c < n; c++) {
		int p = c;
		for (int i = c + 1; i < n; i++)
			if (fabs(a[i][c]) > fabs(a[p][c])) p = i;
		for (int j = 0; j <= n; j++) {
			double t = a[c][j];
			a[c][j]  = a[p][j];
			a[p][j]  = t;
		}
		for (int i = c + 1; i < n; i++) {
			double f = a[i][c] / a[c][c];
			for (int j = c; j <= n; j++) a[i][j] -= f * a[c][j];
		}
	}
	memset(x, 0, N * sizeof(double));
	for (int i = n - 1; i >= 0; i--) {
		double s = a[i][n];
		for (int j = i + 1; j < n; j++) s -= a[i][j] * x[idx[j]];
		x[idx[i]] = s / a[i][i];
	}
}

static double wnorm(const double *v, const unsigned char *mask)
{
	return weightedNorm(N, v, mask, w);
}

static void residual(const double *x, double *r)
{
	product(NULL, x, r);
	for (int i = 0; i < N; i++) r[i] = b[i] - r[i];
}

static void setup(void)
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) A[i][j] = (urand() - 0.5) / sqrt((double) N);
		A[i][i] += 1.0 + 0.5 * urand();
		b[i] = urand() - 0.5;
		w[i] = 0.2 + urand();
	}
}

static void testFull(void)
{
	unsigned char mask[N];
	double        x[N], xd[N];
	Gmres        *gm = gmresInit(N, N);

	memset(mask, 1, N);
	GmresResult res = gmresSolve(gm, product, NULL, b, mask, w, N, 0.0, x);
	denseSolve(mask, xd);
	double err = 0.0, scale = 0.0;
	for (int i = 0; i < N; i++) {
		err   = fmax(err, fabs(x[i] - xd[i]));
		scale = fmax(scale, fabs(xd[i]));
	}
	CHECK(err <= 1.0E-12 * scale, "k = n: %.2e of %.2e", err, scale);
	CHECK(res.rnorm <= 1.0E-12 * res.rnorm0, "k = n: residual %.2e of %.2e", res.rnorm,
		  res.rnorm0);
	gmresFree(gm);
}

static void testPartial(void)
{
	unsigned char mask[N];
	double        x[N], y[N], r[N], av[N];
	Gmres        *gm = gmresInit(N, N);

	memset(mask, 1, N);
	for (int k = 1; k <= 12; k++) {
		GmresResult res = gmresSolve(gm, product, NULL, b, mask, w, k, 0.0, x);
		CHECK(res.Niter == k && !res.breakdown, "k = %d: %d products", k, res.Niter);
		residual(x, r);
		CHECK_CLOSE(wnorm(r, mask), res.rnorm, 1.0E-10, 1.0E-14 * res.rnorm0,
					"k = %d: rotation norm vs true", k);
		for (int l = 0; l < k; l++) {
			product(NULL, gm->V[l], av);
			double d = 0.0;
			for (int i = 0; i < N; i++) d += w[i] * w[i] * av[i] * r[i];
			CHECK(fabs(d) <= 1.0E-12 * wnorm(av, mask) * res.rnorm0,
				  "k = %d: <A v_%d, r>_W = %.2e", k, l, d);
		}
		for (int j = 0; j < k; j++)
			CHECK(gm->resid[j + 1] <= gm->resid[j] * (1.0 + 1.0E-14),
				  "k = %d: residual grew at %d", k, j);
	}

	/* --- A tolerance stops at the first product that meets it -- */

	GmresResult full = gmresSolve(gm, product, NULL, b, mask, w, N, 0.0, x);
	double      tol  = 1.0E-4 * full.rnorm0;
	int         want = 0;
	while (gm->resid[want] > tol) want++;
	GmresResult res = gmresSolve(gm, product, NULL, b, mask, w, N, tol, x);
	CHECK(res.Niter == want && res.rnorm <= tol && res.Niter < N,
		  "tolerance: %d products (want %d), %.2e vs %.2e", res.Niter, want, res.rnorm,
		  tol);

	/* --- Two runs are bit-identical --                 ---------- */

	gmresSolve(gm, product, NULL, b, mask, w, 7, 0.0, x);
	gmresSolve(gm, product, NULL, b, mask, w, 7, 0.0, y);
	CHECK(!memcmp(x, y, sizeof(x)), "two runs differ");
	gmresFree(gm);
}

static void testMask(void)
{
	unsigned char mask[N];
	double        x[N], xd[N];
	Gmres        *gm = gmresInit(N, N);

	for (int i = 0; i < N; i++) mask[i] = (i % 4 != 1);
	gmresSolve(gm, product, NULL, b, mask, w, N, 0.0, x);
	denseSolve(mask, xd);
	double err = 0.0, scale = 0.0;
	int    outside = 0;
	for (int i = 0; i < N; i++) {
		err   = fmax(err, fabs(x[i] - xd[i]));
		scale = fmax(scale, fabs(xd[i]));
		if (!mask[i] && x[i] != 0.0) outside++;
	}
	CHECK(err <= 1.0E-12 * scale, "mask: %.2e of %.2e", err, scale);
	CHECK(outside == 0, "mask: %d entries outside are not zero", outside);
	gmresFree(gm);
}

static void testIdentity(void)
{
	unsigned char mask[N];
	double        x[N];
	Gmres        *gm = gmresInit(N, 10);

	memset(mask, 1, N);
	Nproducts       = 0;
	GmresResult res = gmresSolve(gm, identity, NULL, b, mask, w, 10, 0.0, x);
	CHECK(res.Niter == 1 && Nproducts == 1 && res.breakdown, "A = I: %d products",
		  res.Niter);
	for (int i = 0; i < N; i++) CHECK_CLOSE(x[i], b[i], 1.0E-15, 0.0, "A = I: x[%d]", i);
	gmresFree(gm);
}

static void testForcing(void)
{
	CHECK(newtonForcing(1, 0.0, 1.0, 0.0, 0.9, 0.0) == 0.9, "%s", "first step");
	CHECK_CLOSE(newtonForcing(0, 0.01, 0.5, 1.0, 0.9, 0.0), 0.9 * 0.25, 1.0E-15, 0.0,
				"%s", "choice 2");
	CHECK_CLOSE(newtonForcing(0, 0.9, 0.1, 1.0, 0.9, 0.0), 0.9 * 0.81, 1.0E-15, 0.0,
				"%s", "safeguard");
	CHECK_CLOSE(newtonForcing(0, 0.3, 0.1, 1.0, 0.9, 0.0), 0.9 * 0.01, 1.0E-15, 0.0,
				"%s", "safeguard below 0.1 not applied");
	CHECK(newtonForcing(0, 0.01, 2.0, 1.0, 0.9, 0.0) == 0.9, "%s", "cap");
	CHECK_CLOSE(newtonForcing(0, 0.01, 1.0E-6, 1.0E-3, 0.9, 1.0E-6), 0.5, 1.0E-15, 0.0,
				"%s", "floor");
	CHECK(newtonForcing(0, 0.01, 1.0E-7, 1.0E-3, 0.9, 1.0E-6) == 0.9, "%s", "floor capped");
}

int main(void)
{
	setup();
	testFull();
	testPartial();
	testMask();
	testIdentity();
	testForcing();
	return check_summary("test_newton");
}

/* ------- end ---------------------------- test_newton.c ----------- */
