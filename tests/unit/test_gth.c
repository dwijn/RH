/* ------- file: -------------------------- test_gth.c --------------

	   Rate equations by GTH elimination (SolveGTH in gth.c, and
	   statEquil in statequil.c, N-31).

	   Rate matrices built by detailed balance have a known solution:
	   with populations p_i and symmetric rates c_ij = c_ji >= 0, the
	   rates R_ij = c_ij / p_i (from level i to level j) satisfy
	   p_i R_ij = p_j R_ji, so p is the stationary state. Rounding R_ij
	   moves the exact solution by a few times eps per population (GTH
	   is componentwise accurate, and so is the problem: a relative
	   change d in each rate changes each population by at most about
	   2 N d, O'Cinneide 1993). The test checks:

	   -- random rates and populations, N = 1 .. 40: every population to
		  a few eps, relative;
	   -- groups of three levels bound by rates 1e12 times larger than
		  the rates between the groups (the fine-structure multiplets of
		  He_large, where LU loses up to 4e-3): GTH still to a few eps,
		  while LU as statEquil used it (SolveLinearEq) is printed with
		  --table for comparison;
	   -- populations spanning 1e-250 .. 1 (several ionization stages),
		  with the reference level at the most populated level and at
		  the least populated one (p rescaled in the back substitution);
	   -- a negative rate, a non-finite rate and a level without rates
		  to the others: the right status and level indices, n
		  unchanged;
	   -- statEquil on an atom with a negative collisional rate at one
		  of two locations: GTH at the other, LU there, both solving the
		  rate equations, and the work counter.

	   Run with --table to print the measured errors.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "statistics.h"
#include "check.h"

extern Atmosphere atmos;

static bool_t table;

/* --- Deterministic pseudo-random numbers in [0, 1) -- ------------- */

static unsigned long long rng_state = 88172645463325252ULL;

static double uniform(void)
{
	rng_state ^= rng_state << 13;
	rng_state ^= rng_state >> 7;
	rng_state ^= rng_state << 17;
	return (double) (rng_state >> 11) / 9007199254740992.0;
}

/* --- R_ij = c_ij / p_i, with c symmetric -- ----------------------- */

static void detailed_balance(int N, double **c, const double *p, double **R)
{
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) R[i][j] = (i == j) ? 0.0 : c[i][j] / p[i];
}

/* --- Max relative error of n against p normalized to ntotal -- ---- */

static double rel_error(int N, const double *n, const double *p, double ntotal)
{
	double sum = 0.0, err = 0.0;

	for (int i = 0; i < N; i++) sum += p[i];
	for (int i = 0; i < N; i++) {
		double want = p[i] * (ntotal / sum);
		err         = MAX(err, fabs(n[i] - want) / want);
	}
	return err;
}

/* --- LU as statEquil used it before N-31 -- ----------------------- */

static void solve_lu(int N, double **R, int root, double ntotal, double *n)
{
	double **A = matrix_double(N, N);

	for (int i = 0; i < N; i++) {
		double diag = 0.0;

		for (int j = 0; j < N; j++) {
			A[i][j] = (i == j) ? 0.0 : R[j][i]; /* rate j -> i */
			diag += (i == j) ? 0.0 : R[i][j];
		}
		A[i][i] = -diag;
		n[i]    = 0.0;
	}
	n[root] = ntotal;
	for (int j = 0; j < N; j++) A[root][j] = 1.0;
	SolveLinearEq(N, A, n, TRUE);
	freeMatrix((void **) A);
}

static int argmax(int N, const double *p)
{
	int m = 0;

	for (int i = 1; i < N; i++)
		if (p[i] > p[m]) m = i;
	return m;
}

/* --- Random rates and populations -- ------------------------------ */

/* --- Measured with gcc on x86-64: max relative error 7e-16 (N = 20, 40);
	   bound 4e-15 --                                  -------------- */

static void test_random(void)
{
	static const int sizes[] = {1, 2, 3, 5, 10, 20, 40};

	for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
		int      N = sizes[s], ib, jb;
		double **c = matrix_double(N, N), **R = matrix_double(N, N);
		double  *p = malloc(N * sizeof(double)), *n = malloc(N * sizeof(double));

		for (int i = 0; i < N; i++) {
			p[i] = pow(10.0, -6.0 * uniform());
			for (int j = 0; j < i; j++)
				c[i][j] = c[j][i] = (uniform() < 0.3) ? 0.0 : pow(10.0, 4.0 * uniform());
			if (i > 0) c[i][i - 1] = c[i - 1][i] = 1.0 + uniform(); /* connected */
		}
		detailed_balance(N, c, p, R);
		enum GTHStatus st = SolveGTH(N, R, argmax(N, p), 2.0, n, &ib, &jb);
		double         e  = rel_error(N, n, p, 2.0);

		CHECK(st == GTH_OK, "random N = %d: status %d", N, st);
		CHECK(e <= 4.0e-15, "random N = %d: relative error %.2e", N, e);
		if (table) printf("  random       N = %2d   GTH %.2e\n", N, e);

		freeMatrix((void **) c);
		freeMatrix((void **) R);
		free(p);
		free(n);
	}
}

/* --- Tightly coupled groups (He_large fine structure) -- ---------- */

/* --- Measured: GTH 5e-16, LU 6e-5;
	   bound for GTH 4e-15 --                          -------------- */

static void test_groups(void)
{
	const int N = 30; /* ten groups of three */
	int       ib, jb, root;
	double  **c = matrix_double(N, N), **R = matrix_double(N, N);
	double   *p = malloc(N * sizeof(double)), *n = malloc(N * sizeof(double));
	double    e_gth, e_lu;

	for (int i = 0; i < N; i++) {
		p[i] = pow(10.0, -8.0 * uniform());
		for (int j = 0; j < i; j++) {
			if (i / 3 == j / 3)
				c[i][j] = 1.0e12 * (1.0 + uniform());
			else
				c[i][j] = (uniform() < 0.5) ? 0.0 : 1.0 + uniform();
			c[j][i] = c[i][j];
		}
	}
	for (int i = 3; i < N; i += 3) c[i][i - 3] = c[i - 3][i] = 1.0; /* connected */
	detailed_balance(N, c, p, R);
	root = argmax(N, p);

	solve_lu(N, R, root, 1.0, n);
	e_lu = rel_error(N, n, p, 1.0);

	enum GTHStatus st = SolveGTH(N, R, root, 1.0, n, &ib, &jb);
	e_gth             = rel_error(N, n, p, 1.0);

	CHECK(st == GTH_OK, "groups: status %d", st);
	CHECK(e_gth <= 4.0e-15, "groups: relative error %.2e", e_gth);
	if (table) printf("  groups       N = %2d   GTH %.2e   LU %.2e\n", N, e_gth, e_lu);

	freeMatrix((void **) c);
	freeMatrix((void **) R);
	free(p);
	free(n);
}

/* --- Populations over 250 orders of magnitude -- ------------------ */

/* --- Measured: 8e-16 (reference level 0) and 1e-15 (level 29,
	   rescaled); bound 2e-14 -- */

static void test_range(void)
{
	const int N = 30;
	int       ib, jb;
	double  **c = matrix_double(N, N), **R = matrix_double(N, N);
	double   *p = malloc(N * sizeof(double)), *n = malloc(N * sizeof(double));

	for (int root = 0; root < N; root += N - 1) {
		for (int i = 0; i < N; i++) {
			p[i] = pow(10.0, -250.0 * i / (N - 1));
			for (int j = 0; j < N; j++) c[i][j] = 0.0;
		}
		for (int i = 1; i < N; i++) c[i][i - 1] = c[i - 1][i] = 1.0 + uniform();
		detailed_balance(N, c, p, R);

		enum GTHStatus st = SolveGTH(N, R, root, 1.0, n, &ib, &jb);
		double         e  = rel_error(N, n, p, 1.0);

		CHECK(st == GTH_OK, "range, root %d: status %d", root, st);
		CHECK(e <= 2.0e-14, "range, root %d: relative error %.2e", root, e);
		if (table) printf("  range 1e-250 root = %2d   GTH %.2e\n", root, e);
	}
	freeMatrix((void **) c);
	freeMatrix((void **) R);
	free(p);
	free(n);
}

/* --- Rates GTH cannot use -- -------------------------------------- */

static void test_status(void)
{
	const int N = 6;
	int       ib, jb;
	double  **R = matrix_double(N, N), n[6], n0[6];

	for (int i = 0; i < N; i++) n0[i] = n[i] = i + 1.0;

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) R[i][j] = (i == j) ? 0.0 : 1.0;
	R[2][5]           = -1.0e-3;
	enum GTHStatus st = SolveGTH(N, R, 0, 1.0, n, &ib, &jb);
	CHECK(st == GTH_NEGATIVE && ib == 2 && jb == 5,
		  "negative rate: status %d at %d -> %d", st, ib, jb);
	CHECK(!memcmp(n, n0, sizeof(n)), "negative rate: n changed");

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) R[i][j] = (i == j) ? 0.0 : 1.0;
	R[4][1] = NAN;
	st      = SolveGTH(N, R, 0, 1.0, n, &ib, &jb);
	CHECK(st == GTH_NONFINITE && ib == 4 && jb == 1,
		  "NaN rate: status %d at %d -> %d", st, ib, jb);
	CHECK(!memcmp(n, n0, sizeof(n)), "NaN rate: n changed");

	/* --- Level 5 has no rates out: the chain is reducible -- ------ */

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) R[i][j] = (i == j || i == 5) ? 0.0 : 1.0;
	st = SolveGTH(N, R, 0, 1.0, n, &ib, &jb);
	CHECK(st == GTH_REDUCIBLE && ib == 5, "reducible: status %d at %d", st, ib);
	CHECK(!memcmp(n, n0, sizeof(n)), "reducible: n changed");

	st = SolveGTH(1, R, 0, 3.0, n, &ib, &jb);
	CHECK(st == GTH_OK && n[0] == 3.0, "N = 1: status %d, n = %g", st, n[0]);

	freeMatrix((void **) R);
}

/* --- statEquil with a negative collisional rate -- ---------------- */

/* --- Measured: residual 2e-16 of the largest rate term, GTH and
	   LU; bound 1e-13 --                      -------------- */

static void test_statequil(void)
{
	const int N = 5, Nspace = 2;
	Atom      atom;
	double  **c = matrix_double(N, N), **R = matrix_double(N, N);
	double    p[5] = {1.0, 0.3, 0.05, 2e-3, 1e-6};

	memset(&atom, 0, sizeof(atom));
	strcpy(atom.ID, "XX");
	atom.Nlevel  = N;
	atom.Gamma   = matrix_double(N * N, Nspace);
	atom.C       = matrix_double(N * N, Nspace);
	atom.n       = matrix_double(N, Nspace);
	atom.ntotal  = malloc(Nspace * sizeof(double));
	atmos.Nspace = Nspace;

	for (int i = 0; i < N; i++)
		for (int j = 0; j < i; j++) c[i][j] = c[j][i] = 1.0 + uniform();
	detailed_balance(N, c, p, R);

	for (int k = 0; k < Nspace; k++) {
		atom.ntotal[k] = 1.0;
		for (int i = 0; i < N; i++) {
			atom.n[i][k] = p[i];
			for (int j = 0; j < N; j++) {
				atom.Gamma[j * N + i][k] = (i == j) ? 0.0 : R[i][j];
				atom.C[j * N + i][k]     = 0.0;
			}
		}
	}
	/* --- At k = 1 the collisional rate 1 -> 3 is negative -- ------ */

	atom.C[3 * N + 1][1] = -0.1;
	atom.Gamma[3 * N + 1][1] += -0.1;

	statEquil(&atom, -1);

	CHECK_CLOSE(atom.n[0][0], p[0] / (p[0] + p[1] + p[2] + p[3] + p[4]),
				1.0e-14, 0.0, "statEquil GTH, n[0] at k = %d", 0);
	for (int k = 0; k < Nspace; k++) {
		double sum = 0.0, res = 0.0, scale = 0.0;

		for (int i = 0; i < N; i++) {
			double in = 0.0, out = 0.0;

			sum += atom.n[i][k];
			for (int j = 0; j < N; j++) {
				if (j == i) continue;
				in += atom.n[j][k] * atom.Gamma[i * N + j][k];
				out += atom.n[i][k] * atom.Gamma[j * N + i][k];
			}
			res   = MAX(res, fabs(in - out));
			scale = MAX(scale, MAX(fabs(in), fabs(out)));
		}
		CHECK_CLOSE(sum, 1.0, 1.0e-14, 0.0, "statEquil, sum at k = %d", k);
		CHECK(res <= 1.0e-13 * scale, "statEquil, residual %.2e at k = %d",
			  res / scale, k);
		if (table)
			printf("  statEquil    k = %d   residual %.2e (%s)\n", k, res / scale,
				   (k == 0) ? "GTH" : "LU");
	}
	freeMatrix((void **) c);
	freeMatrix((void **) R);
	freeMatrix((void **) atom.Gamma);
	freeMatrix((void **) atom.C);
	freeMatrix((void **) atom.n);
	free(atom.ntotal);
}

int main(int argc, char *argv[])
{
	table = (argc > 1 && !strcmp(argv[1], "--table"));

	test_random();
	test_groups();
	test_range();
	test_status();
	test_statequil();

	return check_summary("test_gth");
}
