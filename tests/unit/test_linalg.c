/* ------- file: -------------------------- test_linalg.c -----------

	   Linear solver of ludcmp.c (LU decomposition with implicit
	   scaling and partial pivoting) on systems with known solutions:

	   -- systems built as A x_known with random, diagonally dominant A
		  (N = 1 .. 60), where the solution is accurate to a few times
		  the machine precision;
	   -- a matrix with zeros on the diagonal that requires pivoting,
		  with integer entries so that the answer and the determinant
		  (from the sign d and the diagonal of the LU decomposition)
		  are exact;
	   -- the Hilbert matrix scaled to integer entries (lcm(1..2N-1)/
		  (i+j+1)), condition number 1.5e10 at N = 8, with and without
		  the `improve' option (iterative refinement); the error is
		  ~cond * eps either way;
	   -- the inverse by LUbacksubst on the columns of the identity,
		  A A^-1 = 1 (one decomposition for several right-hand sides).

	   Run with --table to print the measured errors.
	   --                                              -------------- */

#include <math.h>
#include <string.h>

#include "rh.h"
#include "check.h"

void LUdecomp(int N, double **A, int *index, double *d);
void LUbacksubst(int N, double **A, int *index, double *b);

static bool_t table;

/* --- Deterministic pseudo-random numbers in [-1, 1) -- ------------ */

static unsigned long long rng_state = 88172645463325252ULL;

static double uniform(void)
{
	rng_state ^= rng_state << 13;
	rng_state ^= rng_state >> 7;
	rng_state ^= rng_state << 17;
	return 2.0 * (double) (rng_state >> 11) / 9007199254740992.0 - 1.0;
}

static double max_abs_diff(int N, const double *a, const double *b)
{
	double m = 0.0;

	for (int i = 0; i < N; i++) m = MAX(m, fabs(a[i] - b[i]));
	return m;
}

/* --- Random diagonally dominant systems -- ------------------------ */

/* --- Measured with gcc on x86-64: max |dx| 9e-16 at N = 60 (x is of
	   order 1), less for smaller N; bound 2e-15 --                 - */

static void test_random(void)
{
	static const int sizes[] = {1, 2, 3, 5, 10, 20, 40, 60};

	for (size_t s = 0; s < sizeof(sizes) / sizeof(sizes[0]); s++) {
		int      N = sizes[s];
		double **A = matrix_double(N, N), **A0 = matrix_double(N, N);
		double  *x = malloc(N * sizeof(double)), *b = malloc(N * sizeof(double));
		double   err[2];

		for (int improve = 0; improve < 2; improve++) {
			for (int i = 0; i < N; i++) {
				x[i] = uniform();
				for (int j = 0; j < N; j++) A0[i][j] = uniform();
				A0[i][i] += (i % 2) ? N : -N; /* dominant, either sign */
			}
			for (int i = 0; i < N; i++) {
				b[i] = 0.0;
				for (int j = 0; j < N; j++) b[i] += A0[i][j] * x[j];
			}
			memcpy(A[0], A0[0], N * N * sizeof(double));
			SolveLinearEq(N, A, b, improve);
			err[improve] = max_abs_diff(N, b, x);
			CHECK(err[improve] < 2e-15, "N = %d, improve = %d: max |dx| = %.3g",
				  N, improve, err[improve]);
		}
		if (table)
			printf("random N = %2d: error %.2e (improve %.2e)\n", N, err[0], err[1]);
		free(x);
		free(b);
		freeMatrix((void **) A);
		freeMatrix((void **) A0);
	}
}

/* --- Pivoting: zero diagonal elements, exact integer arithmetic -- - */

static void test_pivot(void)
{
	/* det = 0*(3-0) - 2*(3-0) + 1*(0-2) = -8 */
	static const double Ainit[3][3] = {{0, 2, 1}, {1, 1, 0}, {2, 0, 3}};
	static const double xk[3]       = {3, -2, 5};
	double              b[3]        = {1, 1, 1}; /* A xk, below */
	double              d, det = 1.0, err;
	int                 index[3];
	double            **A = matrix_double(3, 3);

	for (int i = 0; i < 3; i++) {
		b[i] = 0.0;
		for (int j = 0; j < 3; j++) b[i] += Ainit[i][j] * xk[j];
	}
	CHECK(b[0] == 1.0 && b[1] == 1.0 && b[2] == 21.0, "b = %g %g %g", b[0], b[1], b[2]);

	for (int improve = 0; improve < 2; improve++) {
		double bb[3] = {b[0], b[1], b[2]};

		memcpy(A[0], Ainit, sizeof(Ainit));
		SolveLinearEq(3, A, bb, improve);
		err = max_abs_diff(3, bb, xk);
		CHECK(err < 1e-14, "improve = %d: max |dx| = %.3g", improve, err);
	}

	/* --- The determinant from the decomposition -- ------------- */

	memcpy(A[0], Ainit, sizeof(Ainit));
	LUdecomp(3, A, index, &d);
	for (int i = 0; i < 3; i++) det *= A[i][i];
	CHECK_CLOSE(d * det, -8.0, 1e-14, 0.0, "determinant, sign d = %g", d);

	/* --- A permutation matrix: det = -1, solution is a permutation -- */

	{
		static const double P[3][3] = {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}};
		double              bb[3]   = {7, 8, 9}; /* x = (8, 9, 7) */

		memcpy(A[0], P, sizeof(P));
		SolveLinearEq(3, A, bb, FALSE);
		CHECK(bb[0] == 8.0 && bb[1] == 9.0 && bb[2] == 7.0,
			  "permutation matrix: x = %g %g %g", bb[0], bb[1], bb[2]);
	}

	/* --- A very badly scaled row: implicit scaling must not be
		   fooled into pivoting on the big row -- ------------------ */

	{
		static const double S[2][2] = {{1e12, 1e12}, {1.0, 2.0}};
		double              bb[2]   = {3e12, 5.0};
		double            **A2      = matrix_double(2, 2);

		/* 1e12 (x + y) = 3e12 -> x + y = 3; x + 2y = 5 -> y = 2, x = 1 */
		memcpy(A2[0], S, sizeof(S));
		SolveLinearEq(2, A2, bb, FALSE);
		CHECK_CLOSE(bb[0], 1.0, 1e-13, 0.0, "scaled rows: n = %d", 0);
		CHECK_CLOSE(bb[1], 2.0, 1e-13, 0.0, "scaled rows: n = %d", 1);
		freeMatrix((void **) A2);
	}
	freeMatrix((void **) A);
}

/* --- Scaled Hilbert matrix and the `improve' flag -- -------------- */

/* --- Measured (gcc, x86-64): max |dx| 1.4e-7 without and 8.9e-7 with
	   the refinement. The error is ~cond * eps and depends on the
	   rounding (FMA), so the bound is 1e-5. The refinement does not
	   help: the residual is evaluated in plain double precision, and
	   for the computed (non-integer) solution it carries an error
	   ~ eps |A| |x| that cond(A) amplifies by as much as the first
	   solution had -- the flag is not an accuracy improvement for
	   ill-conditioned systems --                                   - */

static void test_hilbert(void)
{
	enum { N = 8 };
	const double lcm = 360360.0; /* lcm(1, ..., 2N-1) */
	double       x[N], b[N], bb[N], err[2];
	double     **A = matrix_double(N, N);

	for (int i = 0; i < N; i++) x[i] = (i % 2) ? -(i + 1) : (i + 1);

	for (int improve = 0; improve < 2; improve++) {
		for (int i = 0; i < N; i++) {
			b[i] = 0.0;
			for (int j = 0; j < N; j++) {
				A[i][j] = lcm / (i + j + 1); /* exact integer */
				b[i] += A[i][j] * x[j];      /* exact: integers < 2^53 */
			}
			bb[i] = b[i];
		}
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				CHECK(A[i][j] == floor(A[i][j]), "A[%d][%d] = %g is not an integer", i, j, A[i][j]);
		SolveLinearEq(N, A, bb, improve);
		err[improve] = max_abs_diff(N, bb, x);
	}
	CHECK(err[0] < 1e-5, "without refinement: max |dx| = %.3g", err[0]);
	CHECK(err[1] < 1e-5, "with refinement: max |dx| = %.3g", err[1]);
	if (table)
		printf("Hilbert N = %d: error %.2e, improved %.2e\n", N, err[0], err[1]);
	freeMatrix((void **) A);
}

/* --- Inverse by repeated back substitution -- --------------------- */

/* --- Measured: max |A A^-1 - 1| = 4.4e-16 for N = 12; bound 1e-15 -- */

static void test_inverse(void)
{
	enum { N = 12 };
	double **A = matrix_double(N, N), **A0 = matrix_double(N, N),
		   **Ainv = matrix_double(N, N);
	double d, maxdev = 0.0, col[N];
	int    index[N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) A0[i][j] = uniform();
		A0[i][i] += (i % 2) ? N : -N;
	}
	memcpy(A[0], A0[0], N * N * sizeof(double));
	LUdecomp(N, A, index, &d);
	for (int j = 0; j < N; j++) {
		for (int i = 0; i < N; i++) col[i] = (i == j) ? 1.0 : 0.0;
		LUbacksubst(N, A, index, col);
		for (int i = 0; i < N; i++) Ainv[i][j] = col[i];
	}
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++) {
			double s = 0.0;
			for (int k = 0; k < N; k++) s += A0[i][k] * Ainv[k][j];
			maxdev = MAX(maxdev, fabs(s - ((i == j) ? 1.0 : 0.0)));
		}
	CHECK(maxdev < 1e-15, "max |A A^-1 - 1| = %.3g", maxdev);
	if (table) printf("inverse N = %d: max |A A^-1 - 1| = %.2e\n", N, maxdev);
	freeMatrix((void **) A);
	freeMatrix((void **) A0);
	freeMatrix((void **) Ainv);
}

int main(int argc, char *argv[])
{
	table = (argc > 1 && strcmp(argv[1], "--table") == 0);

	test_random();
	test_pivot();
	test_hilbert();
	test_inverse();

	return check_summary("test_linalg");
}
