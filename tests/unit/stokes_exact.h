/* ------- file: -------------------------- stokes_exact.h ----------

	   Exact Stokes vector in a slab with constant propagation matrix
	   K = 1 + K' (in units of the Stokes I opacity) and an unpolarized
	   source vector S = (S(tau), 0, 0, 0) with S(tau) a polynomial
	   (slab_exact.h), along a ray with direction cosine mu:

		 outgoing:  I+(tau) = sum_j (mu K^-1)^j S^(j)(tau),
		 incoming, I-(tau0) = 0 at the top:
				   I-(tau) = P-(tau) - exp(-K (tau - tau0) / mu) P-(tau0),
				   P-(tau) = sum_j (-mu K^-1)^j S^(j)(tau),

	   the matrix form of the scalar solution, from mu dI/dtau = K (I - S).
	   Used by the polarized formal-solver unit tests.
	   --                                              -------------- */

#ifndef STOKES_EXACT_H
#define STOKES_EXACT_H

#include <math.h>
#include <string.h>

#include "slab_exact.h"

/* --- Reduced propagation matrix: eta_Q, eta_U, eta_V, rho_Q, rho_U,
	   rho_V over chi_I, typical of a strong Zeeman-split line core --- */

static const double K_ELEMENTS[6] = {0.12, -0.05, 0.25, 0.04, 0.09, 0.3};

typedef double Mat[4][4];
typedef double Vec[4];

static inline void mat_mul(const Mat A, const Mat B, Mat C)
{
	Mat T;
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) {
			T[i][j] = 0.0;
			for (int l = 0; l < 4; l++) T[i][j] += A[i][l] * B[l][j];
		}
	memcpy(C, T, sizeof(Mat));
}

static inline void mat_vec(const Mat A, const Vec x, Vec y)
{
	Vec t;
	for (int i = 0; i < 4; i++) {
		t[i] = 0.0;
		for (int l = 0; l < 4; l++) t[i] += A[i][l] * x[l];
	}
	memcpy(y, t, sizeof(Vec));
}

/* --- K = 1 + K' with the sign pattern of StokesK() -- ----------- */

static inline void full_K(const double e[6], double f, Mat K)
{
	double eQ = f * e[0], eU = f * e[1], eV = f * e[2], rQ = f * e[3],
		   rU = f * e[4], rV = f * e[5];
	Mat M = {{1.0, eQ, eU, eV}, {eQ, 1.0, rV, -rU}, {eU, -rV, 1.0, rQ}, {eV, rU, -rQ, 1.0}};
	memcpy(K, M, sizeof(Mat));
}

/* --- Inverse by Gauss-Jordan elimination with partial pivoting -- */

static inline void mat_inv(const Mat A, Mat Ainv)
{
	double a[4][8];
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 8; j++)
			a[i][j] = (j < 4) ? A[i][j] : (j - 4 == i);
	for (int c = 0; c < 4; c++) {
		int p = c;
		for (int r = c + 1; r < 4; r++)
			if (fabs(a[r][c]) > fabs(a[p][c])) p = r;
		for (int j = 0; j < 8; j++) {
			double t = a[c][j];
			a[c][j]  = a[p][j];
			a[p][j]  = t;
		}
		double d = a[c][c];
		for (int j = 0; j < 8; j++) a[c][j] /= d;
		for (int r = 0; r < 4; r++)
			if (r != c) {
				double f = a[r][c];
				for (int j = 0; j < 8; j++) a[r][j] -= f * a[c][j];
			}
	}
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) Ainv[i][j] = a[i][j + 4];
}

/* --- exp(X) - 1 without cancellation for small X: Taylor series
	   after scaling X by 2^-s, then (E + 1)^2 - 1 = E (E + 2) s times */

static inline void mat_expm1(const Mat X, Mat E)
{
	double norm = 0.0;
	int    s    = 0;
	Mat    Y, T;

	for (int i = 0; i < 4; i++) {
		double r = 0.0;
		for (int j = 0; j < 4; j++) r += fabs(X[i][j]);
		norm = MAX(norm, r);
	}
	while (norm > 0.5) {
		norm *= 0.5;
		s++;
	}
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) Y[i][j] = ldexp(X[i][j], -s);

	memcpy(E, Y, sizeof(Mat)); /* E = Y + Y^2/2! + ... */
	memcpy(T, Y, sizeof(Mat));
	for (int n = 2; n < 30; n++) {
		mat_mul(T, Y, T);
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++) {
				T[i][j] /= n;
				E[i][j] += T[i][j];
			}
	}
	for (; s > 0; s--) {
		Mat E2;
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++) E2[i][j] = E[i][j] + 2.0 * (i == j);
		mat_mul(E, E2, E);
	}
}

/* --- d^j/dtau^j of S = sum_n s_n tau^n, and the difference of that
	   between tau and tau0 without cancellation --   -------------- */

static inline double S_deriv(const double s[NPOLY], int j, double tau)
{
	double r = 0.0;
	for (int n = NPOLY - 1; n >= j; n--) {
		double f = 1.0;
		for (int i = 0; i < j; i++) f *= n - i;
		r = r * tau + f * s[n];
	}
	return r;
}

static inline double S_deriv_diff(const double s[NPOLY], int j, double tau, double tau0)
{
	double r = 0.0;
	for (int n = j + 1; n < NPOLY; n++) {
		double f = 1.0, d = 0.0;
		int    p = n - j;
		for (int i = 0; i < j; i++) f *= n - i;
		for (int i = 0; i < p; i++) d += pow(tau, p - 1 - i) * pow(tau0, i);
		r += f * s[n] * d;
	}
	return (tau - tau0) * r;
}

/* --- sum_j (m Kinv)^j v_j applied to e_I, v_j = d^j S --- ------- */

static inline void series(const Mat Kinv, double m, const double v[NPOLY], Vec out)
{
	Vec acc = {v[NPOLY - 1], 0.0, 0.0, 0.0};
	for (int j = NPOLY - 2; j >= 0; j--) {
		mat_vec(Kinv, acc, acc);
		for (int n = 0; n < 4; n++) acc[n] *= m;
		acc[0] += v[j];
	}
	memcpy(out, acc, sizeof(Vec));
}

static inline void I_exact_pol(const double s[NPOLY], const Mat K, double tau,
							   double tau0, double mu, bool_t to_obs, Vec I)
{
	Mat    Kinv;
	double v[NPOLY];

	mat_inv(K, Kinv);
	if (to_obs) {
		for (int j = 0; j < NPOLY; j++) v[j] = S_deriv(s, j, tau);
		series(Kinv, mu, v, I);
	} else {
		Vec Pdiff, Ptop;
		Mat X, E;
		for (int j = 0; j < NPOLY; j++) v[j] = S_deriv_diff(s, j, tau, tau0);
		series(Kinv, -mu, v, Pdiff);
		for (int j = 0; j < NPOLY; j++) v[j] = S_deriv(s, j, tau0);
		series(Kinv, -mu, v, Ptop);
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++) X[i][j] = -K[i][j] * (tau - tau0) / mu;
		mat_expm1(X, E);
		mat_vec(E, Ptop, Ptop);
		for (int n = 0; n < 4; n++) I[n] = Pdiff[n] - Ptop[n];
	}
}

#endif /* !STOKES_EXACT_H */
