/* ------- file: -------------------------- mpe.c -------------------

	   Minimal polynomial extrapolation (MPE) of a vector sequence
	   V_0, ..., V_K: the estimate of its limit used by ERROR_LIMIT
	   (N-3, docs/SOLVER.md).

	   --------------------------                      ----------RH-- */

/* --- With the differences D_j = V_(j+1) - V_j, j = 0 .. K-1, MPE
	   takes the coefficients c_0 .. c_(K-2) that minimize

		 | sum_(j < K-1) c_j D_j + D_(K-1) |   (Euclidean),

	   sets c_(K-1) = 1, gamma_j = c_j / sum_l c_l, and estimates the
	   limit as s = sum_(j < K) gamma_j V_j. For V_j = s + M^j e
	   (a linear iteration) s is exact once K - 1 is at least the degree
	   of the minimal polynomial of M on e and 1 is not an eigenvalue
	   of M, whatever the others (also > 1 in modulus, where the
	   iteration diverges).

	   The least squares is solved by modified Gram-Schmidt with one
	   reorthogonalization: near the solution the differences are
	   close to parallel. A column that is a combination of the earlier
	   ones (to 1e-12 of its norm) gets c = 0.

	   Only the differences are kept, K vectors of N:

		 s - V_0 = sum_(l = 1 .. K-1) gamma_l (V_l - V_0)
				 = sum_(j < K-1) beta_j D_j,  beta_j = sum_(l > j) gamma_l,

	   and with D_j = sum_(q <= j) Q_q R_(qj) from the QR, s - V_0 =
	   sum_q Q_q a_q with a_q = sum_(j >= q) R_(qj) beta_j. --------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "mpe.h"

#define MPE_DEPENDENT 1.0E-12

static double dot(long N, const double *a, const double *b)
{
	register long i;
	double        s = 0.0;

	for (i = 0; i < N; i++) s += a[i] * b[i];
	return s;
}

/* ------- begin -------------------------- mpeExtrapolate.c -------- */

bool_t mpeExtrapolate(int K, long N, double **D)
{
	register long i;
	register int  l, q, pass;

	int     Kc = K - 1;
	bool_t  ok = TRUE;
	double *R, *c, *a, *beta, s, nrm0, nrm, sum;
	bool_t *alive;

	/* --- D[0 .. K-1]: the differences D_j, overwritten. Returns TRUE
		   with s - V_0 in D[K-1], or FALSE (D[K-1] unchanged) when
		   sum_l c_l vanishes (no finite limit) or a coefficient is not
		   finite. K >= 2 --                            -------------- */

	if (K < 2) return FALSE;

	R     = (double *) calloc((size_t) Kc * Kc, sizeof(double));
	c     = (double *) calloc(Kc + 1, sizeof(double));
	a     = (double *) calloc(Kc, sizeof(double));
	beta  = (double *) calloc(Kc, sizeof(double));
	alive = (bool_t *) calloc(Kc, sizeof(bool_t));

	/* --- QR of the columns D_0 .. D_(K-2) in place -- ------------- */

	for (l = 0; l < Kc; l++) {
		nrm0 = sqrt(dot(N, D[l], D[l]));
		for (pass = 0; pass < 2; pass++) {
			for (q = 0; q < l; q++) {
				if (!alive[q]) continue;
				s = dot(N, D[q], D[l]);
				R[q * Kc + l] += s;
				for (i = 0; i < N; i++) D[l][i] -= s * D[q][i];
			}
		}
		nrm = sqrt(dot(N, D[l], D[l]));
		if (nrm > MPE_DEPENDENT * nrm0 && nrm > 0.0) {
			alive[l]      = TRUE;
			R[l * Kc + l] = nrm;
			for (i = 0; i < N; i++) D[l][i] /= nrm;
		} else {
			memset(D[l], 0, N * sizeof(double));
		}
	}
	/* --- R c = -Q^T D_(K-1), back substitution -- ----------------- */

	for (l = Kc - 1; l >= 0; l--) {
		if (!alive[l]) continue;
		s = -dot(N, D[l], D[Kc]);
		for (q = l + 1; q < Kc; q++) s -= R[l * Kc + q] * c[q];
		c[l] = s / R[l * Kc + l];
	}
	c[Kc] = 1.0;
	for (sum = 0.0, l = 0; l <= Kc; l++) sum += c[l];
	if (!(fabs(sum) > 1.0E-14) || !isfinite(sum)) ok = FALSE;

	if (ok) {
		/* --- beta_j = sum_(l > j) gamma_l, a_q = sum_(j >= q) R_qj
			   beta_j --                                  -------------- */

		for (s = 0.0, l = Kc; l >= 1; l--) {
			s += c[l] / sum;
			beta[l - 1] = s;
		}
		for (q = 0; q < Kc; q++) {
			for (a[q] = 0.0, l = q; l < Kc; l++) a[q] += R[q * Kc + l] * beta[l];
			if (!isfinite(a[q])) ok = FALSE;
		}
	}
	if (ok) {
		for (i = 0; i < N; i++) {
			for (s = 0.0, q = 0; q < Kc; q++) s += D[q][i] * a[q];
			D[Kc][i] = s;
		}
	}
	free(R);
	free(c);
	free(a);
	free(beta);
	free(alive);

	return ok;
}
/* ------- end ---------------------------- mpeExtrapolate.c -------- */
