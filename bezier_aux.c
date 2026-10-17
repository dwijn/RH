/* ------- file: -------------------------- bezier_aux.c ------------

   Auxiliary routines for Cubic DELO-Bezier (polarized) and
   cubic short-char Bezier solvers.

   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
			   (Derivatives) Fritsch & Butland (1984),

   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

   Modifications:
		   2017-03-12, JdlCR: Created!

   Adapted for RH:
		   Han Uitenbroek

	   Last modified: Tue Apr 28 14:30:54 2020 --

	   --------------------------                      ----------RH-- */

#include <math.h>
#include <string.h> // memcpy, memset

#include "bezier.h"

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- solveLiearFast.c -------- */

void solveLinearFast(double A[4][4], double B[4])
{

	// --- Simple Gaussian elimination with partial pivoting --- //

	// A is the matrix with coeffs that multiply X. B is the right hand
	// term (on input). The result is overwritten into B. All operations
	// are done in-place
	// It uses the macro "swap" to swap the values of two rows when a pivot
	// is found

	// The algorithm is a pretty standard Gauss-Jordan algorithm, here optimized a
	// bit for performance.

	register int i, j, k;
	int          maxrow, swapme;
	double       maxel, tmp;

	for (i = 0; i < 4; i++) {

		// Find pivot

		maxel  = fabs(A[i][i]);
		maxrow = i, swapme = 0;

		for (k = i + 1; k < 4; k++) {
			tmp = fabs(A[k][i]);
			if (tmp > maxel) {
				maxel  = tmp;
				maxrow = k, swapme = 1;
			}
		}

		// swap
		if (swapme) {
			for (k = i; k < 4; k++) SWAP(A[maxrow][k], A[i][k], tmp);
			SWAP(B[maxrow], B[i], tmp);
		}

		// Set to zero relevant columns
		for (k = i + 1; k < 4; k++) {
			tmp     = -A[k][i] / A[i][i];
			A[k][i] = 0.0;
			for (j = i + 1; j < 4; j++) {
				A[k][j] += tmp * A[i][j];
			}
			B[k] += tmp * B[i];
		}
	}

	// Solve upper triagonal system and store in B

	for (i = 3; i >= 0; i--) {
		B[i] /= A[i][i];
		for (k = i - 1; k >= 0; k--) {
			B[k] -= A[k][i] * B[i];
		}
	}
}

/* ------- end ---------------------------- solveLinearFast.c ------- */

/* ------- begin -------------------------- cent_deriv.c ------------ */

inline double cent_deriv(double odx, double dx,
						 double yu, double y0, double yd)
{
	/* --- Derivatives from Steffen (1990) --- */

	const double S0 = (yd - y0) / dx;
	const double Su = (y0 - yu) / odx;
	const double P0 = fabs((Su * dx + S0 * odx) / (odx + dx)) * 0.5;

	return (signFortran(S0) + signFortran(Su)) *
		   fmin(fabs(Su), fmin(fabs(S0), P0));
}
/* ------- end ---------------------------- cent_deriv.c ------------ */

/* ------- begin -------------------------- cent_deriv_mat.c -------- */

inline void cent_deriv_mat(double wprime[4][4], double dsup, double dsdn,
						   double chiup[4][4], double chic[4][4],
						   double chidn[4][4])
{
	register int i, j;

	for (j = 0; j < 4; j++)
		for (i = 0; i < 4; i++)
			wprime[j][i] = cent_deriv(dsup, dsdn, chiup[j][i], chic[j][i],
									  chidn[j][i]);
}
/* ------- end ---------------------------- cent_deriv_mat.c -------- */

/* ------- begin -------------------------- cent_deriv_vec.c -------- */

inline void cent_deriv_vec(double wprime[4], double dsup, double dsdn,
						   double chiup[4], double chic[4], double chidn[4])
{
	register int i;

	for (i = 0; i < 4; i++)
		wprime[i] = cent_deriv(dsup, dsdn, chiup[i], chic[i], chidn[i]);
}
/* ------- end ---------------------------- cent_deriv_vec.c -------- */

/* ------- begin -------------------------- m4m.c ------------------- */

inline void m4m(double a[4][4], double b[4][4], double c[4][4])
{

	/* --- Matrix multiplication --                      -------------- */

	register int i, j, k;
	memset(&c[0][0], 0, sizeof(double) * 16);
	for (j = 0; j < 4; j++)
		for (i = 0; i < 4; i++)
			for (k = 0; k < 4; k++)
				c[j][i] += a[k][i] * b[j][k];
}
/* ------- end ---------------------------- m4m.c ------------------- */

/* ------- begin -------------------------- m4v.c ------------------- */

/* --- Matrix/vector multiplication.
	   We use matrix as float to be able to use Intel's
	   matrix inversion --                         ------------------ */

inline void m4v(float a[4][4], double b[4], double c[4])
{
	register int k, i;
	memset(&c[0], 0, sizeof(double) * 4);
	for (i = 0; i < 4; i++)
		for (k = 0; k < 4; k++)
			c[i] += ((double) a[i][k]) * b[k];
}
/* ------- end ---------------------------- m4v.c ------------------- */

/* ------- begin -------------------------- Svec.c ------------------ */

inline void Svec(int k, double **S, double *Sf)
{
	/* --- Extracts the Source vector at depth-points k -- ------------ */

	Sf[0] = S[0][k], Sf[1] = S[1][k], Sf[2] = S[2][k], Sf[3] = S[3][k];
}
/* ------- end ---------------------------- Svec.c ------------------ */

/* ------- begin -------------------------- Bezier3_coeffs ---------- */

/* --- Taylor coefficients of alpha, beta, gamma and theta: the
	   coefficients of dt^1 to dt^18, from the exact expansion of the
	   closed forms below. For dt < 1 the series is accurate to
	   rounding, while the closed forms lose digits to cancellation
	   (1e-9 at dt = 0.05, 7e-15 at dt = 1); the series also serves
	   small negative dt, which interpolated opacities can give in 2-D
	   and 3-D. For dt > 50 the closed forms are expanded in 1/dt, so
	   that a huge or infinite dt gives the limits instead of NaN.
	   Checked against mpmath in tests/unit/test_formal_1d.c -- ---- */

#define BEZIER3_NTAYLOR 18

static const double bezier3_taylor[4][BEZIER3_NTAYLOR] = {
	/* alpha: 1/4, -1/5, 1/12, -1/42, 1/192, -1/1080, ... */
	{0.25, -0.20000000000000001, 0.083333333333333329, -0.023809523809523808,
	 0.005208333333333333, -0.00092592592592592596, 0.00013888888888888889,
	 -1.8037518037518038e-05, 2.0667989417989419e-06, -2.1197937864604531e-07,
	 1.9683799445704209e-08, -1.6701405590294479e-09, 1.3047973117417563e-10,
	 -9.4464963746009497e-12, 6.3726364431831805e-13, -4.0248230167472717e-14,
	 2.3897386661936928e-15, -1.3387891687359624e-16},
	/* beta: 1/4, -1/20, 1/120, -1/840, 1/6720, -1/60480, ... */
	{0.25, -0.050000000000000003, 0.0083333333333333332, -0.0011904761904761906,
	 0.00014880952380952382, -1.6534391534391536e-05, 1.6534391534391535e-06,
	 -1.5031265031265032e-07, 1.252605419272086e-08, -9.6354263020929685e-10,
	 6.882447358637835e-11, -4.5882982390918896e-12, 2.867686399432431e-13,
	 -1.6868743526073125e-14, 9.3715241811517366e-16, -4.932381147974598e-17,
	 2.466190573987299e-18, -1.1743764638034756e-19},
	/* gamma: 1/4, -3/20, 1/20, -1/84, 1/448, -1/2880, ... */
	{0.25, -0.14999999999999999, 0.050000000000000003, -0.011904761904761904,
	 0.002232142857142857, -0.00034722222222222224, 4.6296296296296294e-05,
	 -5.4112554112554113e-06, 5.6367243867243867e-07, -5.2994844661511327e-08,
	 4.5424152567009714e-09, -3.5788726264916739e-10, 2.6095946234835123e-11,
	 -1.7712180702376782e-12, 1.1245829017382083e-13, -6.7080383612454531e-15,
	 3.7732715782005672e-16, -2.0081837531039433e-17},
	/* theta: 1/4, -1/10, 1/40, -1/210, 1/1344, -1/10080, ... */
	{0.25, -0.10000000000000001, 0.025000000000000001, -0.0047619047619047623,
	 0.00074404761904761901, -9.9206349206349206e-05, 1.1574074074074073e-05,
	 -1.2025012025012026e-06, 1.1273448773448773e-07, -9.6354263020929689e-09,
	 7.5706920945016183e-10, -5.5059578869102679e-11, 3.7279923192621606e-12,
	 -2.3616240936502375e-13, 1.4057286271727604e-14, -7.8918098367593568e-16,
	 4.1925239757784085e-17, -2.1138776348462563e-18}};

inline void Bezier3_coeffs(double dt, double *alpha, double *beta,
						   double *gamma, double *theta, double *eps)
{
	/* ---
	   Integration coeffs. for cubic Bezier interpolants
	   Use Taylor expansion if dtau is small
	   Note:
		  alpha -> Su
	  beta  -> Sc
	  Gamma -> Cc
	  theta -> Cu

	   --- */

	if (fabs(dt) < 1.0) {
		double a = 0.0, b = 0.0, g = 0.0, t = 0.0;
		int    nterm = BEZIER3_NTAYLOR;

		/* --- Terms needed for rounding accuracy: up to dt^9 below
			   0.05, dt^13 below 0.25, all 18 below 1. Below 0.05 also
			   exp(-dt) from its series, to dt^9 --         -------------- */

		if (fabs(dt) < 0.05) {
			static const double inv_n[10] = {0.0, 1.0, 1.0 / 2, 1.0 / 3,
											 1.0 / 4, 1.0 / 5, 1.0 / 6,
											 1.0 / 7, 1.0 / 8, 1.0 / 9};
			double              e         = 1.0;

			for (int n = 9; n >= 1; n--)
				e = 1.0 - dt * inv_n[n] * e;
			*eps  = e;
			nterm = 9;
		} else {
			*eps = exp(-dt);
			if (fabs(dt) < 0.25) nterm = 13;
		}

		for (int n = nterm - 1; n >= 0; n--) {
			a = a * dt + bezier3_taylor[0][n];
			b = b * dt + bezier3_taylor[1][n];
			g = g * dt + bezier3_taylor[2][n];
			t = t * dt + bezier3_taylor[3][n];
		}
		*alpha = a * dt;
		*beta  = b * dt;
		*gamma = g * dt;
		*theta = t * dt;
	} else if (dt > 50.0) {
		double r = 1.0 / dt; /* exp(-dt) terms below 1e-17 */

		*eps   = exp(-dt);
		*alpha = 6.0 * r * r * r;
		*beta  = 1.0 + r * (-3.0 + r * (6.0 - 6.0 * r));
		*gamma = r * r * (6.0 - 18.0 * r);
		*theta = r * (3.0 + r * (-12.0 + 18.0 * r));
	} else {
		double dt2 = dt * dt, dt3 = dt2 * dt;

		*eps   = exp(-dt);
		*alpha = -(-6.0 + (6.0 + 6.0 * dt + 3.0 * dt2 + dt3) * (*eps)) / dt3;
		*beta  = (-6.0 + dt * (6.0 + dt * (dt - 3.0)) + 6.0 * (*eps)) / dt3;
		*gamma = 3.0 * (2.0 * dt - 6.0 + *eps * (6.0 + dt * (dt + 4.0))) / dt3;
		*theta = 3.0 * (-2.0 * (*eps) * (dt + 3.0) + 6.0 + (dt - 4.0) * dt) / dt3;
	}
}
/* ------- end ---------------------------- Bezier3_coeffs ---------- */

/* ------- begin -------------------------- m4inv.c ----------------- */

void m4inv(double MI[4][4])
{

	/* --- In-place Shipley-Coleman matrix inversion
		   Fast, but ... how accurate??
		   Pivoting is always done in the diagonal.
		   Copied here just in case the SIMD matrix inversion
		   gives troubles. --                            -------------- */

	register int k, i, j;

	for (k = 0; k < 4; k++) {

		/* --- The pivot element --                        -------------- */

		MI[k][k] = -1.0 / MI[k][k];

		/* --- The pivot column --                         -------------- */

		for (i = 0; i < 4; ++i)
			if (i != k) MI[i][k] *= MI[k][k];

		/* --- Elements not in a pivot row or column --    -------------- */

		for (i = 0; i < 4; ++i) {
			if (i != k)
				for (j = 0; j < 4; ++j)
					if (j != k)
						MI[i][j] += MI[i][k] * MI[k][j];
		}
		/* --- Elements in a pivot row --                  -------------- */

		for (i = 0; i < 4; ++i) {
			if (i != k)
				MI[k][i] *= MI[k][k];
		}
	}

	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j) MI[i][j] = -MI[i][j];
	}
	return;
}
/* ------- end ---------------------------- m4inv.c ----------------- */
