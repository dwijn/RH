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

/* ------- begin -------------------------- cent_deriv_pair.c ------- */

static inline void cent_deriv_pair(double odx, double dx, double yu[2],
								   double y0[2], double yd[2], double d[2])
{
	/* --- Steffen (1990) derivative of a 2-vector that turns with the
		   reference frame of linear polarization, such as (S_Q, S_U)
		   or (eta_Q, eta_U) (N-25). Limiting each component on its own
		   does not commute with turning the frame; this uses only
		   lengths and a dot product, so it does. The derivative is
		   zero where the one-sided slopes a, b point in opposite
		   directions (a . b <= 0, an extremum), and otherwise the
		   weighted centred slope p, shortened to at most twice the
		   shorter of a and b:

			 d = p min(1, 2 |a| / |p|, 2 |b| / |p|).

		   For a single component this is cent_deriv. --  -------------- */

	double a[2], b[2], p[2], na, nb, np, f;
	int    i;

	for (i = 0; i < 2; i++) {
		a[i] = (y0[i] - yu[i]) / odx;
		b[i] = (yd[i] - y0[i]) / dx;
		p[i] = (a[i] * dx + b[i] * odx) / (odx + dx);
	}
	na = hypot(a[0], a[1]);
	nb = hypot(b[0], b[1]);
	np = hypot(p[0], p[1]);

	if (a[0] * b[0] + a[1] * b[1] <= 0.0 || np == 0.0) {
		d[0] = d[1] = 0.0;
		return;
	}
	f    = fmin(1.0, 2.0 * fmin(na, nb) / np);
	d[0] = f * p[0];
	d[1] = f * p[1];
}
/* ------- end ---------------------------- cent_deriv_pair.c ------- */

/* ------- begin -------------------------- cent_deriv_mat.c -------- */

inline void cent_deriv_mat(double wprime[4][4], double dsup, double dsdn,
						   double chiup[4][4], double chic[4][4],
						   double chidn[4][4])
{
	/* --- The pairs (K[0][1], K[0][2]) and (K[1][0], K[2][0]) hold
		   (eta_Q, eta_U), (K[1][3], K[2][3]) and (K[3][1], K[3][2])
		   hold (rho_Q, rho_U) up to signs: they turn with the frame and
		   are limited as vectors. The other elements do not change
		   when the frame turns --                       -------------- */

	static const int pair[4][2][2] = {{{0, 1}, {0, 2}}, {{1, 0}, {2, 0}}, {{1, 3}, {2, 3}}, {{3, 1}, {3, 2}}};
	register int     i, j, n;
	double           yu[2], y0[2], yd[2], d[2];

	for (j = 0; j < 4; j++)
		for (i = 0; i < 4; i++)
			wprime[j][i] = cent_deriv(dsup, dsdn, chiup[j][i], chic[j][i],
									  chidn[j][i]);

	for (n = 0; n < 4; n++) {
		for (i = 0; i < 2; i++) {
			yu[i] = chiup[pair[n][i][0]][pair[n][i][1]];
			y0[i] = chic[pair[n][i][0]][pair[n][i][1]];
			yd[i] = chidn[pair[n][i][0]][pair[n][i][1]];
		}
		cent_deriv_pair(dsup, dsdn, yu, y0, yd, d);
		for (i = 0; i < 2; i++) wprime[pair[n][i][0]][pair[n][i][1]] = d[i];
	}
}
/* ------- end ---------------------------- cent_deriv_mat.c -------- */

/* ------- begin -------------------------- cent_deriv_vec.c -------- */

inline void cent_deriv_vec(double wprime[4], double dsup, double dsdn,
						   double chiup[4], double chic[4], double chidn[4])
{
	/* --- I and V on their own, (Q, U) as a vector (N-25) -- -------- */

	wprime[0] = cent_deriv(dsup, dsdn, chiup[0], chic[0], chidn[0]);
	wprime[3] = cent_deriv(dsup, dsdn, chiup[3], chic[3], chidn[3]);
	cent_deriv_pair(dsup, dsdn, chiup + 1, chic + 1, chidn + 1, wprime + 1);
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
	   coefficients of dt^1 to dt^13, from the exact expansion of the
	   closed forms below. The closed forms lose digits to cancellation
	   at small dt (1e-9 at dt = 0.05, 5e-12 at dt = 0.2, 7e-15 at
	   dt = 1), so below |dt| = 0.2 the series is used, with as many
	   terms as a relative error of about 1e-12 needs. It also serves
	   small negative dt, which interpolated opacities can give in 2-D
	   and 3-D. For dt > 40 the closed forms are expanded in 1/dt
	   (exp(-dt) negligible), so that a huge or infinite dt gives the
	   limits instead of NaN. Checked against mpmath in
	   tests/unit/test_formal_1d.c --                 -------------- */

#define BEZIER3_NTAYLOR 9

static const double bezier3_taylor[4][BEZIER3_NTAYLOR] = {
	/* alpha: 1/4, -1/5, 1/12, -1/42, 1/192, -1/1080, ... */
	{0.25, -0.20000000000000001, 0.083333333333333329, -0.023809523809523808,
	 0.005208333333333333, -0.00092592592592592596, 0.00013888888888888889, -1.8037518037518038e-05,
	 2.0667989417989419e-06},
	/* beta: 1/4, -1/20, 1/120, -1/840, 1/6720, -1/60480, ... */
	{0.25, -0.050000000000000003, 0.0083333333333333332, -0.0011904761904761906,
	 0.00014880952380952382, -1.6534391534391536e-05, 1.6534391534391535e-06, -1.5031265031265032e-07,
	 1.252605419272086e-08},
	/* gamma: 1/4, -3/20, 1/20, -1/84, 1/448, -1/2880, ... */
	{0.25, -0.14999999999999999, 0.050000000000000003, -0.011904761904761904,
	 0.002232142857142857, -0.00034722222222222224, 4.6296296296296294e-05, -5.4112554112554113e-06,
	 5.6367243867243867e-07},
	/* theta: 1/4, -1/10, 1/40, -1/210, 1/1344, -1/10080, ... */
	{0.25, -0.10000000000000001, 0.025000000000000001, -0.0047619047619047623,
	 0.00074404761904761901, -9.9206349206349206e-05, 1.1574074074074073e-05, -1.2025012025012026e-06,
	 1.1273448773448773e-07}};

/* --- c[0] + c[1] x + ... + c[6] x^6 and to c[8] x^8, with x2 = x^2
	   and x4 = x^4 (Estrin's scheme: short dependency chains) -- -- */

#define BEZIER3_ESTRIN7(c, x, x2, x4)                           \
	(((c)[0] + (c)[1] * (x)) + ((c)[2] + (c)[3] * (x)) * (x2) + \
	 (((c)[4] + (c)[5] * (x)) + (c)[6] * (x2)) * (x4))

#define BEZIER3_ESTRIN9(c, x, x2, x4)                                    \
	(((c)[0] + (c)[1] * (x)) + ((c)[2] + (c)[3] * (x)) * (x2) +          \
	 (((c)[4] + (c)[5] * (x)) + ((c)[6] + (c)[7] * (x)) * (x2)) * (x4) + \
	 (c)[8] * (x4) * (x4))

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

	const double *c0 = bezier3_taylor[0], *c1 = bezier3_taylor[1],
				 *c2 = bezier3_taylor[2], *c3 = bezier3_taylor[3];

	/* --- Below |dt| = 0.05 eps follows from the sum of the weights,
		   which is 1 - exp(-dt) (constant S). The number of terms
		   keeps the relative error below about 1e-12: in 2-D and 3-D
		   most intervals are optically very thin --     -------------- */

	if (fabs(dt) < 2.0e-6) {
		*alpha = dt * (c0[0] + c0[1] * dt);
		*beta  = dt * (c1[0] + c1[1] * dt);
		*gamma = dt * (c2[0] + c2[1] * dt);
		*theta = dt * (c3[0] + c3[1] * dt);
		*eps   = 1.0 - (*alpha + *beta + *gamma + *theta);
	} else if (fabs(dt) < 1.0e-3) {
		*alpha = dt * (c0[0] + dt * (c0[1] + dt * (c0[2] + dt * c0[3])));
		*beta  = dt * (c1[0] + dt * (c1[1] + dt * (c1[2] + dt * c1[3])));
		*gamma = dt * (c2[0] + dt * (c2[1] + dt * (c2[2] + dt * c2[3])));
		*theta = dt * (c3[0] + dt * (c3[1] + dt * (c3[2] + dt * c3[3])));
		*eps   = 1.0 - (*alpha + *beta + *gamma + *theta);
	} else if (fabs(dt) < 0.05) {
		double d2 = dt * dt, d4 = d2 * d2;

		*alpha = dt * BEZIER3_ESTRIN7(c0, dt, d2, d4);
		*beta  = dt * BEZIER3_ESTRIN7(c1, dt, d2, d4);
		*gamma = dt * BEZIER3_ESTRIN7(c2, dt, d2, d4);
		*theta = dt * BEZIER3_ESTRIN7(c3, dt, d2, d4);
		*eps   = 1.0 - (*alpha + *beta + *gamma + *theta);
	} else if (fabs(dt) < 0.2) {

		/* --- Terms to dt^9 --                           -------------- */

		double d2 = dt * dt, d4 = d2 * d2;

		*alpha = dt * BEZIER3_ESTRIN9(c0, dt, d2, d4);
		*beta  = dt * BEZIER3_ESTRIN9(c1, dt, d2, d4);
		*gamma = dt * BEZIER3_ESTRIN9(c2, dt, d2, d4);
		*theta = dt * BEZIER3_ESTRIN9(c3, dt, d2, d4);
		*eps   = exp(-dt);
	} else if (dt > 40.0) {
		double r = 1.0 / dt; /* exp(-dt) terms below 1e-17 relative */

		*eps   = 0.0;
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
