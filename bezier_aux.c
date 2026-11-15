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

#include <complex.h>
#undef I /* complex.h's imaginary unit; RH uses I for intensities */
#include <math.h>
#include <string.h> // memcpy, memset

#include "rh.h"
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

	double a[2], b[2], p[2], na2, nb2, np2, f;
	int    i;

	for (i = 0; i < 2; i++) {
		a[i] = (y0[i] - yu[i]) / odx;
		b[i] = (yd[i] - y0[i]) / dx;
		p[i] = (a[i] * dx + b[i] * odx) / (odx + dx);
	}
	na2 = a[0] * a[0] + a[1] * a[1];
	nb2 = b[0] * b[0] + b[1] * b[1];
	np2 = p[0] * p[0] + p[1] * p[1];

	if (a[0] * b[0] + a[1] * b[1] <= 0.0 || np2 == 0.0) {
		d[0] = d[1] = 0.0;
		return;
	}
	f    = fmin(1.0, 2.0 * sqrt(fmin(na2, nb2) / np2));
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
	static const int single[8][2]  = {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {0, 3}, {3, 0}, {1, 2}, {2, 1}};
	register int     i, j, n;
	double           yu[2], y0[2], yd[2], d[2];

	for (n = 0; n < 8; n++) {
		i            = single[n][0];
		j            = single[n][1];
		wprime[i][j] = cent_deriv(dsup, dsdn, chiup[i][j], chic[i][j], chidn[i][j]);
	}

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

/* --- Largest |K'| dtau below which a DELO step does not amplify the
	   Stokes vector (B-55). With S = 0 the norm of I cannot grow for a
	   physical K (|eta| <= 1); measured worst |I_loc| / |I_uw| over
	   400 000 steps with random physical K' at each end (|rho| up to
	   30, |eta| < 0.95) and dtau 1e-3 to 1e2, binned by
	   x = max |K'| dtau: the Bezier step 0.986 for x <= 1, above 1 up
	   to 1.4 for x <= 3 and 1.4e4 beyond; DELO-linear 0.998 for
	   x <= 0.25, 1.002 at 0.3, 8 beyond; the exponential step never
	   above 0.999 --                                  -------------- */

#define STOKES_XMAX_BEZIER3 1.0
#define STOKES_XMAX_LINEAR  0.2

static double maxKdtau(double dtau, double K_uw[4][4], double K_loc[4][4])
{
	double m = 0.0;

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) m = fmax(m, fmax(fabs(K_uw[i][j]), fabs(K_loc[i][j])));
	return m * dtau;
}

/* --- Use the exponential step: where max |K'| dtau exceeds xmax and
	   K' is physical at both ends, |eta| <= 1 (dichroism not larger
	   than the absorption in I). Where |eta| > 1 a Zeeman component
	   has negative opacity (population inversion): the transfer
	   equation itself amplifies there, and the exponential step would
	   reproduce that exactly, so the DELO step is kept, as before
	   B-55 --                                          -------------- */

static bool_t use_exp_step(double dtau, double K_uw[4][4], double K_loc[4][4],
						   double xmax)
{
	if (maxKdtau(dtau, K_uw, K_loc) <= xmax) return FALSE;
	for (int e = 0; e < 2; e++) {
		double(*K)[4] = e ? K_loc : K_uw;
		if (K[0][1] * K[0][1] + K[0][2] * K[0][2] + K[0][3] * K[0][3] > 1.0) return FALSE;
	}
	return TRUE;
}

/* ------- begin -------------------------- invert4.c --------------- */

static void invert4(double A[4][4], double Ainv[4][4])
{
	/* --- Ainv = A^-1 by Gauss-Jordan elimination with partial
		   pivoting; for the thick steps of exp_step_const, where A is
		   well conditioned (its eigenvalues have real parts below -40)
		   --                                              -------------- */

	double M[4][8];

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) {
			M[i][j]     = A[i][j];
			M[i][4 + j] = (i == j);
		}
	for (int c = 0; c < 4; c++) {
		int p = c;
		for (int i = c + 1; i < 4; i++)
			if (fabs(M[i][c]) > fabs(M[p][c])) p = i;
		if (p != c)
			for (int j = 0; j < 8; j++) {
				double t;
				SWAP(M[c][j], M[p][j], t);
			}
		double piv = 1.0 / M[c][c];
		for (int j = 0; j < 8; j++) M[c][j] *= piv;
		for (int i = 0; i < 4; i++) {
			if (i == c || M[i][c] == 0.0) continue;
			double f = M[i][c];
			for (int j = 0; j < 8; j++) M[i][j] -= f * M[c][j];
		}
	}
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) Ainv[i][j] = M[i][4 + j];
}
/* ------- end ---------------------------- invert4.c --------------- */

/* ------- begin -------------------------- exp_phi4.c -------------- */

static void exp_phi4(double A[4][4], double E[4][4], double F[4][4][4])
{
	/* --- E = exp(A) and F[k-1] = phi_k(A), k = 1 .. 4, of a 4 x 4
		   matrix, phi_k(A) = sum_n A^n / (n + k)!. With E(t) = exp(t A)
		   and G_k(t) = t^k phi_k(t A) (the first block row of the
		   exponential of the block matrix [[A, 1, 0, ..], [0, 0, 1, ..],
		   ..] at time t):

			 E(2t)   = E(t)^2,
			 G_k(2t) = E(t) G_k(t) + sum_(j = 1 .. k) G_j(t) t^(k-j) / (k-j)!.

		   Start from t = 2^-s with |t A| <= 1/2 (infinity norm), where 18
		   Taylor terms give E and G_k to rounding, and double s times.
		   Only 4 x 4 products: about 20 times cheaper than the 20 x 20
		   block exponential --                            -------------- */

	enum { NTERM = 18 };
	double B[4][4], P[4][4], T[4][4], G[4][4][4], Gn[4][4][4], norm = 0.0, t;
	int    s = 0;

	for (int i = 0; i < 4; i++) {
		double row = 0.0;
		for (int j = 0; j < 4; j++) row += fabs(A[i][j]);
		norm = fmax(norm, row);
	}
	if (norm > 0.5) s = (int) ceil(log2(norm / 0.5));
	t = ldexp(1.0, -s);

	/* --- Taylor series at t: E = sum B^n / n!, G_k = t^k sum B^n / (n + k)! */

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++) {
			B[i][j] = t * A[i][j];
			P[i][j] = (i == j) ? 1.0 : 0.0;
			E[i][j] = P[i][j];
			for (int k = 0; k < 4; k++) G[k][i][j] = 0.0;
		}
	{
		double fact[NTERM + 5];

		double tk[4] = {t, t * t, t * t * t, t * t * t * t};

		fact[0] = 1.0;
		for (int n = 1; n < NTERM + 5; n++) fact[n] = fact[n - 1] * n;
		for (int n = 0; n <= NTERM; n++) {
			if (n > 0) {
				for (int i = 0; i < 4; i++)
					for (int j = 0; j < 4; j++) {
						double sum = 0.0;
						for (int m = 0; m < 4; m++) sum += P[i][m] * B[m][j];
						T[i][j] = sum;
					}
				memcpy(P, T, sizeof(P));
				for (int i = 0; i < 4; i++)
					for (int j = 0; j < 4; j++) E[i][j] += P[i][j] / fact[n];
			}
			for (int k = 0; k < 4; k++) {
				double c = tk[k] / fact[n + k + 1];
				for (int i = 0; i < 4; i++)
					for (int j = 0; j < 4; j++) G[k][i][j] += c * P[i][j];
			}
		}
	}

	/* --- s doublings --                                -------------- */

	for (int d = 0; d < s; d++, t *= 2.0) {
		for (int k = 0; k < 4; k++)
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++) {
					double sum = 0.0, c = 1.0;
					for (int m = 0; m < 4; m++) sum += E[i][m] * G[k][m][j];
					/* sum over j' = k .. 1 of G_j' t^(k-j') / (k-j')! */
					for (int q = k; q >= 0; q--) {
						sum += G[q][i][j] * c;
						c *= t / (k - q + 1);
					}
					Gn[k][i][j] = sum;
				}
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++) {
				double sum = 0.0;
				for (int m = 0; m < 4; m++) sum += E[i][m] * E[m][j];
				T[i][j] = sum;
			}
		memcpy(E, T, sizeof(T));
		memcpy(G, Gn, sizeof(G));
	}
	memcpy(F, G, sizeof(G));
}
/* ------- end ---------------------------- exp_phi4.c -------------- */

/* ------- begin -------------------------- exp_phis.c -------------- */

/* --- 1 / j!, j = 0 .. 39 (correctly rounded) --      -------------- */

static const double ifact[40] = {
	1, 1, 0.5, 0.16666666666666666,
	0.041666666666666664, 0.0083333333333333332, 0.0013888888888888889, 0.00019841269841269841,
	2.4801587301587302e-05, 2.7557319223985893e-06, 2.7557319223985888e-07, 2.505210838544172e-08,
	2.08767569878681e-09, 1.6059043836821613e-10, 1.1470745597729725e-11, 7.6471637318198164e-13,
	4.7794773323873853e-14, 2.8114572543455206e-15, 1.5619206968586225e-16, 8.2206352466243295e-18,
	4.1103176233121648e-19, 1.9572941063391263e-20, 8.8967913924505741e-22, 3.8681701706306841e-23,
	1.6117375710961184e-24, 6.4469502843844736e-26, 2.4795962632247976e-27, 9.183689863795546e-29,
	3.2798892370698378e-30, 1.1309962886447716e-31, 3.7699876288159054e-33, 1.2161250415535179e-34,
	3.8003907548547434e-36, 1.1516335620771951e-37, 3.3871575355211618e-39, 9.6775929586318907e-41,
	2.6882202662866363e-42, 7.2654601791530714e-44, 1.911963205040282e-45, 4.9024697565135435e-47};

static void exp_phis_real(double z, double f[5])
{
	/* --- f[0] = e^z and f[k] = phi_k(z), k = 1 .. 4, real z: the
		   power series of phi_4 and downward phi_k = 1/k! + z phi_(k+1)
		   where |z| < 2, else e^z and upward
		   phi_(k+1) = (phi_k - 1/k!) / z --              -------------- */

	if (fabs(z) < 2.0) {
		double term = 1.0 / 24.0, sum = term;

		for (int n = 1; n <= 30; n++) {
			term *= z / (n + 4);
			sum += term;
		}
		f[4] = sum;
		f[3] = 1.0 / 6.0 + z * f[4];
		f[2] = 0.5 + z * f[3];
		f[1] = 1.0 + z * f[2];
		f[0] = 1.0 + z * f[1];
	} else {
		f[0] = exp(z);
		f[1] = (f[0] - 1.0) / z;
		f[2] = (f[1] - 1.0) / z;
		f[3] = (f[2] - 0.5) / z;
		f[4] = (f[3] - 1.0 / 6.0) / z;
	}
}

static void exp_phis(double complex z, double complex f[5])
{
	/* --- The same for complex z --                    -------------- */

	if (cabs(z) < 2.0) {
		double complex term = 1.0 / 24.0, sum = term;

		for (int n = 1; n <= 30; n++) {
			term *= z / (n + 4);
			sum += term;
		}
		f[4] = sum;
		f[3] = 1.0 / 6.0 + z * f[4];
		f[2] = 0.5 + z * f[3];
		f[1] = 1.0 + z * f[2];
		f[0] = 1.0 + z * f[1];
	} else {
		f[0] = cexp(z);
		f[1] = (f[0] - 1.0) / z;
		f[2] = (f[1] - 1.0) / z;
		f[3] = (f[2] - 0.5) / z;
		f[4] = (f[3] - 1.0 / 6.0) / z;
	}
}
/* ------- end ---------------------------- exp_phis.c -------------- */

/* ------- begin -------------------------- exp_pair.c -------------- */

static void exp_pair(double z, double s, double ge[5], double go[5])
{
	/* --- Even and odd parts of g_k(x) = f_k(z + x) (f_0 = exp,
		   f_k = phi_k) on the pair x = +-r, r = sqrt(s), s of either
		   sign (r imaginary for s < 0):

			 ge_k = [g_k(r) + g_k(-r)] / 2,  go_k = [g_k(r) - g_k(-r)] / (2 r),

		   both entire in s. Where |r| >= 1 directly. Else without the
		   division by r: ge_0 = e^z cosh r and go_0 = e^z sinh(r)/r, and
		   for k >= 1, where z <= -3, upward with
		   phi_(k+1)(w) = (phi_k(w) - 1/k!) / w written for the pair
		   (den = z^2 - s >= 8), else from the power series, which for
		   k >= 1 has no cancellation to speak of (|z| + |r| < 4; for
		   k = 0 it lost up to 1e2 ulp near z = -3) --     -------------- */

	static const double c[4] = {1.0, 1.0, 0.5, 1.0 / 6.0};
	double              m    = sqrt(fabs(s));

	if (m >= 1.0) {
		if (s > 0.0) {
			double fp[5], fm[5];

			exp_phis_real(z + m, fp);
			exp_phis_real(z - m, fm);
			for (int k = 0; k < 5; k++) {
				ge[k] = 0.5 * (fp[k] + fm[k]);
				go[k] = 0.5 * (fp[k] - fm[k]) / m;
			}
		} else {
			double complex fc[5];

			exp_phis(z + _Complex_I * m, fc);
			for (int k = 0; k < 5; k++) {
				ge[k] = creal(fc[k]);
				go[k] = cimag(fc[k]) / m;
			}
		}
		return;
	}
	double C, Sh, ez = exp(z);

	if (m < 1.0e-3) {
		C  = 1.0 + s * (0.5 + s / 24.0);
		Sh = 1.0 + s * (1.0 / 6.0 + s / 120.0);
	} else if (s > 0.0) {
		C  = cosh(m);
		Sh = sinh(m) / m;
	} else {
		C  = cos(m);
		Sh = sin(m) / m;
	}
	ge[0] = ez * C;
	go[0] = ez * Sh;

	if (z <= -3.0) {
		double den = z * z - s;

		for (int k = 0; k < 4; k++) {
			ge[k + 1] = (z * (ge[k] - c[k]) - s * go[k]) / den;
			go[k + 1] = (z * go[k] - (ge[k] - c[k])) / den;
		}
	} else {
		/* --- g_k(x) = sum_n (z + x)^n / (n + k)!, k = 1 .. 4, with
			   E_n = [(z + r)^n + (z - r)^n] / 2 and
			   O_n = [(z + r)^n - (z - r)^n] / (2 r); the terms fall
			   below 1e-18 of the first by n = 34 at the latest -- */

		const double *f = ifact;
		double        E = 1.0, O = 0.0;

		for (int k = 1; k < 5; k++) ge[k] = go[k] = 0.0;
		for (int n = 0; n < 35; n++) {
			for (int k = 1; k < 5; k++) {
				ge[k] += E * f[n + k];
				go[k] += O * f[n + k];
			}
			if ((fabs(E) + fabs(O)) * f[n + 1] < 1.0e-18) break;

			double En = z * E + s * O;

			O = z * O + E;
			E = En;
		}
	}
}
/* ------- end ---------------------------- exp_pair.c -------------- */

/* ------- begin -------------------------- exp_closed_form.c ------- */

/* --- Smallest (a^2 + b^2) / (dtau^2 (eta^2 + rho^2)) for the closed
	   form: its error grows as the inverse of that ratio (measured 8,
	   17 and 37 ulp at 0.59, 0.07 and 0.02) --         -------------- */

#define EXP_CLOSED_MIN_RATIO 0.1

static bool_t exp_closed_form(double dtau, double K[4][4], double v[5][4],
							  double I_loc[4])
{
	/* --- The sum of f_k(A) v_k, k = 0 .. 4 (f_0 = exp, f_k = phi_k),
		   A = -dtau (1 + K'), in closed form. K' = [[0, eta], [eta, R]]
		   with R antisymmetric (rho) generates a Lorentz transformation:
		   B = -dtau K' has eigenvalues +-a and +-i b, a^2 and -b^2 the
		   roots of s^2 - p s - q^2, p = dtau^2 (eta^2 - rho^2),
		   q = dtau^2 eta.rho, and (B^2 - a^2)(B^2 + b^2) = 0 (also
		   where B is not diagonalizable: B^2 has the distinct
		   eigenvalues a^2 and -b^2). With g_k(x) = f_k(-dtau + x) and
		   its even and odd parts on the two pairs (exp_pair),
		   interpolation on the spectrum gives

			 f_k(A) = al0_k + be0_k B + al1_k B^2 + be1_k B^3,

		   so that sum_k f_k(A) v_k = u0 + B (u1 + B (u2 + B u3)) with
		   u0 = sum_k al0_k v_k etc. Returns FALSE where K' is not of
		   that form or close to degenerate (|eta| = |rho|,
		   eta _|_ rho: a^2 + b^2 small against |B|^2) -- -------------- */

	double e2 = SQ(K[0][1]) + SQ(K[0][2]) + SQ(K[0][3]);
	double r2 = SQ(K[2][3]) + SQ(K[1][3]) + SQ(K[1][2]);
	double er = K[0][1] * K[2][3] - K[0][2] * K[1][3] + K[0][3] * K[1][2];

	for (int i = 0; i < 4; i++) {
		if (K[i][i] != 0.0) return FALSE;
		for (int j = i + 1; j < 4; j++)
			if (K[j][i] != ((i == 0) ? K[i][j] : -K[i][j])) return FALSE;
	}

	double p = e2 - r2, q2 = er * er, R = sqrt(0.25 * p * p + q2), A2, B2;

	if (p >= 0.0) {
		A2 = 0.5 * p + R;
		B2 = (A2 > 0.0) ? q2 / A2 : 0.0;
	} else {
		B2 = -0.5 * p + R;
		A2 = q2 / B2;
	}
	double a2 = dtau * dtau * A2, b2 = dtau * dtau * B2, D = a2 + b2;

	if (!(D >= EXP_CLOSED_MIN_RATIO * dtau * dtau * (e2 + r2)) || D == 0.0) return FALSE;

	double gea[5], goa[5], geb[5], gob[5], u[4][4] = {{0.0}}, w[4], x[4];

	exp_pair(-dtau, a2, gea, goa);
	exp_pair(-dtau, -b2, geb, gob);
	for (int k = 0; k < 5; k++) {
		double c[4] = {(gea[k] * b2 + geb[k] * a2) / D, (goa[k] * b2 + gob[k] * a2) / D,
					   (gea[k] - geb[k]) / D, (goa[k] - gob[k]) / D};

		for (int p = 0; p < 4; p++)
			for (int i = 0; i < 4; i++) u[p][i] += c[p] * v[k][i];
	}
	/* --- Horner: w = u2 + B u3, w = u1 + B w, I = u0 + B w,
		   B = -dtau K' --                                 -------------- */

	memcpy(w, u[3], sizeof(w));
	for (int p = 2; p >= 0; p--) {
		for (int i = 0; i < 4; i++)
			x[i] = u[p][i] - dtau * (K[i][0] * w[0] + K[i][1] * w[1] + K[i][2] * w[2] +
									 K[i][3] * w[3]);
		memcpy(w, x, sizeof(w));
	}
	memcpy(I_loc, w, sizeof(w));
	return TRUE;
}
/* ------- end ---------------------------- exp_closed_form.c ------- */

/* ------- begin -------------------------- exp_step_const.c -------- */

enum exp_mode { EXP_ANY,
				EXP_TAYLOR,
				EXP_CLOSED };

static bool_t exp_step_const(double dtau, double K[4][4], double S_uw[4],
							 double S_loc[4], double dS_uw[4], double dS_loc[4],
							 double I_uw[4], double I_loc[4], enum exp_mode mode)
{
	/* --- dI/dtau = S - K I, K = 1 + K', over dtau with constant K'
		   and S the cubic Hermite polynomial through S_uw, S_loc with
		   the derivatives dS_uw, dS_loc (in tau). With A = -(1 + K') dtau,
		   u = t / dtau and S = a0 + a1 u + a2 u^2 + a3 u^3,

			 I_loc = e^A I_uw + dtau sum_m a_m m! phi_(m+1)(A),

		   phi_k(A) = sum_n A^n / (n + k)!, so that K need not be
		   inverted. e^A is the exact propagator of a constant physical
		   K and never increases the norm of the Stokes vector. Exact
		   for constant K' and cubic S. By the closed form
		   (exp_closed_form) where it applies, else by the Taylor series
		   and squarings (exp_phi4), or for a thick step by one inverse;
		   mode restricts that to one of the two (for the unit test).
		   Returns FALSE only for EXP_CLOSED where the closed form does
		   not apply. I_uw and I_loc must not overlap --  -------------- */

	double A[4][4], E[4][4], F[4][4][4], a[4][4], eta2 = 0.0;

	/* --- Monomial coefficients of the Hermite cubic, already times
		   m! (0!, 1!, 2!, 3!) --                         -------------- */

	for (int n = 0; n < 4; n++) {
		double m0 = dtau * dS_uw[n], m1 = dtau * dS_loc[n];

		a[0][n] = S_uw[n];
		a[1][n] = m0;
		a[2][n] = 2.0 * (3.0 * (S_loc[n] - S_uw[n]) - 2.0 * m0 - m1);
		a[3][n] = 6.0 * (2.0 * (S_uw[n] - S_loc[n]) + m0 + m1);
	}
	for (int j = 1; j < 4; j++) eta2 += SQ(K[0][j]);
	bool_t thick = (eta2 < 1.0 && dtau * (1.0 - sqrt(eta2)) > 40.0);

	if (mode == EXP_CLOSED || (mode == EXP_ANY && !thick)) {
		double v[5][4];

		for (int n = 0; n < 4; n++) {
			v[0][n] = I_uw[n];
			for (int m = 0; m < 4; m++) v[m + 1][n] = dtau * a[m][n];
		}
		if (exp_closed_form(dtau, K, v, I_loc)) return TRUE;
		if (mode == EXP_CLOSED) return FALSE;
	}
	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			A[i][j] = -dtau * (K[i][j] + (i == j));

	/* --- Thick step: the eigenvalues of A have real parts at most
		   -dtau (1 - |eta|), so where that is below -40, e^A < 4e-18
		   and phi_k(A) = A^-k (e^A - sum_(j<k) A^j / j!) is
		   -sum_(j<k) A^(j-k) / j! to rounding: one inverse instead of
		   the series and the squarings --              -------------- */

	if (thick) {
		double              Ainv[4][4], Pw[4][4][4];
		static const double ifact4[4] = {1.0, 1.0, 0.5, 1.0 / 6.0};

		invert4(A, Ainv);
		memcpy(Pw[0], Ainv, sizeof(Ainv)); /* Pw[p] = A^-(p+1) */
		for (int q = 1; q < 4; q++)
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++) {
					double sum = 0.0;
					for (int m = 0; m < 4; m++) sum += Pw[q - 1][i][m] * Ainv[m][j];
					Pw[q][i][j] = sum;
				}
		for (int k = 1; k <= 4; k++) /* phi_k = -sum_(j<k) A^(j-k) / j! */
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++) {
					double sum = 0.0;
					for (int jj = 0; jj < k; jj++) sum -= Pw[k - jj - 1][i][j] * ifact4[jj];
					F[k - 1][i][j] = sum;
				}
		memset(E, 0, sizeof(E));
	} else
		exp_phi4(A, E, F);

	for (int i = 0; i < 4; i++) {
		I_loc[i] = 0.0;
		for (int j = 0; j < 4; j++) {
			I_loc[i] += E[i][j] * I_uw[j];
			for (int m = 0; m < 4; m++) I_loc[i] += dtau * F[m][i][j] * a[m][j];
		}
	}
	return TRUE;
}

/* --- The two paths separately, for the unit test --  -------------- */

bool_t Stokes_Exp_const_closed(double dtau, double K[4][4], double S_uw[4],
							   double S_loc[4], double dS_uw[4], double dS_loc[4],
							   double I_uw[4], double I_loc[4])
{
	return exp_step_const(dtau, K, S_uw, S_loc, dS_uw, dS_loc, I_uw, I_loc, EXP_CLOSED);
}

void Stokes_Exp_const_taylor(double dtau, double K[4][4], double S_uw[4],
							 double S_loc[4], double dS_uw[4], double dS_loc[4],
							 double I_uw[4], double I_loc[4])
{
	exp_step_const(dtau, K, S_uw, S_loc, dS_uw, dS_loc, I_uw, I_loc, EXP_TAYLOR);
}
/* ------- end ---------------------------- exp_step_const.c -------- */

/* ------- begin -------------------------- Stokes_Exp_step.c ------- */

void Stokes_Exp_step(double dtau, double K_uw[4][4], double K_loc[4][4],
					 double S_uw[4], double S_loc[4], double dS_uw[4],
					 double dS_loc[4], double I_uw[4], double I_loc[4])
{
	/* --- One step of the polarized transfer equation dI/dtau = S - K I,
		   K = 1 + K', from the upwind to the local point over dtau, for
		   the DELO steps' fallback (B-55): two exact half steps
		   (exp_step_const), K' = K'_uw on the first and K'_loc on the
		   second, with S the cubic Hermite polynomial through S_uw,
		   S_loc and the derivatives dS_uw, dS_loc, split at the
		   midpoint (exactly). The propagators never increase the norm
		   of the Stokes vector, however large the magneto-optical terms
		   are, whereas the DELO steps amplify it where max |K'| dtau is
		   large. For a thin step the product of the two propagators is
		   e^(-(1 + K'_mean) dtau) to second order; in a thick step the
		   second half gives the local equilibrium, I -> (1 + K'_loc)^-1
		   S_loc. A single step with the mean K' does not: where the
		   local dichroism exceeds the upwind one it leaves I several
		   per cent below the local source (5 % at |eta| 0.46 vs 0.12),
		   and the near-cancelled preconditioned rates of a thick line
		   core turned negative (B-55, found by an independent analysis).
		   Exact for constant K' and cubic S. I_uw and I_loc must not
		   overlap --                                       -------------- */

	double h = 0.5 * dtau, S_mid[4], dS_mid[4], I_mid[4];

	for (int n = 0; n < 4; n++) {
		S_mid[n]  = 0.5 * (S_uw[n] + S_loc[n]) + 0.125 * dtau * (dS_uw[n] - dS_loc[n]);
		dS_mid[n] = 1.5 * (S_loc[n] - S_uw[n]) / dtau - 0.25 * (dS_uw[n] + dS_loc[n]);
	}
	/* --- With the same K' at both ends one exponential step is the
		   same thing --                                   -------------- */

	if (!memcmp(K_uw, K_loc, 16 * sizeof(double))) {
		exp_step_const(dtau, K_uw, S_uw, S_loc, dS_uw, dS_loc, I_uw, I_loc, EXP_ANY);
		return;
	}

	/* --- Where the second half is thick (its e^A is dropped, the
		   condition of exp_step_const), I_mid does not matter --  ---- */

	double eta2 = SQ(K_loc[0][1]) + SQ(K_loc[0][2]) + SQ(K_loc[0][3]);

	if (eta2 < 1.0 && h * (1.0 - sqrt(eta2)) > 40.0)
		memset(I_mid, 0, sizeof(I_mid));
	else
		exp_step_const(h, K_uw, S_uw, S_mid, dS_uw, dS_mid, I_uw, I_mid, EXP_ANY);
	exp_step_const(h, K_loc, S_mid, S_loc, dS_mid, dS_loc, I_mid, I_loc, EXP_ANY);
}
/* ------- end ---------------------------- Stokes_Exp_step.c ------- */

/* ------- begin -------------------------- Bezier3_Stokes_step.c --- */

double Bezier3_Stokes_step(double dtau_uw, double K_uw[4][4],
						   double K_loc[4][4], double dK_uw[4][4],
						   double dK_loc[4][4], double S_uw[4],
						   double S_loc[4], double dS_uw[4],
						   double dS_loc[4], double I_uw[4],
						   double I_loc[4])
{
	/* --- One step of the cubic DELO-Bezier solver for polarized light
		   (de la Cruz Rodriguez & Piskunov 2013), from the upwind to the
		   local point over the optical path dtau_uw: the Stokes vector
		   I_loc at the local point from I_uw at the upwind point, given
		   the reduced propagation matrix K' (StokesK), the source vector
		   S and their derivatives in tau at both points. Returns the
		   diagonal of the approximate operator, beta + theta.

		   Along the ray dI/dtau = S_eff - I, with S_eff = S - K' I, so
		   S_eff' = S' - dK' I - K' S + (K'^2 + K') I. The control points
		   are S_eff(u) + dt/3 S_eff'(u) (weight gamma) and
		   S_eff(0) - dt/3 S_eff'(0) (weight theta); the terms in I(0) go
		   to the left-hand side, M_loc.

		   The callers differ only in the derivatives they pass. At the
		   local point: centred and Steffen-limited (cent_deriv_vec,
		   cent_deriv_mat), or one-sided where the ray ends. At the
		   upwind point: in 1-D the one of the previous point (one-sided
		   in the first interval), in 2-D and 3-D, where the upwind point
		   is interpolated, the slope of the quadratic Bezier curve
		   (N-24). One function for all three since Q-25. Where
		   max |K'| dtau exceeds STOKES_XMAX_BEZIER3 and K' is physical
		   the step is the exponential one (Stokes_Exp_step,
		   use_exp_step, B-55). I_uw and I_loc
		   must not overlap: I_loc is built while I_uw is read
		   --                                              -------------- */

	static const double ident[4][4] =
		{{1.0, 0.0, 0.0, 0.0},
		 {0.0, 1.0, 0.0, 0.0},
		 {0.0, 0.0, 1.0, 0.0},
		 {0.0, 0.0, 0.0, 1.0}};

	register int i, j;
	double       eps = 0.0, alpha = 0.0, beta = 0.0, gamma = 0.0, theta = 0.0,
		   dt03;
	double A[4][4], Ma[4][4], Mb[4][4], Mc[4][4], M_loc[4][4];

	/* --- Where K' is large over the step the DELO-Bezier control
		   points amplify the polarization (B-55): the exponential step
		   instead --                                     -------------- */

	/* --- The approximate operator stays the scalar beta + theta of the
		   DELO-Bezier step, which ignores K', as everywhere else. The
		   exponential step's own coefficient of S_loc in I is the
		   polarized one, (K^-1)[0][0] in thick layers, above 1 with
		   dichroism (1.19 measured), and with it the first population
		   update of the 2-D test case inverted Ca II 866 nm (B-55) -- */

	Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);

	if (use_exp_step(dtau_uw, K_uw, K_loc, STOKES_XMAX_BEZIER3)) {
		Stokes_Exp_step(dtau_uw, K_uw, K_loc, S_uw, S_loc, dS_uw, dS_loc, I_uw,
						I_loc);
		return beta + theta;
	}

	m4m(K_uw, K_uw, Ma);  // Ku # Ku
	m4m(K_loc, K_loc, A); // K0 # K0

	dt03 = dtau_uw * 0.3333333333333333333333333;
	for (j = 0; j < 4; j++) {
		for (i = 0; i < 4; i++) {
			M_loc[j][i] = ident[j][i] + beta * K_loc[j][i] + theta * (dt03 * (A[j][i] - dK_loc[j][i] + K_loc[j][i]) + K_loc[j][i]);

			Ma[j][i] = eps * ident[j][i] - alpha * K_uw[j][i] + gamma * (dt03 * (Ma[j][i] - dK_uw[j][i] + K_uw[j][i]) - K_uw[j][i]);

			Mb[j][i] = alpha * ident[j][i] + gamma * (ident[j][i] -
													  dt03 * K_uw[j][i]);
			Mc[j][i] = beta * ident[j][i] + theta * (ident[j][i] +
													 dt03 * K_loc[j][i]);
		}
	}

	/* --- Right-hand side Ma I_uw + Mb S_uw + Mc S_loc +
		   (gamma dS_uw - theta dS_loc) dtau / 3 --       -------------- */

	for (i = 0; i < 4; i++) {
		I_loc[i] = 0.0;
		for (j = 0; j < 4; j++) {
			I_loc[i] += Ma[i][j] * I_uw[j] + Mb[i][j] * S_uw[j] +
						Mc[i][j] * S_loc[j];
		}
		I_loc[i] += dt03 * (gamma * dS_uw[i] - theta * dS_loc[i]);
	}
	solveLinearFast(M_loc, I_loc);

	return beta + theta;
}
/* ------- end ---------------------------- Bezier3_Stokes_step.c --- */

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

/* ------- begin -------------------------- longchar_dtau.c --------- */

static void longchar_dtau(double ds_uw, double ds_dw, double chi_uw,
						  double chi_loc, double chi_dw, double *dtau_uw,
						  double *dtau_dw)
{
	/* --- Optical paths of the sections upwind and downwind of the local
		   point of a long characteristic, chi integrated with the
		   Bezier curve of the short characteristics: the Steffen-limited
		   centred derivative at the local point, the slope of the
		   quadratic Bezier curve (N-24) at the upwind point and the
		   secant at the downwind point --              -------------- */

	double dchi_c, dchi_uw, dchi_dw, ds03, c1, c2;

	dchi_c  = cent_deriv(ds_uw, ds_dw, chi_uw, chi_loc, chi_dw);
	dchi_uw = 2.0 * (chi_loc - chi_uw) / ds_uw - dchi_c;
	dchi_dw = (chi_dw - chi_loc) / ds_dw;

	ds03     = ds_uw * 0.333333333333333333333333333;
	c1       = chi_loc - ds03 * dchi_c;
	c2       = chi_uw + ds03 * dchi_uw;
	*dtau_uw = ds_uw * (chi_loc + chi_uw + c1 + c2) * 0.25;

	ds03     = ds_dw * 0.333333333333333333333333333;
	c1       = chi_loc + ds03 * dchi_c;
	c2       = chi_dw - ds03 * dchi_dw;
	*dtau_dw = ds_dw * (chi_loc + chi_dw + c1 + c2) * 0.25;
}
/* ------- end ---------------------------- longchar_dtau.c --------- */

/* ------- begin -------------------------- LongChar_Bezier3.c ------ */

double LongChar_Bezier3(double ds_uw, double ds_dw, double chi_uw,
						double chi_loc, double chi_dw, double S_uw,
						double S_loc, double S_dw, double I_uw)
{
	/* --- One section of a long characteristic in 2-D or 3-D (N-39):
		   the intensity at the local (downwind) end of the section from
		   the upwind point, given chi and S at the upwind, local and
		   next downwind points along the characteristic and the
		   lengths ds_uw (upwind to local) and ds_dw (local to downwind).

		   The cubic Bezier step of the short characteristics
		   (Piecewise_Bezier3_2D and _3D): Steffen-limited centred
		   derivatives at the local point, the upwind slope of the
		   quadratic Bezier curve through the upwind and local points
		   (N-24), chi integrated with the same Bezier curve. The control
		   points lie between the values at the two ends, so that I stays
		   between the smallest and the largest of I_uw, S_uw, S_loc and
		   S_dw, however thin or thick the sections are. S_dw enters only
		   through the Steffen derivative at the local end, and I can
		   decrease when it grows: in the last section S_dw is the local
		   point, whose weight is then small and negative, down to about
		   -3 where S_loc ~ S_dw and the downwind section is thin
		   (test_longchar_3d). LongChar_Linear does not depend on S_dw
		   --                                              -------------- */

	double c1, c2, dtau_uw, dtau_dw, dS_c, dS_uw, dt03;
	double alpha = 0.0, beta = 0.0, gamma = 0.0, theta = 0.0, eps = 0.0;

	longchar_dtau(ds_uw, ds_dw, chi_uw, chi_loc, chi_dw, &dtau_uw, &dtau_dw);

	dS_c  = cent_deriv(dtau_uw, dtau_dw, S_uw, S_loc, S_dw);
	dt03  = dtau_uw * 0.333333333333333333333333333;
	dS_uw = 2.0 * (S_loc - S_uw) / dtau_uw - dS_c;
	c1    = S_loc - dt03 * dS_c;
	c2    = S_uw + dt03 * dS_uw;

	Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);

	return I_uw * eps + beta * S_loc + alpha * S_uw + theta * c1 + gamma * c2;
}
/* ------- end ---------------------------- LongChar_Bezier3.c ------ */

/* ------- begin -------------------------- LongChar_Linear.c ------- */

double LongChar_Linear(double ds_uw, double chi_uw, double chi_loc,
					   double S_uw, double S_loc, double I_uw)
{
	/* --- One section of a long characteristic with S linear in tau
		   and the trapezoidal optical path, as Piecewise_Linear_2D and
		   _3D (N-39) --                                  -------------- */

	double dtau_uw, w[2];

	dtau_uw = 0.5 * (chi_uw + chi_loc) * ds_uw;
	w2(dtau_uw, w);

	return (1.0 - w[0]) * I_uw + w[0] * S_loc + w[1] * (S_uw - S_loc) / dtau_uw;
}
/* ------- end ---------------------------- LongChar_Linear.c ------- */

/* ------- begin -------------------------- LongChar_Stokes_Bezier3.c */

void LongChar_Stokes_Bezier3(double ds_uw, double ds_dw, double chi_uw,
							 double chi_loc, double chi_dw,
							 double K_uw[4][4], double K_loc[4][4],
							 double K_dw[4][4], double S_uw[4],
							 double S_loc[4], double S_dw[4],
							 double I_uw[4], double I_loc[4])
{
	/* --- One section of a polarized long characteristic in 2-D or 3-D
		   (N-41): the Stokes vector I_loc at the local (downwind) end of
		   the section from I_uw at the upwind point, given chi, the
		   reduced propagation matrix K' (StokesK) and the source vector
		   S at the upwind, local and next downwind points and the
		   lengths ds_uw and ds_dw, as LongChar_Bezier3: the DELO-Bezier
		   step of the short characteristics (Bezier3_Stokes_step) with
		   the Steffen-limited centred derivatives at the local point and
		   the slopes of the quadratic Bezier curves (N-24) at the
		   upwind point. With K' = 0 it is LongChar_Bezier3 for I, and
		   Q = U = V stay 0 for an unpolarized S. I_uw and I_loc must
		   not overlap --                                 -------------- */

	register int n, m;
	double       dtau_uw, dtau_dw, dS_c[4], dS_uw[4], dK_c[4][4],
		dK_uw[4][4];

	longchar_dtau(ds_uw, ds_dw, chi_uw, chi_loc, chi_dw, &dtau_uw, &dtau_dw);

	cent_deriv_vec(dS_c, dtau_uw, dtau_dw, S_uw, S_loc, S_dw);
	cent_deriv_mat(dK_c, dtau_uw, dtau_dw, K_uw, K_loc, K_dw);
	for (n = 0; n < 4; n++) {
		dS_uw[n] = 2.0 * (S_loc[n] - S_uw[n]) / dtau_uw - dS_c[n];
		for (m = 0; m < 4; m++)
			dK_uw[n][m] = 2.0 * (K_loc[n][m] - K_uw[n][m]) / dtau_uw - dK_c[n][m];
	}
	Bezier3_Stokes_step(dtau_uw, K_uw, K_loc, dK_uw, dK_c, S_uw, S_loc,
						dS_uw, dS_c, I_uw, I_loc);
}
/* ------- end ---------------------------- LongChar_Stokes_Bezier3.c */

/* ------- begin -------------------------- LongChar_Stokes_Linear.c - */

void LongChar_Stokes_Linear(double ds_uw, double chi_uw, double chi_loc,
							double K_uw[4][4], double K_loc[4][4],
							double S_uw[4], double S_loc[4],
							double I_uw[4], double I_loc[4])
{
	/* --- One section of a polarized long characteristic with the
		   DELO method (Rees, Murphy & Durrant 1989, ApJ 339, 1093):
		   along the section dI/dtau = S_eff - I with S_eff = S - K' I,
		   here linear in tau, and the trapezoidal optical path, as
		   LongChar_Linear (N-41):

			 (1 + (w0 - w1 / dtau) K'_loc) I_loc =
			   ((1 - w0) - w1 / dtau K'_uw) I_uw + w0 S_loc
				 + w1 (S_uw - S_loc) / dtau.

		   With K' = 0 it is LongChar_Linear for I. Where max |K'| dtau
		   exceeds STOKES_XMAX_LINEAR and K' is physical the step is the
		   exponential one instead (Stokes_Exp_step, use_exp_step,
		   B-55). I_uw and I_loc must not
		   overlap --                                     -------------- */

	register int n, m;
	double       dtau_uw, w[2], c_uw, M_loc[4][4];

	dtau_uw = 0.5 * (chi_uw + chi_loc) * ds_uw;
	if (use_exp_step(dtau_uw, K_uw, K_loc, STOKES_XMAX_LINEAR)) {
		double sec[4];

		for (n = 0; n < 4; n++) sec[n] = (S_loc[n] - S_uw[n]) / dtau_uw;
		Stokes_Exp_step(dtau_uw, K_uw, K_loc, S_uw, S_loc, sec, sec, I_uw, I_loc);
		return;
	}
	w2(dtau_uw, w);
	c_uw = w[1] / dtau_uw;

	for (n = 0; n < 4; n++) {
		I_loc[n] = w[0] * S_loc[n] + c_uw * (S_uw[n] - S_loc[n]);
		for (m = 0; m < 4; m++) {
			M_loc[n][m] = (w[0] - c_uw) * K_loc[n][m];
			I_loc[n] -= c_uw * K_uw[n][m] * I_uw[m];
		}
		M_loc[n][n] += 1.0;
		I_loc[n] += (1.0 - w[0]) * I_uw[n];
	}
	solveLinearFast(M_loc, I_loc);
}
/* ------- end ---------------------------- LongChar_Stokes_Linear.c - */
