/* ------- file: -------------------------- test_stokes_step.c ------

	   One step of the polarized transfer equation dI/dtau = S - K I,
	   K = 1 + K' (bezier_aux.c): the exponential step Stokes_Exp_step,
	   and the DELO steps Bezier3_Stokes_step and LongChar_Stokes_Linear,
	   which fall back to it where max |K'| dtau is large (B-55).

	   -- Stokes_Exp_step against closed forms: with K' = 0 the scalar
		  solution for a cubic S (e^-h I + h sum a_m m! phi_(m+1)(-h),
		  phi from their series or closed forms), Q = U = V = 0; a pure
		  magneto-optical rho_V turns (Q, U) by rho_V h and damps them
		  by e^-h, for rho_V h up to 1e4;
	   -- local equilibrium: in a thick step with a line-only source
		  S_loc = S0 (1 + K'_loc) e0, I_loc -> S0 e0 whatever K'_uw is;
	   -- exact for constant K' and cubic S: one step equals two half
		  steps, for random physical K' (|eta| < 0.9, |rho| up to 10)
		  and h from 1e-2 to 1e2;
	   -- no amplification: with S = 0 the norm of the Stokes vector
		  cannot grow for a physical K' (|eta| <= 1). Over random
		  physical K' at either end (|eta| < 0.95, |rho| up to 30),
		  h from 1e-3 to 1e2 and both derivative choices of the 2-D and
		  3-D solvers, none of the three steps raises it. Before B-55
		  the Bezier step did by up to 1.4e4, DELO-linear by up to 8,
		  where K' fell steeply over an optically long step; the
		  2-D and 3-D FULL_STOKES runs diverged on that.

	   -- the cone: the Bezier control points of S are kept in it
		  (N-47), so with S and I_uw in the cone I_loc is too, with
		  K' = 0 and in the exponential step.

	   -v prints the measured values.
	   --                                              -------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "rh.h"
#include "bezier.h"
#include "check.h"

/* --- One step against 64 substeps with K' changing by 5 % across a
	   step of h = 1 .. 300 (rho up to 10): measured 1.6e-2, the
	   first-order error of the two-half-step approximation, largest
	   where the Faraday rotation rho h is large --         -------------- */
#define VAR_BOUND 3.0e-2

static unsigned seed = 2026;

static double urand(void)
{
	seed = seed * 1103515245u + 12345u;
	return ((seed >> 8) & 0xffffff) / 16777215.0;
}

/* --- Random physical K': |eta| = eta in a random direction, rho
	   components uniform in [-rho, rho] --              -------------- */

static void random_K(double K[4][4], double eta, double rho)
{
	double th = acos(2.0 * urand() - 1.0), ph = 2.0 * M_PI * urand(), e[3], r[3];

	e[0] = eta * sin(th) * cos(ph);
	e[1] = eta * sin(th) * sin(ph);
	e[2] = eta * cos(th);
	for (int n = 0; n < 3; n++) r[n] = (2.0 * urand() - 1.0) * rho;

	double k[4][4] = {{0.0, e[0], e[1], e[2]},
					  {e[0], 0.0, r[2], -r[1]},
					  {e[1], -r[2], 0.0, r[0]},
					  {e[2], r[1], -r[0], 0.0}};
	memcpy(K, k, sizeof(k));
}

static double norm4(const double v[4])
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
}

/* --- phi_k(-h) for scalar h: series below 1, the recurrence above -- */

static double phi_scalar(int k, double h)
{
	if (h < 1.0) {
		double s = 0.0, t = 1.0;

		for (int j = 1; j <= k; j++) t /= j;
		for (int n = 0; n < 40; n++) {
			s += t;
			t *= -h / (n + k + 1);
		}
		return s;
	}
	double p = exp(-h); /* phi_0 */
	for (int j = 1; j <= k; j++) {
		double fact = 1.0;
		for (int i = 2; i < j; i++) fact *= i;
		p = (1.0 / fact - p) / h; /* phi_j = (1/(j-1)! - phi_(j-1)) / h */
	}
	return p;
}

static void check_closed_forms(bool_t verbose)
{
	const double hs[] = {1.0e-6, 1.0e-3, 0.1, 1.0, 10.0, 100.0};
	double       err0 = 0.0, quv = 0.0, errrot = 0.0, thick = 0.0;

	for (int i = 0; i < 6; i++) {
		double h = hs[i], Z[4][4] = {{0.0}}, Su[4] = {0.3, 0, 0, 0}, Sl[4] = {0.7, 0, 0, 0},
			   du[4] = {0.2, 0, 0, 0}, dl[4] = {-0.1, 0, 0, 0}, Iu[4] = {0.5, 0, 0, 0}, Il[4];
		double m0 = h * du[0], m1 = h * dl[0], a[4];

		a[0]        = Su[0];
		a[1]        = m0;
		a[2]        = 2.0 * (3.0 * (Sl[0] - Su[0]) - 2.0 * m0 - m1);
		a[3]        = 6.0 * (2.0 * (Su[0] - Sl[0]) + m0 + m1);
		double want = exp(-h) * Iu[0];
		for (int m = 0; m < 4; m++) want += h * a[m] * phi_scalar(m + 1, h);
		Stokes_Exp_step(h, Z, Z, Su, Sl, du, dl, Iu, Il);
		err0 = MAX(err0, fabs(Il[0] - want) / fabs(want));
		quv  = MAX(quv, MAX(fabs(Il[1]), MAX(fabs(Il[2]), fabs(Il[3]))));

		for (double rh = 1.0; rh <= 1.0e4; rh *= 10.0) {
			double K[4][4] = {{0.0}}, S[4] = {0.0}, J[4] = {1.0, 1.0, 0.0, 0.0}, O[4];

			K[1][2] = rh / h;
			K[2][1] = -rh / h;
			Stokes_Exp_step(h, K, K, S, S, S, S, J, O);
			if (h > 40.0) {
				/* --- Thick-step branch: e^A dropped, below 4e-18 -- */
				thick = MAX(thick, norm4(O) / norm4(J));
				continue;
			}
			double q = O[1] / exp(-h), u = O[2] / exp(-h);
			errrot = MAX(errrot, MIN(hypot(q - cos(rh), u - sin(rh)),
									 hypot(q - cos(rh), u + sin(rh))) /
									 MAX(1.0, rh * 1.0e-4));
		}
	}
	if (verbose)
		printf("Exp step, K' = 0 cubic S: rel error %.2e, QUV %.1e; rho_V rotation: "
			   "error %.2e (per 1e4 of rho h), thick |I_loc| / |I_uw| %.1e\n",
			   err0, quv, errrot, thick);
	CHECK(err0 <= 1.0e-13 && quv == 0.0, "Exp step, K' = 0: error %.3g, QUV %.3g", err0, quv);
	CHECK(errrot <= 1.0e-11, "Exp step, rho_V rotation: error %.3g", errrot);
	CHECK(thick <= 4.0e-18, "Exp step, thick: |I_loc| / |I_uw| = %.3g with S = 0", thick);
}

static void check_cubic_exact(bool_t verbose)
{
	double worst = 0.0;

	for (int t = 0; t < 5000; t++) {
		double h     = pow(10.0, 4.0 * urand() - 2.0), K[4][4], c[4][4], S[3][4], d[3][4],
			   Iu[4] = {1.0, 0.1, -0.2, 0.05}, A[4], B[4], C[4], scale = 0.0;

		random_K(K, 0.9 * urand(), pow(10.0, 2.0 * urand() - 1.0));
		for (int n = 0; n < 4; n++)
			for (int k = 0; k < 4; k++) c[n][k] = (urand() - 0.5) / pow(h, k);
		for (int p = 0; p < 3; p++) {
			double tt = 0.5 * h * p;
			for (int n = 0; n < 4; n++) {
				S[p][n] = c[n][0] + tt * (c[n][1] + tt * (c[n][2] + tt * c[n][3]));
				d[p][n] = c[n][1] + tt * (2.0 * c[n][2] + 3.0 * tt * c[n][3]);
			}
		}
		Stokes_Exp_step(h, K, K, S[0], S[2], d[0], d[2], Iu, A);
		Stokes_Exp_step(0.5 * h, K, K, S[0], S[1], d[0], d[1], Iu, B);
		Stokes_Exp_step(0.5 * h, K, K, S[1], S[2], d[1], d[2], B, C);
		scale = norm4(C);
		for (int n = 0; n < 4; n++) worst = MAX(worst, fabs(A[n] - C[n]) / scale);
	}
	if (verbose) printf("Exp step, cubic S and constant K': one step vs two halves %.2e\n", worst);
	CHECK(worst <= 1.0e-12, "Exp step, cubic S: one step and two halves differ by %.3g", worst);
}

/* --- Local equilibrium in a thick step (B-55): with a line-only
	   source S_loc = S0 (1 + K'_loc) e0 the solution tends to S0 e0
	   whatever K' upwind is. A step with the mean K' over the step
	   left I 5 % low where the local dichroism exceeded the upwind
	   one; the exponential fallback uses K'_uw on the first half and
	   K'_loc on the second, and the Bezier step gets it right
	   itself --                                         -------------- */

static void check_thick_limit(bool_t verbose)
{
	double worst = 0.0;

	for (int t = 0; t < 2000; t++) {
		double h = 1.0e4, Ku[4][4], Kl[4][4], Su[4], Sl[4], Iu[4] = {1.0, 0.0, 0.0, 0.0},
			   Il[4], dKu[4][4], dKc[4][4], Z[4] = {0.0};

		random_K(Ku, 0.9 * urand(), pow(10.0, 2.0 * urand() - 1.0));
		random_K(Kl, 0.9 * urand(), pow(10.0, 2.0 * urand() - 1.0));
		for (int i = 0; i < 4; i++) {
			Su[i] = (i == 0) + Ku[i][0];
			Sl[i] = (i == 0) + Kl[i][0];
			for (int j = 0; j < 4; j++) {
				dKu[i][j] = dKc[i][j] = (Kl[i][j] - Ku[i][j]) / h;
			}
		}
		Bezier3_Stokes_step(h, Ku, Kl, dKu, dKc, Su, Sl, Z, Z, Iu, Il);
		worst = MAX(worst, fabs(Il[0] - 1.0));
		for (int n = 1; n < 4; n++) worst = MAX(worst, fabs(Il[n]));
	}
	if (verbose) printf("thick step (dtau 1e4), line-only source: |I - S0 e0| <= %.2e\n", worst);
	CHECK(worst <= 1.0e-4, "thick step: I differs from the local equilibrium by %.3g", worst);
}

/* --- Every thickness (code review): one step against 64 substeps.
	   With constant K' and cubic S the substeps (thin, series branch)
	   are exact, so one step over h = 1 .. 300, which crosses into the
	   thick branch, must agree to rounding. With K' changing linearly
	   from K'_uw to K'_loc the 64 substeps (K' at their own ends) are
	   the accurate reference for the two-half-step approximation;
	   its error is first order in the change of K' (measured, printed
	   with -v) --                                       -------------- */

static void substeps(double h, double Ku[4][4], double Kl[4][4], double c[4][4],
					 double I[4], int N)
{
	for (int s = 0; s < N; s++) {
		double t0 = h * s / N, t1 = h * (s + 1) / N, K0[4][4], K1[4][4], S0[4], S1[4], d0[4],
			   d1[4], In[4];

		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++) {
				K0[i][j] = Ku[i][j] + (Kl[i][j] - Ku[i][j]) * t0 / h;
				K1[i][j] = Ku[i][j] + (Kl[i][j] - Ku[i][j]) * t1 / h;
			}
		for (int n = 0; n < 4; n++) {
			S0[n] = c[n][0] + t0 * (c[n][1] + t0 * (c[n][2] + t0 * c[n][3]));
			S1[n] = c[n][0] + t1 * (c[n][1] + t1 * (c[n][2] + t1 * c[n][3]));
			d0[n] = c[n][1] + t0 * (2.0 * c[n][2] + 3.0 * t0 * c[n][3]);
			d1[n] = c[n][1] + t1 * (2.0 * c[n][2] + 3.0 * t1 * c[n][3]);
		}
		Stokes_Exp_step(h / N, K0, K1, S0, S1, d0, d1, I, In);
		memcpy(I, In, sizeof(In));
	}
}

static void check_thickness_range(bool_t verbose)
{
	double worst_const = 0.0, worst_var = 0.0;

	for (int t = 0; t < 400; t++) {
		double h = pow(10.0, 2.48 * urand()), Ku[4][4], Kl[4][4], c[4][4], A[4] = {1, 0.1, -0.05, 0.02},
			   B[4], C[4];

		random_K(Ku, 0.8 * urand(), pow(10.0, 2.0 * urand() - 1.0));
		for (int n = 0; n < 4; n++)
			for (int k = 0; k < 4; k++) c[n][k] = ((n == 0) + 0.5 * (urand() - 0.5)) / pow(h, k);
		memcpy(B, A, sizeof(A));
		substeps(h, Ku, Ku, c, B, 1);
		memcpy(C, A, sizeof(A));
		substeps(h, Ku, Ku, c, C, 64);
		for (int n = 0; n < 4; n++) worst_const = MAX(worst_const, fabs(B[n] - C[n]) / norm4(C));

		random_K(Kl, 0.8 * urand(), pow(10.0, 2.0 * urand() - 1.0));
		for (int i = 0; i < 4; i++)
			for (int j = 0; j < 4; j++) Kl[i][j] = Ku[i][j] + 0.05 * (Kl[i][j] - Ku[i][j]);
		memcpy(B, A, sizeof(A));
		substeps(h, Ku, Kl, c, B, 1);
		memcpy(C, A, sizeof(A));
		substeps(h, Ku, Kl, c, C, 64);
		for (int n = 0; n < 4; n++) worst_var = MAX(worst_var, fabs(B[n] - C[n]) / norm4(C));
	}
	if (verbose)
		printf("h 1 .. 300, one step vs 64 substeps: constant K' %.2e, K' changing by 5 %% %.2e\n",
			   worst_const, worst_var);
	CHECK(worst_const <= 1.0e-12, "constant K': one step and 64 substeps differ by %.3g", worst_const);
	CHECK(worst_var <= VAR_BOUND, "K' changing: one step and 64 substeps differ by %.3g", worst_var);
}

static void check_no_amplification(bool_t verbose)
{
	double      worst[4] = {0.0, 0.0, 0.0, 0.0};
	const char *name[4]  = {"Bezier3 (dK_c = 0, 2 x secant upwind)", "Bezier3 (secant)",
							"DELO-linear", "Exp"};

	for (int t = 0; t < 100000; t++) {
		double h = pow(10.0, 5.0 * urand() - 3.0), Ku[4][4], Kl[4][4], Z[4] = {0.0},
			   p = urand(), th = acos(2.0 * urand() - 1.0), ph = 2.0 * M_PI * urand(), Iu[4],
			   Il[4];

		random_K(Ku, 0.95 * urand(), pow(10.0, 3.0 * urand() - 1.5));
		random_K(Kl, 0.95 * urand(), pow(10.0, 3.0 * urand() - 1.5));
		if (t % 2) {
			/* --- Neighbouring K' close to each other too -- --- */
			for (int i = 0; i < 4; i++)
				for (int j = 0; j < 4; j++) Kl[i][j] = Ku[i][j] * (1.0 + 0.01 * (urand() - 0.5));
		}
		Iu[0] = 1.0;
		Iu[1] = p * sin(th) * cos(ph);
		Iu[2] = p * sin(th) * sin(ph);
		Iu[3] = p * cos(th);

		for (int s = 0; s < 4; s++) {
			if (s < 2) {
				double dKu[4][4], dKc[4][4];
				for (int i = 0; i < 4; i++)
					for (int j = 0; j < 4; j++) {
						double sec = (Kl[i][j] - Ku[i][j]) / h;
						dKu[i][j]  = s ? sec : 2.0 * sec;
						dKc[i][j]  = s ? sec : 0.0;
					}
				Bezier3_Stokes_step(h, Ku, Kl, dKu, dKc, Z, Z, Z, Z, Iu, Il);
			} else if (s == 2)
				LongChar_Stokes_Linear(h, 1.0, 1.0, Ku, Kl, Z, Z, Iu, Il);
			else
				Stokes_Exp_step(h, Ku, Kl, Z, Z, Z, Z, Iu, Il);
			double a = norm4(Il) / norm4(Iu);
			worst[s] = MAX(worst[s], (a == a) ? a : HUGE_VAL);
		}
	}
	for (int s = 0; s < 4; s++) {
		if (verbose) printf("S = 0, largest |I_loc| / |I_uw|: %-38s %.6f\n", name[s], worst[s]);
		CHECK(worst[s] <= 1.0, "%s amplifies the Stokes vector by %.4g", name[s], worst[s]);
	}
}

/* --- The closed form of the constant-K' step against the Taylor
	   series and squarings (Stokes_Exp_const_closed and
	   Stokes_Exp_const_taylor). B = -dtau K' has eigenvalues +-a,
	   +-i b; with eta parallel to rho, a = dtau |eta| and
	   b = dtau |rho|. Cases: a grid of a from 1e-6 to 1e2 and b from
	   1e-6 to 1e4; dtau around 3, where exp_pair changes from the power
	   series to the recurrence, with small a and b; random K'; K' with
	   |eta| = |rho| at an angle whose cosine is the degeneracy ratio
	   (a^2 + b^2) / (dtau^2 (eta^2 + rho^2)), just above and below
	   the cut of 0.1. The closed form must be taken wherever the ratio
	   is above the cut and not below it, and agree with the Taylor
	   path within a bound times max(10, b), since both are limited by
	   the rounding of the phase b where b is large. Measured (gcc,
	   x86-64), bounds 3 to 5 times above: grid 4.2e-15, dtau ~ 3 2e-16
	   (2.4e-16 with the series of exp_pair also for k = 0, which
	   lost 1e2 ulp in e^z cosh r near z = -3: too little in a whole
	   step to be told from the Taylor path's own error), near the cut
	   3e-15, random 1.4e-14. The random K' include strongly dichroic steps with
	   dtau of order 50 to 70 (|eta| 0.8), where e^-dtau(1 - a/dtau)
	   and the polynomial terms of phi_k nearly cancel: there both
	   paths are 1e-13 off long double (closed 1.2e-13, Taylor
	   1.5e-13), a limit of the problem, not of either method -- --- */

/* --- The Bezier control points of S stay in the cone (N-47):
	   source vectors S_uw, S_loc, S_dw and I_uw in the cone, the
	   derivatives of the 2-D and 3-D solvers (cent_deriv_vec at the
	   local point, 2 x secant minus it upwind, N-24) and of 1-D (the
	   secant upwind), K' = 0 (the DELO step a combination of nodes and
	   control points with weights >= 0) and a random physical K' in a
	   step thick enough for the exponential one: I_loc stays in the
	   cone. First the case of the interpolation audit (S_I 30 times
	   higher upwind, an extremum of S_I at the local point while S_Q
	   falls monotonically, |Q| / I = 0.9 at the nodes): |Q| / I of
	   I_loc was 2.08 before N-47; then random nodes with |pol| / I up
	   to 0.99 (up to 6 before). A linear S (no extremum) is not
	   changed (check_cubic_exact, test_formal_stokes) -- ---------- */

static double cone_ratio(const double v[4])
{
	double p = sqrt(SQ(v[1]) + SQ(v[2]) + SQ(v[3]));

	return (v[0] > 0.0) ? p / v[0] : ((p > 0.0 || v[0] < 0.0) ? HUGE_VAL : 0.0);
}

static void random_in_cone(double v[4], double I, double pmax)
{
	double p = pmax * urand(), th = acos(2.0 * urand() - 1.0), ph = 2.0 * M_PI * urand();

	v[0] = I;
	v[1] = I * p * sin(th) * cos(ph);
	v[2] = I * p * sin(th) * sin(ph);
	v[3] = I * p * cos(th);
}

static double cone_step(double dtu, double dtd, double Ku[4][4], double Kl[4][4],
						double Su[4], double Sl[4], double Sd[4], double Iu[4], int upwind)
{
	double dSc[4], dSu[4], Il[4], dKu[4][4] = {{0.0}}, dKc[4][4] = {{0.0}};

	cent_deriv_vec(dSc, dtu, dtd, Su, Sl, Sd);
	for (int n = 0; n < 4; n++) {
		double sec = (Sl[n] - Su[n]) / dtu;
		dSu[n]     = upwind ? 2.0 * sec - dSc[n] : sec;
	}
	Bezier3_Stokes_step(dtu, Ku, Kl, dKu, dKc, Su, Sl, dSu, dSc, Iu, Il);
	return cone_ratio(Il);
}

static void check_cone(bool_t verbose)
{
	double Z[4][4] = {{0.0}}, worst0 = 0.0, worstK = 0.0;

	/* --- The case of the audit --                    -------------- */

	double Su[4] = {30.0, 27.0, 0.0, 0.0}, Sl[4] = {1.0, 0.9, 0.0, 0.0},
		   Sd[4] = {1.5, 0.0, 0.0, 0.0}, Iu[4] = {30.0, 27.0, 0.0, 0.0};
	double audit = cone_step(3.0, 0.01, Z, Z, Su, Sl, Sd, Iu, 1);

	for (int t = 0; t < 200000; t++) {
		double dtu = pow(10.0, 3.0 * urand() - 1.5), dtd = pow(10.0, 3.0 * urand() - 2.0),
			   Ku[4][4], Kl[4][4];

		random_in_cone(Su, pow(30.0, 2.0 * urand() - 1.0), 0.99);
		random_in_cone(Sl, 1.0, 0.99);
		random_in_cone(Sd, pow(30.0, 2.0 * urand() - 1.0), 0.99);
		random_in_cone(Iu, pow(30.0, 2.0 * urand() - 1.0), 0.99);
		worst0 = MAX(worst0, cone_step(dtu, dtd, Z, Z, Su, Sl, Sd, Iu, t % 2));

		/* --- A physical K' and a step thick enough for the
			   exponential step (max |K'| dtau > 1) -- ------------ */

		random_K(Ku, 0.95 * urand(), pow(10.0, 2.0 * urand() - 1.0));
		random_K(Kl, 0.95 * urand(), pow(10.0, 2.0 * urand() - 1.0));
		dtu    = 2.0 + 30.0 * urand();
		worstK = MAX(worstK, cone_step(dtu, dtd, Ku, Kl, Su, Sl, Sd, Iu, t % 2));
	}
	if (verbose)
		printf("cone: audit case |pol| / I %.6f; random, K' = 0: %.6f; exponential step: %.6f\n",
			   audit, worst0, worstK);
	CHECK(audit <= 1.0 && worst0 <= 1.0 && worstK <= 1.0,
		  "I_loc leaves the cone: audit case %.4g, random K' = 0 %.4g, exponential step %.4g",
		  audit, worst0, worstK);
}

#define CLOSED_BOUND_GRID   1.5e-14
#define CLOSED_BOUND_DTAU3  1.0e-15
#define CLOSED_BOUND_CUT    1.0e-14
#define CLOSED_BOUND_RANDOM 5.0e-14

static int const_step_diff(double dtau, double K[4][4], double *diff)
{
	double S0[4], S1[4], d0[4], d1[4], I0[4], Ic[4], It[4], scale = 0.0, err = 0.0;

	for (int i = 0; i < 4; i++) {
		double w = (i == 0) ? 1.0 : 0.3, o = (i == 0) ? 1.5 : 0.0;

		S0[i] = o + w * (2.0 * urand() - 1.0);
		S1[i] = o + w * (2.0 * urand() - 1.0);
		d0[i] = w * (2.0 * urand() - 1.0);
		d1[i] = w * (2.0 * urand() - 1.0);
		I0[i] = o + w * (2.0 * urand() - 1.0);
	}
	if (!Stokes_Exp_const_closed(dtau, K, S0, S1, d0, d1, I0, Ic)) return 0;
	Stokes_Exp_const_taylor(dtau, K, S0, S1, d0, d1, I0, It);
	for (int i = 0; i < 4; i++) {
		scale = fmax(scale, fmax(fabs(It[i]), fmax(fabs(I0[i]), fabs(S1[i]))));
		err   = fmax(err, fabs(Ic[i] - It[i]));
	}
	*diff = err / scale;
	return 1;
}

/* --- K' with eta = ue |eta|, rho = ur |rho| (unit vectors) --  ---- */

static void make_K(double K[4][4], const double ue[3], double eta, const double ur[3],
				   double rho)
{
	double e[3]    = {ue[0] * eta, ue[1] * eta, ue[2] * eta};
	double q[3]    = {ur[0] * rho, ur[1] * rho, ur[2] * rho};
	double k[4][4] = {{0.0, e[0], e[1], e[2]},
					  {e[0], 0.0, q[2], -q[1]},
					  {e[1], -q[2], 0.0, q[0]},
					  {e[2], q[1], -q[0], 0.0}};

	memcpy(K, k, sizeof(k));
}

static void random_unit(double u[3])
{
	double th = acos(2.0 * urand() - 1.0), ph = 2.0 * M_PI * urand();

	u[0] = sin(th) * cos(ph);
	u[1] = sin(th) * sin(ph);
	u[2] = cos(th);
}

static void check_closed_form(bool_t verbose)
{
	double worst[4] = {0.0}, d, K[4][4], u[3], w[3];
	int    nmiss[4] = {0}, nwrong_cut = 0;

	/* --- 0: (a, b) grid, 1: dtau around 3 --         -------------- */

	for (int la = -6; la <= 2; la++)
		for (int lb = -6; lb <= 4; lb++) {
			double a = pow(10.0, la), b = pow(10.0, lb), dtau = fmax(1.0, a / 0.9);

			random_unit(u);
			make_K(K, u, a / dtau, u, b / dtau);
			if (const_step_diff(dtau, K, &d))
				worst[0] = fmax(worst[0], d / fmax(10.0, b));
			else
				nmiss[0]++;
		}
	static const double dt3[] = {2.0, 2.5, 2.9, 2.999, 3.0, 3.001, 3.1, 4.0};
	static const double ab3[] = {1.0e-4, 1.0e-2, 0.3, 0.9};

	for (int t = 0; t < 8; t++)
		for (int ia = 0; ia < 4; ia++)
			for (int ib = 0; ib < 4; ib++) {
				double dtau = dt3[t];

				random_unit(u);
				make_K(K, u, ab3[ia] / dtau, u, ab3[ib] / dtau);
				if (const_step_diff(dtau, K, &d))
					worst[1] = fmax(worst[1], d / 10.0);
				else
					nmiss[1]++;
			}

	/* --- 2: random K' --                              -------------- */

	for (int n = 0; n < 20000; n++) {
		double rho = pow(10.0, -3.0 + 6.0 * urand()), dtau = pow(10.0, -3.0 + 5.0 * urand());

		random_K(K, 0.95 * urand(), rho);
		if (const_step_diff(dtau, K, &d))
			worst[2] = fmax(worst[2], d / fmax(10.0, 2.0 * rho * dtau));
		else
			nmiss[2]++;
	}

	/* --- 3: |eta| = |rho| at an angle with cosine = the ratio, just
		   above (closed form) and below (Taylor) the cut 0.1 --    -- */

	for (int n = 0; n < 2000; n++) {
		double m = 0.05 + 0.9 * urand(), dtau = pow(10.0, -1.0 + 2.0 * urand());
		double cosang = (n % 2) ? 0.1 * (1.0 + 1.0e-3 + urand()) : 0.1 * (1.0 - 1.0e-3) * urand();
		double t[3], sn = sqrt(1.0 - cosang * cosang);

		random_unit(u);
		random_unit(t); /* w = cosang u + sn (unit vector _|_ u) */
		double dot = t[0] * u[0] + t[1] * u[1] + t[2] * u[2];
		for (int i = 0; i < 3; i++) t[i] -= dot * u[i];
		double tn = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
		for (int i = 0; i < 3; i++) w[i] = cosang * u[i] + sn * t[i] / tn;
		make_K(K, u, m, w, m);
		int closed = const_step_diff(dtau, K, &d);

		if (closed != (n % 2)) nwrong_cut++;
		if (closed) worst[3] = fmax(worst[3], d / fmax(10.0, m * dtau));
	}
	if (verbose)
		printf("closed form vs Taylor, / max(10, b): grid %.3g, dtau ~ 3 %.3g, random %.3g,"
			   " near the cut %.3g; not taken %d, %d, %d; wrong side of the cut %d\n",
			   worst[0], worst[1], worst[2], worst[3], nmiss[0], nmiss[1], nmiss[2], nwrong_cut);
	CHECK(nmiss[0] == 0 && nmiss[1] == 0, "closed form not taken at %d grid points, %d near dtau 3",
		  nmiss[0], nmiss[1]);
	CHECK(nwrong_cut == 0, "closed form on the wrong side of the degeneracy cut %d times", nwrong_cut);
	CHECK(worst[0] <= CLOSED_BOUND_GRID, "closed form, (a, b) grid: %.3g max(10, b)", worst[0]);
	CHECK(worst[1] <= CLOSED_BOUND_DTAU3, "closed form, dtau ~ 3: %.3g max(10, b)", worst[1]);
	CHECK(worst[2] <= CLOSED_BOUND_RANDOM, "closed form, random K': %.3g max(10, b)", worst[2]);
	CHECK(worst[3] <= CLOSED_BOUND_CUT, "closed form, near the cut: %.3g max(10, b)", worst[3]);
}

int main(int argc, char *argv[])
{
	bool_t verbose = (argc > 1 && !strcmp(argv[1], "-v"));

	check_closed_forms(verbose);
	check_closed_form(verbose);
	check_cubic_exact(verbose);
	check_thick_limit(verbose);
	check_thickness_range(verbose);
	check_no_amplification(verbose);
	check_cone(verbose);
	return check_summary("test_stokes_step");
}
/* ------- end ---------------------------- test_stokes_step.c ------ */
