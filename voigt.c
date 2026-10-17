/* ------- file: -------------------------- voigt.c -----------------

	   Voigt function H(a, v) and Faraday-Voigt function F(a, v):

		 H + i F = w(v + ia),  w(z) = exp(-z^2) erfc(-iz),

	   the Faddeeva function, with damping parameter a >= 0 and
	   frequency v, both in Doppler widths. H(a, 0) = 1 for a = 0 (not
	   normalized to unit area; divide by sqrt(pi) for that).

	   A hybrid of two approximations, after Schreier (2018, JQSRT 213,
	   13), who combines Humlicek's asymptotic region with Weideman's
	   rational approximation:

	   -- |v| + a >= 8: the Laplace continued fraction
			w(z) = (i/sqrt(pi)) / (z - (1/2) / (z - 1 / (z - (3/2) / ...)))
		  to 6 levels (Humlicek's region I is its first level);
	   -- inside: Weideman's (1994, SIAM J. Numer. Anal. 31, 1497)
		  rational approximation with N = 32 terms;
	   -- a = 0: H = exp(-v^2) exactly (the approximations have an
		  absolute error, which no longer is small compared with the
		  Gaussian far from line center).

	   Against mpmath (tests/unit/test_voigt.c), for a from 0 to 100 and
	   v from 0 to 1e4: relative error of H below 1.6e-6 where H >=
	   1e-12 H(a, 0) (largest for a ~ 1e-6 just outside the Doppler
	   core), 2.7e-10 relative to H(a, 0), and F within 1e-9 of max |F|.
	   About 50 ns per call. This replaces the five algorithms of RH 2
	   (Armstrong, Rybicki, Hui et al., Humlicek, lookup table; errors
	   up to 1.2e-4), so that profiles with and without polarization are
	   the same function.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>

#include "rh.h"

#define INVSQRTPI 0.56418958354775628695

/* --- Asymptotic region and depth of the continued fraction: |v| + a
	   >= 8 with 6 levels is accurate to 3e-10 of H(a, 0) -- ---------- */

#define VOIGT_ASYMPTOTIC 8.0
#define VOIGT_CF_LEVELS  6

/* --- Weideman's coefficients for N = 32, highest power first, from
	   tests/unit/make_weideman_coefficients.py --   -------------- */

#define WEIDEMAN_N 32
#define WEIDEMAN_L 4.75682846001088426687

static const double weideman_a[WEIDEMAN_N] = {
	-1.3033481245090925e-12,
	3.7410343579709147e-12,
	8.0303941233429979e-12,
	-2.1543623534540259e-11,
	-5.544247434751948e-11,
	1.1658238505513269e-10,
	4.1537422837681178e-10,
	-5.2310221157958813e-10,
	-3.2080152206021089e-9,
	8.1248891740501905e-10,
	2.3797556704926332e-8,
	2.2930439030891306e-8,
	-1.481307891793793e-7,
	-4.1840763711855293e-7,
	4.2558331373556316e-7,
	4.4015317315303742e-6,
	6.821031944000645e-6,
	-2.1409619201818263e-5,
	-1.3075449254609854e-4,
	-2.4532980270018119e-4,
	3.9259136070078965e-4,
	4.5195411053492856e-3,
	1.9006155784845476e-2,
	5.7304403529837191e-2,
	1.4060716226893785e-1,
	2.9544451071508732e-1,
	5.4601397206393423e-1,
	9.0192548936479998e-1,
	1.3455441692345451,
	1.8256696296324813,
	2.2635372999002678,
	2.5722534081245693,
};

/* ------- begin -------------------------- Voigt.c ----------------- */

double Voigt(double a, double v, double *F)
{
	double x = fabs(v), H, Fx;

	if (x + a >= VOIGT_ASYMPTOTIC) {

		/* --- Continued fraction as the ratio A/B of its convergents:
			   A_k = z A_(k-1) - (k/2) A_(k-2), and the same for B, so
			   that there is a single division --          -------------- */

		double Ar = 1.0, Ai = 0.0, Br = x, Bi = a;
		double Apr = 0.0, Api = 0.0, Bpr = 1.0, Bpi = 0.0;

		for (int k = 1; k <= VOIGT_CF_LEVELS; k++) {
			double c   = 0.5 * k;
			double nAr = x * Ar - a * Ai - c * Apr, nAi = x * Ai + a * Ar - c * Api;
			double nBr = x * Br - a * Bi - c * Bpr, nBi = x * Bi + a * Br - c * Bpi;

			Apr = Ar;
			Api = Ai;
			Ar  = nAr;
			Ai  = nAi;
			Bpr = Br;
			Bpi = Bi;
			Br  = nBr;
			Bi  = nBi;
		}
		double q = 1.0 / (Br * Br + Bi * Bi);

		/* --- w = (i / sqrt(pi)) A / B --              -------------- */

		H  = INVSQRTPI * (Ar * Bi - Ai * Br) * q;
		Fx = INVSQRTPI * (Ar * Br + Ai * Bi) * q;
	} else {

		/* --- Weideman: with r = 1 / (L - iz) and Z = (L + iz) r,
			   w = (2 p(Z) r + 1/sqrt(pi)) r. p is summed as its odd and
			   even powers in Z^2, two independent Horner chains -- --- */

		double dr = WEIDEMAN_L + a, q = 1.0 / (dr * dr + x * x);
		double rr = dr * q, ri = x * q;
		double nr = WEIDEMAN_L - a;
		double Zr = nr * rr - x * ri, Zi = nr * ri + x * rr;
		double Z2r = Zr * Zr - Zi * Zi, Z2i = 2.0 * Zr * Zi;
		double odr = weideman_a[0], odi = 0.0, evr = weideman_a[1], evi = 0.0;

		for (int n = 2; n < WEIDEMAN_N; n += 2) {
			double t = odr * Z2r - odi * Z2i + weideman_a[n];

			odi = odr * Z2i + odi * Z2r;
			odr = t;
			t   = evr * Z2r - evi * Z2i + weideman_a[n + 1];
			evi = evr * Z2i + evi * Z2r;
			evr = t;
		}
		double pr = odr * Zr - odi * Zi + evr, pi = odr * Zi + odi * Zr + evi;
		double sr = 2.0 * (pr * rr - pi * ri) + INVSQRTPI,
			   si = 2.0 * (pr * ri + pi * rr);

		H  = sr * rr - si * ri;
		Fx = sr * ri + si * rr;
	}
	if (a == 0.0) H = exp(-v * v);

	/* --- H is even in v, F odd --                     -------------- */

	if (F) *F = (v < 0.0) ? -Fx : Fx;
	return H;
}
/* ------- end ---------------------------- Voigt.c ----------------- */
