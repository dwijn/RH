/* ------- file: -------------------------- test_hminus_ff.c --------

	   H-minus free-free absorption coefficient (hydrogen.c:
	   Hminus_ff_kappa, Hminus_ff; N-45), the fits of John (1988,
	   A&A 193, 189, Tables 3a and 3b):

	   -- they reproduce Table 1 of Bell & Berrington (1987, J. Phys. B
		  20, 801) at theta = 0.5, 1, 2 and 3.6 from 182.3 nm to
		  15.2 micron to 1.5 %;
	   -- the lambda^2 tail joins Table 3b at 182.3 nm, and Tables 3a
		  and 3b agree to 0.6 % at 364.5 nm for theta >= 0.5;
	   -- theta above 3.6 gives the values at 3.6;
	   -- positive and rising with theta for theta = 0.05-3.6 from
		  50 nm to 1 mm; zero for lambda <= 0;
	   -- Hminus_ff = n(H I) n_e k T kappa.
	   --                                              -------------- */

#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "background.h"
#include "constant.h"
#include "check.h"
#include "state.h"

/* --- 1.0E-26 cm^4/dyn in m^5/J --                    -------------- */

#define K26 1.0E-29

#define NBB 22

/* --- Bell & Berrington (1987), Table 1: wavelengths [Angstrom] and
	   k [1e-26 cm^4/dyn] at theta = 0.5, 1, 2, 3.6 --  -------------- */

static const double BB_lambda[NBB] = {
	151890, 113918, 91134, 45567, 30378, 22784, 18227, 15189,
	13019, 11392, 10126, 9113, 7595, 6510, 5696, 5063,
	4557, 3645, 3038, 2604, 2278, 1823};

static const double BB_theta[4] = {0.5, 1.0, 2.0, 3.6};

static const double BB_k[4][NBB] = {
	{75.1, 42.3, 27.0, 6.79, 3.04, 1.73, 1.11, .781, .579, .448, .358,
	 .293, .208, .154, .121, .0965, .0791, .0520, .0364, .0277, .0228,
	 .0178},
	{144, 80.8, 51.9, 13.1, 5.86, 3.37, 2.17, 1.52, 1.13, .871, .702,
	 .576, .409, .309, .241, .195, .161, .108, .0789, .0615, .0499,
	 .0402},
	{253, 142, 91.4, 23.0, 10.3, 5.95, 3.80, 2.69, 2.02, 1.57, 1.28,
	 1.06, .774, .597, .475, .390, .327, .225, .169, .136, .110, .0896},
	{388, 219, 140, 35.0, 15.6, 9.06, 5.97, 4.31, 3.31, 2.60, 2.17,
	 1.83, 1.36, 1.07, .861, .711, .602, .418, .318, .262, .211, .172}};

int main(void)
{
	int    i, j, k;
	double lambda, theta, r, prev;

	/* --- Bell & Berrington, Table 1 --              -------------- */

	for (j = 0; j < 4; j++) {
		for (i = 0; i < NBB; i++) {
			lambda = BB_lambda[i] / 10.0;
			r      = Hminus_ff_kappa(lambda, BB_theta[j]) / (BB_k[j][i] * K26);
			CHECK(r > 0.985 && r < 1.016,
				  "John/BB = %.4f at %.1f nm, theta = %.1f", r, lambda,
				  BB_theta[j]);
		}
	}

	/* --- Joins at 182.3 and 364.5 nm --              -------------- */

	for (theta = 0.1; theta <= 3.6; theta += 0.1) {
		CHECK_CLOSE(Hminus_ff_kappa(182.3 * (1.0 - 1.0E-12), theta),
					Hminus_ff_kappa(182.3, theta), 1.0E-9, 0.0,
					"join at 182.3 nm, theta = %.1f", theta);
		CHECK_CLOSE(Hminus_ff_kappa(100.0, theta),
					Hminus_ff_kappa(182.3, theta) * SQ(100.0 / 182.3), 1.0E-12,
					0.0, "lambda^2 tail at 100 nm, theta = %.1f", theta);
		if (theta >= 0.5 - 1.0E-9) {
			r = Hminus_ff_kappa(364.5, theta) /
				Hminus_ff_kappa(364.5 * (1.0 - 1.0E-12), theta);
			CHECK(fabs(r - 1.0) < 0.006, "3a/3b = %.4f at theta = %.1f", r,
				  theta);
		}
	}

	/* --- theta clamped at 3.6; positive and rising below --  ------ */

	for (lambda = 50.0; lambda <= 1.0E6; lambda *= 1.1) {
		CHECK(Hminus_ff_kappa(lambda, 5.0) == Hminus_ff_kappa(lambda, 3.6),
			  "theta = 5 not clamped at %.1f nm", lambda);
		CHECK(Hminus_ff_kappa(lambda, 20.0) == Hminus_ff_kappa(lambda, 3.6),
			  "theta = 20 not clamped at %.1f nm", lambda);
		prev = 0.0;
		for (theta = 0.05; theta <= 3.6 + 1.0E-9; theta += 0.05) {
			r = Hminus_ff_kappa(lambda, theta);
			CHECK(r > prev, "not rising at %.1f nm, theta = %.2f", lambda,
				  theta);
			prev = r;
		}
	}
	CHECK(Hminus_ff_kappa(0.0, 1.0) == 0.0, "%s", "lambda = 0");
	CHECK(Hminus_ff_kappa(-1.0, 1.0) == 0.0, "%s", "lambda < 0");

	/* --- Hminus_ff on a small atmosphere --          -------------- */

	{
		enum { NSPACE = 5 };
		static double T[NSPACE]  = {1000.0, 1400.0, 5000.0, 10000.0, 20000.0},
					  ne[NSPACE] = {1.0E15, 1.0E17, 1.0E19, 1.0E21, 1.0E22},
					  nH[NSPACE] = {1.0E23, 1.0E22, 1.0E21, 1.0E20, 1.0E19},
					  chi[NSPACE];
		static double *n[1] = {nH};
		static Atom    H;
		double         want, lambdas[] = {150.0, 300.0, 500.0, 2000.0, 2.0E5};

		H.n              = n;
		rh->atmos.H      = &H;
		rh->atmos.Nspace = NSPACE;
		rh->atmos.T      = T;
		rh->atmos.ne     = ne;
		for (i = 0; i < 5; i++) {
			CHECK(Hminus_ff(lambdas[i], chi), "%s", "Hminus_ff returns FALSE");
			for (k = 0; k < NSPACE; k++) {
				want = nH[k] * ne[k] * KBOLTZMANN * T[k] *
					   Hminus_ff_kappa(lambdas[i], THETA0 / T[k]);
				CHECK_CLOSE(chi[k], want, 1.0E-13, 0.0,
							"Hminus_ff at %.0f nm, T = %.0f K", lambdas[i],
							T[k]);
			}
		}
	}
	return check_summary("test_hminus_ff");
}
/* ------- end ---------------------------- test_hminus_ff.c -------- */
