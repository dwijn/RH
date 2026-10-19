/* ------- file: -------------------------- test_formal_stokes.c ----

	   Polarized 1-D formal solver (rhf1d/bezier_1D.c,
	   Piece_Stokes_Bezier3_1D, DELO-Bezier) on slabs with known
	   solutions.

	   Along a ray with direction cosine mu, in optical depth tau of
	   Stokes I (increasing inward), the transfer equation for the
	   Stokes vector I is

		 mu dI/dtau = K (I - S),   K = 1 + K',

	   with K' the reduced propagation matrix of StokesK() (absorption
	   eta_Q,U,V / chi_I, magneto-optical rho_Q,U,V / chi_I) and S the
	   source vector of Unno-Rachkovsky form, S = (S(tau), 0, 0, 0)
	   (an unpolarized source function, as in LTE). RH is passed the
	   emission vector over chi_I, K S.

	   For constant K and polynomial S(tau) the solution is the matrix
	   form of the scalar one in slab_exact.h, with mu -> mu K^-1:

		 outgoing:  I+(tau) = sum_j (mu K^-1)^j S^(j)(tau),
		 incoming, I-(tau0) = 0 at the top:
				   I-(tau) = P-(tau) - exp(-K (tau - tau0) / mu) P-(tau0),
				   P-(tau) = sum_j (-mu K^-1)^j S^(j)(tau).

	   Checks:
	   -- with K' = 0 the polarized solver must give the unpolarized
		  Bezier3 solution of Piecewise_Bezier3_1D, and Q = U = V = 0;
	   -- I+ exact (to rounding) for linear S and constant K, where the
		  solution is linear in tau (away from the lower boundary,
		  where the solver sets Q = U = V = 0);
	   -- convergence order of I+ and I- for cubic S, against the exact
		  solution;
	   -- convergence order with depth-dependent K, from three grids
		  (no closed form; Richardson).
	   -v prints the measured errors.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "inputs.h"
#include "bezier.h"
#include "check.h"
#include "stokes_exact.h"

extern Geometry   geometry;
extern Atmosphere atmos;
extern Spectrum   spectrum;
extern InputData  input;

void Piecewise_Bezier3_1D(int nspect, int mu, bool_t to_obs, double *chi,
						  double *S, double *I, double *Psi);
void Piece_Stokes_Bezier3_1D(int nspect, int mu, bool_t to_obs, double *chi,
							 double **S, double **I, double *Psi);

#define NMU 3
static double mus[NMU] = {0.1, 0.5, 1.0};

/* --- Slab of Ndep points logarithmic in tau, constant chi = 1
	   (height = -tau). K'(tau) = f(tau) K_ELEMENTS, f = 1 or varying -- */

typedef struct {
	int     Ndep;
	double *tau, *chi, **S, **I, *Is;
	bool_t  vary;
} Slab;

#define TAU_MIN 1.0e-4
#define TAU_MAX 1.0e2

static double K_scale(const Slab *slab, double tau)
{
	return slab->vary ? 1.0 / (1.0 + 0.5 * tau) : 1.0;
}

static void slab_init(Slab *slab, int Ndep, bool_t vary, double kscale)
{
	slab->Ndep = Ndep;
	slab->vary = vary;
	slab->tau  = malloc(Ndep * sizeof(double));
	slab->chi  = malloc(Ndep * sizeof(double));
	slab->Is   = malloc(Ndep * sizeof(double));
	slab->S    = matrix_double(4, Ndep);
	slab->I    = matrix_double(4, Ndep);

	geometry.Ndep   = Ndep;
	atmos.Nspace    = Ndep;
	geometry.height = realloc(geometry.height, Ndep * sizeof(double));

	ActiveSet *as = &spectrum.as[0];
	as->chi_c     = realloc(as->chi_c, 4 * Ndep * sizeof(double));
	as->chip_c    = realloc(as->chip_c, 3 * Ndep * sizeof(double));

	for (int k = 0; k < Ndep; k++) {
		double t = TAU_MIN * pow(TAU_MAX / TAU_MIN, (double) k / (Ndep - 1)),
			   f = kscale * K_scale(slab, t);

		slab->tau[k]       = t;
		slab->chi[k]       = 1.0;
		geometry.height[k] = -t;

		/* --- StokesK: K[0][1..3] = chi_c[Q,U,V], K[1][2] = chip_c[2N]
			   (rho_V), K[1][3] = -chip_c[N] (rho_U), K[2][3] = chip_c[0]
			   (rho_Q), all over chi_I = 1 --          -------------- */

		as->chi_c[k]             = 1.0;
		as->chi_c[Ndep + k]      = f * K_ELEMENTS[0];
		as->chi_c[2 * Ndep + k]  = f * K_ELEMENTS[1];
		as->chi_c[3 * Ndep + k]  = f * K_ELEMENTS[2];
		as->chip_c[k]            = f * K_ELEMENTS[3];
		as->chip_c[Ndep + k]     = f * K_ELEMENTS[4];
		as->chip_c[2 * Ndep + k] = f * K_ELEMENTS[5];
	}
}

static void slab_free(Slab *slab)
{
	free(slab->tau);
	free(slab->chi);
	free(slab->Is);
	freeMatrix((void **) slab->S);
	freeMatrix((void **) slab->I);
}

/* --- Emission vector over chi_I, K S for S = (S(tau), 0, 0, 0) -- */

static void set_source(Slab *slab, const double s[NPOLY], double kscale)
{
	for (int k = 0; k < slab->Ndep; k++) {
		Mat K;
		full_K(K_ELEMENTS, kscale * K_scale(slab, slab->tau[k]), K);
		double St = S_deriv(s, 0, slab->tau[k]);
		for (int n = 0; n < 4; n++) slab->S[n][k] = K[n][0] * St;
	}
	for (int mu = 0; mu < NMU; mu++) {
		geometry.Itop[0][mu]    = 0.0;
		geometry.Ibottom[0][mu] = S_deriv(s, 0, TAU_MAX) + mus[mu] * s[1] +
								  2.0 * mus[mu] * mus[mu] * s[2] +
								  6.0 * pow(mus[mu], 3) * s[3];
	}
}

/* --- Largest error of any Stokes component, relative to Stokes I,
	   at the emergent point and at interior points with tau below
	   TAU_CHECK (deeper, the lower boundary condition Q = U = V = 0
	   differs from the semi-infinite solution) --    -------------- */

#define TAU_CHECK 10.0

typedef struct {
	double emergent, interior;
} Errors;

static void errors_exact(Errors *e, Slab *slab, const double s[NPOLY],
						 double kscale, int mu, bool_t to_obs)
{
	Mat K;
	full_K(K_ELEMENTS, kscale, K);
	Piece_Stokes_Bezier3_1D(0, mu, to_obs, slab->chi, slab->S, slab->I, NULL);
	for (int k = 0; k < slab->Ndep && slab->tau[k] <= TAU_CHECK; k++) {
		Vec    Ie;
		double err = 0.0;

		I_exact_pol(s, K, slab->tau[k], slab->tau[0], mus[mu], to_obs, Ie);
		if (k == 0 && !to_obs) continue; /* boundary value */
		for (int n = 0; n < 4; n++)
			err = MAX(err, fabs(slab->I[n][k] - Ie[n]) / fabs(Ie[0]));
		if (k == 0)
			e->emergent = MAX(e->emergent, err);
		else
			e->interior = MAX(e->interior, err);
	}
}

static const double s_linear[NPOLY] = {1.0, 1.0, 0.0, 0.0};
static const double s_cubic[NPOLY]  = {1.0, 1.0, 0.5, 0.1};

/* --- K' = 0: same as the unpolarized solver -- ----------------- */

static void check_unpolarized(bool_t verbose)
{
	Slab   slab;
	double dI = 0.0, dQUV = 0.0;

	slab_init(&slab, 81, FALSE, 0.0);
	set_source(&slab, s_cubic, 0.0);
	for (int mu = 0; mu < NMU; mu++)
		for (int to_obs = 0; to_obs <= 1; to_obs++) {
			Piece_Stokes_Bezier3_1D(0, mu, to_obs, slab.chi, slab.S, slab.I, NULL);
			Piecewise_Bezier3_1D(0, mu, to_obs, slab.chi, slab.S[0], slab.Is, NULL);
			for (int k = 0; k < slab.Ndep; k++) {
				dI = MAX(dI, fabs(slab.I[0][k] - slab.Is[k]) / MAX(fabs(slab.Is[k]), 1e-300));
				for (int n = 1; n < 4; n++) dQUV = MAX(dQUV, fabs(slab.I[n][k]));
			}
		}
	if (verbose)
		printf("K' = 0: polarized vs unpolarized Bezier3 I %.2e, max |QUV| %.2e\n",
			   dI, dQUV);
	CHECK(dI <= 1.0e-13, "K' = 0: polarized and unpolarized I differ by %.3g", dI);
	CHECK(dQUV == 0.0, "K' = 0: Q, U, V not zero (%.3g)", dQUV);
	slab_free(&slab);
}

/* --- Linear S, constant K: I+ is linear in tau and exact -- ---- */

static void check_linear(bool_t verbose)
{
	Slab   slab;
	Errors e = {0.0, 0.0};

	slab_init(&slab, 41, FALSE, 1.0);
	set_source(&slab, s_linear, 1.0);
	for (int mu = 0; mu < NMU; mu++) errors_exact(&e, &slab, s_linear, 1.0, mu, TRUE);
	if (verbose)
		printf("linear S: I+ emergent %.2e, interior %.2e\n", e.emergent, e.interior);
	CHECK(e.emergent <= 1.0e-13 && e.interior <= 1.0e-13,
		  "linear S: I+ not exact (emergent %.3g, interior %.3g)", e.emergent,
		  e.interior);
	slab_free(&slab);
}

/* --- Cubic S, constant K: errors and order against the exact
	   solution, both directions. Measured (gcc, x86-64): I+ emergent
	   order 4.00, 9.9e-7 at Ndep = 161 (the unpolarized solver: 9.3e-7),
	   interior 2.1e-5; I- interior order 2.92, 2.3e-5 (the one-sided
	   derivatives at the ends, carried along the ray, as unpolarized).
	   Bounds about twice the errors --                -------------- */

enum { NREFINE = 4 };

static void check_cubic(bool_t verbose)
{
	const double min_order[2] = {2.6, 3.7}, max_err[2] = {5.0e-5, 2.0e-6};

	for (int to_obs = 0; to_obs <= 1; to_obs++) {
		Errors e[NREFINE];
		double order = HUGE_VAL;

		for (int r = 0; r < NREFINE; r++) {
			Slab slab;
			int  Ndep = 40 * (1 << r) + 1;

			e[r].emergent = e[r].interior = 0.0;
			slab_init(&slab, Ndep, FALSE, 1.0);
			set_source(&slab, s_cubic, 1.0);
			for (int mu = 0; mu < NMU; mu++)
				errors_exact(&e[r], &slab, s_cubic, 1.0, mu, to_obs);
			slab_free(&slab);
			if (verbose)
				printf("  cubic S, %s Ndep %4d  emergent %.3e  interior %.3e\n",
					   to_obs ? "I+" : "I-", Ndep, e[r].emergent, e[r].interior);
			if (r > 0) {
				double ee = to_obs ? e[r].emergent : e[r].interior,
					   ep = to_obs ? e[r - 1].emergent : e[r - 1].interior;
				order     = MIN(order, log2(ep / ee));
			}
		}
		double err161 = to_obs ? e[2].emergent : e[2].interior;
		CHECK(err161 <= max_err[to_obs], "cubic S, %s: error %.3g > %.3g at Ndep = 161",
			  to_obs ? "I+ emergent" : "I- interior", err161, max_err[to_obs]);
		if (verbose)
			printf("cubic S, %s: order %.2f\n", to_obs ? "I+ emergent" : "I- interior",
				   order);
		CHECK(order >= min_order[to_obs], "cubic S, %s: order %.2f < %.2f",
			  to_obs ? "I+" : "I-", order, min_order[to_obs]);
	}
}

/* --- Depth-dependent K: order of the emergent Stokes vector from
	   grids N, 2N, 4N. Measured 4.00 (1.02 with the sign error of dK
	   at the central point, B-25) --                  -------------- */

static void check_varying_K(bool_t verbose, double min_order)
{
	enum { NR = 5 };
	double Iem[NR][NMU][4], order = HUGE_VAL;

	for (int r = 0; r < NR; r++) {
		Slab slab;
		slab_init(&slab, 40 * (1 << r) + 1, TRUE, 1.0);
		set_source(&slab, s_cubic, 1.0);
		for (int mu = 0; mu < NMU; mu++) {
			Piece_Stokes_Bezier3_1D(0, mu, TRUE, slab.chi, slab.S, slab.I, NULL);
			for (int n = 0; n < 4; n++) Iem[r][mu][n] = slab.I[n][0];
		}
		slab_free(&slab);
	}
	for (int r = 2; r < NR; r++) {
		double d1 = 0.0, d2 = 0.0;
		for (int mu = 0; mu < NMU; mu++)
			for (int n = 0; n < 4; n++) {
				d1 = MAX(d1, fabs(Iem[r - 1][mu][n] - Iem[r - 2][mu][n]) / Iem[r][mu][0]);
				d2 = MAX(d2, fabs(Iem[r][mu][n] - Iem[r - 1][mu][n]) / Iem[r][mu][0]);
			}
		if (verbose)
			printf("  varying K: Ndep %4d  change %.3e  order %.2f\n", 40 * (1 << r) + 1,
				   d2, log2(d1 / d2));
		order = MIN(order, log2(d1 / d2));
	}
	if (verbose) printf("varying K: order %.2f\n", order);
	CHECK(order >= min_order, "varying K: order %.2f < %.2f", order, min_order);
}

int main(int argc, char *argv[])
{
	bool_t verbose = (argc > 1 && strcmp(argv[1], "-v") == 0);

	geometry.Nrays             = NMU;
	geometry.muz               = mus;
	geometry.Itop              = matrix_double(1, NMU);
	geometry.Ibottom           = matrix_double(1, NMU);
	geometry.vboundary[TOP]    = ZERO;
	geometry.vboundary[BOTTOM] = IRRADIATED;

	spectrum.as                      = calloc(1, sizeof(ActiveSet));
	atmos.backgrflags                = calloc(1, sizeof(flags));
	atmos.backgrflags[0].ispolarized = TRUE;
	input.magneto_optical            = TRUE;

	check_unpolarized(verbose);
	check_linear(verbose);
	check_cubic(verbose);
	check_varying_K(verbose, 3.7);

	return check_summary("test_formal_stokes");
}
