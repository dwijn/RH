/* ------- file: -------------------------- test_formal_1d.c --------

	   1-D formal solvers on slabs with analytic solutions.

	   The source function is a polynomial in optical depth, where the
	   intensities in both directions are known in closed form
	   (slab_exact.h). Feautrier solves for P = (I+ + I-) / 2.

	   Checks:
	   -- exact (to rounding) where the method's interpolant contains
		  the solution: linear S for the linear and cubic Bezier
		  short characteristics;
	   -- convergence order as the grid is refined, on a cubic S and
		  on an opacity that varies exponentially with height;
	   -- the diagonal operator Psi of the linear solvers is the exact
		  derivative dI(k)/dS(k).
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "bezier.h"
#include "check.h"
#include "slab_exact.h"
#include "bezier_reference.h"

extern Geometry   geometry;
extern Atmosphere atmos;
extern Spectrum   spectrum;

/* --- The solvers under test (rhf1d/) --              -------------- */

void   Piecewise_Linear_1D(int nspect, int mu, bool_t to_obs, double *chi,
						   double *S, double *I, double *Psi);
void   Piecewise_Bezier3_1D(int nspect, int mu, bool_t to_obs, double *chi,
							double *S, double *I, double *Psi);
double Feautrier(int nspect, int mu, double *chi, double *S,
				 enum FeautrierOrder F_order, double *P, double *Psi);

#define NMU 3
static double mus[NMU] = {0.1, 0.5, 1.0};

/* --- An atmosphere of Ndep points, logarithmic in optical depth as
	   in a stellar atmosphere. Opacity either constant (height = -tau)
	   or exponential in height (chi = tau = exp(-height)), which the
	   solvers must integrate. With IRRADIATED lower boundary the exact
	   I+ enters at the bottom, so the finite slab has the
	   semi-infinite solution --                       -------------- */

typedef struct {
	int     Ndep;
	double *tau, *chi, *S, *I, *Psi, *P, *dI;
} Slab;

static void slab_init(Slab *slab, int Ndep, bool_t exp_chi)
{
	const double tau_min = 1.0e-4, tau_max = 1.0e2;

	slab->Ndep = Ndep;
	slab->tau  = malloc(Ndep * sizeof(double));
	slab->chi  = malloc(Ndep * sizeof(double));
	slab->S    = malloc(Ndep * sizeof(double));
	slab->I    = malloc(Ndep * sizeof(double));
	slab->Psi  = malloc(Ndep * sizeof(double));
	slab->P    = malloc(Ndep * sizeof(double));
	slab->dI   = malloc(Ndep * sizeof(double));

	geometry.Ndep   = Ndep;
	geometry.height = realloc(geometry.height, Ndep * sizeof(double));
	for (int k = 0; k < Ndep; k++) {
		slab->tau[k] = tau_min * pow(tau_max / tau_min, (double) k / (Ndep - 1));
		if (exp_chi) {
			slab->chi[k]       = slab->tau[k];
			geometry.height[k] = -log(slab->tau[k]);
		} else {
			slab->chi[k]       = 1.0;
			geometry.height[k] = -slab->tau[k];
		}
	}
}

static void slab_free(Slab *slab)
{
	free(slab->tau);
	free(slab->chi);
	free(slab->S);
	free(slab->I);
	free(slab->Psi);
	free(slab->P);
	free(slab->dI);
}

static void set_source(Slab *slab, const double s[NPOLY])
{
	for (int k = 0; k < slab->Ndep; k++) {
		double t   = slab->tau[k];
		slab->S[k] = s[0] + t * (s[1] + t * (s[2] + t * s[3]));
	}
	for (int mu = 0; mu < NMU; mu++) {
		geometry.Itop[0][mu] = 0.0;
		geometry.Ibottom[0][mu] =
			I_exact(s, slab->tau[slab->Ndep - 1], 0.0, mus[mu], TRUE);
	}
}

/* --- Errors relative to the exact solution: of the emergent
	   intensity (I+ at the top; for Feautrier its return value), at
	   the interior points, and at the two points nearest to either
	   end, where the solvers use one-sided differences -- --------- */

typedef struct {
	double emergent, interior, ends;
	int    worst_k; /* largest interior error */
} Errors;

static double rel_err(double got, double want)
{
	return fabs(got - want) / MAX(fabs(want), 1.0e-300);
}

static void accumulate(Errors *e, const double *got, const double *want,
					   int N)
{
	for (int k = 0; k < N; k++) {
		double err = rel_err(got[k], want[k]);

		if (k <= 1 || k >= N - 2)
			e->ends = MAX(e->ends, err);
		else if (err > e->interior) {
			e->interior = err;
			e->worst_k  = k;
		}
	}
}

typedef void (*ShortChar)(int, int, bool_t, double *, double *, double *,
						  double *);

static void sc_errors(Errors *e, ShortChar solver, Slab *slab,
					  const double s[NPOLY], int mu, bool_t to_obs)
{
	solver(0, mu, to_obs, slab->chi, slab->S, slab->I, NULL);
	for (int k = 0; k < slab->Ndep; k++)
		slab->dI[k] = I_exact(s, slab->tau[k], slab->tau[0], mus[mu], to_obs);
	accumulate(e, slab->I, slab->dI, slab->Ndep);
	if (to_obs) e->emergent = MAX(e->emergent, rel_err(slab->I[0], slab->dI[0]));
}

/* --- Feautrier: P = (I+ + I-) / 2 and the emergent I+ -- ------- */

static void feautrier_errors(Errors *e, Slab *slab, const double s[NPOLY],
							 int mu, enum FeautrierOrder order)
{
	double Iplus = Feautrier(0, mu, slab->chi, slab->S, order, slab->P, NULL);

	for (int k = 0; k < slab->Ndep; k++)
		slab->dI[k] = 0.5 * (I_exact(s, slab->tau[k], 0.0, mus[mu], TRUE) +
							 I_exact(s, slab->tau[k], slab->tau[0], mus[mu], FALSE));
	accumulate(e, slab->P, slab->dI, slab->Ndep);
	e->emergent = MAX(e->emergent,
					  rel_err(Iplus, I_exact(s, slab->tau[0], 0.0, mus[mu], TRUE)));
}

/* --- The methods. Bounds from the convergence on a cubic S at
	   Ndep = 41 to 321 (5 to 40 points per decade): the order 0.3
	   below the measured one, the error at Ndep = 161 twice the
	   measured one; [0] constant, [1] exponential opacity -- ------ */

typedef struct {
	const char         *name;
	ShortChar           solver; /* NULL: Feautrier */
	enum FeautrierOrder order;
	bool_t              exact_linear;
	double              order_emergent[2], order_interior[2];
	double              err_emergent[2], err_interior[2], err_ends[2];
} Method;

/* --- Measured (gcc, x86-64), constant / exponential opacity:
	   LINEAR     order 2.00/2.00 emergent, 1.91/1.89 interior;
				  error 9.7e-4/1.4e-4 emergent, 3.2e-3/3.1e-3 interior,
				  2.2e-3/2.1e-3 ends
	   BEZIER3    order 4.00/4.00 emergent, 2.16/1.84 interior (the
				  interior error comes from the one-sided derivatives at
				  the ends, carried along the ray);
				  error 9.3e-7/3.4e-7 emergent, 3.6e-5/1.4e-4 interior,
				  1.5e-3/1.5e-3 ends
	   FEAUTRIER  order 2.00/1.99; error 4.4e-4/3.9e-4 everywhere
	   HERMITE    order 3.97/1.98 emergent, 4.00/1.98 interior (with
				  varying opacity the trapezoidal optical depth limits it
				  to second order); error 3.6e-7/8.3e-4 emergent,
				  7.8e-7/8.3e-4 interior, 3.6e-7/8.3e-4 ends -- -------- */

/* clang-format off */
static const Method methods[] = {
	{"LINEAR", Piecewise_Linear_1D, STANDARD, TRUE,
	 {1.7, 1.7}, {1.6, 1.6}, {2e-3, 3e-4}, {7e-3, 7e-3}, {5e-3, 5e-3}},
	{"BEZIER3", Piecewise_Bezier3_1D, STANDARD, TRUE,
	 {3.7, 3.7}, {1.8, 1.5}, {2e-6, 7e-7}, {8e-5, 3e-4}, {3e-3, 3e-3}},
	{"FEAUTRIER", NULL, STANDARD, FALSE,
	 {1.7, 1.7}, {1.7, 1.7}, {9e-4, 8e-4}, {9e-4, 8e-4}, {9e-4, 8e-4}},
	{"FEAUTRIER_HERMITE", NULL, FEAUTRIER_HERMITE, FALSE,
	 {3.6, 1.7}, {3.7, 1.7}, {8e-7, 1.7e-3}, {1.6e-6, 1.7e-3}, {8e-7, 1.7e-3}}};
/* clang-format on */

#define N_METHODS (int) (sizeof(methods) / sizeof(methods[0]))

/* --- Worst errors over directions --                -------------- */

static Errors method_errors(const Method *m, Slab *slab, const double s[NPOLY])
{
	Errors e = {0.0, 0.0, 0.0, 0};

	for (int mu = 0; mu < NMU; mu++) {
		if (m->solver) {
			for (int to_obs = 0; to_obs <= 1; to_obs++)
				sc_errors(&e, m->solver, slab, s, mu, to_obs);
		} else
			feautrier_errors(&e, slab, s, mu, m->order);
	}
	return e;
}

/* --- Linear S at constant opacity lies in the interpolants of both
	   short-characteristics solvers: exact up to rounding (LINEAR
	   1.3e-15) or the accuracy of the Bezier weights (1.1e-12). With
	   exponential opacity the optical depth is not exact,
	   so neither is I --                              -------------- */

static void check_exact(const Method *m, bool_t verbose)
{
	const double s[NPOLY] = {1.0, 2.0, 0.0, 0.0};
	Slab         slab;
	Errors       e;

	slab_init(&slab, 61, FALSE);
	set_source(&slab, s);
	e = method_errors(m, &slab, s);
	if (verbose)
		printf("%-18s linear S: error %.2e\n", m->name, MAX(e.interior, e.ends));
	CHECK(MAX(e.interior, e.ends) < 1.0e-11,
		  "%s, linear S: error %.3g in the interior, %.3g at the ends",
		  m->name, e.interior, e.ends);
	slab_free(&slab);
}

/* --- Order of convergence on a cubic S --          -------------- */

static void check_convergence(const Method *m, bool_t exp_chi, bool_t verbose)
{
	const double s[NPOLY] = {1.0, 1.0, 0.5, 0.1};
	const char  *chi_name = exp_chi ? "exp chi" : "const chi";
	enum { NREFINE = 4 };
	Errors e[NREFINE];
	double order_emergent = HUGE_VAL, order_interior = HUGE_VAL;

	for (int r = 0; r < NREFINE; r++) {
		Slab slab;
		int  Ndep = 40 * (1 << r) + 1;

		slab_init(&slab, Ndep, exp_chi);
		set_source(&slab, s);
		e[r] = method_errors(m, &slab, s);
		if (verbose)
			printf("  %-18s %-9s Ndep %4d  emergent %.3e  interior %.3e "
				   "(tau %.3g)  ends %.3e\n",
				   m->name, chi_name, Ndep, e[r].emergent, e[r].interior,
				   slab.tau[e[r].worst_k], e[r].ends);
		slab_free(&slab);
		if (r > 0) {
			order_emergent = MIN(order_emergent, log2(e[r - 1].emergent / e[r].emergent));
			order_interior = MIN(order_interior, log2(e[r - 1].interior / e[r].interior));
		}
	}
	if (verbose)
		printf("%-18s %-9s order emergent %.2f, interior %.2f\n", m->name,
			   chi_name, order_emergent, order_interior);

	CHECK(order_emergent >= m->order_emergent[exp_chi],
		  "%s, %s: order of emergent I %.2f < %.2f", m->name, chi_name,
		  order_emergent, m->order_emergent[exp_chi]);
	CHECK(order_interior >= m->order_interior[exp_chi],
		  "%s, %s: order in the interior %.2f < %.2f", m->name, chi_name,
		  order_interior, m->order_interior[exp_chi]);
	CHECK(e[2].emergent <= m->err_emergent[exp_chi],
		  "%s, %s: error of emergent I %.3g > %.3g", m->name, chi_name,
		  e[2].emergent, m->err_emergent[exp_chi]);
	CHECK(e[2].interior <= m->err_interior[exp_chi],
		  "%s, %s: error in the interior %.3g > %.3g", m->name, chi_name,
		  e[2].interior, m->err_interior[exp_chi]);
	CHECK(e[2].ends <= m->err_ends[exp_chi],
		  "%s, %s: error at the end points %.3g > %.3g", m->name, chi_name,
		  e[2].ends, m->err_ends[exp_chi]);
}

/* --- Diagonal operator. The linear and Feautrier solvers are linear
	   in S, so Psi(k) must equal I(k) (P(k) for Feautrier) for
	   S = unit vector e_k and zero incident radiation. Bezier3 is not
	   linear in S (Steffen-limited derivatives), and its Psi = beta +
	   theta leaves out how the control point moves with S(k): there
	   Psi is compared with a finite-difference dI(k)/dS(k) and only
	   reported --                                     -------------- */

static double psi_error(const Method *m, Slab *slab, int mu, bool_t to_obs)
{
	double  err = 0.0, *Psi = malloc(slab->Ndep * sizeof(double));
	double *S1 = malloc(slab->Ndep * sizeof(double));

	if (m->solver)
		m->solver(0, mu, to_obs, slab->chi, slab->S, slab->I, Psi);
	else
		Feautrier(0, mu, slab->chi, slab->S, m->order, slab->P, Psi);

	for (int k = 0; k < slab->Ndep; k++) {
		double dIdS;

		if (m->solver == Piecewise_Bezier3_1D) {
			double dS = 1.0e-6 * slab->S[k], I0;

			m->solver(0, mu, to_obs, slab->chi, slab->S, slab->I, NULL);
			I0 = slab->I[k];
			slab->S[k] += dS;
			m->solver(0, mu, to_obs, slab->chi, slab->S, slab->I, NULL);
			slab->S[k] -= dS;
			dIdS = (slab->I[k] - I0) / dS;
		} else {
			double Itop = geometry.Itop[0][mu], Ibottom = geometry.Ibottom[0][mu];

			for (int l = 0; l < slab->Ndep; l++)
				S1[l] = (l == k) ? 1.0 : 0.0;
			geometry.Itop[0][mu] = geometry.Ibottom[0][mu] = 0.0;
			if (m->solver) {
				m->solver(0, mu, to_obs, slab->chi, S1, slab->I, NULL);
				dIdS = slab->I[k];
			} else {
				Feautrier(0, mu, slab->chi, S1, m->order, slab->P, NULL);
				dIdS = slab->P[k];
			}
			geometry.Itop[0][mu]    = Itop;
			geometry.Ibottom[0][mu] = Ibottom;
		}
		/* --- No contribution from the upwind boundary point -- --- */

		if (Psi[k] == 0.0 && dIdS == 0.0) continue;
		err = MAX(err, rel_err(Psi[k], dIdS));
	}
	free(Psi);
	free(S1);
	return err;
}

static void check_psi(const Method *m, double bound, bool_t verbose)
{
	const double s[NPOLY] = {1.0, 1.0, 0.5, 0.1};
	double       err      = 0.0;
	Slab         slab;

	slab_init(&slab, 61, TRUE);
	set_source(&slab, s);
	for (int mu = 0; mu < NMU; mu++) {
		if (m->solver) {
			for (int to_obs = 0; to_obs <= 1; to_obs++)
				err = MAX(err, psi_error(m, &slab, mu, to_obs));
		} else
			err = MAX(err, psi_error(m, &slab, mu, FALSE));
	}
	slab_free(&slab);
	if (verbose) printf("%-18s Psi vs dI/dS: %.2e\n", m->name, err);
	CHECK(err <= bound, "%s: Psi differs from dI/dS by %.3g > %.3g", m->name,
		  err, bound);
}

/* --- Cubic Bezier integration weights against mpmath, for dt from
	   1e-8 to 1e3: relative error, absolute for eps = exp(-dt) --- */

static void check_bezier_weights(double bound, bool_t verbose)
{
	static const char *names[5] = {"eps", "alpha", "beta", "gamma", "theta"};

	for (int n = 0; n < BEZIER_REF_N; n++) {
		double dt = bezier_ref[n][0], w[5], worst = 0.0;
		int    iworst = 0;

		Bezier3_coeffs(dt, &w[1], &w[2], &w[3], &w[4], &w[0]);
		for (int i = 0; i < 5; i++) {
			double want = bezier_ref[n][i + 1],
				   err  = (i == 0) ? fabs(w[i] - want) : fabs(w[i] - want) / want;
			if (err > worst) {
				worst  = err;
				iworst = i;
			}
		}
		if (verbose)
			printf("  Bezier3_coeffs dt = %-10.3g worst %-5s %.2e\n", dt,
				   names[iworst], worst);
		CHECK(worst <= bound, "Bezier3_coeffs(%g): %s relative error %.3g > %.3g",
			  dt, names[iworst], worst, bound);
	}
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

	check_bezier_weights(1.0e-11, verbose); /* 8.5e-13 at dt = 0.32 */

	for (int n = 0; n < N_METHODS; n++) {
		if (methods[n].exact_linear) check_exact(&methods[n], verbose);
		check_convergence(&methods[n], FALSE, verbose);
		check_convergence(&methods[n], TRUE, verbose);
	}
	check_psi(&methods[0], 1.0e-14, verbose); /* measured 2e-16 */
	check_psi(&methods[1], 0.6, verbose);     /* 0.28: approximate */
	check_psi(&methods[2], 1.0e-13, verbose); /* 1.9e-15 */

	return check_summary("test_formal_1d");
}
