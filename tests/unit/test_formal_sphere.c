/* ------- file: -------------------------- test_formal_sphere.c ----

	   Spherical formal solution (rhsphere/feautrier.c along the rays of
	   rhsphere/geometry.c getRays) in a shell of constant opacity chi
	   and constant source function S = B(T), where the intensities
	   along every ray are known in closed form. A ray with impact
	   parameter p has path coordinate s = sqrt(r^2 - p^2) (r from the
	   center), s0 at the outer radius:

	   -- shell rays (p above the core; Feautrier reflects at the
		  tangent point s = 0):
			I-(s) = S (1 - exp(-chi (s0 - s))),
			I+(s) = S (1 - exp(-chi (s0 + s)));
	   -- core rays (thermalized at the core, I+ = B there, so with
		  S = B everywhere): I+(s) = S, I-(s) as above.

	   Feautrier returns P = (I+ + I-)/2 and the emergent I+(s0).
	   Checks the order of convergence and the errors as the shell is
	   refined, for STANDARD and FEAUTRIER_HERMITE. -v prints the
	   measured errors.
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
#include "check.h"

extern Geometry    geometry;
extern Atmosphere  atmos;
extern Spectrum    spectrum;
extern CommandLine commandline;

double Feautrier(int nspect, int mu, double *chi, double *S,
				 enum FeautrierOrder F_order, double *P, double *Psi);

#define RADIUS 1.0 /* core radius */
#define HEIGHT 1.0 /* shell thickness */
#define CHI    10.0
#define TEMP   6000.0
#define LAMBDA 500.0

/* --- Shell of Nr radii, logarithmic in radial optical depth from
	   1e-3 to CHI * HEIGHT, r[k] the height above the core -- ------- */

static void shell_init(int Nr)
{
	geometry.Nradius        = Nr;
	geometry.Radius         = RADIUS;
	geometry.Ncore          = 5;
	geometry.rboundary[TOP] = ZERO;
	geometry.r              = malloc(Nr * sizeof(double));
	for (int k = 0; k < Nr; k++) {
		double tau    = (k == 0) ? 0.0 : 1.0e-3 * pow(CHI * HEIGHT / 1.0e-3, (double) (k - 1) / (Nr - 2));
		geometry.r[k] = HEIGHT - tau / CHI;
	}
	geometry.r[Nr - 1] = 0.0;

	atmos.T = malloc(Nr * sizeof(double));
	for (int k = 0; k < Nr; k++) atmos.T[k] = TEMP;

	getRays(&geometry);
}

static void shell_free(void)
{
	for (int mu = 0; mu < geometry.Nrays; mu++) {
		free(geometry.rays[mu].s);
		free(geometry.rays[mu].xmu);
		free(geometry.rays[mu].wmu);
	}
	free(geometry.rays);
	free(geometry.r);
	free(atmos.T);
}

/* --- Shell rays whose tangent point lies above radial optical depth
	   TAU_GRAZING are grazing rays near the top: refining the radial
	   grid moves them closer to the surface instead of sampling them
	   more finely (the outermost always has 3 points, PTOP = 2), so
	   their errors do not converge and are bounded separately -- --- */

#define TAU_GRAZING 0.1

static bool_t grazing(const Ray *ray)
{
	return ray->type == SHELL_RAY && (HEIGHT - ray->p) * CHI < TAU_GRAZING;
}

typedef struct {
	double P, emergent, grazing;
} Errors;

/* --- Worst relative errors of P and the emergent intensity over all
	   rays --                                          -------------- */

static Errors shell_errors(int Nr, enum FeautrierOrder order, double S0)
{
	Errors  e   = {0.0, 0.0, 0.0};
	double *chi = malloc(Nr * sizeof(double)), *S = malloc(Nr * sizeof(double)),
		   *P = malloc(Nr * sizeof(double));

	shell_init(Nr);
	for (int k = 0; k < Nr; k++) {
		chi[k] = CHI;
		S[k]   = S0;
	}
	for (int mu = 0; mu < geometry.Nrays; mu++) {
		Ray   *ray   = &geometry.rays[mu];
		double s0    = ray->s[0],
			   Iplus = Feautrier(0, mu, chi, S, order, P, NULL),
			   Iemit = (ray->type == SHELL_RAY) ? S0 * (1.0 - exp(-2.0 * CHI * s0)) : S0;

		for (int k = 0; k < ray->Ns; k++) {
			double s      = ray->s[k],
				   Iminus = S0 * (1.0 - exp(-CHI * (s0 - s))),
				   Ip     = (ray->type == SHELL_RAY) ? S0 * (1.0 - exp(-CHI * (s0 + s))) : S0,
				   Pe     = 0.5 * (Ip + Iminus);

			if (grazing(ray))
				e.grazing = MAX(e.grazing, fabs(P[k] - Pe) / Pe);
			else
				e.P = MAX(e.P, fabs(P[k] - Pe) / Pe);
		}
		if (grazing(ray))
			e.grazing = MAX(e.grazing, fabs(Iplus - Iemit) / Iemit);
		else
			e.emergent = MAX(e.emergent, fabs(Iplus - Iemit) / Iemit);
	}
	shell_free();
	free(chi);
	free(S);
	free(P);
	return e;
}

/* --- Measured (gcc, x86-64), Nr = 41 to 321: order 1.39 (P) and
	   1.60 (emergent) for STANDARD, 1.45 and 1.18 for HERMITE, below
	   their formal 2 and 4: the rays are sampled where they cross the
	   radii, s ~ sqrt(r - p), which is coarse near the tangent points.
	   Nr = 161: P 1.1e-3 / 7.9e-4, emergent 5.4e-4 / 1.6e-4; grazing
	   rays up to 7.2e-3 / 4.9e-3 (at Nr = 41). Bounds about twice
	   the errors, 0.3 below the orders --             -------------- */

typedef struct {
	const char         *name;
	enum FeautrierOrder order;
	double              min_order, max_P, max_emergent; /* Nr = 161 */
	double              max_grazing;                    /* all Nr */
} Method;

/* clang-format off */
static const Method methods[] = {
	{"FEAUTRIER", STANDARD, 1.1, 2.2e-3, 1.1e-3, 1.5e-2},
	{"FEAUTRIER_HERMITE", FEAUTRIER_HERMITE, 0.9, 1.6e-3, 3.2e-4, 1.0e-2}};
/* clang-format on */

#define N_METHODS (int) (sizeof(methods) / sizeof(methods[0]))

int main(int argc, char *argv[])
{
	bool_t verbose = (argc > 1 && strcmp(argv[1], "-v") == 0);
	double S0, T = TEMP;

	commandline.logfile = stdout;
	commandline.quiet   = TRUE;

	spectrum.lambda    = malloc(sizeof(double));
	spectrum.lambda[0] = LAMBDA;
	Planck(1, &T, LAMBDA, &S0);

	for (int n = 0; n < N_METHODS; n++) {
		const Method *m = &methods[n];
		enum { NREFINE = 4 };
		Errors e[NREFINE];
		double order_P = HUGE_VAL, order_em = HUGE_VAL, grazing_max = 0.0;

		for (int r = 0; r < NREFINE; r++) {
			int Nr = 40 * (1 << r) + 1;

			e[r]        = shell_errors(Nr, m->order, S0);
			grazing_max = MAX(grazing_max, e[r].grazing);
			if (verbose)
				printf("  %-18s Nr %4d  P %.3e  emergent %.3e  grazing rays %.3e\n",
					   m->name, Nr, e[r].P, e[r].emergent, e[r].grazing);
			if (r > 0) {
				order_P  = MIN(order_P, log2(e[r - 1].P / e[r].P));
				order_em = MIN(order_em, log2(e[r - 1].emergent / e[r].emergent));
			}
		}
		if (verbose)
			printf("%-18s order P %.2f, emergent %.2f\n", m->name, order_P, order_em);

		CHECK(order_P >= m->min_order, "%s: order of P %.2f < %.2f", m->name,
			  order_P, m->min_order);
		CHECK(order_em >= m->min_order, "%s: order of emergent I %.2f < %.2f",
			  m->name, order_em, m->min_order);
		CHECK(e[2].P <= m->max_P, "%s: error of P %.3g > %.3g", m->name, e[2].P,
			  m->max_P);
		CHECK(e[2].emergent <= m->max_emergent, "%s: emergent error %.3g > %.3g",
			  m->name, e[2].emergent, m->max_emergent);
		CHECK(grazing_max <= m->max_grazing, "%s: grazing-ray error %.3g > %.3g",
			  m->name, grazing_max, m->max_grazing);
	}
	return check_summary("test_formal_sphere");
}
