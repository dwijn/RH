/* ------- file: -------------------------- test_formal_2d.c --------

	   2-D formal solvers on a horizontally uniform, periodic slab, where
	   every column must reproduce the 1-D analytic intensities along the
	   ray's vertical direction cosine muz (slab_exact.h).

	   Two meshes:
	   -- "layers": columns so wide that every upwind and downwind point
		  of the short characteristics lies on the layer above or below,
		  where the interpolation is exact for uniform data: the 2-D
		  solvers then reduce to the 1-D ones;
	   -- "columns": columns narrow enough that deep layers are crossed
		  at column lines, where S, chi and I are interpolated
		  vertically between layers (quadratic, monotonic).

	   Checks: all columns agree; I+ is exact for linear S (a linear
	   function of height survives any of the interpolations), and so is
	   I- on the "layers" mesh; convergence order on a cubic S; Psi of
	   the linear solver against the response to a unit S. -v prints the
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
#include "slab_exact.h"

extern CommandLine commandline;

/* --- Rays: (muz, mux), both signs of mux and a vertical one -- ---- */

#define NRAYS 5
static const double ray_muz[NRAYS] = {0.5, 0.5, 0.15, 0.9, 1.0};
static const double ray_mux[NRAYS] = {0.6, -0.6, 0.7, -0.3, 0.0};

#define NX 6

typedef void (*Solver2D)(Geometry *, int, int, bool_t, double *, double *,
						 double *, double *);

/* --- A periodic slab of Nz layers, logarithmic in optical depth
	   from 1e-4 to 1e2, unit opacity (height = -tau), columns DX
	   apart. Meshes stay allocated (and reachable) to the end -- ---- */

typedef struct {
	Geometry geo;
	int      Nz;
	double  *tau, *chi, *S, *I, *Psi, *S1;
} Slab;

static Slab *slab_new(int Nz, double DX)
{
	Slab     *slab   = calloc(1, sizeof(Slab));
	Geometry *geo    = &slab->geo;
	int       Nspace = NX * Nz;

	slab->Nz  = Nz;
	slab->tau = malloc(Nz * sizeof(double));
	slab->chi = malloc(Nspace * sizeof(double));
	slab->S   = malloc(Nspace * sizeof(double));
	slab->S1  = malloc(Nspace * sizeof(double));
	slab->I   = malloc(Nspace * sizeof(double));
	slab->Psi = malloc(Nspace * sizeof(double));

	geo->hboundary      = PERIODIC;
	geo->bvalue[TOP]    = ZERO;
	geo->bvalue[BOTTOM] = IRRADIATED;
	geo->Nx             = NX;
	geo->Nz             = Nz;
	geo->Nrays          = NRAYS;

	geo->x  = malloc(NX * sizeof(double));
	geo->dx = malloc(NX * sizeof(double));
	for (int l = 0; l < NX; l++) {
		geo->x[l]  = l * DX;
		geo->dx[l] = DX;
	}
	geo->z  = malloc(Nz * sizeof(double));
	geo->dz = malloc(Nz * sizeof(double));
	for (int k = 0; k < Nz; k++) {
		slab->tau[k] = 1.0e-4 * pow(1.0e6, (double) k / (Nz - 1));
		geo->z[k]    = -slab->tau[k];
	}
	for (int k = 0; k < Nz - 1; k++) geo->dz[k] = geo->z[k] - geo->z[k + 1];
	geo->dz[Nz - 1] = 0.0;
	for (int n = 0; n < Nspace; n++) slab->chi[n] = 1.0;

	geo->mux = malloc(NRAYS * sizeof(double));
	geo->muy = malloc(NRAYS * sizeof(double));
	geo->muz = malloc(NRAYS * sizeof(double));
	geo->wmu = malloc(NRAYS * sizeof(double));
	for (int mu = 0; mu < NRAYS; mu++) {
		geo->mux[mu] = ray_mux[mu];
		geo->muz[mu] = ray_muz[mu];
		geo->muy[mu] = sqrt(1.0 - SQ(ray_mux[mu]) - SQ(ray_muz[mu]));
		geo->wmu[mu] = 1.0 / NRAYS;
	}
	geo->Itop    = matrix_double(1, NX);
	geo->Ibottom = matrix_double(1, NX);

	fillMesh(geo);
	return slab;
}

static void set_source(Slab *slab, const double s[NPOLY])
{
	for (int k = 0; k < slab->Nz; k++) {
		double t = slab->tau[k];
		for (int l = 0; l < NX; l++)
			slab->S[k * NX + l] = s[0] + t * (s[1] + t * (s[2] + t * s[3]));
	}
}

/* --- Solve one ray; the irradiated lower boundary has one value per
	   column for all rays, so it is set per ray --   -------------- */

static void solve(Solver2D solver, Slab *slab, const double s[NPOLY], int mu,
				  bool_t to_obs, double *Psi)
{
	double Ib = I_exact(s, slab->tau[slab->Nz - 1], 0.0, ray_muz[mu], TRUE);

	for (int l = 0; l < NX; l++) {
		slab->geo.Itop[0][l]    = 0.0;
		slab->geo.Ibottom[0][l] = Ib;
	}
	solver(&slab->geo, 0, mu, to_obs, slab->chi, slab->S, slab->I, Psi);
}

/* --- Errors relative to the exact solution: emergent I+ (top
	   layer), the interior, and the two layers nearest either end;
	   and the largest difference between columns -- -------------- */

typedef struct {
	double emergent, interior, ends, columns;
} Errors;

static void accumulate(Errors *e, Slab *slab, const double s[NPOLY], int mu,
					   bool_t to_obs)
{
	int Nz = slab->Nz;

	for (int k = 0; k < Nz; k++) {
		double want = I_exact(s, slab->tau[k], slab->tau[0], ray_muz[mu], to_obs);

		for (int l = 0; l < NX; l++) {
			double got = slab->I[k * NX + l],
				   err = fabs(got - want) / MAX(fabs(want), 1.0e-300);

			if (k <= 1 || k >= Nz - 2)
				e->ends = MAX(e->ends, err);
			else
				e->interior = MAX(e->interior, err);
			if (to_obs && k == 0) e->emergent = MAX(e->emergent, err);
			e->columns = MAX(e->columns, fabs(got - slab->I[k * NX]) /
											 MAX(fabs(slab->I[k * NX]), 1.0e-300));
		}
	}
}

static Errors ray_errors(Solver2D solver, Slab *slab, const double s[NPOLY],
						 bool_t only_to_obs)
{
	Errors e = {0.0, 0.0, 0.0, 0.0};

	set_source(slab, s);
	for (int mu = 0; mu < NRAYS; mu++)
		for (int to_obs = only_to_obs ? 1 : 0; to_obs <= 1; to_obs++) {
			solve(solver, slab, s, mu, to_obs, NULL);
			accumulate(&e, slab, s, mu, to_obs);
		}
	return e;
}

/* --- The meshes: [r][m] for refinement r (Nz = 41 ... 321) and
	   m = 0 "layers" (DX = 1e8), 1 "columns" (DX = 0.05) -- ------- */

enum { NREFINE = 4,
	   NMESH   = 2 };

static Slab        *slabs[NREFINE][NMESH];
static const char  *mesh_name[NMESH] = {"layers", "columns"};
static const double mesh_dx[NMESH]   = {1.0e8, 0.05};

/* --- Bounds, about twice the measured values (gcc, x86-64), per
	   solver and mesh --                              -------------- */

typedef struct {
	const char *name;
	Solver2D    solver;
	double      exact_up[NMESH], exact_down[NMESH];
	double      order_emergent[NMESH], err_emergent[NMESH], err_interior[NMESH];
	double      columns;
} Method;

/* --- Measured (gcc, x86-64), "layers" / "columns":
	   LINEAR   linear S up 1.3e-15 / 2.6e-15, down 1.3e-15 / 7.6e-4;
			    order 2.00; Nz = 161: emergent 9.7e-4, interior 3.2e-3
			    on both meshes (the same as in 1-D)
	   BEZIER3  linear S up 2.8e-13 / 1.0e-13, down 4.0e-13 / 7.6e-4;
			    order 1.99; emergent 4.8e-4, interior 1.4e-3 on both:
			    second order where the 1-D solver is fourth, because the
			    upwind control point uses the secant slope (N-24)
	   Columns of the uniform slab agree to 8.0e-3 on "columns" (in
	   the two layers nearest the bottom), exactly on "layers" -- --- */

/* clang-format off */
static const Method methods[] = {
	{"LINEAR", Piecewise_Linear_2D,
	 {1e-14, 1e-14}, {1e-14, 1.6e-3}, {1.7, 1.7}, {2e-3, 2e-3}, {7e-3, 7e-3}, 1.6e-2},
	{"BEZIER3", Piecewise_Bezier3_2D,
	 {1e-12, 1e-12}, {1e-12, 1.6e-3}, {1.7, 1.7}, {1e-3, 1e-3}, {3e-3, 3e-3}, 1.6e-2}};
/* clang-format on */

#define N_METHODS (int) (sizeof(methods) / sizeof(methods[0]))

static void check_method(const Method *m, bool_t verbose)
{
	const double lin[NPOLY] = {1.0, 2.0, 0.0, 0.0}, cub[NPOLY] = {1.0, 1.0, 0.5, 0.1};

	for (int mesh = 0; mesh < NMESH; mesh++) {
		const char *mn = mesh_name[mesh];
		Errors      up = ray_errors(m->solver, slabs[1][mesh], lin, TRUE), all;
		Errors      e[NREFINE];
		double      down = 0.0, order = HUGE_VAL, columns;

		/* --- Linear S: exact upward, and downward on "layers" -- -- */

		all  = ray_errors(m->solver, slabs[1][mesh], lin, FALSE);
		down = MAX(all.interior, all.ends);
		if (verbose)
			printf("%-8s %-8s linear S: up %.2e, down %.2e\n", m->name, mn,
				   MAX(up.interior, up.ends), down);
		CHECK(MAX(up.interior, up.ends) <= m->exact_up[mesh],
			  "%s, %s, linear S, upward: error %.3g > %.3g", m->name, mn,
			  MAX(up.interior, up.ends), m->exact_up[mesh]);
		CHECK(down <= m->exact_down[mesh],
			  "%s, %s, linear S, both directions: error %.3g > %.3g", m->name,
			  mn, down, m->exact_down[mesh]);

		/* --- Cubic S: convergence --                   -------------- */

		columns = 0.0;
		for (int r = 0; r < NREFINE; r++) {
			e[r]    = ray_errors(m->solver, slabs[r][mesh], cub, FALSE);
			columns = MAX(columns, e[r].columns);
			if (verbose)
				printf("  %-8s %-8s Nz %4d  emergent %.3e  interior %.3e  ends %.3e\n",
					   m->name, mn, slabs[r][mesh]->Nz, e[r].emergent,
					   e[r].interior, e[r].ends);
			if (r > 0) order = MIN(order, log2(e[r - 1].emergent / e[r].emergent));
		}
		if (verbose)
			printf("%-8s %-8s order emergent %.2f; columns agree to %.2e\n",
				   m->name, mn, order, columns);

		CHECK(order >= m->order_emergent[mesh], "%s, %s: order %.2f < %.2f",
			  m->name, mn, order, m->order_emergent[mesh]);
		CHECK(e[2].emergent <= m->err_emergent[mesh],
			  "%s, %s: emergent error %.3g > %.3g", m->name, mn, e[2].emergent,
			  m->err_emergent[mesh]);
		CHECK(e[2].interior <= m->err_interior[mesh],
			  "%s, %s: interior error %.3g > %.3g", m->name, mn, e[2].interior,
			  m->err_interior[mesh]);
		CHECK(columns <= m->columns, "%s, %s: columns differ by %.3g > %.3g",
			  m->name, mn, columns, m->columns);
	}
}

/* --- Psi of the linear solver against the response to a unit S at
	   (k, l) with zero incident radiation (its I is linear in S). On
	   the "layers" mesh Psi is that response exactly. On "columns" the
	   periodic row couples S(k, l) back into the upwind intensity of
	   (k, l) through the column-line interpolation, which Psi leaves
	   out by design: reported, not checked --          -------------- */

static double psi_error(Slab *slab)
{
	int    Nspace = NX * slab->Nz;
	double err    = 0.0;

	for (int mu = 0; mu < NRAYS; mu++)
		for (int to_obs = 0; to_obs <= 1; to_obs++) {
			for (int n = 0; n < Nspace; n++) slab->S[n] = 1.0 + slab->tau[n / NX];
			for (int l = 0; l < NX; l++)
				slab->geo.Itop[0][l] = slab->geo.Ibottom[0][l] = 0.0;
			Piecewise_Linear_2D(&slab->geo, 0, mu, to_obs, slab->chi, slab->S,
								slab->I, slab->Psi);

			for (int n = 0; n < Nspace; n += 7) {
				memset(slab->S1, 0, Nspace * sizeof(double));
				slab->S1[n] = 1.0;
				Piecewise_Linear_2D(&slab->geo, 0, mu, to_obs, slab->chi,
									slab->S1, slab->I, NULL);
				if (slab->Psi[n] == 0.0 && slab->I[n] == 0.0) continue;
				err = MAX(err, fabs(slab->Psi[n] - slab->I[n]) /
								   MAX(fabs(slab->I[n]), 1.0e-300));
			}
		}
	return err;
}

static void check_psi_linear(double bound, bool_t verbose)
{
	double layers = psi_error(slabs[0][0]), columns = psi_error(slabs[0][1]);

	if (verbose)
		printf("LINEAR   Psi vs response to unit S: layers %.2e, columns %.2e\n",
			   layers, columns);
	CHECK(layers <= bound, "LINEAR, layers: Psi differs from dI/dS by %.3g > %.3g",
		  layers, bound);
}

int main(int argc, char *argv[])
{
	bool_t verbose = (argc > 1 && strcmp(argv[1], "-v") == 0);

	/* --- fillMesh reports the number of long characteristics -- */

	commandline.logfile = stdout;
	commandline.quiet   = !verbose;

	for (int r = 0; r < NREFINE; r++)
		for (int m = 0; m < NMESH; m++)
			slabs[r][m] = slab_new(40 * (1 << r) + 1, mesh_dx[m]);

	for (int n = 0; n < N_METHODS; n++) check_method(&methods[n], verbose);
	check_psi_linear(1.0e-14, verbose); /* layers: 2.2e-16 */

	return check_summary("test_formal_2d");
}
