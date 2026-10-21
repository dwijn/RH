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
	   the linear solver against the response to a unit S.

	   Polarized (Piece_Stokes_Bezier3_2D) with a constant propagation
	   matrix and an unpolarized source vector (stokes_exact.h): with
	   K' = 0 it must equal the unpolarized Bezier3 solver; I+ is exact
	   for linear S; the emergent Stokes vector converges for cubic S.
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
#include "check.h"
#include "stokes_exact.h"

extern CommandLine commandline;
extern Atmosphere  atmos;
extern Spectrum    spectrum;
extern InputData   input;

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

static Slab *slab_new(int Nz, double DX, bool_t uneven)
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
		geo->dx[l] = uneven ? DX * ((l % 2) ? 1.4 : 0.6) : DX;
		geo->x[l]  = (l == 0) ? 0.0 : geo->x[l - 1] + geo->dx[l - 1];
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

/* --- The meshes: [r][m] for refinement r (Nz = 41 ... 321, 7 to 53
	   points per decade) and m = 0 "layers" (DX = 1e8), 1 "columns"
	   (DX = 0.05), 2 "uneven" (columns alternately 0.03 and 0.07
	   wide) --                                        -------------- */

enum { NREFINE = 4,
	   NMESH   = 3 };

static Slab        *slabs[NREFINE][NMESH];
static const char  *mesh_name[NMESH]   = {"layers", "columns", "uneven"};
static const double mesh_dx[NMESH]     = {1.0e8, 0.05, 0.05};
static const bool_t mesh_uneven[NMESH] = {FALSE, FALSE, TRUE};

/* --- Bounds, about twice the measured values (gcc, x86-64), per
	   solver and mesh --                              -------------- */

typedef struct {
	const char *name;
	Solver2D    solver;
	double      exact_up[NMESH], exact_down[NMESH];
	double      order_emergent[NMESH], err_emergent[NMESH], err_interior[NMESH];
	double      columns;
} Method;

/* --- Measured (gcc, x86-64), "layers" / "columns" ("uneven" as
	   "columns", down 7.1e-4, columns agree to 7.5e-3):
	   LINEAR   linear S up 1.3e-15 / 2.6e-15, down 1.3e-15 / 7.6e-4;
				order 2.00; Nz = 161: emergent 9.7e-4, interior 3.2e-3
				on both meshes (the same as in 1-D)
	   BEZIER3  linear S up 2.8e-13 / 1.0e-13, down 4.0e-13 / 7.6e-4;
				"layers": order 2.82, emergent 2.5e-5, interior 1.2e-4;
				"columns": order 2.58, emergent 1.0e-4, interior 4.6e-4
				("uneven" 2.58, 1.1e-4, 4.3e-4). The upwind control
				point uses the slope of the quadratic Bezier curve,
				2 x secant - local derivative (N-24; with the secant
				alone: order 1.99, emergent 4.8e-4 on all meshes). On
				"columns" the vertical interpolation of I limits the
				order
	   Columns of the uniform slab agree to 8.0e-3 on "columns" (in
	   the two layers nearest the bottom), exactly on "layers" -- --- */

/* clang-format off */
static const Method methods[] = {
	{"LINEAR", Piecewise_Linear_2D,
	 {1e-14, 1e-14, 1e-14}, {1e-14, 1.6e-3, 1.5e-3}, {1.7, 1.7, 1.7},
	 {2e-3, 2e-3, 2e-3}, {7e-3, 7e-3, 7e-3}, 1.6e-2},
	{"BEZIER3", Piecewise_Bezier3_2D,
	 {1e-12, 1e-12, 1e-12}, {1e-12, 1.6e-3, 1.5e-3}, {2.6, 2.4, 2.4},
	 {5e-5, 2e-4, 2e-4}, {2.5e-4, 1e-3, 1e-3}, 1.6e-2}};
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

/* --- Polarized: constant K' = kscale K_ELEMENTS in every point,
	   source vector K (S, 0, 0, 0). Arrays for the Stokes vectors are
	   per slab and stay allocated --                  -------------- */

static double **stokes_S[NREFINE][NMESH], **stokes_I[NREFINE][NMESH];

static void stokes_set(Slab *slab, double **Sv, const double s[NPOLY], double kscale)
{
	int        Nspace = NX * slab->Nz;
	ActiveSet *as     = &spectrum.as[0];
	Mat        K;

	atmos.Nspace = Nspace;
	as->chi_c    = realloc(as->chi_c, 4 * Nspace * sizeof(double));
	as->chip_c   = realloc(as->chip_c, 3 * Nspace * sizeof(double));
	full_K(K_ELEMENTS, kscale, K);
	set_source(slab, s);
	for (int n = 0; n < Nspace; n++) {
		as->chi_c[n]               = 1.0;
		as->chi_c[Nspace + n]      = K[0][1];
		as->chi_c[2 * Nspace + n]  = K[0][2];
		as->chi_c[3 * Nspace + n]  = K[0][3];
		as->chip_c[n]              = K[2][3]; /* rho_Q */
		as->chip_c[Nspace + n]     = K[3][1]; /* rho_U */
		as->chip_c[2 * Nspace + n] = K[1][2]; /* rho_V */
		for (int i = 0; i < 4; i++) Sv[i][n] = K[i][0] * slab->S[n];
	}
}

/* --- Largest error of any Stokes component relative to I, over the
	   layers with tau <= 10 (deeper, the boundary condition Q = U = V
	   = 0 at the bottom differs from the semi-infinite solution) and
	   all columns; and of the emergent vector --     -------------- */

static void stokes_errors(Slab *slab, double **I, const double s[NPOLY], int mu,
						  bool_t to_obs, double *interior, double *emergent)
{
	Mat K;

	full_K(K_ELEMENTS, 1.0, K);
	for (int k = (to_obs ? 0 : 1); k < slab->Nz && slab->tau[k] <= 10.0; k++) {
		Vec want;

		I_exact_pol(s, K, slab->tau[k], slab->tau[0], ray_muz[mu], to_obs, want);
		for (int l = 0; l < NX; l++)
			for (int n = 0; n < 4; n++) {
				double err = fabs(I[n][k * NX + l] - want[n]) / fabs(want[0]);
				if (k == 0)
					*emergent = MAX(*emergent, err);
				else
					*interior = MAX(*interior, err);
			}
	}
}

static void stokes_solve(Slab *slab, double **Sv, double **I, const double s[NPOLY],
						 int mu, bool_t to_obs)
{
	double Ib = I_exact(s, slab->tau[slab->Nz - 1], 0.0, ray_muz[mu], TRUE);

	for (int l = 0; l < NX; l++) {
		slab->geo.Itop[0][l]    = 0.0;
		slab->geo.Ibottom[0][l] = Ib;
	}
	Piece_Stokes_Bezier3_2D(&slab->geo, 0, mu, to_obs, slab->chi, Sv, I, NULL);
}

/* --- Measured (gcc, x86-64): K' = 0 equal to the unpolarized solver
	   to 4.6e-16 on all meshes; linear S, I+ exact to 2.9e-13 (layers)
	   / 1.1e-13 (before the fix of B-26 the emergent vector was
	   computed with the source vector of a point in the layer below:
	   error 1.8e-9, and 38% at K' = 0 in the bottom layer of I-);
	   cubic S on "layers": order 2.80, emergent 3.3e-5 at Nz = 161,
	   as unpolarized (N-24; with the secant upwind slope order 1.96,
	   emergent 5.8e-4) --                            -------------- */

static void check_stokes(bool_t verbose)
{
	const double lin[NPOLY] = {1.0, 2.0, 0.0, 0.0}, cub[NPOLY] = {1.0, 1.0, 0.5, 0.1};

	spectrum.as                      = calloc(1, sizeof(ActiveSet));
	atmos.backgrflags                = calloc(1, sizeof(flags));
	atmos.backgrflags[0].ispolarized = TRUE;
	input.magneto_optical            = TRUE;

	for (int r = 0; r < NREFINE; r++)
		for (int m = 0; m < NMESH; m++) {
			stokes_S[r][m] = matrix_double(4, NX * slabs[r][m]->Nz);
			stokes_I[r][m] = matrix_double(4, NX * slabs[r][m]->Nz);
		}

	for (int mesh = 0; mesh < NMESH; mesh++) {
		const char *mn   = mesh_name[mesh];
		Slab       *slab = slabs[1][mesh];
		double    **Sv = stokes_S[1][mesh], **I = stokes_I[1][mesh];
		double      unpol = 0.0, quv = 0.0, lin_err = 0.0, lin_em = 0.0;

		/* --- K' = 0: the unpolarized solution, Q = U = V = 0 -- --- */

		stokes_set(slab, Sv, cub, 0.0);
		for (int mu = 0; mu < NRAYS; mu++)
			for (int to_obs = 0; to_obs <= 1; to_obs++) {
				stokes_solve(slab, Sv, I, cub, mu, to_obs);
				solve(Piecewise_Bezier3_2D, slab, cub, mu, to_obs, NULL);
				for (int n = 0; n < NX * slab->Nz; n++) {
					unpol = MAX(unpol, fabs(I[0][n] - slab->I[n]) / fabs(slab->I[n]));
					for (int i = 1; i < 4; i++) quv = MAX(quv, fabs(I[i][n]));
				}
			}

		/* --- Linear S: I+ exact --                     -------------- */

		stokes_set(slab, Sv, lin, 1.0);
		for (int mu = 0; mu < NRAYS; mu++) {
			stokes_solve(slab, Sv, I, lin, mu, TRUE);
			stokes_errors(slab, I, lin, mu, TRUE, &lin_err, &lin_em);
		}
		if (verbose)
			printf("STOKES   %-8s K' = 0 vs unpolarized %.2e (QUV %.1e); linear S up "
				   "%.2e, emergent %.2e\n",
				   mn, unpol, quv, lin_err, lin_em);
		CHECK(unpol <= 1.0e-13 && quv == 0.0,
			  "STOKES, %s: K' = 0 differs from unpolarized by %.3g (QUV %.3g)", mn,
			  unpol, quv);
		CHECK(MAX(lin_err, lin_em) <= 1.0e-12, "STOKES, %s: linear S, I+ error %.3g",
			  mn, MAX(lin_err, lin_em));
	}

	/* --- Cubic S on "layers": emergent Stokes vector converges -- */

	double em[NREFINE], order = HUGE_VAL;
	for (int r = 0; r < NREFINE; r++) {
		double interior = 0.0;
		em[r]           = 0.0;
		stokes_set(slabs[r][0], stokes_S[r][0], cub, 1.0);
		for (int mu = 0; mu < NRAYS; mu++) {
			stokes_solve(slabs[r][0], stokes_S[r][0], stokes_I[r][0], cub, mu, TRUE);
			stokes_errors(slabs[r][0], stokes_I[r][0], cub, mu, TRUE, &interior, &em[r]);
		}
		if (verbose)
			printf("  STOKES   layers   Nz %4d  emergent %.3e  interior %.3e\n",
				   slabs[r][0]->Nz, em[r], interior);
		if (r > 0) order = MIN(order, log2(em[r - 1] / em[r]));
	}
	if (verbose) printf("STOKES   layers   order emergent %.2f\n", order);
	CHECK(order >= 2.6, "STOKES, layers: order %.2f < 2.6", order);
	CHECK(em[2] <= 7.0e-5, "STOKES, layers: emergent error %.3g > 7e-5", em[2]);
}

int main(int argc, char *argv[])
{
	bool_t verbose = (argc > 1 && strcmp(argv[1], "-v") == 0);

	/* --- fillMesh reports the number of long characteristics -- */

	commandline.logfile = stdout;
	commandline.quiet   = !verbose;

	for (int r = 0; r < NREFINE; r++)
		for (int m = 0; m < NMESH; m++)
			slabs[r][m] = slab_new(40 * (1 << r) + 1, mesh_dx[m], mesh_uneven[m]);

	for (int n = 0; n < N_METHODS; n++) check_method(&methods[n], verbose);
	check_psi_linear(1.0e-14, verbose); /* layers: 2.2e-16 */
	check_stokes(verbose);

	return check_summary("test_formal_2d");
}
