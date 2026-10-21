/* ------- file: -------------------------- test_formal_3d.c --------

	   3-D formal solvers (ShortChar: linear and cubic Bezier short
	   characteristics, bilinear and bicubic horizontal interpolation)
	   on a horizontally uniform, periodic slab, where every column must
	   reproduce the 1-D analytic intensities along the ray's vertical
	   direction cosine muz (slab_exact.h).

	   As in test_formal_2d.c, two meshes: "layers" (columns so wide
	   that every upwind point lies on the layer above or below) and
	   "columns" (narrow columns: deep layers are crossed at vertical
	   planes, with long characteristics back to the next layer).

	   Checks: all columns agree; I+ is exact for linear S; convergence
	   order on a cubic S. -v prints the measured errors.
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

extern InputData   input;
extern CommandLine commandline;
extern Atmosphere  atmos;
extern Spectrum    spectrum;

/* --- Rays (mux, muy, muz) in different octants and a vertical one - */

#define NRAYS 5
static const double ray_mu[NRAYS][3] = {{0.6, 0.3, 0.7416198487095663},
										{-0.5, 0.5, 0.7071067811865476},
										{0.7, -0.6, 0.3872983346207417},
										{-0.2, -0.3, 0.9327379053088815},
										{0.0, 0.0, 1.0}};

#define NX 4
#define NY 4

/* --- A periodic slab of Nz layers, logarithmic in optical depth
	   from 1e-4 to 1e2, unit opacity (height = -tau), columns DX
	   apart in x and y. Meshes stay allocated (and reachable) to the
	   end --                                          -------------- */

typedef struct {
	Geometry geo;
	int      Nz;
	double  *tau, *chi, *S, *I;
} Slab;

static Slab *slab_new(int Nz, double DX)
{
	Slab     *slab   = calloc(1, sizeof(Slab));
	Geometry *geo    = &slab->geo;
	int       Nplane = NX * NY, Nspace = Nplane * Nz;

	slab->Nz  = Nz;
	slab->tau = malloc(Nz * sizeof(double));
	slab->chi = malloc(Nspace * sizeof(double));
	slab->S   = malloc(Nspace * sizeof(double));
	slab->I   = malloc(Nspace * sizeof(double));

	geo->x_boundary = geo->y_boundary = PERIODIC;
	geo->z_boundary_value[TOP]        = ZERO;
	geo->z_boundary_value[BOTTOM]     = IRRADIATED;
	geo->Nx                           = NX;
	geo->Ny                           = NY;
	geo->Nz                           = Nz;
	geo->Nplane                       = Nplane;
	geo->Nrays                        = NRAYS;
	geo->dx = geo->dy = DX;

	geo->z = malloc(Nz * sizeof(double));
	for (int k = 0; k < Nz; k++) {
		slab->tau[k] = 1.0e-4 * pow(1.0e6, (double) k / (Nz - 1));
		geo->z[k]    = -slab->tau[k];
	}
	for (int n = 0; n < Nspace; n++) slab->chi[n] = 1.0;

	geo->mux = malloc(NRAYS * sizeof(double));
	geo->muy = malloc(NRAYS * sizeof(double));
	geo->muz = malloc(NRAYS * sizeof(double));
	geo->wmu = malloc(NRAYS * sizeof(double));
	for (int mu = 0; mu < NRAYS; mu++) {
		geo->mux[mu] = ray_mu[mu][0];
		geo->muy[mu] = ray_mu[mu][1];
		geo->muz[mu] = ray_mu[mu][2];
		geo->wmu[mu] = 1.0 / NRAYS;
	}
	geo->Itop    = matrix_double(1, Nplane);
	geo->Ibottom = matrix_double(1, Nplane);

	fillMesh(geo);
	return slab;
}

static void set_source(Slab *slab, const double s[NPOLY])
{
	int Nplane = NX * NY;

	for (int k = 0; k < slab->Nz; k++) {
		double t = slab->tau[k];
		for (int n = 0; n < Nplane; n++)
			slab->S[k * Nplane + n] = s[0] + t * (s[1] + t * (s[2] + t * s[3]));
	}
}

/* --- Errors relative to the exact solution: emergent I+ (top
	   layer), the interior, the two layers nearest either end, and
	   the largest difference between columns --      -------------- */

typedef struct {
	double emergent, interior, ends, columns;
} Errors;

static Errors ray_errors(Slab *slab, const double s[NPOLY], bool_t only_to_obs)
{
	int    Nz = slab->Nz, Nplane = NX * NY;
	Errors e = {0.0, 0.0, 0.0, 0.0};

	set_source(slab, s);
	for (int mu = 0; mu < NRAYS; mu++)
		for (int to_obs = only_to_obs ? 1 : 0; to_obs <= 1; to_obs++) {
			double muz = ray_mu[mu][2],
				   Ib  = I_exact(s, slab->tau[Nz - 1], 0.0, muz, TRUE);

			for (int n = 0; n < Nplane; n++) {
				slab->geo.Itop[0][n]    = 0.0;
				slab->geo.Ibottom[0][n] = Ib;
			}
			ShortChar(&slab->geo, 0, mu, to_obs, slab->chi, slab->S, slab->I,
					  NULL);

			for (int k = 0; k < Nz; k++) {
				double want = I_exact(s, slab->tau[k], slab->tau[0], muz, to_obs);

				for (int n = 0; n < Nplane; n++) {
					double got = slab->I[k * Nplane + n],
						   ref = slab->I[k * Nplane],
						   err = fabs(got - want) / MAX(fabs(want), 1.0e-300);

					if (k <= 1 || k >= Nz - 2)
						e.ends = MAX(e.ends, err);
					else
						e.interior = MAX(e.interior, err);
					if (to_obs && k == 0) e.emergent = MAX(e.emergent, err);
					e.columns = MAX(e.columns, fabs(got - ref) / MAX(fabs(ref), 1.0e-300));
				}
			}
		}
	return e;
}

/* --- Meshes [i][r][m]: horizontal interpolation i (the stencils are
	   built for it by fillMesh), refinement r (Nz = 41 ... 321, 7 to 53
	   points per decade), m = 0 "layers" (DX = 1e8), 1 "columns"
	   (DX = 0.05) --                                  -------------- */

enum { NINTERP = 2,
	   NREFINE = 4,
	   NMESH   = 2 };

static Slab        *slabs[NINTERP][NREFINE][NMESH];
static const char  *mesh_name[NMESH] = {"layers", "columns"};
static const double mesh_dx[NMESH]   = {1.0e8, 0.05};

/* --- Bounds, about twice the measured values (gcc, x86-64), per
	   solver, horizontal interpolation and mesh --    -------------- */

typedef struct {
	const char     *name;
	enum S_interpol S_interpolation;
	enum order_3D   interpolate_3D;
	double          exact_up[NMESH], exact_down[NMESH], order_emergent[NMESH];
	double          err_emergent[NMESH], err_interior[NMESH], columns;
} Method;

/* --- Measured (gcc, x86-64), "layers" / "columns", Nz = 161:
	   linear S, upward 1.7e-15 / 2.0e-15 (LINEAR), 1.2e-13 / 1.0e-13
	   (BEZIER3); both directions exact on "layers", 1.5e-3 on
	   "columns" (I- is not linear in height). On "layers" the same as
	   2-D and, for LINEAR, 1-D (BEZIER3: order 2.82, emergent 2.5e-5,
	   interior 1.2e-4 with the quadratic Bezier upwind slope of N-24;
	   order 1.99, emergent 4.8e-4 with the secant). On
	   "columns" the vertical interpolation between two layers and the
	   long characteristics limit the order to ~1.5 at these
	   resolutions: emergent 3.8e-3 to 3.9e-3, interior 9.5e-3 to
	   9.8e-3; columns of the uniform slab agree to 3.0e-2 at Nz = 41,
	   5.6e-3 at 161 --                                -------------- */

/* clang-format off */
static const Method methods[] = {
	{"LINEAR/bilinear", S_LINEAR, LINEAR_3D,
	 {1e-14, 1e-14}, {1e-14, 3.1e-3}, {1.7, 1.2}, {2e-3, 8e-3}, {7e-3, 2e-2}, 6.1e-2},
	{"LINEAR/bicubic", S_LINEAR, BICUBIC_3D,
	 {1e-14, 1e-14}, {1e-14, 3.1e-3}, {1.7, 1.2}, {2e-3, 8e-3}, {7e-3, 2e-2}, 6.1e-2},
	{"BEZIER3/bilinear", S_BEZIER3, LINEAR_3D,
	 {1e-12, 1e-12}, {1e-12, 3.1e-3}, {2.6, 1.2}, {5e-5, 8e-3}, {2.5e-4, 2e-2}, 6.1e-2},
	{"BEZIER3/bicubic", S_BEZIER3, BICUBIC_3D,
	 {1e-12, 1e-12}, {1e-12, 3.1e-3}, {2.6, 1.2}, {5e-5, 8e-3}, {2.5e-4, 2e-2}, 6.1e-2}};
/* clang-format on */

#define N_METHODS (int) (sizeof(methods) / sizeof(methods[0]))

static void check_method(const Method *m, bool_t verbose)
{
	const double lin[NPOLY] = {1.0, 2.0, 0.0, 0.0}, cub[NPOLY] = {1.0, 1.0, 0.5, 0.1};

	input.S_interpolation = m->S_interpolation;
	input.interpolate_3D  = m->interpolate_3D;

	for (int mesh = 0; mesh < NMESH; mesh++) {
		Slab      **sl = &slabs[m->interpolate_3D][0][0];
		const char *mn = mesh_name[mesh];
		Errors      up = ray_errors(sl[1 * NMESH + mesh], lin, TRUE), e[NREFINE],
			   all     = ray_errors(sl[1 * NMESH + mesh], lin, FALSE);
		double exact = MAX(up.interior, up.ends), down = MAX(all.interior, all.ends),
			   order = HUGE_VAL, columns = 0.0;

		for (int r = 0; r < NREFINE; r++) {
			e[r]    = ray_errors(sl[r * NMESH + mesh], cub, FALSE);
			columns = MAX(columns, e[r].columns);
			if (verbose)
				printf("  %-16s %-8s Nz %4d  emergent %.3e  interior %.3e  ends %.3e  columns %.2e\n",
					   m->name, mn, sl[r * NMESH + mesh]->Nz, e[r].emergent,
					   e[r].interior, e[r].ends, e[r].columns);
			if (r > 0) order = MIN(order, log2(e[r - 1].emergent / e[r].emergent));
		}
		if (verbose)
			printf("%-16s %-8s linear S up %.2e, down %.2e; order emergent %.2f; columns agree to %.2e\n",
				   m->name, mn, exact, down, order, columns);

		CHECK(exact <= m->exact_up[mesh], "%s, %s, linear S, upward: error %.3g > %.3g",
			  m->name, mn, exact, m->exact_up[mesh]);
		CHECK(down <= m->exact_down[mesh], "%s, %s, linear S, both directions: error %.3g > %.3g",
			  m->name, mn, down, m->exact_down[mesh]);
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

/* --- Polarized (ShortChar_Stokes, DELO-Bezier), as in
	   test_formal_2d: constant K' = kscale K_ELEMENTS, source vector
	   K (S, 0, 0, 0) --                               -------------- */

static double **stokes_S[NINTERP][NMESH], **stokes_I[NINTERP][NMESH];

static void stokes_set(Slab *slab, double **Sv, const double s[NPOLY], double kscale)
{
	int        Nspace = NX * NY * slab->Nz;
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

static void set_boundary(Slab *slab, const double s[NPOLY], int mu)
{
	double Ib = I_exact(s, slab->tau[slab->Nz - 1], 0.0, ray_mu[mu][2], TRUE);

	for (int n = 0; n < NX * NY; n++) {
		slab->geo.Itop[0][n]    = 0.0;
		slab->geo.Ibottom[0][n] = Ib;
	}
}

/* --- Measured (gcc, x86-64), both interpolations and meshes, Nz =
	   81: K' = 0 equal to the unpolarized Bezier3 solver to 4.6e-16;
	   linear S, I+ (tau <= 10, all columns, any Stokes component
	   relative to I) exact to 1.2e-13 --              -------------- */

static void check_stokes(bool_t verbose)
{
	const double lin[NPOLY] = {1.0, 2.0, 0.0, 0.0}, cub[NPOLY] = {1.0, 1.0, 0.5, 0.1};
	int          Nplane = NX * NY;

	spectrum.as                      = calloc(1, sizeof(ActiveSet));
	atmos.backgrflags                = calloc(1, sizeof(flags));
	atmos.backgrflags[0].ispolarized = TRUE;
	input.magneto_optical            = TRUE;
	input.S_interpolation            = S_BEZIER3;
	input.S_interpolation_stokes     = DELO_BEZIER3;

	for (int i = 0; i < NINTERP; i++) {
		input.interpolate_3D = (i == 0) ? LINEAR_3D : BICUBIC_3D;
		for (int mesh = 0; mesh < NMESH; mesh++) {
			Slab    *slab   = slabs[i][1][mesh];
			int      Nspace = Nplane * slab->Nz;
			double **Sv = stokes_S[i][mesh] = matrix_double(4, Nspace),
				   **I = stokes_I[i][mesh] = matrix_double(4, Nspace);
			double unpol = 0.0, quv = 0.0, lin_err = 0.0;
			Mat    K;

			/* --- K' = 0: the unpolarized solution --     -------------- */

			stokes_set(slab, Sv, cub, 0.0);
			for (int mu = 0; mu < NRAYS; mu++)
				for (int to_obs = 0; to_obs <= 1; to_obs++) {
					set_boundary(slab, cub, mu);
					ShortChar_Stokes(&slab->geo, 0, mu, to_obs, slab->chi, Sv, I, NULL);
					ShortChar(&slab->geo, 0, mu, to_obs, slab->chi, slab->S, slab->I, NULL);
					for (int n = 0; n < Nspace; n++) {
						unpol = MAX(unpol, fabs(I[0][n] - slab->I[n]) / fabs(slab->I[n]));
						for (int j = 1; j < 4; j++) quv = MAX(quv, fabs(I[j][n]));
					}
				}

			/* --- Linear S: I+ exact --                   -------------- */

			stokes_set(slab, Sv, lin, 1.0);
			full_K(K_ELEMENTS, 1.0, K);
			for (int mu = 0; mu < NRAYS; mu++) {
				set_boundary(slab, lin, mu);
				ShortChar_Stokes(&slab->geo, 0, mu, TRUE, slab->chi, Sv, I, NULL);
				for (int k = 0; k < slab->Nz && slab->tau[k] <= 10.0; k++) {
					Vec want;
					I_exact_pol(lin, K, slab->tau[k], slab->tau[0], ray_mu[mu][2], TRUE, want);
					for (int n = 0; n < Nplane; n++)
						for (int j = 0; j < 4; j++)
							lin_err = MAX(lin_err, fabs(I[j][k * Nplane + n] - want[j]) / want[0]);
				}
			}
			if (verbose)
				printf("STOKES %-7s %-8s K' = 0 vs unpolarized %.2e (QUV %.1e); linear S up %.2e\n",
					   i ? "bicubic" : "linear", mesh_name[mesh], unpol, quv, lin_err);
			CHECK(unpol <= 1.0e-13 && quv == 0.0,
				  "STOKES, %s: K' = 0 differs from unpolarized by %.3g (QUV %.3g)",
				  mesh_name[mesh], unpol, quv);
			CHECK(lin_err <= 1.0e-12, "STOKES, %s: linear S, I+ error %.3g",
				  mesh_name[mesh], lin_err);
		}
	}
}

int main(int argc, char *argv[])
{
	bool_t verbose = (argc > 1 && strcmp(argv[1], "-v") == 0);

	commandline.logfile = stdout;
	commandline.quiet   = !verbose;

	for (int i = 0; i < NINTERP; i++) {
		input.interpolate_3D = (i == 0) ? LINEAR_3D : BICUBIC_3D;
		for (int r = 0; r < NREFINE; r++)
			for (int m = 0; m < NMESH; m++)
				slabs[i][r][m] = slab_new(40 * (1 << r) + 1, mesh_dx[m]);
	}

	for (int n = 0; n < N_METHODS; n++) check_method(&methods[n], verbose);
	check_stokes(verbose);

	return check_summary("test_formal_3d");
}
