/* ------- file: -------------------------- test_longchar_3d.c ------

	   3-D formal solution with long characteristics (N-39, B-52, B-53)
	   on a periodic slab whose source function varies horizontally:

		 S(x, y, tau) = s0 + s1 tau + A sin(2 pi (x + y)),

	   constant opacity chi, period 1 in x and y, so that the optical
	   depth is tau = chi * depth. The vertical part is linear, which the
	   vertical interpolation and every integration reproduce; the errors
	   come from the horizontal structure. Layers are Q dx apart, so the
	   oblique rays (Q |mux| / muz or Q |muy| / muz > 1) cross vertical
	   planes before the next layer and use long characteristics;
	   refining dx and dz together keeps the same crossings. chi = 1
	   gives optically thin sections along the long characteristics
	   (0.03 to 0.2 at N = 32), chi = 30 thick ones.

	   Exact intensity along direction d (d = (mux, muy, muz) for the
	   upward sweep, -(mux, muy, muz) for the downward one), phase
	   phi = 2 pi (x + y), kappa = 2 pi (d_x + d_y) / chi per unit
	   optical path, optical path L to the boundary:

		 I = I_lin + A Im[exp(i phi) (1 - exp(-(1 + i kappa) L)) / (1 + i kappa)],

	   with I- = 0 at the top and the bottom irradiated with the
	   semi-infinite value (L -> infinity), which is then exact for I+.

	   Checks, for linear and Bezier short characteristics, both
	   horizontal interpolations and both long-characteristic schemes
	   (S_INTERPOLATION_LONGCHAR): the errors stay below bounds and
	   converge, and the two long-characteristic schemes agree (they
	   carry the intensity only to the first point of each row, the
	   others interpolate it in the vertical planes). -v prints the
	   errors.
	   --                                              -------------- */

#include <complex.h>
#undef I /* the meshes have a member I */
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
#include "bezier.h"
#include "slab_exact.h"

extern InputData   input;
extern CommandLine commandline;

#define NRAYS 5
static const double ray_mu[NRAYS][3] = {{0.6, 0.3, 0.7416198487095663},
										{-0.5, 0.5, 0.7071067811865476},
										{0.7, -0.6, 0.3872983346207417},
										{-0.2, -0.3, 0.9327379053088815},
										{0.0, 0.0, 1.0}};

#define Q       4.0 /* layer spacing / dx */
#define DEPTH   4.0 /* geometric depth of the slab */
#define AMP     0.5
#define NREFINE 3
#define NCHI    2

static const double s_lin[NPOLY]    = {1.0, 0.5, 0.0, 0.0};
static const int    mesh_N[NREFINE] = {8, 16, 32};
static const double chi_value[NCHI] = {1.0, 30.0};

typedef struct {
	Geometry geo;
	int      N, Nz;
	double   dx, chi0, *tau, *chi, *S, *I;
} Mesh;

static Mesh *mesh_new(int N, double chi0)
{
	Mesh     *ms  = calloc(1, sizeof(Mesh));
	Geometry *geo = &ms->geo;
	double    dx  = 1.0 / N;
	int       Nz = (int) lround(DEPTH / (Q * dx)) + 1, Nplane = N * N,
		Nspace = Nplane * Nz;

	ms->N    = N;
	ms->Nz   = Nz;
	ms->dx   = dx;
	ms->chi0 = chi0;
	ms->tau  = malloc(Nz * sizeof(double));
	ms->chi  = malloc(Nspace * sizeof(double));
	ms->S    = malloc(Nspace * sizeof(double));
	ms->I    = malloc(Nspace * sizeof(double));

	geo->x_boundary = geo->y_boundary = PERIODIC;
	geo->z_boundary_value[TOP]        = ZERO;
	geo->z_boundary_value[BOTTOM]     = IRRADIATED;
	geo->Nx                           = N;
	geo->Ny                           = N;
	geo->Nz                           = Nz;
	geo->Nplane                       = Nplane;
	geo->Nrays                        = NRAYS;
	geo->dx = geo->dy = dx;

	geo->z = malloc(Nz * sizeof(double));
	for (int k = 0; k < Nz; k++) {
		geo->z[k]  = -k * Q * dx;
		ms->tau[k] = -chi0 * geo->z[k];
	}
	for (int k = 0; k < Nz; k++)
		for (int m = 0; m < N; m++)
			for (int l = 0; l < N; l++) {
				int n      = k * Nplane + m * N + l;
				ms->chi[n] = chi0;
				ms->S[n]   = s_lin[0] + s_lin[1] * ms->tau[k] +
						   AMP * sin(2.0 * M_PI * (l + m) * dx);
			}

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
	return ms;
}

static void mesh_free(Mesh *ms)
{
	/* --- fillMesh sets the stencil_uw of layers 1 .. Nz - 1 and the
		   stencil_dw of layers 0 .. Nz - 2; the others are not
		   initialized (code review) --                  -------------- */

	Geometry *geo = &ms->geo;

	for (int mu = 0; mu < NRAYS; mu++) {
		for (int k = 0; k < ms->Nz; k++) {
			Stencil *st[2] = {(k > 0) ? &geo->stencil_uw[mu][k] : NULL,
							  (k < ms->Nz - 1) ? &geo->stencil_dw[mu][k] : NULL};

			for (int w = 0; w < 2; w++)
				if (st[w] && st[w]->longchar) {
					free(st[w]->longchar->stencil);
					free(st[w]->longchar);
				}
		}
		free(geo->stencil_uw[mu]);
		free(geo->stencil_dw[mu]);
	}
	free(geo->stencil_uw);
	free(geo->stencil_dw);
	freeMatrix((void **) geo->Itop);
	freeMatrix((void **) geo->Ibottom);
	free(geo->z);
	free(geo->mux);
	free(geo->muy);
	free(geo->muz);
	free(geo->wmu);
	free(ms->tau);
	free(ms->chi);
	free(ms->S);
	free(ms->I);
	free(ms);
}

/* --- Horizontal part of the exact intensity, optical path L -- --- */

static double I_horizontal(double phi, double kappa, double L)
{
	double complex z = 1.0 + _Complex_I * kappa, f;

	f = isinf(L) ? 1.0 / z : (1.0 - cexp(-z * L)) / z;
	return AMP * cimag(cexp(_Complex_I * phi) * f);
}

/* --- Largest error relative to S0 = s0 (the scale of S at the top,
	   where I- is small, so that a relative error in I- would only
	   measure how small it is): I+ (upward) in the interior and at
	   the top, I- (downward) in the interior --       -------------- */

typedef struct {
	double up, emergent, down;
} Errors;

static Errors mesh_errors(Mesh *ms, double *Isave)
{
	/* --- Isave, if not NULL, receives I of every ray and direction
		   (2 NRAYS Nspace) --                            -------------- */

	int    N = ms->N, Nz = ms->Nz, Nplane = N * N;
	Errors e = {0.0, 0.0, 0.0};

	for (int mu = 0; mu < NRAYS; mu++)
		for (int to_obs = 0; to_obs <= 1; to_obs++) {
			double sgn = to_obs ? 1.0 : -1.0, muz = ray_mu[mu][2],
				   kappa = 2.0 * M_PI * sgn * (ray_mu[mu][0] + ray_mu[mu][1]) / ms->chi0;

			for (int n = 0; n < Nplane; n++) {
				double phi = 2.0 * M_PI * (n % N + n / N) * ms->dx;

				ms->geo.Itop[0][n]    = 0.0;
				ms->geo.Ibottom[0][n] = I_exact(s_lin, ms->tau[Nz - 1], 0.0, muz, TRUE) +
										I_horizontal(phi, kappa, INFINITY);
			}
			ShortChar(&ms->geo, 0, mu, to_obs, ms->chi, ms->S, ms->I, NULL);
			if (Isave)
				memcpy(Isave + (2 * mu + to_obs) * (size_t) (Nz * Nplane), ms->I,
					   Nz * Nplane * sizeof(double));

			for (int k = 0; k < Nz - 1; k++)
				for (int n = 0; n < Nplane; n++) {
					double phi  = 2.0 * M_PI * (n % N + n / N) * ms->dx,
						   L    = to_obs ? INFINITY : ms->tau[k] / muz,
						   want = I_exact(s_lin, ms->tau[k], 0.0, muz, to_obs) +
								  I_horizontal(phi, kappa, L),
						   err = fabs(ms->I[k * Nplane + n] - want) / s_lin[0];

					if (!to_obs)
						e.down = MAX(e.down, err);
					else if (k == 0)
						e.emergent = MAX(e.emergent, err);
					else
						e.up = MAX(e.up, err);
				}
		}
	return e;
}

/* --- Measured (gcc, x86-64), N = 8, 16, 32, both horizontal
	   interpolations. chi = 1: I+ and emergent 1.7e-2 to 1.8e-2 at
	   N = 32, order 0.85 to 0.95 (the diffusion of the repeated
	   interpolation of I: h^2 per layer, 1 / dz layers per unit
	   depth); chi = 30: I+ 4.9e-3 (linear short characteristics),
	   1.4e-3 to 1.6e-3 (Bezier), emergent 4.6e-3 to 4.9e-3, order 1.4
	   to 1.7. I- is worst in the first layer below the top: in one
	   optically thick layer under I = 0 it rises to about S within
	   tau ~ 1, which the linear interpolation of I between the two
	   layers in the vertical planes cannot follow (6.7e-2 at chi = 1,
	   8.2e-2 at chi = 30, N = 32), the same for both long-characteristic
	   schemes. The intensities of the two long-characteristic schemes
	   differ, point by point, by at most 1.3e-2, 2.8e-3, 4.1e-4 of S0
	   at chi = 1 (N = 8, 16, 32) and 5.3e-5, 4.7e-4, 6.4e-4 at
	   chi = 30 --                                     -------------- */

typedef struct {
	double up_max[NCHI], down_max[NCHI], order_min, lc_agree[NCHI][NREFINE];
} Bounds;

/* --- One section, LongChar_Linear and LongChar_Bezier3 directly, on
	   sections of very different optical thickness (code review): the
	   intensity stays between the smallest and the largest of I_uw,
	   S_uw, S_loc, S_dw, and does not decrease when S_loc grows. The
	   weight of S_dw, which in the last section is the local point, is
	   0 for the linear step and is checked. For the Bezier step it
	   enters through the Steffen derivative at the local end and is
	   not bounded here: small (-0.07 at worst for constant chi and
	   S_dw away from S_loc) but down to about -3 where S_loc ~ S_dw
	   and the downwind section is thin, where the limiter sets the
	   derivative to twice the downwind slope. The parabolic step of
	   N-39 had about -1/dtau_dw there for any S --     -------------- */

static bool_t verbose_section = FALSE;

static double section(int scheme, const double *v)
{
	/* --- v: ds_uw, ds_dw, chi_uw, chi_loc, chi_dw, S_uw, S_loc, S_dw,
		   I_uw --                                       -------------- */

	return (scheme == S_BEZIER3)
			   ? LongChar_Bezier3(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8])
			   : LongChar_Linear(v[0], v[2], v[3], v[5], v[6], v[8]);
}

static void check_section(void)
{
	const double dtau[]    = {1.0e-6, 1.0e-3, 0.1, 1.0, 3.0, 30.0, 1.0e3};
	const int    Nd        = sizeof(dtau) / sizeof(dtau[0]);
	unsigned     seed      = 12345;
	int          nrange[2] = {0, 0}, nloc[2] = {0, 0}, ndw[2] = {0, 0}, ntot = 0;
	double       wdw_min[2] = {0.0, 0.0};

	for (int iu = 0; iu < Nd; iu++)
		for (int id = 0; id < Nd; id++)
			for (int t = 0; t < 40; t++) {
				double v[9], lo, hi;

				/* --- chi varies by up to 10 along the section; ds sets
					   dtau to about the target --            ---------- */

				for (int j = 2; j < 5; j++) {
					seed = seed * 1103515245u + 12345u;
					v[j] = pow(10.0, ((seed >> 8) & 0xffff) / 65535.0);
				}
				v[0] = dtau[iu] / v[3];
				v[1] = dtau[id] / v[3];
				for (int j = 5; j < 9; j++) {
					seed = seed * 1103515245u + 12345u;
					v[j] = ((seed >> 8) & 0xffff) / 65535.0;
				}
				lo = MIN(MIN(v[5], v[6]), MIN(v[7], v[8]));
				hi = MAX(MAX(v[5], v[6]), MAX(v[7], v[8]));
				ntot++;
				for (int sc = 0; sc < 2; sc++) {
					double I0 = section(sc, v), w[9], h = 1.0e-6, tol = 1.0e-12, wdw;

					if (I0 < lo - tol || I0 > hi + tol) nrange[sc]++;
					memcpy(w, v, sizeof(w));
					w[6] += h;
					if (section(sc, w) < I0 - tol) nloc[sc]++;
					memcpy(w, v, sizeof(w));
					w[7] += h;
					wdw         = (section(sc, w) - I0) / h;
					wdw_min[sc] = MIN(wdw_min[sc], wdw);
					if (sc == S_LINEAR && wdw != 0.0) ndw[sc]++;
				}
			}
	for (int sc = 0; sc < 2; sc++) {
		const char *name = sc ? "LongChar_Bezier3" : "LongChar_Linear";

		if (verbose_section)
			printf("%s: smallest weight of S_dw %.3g\n", name, wdw_min[sc]);
		CHECK(nrange[sc] == 0, "%s: %d of %d sections outside the range of I_uw and S", name,
			  nrange[sc], ntot);
		CHECK(nloc[sc] == 0, "%s: %d of %d sections decrease with S_loc", name, nloc[sc], ntot);
		if (sc == S_LINEAR)
			CHECK(ndw[sc] == 0, "%s: %d of %d sections depend on S_dw", name, ndw[sc], ntot);
	}
}

static const Bounds bounds = {{3.0e-2, 1.0e-2},
							  {0.1, 0.15},
							  0.75,
							  {{2.0e-2, 4.0e-3, 6.0e-4}, {1.0e-4, 7.0e-4, 1.0e-3}}};

int main(int argc, char *argv[])
{
	bool_t verbose = (argc > 1 && strcmp(argv[1], "-v") == 0);

	const enum S_interpol sc[2]    = {S_LINEAR, S_BEZIER3};
	const char           *sname[2] = {"LINEAR", "BEZIER3"};

	commandline.logfile = stdout;
	commandline.quiet   = TRUE;

	for (int c = 0; c < NCHI; c++)
		for (int i = 0; i < 2; i++) {
			Mesh       *ms[NREFINE];
			const char *iname = i ? "bicubic" : "linear";

			input.interpolate_3D = (i == 0) ? LINEAR_3D : BICUBIC_3D;
			for (int r = 0; r < NREFINE; r++) ms[r] = mesh_new(mesh_N[r], chi_value[c]);

			for (int a = 0; a < 2; a++) {
				Errors  e[2][NREFINE];
				double *Is[2][NREFINE];

				input.S_interpolation = sc[a];
				for (int b = 0; b < 2; b++) {
					input.S_interpolation_longchar = sc[b];
					for (int r = 0; r < NREFINE; r++) {
						Is[b][r] = malloc(2 * NRAYS * (size_t) (ms[r]->Nz * ms[r]->N * ms[r]->N) *
										  sizeof(double));
						e[b][r]  = mesh_errors(ms[r], Is[b][r]);
						if (verbose)
							printf("  chi %4.1f %-7s SC %-7s LC %-7s N %2d  I+ %.3e  emergent %.3e  I- %.3e\n",
								   chi_value[c], iname, sname[a], sname[b], mesh_N[r],
								   e[b][r].up, e[b][r].emergent, e[b][r].down);
					}
					for (int r = 1; r < NREFINE; r++) {
						double order = log2(e[b][r - 1].emergent / e[b][r].emergent);

						CHECK(order >= bounds.order_min,
							  "chi %g, %s, SC %s, LC %s, N %d: emergent order %.2f < %.2f",
							  chi_value[c], iname, sname[a], sname[b], mesh_N[r], order,
							  bounds.order_min);
					}
					{
						Errors *f = &e[b][NREFINE - 1];
						double  m = MAX(f->up, f->emergent);

						CHECK(m <= bounds.up_max[c], "chi %g, %s, SC %s, LC %s: I+ error %.3g > %.3g",
							  chi_value[c], iname, sname[a], sname[b], m, bounds.up_max[c]);
						CHECK(f->down <= bounds.down_max[c],
							  "chi %g, %s, SC %s, LC %s: I- error %.3g > %.3g", chi_value[c],
							  iname, sname[a], sname[b], f->down, bounds.down_max[c]);
					}
				}

				/* --- Linear and Bezier long characteristics give nearly
					   the same intensities, point by point, relative to
					   S0 (code review: not only similar errors) -- ---- */

				for (int r = 0; r < NREFINE; r++) {
					size_t n_all = 2 * NRAYS * (size_t) (ms[r]->Nz * ms[r]->N * ms[r]->N);
					double d     = 0.0;

					for (size_t n = 0; n < n_all; n++)
						d = MAX(d, fabs(Is[0][r][n] - Is[1][r][n]));
					d /= s_lin[0];
					if (verbose)
						printf("chi %4.1f %-7s SC %-7s N %2d: LC LINEAR and BEZIER3 differ by %.1e\n",
							   chi_value[c], iname, sname[a], mesh_N[r], d);
					CHECK(d <= bounds.lc_agree[c][r],
						  "chi %g, %s, SC %s, N %d: long-characteristic schemes differ by %.3g",
						  chi_value[c], iname, sname[a], mesh_N[r], d);
					free(Is[0][r]);
					free(Is[1][r]);
				}
			}
			for (int r = 0; r < NREFINE; r++) mesh_free(ms[r]);
		}
	verbose_section = verbose;
	check_section();
	return check_summary("test_longchar_3d");
}
/* ------- end ---------------------------- test_longchar_3d.c ------ */
