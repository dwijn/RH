/* ------- file: -------------------------- test_zeeman.c -----------

	   Zeeman splitting (zeeman.c) and the polarized absorption and
	   magneto-optical profiles (profile.c, with rhf1d/project.c for the
	   angles) of fictional lines, against expressions derived from
	   first principles (docs/unit-tests, section "Zeeman effect").

	   Line configurations (lower -> upper level):
	   -- 1S0 -> 1P1, the normal triplet (g = 1), and the same line with
		  an explicit effective Lande factor (three components);
	   -- 2S1/2 -> 2P1/2 and 2S1/2 -> 2P3/2 (anomalous, as Na D1, D2);
	   -- 5P2 -> 5D2 (dJ = 0, no pi component from M = 0 to M = 0) and
		  5P1 -> 5D0 (upper J smaller than lower);
	   -- 5D0 -> 5F1: g(5F1) = 0, no splitting, no polarization at any
		  field strength;
	   -- 1S0 -> 3P0: J = 0 -> 0, which zeemanPolarizable() must reject
		  (B-27: its Zeeman strengths were normalized by zero).

	   Checks:
	   -- zeemanPolarizable(): every line except J = 0 -> 0;
	   -- Zeeman(): components, shifts g_l M_l - g_u M_u (Larmor units),
		  strengths 3 (Ju Jl 1; -Mu Ml -q)^2 normalized per q, g from LS
		  coupling, the effective Lande factor against its formula and
		  against the barycentre of the q = +1 components;
	   -- Profile(): phi_I, phi_Q, phi_U, phi_V and psi_Q, psi_U, psi_V
		  for a set of field strengths and directions, two rays and both
		  directions along them, against the propagation matrix derived
		  from the polarization of the circular and linear normal modes,
		  in the IAU convention, with RH's sign of Stokes V applied (RH:
		  positive V is negative helicity, opposite to IAU 1973);
	   -- Milne-Eddington: RH's profiles at one point turned into the
		  opacity and emissivity of a slab with constant field, constant
		  line-to-continuum ratio and S linear in continuum optical
		  depth, solved with the polarized 1-D solver: the emergent
		  Stokes vector across the line against I = S0 + mu S1 K^-1 e_I
		  with K from the derived profiles (Unno-Rachkovsky), and in a
		  weak field against V = -g_eff vB cos(theta) dI/dv (RH's sign).
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
#include "constant.h"
#include "check.h"
#include "stokes_exact.h"

extern Atmosphere  atmos;
extern Geometry    geometry;
extern InputData   input;
extern CommandLine commandline;
extern Spectrum    spectrum;

void Bproject(void);
void Piece_Stokes_Bezier3_1D(int nspect, int mu, bool_t to_obs, double *chi,
							 double **S, double **I, double *Psi);

/* --- Line configurations: labels (parsed by determinate() for S
	   and L), statistical weights (giving J), and the quantum numbers
	   the test uses independently --                  -------------- */

typedef struct {
	const char *name, *label_l, *label_u;
	double      S_l, L_l, J_l, S_u, L_u, J_u;
	double      g_eff_explicit; /* 0: from the labels */
} Config;

/* clang-format off */
static const Config configs[] = {
	{"1S0 -> 1P1",           "TEST 3S 1SE", "TEST 3P 1PO", 0.0, 0, 0.0, 0.0, 1, 1.0, 0.0},
	{"triplet, g_eff = 1.5", "TEST 3S 1SE", "TEST 3P 1PO", 0.0, 0, 0.0, 0.0, 1, 1.0, 1.5},
	{"2S1/2 -> 2P1/2",       "TEST 3S 2SE", "TEST 3P 2PO", 0.5, 0, 0.5, 0.5, 1, 0.5, 0.0},
	{"2S1/2 -> 2P3/2",       "TEST 3S 2SE", "TEST 3P 2PO", 0.5, 0, 0.5, 0.5, 1, 1.5, 0.0},
	{"5P2 -> 5D2",           "TEST 4P 5PE", "TEST 4D 5DO", 2.0, 1, 2.0, 2.0, 2, 2.0, 0.0},
	{"5P1 -> 5D0",           "TEST 4P 5PE", "TEST 4D 5DO", 2.0, 1, 1.0, 2.0, 2, 0.0, 0.0},
	{"5D0 -> 5F1 (null)",    "TEST 4D 5DE", "TEST 4F 5FO", 2.0, 2, 0.0, 2.0, 3, 1.0, 0.0},
	{"1S0 -> 3P0 (J 0-0)",   "TEST 3S 1SE", "TEST 3P 3PO", 0.0, 0, 0.0, 1.0, 1, 0.0, 0.0}};
/* clang-format on */

#define N_CONFIG (int) (sizeof(configs) / sizeof(configs[0]))

/* --- Lande factor in LS coupling --                 -------------- */

static double lande_LS(double S, double L, double J)
{
	if (J == 0.0) return 0.0; /* no M != 0 sublevels: irrelevant */
	return 1.0 + (J * (J + 1.0) + S * (S + 1.0) - L * (L + 1.0)) / (2.0 * J * (J + 1.0));
}

/* --- The expected components: all (Mu, Ml) with |Ml - Mu| <= 1, in
	   the loop order of Zeeman() (Ml outer, Mu inner) -- ----------- */

#define MAX_COMP 64

typedef struct {
	int    N, q[MAX_COMP];
	double shift[MAX_COMP], strength[MAX_COMP];
} Pattern;

static void expected_pattern(const Config *c, Pattern *p)
{
	double norm[3] = {0.0, 0.0, 0.0};
	double g_l = lande_LS(c->S_l, c->L_l, c->J_l), g_u = lande_LS(c->S_u, c->L_u, c->J_u);

	p->N = 0;
	if (c->g_eff_explicit != 0.0) {
		for (int q = -1; q <= 1; q++) {
			p->q[p->N]        = q;
			p->shift[p->N]    = q * c->g_eff_explicit;
			p->strength[p->N] = 1.0;
			p->N++;
		}
		return;
	}
	for (double Ml = -c->J_l; Ml <= c->J_l; Ml++)
		for (double Mu = -c->J_u; Mu <= c->J_u; Mu++) {
			if (fabs(Mu - Ml) > 1.0) continue;
			int q             = (int) (Ml - Mu);
			p->q[p->N]        = q;
			p->shift[p->N]    = g_l * Ml - g_u * Mu;
			p->strength[p->N] = 3.0 * SQ(w3js(c->J_u, c->J_l, 1.0, -Mu, Ml, -q));
			norm[q + 1] += p->strength[p->N];
			p->N++;
		}
	for (int n = 0; n < p->N; n++)
		p->strength[n] = (norm[p->q[n] + 1] > 0.0) ? p->strength[n] / norm[p->q[n] + 1] : 0.0;
}

/* --- A fictional atom with one line --             -------------- */

#define LAMBDA0 500.0 /* nm */
#define NLAMBDA 81

static Atom       atom;
static AtomicLine line;
static char      *labels[2];
static double     weights[2];

static void make_line(const Config *c)
{
	memset(&line, 0, sizeof(line));
	labels[0]  = strdup(c->label_l);
	labels[1]  = strdup(c->label_u);
	weights[0] = 2.0 * c->J_l + 1.0;
	weights[1] = 2.0 * c->J_u + 1.0;
	atom.label = labels;
	atom.g     = weights;
	strcpy(atom.ID, "TE");
	line.atom          = &atom;
	line.i             = 0;
	line.j             = 1;
	line.lambda0       = LAMBDA0;
	line.g_Lande_eff   = c->g_eff_explicit;
	line.polarizable   = TRUE;
	line.symmetric     = FALSE;
	line.Voigt         = TRUE;
	line.Grad          = 4.0e9; /* a ~ 0.05 */
	line.Ncomponent    = 1;
	line.c_shift       = calloc(1, sizeof(double));
	line.c_fraction    = malloc(sizeof(double));
	line.c_fraction[0] = 1.0;
	line.Nlambda       = NLAMBDA;
	line.lambda        = malloc(NLAMBDA * sizeof(double));
	for (int la = 0; la < NLAMBDA; la++)
		line.lambda[la] = LAMBDA0 + 0.04 * (2.0 * la / (NLAMBDA - 1) - 1.0);
}

static void free_line(void)
{
	if (line.zm) {
		free(line.zm->q);
		free(line.zm->strength);
		free(line.zm->shift);
		free(line.zm);
	}
	if (line.phi) {
		freeMatrix((void **) line.phi);
		freeMatrix((void **) line.phi_Q);
		freeMatrix((void **) line.phi_U);
		freeMatrix((void **) line.phi_V);
		freeMatrix((void **) line.psi_Q);
		freeMatrix((void **) line.psi_U);
		freeMatrix((void **) line.psi_V);
	}
	free(line.wphi);
	free(line.c_shift);
	free(line.c_fraction);
	free(line.lambda);
	free(labels[0]);
	free(labels[1]);
}

/* --- Zeeman(): pattern and effective Lande factor -- ------------ */

/* --- Lines with an explicit g_eff, or |Ju - Jl| <= 1 except
	   J = 0 -> 0 --                                   -------------- */

static bool_t expect_polarizable(const Config *c)
{
	return c->g_eff_explicit != 0.0 ||
		   (fabs(c->J_u - c->J_l) <= 1.0 && !(c->J_u == 0.0 && c->J_l == 0.0));
}

static void check_pattern(const Config *c, bool_t verbose)
{
	Pattern p;
	double  err_s = 0.0, err_sh = 0.0;

	make_line(c);
	bool_t pol = zeemanPolarizable(&atom, 0, 1, c->g_eff_explicit);
	CHECK(pol == expect_polarizable(c), "%s: zeemanPolarizable() %d, expected %d",
		  c->name, pol, expect_polarizable(c));
	if (!pol) {
		if (verbose) printf("%-22s not polarizable\n", c->name);
		free_line();
		return;
	}
	expected_pattern(c, &p);
	Zeeman(&line);

	CHECK(line.zm->Ncomponent == p.N, "%s: %d components, expected %d", c->name,
		  line.zm->Ncomponent, p.N);
	if (line.zm->Ncomponent == p.N) {
		for (int n = 0; n < p.N; n++) {
			CHECK(line.zm->q[n] == p.q[n], "%s: component %d has q = %d, expected %d",
				  c->name, n, line.zm->q[n], p.q[n]);
			err_sh = MAX(err_sh, fabs(line.zm->shift[n] - p.shift[n]));
			err_s  = isfinite(line.zm->strength[n])
						 ? MAX(err_s, fabs(line.zm->strength[n] - p.strength[n]))
						 : HUGE_VAL;
		}
	}

	/* --- Effective Lande factor: the formula, and the barycentre of
		   the q = +1 components (where there are any) -- ----------- */

	double g_l = lande_LS(c->S_l, c->L_l, c->J_l), g_u = lande_LS(c->S_u, c->L_u, c->J_u);
	double g_eff = (c->g_eff_explicit != 0.0)
					   ? c->g_eff_explicit
					   : 0.5 * (g_u + g_l) + 0.25 * (g_u - g_l) *
												 (c->J_u * (c->J_u + 1.0) - c->J_l * (c->J_l + 1.0));
	double bary = 0.0, wsum = 0.0;
	for (int n = 0; n < p.N; n++)
		if (p.q[n] == 1) {
			bary += p.strength[n] * p.shift[n];
			wsum += p.strength[n];
		}

	if (verbose)
		printf("%-22s %2d components, strengths %.1e, shifts %.1e, g_eff %.4f "
			   "(formula %.4f, barycentre %.4f)\n",
			   c->name, line.zm->Ncomponent, err_s, err_sh, line.zm->g_eff, g_eff,
			   (wsum > 0.0) ? bary : 0.0);
	CHECK(err_s <= 1.0e-14, "%s: strengths off by %.3g", c->name, err_s);
	CHECK(err_sh <= 1.0e-14, "%s: shifts off by %.3g", c->name, err_sh);
	CHECK(fabs(line.zm->g_eff - g_eff) <= 1.0e-14, "%s: g_eff %.15g, expected %.15g",
		  c->name, line.zm->g_eff, g_eff);
	if (wsum > 0.0)
		CHECK(fabs(bary - g_eff) <= 1.0e-14,
			  "%s: barycentre of the q = +1 components %.15g, g_eff %.15g", c->name,
			  bary, g_eff);
	free_line();
}

/* --- Atmosphere points: field strength [T], inclination and azimuth
	   with respect to the vertical and x, Doppler width -- --------- */

typedef struct {
	double B, gamma, chi, vbroad;
} Point;

/* clang-format off */
static Point points[] = {
	{0.0,  0.0,  0.0,  3.0e3},
	{0.01, 0.3,  0.4,  3.0e3},
	{0.1,  1.2,  2.0,  2.0e3},
	{0.5,  2.5, -0.7,  4.0e3},  /* pointing away from the observer */
	{0.3,  M_PI / 2.0, 0.25, 3.0e3},
	{0.2,  0.0,  0.0,  3.0e3},
	{0.2,  0.0,  0.0,  3.0e3}}; /* set to point along the oblique ray */
/* clang-format on */

#define NPOINT (int) (sizeof(points) / sizeof(points[0]))

#define NRAYS 2
static double ray_mux[NRAYS] = {0.0, 0.5}, ray_muy[NRAYS] = {0.0, 0.3},
			  ray_muz[NRAYS] = {1.0, 0.8124038404635961}, ray_wmu[NRAYS] = {0.5, 0.5};

static void setup_atmosphere(void)
{
	points[NPOINT - 1].gamma = acos(ray_muz[1]);
	points[NPOINT - 1].chi   = atan2(ray_muy[1], ray_mux[1]);

	atmos.Nspace   = NPOINT;
	atmos.Nrays    = NRAYS;
	atmos.Stokes   = TRUE;
	atmos.moving   = FALSE;
	atmos.wmu      = ray_wmu;
	atmos.B        = malloc(NPOINT * sizeof(double));
	atmos.gamma_B  = malloc(NPOINT * sizeof(double));
	atmos.chi_B    = malloc(NPOINT * sizeof(double));
	atom.vbroad    = malloc(NPOINT * sizeof(double));
	geometry.Nrays = NRAYS;
	geometry.mux   = ray_mux;
	geometry.muy   = ray_muy;
	geometry.muz   = ray_muz;
	geometry.wmu   = ray_wmu;
	for (int k = 0; k < NPOINT; k++) {
		atmos.B[k]       = points[k].B;
		atmos.gamma_B[k] = points[k].gamma;
		atmos.chi_B[k]   = points[k].chi;
		atom.vbroad[k]   = points[k].vbroad;
	}
	Bproject();

	input.StokesMode      = FULL_STOKES;
	input.magneto_optical = TRUE;
	input.limit_memory    = FALSE;
	input.Nthreads        = 1;
}

/* --- Expected profiles. For propagation direction n with the
	   polarization frame (e1, e2) (e1, e2, n right-handed; Q > 0 along
	   e1, U > 0 at 45 degrees from e1 toward e2) and theta, chi the
	   angles of B with n and of its projection with e1, the IAU
	   propagation matrix (dI/ds = -K I + ...) has

		 eta_I = [H0 sin^2 + (H+ + H-) (1 + cos^2) / 2] / 2
		 eta_Q = [H0 - (H+ + H-) / 2] sin^2 cos 2chi / 2   (U: sin 2chi)
		 eta_V = (H- - H+) cos / 2
		 rho_Q = -[F0 - (F+ + F-) / 2] sin^2 cos 2chi / 2  (U: sin 2chi)
		 rho_V = (F+ - F-) cos / 2

	   with Hq, Fq the strength-weighted sums of H and F over the
	   components q = Ml - Mu, v = (lambda - lambda0) / Doppler width:
	   the q = -1 (blue) components absorb positive helicity (IAU
	   V > 0), and F is the real part of the refractive index
	   (Kramers-Kronig). RH writes V with the opposite sign: eta_V,
	   rho_Q and rho_U change sign, rho_V does not --  -------------- */

typedef struct {
	double phi, phi_Q, phi_U, phi_V, psi_Q, psi_U, psi_V;
} Profiles;

static void expected_profiles(const Pattern *p, double v, double a, double vB,
							  double cos_t, double c2x, double s2x, Profiles *e)
{
	double H[3] = {0.0, 0.0, 0.0}, F[3] = {0.0, 0.0, 0.0};

	for (int n = 0; n < p->N; n++) {
		double Fn, Hn = Voigt(a, v - p->shift[n] * vB, &Fn);
		H[p->q[n] + 1] += p->strength[n] * Hn;
		F[p->q[n] + 1] += p->strength[n] * Fn;
	}
	double sin2 = 1.0 - cos_t * cos_t;
	double eI   = 0.5 * (H[1] * sin2 + 0.5 * (H[2] + H[0]) * (1.0 + cos_t * cos_t));
	double eQU  = 0.5 * (H[1] - 0.5 * (H[2] + H[0])) * sin2;
	double eV   = 0.5 * (H[0] - H[2]) * cos_t;
	double rQU  = -0.5 * (F[1] - 0.5 * (F[2] + F[0])) * sin2;
	double rV   = 0.5 * (F[2] - F[0]) * cos_t;

	/* --- RH's sign of V: D K D with D = diag(1, 1, 1, -1) -- ------ */

	e->phi   = eI;
	e->phi_Q = eQU * c2x;
	e->phi_U = eQU * s2x;
	e->phi_V = -eV;
	e->psi_Q = -rQU * c2x;
	e->psi_U = -rQU * s2x;
	e->psi_V = rV;
}

/* --- Angles of B for ray mu, direction to_obs. Upward (to the
	   observer): n = (mux, muy, muz), e1 the projection of the vertical
	   on the plane of the sky (x for a vertical ray), e2 = n x e1.
	   Downward: n' = -n with RH's frame e1' = e2, e2' = e1 -- ------ */

static void field_angles(int k, int mu, bool_t to_obs, double *cos_t, double *c2x,
						 double *s2x)
{
	double b[3] = {sin(points[k].gamma) * cos(points[k].chi),
				   sin(points[k].gamma) * sin(points[k].chi), cos(points[k].gamma)};
	double n[3] = {ray_mux[mu], ray_muy[mu], ray_muz[mu]}, e1[3], e2[3];

	if (ray_muz[mu] == 1.0) {
		e1[0] = 1.0, e1[1] = 0.0, e1[2] = 0.0;
	} else {
		double s = sqrt(1.0 - SQ(ray_muz[mu]));
		for (int i = 0; i < 3; i++) e1[i] = ((i == 2) - ray_muz[mu] * n[i]) / s;
	}
	e2[0] = n[1] * e1[2] - n[2] * e1[1];
	e2[1] = n[2] * e1[0] - n[0] * e1[2];
	e2[2] = n[0] * e1[1] - n[1] * e1[0];

	double bn = b[0] * n[0] + b[1] * n[1] + b[2] * n[2],
		   b1 = b[0] * e1[0] + b[1] * e1[1] + b[2] * e1[2],
		   b2 = b[0] * e2[0] + b[1] * e2[1] + b[2] * e2[2];
	if (!to_obs) {
		double t = b1;
		bn       = -bn;
		b1       = b2;
		b2       = t;
	}
	double bt2 = b1 * b1 + b2 * b2;

	*cos_t = bn;
	*c2x   = (bt2 > 0.0) ? (b1 * b1 - b2 * b2) / bt2 : 1.0;
	*s2x   = (bt2 > 0.0) ? 2.0 * b1 * b2 / bt2 : 0.0;
}

static void check_profiles(const Config *c, bool_t verbose)
{
	const double Larmor = (Q_ELECTRON / (4.0 * PI * M_ELECTRON)) * (LAMBDA0 * NM_TO_M);
	Pattern      p;
	double       err[7] = {0}, scale = 0.0, quv = 0.0;

	make_line(c);
	expected_pattern(c, &p);
	Profile(&line);

	for (int la = 0; la < NLAMBDA; la++)
		for (int mu = 0; mu < NRAYS; mu++)
			for (int to_obs = 0; to_obs <= 1; to_obs++) {
				int lamu = 2 * (NRAYS * la + mu) + to_obs;

				for (int k = 0; k < NPOINT; k++) {
					double vb = atom.vbroad[k], sv = 1.0 / (SQRTPI * vb),
						   v = (line.lambda[la] - LAMBDA0) * CLIGHT / (vb * LAMBDA0),
						   a = line.Grad * (NM_TO_M * LAMBDA0) / (4.0 * PI) / vb,
						   ct, c2x, s2x;
					Profiles e;

					field_angles(k, mu, to_obs, &ct, &c2x, &s2x);
					expected_profiles(&p, v, a, Larmor * points[k].B / vb, ct, c2x, s2x, &e);

					double got[7]  = {line.phi[lamu][k], line.phi_Q[lamu][k],
									  line.phi_U[lamu][k], line.phi_V[lamu][k],
									  line.psi_Q[lamu][k], line.psi_U[lamu][k],
									  line.psi_V[lamu][k]},
						   want[7] = {e.phi, e.phi_Q, e.phi_U, e.phi_V,
									  e.psi_Q, e.psi_U, e.psi_V};

					scale = MAX(scale, want[0] * sv);
					for (int i = 0; i < 7; i++) {
						double d = isfinite(got[i]) ? fabs(got[i] - want[i] * sv) : HUGE_VAL;
						err[i]   = MAX(err[i], d);
						if (i > 0) quv = MAX(quv, fabs(got[i]));
					}
				}
			}
	for (int i = 0; i < 7; i++) err[i] /= scale;

	if (verbose)
		printf("%-22s phi I %.1e Q %.1e U %.1e V %.1e  psi Q %.1e U %.1e V %.1e  "
			   "(max |QUV| %.2g of phi_I)\n",
			   c->name, err[0], err[1], err[2], err[3], err[4], err[5], err[6],
			   quv / scale);

	const char *what[7] = {"phi_I", "phi_Q", "phi_U", "phi_V", "psi_Q", "psi_U", "psi_V"};
	for (int i = 0; i < 7; i++)
		CHECK(err[i] <= 1.0e-13, "%s: %s off by %.3g (relative to max phi_I)", c->name,
			  what[i], err[i]);
	free_line();
}

/* --- Milne-Eddington slab. The profiles of point k, ray mu (to the
	   observer) are taken from Profile() for all wavelengths; then a
	   slab of NDEP points, continuum opacity 1 per unit height (height
	   = -tau_c), line opacity chi_l phi with chi_l = ETA0 sqrt(pi)
	   vbroad (line centre about ETA0 times the continuum), S = S0 +
	   S1 tau_c for line and continuum, is solved per wavelength. In
	   continuum optical depth the propagation matrix is
	   K = 1 + chi_l Phi, Phi the profile matrix in RH's convention, and
	   the emergent Stokes vector of the semi-infinite atmosphere is
	   S0 e_I + mu S1 K^-1 e_I --                       -------------- */

#define NDEP 41
#define ETA0 5.0
#define S0   1.0
#define S1   2.0

typedef struct {
	int    config, point, mu;
	double tol;
} MECase;

static double me_profile[NLAMBDA][7];

static void solve_me_slab(int mu, double chi_l, double Iem[NLAMBDA][4])
{
	ActiveSet *as  = &spectrum.as[0];
	double    *chi = malloc(NDEP * sizeof(double)), **Sv = matrix_double(4, NDEP),
		   **I = matrix_double(4, NDEP), *tau = malloc(NDEP * sizeof(double));

	atmos.Nspace  = NDEP;
	geometry.Ndep = NDEP;
	free(geometry.height);
	geometry.height = malloc(NDEP * sizeof(double));
	as->chi_c       = realloc(as->chi_c, 4 * NDEP * sizeof(double));
	as->chip_c      = realloc(as->chip_c, 3 * NDEP * sizeof(double));
	for (int k = 0; k < NDEP; k++) {
		tau[k]             = (k == 0) ? 0.0 : 1.0e-4 * pow(1.0e6, (double) (k - 1) / (NDEP - 2));
		geometry.height[k] = -tau[k];
	}
	for (int la = 0; la < NLAMBDA; la++) {
		const double *f = me_profile[la];

		for (int k = 0; k < NDEP; k++) {
			double S = S0 + S1 * tau[k];

			chi[k]                   = 1.0 + chi_l * f[0];
			as->chi_c[k]             = chi[k];
			as->chi_c[NDEP + k]      = chi_l * f[1];
			as->chi_c[2 * NDEP + k]  = chi_l * f[2];
			as->chi_c[3 * NDEP + k]  = chi_l * f[3];
			as->chip_c[k]            = chi_l * f[4];
			as->chip_c[NDEP + k]     = chi_l * f[5];
			as->chip_c[2 * NDEP + k] = chi_l * f[6];
			Sv[0][k]                 = S;
			for (int n = 1; n < 4; n++) Sv[n][k] = S * chi_l * f[n] / chi[k];
		}
		geometry.Ibottom[0][mu] = S0 + S1 * tau[NDEP - 1] + ray_muz[mu] * S1;
		Piece_Stokes_Bezier3_1D(0, mu, TRUE, chi, Sv, I, NULL);
		for (int n = 0; n < 4; n++) Iem[la][n] = I[n][0];
	}
	free(chi);
	free(tau);
	freeMatrix((void **) Sv);
	freeMatrix((void **) I);
}

static void check_milne_eddington(const MECase *mc, bool_t verbose)
{
	const Config *c      = &configs[mc->config];
	const double  Larmor = (Q_ELECTRON / (4.0 * PI * M_ELECTRON)) * (LAMBDA0 * NM_TO_M);
	int           k = mc->point, mu = mc->mu;
	double        vb = points[k].vbroad, sv = 1.0 / (SQRTPI * vb), chi_l = ETA0 / sv,
		   a = 4.0e9 * (NM_TO_M * LAMBDA0) / (4.0 * PI) / vb, vB = Larmor * points[k].B / vb;
	double  Iem[NLAMBDA][4], err                                 = 0.0, ct, c2x, s2x;
	Pattern p;

	/* --- RH's profiles of point k, ray mu, upward --  -------------- */

	atmos.Nspace = NPOINT;
	make_line(c);
	expected_pattern(c, &p);
	Profile(&line);
	for (int la = 0; la < NLAMBDA; la++) {
		int lamu          = 2 * (NRAYS * la + mu) + 1;
		me_profile[la][0] = line.phi[lamu][k];
		me_profile[la][1] = line.phi_Q[lamu][k];
		me_profile[la][2] = line.phi_U[lamu][k];
		me_profile[la][3] = line.phi_V[lamu][k];
		me_profile[la][4] = line.psi_Q[lamu][k];
		me_profile[la][5] = line.psi_U[lamu][k];
		me_profile[la][6] = line.psi_V[lamu][k];
	}
	solve_me_slab(mu, chi_l, Iem);

	/* --- Expected: K from the derived profiles --    -------------- */

	field_angles(k, mu, TRUE, &ct, &c2x, &s2x);
	double v[NLAMBDA], Iexp[NLAMBDA][4];
	for (int la = 0; la < NLAMBDA; la++) {
		Profiles e;
		Mat      K, Kinv;

		v[la] = (line.lambda[la] - LAMBDA0) * CLIGHT / (vb * LAMBDA0);
		expected_profiles(&p, v[la], a, vB, ct, c2x, s2x, &e);
		double eI = ETA0 * e.phi, eQ = ETA0 * e.phi_Q, eU = ETA0 * e.phi_U,
			   eV = ETA0 * e.phi_V, rQ = ETA0 * e.psi_Q, rU = ETA0 * e.psi_U,
			   rV = ETA0 * e.psi_V;
		Mat M     = {{1.0 + eI, eQ, eU, eV},
					 {eQ, 1.0 + eI, rV, -rU},
					 {eU, -rV, 1.0 + eI, rQ},
					 {eV, rU, -rQ, 1.0 + eI}};
		memcpy(K, M, sizeof(Mat));
		mat_inv(K, Kinv);
		for (int n = 0; n < 4; n++) {
			Iexp[la][n] = ray_muz[mu] * S1 * Kinv[n][0] + ((n == 0) ? S0 : 0.0);
			err         = MAX(err, fabs(Iem[la][n] - Iexp[la][n]) / Iexp[la][0]);
		}
	}

	/* --- Weak field: to first order in vB, V = -g_eff vB cos(theta)
		   dI0/dv in RH's sign, I0 the field-free emergent intensity
		   (here mu S1 / (1 + ETA0 phi) + S0), differentiated with a
		   small step --                                -------------- */

	double wfa = 0.0, vmax = 0.0;
	if (vB < 0.05) {
		for (int la = 0; la < NLAMBDA; la++) {
			double Ipm[2];

			for (int sgn = 0; sgn < 2; sgn++) {
				Profiles e0;
				expected_profiles(&p, v[la] + (sgn ? 1e-4 : -1e-4), a, 0.0, ct, c2x, s2x, &e0);
				Ipm[sgn] = S0 + ray_muz[mu] * S1 / (1.0 + ETA0 * e0.phi);
			}
			double Vw = -line.zm->g_eff * vB * ct * (Ipm[1] - Ipm[0]) / 2.0e-4;
			wfa       = MAX(wfa, fabs(Iem[la][3] - Vw));
			vmax      = MAX(vmax, fabs(Iem[la][3]));
		}
		wfa /= vmax;
	}
	if (verbose)
		printf("ME %-22s point %d ray %d: Stokes vector %.1e of I%s", c->name, k, mu, err,
			   (vB < 0.05) ? "" : "\n");
	if (vB < 0.05 && verbose)
		printf("; weak field: max |V| %.2e, V - WFA %.1e of max |V|\n", vmax, wfa);
	CHECK(err <= mc->tol, "ME %s, point %d, ray %d: Stokes vector off by %.3g", c->name,
		  k, mu, err);
	if (vB < 0.05)
		CHECK(wfa <= 2.0e-3, "ME %s, point %d: V differs from the weak-field limit by %.3g",
			  c->name, k, wfa);
	free_line();
}

int main(int argc, char *argv[])
{
	bool_t verbose = (argc > 1 && strcmp(argv[1], "-v") == 0);

	commandline.logfile = stdout;
	commandline.quiet   = TRUE;

	setup_atmosphere();
	for (int n = 0; n < N_CONFIG; n++) check_pattern(&configs[n], verbose);
	for (int n = 0; n < N_CONFIG; n++)
		if (expect_polarizable(&configs[n])) check_profiles(&configs[n], verbose);

	/* --- Milne-Eddington: (config, point, ray, bound) -- ---------- */

	static const MECase me[]         = {{0, 2, 0, 1e-12}, {0, 2, 1, 1e-12}, {0, 3, 1, 1e-12}, {2, 4, 1, 1e-12}, {3, 1, 0, 1e-12}, {4, 3, 0, 1e-12}, {4, 6, 1, 1e-12}, {6, 3, 1, 1e-12}, {0, 1, 0, 1e-12}};
	spectrum.as                      = calloc(1, sizeof(ActiveSet));
	atmos.backgrflags                = calloc(1, sizeof(flags));
	atmos.backgrflags[0].ispolarized = TRUE;
	geometry.vboundary[TOP]          = ZERO;
	geometry.vboundary[BOTTOM]       = IRRADIATED;
	geometry.Itop                    = matrix_double(1, NRAYS);
	geometry.Ibottom                 = matrix_double(1, NRAYS);
	for (int n = 0; n < (int) (sizeof(me) / sizeof(me[0])); n++)
		check_milne_eddington(&me[n], verbose);

	return check_summary("test_zeeman");
}
