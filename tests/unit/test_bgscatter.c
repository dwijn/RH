/* ------- file: -------------------------- test_bgscatter.c --------

	   The background scattering iteration (bgscatter.c, N-51) on a
	   semi-infinite two-level scattering slab, S = (1 - eps) J + eps B
	   with B = 1, solved with the 1-D linear short characteristics.

	   Checks:
	   -- the row of the operator from Piecewise_Linear_1D_tri is the
		  exact tridiagonal of dJ/du (columns by linearity: one formal
		  solution per point);
	   -- the row from Piecewise_Bezier3_1D_tri is the derivative of
		  the Bezier solution where no limiter switches (finite
		  differences on a smooth S);
	   -- bgScatterIterate from J = 0 reaches the solution of the
		  linear system (I - Lambda sigma) J = J(0) (dense solve) with
		  the tridiagonal and the diagonal operator, eps = 1e-4 and
		  1e-6, within the solves recorded below;
	   -- BACKGROUND_N_MAX_ITER bounds the solves, and the result then
		  is still positive;
	   -- BACKGROUND_ITER_LIMIT bounds the error of J, not only the
		  last step (N-66): from J = 0, limits 1e-3 to 1e-7, eps 1e-4
		  and 1e-6, every operator.
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
#include "bgscatter.h"
#include "geometryops.h"
#include "check.h"
#include "state.h"

/* --- bgscatter.c asks the geometry (geometryOps) whether it
	   supports the iteration: here it does, with the linear fallback of
	   rhf1d (N-52) --                                -------------- */

static const char *bgSupported(RHState *rh)
{
	(void) rh;
	return NULL;
}

static bool_t bgLinearFallback(RHState *rh)
{
	return rh->input.S_interpolation != S_LINEAR;
}

const GeometryOps geometryOps = {.name               = "unit test",
								 .bg_supported       = bgSupported,
								 .bg_linear_fallback = bgLinearFallback};

#define NMU 3
static double mus[NMU], wmus[NMU];

typedef struct {
	int    N;
	double eps, *S, *I, *Tri;
	bool_t bezier;
} Slab;

/* --- Gauss-Legendre on [0, 1] with 3 points --      -------------- */

static void gauss3(void)
{
	const double x = sqrt(0.6);

	mus[0]  = 0.5 * (1.0 - x);
	mus[1]  = 0.5;
	mus[2]  = 0.5 * (1.0 + x);
	wmus[0] = wmus[2] = 5.0 / 18.0;
	wmus[1]           = 8.0 / 18.0;
}

/* --- An atmosphere logarithmic in optical depth, chi = 1, from
	   1e-4 to beyond the thermalization depth (100 / sqrt(eps)), the
	   bottom irradiated with B = 1 --                -------------- */

static void slab_init(Slab *s, int N, double eps, bool_t bezier)
{
	const double tau_min = 1.0e-4, tau_max = 100.0 / sqrt(eps);

	s->N      = N;
	s->eps    = eps;
	s->bezier = bezier;
	s->S      = malloc(N * sizeof(double));
	s->I      = malloc(N * sizeof(double));
	s->Tri    = malloc(3 * N * sizeof(double));

	rh->geometry->Ndep   = N;
	rh->geometry->height = realloc(rh->geometry->height, N * sizeof(double));
	for (int k = 0; k < N; k++)
		rh->geometry->height[k] =
			-tau_min * pow(tau_max / tau_min, (double) k / (N - 1));
	for (int mu = 0; mu < NMU; mu++) {
		rh->geometry->Itop[0][mu]    = 0.0;
		rh->geometry->Ibottom[0][mu] = 1.0;
	}
}

static void slab_free(Slab *s)
{
	free(s->S);
	free(s->I);
	free(s->Tri);
}

/* --- The sweep of bgScatterIterate: S = eps + (1 - eps) u, J and
	   A[3k + i] = dJ_k / du_k-1+i, as in src/geometry/plane1d/formal.c
	   --                                              -------------- */

static void sweep(RHState *rh, void *ctx, const double *u, double *J, double *A,
				  double *albedo)
{
	Slab  *s   = (Slab *) ctx;
	int    N   = s->N;
	double sca = 1.0 - s->eps, *chi = malloc(N * sizeof(double)), dm, dp;

	(void) rh;
	for (int k = 0; k < N; k++) {
		chi[k]  = 1.0;
		s->S[k] = s->eps + sca * u[k];
		J[k]    = 0.0;
		if (albedo) albedo[k] = sca;
	}
	for (int k = 0; k < 3 * N; k++) A[k] = 0.0;
	for (int mu = 0; mu < NMU; mu++) {
		for (int to_obs = 0; to_obs <= 1; to_obs++) {
			double w = 0.5 * wmus[mu];

			if (s->bezier)
				Piecewise_Bezier3_1D_tri(rh, 0, mu, to_obs, chi, s->S, s->I, NULL, s->Tri);
			else
				Piecewise_Linear_1D_tri(rh, 0, mu, to_obs, chi, s->S, s->I, NULL, s->Tri);
			for (int k = 0; k < N; k++) {
				J[k] += w * s->I[k];
				dm = (to_obs) ? s->Tri[3 * k + 2] : s->Tri[3 * k];
				dp = (to_obs) ? s->Tri[3 * k] : s->Tri[3 * k + 2];
				if (k > 0) A[3 * k] += w * dm * sca;
				A[3 * k + 1] += w * s->Tri[3 * k + 1] * sca;
				if (k < N - 1) A[3 * k + 2] += w * dp * sca;
			}
		}
	}
	free(chi);
}

/* --- Dense Gaussian elimination with partial pivoting, M x = b -- */

static void dense_solve(int N, double *M, double *b)
{
	for (int c = 0; c < N; c++) {
		int p = c;
		for (int r = c + 1; r < N; r++)
			if (fabs(M[r * N + c]) > fabs(M[p * N + c])) p = r;
		for (int j = 0; j < N; j++) {
			double t     = M[c * N + j];
			M[c * N + j] = M[p * N + j];
			M[p * N + j] = t;
		}
		double t = b[c];
		b[c]     = b[p];
		b[p]     = t;
		for (int r = c + 1; r < N; r++) {
			double f = M[r * N + c] / M[c * N + c];
			for (int j = c; j < N; j++) M[r * N + j] -= f * M[c * N + j];
			b[r] -= f * b[c];
		}
	}
	for (int r = N - 1; r >= 0; r--) {
		for (int j = r + 1; j < N; j++) b[r] -= M[r * N + j] * b[j];
		b[r] /= M[r * N + r];
	}
}

/* --- The exact solution and the exact Jacobian (J is affine in u:
	   column j = J(e_j) - J(0)) --                    -------------- */

static void exact(Slab *s, double *Jx, double *Lam)
{
	int     N = s->N;
	double *u = calloc(N, sizeof(double)), *J0 = malloc(N * sizeof(double)),
		   *J1 = malloc(N * sizeof(double)), *A = malloc(3 * N * sizeof(double)),
		   *M = malloc(N * N * sizeof(double));

	sweep(rh, s, u, J0, A, NULL);
	for (int j = 0; j < N; j++) {
		u[j] = 1.0;
		sweep(rh, s, u, J1, A, NULL);
		u[j] = 0.0;
		for (int k = 0; k < N; k++) Lam[k * N + j] = J1[k] - J0[k];
	}
	for (int k = 0; k < N; k++)
		for (int j = 0; j < N; j++) M[k * N + j] = (k == j) - Lam[k * N + j];
	memcpy(Jx, J0, N * sizeof(double));
	dense_solve(N, M, Jx);
	free(u);
	free(J0);
	free(J1);
	free(A);
	free(M);
}

/* --- The tridiagonal of the linear solver is that of dJ/du -- --- */

static void check_tri_linear(double eps)
{
	Slab    s;
	int     N   = 48;
	double *Lam = malloc(N * N * sizeof(double)), *Jx = malloc(N * sizeof(double)),
		   *u = calloc(N, sizeof(double)), *J = malloc(N * sizeof(double)),
		   *A = malloc(3 * N * sizeof(double)), err = 0.0;

	slab_init(&s, N, eps, FALSE);
	exact(&s, Jx, Lam);
	sweep(rh, &s, u, J, A, NULL);
	for (int k = 0; k < N; k++)
		for (int i = 0; i < 3; i++) {
			int j = k - 1 + i;
			if (j < 0 || j >= N) continue;
			err = fmax(err, fabs(A[3 * k + i] - Lam[k * N + j]));
		}
	CHECK(err < 1.0e-13, "linear: tridiagonal of dJ/du off by %.2e (eps %.0e)", err, eps);
	slab_free(&s);
	free(Lam);
	free(Jx);
	free(u);
	free(J);
	free(A);
}

/* --- The Bezier row against finite differences on a smooth S
	   (u = J of the exact linear solution, no limiter switch
	   crossed by a relative step of 1e-7) --          -------------- */

static void check_tri_bezier(void)
{
	Slab    s;
	int     N   = 48;
	double *Lam = malloc(N * N * sizeof(double)), *u = malloc(N * sizeof(double)),
		   *J0 = malloc(N * sizeof(double)), *J1 = malloc(N * sizeof(double)),
		   *A = malloc(3 * N * sizeof(double)), *A1 = malloc(3 * N * sizeof(double)),
		   err = 0.0, du;

	slab_init(&s, N, 1.0e-2, FALSE);
	exact(&s, u, Lam); /* a smooth, monotonic u */
	s.bezier = TRUE;
	sweep(rh, &s, u, J0, A, NULL);
	for (int j = 0; j < N; j++) {
		du = 1.0e-7 * u[j];
		u[j] += du;
		sweep(rh, &s, u, J1, A1, NULL);
		u[j] -= du;
		for (int k = (j > 0) ? j - 1 : 0; k <= j + 1 && k < N; k++)
			err = fmax(err, fabs((J1[k] - J0[k]) / du - A[3 * k + (j - k + 1)]));
	}
	CHECK(err < 1.0e-6, "Bezier: tridiagonal of dJ/du off by %.2e", err);
	slab_free(&s);
	free(Lam);
	free(u);
	free(J0);
	free(J1);
	free(A);
	free(A1);
}

/* --- The iteration from J = 0 --                    -------------- */

static void check_iterate(double eps, enum bg_operator op, int depth,
						  int max_solves)
{
	Slab    s;
	int     N   = 64;
	double *Lam = malloc(N * N * sizeof(double)), *Jx = malloc(N * sizeof(double)),
		   *Jdag = calloc(N, sizeof(double)), *Jfs = malloc(N * sizeof(double)),
		   *A = malloc(3 * N * sizeof(double)), *alb = malloc(N * sizeof(double)),
		   *J = malloc(N * sizeof(double)), err = 0.0, kap;
	BgScatterStats st;

	slab_init(&s, N, eps, FALSE);
	exact(&s, Jx, Lam);
	rh->input.bg_operator  = op;
	rh->input.bg_depth     = depth;
	rh->input.bg_iterLimit = 1.0e-10;
	rh->input.bg_NmaxIter  = 1000;
	for (int k = 0; k < N; k++) alb[k] = 1.0 - eps;

	bgScatterStats(rh, &st, TRUE);
	sweep(rh, &s, Jdag, Jfs, A, NULL);
	kap = 0.0;
	bgScatterIterate(rh, N, Jdag, Jfs, A, alb, J, &kap, sweep, &s);
	bgScatterStats(rh, &st, TRUE);
	for (int k = 0; k < N; k++) err = fmax(err, fabs(J[k] / Jx[k] - 1.0));

	CHECK(st.Nunconverged == 0 && err < 1.0e-8,
		  "%s depth %d eps %.0e: relative error %.2e after %ld solves%s",
		  (op == BG_OPERATOR_TRIDIAGONAL) ? "TRIDIAGONAL" : "DIAGONAL", depth, eps,
		  err, st.Nsweeps, st.Nunconverged ? " (not converged)" : "");
	CHECK(st.Nsweeps <= max_solves, "%s depth %d eps %.0e: %ld solves > %d",
		  (op == BG_OPERATOR_TRIDIAGONAL) ? "TRIDIAGONAL" : "DIAGONAL", depth, eps,
		  st.Nsweeps, max_solves);
	if (getenv("RH_TEST_VERBOSE"))
		fprintf(stderr, "  %s depth %d eps %.0e: %ld solves, error %.2e\n",
				(op == BG_OPERATOR_TRIDIAGONAL) ? "TRIDIAGONAL" : "DIAGONAL", depth,
				eps, st.Nsweeps, err);

	/* --- BACKGROUND_N_MAX_ITER bounds the work; J stays positive -- */

	rh->input.bg_NmaxIter = 3;
	kap                   = 0.0;
	bgScatterIterate(rh, N, Jdag, Jfs, A, alb, J, &kap, sweep, &s);
	bgScatterStats(rh, &st, TRUE);
	CHECK(st.Nsweeps == 3 && st.Nunconverged == 1, "N_MAX_ITER = 3: %ld solves",
		  st.Nsweeps);
	for (int k = 0; k < N; k++) CHECK(J[k] > 0.0, "J[%d] = %g", k, J[k]);

	slab_free(&s);
	free(Lam);
	free(Jx);
	free(Jdag);
	free(Jfs);
	free(A);
	free(alb);
	free(J);
}

/* --- The limit bounds the error (N-66). The step alone did not:
	   stopping on the step, 10 of these 12 cases failed, the diagonal
	   operator at eps = 1e-6 and limit 1e-3 with an error 41 times the
	   limit --                                        -------------- */

static void check_error_bound(double eps, enum bg_operator op, int depth)
{
	Slab    s;
	int     N   = 64;
	double *Lam = malloc(N * N * sizeof(double)), *Jx = malloc(N * sizeof(double)),
		   *Jdag = calloc(N, sizeof(double)), *Jfs = malloc(N * sizeof(double)),
		   *A = malloc(3 * N * sizeof(double)), *alb = malloc(N * sizeof(double)),
		   *J = malloc(N * sizeof(double)), err, kap;
	BgScatterStats st;

	slab_init(&s, N, eps, FALSE);
	exact(&s, Jx, Lam);
	rh->input.bg_operator = op;
	rh->input.bg_depth    = depth;
	rh->input.bg_NmaxIter = 5000;
	for (int k = 0; k < N; k++) alb[k] = 1.0 - eps;

	for (double tol = 1.0e-3; tol > 1.0e-8; tol *= 0.01) {
		rh->input.bg_iterLimit = tol;
		kap                    = 0.0;
		bgScatterStats(rh, &st, TRUE);
		sweep(rh, &s, Jdag, Jfs, A, NULL);
		bgScatterIterate(rh, N, Jdag, Jfs, A, alb, J, &kap, sweep, &s);
		bgScatterStats(rh, &st, TRUE);
		err = 0.0;
		for (int k = 0; k < N; k++) err = fmax(err, fabs(J[k] / Jx[k] - 1.0));
		CHECK(st.Nunconverged == 0 && err <= tol,
			  "%s depth %d eps %.0e limit %.0e: relative error %.2e after %ld"
			  " solves%s",
			  (op == BG_OPERATOR_TRIDIAGONAL) ? "TRIDIAGONAL" : "DIAGONAL", depth,
			  eps, tol, err, st.Nsweeps, st.Nunconverged ? " (not converged)" : "");
		if (getenv("RH_TEST_VERBOSE"))
			fprintf(stderr, "  %s depth %d eps %.0e limit %.0e: %ld solves, error"
							" %.2e, kappa %.3g\n",
					(op == BG_OPERATOR_TRIDIAGONAL) ? "TRIDIAGONAL" : "DIAGONAL",
					depth, eps, tol, st.Nsweeps, err, kap);
	}
	slab_free(&s);
	free(Lam);
	free(Jx);
	free(Jdag);
	free(Jfs);
	free(A);
	free(alb);
	free(J);
}

/* --- From the solution of a nearby problem (the populations of the
	   next iteration: eps 5 % larger) the iteration takes a few steps,
	   and the error stays below the limit. kappa from those steps
	   alone is 3-5 times too small with the tridiagonal operator
	   (RH_TEST_VERBOSE); the slab is forgiving enough that the error
	   still meets the limit without the last kappa, unlike the MURaM
	   column of N-66 (Mn II 257.68 nm off by 3e-2 at limit 1e-4) -- */

static double warm_error(double eps, enum bg_operator op, double tol,
						 int keep, double *kap_out)
{
	Slab    s;
	int     N   = 64;
	double *Lam = malloc(N * N * sizeof(double)), *Jx = malloc(N * sizeof(double)),
		   *J0 = calloc(N, sizeof(double)), *J1 = malloc(N * sizeof(double)),
		   *Jfs = malloc(N * sizeof(double)), *A = malloc(3 * N * sizeof(double)),
		   *alb = malloc(N * sizeof(double)), *J = malloc(N * sizeof(double)),
		   err = 0.0, kap = 0.0;

	slab_init(&s, N, eps, FALSE);
	rh->input.bg_operator  = op;
	rh->input.bg_depth     = 3;
	rh->input.bg_NmaxIter  = 5000;
	rh->input.bg_iterLimit = tol;
	for (int k = 0; k < N; k++) alb[k] = 1.0 - eps;
	sweep(rh, &s, J0, Jfs, A, NULL);
	bgScatterIterate(rh, N, J0, Jfs, A, alb, J1, &kap, sweep, &s);

	s.eps *= 1.05;
	exact(&s, Jx, Lam);
	for (int k = 0; k < N; k++) alb[k] = 1.0 - s.eps;
	if (!keep) kap = 0.0;
	sweep(rh, &s, J1, Jfs, A, NULL);
	bgScatterIterate(rh, N, J1, Jfs, A, alb, J, &kap, sweep, &s);
	for (int k = 0; k < N; k++) err = fmax(err, fabs(J[k] / Jx[k] - 1.0));
	*kap_out = kap;

	slab_free(&s);
	free(Lam);
	free(Jx);
	free(J0);
	free(J1);
	free(Jfs);
	free(A);
	free(alb);
	free(J);
	return err;
}

static void check_warm_start(double eps, enum bg_operator op)
{
	double kap, kap0, err, err0;

	for (double tol = 1.0e-3; tol > 1.0e-8; tol *= 0.01) {
		err  = warm_error(eps, op, tol, TRUE, &kap);
		err0 = warm_error(eps, op, tol, FALSE, &kap0);
		CHECK(err <= tol, "%s eps %.0e limit %.0e: warm start off by %.2e",
			  (op == BG_OPERATOR_TRIDIAGONAL) ? "TRIDIAGONAL" : "DIAGONAL", eps,
			  tol, err);
		if (getenv("RH_TEST_VERBOSE"))
			fprintf(stderr,
					"  warm %s eps %.0e limit %.0e: error %.2e (kappa %.3g);"
					" without the last kappa %.2e (kappa %.3g)\n",
					(op == BG_OPERATOR_TRIDIAGONAL) ? "TRIDIAGONAL" : "DIAGONAL", eps,
					tol, err, kap, err0, kap0);
	}
}

/* --- Iterate before the main sweep: the first sweep iterates every
	   wavelength first; one that needs no iteration from there is
	   dropped, and bgScatterAfter picks it up again when the main
	   sweep finds the scattering changes J -- --       -------------- */

static void check_before_after(void)
{
	Slab    s;
	int     N   = 64;
	double  eps = 1.0e-4, err = 0.0;
	double *Lam = malloc(N * N * sizeof(double)), *Jx = malloc(N * sizeof(double)),
		   *Jdag = calloc(N, sizeof(double)), *Js = malloc(N * sizeof(double)),
		   *Jfs = malloc(N * sizeof(double)), *A = malloc(3 * N * sizeof(double)),
		   *alb = malloc(N * sizeof(double)), *u = malloc(N * sizeof(double));

	slab_init(&s, N, eps, FALSE);
	exact(&s, Jx, Lam);
	rh->commandline.logfile   = stderr; /* bgScatterSetup's message */
	rh->spectrum.Nspect       = 2;
	rh->input.S_interpolation = S_LINEAR;
	rh->input.bg_operator     = BG_OPERATOR_TRIDIAGONAL;
	rh->input.bg_depth        = 3;
	rh->input.bg_iterLimit    = 1.0e-10;
	rh->input.bg_NmaxIter     = 1000;
	bgScatterSetup(rh);
	CHECK(bgScatterPending(rh, 0) && bgScatterPending(rh, 1),
		  "the first sweep does not iterate first");

	bgScatterBefore(rh, 0, N, Jdag, Js, sweep, &s);
	for (int k = 0; k < N; k++) err = fmax(err, fabs(Js[k] / Jx[k] - 1.0));
	CHECK(err < 1.0e-8, "J* from J = 0 off by %.2e", err);
	CHECK(bgScatterPending(rh, 0), "an iterated wavelength was dropped");

	memcpy(u, Js, N * sizeof(double));
	bgScatterBefore(rh, 0, N, u, Js, sweep, &s);
	CHECK(!bgScatterPending(rh, 0), "a converged wavelength stays pending");
	CHECK(bgScatterPending(rh, 1), "wavelength 1 was dropped");

	sweep(rh, &s, Jx, Jfs, A, alb);
	bgScatterAfter(rh, 0, N, Jx, Jfs, A, alb);
	CHECK(!bgScatterPending(rh, 0), "bgScatterAfter: converged J made it pending");
	for (int k = 0; k < N; k++) u[k] = 0.5 * Jx[k];
	sweep(rh, &s, u, Jfs, A, alb);
	bgScatterAfter(rh, 0, N, u, Jfs, A, alb);
	CHECK(bgScatterPending(rh, 0), "bgScatterAfter: J off by 2 did not make it pending");

	slab_free(&s);
	free(Lam);
	free(Jx);
	free(Jdag);
	free(Js);
	free(Jfs);
	free(A);
	free(alb);
	free(u);
}

/* --- With S_BEZIER3, the choice of the linear solver per wavelength
	   (N-52): above the threshold only, decided afresh, and a query
	   (albedo < 0) keeps the last choice -- --       -------------- */

static void check_linear_choice(void)
{
	BgScatterStats st;

	rh->commandline.logfile   = stderr; /* bgScatterSetup's message */
	rh->spectrum.Nspect       = 4;
	rh->input.S_interpolation = S_BEZIER3;
	rh->input.bg_NmaxIter     = 50;
	rh->input.bg_linearAlbedo = 0.9;
	rh->input.bg_linearTau    = 1.0;
	bgScatterSetup(rh);

	CHECK(!bgScatterLinear(rh, 0, -1.0), "query set the flag");
	CHECK(!bgScatterLinear(rh, 0, 0.9), "albedo = threshold chose linear");
	CHECK(bgScatterLinear(rh, 0, 0.95), "albedo 0.95 > 0.9 did not choose linear");
	CHECK(bgScatterLinear(rh, 0, -1.0), "a query lost the choice");
	CHECK(!bgScatterLinear(rh, 1, -1.0), "wavelength 1 chose linear");
	bgScatterStats(rh, &st, TRUE);
	CHECK(st.Nlinear == 1 && st.Nto == 1 && st.Nback == 0,
		  "linear %ld, +%ld -%ld, not 1, +1 -0", st.Nlinear, st.Nto, st.Nback);
	CHECK(!bgScatterLinear(rh, 0, 0.1), "albedo 0.1 kept the linear solver");
	CHECK(!bgScatterLinear(rh, 0, -1.0), "a query lost the switch back");
	bgScatterStats(rh, &st, TRUE);
	CHECK(st.Nlinear == 0 && st.Nto == 0 && st.Nback == 1,
		  "linear %ld, +%ld -%ld, not 0, +0 -1", st.Nlinear, st.Nto, st.Nback);
}

int main(void)
{
	gauss3();
	rh->geometry->Nrays             = NMU;
	rh->geometry->muz               = mus;
	rh->geometry->wmu               = wmus;
	rh->geometry->Itop              = matrix_double(rh, 1, NMU);
	rh->geometry->Ibottom           = matrix_double(rh, 1, NMU);
	rh->geometry->vboundary[TOP]    = ZERO;
	rh->geometry->vboundary[BOTTOM] = IRRADIATED;

	check_tri_linear(1.0e-4);
	check_tri_bezier();

	/* --- Solves to 1e-10, limit about 1.5 x the measured: 24, 66,
		   520 (the diagonal operator) and 152 (no Anderson); with the
		   test on the step alone (before N-66) 23, 36, 494, 95 -- -- */

	check_iterate(1.0e-4, BG_OPERATOR_TRIDIAGONAL, 3, 35);
	check_iterate(1.0e-6, BG_OPERATOR_TRIDIAGONAL, 3, 100);
	check_iterate(1.0e-4, BG_OPERATOR_DIAGONAL, 3, 780);
	check_iterate(1.0e-6, BG_OPERATOR_TRIDIAGONAL, 0, 230);
	check_error_bound(1.0e-4, BG_OPERATOR_TRIDIAGONAL, 3);
	check_error_bound(1.0e-6, BG_OPERATOR_TRIDIAGONAL, 3);
	check_error_bound(1.0e-6, BG_OPERATOR_DIAGONAL, 3);
	check_error_bound(1.0e-6, BG_OPERATOR_TRIDIAGONAL, 0);
	check_warm_start(1.0e-4, BG_OPERATOR_TRIDIAGONAL);
	check_warm_start(1.0e-6, BG_OPERATOR_TRIDIAGONAL);
	check_warm_start(1.0e-6, BG_OPERATOR_DIAGONAL);
	check_before_after();
	check_linear_choice();

	return check_summary("test_bgscatter");
}
/* ------- end ---------------------------- test_bgscatter.c -------- */
