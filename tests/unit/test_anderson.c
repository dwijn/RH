/* ------- file: -------------------------- test_anderson.c ---------

	   Anderson acceleration (anderson.c, S-3 in docs/SOLVER.md).

	   -- Walker & Ni (2011), Theorem 2.2: for a linear map
		  G(x) = M x + b and full depth, Anderson type II with mixing 1
		  gives x_(k+1) = G(x_k^GMRES), the image of the k-th GMRES
		  iterate for (I - M) x = b from the same x_0. GMRES is computed
		  here by Arnoldi and a small least squares; agreement to 1e-10
		  for k = 1 .. 12 on a nonsymmetric M.
		  This checks the columns, the QR and the update together.
	   -- depth 0, and depth 5 without extrapolation, are the plain
		  iteration bit for bit on a nonlinear map;
	   -- collinear difference columns: the oldest are dropped, the
		  condition number of what is used is at most cond_max, and
		  nothing is NaN;
	   -- the safeguard: an extrapolated point whose residual grew is
		  rejected, the next point is the image of the last accepted
		  pair exactly, the history restarts from that pair, and the
		  counters are right; with the safeguard off it is accepted;
	   -- damping (theta = 1 repeats the step, theta = 0 gives g) and
		  abandoning (plain step, history restarted);
	   -- masked entries take g's value exactly and do not change the
		  others: the same bits as the problem without them;
	   -- coefficients per block equal separate runs on each block, and
		  a block that uses no column takes g exactly (mask 0);
	   -- two runs are bit-identical.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "anderson.h"
#include "check.h"

#define N    30
#define KMAX 12

/* --- Deterministic pseudo-random numbers in [0, 1) -- ------------- */

static unsigned long long rng_state = 88172645463325252ULL;

static double urand(void)
{
	rng_state ^= rng_state << 13;
	rng_state ^= rng_state >> 7;
	rng_state ^= rng_state << 17;
	return (rng_state >> 11) * (1.0 / 9007199254740992.0);
}

static double M[N][N], bvec[N];

static void G_linear(const double *x, double *g)
{
	for (int i = 0; i < N; i++) {
		double s = bvec[i];
		for (int j = 0; j < N; j++) s += M[i][j] * x[j];
		g[i] = s;
	}
}

static void G_nonlinear(const double *x, double *g, long n)
{
	for (long i = 0; i < n; i++)
		g[i] = 0.5 * cos(x[(i + 1) % n]) + 0.3 * x[i] + 0.1;
}

/* --- Run Anderson on the linear map for K steps; x[k] = x_k -- --- */

static void runLinear(int m, int K, double x[][N], double cond_max)
{
	long            start[2] = {0, N};
	AndersonOptions opt      = {1, start, 0, cond_max, 0.0};
	double          w[N], g[N], gamma[KMAX + 1];
	AndersonSolve   info[1];
	Anderson       *a = andersonInit(N, m);

	for (int i = 0; i < N; i++) w[i] = 1.0;
	for (int k = 0; k < K; k++) {
		G_linear(x[k], g);
		andersonStep(a, &opt, x[k], g, NULL, w, k >= 1, x[k + 1], gamma, info);
	}
	andersonFree(a);
}

/* --- GMRES iterate x_k for (I - M) x = b from x0 -- --------------- */

static void gmres(const double *x0, int k, double *xk)
{
	static double V[KMAX + 1][N];
	double        H[KMAX + 1][KMAX], r[N], beta, Q[KMAX + 1][KMAX],
		Rm[KMAX][KMAX], c[KMAX + 1], y[KMAX];

	memset(H, 0, sizeof(H));
	G_linear(x0, r); /* r = b + M x0 */
	for (int i = 0; i < N; i++) r[i] -= x0[i];
	beta = 0.0;
	for (int i = 0; i < N; i++) beta += r[i] * r[i];
	beta = sqrt(beta);
	for (int i = 0; i < N; i++) V[0][i] = r[i] / beta;

	for (int j = 0; j < k; j++) { /* Arnoldi with A = I - M, MGS twice */
		double Av[N];
		for (int i = 0; i < N; i++) {
			double s = V[j][i];
			for (int l = 0; l < N; l++) s -= M[i][l] * V[j][l];
			Av[i] = s;
		}
		for (int pass = 0; pass < 2; pass++) {
			for (int l = 0; l <= j; l++) {
				double h = 0.0;
				for (int i = 0; i < N; i++) h += V[l][i] * Av[i];
				H[l][j] += h;
				for (int i = 0; i < N; i++) Av[i] -= h * V[l][i];
			}
		}
		double h = 0.0;
		for (int i = 0; i < N; i++) h += Av[i] * Av[i];
		H[j + 1][j] = sqrt(h);
		for (int i = 0; i < N; i++) V[j + 1][i] = Av[i] / H[j + 1][j];
	}
	/* --- min || beta e1 - H y || by MGS QR of the (k+1) x k matrix -- */

	memset(Rm, 0, sizeof(Rm));
	for (int j = 0; j < k; j++) {
		double v[KMAX + 1];
		for (int i = 0; i <= k; i++) v[i] = H[i][j];
		for (int pass = 0; pass < 2; pass++) {
			for (int l = 0; l < j; l++) {
				double s = 0.0;
				for (int i = 0; i <= k; i++) s += Q[i][l] * v[i];
				Rm[l][j] += s;
				for (int i = 0; i <= k; i++) v[i] -= s * Q[i][l];
			}
		}
		double s = 0.0;
		for (int i = 0; i <= k; i++) s += v[i] * v[i];
		Rm[j][j] = sqrt(s);
		for (int i = 0; i <= k; i++) Q[i][j] = v[i] / Rm[j][j];
	}
	for (int l = 0; l < k; l++) c[l] = Q[0][l] * beta;
	for (int l = k - 1; l >= 0; l--) {
		double s = c[l];
		for (int j = l + 1; j < k; j++) s -= Rm[l][j] * y[j];
		y[l] = s / Rm[l][l];
	}
	for (int i = 0; i < N; i++) {
		double s = x0[i];
		for (int j = 0; j < k; j++) s += V[j][i] * y[j];
		xk[i] = s;
	}
}

static void testWalkerNi(void)
{
	static double x[KMAX + 2][N], y[KMAX + 2][N];
	double        xg[N], gx[N], err, scale;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) M[i][j] = (2.0 * urand() - 1.0) * 1.1 / sqrt(N);
		bvec[i] = urand();
		x[0][i] = 2.0 * urand() - 1.0;
	}
	runLinear(KMAX + 1, KMAX + 1, x, 1.0E14);

	for (int k = 1; k <= KMAX; k++) {
		gmres(x[0], k, xg);
		G_linear(xg, gx);
		err = scale = 0.0;
		for (int i = 0; i < N; i++) {
			err   = fmax(err, fabs(x[k + 1][i] - gx[i]));
			scale = fmax(scale, fabs(gx[i]));
		}
		CHECK(err <= 1.0E-10 * scale, "x_%d vs G(x_%d^GMRES): %.2e of %.2e", k + 1,
			  k, err, scale);
	}
	/* --- Two runs are bit-identical --                -------------- */

	memcpy(y[0], x[0], sizeof(y[0]));
	runLinear(KMAX + 1, KMAX + 1, y, 1.0E14);
	CHECK(memcmp(x, y, sizeof(x)) == 0, "two runs differ");
}

static void testPlain(void)
{
	long            start[2] = {0, N};
	AndersonOptions opt      = {1, start, 0, 1.0E8, 1.0};
	double          x[N], y[N], z[N], g[N], w[N], gamma[5];
	AndersonSolve   info[1];
	Anderson       *a0 = andersonInit(N, 0), *a5 = andersonInit(N, 5);
	int             same0 = 1, same5 = 1;

	for (int i = 0; i < N; i++) {
		x[i] = y[i] = z[i] = urand();
		w[i]               = 1.0;
	}
	for (int k = 0; k < 20; k++) {
		G_nonlinear(y, g, N);
		CHECK(andersonStep(a0, &opt, y, g, NULL, w, 1, y, gamma, info) ==
				  ANDERSON_PLAIN,
			  "depth 0 extrapolated at %d", k);
		G_nonlinear(z, g, N);
		andersonStep(a5, &opt, z, g, NULL, w, 0, z, gamma, info);
		G_nonlinear(x, g, N);
		memcpy(x, g, sizeof(x));
		same0 &= !memcmp(x, y, sizeof(x));
		same5 &= !memcmp(x, z, sizeof(x));
	}
	CHECK(same0, "depth 0 differs from the plain iteration");
	CHECK(same5, "depth 5 without extrapolation differs from the plain iteration");
	CHECK(a0->Nplain == 20 && a0->Nextrapolated == 0, "counts %ld %ld", a0->Nplain,
		  a0->Nextrapolated);
	andersonFree(a0);
	andersonFree(a5);
}

static void testCollinear(void)
{
	double        u[N], g[N], d[N], w[N], gamma[3];
	AndersonSolve s;
	Anderson     *a = andersonInit(N, 3);

	/* --- u_i = u_0 + i d on a linear map: all columns (I - M) d -- */

	for (int i = 0; i < N; i++) {
		u[i] = urand();
		d[i] = urand() - 0.5;
		w[i] = 1.0;
	}
	for (int p = 0; p < 4; p++) {
		double up[N];
		for (int i = 0; i < N; i++) up[i] = u[i] + p * d[i];
		G_linear(up, g);
		andersonPush(a, up, g, NULL);
	}
	s = andersonSolve(a, 0, N, w, 1.0E8, gamma);
	CHECK(s.Ncols == 1 && s.Ndrop == 2, "columns used %d, dropped %d", s.Ncols,
		  s.Ndrop);
	CHECK(s.cond <= 1.0E8 && isfinite(gamma[0]), "cond %.3e gamma %g", s.cond,
		  gamma[0]);

	/* --- Nearly collinear (1e-12): the same; 1e-3: all kept -- ----- */

	for (int eps = 0; eps < 2; eps++) {
		andersonReset(a, 0);
		for (int p = 0; p < 4; p++) {
			double up[N];
			for (int i = 0; i < N; i++)
				up[i] = u[i] + p * d[i] + (eps ? 1.0E-3 : 1.0E-12) * p * p * urand();
			G_linear(up, g);
			andersonPush(a, up, g, NULL);
		}
		s = andersonSolve(a, 0, N, w, 1.0E8, gamma);
		if (eps)
			CHECK(s.Ncols == 3 && s.cond <= 1.0E8, "1e-3: used %d cond %.3e", s.Ncols,
				  s.cond);
		else
			CHECK(s.Ncols == 1 && s.cond <= 1.0E8, "1e-12: used %d cond %.3e",
				  s.Ncols, s.cond);
		for (int j = 0; j < s.Ncols; j++)
			CHECK(isfinite(gamma[j]), "gamma[%d] = %g", j, gamma[j]);
	}
	andersonFree(a);
}

static void testSafeguard(void)
{
	long            start[2] = {0, N};
	AndersonOptions opt      = {1, start, 0, 1.0E8, 1.0};
	double          u[N], g[N], w[N], un[N], g2[N], gamma[4];
	AndersonSolve   info[1];
	Anderson       *a;

	for (int off = 0; off < 2; off++) {
		opt.safeguard = off ? 0.0 : 1.0;
		a             = andersonInit(N, 4);
		for (int i = 0; i < N; i++) {
			u[i] = urand();
			w[i] = 1.0;
		}
		G_linear(u, g);
		CHECK(andersonStep(a, &opt, u, g, NULL, w, 1, un, gamma, info) ==
				  ANDERSON_PLAIN,
			  "first step");
		memcpy(u, un, sizeof(u));
		G_linear(u, g);
		memcpy(g2, g, sizeof(g2));
		CHECK(andersonStep(a, &opt, u, g, NULL, w, 1, un, gamma, info) ==
				  ANDERSON_EXTRAPOLATED,
			  "second step");
		/* --- A residual 100 times larger at the trial point -- ----- */
		memcpy(u, un, sizeof(u));
		for (int i = 0; i < N; i++) g[i] = u[i] + 100.0 * (g2[i] - u[i]) + 1.0;
		if (off) {
			CHECK(andersonStep(a, &opt, u, g, NULL, w, 1, un, gamma, info) ==
					  ANDERSON_EXTRAPOLATED,
				  "safeguard off: not accepted");
			CHECK(a->Nrejected == 0, "rejected %ld", a->Nrejected);
		} else {
			CHECK(andersonStep(a, &opt, u, g, NULL, w, 1, un, gamma, info) ==
					  ANDERSON_REJECTED,
				  "not rejected");
			CHECK(!memcmp(un, g2, sizeof(un)), "next point is not G(last accepted)");
			CHECK(a->Npairs == 1 && a->Nrejected == 1 && a->Nextrapolated == 1 &&
					  a->Nplain == 1,
				  "pairs %d rejected %ld extrapolated %ld plain %ld", a->Npairs,
				  a->Nrejected, a->Nextrapolated, a->Nplain);
			/* --- The next evaluation extrapolates from 2 pairs -- --- */
			memcpy(u, un, sizeof(u));
			G_linear(u, g);
			CHECK(andersonStep(a, &opt, u, g, NULL, w, 1, un, gamma, info) ==
						  ANDERSON_EXTRAPOLATED &&
					  info[0].Ncols == 1,
				  "after restart: columns %d", info[0].Ncols);
			/* --- Damping and abandoning --              -------------- */
			double ud[N], g_new[N];
			memcpy(g_new, g, sizeof(g_new));
			andersonDamp(a, &opt, gamma, info, 1.0, ud);
			CHECK(!memcmp(ud, un, sizeof(ud)), "theta = 1 differs from the step");
			andersonDamp(a, &opt, gamma, info, 0.0, ud);
			CHECK(!memcmp(ud, g_new, sizeof(ud)), "theta = 0 is not g");
			andersonAbandon(a, ud);
			CHECK(!memcmp(ud, g_new, sizeof(ud)) && a->Npairs == 1 && !a->trial,
				  "abandon: pairs %d", a->Npairs);
		}
		andersonFree(a);
	}
}

static void testMask(void)
{
	enum { NM = 10 };
	double u[3][NM], g[3][NM], ur[3][NM - 2], gr[3][NM - 2], w[NM],
		gamma[2], gammar[2], out[NM], outr[NM - 2];
	unsigned char valid[NM];
	const int     skip[2] = {3, 7};
	Anderson     *a = andersonInit(NM, 2), *ar = andersonInit(NM - 2, 2);
	AndersonSolve s, sr;

	for (int i = 0; i < NM; i++) w[i] = 1.0 + i;
	for (int p = 0; p < 3; p++) {
		for (int i = 0, ir = 0; i < NM; i++) {
			u[p][i] = urand();
			g[p][i] = 0.5 * u[p][i] + urand() * 0.1;
			if (i == skip[0] || i == skip[1]) {
				if (p == 1) u[p][i] = g[p][i] = 1.0E300; /* garbage */
				continue;
			}
			ur[p][ir] = u[p][i];
			gr[p][ir] = g[p][i];
			ir++;
		}
		for (int i = 0; i < NM; i++) valid[i] = !(p == 1 && (i == 3 || i == 7));
		andersonPush(a, u[p], g[p], valid);
		andersonPush(ar, ur[p], gr[p], NULL);
	}
	double wr[NM - 2];
	for (int i = 0, ir = 0; i < NM; i++)
		if (i != 3 && i != 7) wr[ir++] = w[i];
	s  = andersonSolve(a, 0, NM, w, 1.0E8, gamma);
	sr = andersonSolve(ar, 0, NM - 2, wr, 1.0E8, gammar);
	andersonApply(a, 0, NM, gamma, s.Ncols, 1.0, out);
	andersonApply(ar, 0, NM - 2, gammar, sr.Ncols, 1.0, outr);
	CHECK(s.Ncols == 2 && sr.Ncols == 2, "columns %d %d", s.Ncols, sr.Ncols);
	CHECK(!memcmp(gamma, gammar, 2 * sizeof(double)), "gamma differs: %.17g %.17g",
		  gamma[0], gammar[0]);
	CHECK(out[3] == g[2][3] && out[7] == g[2][7], "masked entries are not g");
	int same = 1;
	for (int i = 0, ir = 0; i < NM; i++)
		if (i != 3 && i != 7) same &= (out[i] == outr[ir++]);
	CHECK(same, "unmasked entries differ from the reduced problem");
	andersonFree(a);
	andersonFree(ar);
}

static void testPerBlock(void)
{
	enum { NB = 12,
		   H  = 6 };
	long            start[3] = {0, H, NB}, s1[2] = {0, H};
	AndersonOptions opt = {2, start, 1, 1.0E8, 0.0}, o1 = {1, s1, 0, 1.0E8, 0.0};
	double          x[NB], xa[H], xb[H], g[NB], ga[H], gb[H], w[NB], gamma[6],
		gam1[3];
	AndersonSolve info[2], i1[1];
	Anderson     *a = andersonInit(NB, 3), *aa = andersonInit(H, 3),
			 *ab = andersonInit(H, 3);
	int same     = 1;

	for (int i = 0; i < NB; i++) {
		x[i] = urand();
		w[i] = 1.0;
	}
	memcpy(xa, x, sizeof(xa));
	memcpy(xb, x + H, sizeof(xb));
	for (int k = 0; k < 8; k++) {
		G_nonlinear(x, g, H);
		G_nonlinear(x + H, g + H, H);
		andersonStep(a, &opt, x, g, NULL, w, k >= 1, x, gamma, info);
		G_nonlinear(xa, ga, H);
		andersonStep(aa, &o1, xa, ga, NULL, w, k >= 1, xa, gam1, i1);
		G_nonlinear(xb, gb, H);
		andersonStep(ab, &o1, xb, gb, NULL, w, k >= 1, xb, gam1, i1);
		same &= !memcmp(x, xa, sizeof(xa)) && !memcmp(x + H, xb, sizeof(xb));
	}
	CHECK(same, "per-block coefficients differ from separate runs");
	andersonFree(a);
	andersonFree(aa);
	andersonFree(ab);

	/* --- A block that uses no column (its pairs all equal, so every
		   column is zero) is not extrapolated: mask 0 and g exactly,
		   while the other block extrapolates -- -------------------- */

	a = andersonInit(NB, 3);
	for (int i = 0; i < NB; i++) x[i] = urand();
	for (int k = 0; k < 3; k++) {
		G_nonlinear(x, g, H);
		for (int i = H; i < NB; i++) g[i] = x[i] = 0.25 + i; /* fixed */
		andersonStep(a, &opt, x, g, NULL, w, k >= 1, x, gamma, info);
	}
	const unsigned char *mask = andersonMask(a);
	int                  zero = 1, kept = 1;
	for (int i = H; i < NB; i++) zero &= !mask[i] && x[i] == 0.25 + i;
	for (int i = 0; i < H; i++) kept &= mask[i];
	CHECK(info[1].Ncols == 0 && info[0].Ncols > 0 && zero && kept,
		  "block without columns: columns %d / %d, mask or values wrong",
		  info[0].Ncols, info[1].Ncols);
	andersonFree(a);
}

int main(void)
{
	testWalkerNi();
	testPlain();
	testCollinear();
	testSafeguard();
	testMask();
	testPerBlock();
	return check_summary("test_anderson");
}

/* ------- end ---------------------------- test_anderson.c --------- */
