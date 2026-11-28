/* ------- file: -------------------------- test_prdh.c -------------

	   Gas-frame interpolation of hybrid PRD (prdh.c, P-12, P-17):

	   1. prdhBracket selects the same points as a scan of the whole
		  grid with the same comparisons (random grids and shifts).
	   2. The weights are a partition of unity: every grid point gets a
		  total weight of 1 from all wavelengths of a ray, for any
		  shift, also at both ends of the grid; with no shift each
		  wavelength gives weight 1 to itself and nothing to the others
		  (before P-17 the last point got 2 from the table and the
		  first 0 with PRDH_LIMIT_MEM). Also where two grid points 1
		  ulp apart shift to the same wavelength (N-42: before, the
		  upper one took every point below it with weight 1).
	   3. addJgasRay gives the same Jgas to the bit from the table
		  (initPRDhTables) and with PRDH_LIMIT_MEM, with locations at
		  rest among moving ones; at rest Jgas = J.
	   4. prdhRhoWeights interpolates linearly inside the line's grid
		  and is constant beyond it.
	   --                                              -------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "inputs.h"
#include "prdh.h"
#include "check.h"
#include "state.h"

static double urand(void)
{
	return rand() / (RAND_MAX + 1.0);
}

/* --- Sorted grid without duplicates: uneven steps, dense clusters
	   (like the lines of a wavelength table) and a few wide gaps -- */

static void makeGrid(int N)
{
	int    n;
	double l = 250.0;

	rh->spectrum.Nspect = N;
	rh->spectrum.lambda = (double *) realloc(rh->spectrum.lambda, N * sizeof(double));
	for (n = 0; n < N; n++) {
		double u = urand(), step;

		step = (u < 0.6) ? 1e-4 * (0.5 + urand()) : (u < 0.95) ? 2e-3 * urand() + 1e-6
															   : 0.05;
		l += step;
		rh->spectrum.lambda[n] = l;
	}
}

static void check_brackets(void)
{
	int         t, la, idx, N = 400;
	double      fac;
	PRDhBracket b;

	for (t = 0; t < 20; t++) {
		makeGrid(N);
		fac = (t == 0) ? 1.0 : 1.0 + ((t % 2) ? 1 : -1) * 3e-4 * urand() * (t % 5 + 1);
		for (la = 0; la < N; la++) {
			prdhBracket(rh, la, fac, &b);

			/* --- Membership: a point is in the first range exactly
				   when the scan of the old code over the whole grid
				   takes it, with the edge rule of P-17 (<= gas at the
				   lower end, > gas at the upper end; the ends by
				   index, N-42), and likewise for the second -- ----- */

			for (idx = 0; idx < N; idx++) {
				const double *lambda = rh->spectrum.lambda;
				double        prv = b.prv, gas = b.gas, nxt = b.nxt;
				int           in1  = (la != 0) ? (lambda[idx] > prv && lambda[idx] <= gas)
											   : (lambda[idx] <= gas);
				int           in2  = (la != N - 1) ? (lambda[idx] > gas && lambda[idx] < nxt)
												   : (lambda[idx] > gas);
				int           got1 = (idx >= b.r[0] && idx < b.r[1]);
				int           got2 = (idx >= b.r[1] && idx < b.r[2]);

				CHECK(got1 == in1 && got2 == in2,
					  "fac %.9f la %d idx %d: ranges [%ld %ld %ld], scans %d %d", fac, la, idx,
					  b.r[0], b.r[1], b.r[2], in1, in2);
			}
		}
	}
}

static void check_partition(void)
{
	int         t, la, N = 300;
	long        idx;
	double     *W, fac, maxdev;
	PRDhBracket b;

	W = (double *) malloc(N * sizeof(double));
	for (t = 0; t < 12; t++) {
		makeGrid(N);
		fac = (t == 0) ? 1.0 : 1.0 + ((t % 2) ? 1 : -1) * 1e-3 * urand() * t;
		memset(W, 0, N * sizeof(double));
		for (la = 0; la < N; la++) {
			prdhBracket(rh, la, fac, &b);
			for (idx = b.r[0]; idx < b.r[2]; idx++) {
				double w = prdhWeight(rh, &b, idx);

				W[idx] += w;
				if (fac == 1.0)
					CHECK(idx == la && w == 1.0, "no shift: la %d gives %g to %ld", la, w, idx);
			}
		}
		maxdev = 0.0;
		for (idx = 0; idx < N; idx++) maxdev = fmax(maxdev, fabs(W[idx] - 1.0));
		CHECK(maxdev <= 4e-15, "fac %.9f: total weights differ from 1 by %.2e (first %g, last %g)",
			  fac, maxdev, W[0], W[N - 1]);
		if (fac == 1.0)
			CHECK(W[0] == 1.0 && W[N - 1] == 1.0, "no shift: first %g, last %g", W[0], W[N - 1]);
	}
	free(W);
}

/* --- N-42: grid points n0 and n0 + 1 one ulp apart, and a shift
	   that maps both to the same gas-frame wavelength --  --------- */

static void check_collapse(void)
{
	const int   N = 300, n0 = 150;
	int         t, la, found = 0;
	long        idx;
	double     *W, fac = 1.0, maxdev;
	PRDhBracket b;

	W = (double *) malloc(N * sizeof(double));
	makeGrid(N);
	rh->spectrum.lambda[n0 + 1] = nextafter(rh->spectrum.lambda[n0], 1e300);
	CHECK(rh->spectrum.lambda[n0 + 1] < rh->spectrum.lambda[n0 + 2], "grid not sorted");

	for (t = 1; t <= 10000 && !found; t++) {
		fac   = 1.0 - 1e-7 * t;
		found = (rh->spectrum.lambda[n0] * fac == rh->spectrum.lambda[n0 + 1] * fac);
	}
	CHECK(found, "no shift maps the pair to one wavelength");

	memset(W, 0, N * sizeof(double));
	for (la = 0; la < N; la++) {
		prdhBracket(rh, la, fac, &b);
		for (idx = b.r[0]; idx < b.r[2]; idx++) W[idx] += prdhWeight(rh, &b, idx);

		/* --- The collapsed interval is empty from both sides -- --- */

		if (la == n0)
			CHECK(b.r[2] == b.r[1], "fac %.9f: la %d takes [%ld %ld) above", fac, la, b.r[1], b.r[2]);
		if (la == n0 + 1)
			CHECK(b.r[1] == b.r[0], "fac %.9f: la %d takes [%ld %ld) below", fac, la, b.r[0], b.r[1]);
	}
	maxdev = 0.0;
	for (idx = 0; idx < N; idx++) maxdev = fmax(maxdev, fabs(W[idx] - 1.0));
	CHECK(maxdev <= 4e-15, "fac %.9f: total weights differ from 1 by %.2e", fac, maxdev);
	free(W);
}

static void check_jgas(void)
{
	const int Nrays = 3, Nspace = 7, N = 250;
	int       mu, to_obs, k, la, rest;
	double  **Jtab, **Jmem, **J, *I, wmu, maxrel = 0.0;

	makeGrid(N);
	rh->atmos.Nrays       = Nrays;
	rh->atmos.Nspace      = Nspace;
	rh->atmos.Nactiveatom = 0;
	rh->spectrum.v_los    = matrix_double(rh, Nrays, Nspace);
	for (mu = 0; mu < Nrays; mu++)
		for (k = 0; k < Nspace; k++)
			rh->spectrum.v_los[mu][k] = (k % 3 == 0) ? 0.0 : 2.0e4 * (2.0 * urand() - 1.0); /* m/s */

	rh->spectrum.Jgas = matrix_double(rh, N, Nspace);
	initJgasLocks(rh);
	rh->input.prdh_limit_mem = FALSE;
	CHECK(initPRDhTables(rh) == RH_OK, "initPRDhTables failed: %s", rh->error.message);

	Jtab = matrix_double(rh, N, Nspace);
	Jmem = matrix_double(rh, N, Nspace);
	J    = matrix_double(rh, N, Nspace);
	I    = (double *) malloc(Nspace * sizeof(double));

	for (rest = 0; rest < 2; rest++) {
		rh->input.prdh_limit_mem = rest; /* first the table, then on the fly */
		memset(rh->spectrum.Jgas[0], 0, (size_t) N * Nspace * sizeof(double));
		srand(7);
		for (la = 0; la < N; la++) {
			for (mu = 0; mu < Nrays; mu++) {
				wmu = 0.5 * (mu + 1.0) / 6.0;
				for (to_obs = 0; to_obs <= 1; to_obs++) {
					for (k = 0; k < Nspace; k++) I[k] = 1.0 + urand();
					addJgasRay(rh, la, mu, to_obs, wmu, I);
					if (rest == 0)
						for (k = 0; k < Nspace; k++) J[la][k] += wmu * I[k];
				}
			}
		}
		memcpy((rest ? Jmem : Jtab)[0], rh->spectrum.Jgas[0], (size_t) N * Nspace * sizeof(double));
	}
	freeJgasBuffer();

	CHECK(memcmp(Jtab[0], Jmem[0], (size_t) N * Nspace * sizeof(double)) == 0,
		  "Jgas from the table and with PRDH_LIMIT_MEM differ");
	for (la = 0; la < N; la++)
		for (k = 0; k < Nspace; k += 3) /* locations at rest */
			maxrel = fmax(maxrel, fabs(Jtab[la][k] / J[la][k] - 1.0));
	CHECK(maxrel <= 1e-15, "at rest Jgas differs from J by %.2e (first %g, last %g)", maxrel,
		  Jtab[0][0] / J[0][0], Jtab[N - 1][0] / J[N - 1][0]);

	free(I);
	freeMatrix(rh, (void **) J);
	freeMatrix(rh, (void **) Jmem);
	freeMatrix(rh, (void **) Jtab);
}

static void check_rho(void)
{
	AtomicLine line;
	double     lam[5] = {279.5, 279.55, 279.6, 279.62, 279.7}, rho[5], lag, frac, got;
	int        la, id0, id1, t;

	memset(&line, 0, sizeof(line));
	line.Nlambda = 5;
	line.lambda  = lam;
	for (la = 0; la < 5; la++) rho[la] = 2.0 + 3.0 * lam[la]; /* linear in lambda */

	for (t = 0; t < 200; t++) {
		lag = 279.45 + 0.3 * urand();
		prdhRhoWeights(&line, lag, &id0, &id1, &frac);
		got = (1.0 - frac) * rho[id0] + frac * rho[id1];
		if (lag <= lam[0])
			CHECK(got == rho[0], "below the grid: %g, want %g", got, rho[0]);
		else if (lag >= lam[4])
			CHECK(got == rho[4], "above the grid: %g, want %g", got, rho[4]);
		else
			CHECK_CLOSE(got, 2.0 + 3.0 * lag, 1e-14, 0.0, "inside, at %.6f", lag);
	}
}

int main(void)
{
	srand(1);
	check_brackets();
	check_partition();
	check_collapse();
	check_jgas();
	check_rho();
	return check_summary("test_prdh");
}
