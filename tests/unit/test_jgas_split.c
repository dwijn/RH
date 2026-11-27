/* ------- file: -------------------------- test_jgas_split.c --------

	   Jgas in the PRD sub-iterations of hybrid PRD (readj.c, N-27).

	   A sub-iteration re-solves only the wavelengths with a PRD line,
	   but Jgas at those wavelengths also holds contributions of
	   neighbouring wavelengths without one. With beginJgas, endJgas and
	   the split of addJgasRay between Jgas and Jgas_fixed:

	   1. After a sub-iteration (only the wavelengths with PRD line
		  re-solved, with new intensities) Jgas at the wavelengths with a
		  PRD line equals Jgas of a formal solution of all wavelengths with
		  the same intensities, to rounding.
	   2. Without the split (zeroing Jgas, as before N-27) the two
		  differ, so the test sees the bug.
	   3. Several sub-iterations in a row keep it so.
	   4. initJgasRows gives a row of Jgas_fixed only to the PRD
		  wavelengths a wavelength without PRD line can reach (some do
		  not here); a contribution to one that has none stops the run
		  (so the test passing shows the reach holds all contributions).
	   5. Without Jgas_fixed (spectrum.Jgas_fixed == NULL) all goes into
		  Jgas as before.
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

#define NRAYS  2
#define NSPACE 7
#define NWAVE  240

static char isPRD[NWAVE];

static double urand(void)
{
	return rand() / (RAND_MAX + 1.0);
}

/* --- Intensity of a ray: fixed for the wavelengths without PRD line,
	   changing with the generation (the rho of the iteration) for those
	   with --                                         -------------- */

static double Ival(int la, int mu, int to_obs, int k, int gen)
{
	double base = 1.0 + 0.5 * sin(0.37 * la + 1.3 * mu + 2.1 * to_obs + 0.7 * k);

	return (isPRD[la]) ? base * (1.0 + 0.1 * (gen + 1) * cos(0.11 * la + k)) : base;
}

static void sweep(bool_t only_prd, int gen)
{
	int    la, mu, to_obs, k;
	double wmu, I[NSPACE];

	for (la = 0; la < NWAVE; la++) {
		if (only_prd && !isPRD[la]) continue;
		for (mu = 0; mu < NRAYS; mu++) {
			wmu = 0.5 * (mu + 1.0) / 3.0;
			for (to_obs = 0; to_obs <= 1; to_obs++) {
				for (k = 0; k < NSPACE; k++) I[k] = Ival(la, mu, to_obs, k, gen);
				addJgasRay(rh, la, mu, to_obs, wmu, I);
			}
		}
	}
}

static void full(int gen)
{
	beginJgas(rh, FALSE);
	sweep(FALSE, gen);
	endJgas(rh, FALSE);
}

static void subiter(int gen)
{
	beginJgas(rh, TRUE);
	sweep(TRUE, gen);
	endJgas(rh, TRUE);
}

static double maxrel_prd(double **A, double **B)
{
	int    la, k;
	double m = 0.0;

	for (la = 0; la < NWAVE; la++) {
		if (!isPRD[la]) continue;
		for (k = 0; k < NSPACE; k++)
			m = fmax(m, fabs(A[la][k] - B[la][k]) / fmax(fabs(B[la][k]), 1e-300));
	}
	return m;
}

int main(void)
{
	int    la, k, mu, gen, nprd = 0;
	double l = 250.0, **ref, err, naive;
	size_t n = (size_t) NWAVE * NSPACE * sizeof(double);

	srand(3);
	rh->spectrum.Nspect = NWAVE;
	rh->spectrum.lambda = (double *) malloc(NWAVE * sizeof(double));
	for (la = 0; la < NWAVE; la++) {
		double u = urand();

		l += (u < 0.6) ? 1e-4 * (0.5 + urand()) : (u < 0.95) ? 2e-3 * urand() + 1e-6
															 : 0.05;
		rh->spectrum.lambda[la] = l;
	}
	rh->atmos.Nrays       = NRAYS;
	rh->atmos.Nspace      = NSPACE;
	rh->atmos.Nactiveatom = 0;
	rh->spectrum.v_los    = matrix_double(rh, NRAYS, NSPACE);
	for (mu = 0; mu < NRAYS; mu++)
		for (k = 0; k < NSPACE; k++) rh->spectrum.v_los[mu][k] = 2.0e4 * (2.0 * urand() - 1.0);

	/* --- Every third wavelength has a PRD line, in two groups, and
		   a block of 70 in a row, whose inner part no wavelength without
		   PRD line reaches --                           -------------- */

	for (la = 0; la < NWAVE; la++) {
		isPRD[la] = ((la / 7) % 3 == 0 && la % 2 == 0) || (la >= 150 && la < 220);
		nprd += isPRD[la];
	}
	CHECK(nprd > 20 && nprd < NWAVE - 20, "%d wavelengths with a PRD line", nprd);

	rh->spectrum.Jgas        = matrix_double(rh, NWAVE, NSPACE);
	rh->input.prdh_limit_mem = FALSE;
	initPRDhTables(rh);
	initJgasRows(rh, isPRD);
	{
		int nrow = 0, nnone = 0;

		for (la = 0; la < NWAVE; la++) {
			CHECK((rh->spectrum.PRDrow[la] == -1) == !isPRD[la], "row of wavelength %d", la);
			nrow += rh->spectrum.PRDrow[la] >= 0;
			nnone += rh->spectrum.PRDrow[la] == -2;
		}
		CHECK(rh->spectrum.Jgas_fixed != NULL && nrow > 20 && nnone > 10,
			  "%d PRD wavelengths with a row of Jgas_fixed, %d without", nrow, nnone);
	}
	initJgasLocks(rh);
	ref = matrix_double(rh, NWAVE, NSPACE);

	/* --- 1 and 3: sub-iterations with changing intensities -- ------ */

	full(0);
	for (gen = 1; gen <= 3; gen++) {
		subiter(gen);
		memcpy(ref[0], rh->spectrum.Jgas[0], n);

		/* --- The full solution of the same intensities, from a
			   clean slate; then back to the sub-iteration's state -- */

		full(gen);
		err = maxrel_prd(ref, rh->spectrum.Jgas);
		CHECK(err <= 1e-13, "generation %d: sub-iteration Jgas differs from the full solution by %.2e",
			  gen, err);
	}

	/* --- 2: the old behaviour, Jgas zeroed and only the wavelengths with
		   PRD line added, differs from the full solution -- ------- */

	full(0);
	{
		double **save = rh->spectrum.Jgas_fixed;

		rh->spectrum.Jgas_fixed = NULL;
		memset(rh->spectrum.Jgas[0], 0, n);
		sweep(TRUE, 1);
		rh->spectrum.Jgas_fixed = save;
	}
	memcpy(ref[0], rh->spectrum.Jgas[0], n);
	full(1);
	naive = maxrel_prd(ref, rh->spectrum.Jgas);
	CHECK(naive > 1e-3, "zeroing Jgas differs from the full solution by only %.2e", naive);

	/* --- 4: without Jgas_fixed everything is in Jgas --------------- */

	{
		double **save = rh->spectrum.Jgas_fixed;

		full(2);
		memcpy(ref[0], rh->spectrum.Jgas[0], n);
		rh->spectrum.Jgas_fixed = NULL;
		beginJgas(rh, FALSE);
		sweep(FALSE, 2);
		endJgas(rh, FALSE);
		rh->spectrum.Jgas_fixed = save;
		err                     = maxrel_prd(ref, rh->spectrum.Jgas);
		CHECK(err <= 1e-13, "without Jgas_fixed Jgas differs by %.2e", err);
	}

	freeJgasBuffer();
	freeMatrix(rh, (void **) ref);
	return check_summary("test_jgas_split");
}
