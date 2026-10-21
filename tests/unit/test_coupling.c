/* ------- file: -------------------------- test_coupling.c ---------

	   Preconditioned radiative rates (addPreconditionedRates in
	   fillgamma.c, B-6).

	   For the active transitions t = (l, u) at one wavelength and angle,
	   with V_t, g_t, wla_t, U_t = (2h nu^3/c^2) g_t V_t, and Ubar_a the
	   sum of U_t over the transitions with upper level a, the kernel
	   adds to Gamma the rates of t with

		 Ieff_t = I - Psi (n_u^old Ubar_u + n_l^old Ubar_l)

	   and the cross-coupling terms that put the bracket back with the
	   new populations, chi_t^old Psi (n_u Ubar_u + n_l Ubar_l). statEquil
	   builds the diagonal from the column sums, so Gamma represents the
	   net rate into level c as

		 sum_{a != c} Gamma[c][a] n_a - n_c sum_{b != c} Gamma[b][c].

	   For random sets of transitions, including the overlaps that the
	   old form got wrong (shared lower level, ladder, continua into an
	   ion level that is the lower level of an ion line), the test
	   checks:

	   -- at the old populations the net rates are the plain radiative
		  rates with I, to rounding: the converged populations solve
		  statistical equilibrium (fixed-point consistency);
	   -- at other populations they are the preconditioned rates above,
		  i.e. linear in the new populations as derived;
	   -- Gamma has entries only between the two levels of a transition;
	   -- the per-point angle weights (wk, spherical geometry) give the
		  same result as the scalar weight.

	   Run with --table to print the measured errors.
	   --                                              -------------- */

#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "check.h"

#define NLEVEL 8
#define NSPACE 5
#define MAXT   (NLEVEL * (NLEVEL - 1) / 2)

static bool_t table;

/* --- Deterministic pseudo-random numbers in [-1, 1) -- ------------ */

static unsigned long long rng_state = 88172645463325252ULL;

static double uniform(void)
{
	rng_state ^= rng_state << 13;
	rng_state ^= rng_state >> 7;
	rng_state ^= rng_state << 17;
	return 2.0 * (double) (rng_state >> 11) / 9007199254740992.0 - 1.0;
}

/* --- Net rate into each level represented by Gamma -- ------------- */

static void net_from_gamma(double **Gam, double n[NLEVEL][NSPACE], int k,
						   double *net)
{
	for (int c = 0; c < NLEVEL; c++) {
		net[c] = 0.0;
		for (int a = 0; a < NLEVEL; a++) {
			if (a == c) continue;
			net[c] += Gam[c * NLEVEL + a][k] * n[a][k];
			net[c] -= Gam[a * NLEVEL + c][k] * n[c][k];
		}
	}
}

/* --- Check one active set of Nt transitions (lower[t], upper[t]) -- */

static void check_set(const char *name, int Nt, const int *tl, const int *tu)
{
	double V[MAXT][NSPACE], wla[MAXT][NSPACE], gij[MAXT][NSPACE], c2[MAXT];
	double n0[NLEVEL][NSPACE], n1[NLEVEL][NSPACE], I[NSPACE], Psi[NSPACE],
		wk[NSPACE];
	double **U   = matrix_double(NLEVEL, NSPACE),
		   **Gam = matrix_double(NLEVEL * NLEVEL, NSPACE),
		   **Gw  = matrix_double(NLEVEL * NLEVEL, NSPACE);
	double wmu = 0.37, err_fix = 0.0, err_lin = 0.0;
	bool_t connected[NLEVEL][NLEVEL] = {{FALSE}};

	for (int t = 0; t < Nt; t++) {
		c2[t]                   = 0.5 * (1.5 + uniform());
		connected[tl[t]][tu[t]] = connected[tu[t]][tl[t]] = TRUE;
		for (int k = 0; k < NSPACE; k++) {
			V[t][k]   = 0.5 * (1.0 + uniform());
			wla[t][k] = 0.5 * (1.0 + uniform());
			gij[t][k] = 0.6 + 0.5 * uniform(); /* inversions occur */
		}
	}
	for (int k = 0; k < NSPACE; k++) {
		I[k]   = 0.5 * (1.0 + uniform());
		Psi[k] = 0.5 * (1.0 + uniform());
		wk[k]  = 0.5 * (1.0 + uniform());
		for (int c = 0; c < NLEVEL; c++) {
			n0[c][k] = 1.0 + uniform();
			n1[c][k] = 1.0 + uniform();
		}
	}

	/* --- Ubar as addtoCoupling accumulates it -- -------------------- */

	for (int c = 0; c < NLEVEL; c++) {
		for (int k = 0; k < NSPACE; k++) U[c][k] = 0.0;
	}
	for (int t = 0; t < Nt; t++) {
		for (int k = 0; k < NSPACE; k++) U[tu[t]][k] += c2[t] * gij[t][k] * V[t][k];
	}

	/* --- Gamma, with scalar weight wmu and with wmu = 1 times wk -- - */

	for (int ij = 0; ij < NLEVEL * NLEVEL; ij++) {
		for (int k = 0; k < NSPACE; k++) Gam[ij][k] = Gw[ij][k] = 0.0;
	}
	for (int t = 0; t < Nt; t++) {
		addPreconditionedRates(NLEVEL, tl[t], tu[t], c2[t], V[t], wla[t], gij[t],
							   n0[tl[t]], n0[tu[t]], U[tl[t]], U[tu[t]], I, Psi,
							   wmu, NULL, NSPACE, Gam);
		addPreconditionedRates(NLEVEL, tl[t], tu[t], c2[t], V[t], wla[t], gij[t],
							   n0[tl[t]], n0[tu[t]], U[tl[t]], U[tu[t]], I, Psi,
							   1.0, wk, NSPACE, Gw);
	}

	for (int k = 0; k < NSPACE; k++) {
		double net[NLEVEL], plain[NLEVEL] = {0.0}, lin[NLEVEL] = {0.0},
							scale = 0.0;

		/* --- Plain rates with I at the old populations, and the
			   preconditioned rates at the new ones --     -------------- */

		for (int t = 0; t < Nt; t++) {
			int    l = tl[t], u = tu[t];
			double w = V[t][k] * wla[t][k] * wmu, F, Ie, chi;

			F = w * (n0[l][k] * I[k] - n0[u][k] * gij[t][k] * (c2[t] + I[k]));
			plain[u] += F;
			plain[l] -= F;
			scale += fabs(w * n0[l][k] * I[k]) +
					 fabs(w * n0[u][k] * gij[t][k] * (c2[t] + I[k]));

			Ie  = I[k] - Psi[k] * (n0[u][k] * U[u][k] + n0[l][k] * U[l][k]);
			chi = V[t][k] * wla[t][k] * (n0[l][k] - gij[t][k] * n0[u][k]);
			F   = w * (n1[l][k] * Ie - n1[u][k] * gij[t][k] * (c2[t] + Ie)) +
				chi * Psi[k] * wmu * (n1[u][k] * U[u][k] + n1[l][k] * U[l][k]);
			lin[u] += F;
			lin[l] -= F;
		}
		net_from_gamma(Gam, n0, k, net);
		for (int c = 0; c < NLEVEL; c++) {
			CHECK_CLOSE(net[c], plain[c], 0.0, 1e-14 * scale,
						"%s: fixed point, net rate into level %d at k = %d",
						name, c, k);
			err_fix = MAX(err_fix, fabs(net[c] - plain[c]) / scale);
		}
		net_from_gamma(Gam, n1, k, net);
		for (int c = 0; c < NLEVEL; c++) {
			CHECK_CLOSE(net[c], lin[c], 0.0, 1e-14 * scale,
						"%s: linear form, net rate into level %d at k = %d",
						name, c, k);
			err_lin = MAX(err_lin, fabs(net[c] - lin[c]) / scale);
		}
		for (int c = 0; c < NLEVEL; c++) {
			for (int a = 0; a < NLEVEL; a++) {
				int ca = c * NLEVEL + a;

				if (a != c && !connected[c][a]) {
					CHECK(Gam[ca][k] == 0.0,
						  "%s: Gamma[%d][%d] = %g at k = %d, no transition",
						  name, c, a, Gam[ca][k], k);
				}
				CHECK_CLOSE(Gw[ca][k] * wmu, Gam[ca][k] * wk[k], 0.0, 1e-14 * scale,
							"%s: per-point weight, Gamma[%d][%d] at k = %d", name,
							c, a, k);
			}
		}
	}
	if (table) {
		printf("  %-24s %2d transitions  fixed point %.1e  linear %.1e\n", name, Nt,
			   err_fix, err_lin);
	}
	freeMatrix((void **) U);
	freeMatrix((void **) Gam);
	freeMatrix((void **) Gw);
}

int main(int argc, char *argv[])
{
	table = (argc > 1 && !strcmp(argv[1], "--table"));

	{
		const int l[] = {0, 1}, u[] = {3, 3};
		check_set("shared upper level", 2, l, u);
	}
	{
		const int l[] = {0, 0}, u[] = {3, 4};
		check_set("shared lower level", 2, l, u);
	}
	{
		const int l[] = {0, 1}, u[] = {1, 6};
		check_set("ladder", 2, l, u);
	}
	{
		const int l[] = {0, 1, 2, 6}, u[] = {6, 6, 6, 7};
		check_set("continua into ion line", 4, l, u);
	}
	for (int s = 0; s < 20; s++) {
		int  tl[MAXT], tu[MAXT], Nt = 0;
		char name[32];

		for (int i = 0; i < NLEVEL; i++) {
			for (int j = i + 1; j < NLEVEL; j++) {
				if (uniform() > 0.4) {
					tl[Nt]   = i;
					tu[Nt++] = j;
				}
			}
		}
		if (Nt == 0) continue;
		sprintf(name, "random set %d", s);
		check_set(name, Nt, tl, tu);
	}
	return check_summary("test_coupling");
}
