/* ------- file: -------------------------- test_ltepops.c ----------

	   LTE populations in logarithms (ltepops.c, N-7): LTEpops_elem
	   (Saha chain of an element) and LTEpops (Saha-Boltzmann
	   populations of the levels of an atom) against a reference
	   computed in long double (64-bit mantissa, exponent range
	   10^+-4932) in linear space, which neither underflows nor
	   overflows for these parameters.

	   -- Element with 7 stages (C I-VII ionization potentials) and
		  one with 30 stages, T from 3000 K to 10^6 K, electron
		  densities from 10^-10 to 10^24 m^-3 (the low ones only to
		  stretch the chain): the populations are finite and
		  non-negative, sum to abund * nHtot and agree with the
		  reference to 10^-12 relative where it is above 10^-3 of the
		  total, to 10^-12 + 10^-14 |ln(n/total)| above 10^-280 of
		  it (the conditioning of exp), and are negligible below;
	   -- atom with 4 stages and 3 levels each, in the same ranges
		  (n* clipped at DBL_MIN, so the same threshold): the same
		  checks, also with Debeye shielding (sum and finiteness
		  only);
	   -- the old chain (linear, double) overflows or loses the
		  dominant stage in the case documented below, the new code
		  does not.

	   Run with --table to print the cases where the old linear chain
	   fails.
	   --                                              -------------- */

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "check.h"

#define NT     25
#define NNE    18
#define NSPACE (NT * NNE)
#define NPF    3

extern Atmosphere atmos;

static bool_t table;

static const double T_PF[NPF] = {1.0e3, 1.0e5, 1.0e7};

/* --- Ionization potentials of C I-VI [eV] and 30 synthetic stages -- */

static const double CHI_C[6] = {11.2603, 24.3833, 47.8878, 64.4939, 392.087, 489.993};

static void setup_atmos(void)
{
	atmos.Nspace = NSPACE;
	atmos.T      = (double *) malloc(NSPACE * sizeof(double));
	atmos.ne     = (double *) malloc(NSPACE * sizeof(double));
	atmos.nHtot  = (double *) malloc(NSPACE * sizeof(double));
	atmos.Npf    = NPF;
	atmos.Tpf    = (double *) T_PF;

	for (int it = 0; it < NT; it++) {
		double T = 3000.0 * pow(1.0e6 / 3000.0, it / (double) (NT - 1));
		for (int ie = 0; ie < NNE; ie++) {
			int k = it * NNE + ie;

			atmos.T[k]     = T;
			atmos.ne[k]    = pow(10.0, -10.0 + 34.0 * ie / (NNE - 1.0));
			atmos.nHtot[k] = 1.0e18 * (1.0 + 0.1 * k);
		}
	}
}

/* --- ln of the partition function of stage i: linear in T, so that
	   Linear reproduces it exactly -- --------------------------- */

static double lnU(int i, double T)
{
	return log(1.0 + 2.0 * (i % 5)) + 0.01 * (i + 1) * T / 1.0e6;
}

static Element make_element(int Nstage, const double *chi)
{
	Element el;

	memset(&el, 0, sizeof(el));
	strcpy(el.ID, "X");
	el.Nstage = Nstage;
	el.abund  = 3.0e-4;
	el.ionpot = (double *) malloc(Nstage * sizeof(double));
	el.pf     = (double **) malloc(Nstage * sizeof(double *));
	el.n      = (double **) malloc(Nstage * sizeof(double *));
	for (int i = 0; i < Nstage; i++) {
		el.ionpot[i] = chi[i] * EV;
		el.pf[i]     = (double *) malloc(NPF * sizeof(double));
		el.n[i]      = (double *) malloc(NSPACE * sizeof(double));
		for (int j = 0; j < NPF; j++) el.pf[i][j] = lnU(i, T_PF[j]);
	}
	return el;
}

static void free_element(Element *el)
{
	for (int i = 0; i < el->Nstage; i++) {
		free(el->pf[i]);
		free(el->n[i]);
	}
	free(el->pf);
	free(el->n);
	free(el->ionpot);
}

/* --- Reference Saha chain in long double, linear space -- -------- */

static void saha_reference(const Element *el, int k, long double *n)
{
	long double C1 = ((long double) HPLANCK / (2.0L * (long double) PI * (long double) M_ELECTRON)) *
					 ((long double) HPLANCK / (long double) KBOLTZMANN);
	long double T = atmos.T[k], sum = 1.0L;

	n[0] = 1.0L;
	for (int i = 1; i < el->Nstage; i++) {
		long double CT_ne = 2.0L * powl(C1 / T, -1.5L) / (long double) atmos.ne[k];

		n[i] = n[i - 1] * CT_ne *
			   expl((long double) (lnU(i, atmos.T[k]) - lnU(i - 1, atmos.T[k])) -
					(long double) el->ionpot[i - 1] / ((long double) KBOLTZMANN * T));
		sum += n[i];
	}
	for (int i = 0; i < el->Nstage; i++)
		n[i] *= (long double) el->abund * (long double) atmos.nHtot[k] / sum;
}

/* --- The chain as it was before N-7: linear, double -- ----------- */

static void saha_old(const Element *el, int k, double *n)
{
	double C1    = (HPLANCK / (2.0 * PI * M_ELECTRON)) * (HPLANCK / KBOLTZMANN);
	double CT_ne = 2.0 * pow(C1 / atmos.T[k], -1.5) / atmos.ne[k], sum = 1.0;

	n[0] = 1.0;
	for (int i = 1; i < el->Nstage; i++) {
		n[i] = n[i - 1] * CT_ne *
			   exp(lnU(i, atmos.T[k]) - lnU(i - 1, atmos.T[k]) -
				   el->ionpot[i - 1] / (KBOLTZMANN * atmos.T[k]));
		sum += n[i];
	}
	n[0] = el->abund * atmos.nHtot[k] / sum;
	for (int i = 1; i < el->Nstage; i++) n[i] *= n[0];
}

/* --- Compare the populations of all depths with the reference.
	   Returns the largest relative difference of populations above
	   1e-3 of the total -- -------------------------------------- */

static double compare_pops(const char *what, double **got, long double **ref, int Nlev,
						   const double *total)
{
	double worst = 0.0;

	for (int k = 0; k < NSPACE; k++) {
		double sum = 0.0;

		for (int i = 0; i < Nlev; i++) {
			double g = got[i][k], r = (double) ref[i][k];

			CHECK(isfinite(g) && g >= 0.0, "%s: n[%d][%d] = %g not finite or negative", what, i, k, g);
			sum += g;
			if (r > 1.0e-280 * total[k]) {
				/* --- The error of exp(L) is the error of L, which sums terms
					   of up to 10^3: 10^-12 for the populations above 10^-3
					   of the total, growing as 10^-14 |ln(n/total)| --- */

				double rel = fabs(g - r) / r;
				double tol = 1.0e-12 + 1.0e-14 * fabs(log(r / total[k]));

				if (r > 1.0e-3 * total[k] && rel > worst) worst = rel;
				CHECK(rel <= tol, "%s: n[%d][%d] = %.17g, want %.17g (rel. diff %g, T = %g, ne = %g)", what,
					  i, k, g, r, rel, atmos.T[k], atmos.ne[k]);
			} else {
				/* --- Underflowing alone: at most the clip -- */

				CHECK(g <= 1.0e-12 * total[k] + DBL_MIN, "%s: n[%d][%d] = %g should be negligible", what, i, k, g);
			}
		}
		CHECK_CLOSE(sum, total[k], 1.0e-12, 0.0, "%s: sum at k = %d", what, k);
	}
	return worst;
}

static void test_element(const char *what, int Nstage, const double *chi)
{
	Element       el    = make_element(Nstage, chi);
	long double **ref   = (long double **) malloc(Nstage * sizeof(long double *));
	double       *total = (double *) malloc(NSPACE * sizeof(double));
	double       *old   = (double *) malloc(Nstage * sizeof(double)), worst;
	long double  *col   = (long double *) malloc(Nstage * sizeof(long double));
	int           Nold  = 0;

	for (int i = 0; i < Nstage; i++) ref[i] = (long double *) malloc(NSPACE * sizeof(long double));
	LTEpops_elem(&el);
	for (int k = 0; k < NSPACE; k++) {
		saha_reference(&el, k, col);
		for (int i = 0; i < Nstage; i++) {
			CHECK(isfinite((double) col[i]) && col[i] >= 0.0L, "reference at k = %d not usable", k);
			ref[i][k] = col[i];
		}
		total[k] = el.abund * atmos.nHtot[k];

		/* --- Where the old chain is wrong in the dominant stage -- */

		saha_old(&el, k, old);
		for (int i = 0; i < Nstage; i++) {
			double r = (double) col[i];

			if (r > 1.0e-3 * total[k] && !(fabs(old[i] - r) <= 1.0e-9 * r)) {
				if (table && Nold < 6)
					printf("%s old chain wrong: T = %g ne = %g stage %d: old %g, ref %g\n", what, atmos.T[k],
						   atmos.ne[k], i, old[i], r);
				Nold++;
				break;
			}
		}
	}
	worst = compare_pops(what, el.n, ref, Nstage, total);
	if (table) printf("%s: largest relative difference of the dominant stages %.3g; old chain wrong at %d of %d points\n", what, worst, Nold,
					  NSPACE);

	/* --- The old chain must fail somewhere in the 30-stage sweep,
		   which is what this test documents (the 7-stage C chain does
		   not need to) -- --------------------------------------- */

	if (Nstage > 10) CHECK(Nold > 0, "%s: no point where the old chain fails", what);

	for (int i = 0; i < Nstage; i++) free(ref[i]);
	free(ref);
	free(total);
	free(old);
	free(col);
	free_element(&el);
}

/* --- Atom: levels 0-2 in stage 0, 3-5 in stage 1, 6-8 in stage 2,
	   9-11 in stage 3; energies [eV] relative to level 0 -- ------- */

#define NLEV 12

static const double E_LEV[NLEV] = {0.0, 1.3, 2.6, 11.0, 12.5, 14.0, 35.0, 37.0, 39.0, 80.0, 83.0, 90.0};

static void test_atom(bool_t Debeye)
{
	Atom          atom;
	long double **ref   = (long double **) malloc(NLEV * sizeof(long double *));
	double       *total = (double *) malloc(NSPACE * sizeof(double));
	double        worst = 0.0;
	long double   C1    = ((long double) HPLANCK / (2.0L * (long double) PI * (long double) M_ELECTRON)) *
					 ((long double) HPLANCK / (long double) KBOLTZMANN);

	memset(&atom, 0, sizeof(atom));
	strcpy(atom.ID, "X");
	atom.Nlevel = NLEV;
	atom.E      = (double *) malloc(NLEV * sizeof(double));
	atom.g      = (double *) malloc(NLEV * sizeof(double));
	atom.stage  = (int *) malloc(NLEV * sizeof(int));
	atom.ntotal = (double *) malloc(NSPACE * sizeof(double));
	atom.nstar  = (double **) malloc(NLEV * sizeof(double *));
	for (int i = 0; i < NLEV; i++) {
		atom.E[i]     = E_LEV[i] * EV;
		atom.g[i]     = 2.0 * (i % 4) + 1.0;
		atom.stage[i] = i / 3;
		atom.nstar[i] = (double *) malloc(NSPACE * sizeof(double));
		ref[i]        = (long double *) malloc(NSPACE * sizeof(long double));
	}
	for (int k = 0; k < NSPACE; k++) {
		atom.ntotal[k] = 2.0e-5 * atmos.nHtot[k];
		total[k]       = atom.ntotal[k];
	}
	LTEpops(&atom, Debeye);

	if (!Debeye) {
		for (int k = 0; k < NSPACE; k++) {
			long double T = atmos.T[k], cNe_T = 0.5L * atmos.ne[k] * powl(C1 / T, 1.5L), sum = 0.0L;

			for (int i = 0; i < NLEV; i++) {
				ref[i][k] = ((long double) atom.g[i] / atom.g[0]) *
							expl(-((long double) atom.E[i] - atom.E[0]) / ((long double) KBOLTZMANN * T)) /
							powl(cNe_T, atom.stage[i] - atom.stage[0]);
				sum += ref[i][k];
			}
			for (int i = 0; i < NLEV; i++) ref[i][k] *= atom.ntotal[k] / sum;
		}
		worst = compare_pops("atom", atom.nstar, ref, NLEV, total);
		if (table) printf("atom: largest relative difference %.3g\n", worst);
	} else {
		for (int k = 0; k < NSPACE; k++) {
			double sum = 0.0;

			for (int i = 0; i < NLEV; i++) {
				CHECK(isfinite(atom.nstar[i][k]) && atom.nstar[i][k] > 0.0, "Debeye: n*[%d][%d] = %g", i, k,
					  atom.nstar[i][k]);
				sum += atom.nstar[i][k];
			}
			CHECK_CLOSE(sum, total[k], 1.0e-12, 0.0, "Debeye: sum at k = %d", k);
		}
	}

	for (int i = 0; i < NLEV; i++) {
		free(atom.nstar[i]);
		free(ref[i]);
	}
	free(atom.nstar);
	free(ref);
	free(total);
	free(atom.E);
	free(atom.g);
	free(atom.stage);
	free(atom.ntotal);
}

int main(int argc, char *argv[])
{
	static double chi30[30];

	table = (argc > 1 && strcmp(argv[1], "--table") == 0);
	setup_atmos();

	test_element("C I-VII", 7, CHI_C);

	for (int i = 0; i < 30; i++) chi30[i] = 3.0 + 1.5 * i;
	test_element("30 stages", 30, chi30);

	test_atom(FALSE);
	test_atom(TRUE);

	return check_summary("test_ltepops");
}
