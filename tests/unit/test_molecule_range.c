/* ------- file: -------------------------- test_molecule_range.c ---

	   Partition function and equilibrium constant of a molecule outside
	   the range [Tmin, Tmax] of their fits (chemequil.c, N-53). The
	   fits used to return 0 there, so that every molecule vanished
	   below Tmin (1000 K for all shipped molecules) and above Tmax.
	   Now:

	   -- the partition function is held at its value at the nearest
		  edge;
	   -- ln K is extrapolated linearly in theta = 5040/T from the
		  nearer edge, with the value and slope of the fit at the edge:
		  K is continuous at the edge, its slope d ln K / d theta is
		  continuous there, ln K is a straight line in theta outside
		  (three points checked), and K is positive and finite;
	   -- inside the range the result is that of the fit (compared
		  with ln K of the log form that the extrapolation uses, by
		  continuity at the edge, which also checks its units).

	   Molecules: the coefficients of the shipped H2 and CO (KURUCZ_85),
	   H2O and H2+ (KURUCZ_70, the latter charged), OH and CH
	   (SAUVAL_TATUM_84) and a synthetic IRWIN_81 and TSUJI_73. Run with
	   --table to print, for the shipped ones, the extrapolated K over
	   the last 200 and 500 K of the fit range against the fit.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "check.h"

typedef struct {
	const char   *id;
	enum fit_type fit;
	int           charge, Nnuclei;
	double        Ediss_eV, Tmin, Tmax;
	int           Npf, Neqc;
	double        pf[8], eqc[8]; /* as in the files: lowest power first */
	int           shipped;
} Spec;

static const Spec specs[] = {
	{"H2", KURUCZ_85, 0, 2, 4.478, 1000, 9000, 7, 7, {0.582145, 16.3760, -49.4684, 112.049, -149.953, 106.531, -30.9791}, {-46.4584, 16.3660, -49.3992, 111.822, -149.567, 106.206, -30.8720}, 1},
	{"CO", KURUCZ_85, 0, 2, 11.091, 1000, 12000, 7, 6, {4.51349, 18.4221, -50.0599, 102.208, -128.504, 87.8414, -24.8533}, {-49.0414, 14.0306, -26.6341, 35.3827, -26.5424, 8.32385}, 1},
	{"H2O", KURUCZ_70, 0, 3, 9.512, 1000, 6000, 0, 5, {0}, {-9.3179E+01, 2.6725E-03, -5.7830E-07, 8.5268E-11, -5.1311E-15}, 1},
	{"H2+", KURUCZ_70, 1, 2, -10.948, 1000, 9000, 0, 5, {0}, {-9.2120E+00, 1.9386E-03, -4.7676E-07, 7.4004E-11, -4.5518E-15}, 1},
	{"OH", SAUVAL_TATUM_84, 0, 2, 4.392, 1000, 9000, 4, 3, {3.0929, -1.6778, 0.6743, -0.1874}, {10.7881, -0.8762, -0.5436}, 1},
	{"CH", SAUVAL_TATUM_84, 0, 2, 3.465, 1000, 9000, 4, 3, {3.3586, -2.0656, 0.9624, -0.2239}, {10.5263, -0.5284, -0.7239}, 1},
	/* synthetic: IRWIN_81 (ln T polynomial for the partition function,
	   Sauval-Tatum form for K) and TSUJI_73 (K from a polynomial in
	   theta; the partition function is not used) */
	{"XIRW", IRWIN_81, 0, 2, 5.0, 500, 7000, 3, 3, {-3.0, 1.2, -0.04}, {10.5, -0.8, -0.5}, 0},
	{"XTSU", TSUJI_73, 0, 2, 5.0, 1000, 9000, 0, 3, {0}, {-9.0, 3.0, 0.1}, 0},
};
#define NSPEC ((int) (sizeof(specs) / sizeof(specs[0])))

static void build(const Spec *s, Molecule *m)
{
	int n;

	memset(m, 0, sizeof(*m));
	strcpy(m->ID, s->id);
	m->fit     = s->fit;
	m->charge  = s->charge;
	m->Nnuclei = s->Nnuclei;
	m->Ediss   = s->Ediss_eV * EV;
	m->Tmin    = s->Tmin;
	m->Tmax    = s->Tmax;
	m->Npf     = s->Npf;
	m->Neqc    = s->Neqc;

	/* --- readMolecule stores the coefficients highest power first -- */

	m->pf_coef  = (double *) malloc(8 * sizeof(double));
	m->eqc_coef = (double *) malloc(8 * sizeof(double));
	for (n = 0; n < s->Npf; n++) m->pf_coef[s->Npf - 1 - n] = s->pf[n];
	for (n = 0; n < s->Neqc; n++) m->eqc_coef[s->Neqc - 1 - n] = s->eqc[n];
}

static double lnK(Molecule *m, double theta)
{
	return log(equilconstant(rh, m, THETA0 / theta));
}

int main(int argc, char *argv[])
{
	int      i, e, table = (argc > 1 && !strcmp(argv[1], "--table"));
	Molecule m;

	for (i = 0; i < NSPEC; i++) {
		const Spec *s = &specs[i];

		build(s, &m);

		for (e = 0; e < 2; e++) {
			double Tedge = e ? m.Tmax : m.Tmin;
			double sgn   = e ? 1.0 : -1.0; /* outside: T = Tedge + sgn * dT */
			double th0   = THETA0 / Tedge;
			double lo, hi, sl_in, sl_out, d = 1.0E-5 * th0;
			double th1, th2, th3, l1, l2, l3, slope;

			/* --- Continuity at the edge (this also checks the units
				   of the log form against the product form) --  --- */

			CHECK_CLOSE(equilconstant(rh, &m, Tedge + sgn * 1.0E-6),
						equilconstant(rh, &m, Tedge), 1.0E-6, 0.0,
						"%s: K continuous at %s", s->id, e ? "Tmax" : "Tmin");

			/* --- Slope in theta continuous across the edge: one-sided
				   differences of ln K over 1e-5 of theta -- */

			lo = lnK(&m, th0 - d);
			hi = lnK(&m, th0 + d);
			if (e) { /* high T = small theta: inside is larger theta */
				sl_in  = (hi - lnK(&m, th0)) / d;
				sl_out = (lnK(&m, th0) - lo) / d;
			} else {
				sl_in  = (lnK(&m, th0) - lo) / d;
				sl_out = (hi - lnK(&m, th0)) / d;
			}
			CHECK_CLOSE(sl_out, sl_in, 1.0E-3, 0.0,
						"%s: d lnK / d theta continuous at %s", s->id,
						e ? "Tmax" : "Tmin");

			/* --- A straight line in theta outside --      ------- */

			th1   = th0 * (e ? 0.8 : 1.3);
			th2   = th0 * (e ? 0.5 : 2.0);
			th3   = th0 * (e ? 0.2 : 5.0);
			l1    = lnK(&m, th1);
			l2    = lnK(&m, th2);
			l3    = lnK(&m, th3);
			slope = (l2 - l1) / (th2 - th1);
			CHECK_CLOSE(l3, l2 + slope * (th3 - th2), 1.0E-10, 1.0E-10,
						"%s: ln K linear in theta beyond %s", s->id,
						e ? "Tmax" : "Tmin");
			CHECK(isfinite(l1) && isfinite(l2) && isfinite(l3) &&
					  equilconstant(rh, &m, THETA0 / th3) > 0.0,
				  "%s: K positive and finite beyond %s", s->id,
				  e ? "Tmax" : "Tmin");

			/* --- Dissociation: K grows towards low T (e = 0) --  */

			if (!e && s->shipped && m.Ediss > 0.0)
				CHECK(l3 > l2 && l2 > l1, "%s: K increases towards low T", s->id);
		}
		/* --- Partition function held at the edges --      ------- */

		if (s->Npf > 0 && s->fit != TSUJI_73) {
			double pmin = partfunction(rh, &m, m.Tmin);
			double pmax = partfunction(rh, &m, m.Tmax);

			CHECK(pmin > 0.0 && pmax > 0.0, "%s: pf positive at the edges", s->id);
			CHECK(partfunction(rh, &m, 0.5 * m.Tmin) == pmin &&
					  partfunction(rh, &m, 300.0) == pmin,
				  "%s: pf held at Tmin below it", s->id);
			CHECK(partfunction(rh, &m, 3.0 * m.Tmax) == pmax &&
					  partfunction(rh, &m, 1.0E6) == pmax,
				  "%s: pf held at Tmax above it", s->id);
		}
		if (table && s->shipped) {
			/* Accuracy of the extrapolation: raise Tmin of a copy to
			   1000 + dT and compare its K at 1000 K (extrapolated)
			   with the fit, then lower Tmax of a copy by dT and the
			   same at the original Tmax. */
			static const double dTs[] = {200.0, 500.0};
			int                 j;

			for (j = 0; j < 2; j++) {
				Molecule c = m;

				c.Tmin = m.Tmin + dTs[j];
				printf("%-4s Tmin raised to %5.0f K: K(%4.0f K) extrapolated / fit"
					   " = %.5f\n",
					   s->id, c.Tmin, m.Tmin,
					   equilconstant(rh, &c, m.Tmin) / equilconstant(rh, &m, m.Tmin));
				c      = m;
				c.Tmax = m.Tmax - dTs[j];
				printf("%-4s Tmax lowered to %5.0f K: K(%4.0f K) extrapolated / fit"
					   " = %.5f\n",
					   s->id, c.Tmax, m.Tmax,
					   equilconstant(rh, &c, m.Tmax) / equilconstant(rh, &m, m.Tmax));
			}
		}
		free(m.pf_coef);
		free(m.eqc_coef);
	}
	return check_summary("test_molecule_range");
}
