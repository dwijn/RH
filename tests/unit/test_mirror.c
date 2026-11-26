/* ------- file: -------------------------- test_mirror.c -----------

	   Integration across the empty blue half of a symmetric line with
	   its red half mirrored (mirrorQuadrature, getlambda.c, N-62).

	   Merged grid: points at 90 - 96 nm, a symmetric line A with its
	   center at 100 nm and red half to 101 nm, and points at 101.5 -
	   110 nm. Nothing lies in A's blue half [99, 100). For a J that is
	   symmetric about 100 nm within A, the mirrored weights on the
	   merged grid must give the trapezoid on the full grid (the merged
	   grid plus the mirror image of A's red half):

	   -- a continuum spanning A, with alpha linear in lambda: the
		  integrand alpha J / lambda;
	   -- a line spanning A with A in its far wing, profile
		  1 / (lambda - lambda0)^2: the integrand phi J;
	   -- with alpha proportional to lambda the weights add up to the
		  length of the range;
	   -- a continuum that ends inside A's red half: the weights still
		  add up to the length of its range;
	   -- a transition entirely redward of A's center, and one that
		  starts at it, keep the trapezoid (wquad NULL), as does A;
	   -- getwlambda_cont and getwlambda_line return the new weights
		  (the latter in Doppler units).
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "inputs.h"
#include "check.h"
#include "state.h"

#define NG    18
#define NFULL (NG + 4)
#define LAM0  100.0

/* --- Merged grid, A's red half at indices 4 - 8 --   -------------- */

static const double grid[NG] = {90.0, 92.0, 94.0, 96.0, 100.0, 100.1,
								100.3, 100.6, 101.0, 101.5, 102.0, 103.0,
								104.0, 105.0, 106.0, 107.0, 108.5, 110.0};

#define IA0    4 /* A's center */
#define NA_RED 5 /* A's red half: 100 .. 101 */

static double Jsym(double lambda)
{
	/* --- Symmetric about LAM0 within A, smooth elsewhere -- -------- */

	double d = fabs(lambda - LAM0);

	if (d <= 1.0) return 1.0 - 0.8 * exp(-SQ(d / 0.2));
	return 1.0 + 0.01 * (lambda - LAM0);
}

static double alphaLin(double lambda)
{
	return 2.0 + 0.05 * lambda;
}

static double profile(double lambda, double lambda0)
{
	return 1.0 / SQ(lambda - lambda0);
}

/* --- The full grid: merged grid plus the mirror image of A's red
	   half, sorted --                                -------------- */

static int fullGrid(double *full)
{
	int la, n = 0;

	for (la = 0; la < IA0; la++) full[n++] = grid[la];
	for (la = NA_RED - 1; la >= 1; la--) full[n++] = 2.0 * LAM0 - grid[IA0 + la];
	for (la = IA0; la < NG; la++) full[n++] = grid[la];
	return n;
}

static double trapezoid(const double *x, const double *f, int n)
{
	double sum = 0.0;
	int    la;

	for (la = 0; la < n - 1; la++) sum += 0.5 * (x[la + 1] - x[la]) * (f[la] + f[la + 1]);
	return sum;
}

int main(void)
{
	Atom            atom, *active[1];
	AtomicLine      lines[2];
	AtomicContinuum conts[4];
	double          full[NFULL], f[NFULL], ref, sum, alpha[NG], alpha_prop[NG], lambda0_T = 120.0;
	int             la, Nfull, Npair;

	rh->commandline.logfile = stderr;
	rh->commandline.quiet   = TRUE;
	rh->atmos.moving        = FALSE;
	rh->atmos.Stokes        = FALSE;

	memset(&atom, 0, sizeof(atom));
	memset(lines, 0, sizeof(lines));
	memset(conts, 0, sizeof(conts));

	/* --- A: symmetric line, center at index IA0 -- -------------- */

	lines[0].symmetric = TRUE;
	lines[0].lambda0   = LAM0;
	lines[0].Nblue     = IA0;
	lines[0].Nlambda   = NA_RED;
	lines[0].lambda    = (double *) grid + IA0;

	/* --- T: asymmetric line over the whole grid, center at 120 nm - */

	lines[1].symmetric = FALSE;
	lines[1].lambda0   = lambda0_T;
	lines[1].Nblue     = 0;
	lines[1].Nlambda   = NG;
	lines[1].lambda    = (double *) grid;

	/* --- Continua: 0 over the whole grid (alpha linear), 1 over the
		   whole grid (alpha proportional to lambda), 2 to 100.3 nm
		   (inside A's red half), 3 from A's center --  -------------- */

	for (la = 0; la < NG; la++) {
		alpha[la]      = alphaLin(grid[la]);
		alpha_prop[la] = 0.03 * grid[la];
	}
	conts[0] = (AtomicContinuum) {.Nblue = 0, .Nlambda = NG, .lambda = (double *) grid, .alpha = alpha};
	conts[1] = (AtomicContinuum) {.Nblue = 0, .Nlambda = NG, .lambda = (double *) grid, .alpha = alpha_prop};
	conts[2] = (AtomicContinuum) {.Nblue = 0, .Nlambda = IA0 + 3, .lambda = (double *) grid, .alpha = alpha_prop};
	conts[3] = (AtomicContinuum) {.Nblue   = IA0,
								  .Nlambda = NG - IA0,
								  .lambda  = (double *) grid + IA0,
								  .alpha   = alpha_prop + IA0};

	strcpy(atom.ID, "XX");
	atom.active    = TRUE;
	atom.Nline     = 2;
	atom.line      = lines;
	atom.Ncont     = 4;
	atom.continuum = conts;
	active[0]      = &atom;

	rh->atmos.Natom       = 1;
	rh->atmos.atoms       = &atom;
	rh->atmos.Nactiveatom = 1;
	rh->atmos.activeatoms = active;
	rh->atmos.Nmolecule   = 0;
	rh->atmos.Nactivemol  = 0;

	Npair = mirrorQuadrature();
	CHECK(Npair == 4, "%d (transition, line) pairs, want 4", Npair);
	CHECK(lines[0].wquad == NULL, "the symmetric line itself has new weights");
	CHECK(lines[1].wquad != NULL, "the line spanning A keeps the trapezoid");
	CHECK(conts[0].wquad != NULL && conts[1].wquad != NULL && conts[2].wquad != NULL,
		  "a continuum spanning A keeps the trapezoid");
	CHECK(conts[3].wquad == NULL, "the continuum starting at A's center has new weights");

	Nfull = fullGrid(full);
	CHECK(Nfull == NFULL, "full grid has %d points, want %d", Nfull, NFULL);

	/* --- Continuum, alpha linear: alpha J / lambda --  -------------- */

	for (la = 0; la < Nfull; la++) f[la] = alphaLin(full[la]) * Jsym(full[la]) / full[la];
	ref = trapezoid(full, f, Nfull);
	sum = 0.0;
	for (la = 0; la < NG; la++)
		sum += getwlambda_cont(&conts[0], la) * alpha[la] * Jsym(grid[la]) / grid[la];
	CHECK_CLOSE(sum, ref, 1e-13, 0.0, "continuum: mirrored %.15g, full grid %.15g", sum, ref);

	/* --- The old trapezoid on the merged grid is far off -- ------- */

	{
		double fm[NG];

		for (la = 0; la < NG; la++) fm[la] = alpha[la] * Jsym(grid[la]) / grid[la];
		CHECK(fabs(trapezoid(grid, fm, NG) / ref - 1.0) > 1e-3,
			  "test does not discriminate: trapezoid on the merged grid is as good");
	}

	/* --- Line T, A in its far wing: phi J, Doppler units -- -------- */

	for (la = 0; la < Nfull; la++) f[la] = profile(full[la], lambda0_T) * Jsym(full[la]);
	ref = trapezoid(full, f, Nfull) * CLIGHT / lambda0_T;
	sum = 0.0;
	for (la = 0; la < NG; la++)
		sum += getwlambda_line(&lines[1], la) * profile(grid[la], lambda0_T) * Jsym(grid[la]);
	CHECK_CLOSE(sum, ref, 1e-13, 0.0, "line: mirrored %.15g, full grid %.15g", sum, ref);

	/* --- Sums of the weights with alpha / lambda constant -- ------- */

	sum = 0.0;
	for (la = 0; la < NG; la++) sum += conts[1].wquad[la];
	CHECK_CLOSE(sum, grid[NG - 1] - grid[0], 1e-13, 0.0, "weights add up to %.15g", sum);

	sum = 0.0;
	for (la = 0; la < conts[2].Nlambda; la++) sum += conts[2].wquad[la];
	CHECK_CLOSE(sum, grid[IA0 + 2] - grid[0], 1e-13, 0.0,
				"weights of the continuum ending inside A add up to %.15g", sum);

	/* --- The trapezoid where nothing changes --     -------------- */

	CHECK(getwlambda_cont(&conts[3], 0) == 0.5 * (grid[IA0 + 1] - grid[IA0]),
		  "continuum from A's center: first weight changed");

	/* --- Calling again gives the same weights --   -------------- */

	{
		double w0[NG];

		memcpy(w0, conts[0].wquad, NG * sizeof(double));
		CHECK(mirrorQuadrature() == Npair, "second call: different number of pairs");
		CHECK(memcmp(w0, conts[0].wquad, NG * sizeof(double)) == 0,
			  "second call: different weights");
	}

	free(lines[1].wquad);
	for (la = 0; la < 4; la++) free(conts[la].wquad);

	return check_summary("test_mirror");
}
