/* ------- file: -------------------------- pchip.c -----------------

	   Monotone piecewise cubic Hermite interpolation (PCHIP) of tabulated
	   data, optionally in log x and log y (N-32).

	   --------------------------                      ----------RH-- */

/* --- Tabulated atomic data (collision strengths and rate coefficients
	   against temperature, photoionization cross-sections against
	   wavelength) must stay positive between the nodes. A cubic spline
	   does not: it is not local, so one bad or steep value rings through
	   the whole table, and it overshoots wherever the data turn sharply,
	   which makes rates and opacities negative (N-32).

	   PCHIP (Fritsch & Carlson 1980, SIAM J. Numer. Anal. 17, 238) is a
	   cubic Hermite interpolant whose slopes are chosen so that on every
	   interval it stays between the two node values: it adds no extrema,
	   is monotone where the data are, and a node affects only the
	   intervals within two nodes of it. The slopes at interior nodes are
	   the weighted harmonic mean of the neighbouring secants (Fritsch &
	   Butland 1984, SIAM J. Sci. Stat. Comput. 5, 300), zero where the
	   secants change sign; the end slopes come from the one-sided
	   three-point formula, limited so that they keep the sign of the end
	   secant and are at most three times it. This is the pchip of
	   Moler (2004), Numerical Computing with MATLAB, section 3.6, and of
	   SciPy's PchipInterpolator. The interpolant is C1.

	   Collision rates and cross-sections vary roughly as powers of
	   temperature and wavelength, so the tables are interpolated in log
	   x and log y (PCHIP_LOGLOG): a power law between nodes is then
	   reproduced exactly, and log y keeps the result positive whatever
	   the slopes. A value <= 0 cannot be interpolated in log y: the
	   intervals that touch one are interpolated in y, where PCHIP still
	   keeps every interval between its two node values, so
	   non-negative data stay non-negative, and every run of positive
	   nodes stays in log y as a table of its own (N-49; before, one
	   such value put the whole table in y). Outside the table the end
	   values are used, as by Linear and the spline routines.

	   The table is copied (transformed, ascending), so the caller's
	   arrays need not outlive the PchipTable, and evaluation keeps no
	   state, so one table can be evaluated from several threads.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "error.h"

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- endSlope.c -------------- */

/* --- One-sided three-point slope at an end of the table, limited to
	   keep the shape (Moler 2004) --                  -------------- */

static double endSlope(double h0, double h1, double del0, double del1)
{
	double d = ((2.0 * h0 + h1) * del0 - h0 * del1) / (h0 + h1);

	if (d * del0 <= 0.0)
		d = 0.0;
	else if (del0 * del1 < 0.0 && fabs(d) > fabs(3.0 * del0))
		d = 3.0 * del0;

	return d;
}
/* ------- end ---------------------------- endSlope.c -------------- */

/* ------- begin -------------------------- pchipSlopes.c ----------- */

/* --- The PCHIP slopes d of the table x, y (N points, x ascending):
	   the weighted harmonic mean of the secants inside, zero where they
	   change sign, endSlope at the ends --          -------------- */

static void pchipSlopes(int N, const double *x, const double *y, double *d)
{
	register int k;

	double *h, *del, w1, w2;

	if (N < 2) {
		if (N == 1) d[0] = 0.0;
		return;
	}
	h   = (double *) malloc((N - 1) * sizeof(double));
	del = (double *) malloc((N - 1) * sizeof(double));
	for (k = 0; k < N - 1; k++) {
		h[k]   = x[k + 1] - x[k];
		del[k] = (y[k + 1] - y[k]) / h[k];
	}
	if (N == 2) {
		d[0] = d[1] = del[0];
	} else {
		for (k = 1; k < N - 1; k++) {
			if (del[k - 1] * del[k] > 0.0) {
				w1   = 2.0 * h[k] + h[k - 1];
				w2   = h[k] + 2.0 * h[k - 1];
				d[k] = (w1 + w2) / (w1 / del[k - 1] + w2 / del[k]);
			} else
				d[k] = 0.0;
		}
		d[0]     = endSlope(h[0], h[1], del[0], del[1]);
		d[N - 1] = endSlope(h[N - 2], h[N - 3], del[N - 2], del[N - 3]);
	}
	free(h);
	free(del);
}
/* ------- end ---------------------------- pchipSlopes.c ----------- */

/* ------- begin -------------------------- pchipTableInit.c -------- */

void pchipTableInit(PchipTable *pt, int N, const double *x, const double *y,
					int scale)
{
	const char   routineName[] = "pchipTableInit";
	register int k;

	bool_t ascend, mixed = FALSE;
	int    kk, a;

	if (N < 1) {
		Error(ERROR_LEVEL_2, routineName, "Table needs at least one point (N = %d)", N);
	}
	pt->N    = N;
	pt->logx = (scale & PCHIP_LOGX) != 0;
	pt->logy = (scale & PCHIP_LOGY) != 0;
	if (pt->logy) {
		for (k = 0; k < N; k++)
			if (!(y[k] > 0.0)) mixed = TRUE;
	}
	pt->x    = (double *) malloc(N * sizeof(double));
	pt->y    = (double *) malloc(N * sizeof(double));
	pt->d    = (double *) calloc(N, sizeof(double));
	pt->ylin = (mixed) ? (double *) malloc(N * sizeof(double)) : NULL;
	pt->dlin = (mixed) ? (double *) calloc(N, sizeof(double)) : NULL;

	/* --- Store ascending in x, transformed. With log y and a node
		   <= 0 (N-49), y is kept in log y at the positive nodes and as
		   it is in ylin --                               -------------- */

	ascend = (N == 1 || x[1] > x[0]);
	for (k = 0; k < N; k++) {
		kk = (ascend) ? k : N - 1 - k;
		if (pt->logx && !(x[kk] > 0.0)) {
			Error(ERROR_LEVEL_2, routineName, "Interpolation in log x needs x > 0 (x = %g)",
				  x[kk]);
		}
		pt->x[k] = (pt->logx) ? log(x[kk]) : x[kk];
		pt->y[k] = (pt->logy && y[kk] > 0.0) ? log(y[kk]) : y[kk];
		if (mixed) pt->ylin[k] = y[kk];
	}
	for (k = 1; k < N; k++) {
		if (!(pt->x[k] > pt->x[k - 1])) {
			Error(ERROR_LEVEL_2, routineName, "Table abscissae are not strictly monotonic"
											  " (point %d of %d)",
				  k, N);
		}
	}
	if (!mixed) {
		pchipSlopes(N, pt->x, pt->y, pt->d);
		return;
	}

	/* --- Mixed (N-49): the slopes in y over the whole table, for the
		   intervals that touch a node <= 0, and in log y over each run
		   of positive nodes, as a table of its own, for the intervals
		   within it --                                   -------------- */

	pchipSlopes(N, pt->x, pt->ylin, pt->dlin);
	for (a = 0; a < N;) {
		if (!(pt->ylin[a] > 0.0)) {
			a++;
			continue;
		}
		for (k = a; k + 1 < N && pt->ylin[k + 1] > 0.0; k++);
		pchipSlopes(k - a + 1, pt->x + a, pt->y + a, pt->d + a);
		a = k + 1;
	}
}
/* ------- end ---------------------------- pchipTableInit.c -------- */

/* ------- begin -------------------------- pchipTableEval.c -------- */

void pchipTableEval(const PchipTable *pt, int N, const double *x, double *y)
{
	register int n;

	int           lo, hi, mid, N1 = pt->N - 1;
	bool_t        logv;
	double        t, h, s, s2, s3, v;
	const double *yv, *dv;

	for (n = 0; n < N; n++) {
		t = (pt->logx) ? ((x[n] > 0.0) ? log(x[n]) : -HUGE_VAL) : x[n];

		/* --- Outside the table the end values; with log y and a node
			   <= 0, each interval in log y if both its nodes are
			   positive, else in y (N-49) --              -------------- */

		if (t <= pt->x[0]) {
			y[n] = (pt->ylin) ? pt->ylin[0] : ((pt->logy) ? exp(pt->y[0]) : pt->y[0]);
			continue;
		} else if (t >= pt->x[N1]) {
			y[n] = (pt->ylin) ? pt->ylin[N1] : ((pt->logy) ? exp(pt->y[N1]) : pt->y[N1]);
			continue;
		}

		/* --- Interval by bisection: x[lo] <= t < x[lo + 1] -- ---- */

		lo = 0;
		hi = N1;
		while (hi - lo > 1) {
			mid = (lo + hi) / 2;
			if (pt->x[mid] > t)
				hi = mid;
			else
				lo = mid;
		}
		logv = pt->logy && (!pt->ylin || (pt->ylin[lo] > 0.0 && pt->ylin[lo + 1] > 0.0));
		yv   = (logv || !pt->ylin) ? pt->y : pt->ylin;
		dv   = (logv || !pt->ylin) ? pt->d : pt->dlin;

		h  = pt->x[lo + 1] - pt->x[lo];
		s  = (t - pt->x[lo]) / h;
		s2 = s * s;
		s3 = s2 * s;
		v  = (2.0 * s3 - 3.0 * s2 + 1.0) * yv[lo] + (s3 - 2.0 * s2 + s) * h * dv[lo] +
			(-2.0 * s3 + 3.0 * s2) * yv[lo + 1] + (s3 - s2) * h * dv[lo + 1];
		y[n] = (logv) ? exp(v) : v;
	}
}
/* ------- end ---------------------------- pchipTableEval.c -------- */

/* ------- begin -------------------------- pchipTableFree.c -------- */

void pchipTableFree(PchipTable *pt)
{
	free(pt->x);
	free(pt->y);
	free(pt->d);
	free(pt->ylin);
	free(pt->dlin);
	pt->x = pt->y = pt->d = pt->ylin = pt->dlin = NULL;
	pt->N                                       = 0;
}
/* ------- end ---------------------------- pchipTableFree.c -------- */
