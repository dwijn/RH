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
	   the slopes. A table with a value <= 0 cannot be interpolated in
	   log y; it is interpolated in y, where PCHIP still keeps every
	   interval between its two node values, so non-negative data stay
	   non-negative. Outside the table the end values are used, as by
	   Linear and the spline routines.

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

extern char messageStr[];

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

/* ------- begin -------------------------- pchipTableInit.c -------- */

void pchipTableInit(PchipTable *pt, int N, const double *x, const double *y,
					int scale)
{
	const char   routineName[] = "pchipTableInit";
	register int k;

	bool_t  ascend;
	int     kk;
	double *h, *del, w1, w2;

	if (N < 1) {
		sprintf(messageStr, "Table needs at least one point (N = %d)", N);
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	pt->N    = N;
	pt->logx = (scale & PCHIP_LOGX) != 0;
	pt->logy = (scale & PCHIP_LOGY) != 0;
	if (pt->logy) {
		for (k = 0; k < N; k++)
			if (!(y[k] > 0.0)) pt->logy = FALSE;
	}
	pt->x = (double *) malloc(N * sizeof(double));
	pt->y = (double *) malloc(N * sizeof(double));
	pt->d = (double *) calloc(N, sizeof(double));

	/* --- Store ascending in x, transformed --         -------------- */

	ascend = (N == 1 || x[1] > x[0]);
	for (k = 0; k < N; k++) {
		kk = (ascend) ? k : N - 1 - k;
		if (pt->logx && !(x[kk] > 0.0)) {
			sprintf(messageStr, "Interpolation in log x needs x > 0 (x = %g)",
					x[kk]);
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
		pt->x[k] = (pt->logx) ? log(x[kk]) : x[kk];
		pt->y[k] = (pt->logy) ? log(y[kk]) : y[kk];
	}
	for (k = 1; k < N; k++) {
		if (!(pt->x[k] > pt->x[k - 1])) {
			sprintf(messageStr, "Table abscissae are not strictly monotonic"
								" (point %d of %d)",
					k, N);
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
	}
	if (N == 1) return;

	h   = (double *) malloc((N - 1) * sizeof(double));
	del = (double *) malloc((N - 1) * sizeof(double));
	for (k = 0; k < N - 1; k++) {
		h[k]   = pt->x[k + 1] - pt->x[k];
		del[k] = (pt->y[k + 1] - pt->y[k]) / h[k];
	}
	if (N == 2) {
		pt->d[0] = pt->d[1] = del[0];
	} else {
		for (k = 1; k < N - 1; k++) {
			if (del[k - 1] * del[k] > 0.0) {
				w1       = 2.0 * h[k] + h[k - 1];
				w2       = h[k] + 2.0 * h[k - 1];
				pt->d[k] = (w1 + w2) / (w1 / del[k - 1] + w2 / del[k]);
			} else
				pt->d[k] = 0.0;
		}
		pt->d[0]     = endSlope(h[0], h[1], del[0], del[1]);
		pt->d[N - 1] = endSlope(h[N - 2], h[N - 3], del[N - 2], del[N - 3]);
	}
	free(h);
	free(del);
}
/* ------- end ---------------------------- pchipTableInit.c -------- */

/* ------- begin -------------------------- pchipTableEval.c -------- */

void pchipTableEval(const PchipTable *pt, int N, const double *x, double *y)
{
	register int n;

	int    lo, hi, mid, N1 = pt->N - 1;
	double t, h, s, s2, s3, v;

	for (n = 0; n < N; n++) {
		t = (pt->logx) ? ((x[n] > 0.0) ? log(x[n]) : -HUGE_VAL) : x[n];

		if (t <= pt->x[0])
			v = pt->y[0];
		else if (t >= pt->x[N1])
			v = pt->y[N1];
		else {
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
			h  = pt->x[lo + 1] - pt->x[lo];
			s  = (t - pt->x[lo]) / h;
			s2 = s * s;
			s3 = s2 * s;
			v  = (2.0 * s3 - 3.0 * s2 + 1.0) * pt->y[lo] +
				(s3 - 2.0 * s2 + s) * h * pt->d[lo] +
				(-2.0 * s3 + 3.0 * s2) * pt->y[lo + 1] + (s3 - s2) * h * pt->d[lo + 1];
		}
		y[n] = (pt->logy) ? exp(v) : v;
	}
}
/* ------- end ---------------------------- pchipTableEval.c -------- */

/* ------- begin -------------------------- pchipTableFree.c -------- */

void pchipTableFree(PchipTable *pt)
{
	free(pt->x);
	free(pt->y);
	free(pt->d);
	pt->x = pt->y = pt->d = NULL;
	pt->N                 = 0;
}
/* ------- end ---------------------------- pchipTableFree.c -------- */
