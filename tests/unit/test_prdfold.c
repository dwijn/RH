/* ------- file: -------------------------- test_prdfold.c ---------

	   Folding of the linear interpolation into the PRD redistribution
	   weights (prdfold.c, P-8). For a line grid q_abs, values J on it,
	   and weights g on the fine absorption grid of prdGrid:

	   -- sum_m W(m) J(m0 + m), with W = prdFold(g) on the range of
		  prdFoldRange, equals sum_q' g(q') Linear(J)(q'), Linear being
		  RH's interpolation (constant beyond the ends of the table), to
		  rounding: for ascending and descending line grids, for the
		  folded fine grid |q'| of symmetric lines, for emission
		  wavelengths in the core, the near wing and far beyond the end
		  of the line grid, and for waveratio != 1 (cross-redistribution);
	   -- the folded weights sum to the sum of g;
	   -- prdFold writes only inside [m0, m0 + len), every point of the
		  fine grid falls inside that range, and the range is no wider
		  than the touched indices;
	   -- prdFineGrid builds the grid of prdGrid's limits, folded to |q|
		  for symmetric lines.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "prdfold.h"
#include "check.h"

#define NTAB  41
#define GUARD 4

/* --- Line grids: dense in the core, wider in the wings, like RH's
	   line wavelength grids (in Doppler units) --     -------------- */

static void make_table(double *q, bool_t half, bool_t descend)
{
	int    n;
	double x;

	for (n = 0; n < NTAB; n++) {
		x                               = (half) ? (double) n / (NTAB - 1) : 2.0 * n / (NTAB - 1) - 1.0;
		x                               = 30.0 * x * fabs(x) + 1.5 * x; /* finer near the core */
		q[(descend) ? NTAB - 1 - n : n] = x;
	}
}

int main(void)
{
	double q_abs[NTAB], J[NTAB], buf[NTAB + 2 * GUARD], *W                   = buf + GUARD, q0, qN,
														*qp, *gw, *Ji, worst = 0.0, worst_sum = 0.0;
	int    it, half, desc, Np, m0, len, m, lap, j, nq, maxlen = 0;
	double fx, waveratio, q_emit;

	for (it = 0; it < NTAB; it++) J[it] = 1.0 + 0.5 * sin(0.7 * it) + 0.01 * it * it;

	for (half = 0; half <= 1; half++) {
		for (desc = 0; desc <= 1; desc++) {
			make_table(q_abs, half, desc);
			for (it = 0; it < 2; it++) {
				waveratio = (it == 0) ? 1.0 : 1.04;
				for (nq = 0; nq <= 160; nq++) {
					q_emit = -45.0 + 90.0 * nq / 160.0;
					Np     = prdGrid(q_emit, waveratio, &q0, &qN);
					qp     = malloc(Np * sizeof(double));
					gw     = malloc(Np * sizeof(double));
					Ji     = malloc(Np * sizeof(double));

					/* --- The fine grid --              -------------- */

					prdFineGrid(q0, Np, FALSE, qp);
					CHECK(fabs(qp[0] - q0) == 0.0 && fabs(qp[Np - 1] - (q0 + (Np - 1) * PRD_DQ)) < 1e-12,
						  "fine grid from %g, %d points", q0, Np);
					for (lap = 0; lap < Np; lap++)
						gw[lap] = exp(-0.1 * (qp[lap] - q_emit) * (qp[lap] - q_emit)) + 0.01;
					if (half) prdFineGrid(q0, Np, TRUE, qp);
					for (lap = 0; lap < Np; lap++)
						CHECK(!half || qp[lap] >= 0.0, "folded grid point %g", qp[lap]);

					/* --- Reference: interpolate, then sum -- ------- */

					Linear(NTAB, q_abs, J, Np, qp, Ji, TRUE);
					double ref = 0.0, gsum = 0.0;
					for (lap = 0; lap < Np; lap++) {
						ref += gw[lap] * Ji[lap];
						gsum += gw[lap];
					}

					/* --- Folded --                     -------------- */

					len = prdFoldRange(NTAB, q_abs, q0, Np, half, &m0);
					CHECK(m0 >= 0 && len >= 1 && m0 + len <= NTAB, "range [%d, %d)", m0, m0 + len);
					if (m0 < 0 || len < 1 || m0 + len > NTAB) {
						free(qp);
						free(gw);
						free(Ji);
						continue;
					}
					maxlen = MAX(maxlen, len);
					for (m = 0; m < NTAB + 2 * GUARD; m++) buf[m] = 0.0;
					prdFold(NTAB, q_abs, Np, qp, gw, m0, len, W);

					double got = 0.0, wsum = 0.0;
					for (m = 0; m < len; m++) {
						got += W[m] * J[m0 + m];
						wsum += W[m];
					}
					for (m = -GUARD; m < 0; m++) CHECK(W[m] == 0.0, "write below the range (%d)", m);
					for (m = len; m < len + GUARD; m++)
						CHECK(W[m] == 0.0, "write above the range (%d)", m);

					/* --- The range is tight: both ends carry weight
						   or are touched by a grid point -- ------- */

					int touch_lo = 0, touch_hi = 0;
					for (lap = 0; lap < Np; lap++) {
						prdFoldIndex(NTAB, q_abs, qp[lap], &j, &fx);
						if (j == m0 || (fx < 1.0 && j + 1 == m0)) touch_lo = 1;
						if (j == m0 + len - 1 || (fx < 1.0 && j + 1 == m0 + len - 1)) touch_hi = 1;
					}
					CHECK(touch_lo && touch_hi, "range [%d, %d) wider than needed", m0, m0 + len);

					worst     = MAX(worst, fabs(got - ref) / fabs(ref));
					worst_sum = MAX(worst_sum, fabs(wsum - gsum) / gsum);
					CHECK(fabs(got - ref) <= 1.0e-13 * fabs(ref),
						  "half %d desc %d waveratio %g q_emit %g: folded %.17g, interpolated %.17g",
						  half, desc, waveratio, q_emit, got, ref);
					CHECK(fabs(wsum - gsum) <= 1.0e-13 * gsum, "sum of folded weights %.17g, of g %.17g",
						  wsum, gsum);
					free(qp);
					free(gw);
					free(Ji);
				}
			}
		}
	}
	printf("folded vs interpolated: largest relative difference %.1e, weight sums %.1e,"
		   " longest row %d of %d\n",
		   worst, worst_sum, maxlen, NTAB);
	return check_summary("test_prdfold");
}
