/* ------- file: -------------------------- prdfold.c ---------------

	   The fine absorption grid of the PRD scattering integrals and the
	   folding of linear interpolation into the redistribution weights
	   (see the comments in scatter.c, P-8). Separate from scatter.c so
	   that the unit tests can check them (tests/unit/test_prdfold.c).
	   --                                              -------------- */

#include <math.h>

#include "rh.h"
#include "atom.h"
#include "prdfold.h"

/* ------- begin -------------------------- prdGrid.c --------------- */

int prdGrid(double q_emit, double waveratio, double *q0, double *qN)
{
	/* --- Establish integration limits over absorption wavelength,
		   using only regions where the redistribution function is
		   non-zero, and the number of points of the fine grid.
		   (See also function GII.) --                   -------------- */

	if (fabs(q_emit) < PRD_QCORE) {
		*q0 = -PRD_QWING;
		*qN = PRD_QWING;
	} else {
		if (fabs(q_emit) < PRD_QWING) {
			if (q_emit > 0.0) {
				*q0 = -PRD_QWING;
				*qN = waveratio * (q_emit + PRD_QSPREAD);
			} else {
				*q0 = waveratio * (q_emit - PRD_QSPREAD);
				*qN = PRD_QWING;
			}
		} else {
			*q0 = waveratio * (q_emit - PRD_QSPREAD);
			*qN = waveratio * (q_emit + PRD_QSPREAD);
		}
	}
	return (int) ((*qN - *q0) / PRD_DQ) + 1;
}
/* ------- end ---------------------------- prdGrid.c --------------- */

/* ------- begin -------------------------- prdFineGrid.c ----------- */

/* --- Point lap of the fine grid starting at q0, folded to |q| for
	   symmetric lines, whose wavelength grid covers one half of the
	   line. prdFineGrid and prdFoldRange both use it, so that the range
	   of a row is exactly that of its grid points --   -------------- */

static double prdFinePoint(double q0, int lap, bool_t symmetric)
{
	double q = q0 + lap * PRD_DQ;

	return (symmetric) ? fabs(q) : q;
}

void prdFineGrid(double q0, int Np, bool_t symmetric, double *qp)
{
	int lap;

	for (lap = 0; lap < Np; lap++) qp[lap] = prdFinePoint(q0, lap, symmetric);
}
/* ------- end ---------------------------- prdFineGrid.c ----------- */

/* ------- begin -------------------------- prdFoldIndex.c ---------- */

void prdFoldIndex(int N, double *qtab, double x, int *j, double *fx)
{
	/* --- Linear interpolation in the table qtab as in Linear: the
		   value at x is fx * y[j] + (1 - fx) * y[j + 1], and constant
		   (fx = 1) beyond the ends --                     -------------- */

	bool_t ascend = (qtab[1] > qtab[0]);
	double xmin = (ascend) ? qtab[0] : qtab[N - 1], xmax = (ascend) ? qtab[N - 1] : qtab[0];

	if (x <= xmin) {
		*j  = (ascend) ? 0 : N - 1;
		*fx = 1.0;
	} else if (x >= xmax) {
		*j  = (ascend) ? N - 1 : 0;
		*fx = 1.0;
	} else {
		Locate(N, qtab, x, j);
		*fx = (qtab[*j + 1] - x) / (qtab[*j + 1] - qtab[*j]);
	}
}
/* ------- end ---------------------------- prdFoldIndex.c ---------- */

/* ------- begin -------------------------- prdFoldRange.c ---------- */

int prdFoldRange(int N, double *qtab, double q0, int Np, bool_t symmetric, int *m0)
{
	/* --- Range [m0, m0 + len) of the line's grid that the linear
		   interpolation onto the fine grid touches; returns len. The
		   touched indices move monotonically with x, so the smallest
		   and the largest x of the fine grid bound the range. Both
		   follow from the ends of the grid, and for a symmetric line
		   that crosses q = 0 from the points next to 0, without
		   building the grid: this runs for every row on every call -- */

	int    j, lo, hi, i0, i;
	double xa = prdFinePoint(q0, 0, symmetric), xb = prdFinePoint(q0, Np - 1, symmetric),
		   xlo = MIN(xa, xb), xhi = MAX(xa, xb), x, fx;

	if (symmetric && q0 < 0.0 && q0 + (Np - 1) * PRD_DQ > 0.0) {
		i0 = (int) floor(-q0 / PRD_DQ);
		for (i = i0 - 1; i <= i0 + 2; i++) {
			if (i < 0 || i >= Np) continue;
			x = prdFinePoint(q0, i, symmetric);
			if (x < xlo) xlo = x;
		}
	}
	prdFoldIndex(N, qtab, xlo, &j, &fx);
	lo = j;
	hi = (fx < 1.0) ? j + 1 : j;
	prdFoldIndex(N, qtab, xhi, &j, &fx);
	lo  = MIN(lo, j);
	hi  = MAX(hi, (fx < 1.0) ? j + 1 : j);
	*m0 = lo;
	return hi - lo + 1;
}
/* ------- end ---------------------------- prdFoldRange.c ---------- */

/* ------- begin -------------------------- prdFold.c --------------- */

void prdFold(int N, double *qtab, int Np, double *qp, double *gw,
			 int m0, int len, double *W)
{
	/* --- W[m - m0] += sum over the fine grid of gw times the
		   interpolation coefficient of line point m --   -------------- */

	int    lap, j;
	double fx;

	for (lap = 0; lap < Np; lap++) {
		prdFoldIndex(N, qtab, qp[lap], &j, &fx);
		W[j - m0] += gw[lap] * fx;
		if (fx < 1.0) W[j + 1 - m0] += gw[lap] * (1.0 - fx);
	}
	(void) len;
}
/* ------- end ---------------------------- prdFold.c --------------- */
