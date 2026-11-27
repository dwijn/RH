/* ------- file: -------------------------- accelerate.c ------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Mon Apr 17 12:00:40 2000 --

	   --------------------------                      ----------RH-- */

/* --- Implements Ng convergence acceleration with general order.

  See: K. C. Ng 1974, J. Chem. Phys. 61, 2680
 Also: L. Auer 1987, in "Numerical Radiative Transfer",
		 ed. W. Kalkofen, pp. 101-109
	   --                                              -------------- */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "accelerate.h"
#include "error.h"
#include "statistics.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- NgInit.c ---------------- */

struct Ng *NgInit(RHState *rh, int N, int Ndelay, int Norder, int Nperiod,
				  double *solution)
{
	/* --- Initialize data structure and allocate space for previous
		   solutions, coefficient matrix A and correction vector b,
		   copy initial solution into previous solution matrix -- ----- */

	register int k;

	struct Ng *Ngs;

	/* --- The history has Norder + 2 rows and Accelerate takes the
		   count modulo Norder + 2, so a negative order corrupts the
		   allocation or divides by zero (B-45). readInput checks the
		   keywords; this covers callers that do not read them -- -- */
	RH_ASSERT(rh, Norder >= 0 && (Norder == 0 || Nperiod >= 1), "NgInit",
			  "Ng order must be at least 0, and the period at least 1"
			  " when the order is positive (order %d, period %d)",
			  Norder, Nperiod);

	Ngs = (struct Ng *) mallocOrAbort(sizeof(struct Ng));

	Ngs->N       = N;
	Ngs->Norder  = Norder;
	Ngs->Nperiod = Nperiod;
	Ngs->Ndelay  = MAX(Ndelay, Norder + 2);

	if (Norder > 0) {
		Ngs->A = matrix_double(rh, Norder, Norder);
		Ngs->b = (double *) mallocOrAbort(Norder * sizeof(double));
	}
	Ngs->previous   = matrix_double(rh, Norder + 2, N);
	Ngs->theStorage = Ngs->previous[0];
	for (k = 0; k < N; k++)
		Ngs->previous[0][k] = solution[k];
	Ngs->count    = 1;
	Ngs->rejected = 0;

	return Ngs;
}
/* ------- end ---------------------------- NgInit.c ---------------- */

/* ------- begin -------------------------- Accelerate.c ------------ */

rh_status Accelerate(RHState *rh, struct Ng *Ngs, double *solution, bool_t *accelerated)
{
	rh_status    status;
	register int i, j, k;

	int      Norder = Ngs->Norder, ip, ipp, i0;
	double **Delta, *weight, g;

	/* --- Store the current solution --                --------------- */

	i = Ngs->count % (Norder + 2);
	for (k = 0; k < Ngs->N; k++) {
		Ngs->previous[i][k] = solution[k];
	}
	(Ngs->count)++;

	/* --- Accelerate only after we have accumulated enough iterations,
		   and then only every Ngs->period th iteration after iteration
		   Ngs->Ndelay --                               --------------- */

	if ((Norder > 0) && (Ngs->count >= Ngs->Ndelay) &&
		!((Ngs->count - Ngs->Ndelay) % Ngs->Nperiod)) {
		getCPU(rh, 4, TIME_START, NULL);

		Delta  = matrix_double(rh, Norder + 1, Ngs->N);
		weight = (double *) mallocOrAbort(Ngs->N * sizeof(double));

		for (i = 0; i <= Norder; i++) {
			ip  = (Ngs->count - 1 - i) % (Norder + 2);
			ipp = (Ngs->count - 2 - i) % (Norder + 2);
			for (k = 0; k < Ngs->N; k++) {
				Delta[i][k] = Ngs->previous[ip][k] - Ngs->previous[ipp][k];
			}
		}
		/* --- Use weighted acceleration --                -------------- */

		for (k = 0; k < Ngs->N; k++) weight[k] = 1.0 / fabs(solution[k]);
		for (i = 0; i < Norder; i++) {
			Ngs->b[i] = 0.0;
			for (j = 0; j < Norder; j++) Ngs->A[i][j] = 0.0;
		}

		/* --- Fill the coefficient matrix and invert linear set --  ---- */

		for (j = 0; j < Norder; j++) {
			for (k = 0; k < Ngs->N; k++) {
				Ngs->b[j] += weight[k] * Delta[0][k] * (Delta[0][k] - Delta[j + 1][k]);
			}
			for (i = 0; i < Norder; i++) {
				for (k = 0; k < Ngs->N; k++) {
					Ngs->A[i][j] += weight[k] *
									(Delta[j + 1][k] - Delta[0][k]) * (Delta[i + 1][k] - Delta[0][k]);
				}
			}
		}
		if ((status = SolveLinearEq(rh, Norder, Ngs->A, Ngs->b, TRUE)) != RH_OK) {
			free(weight);
			freeMatrix(rh, (void **) Delta);
			*accelerated = FALSE;
			return status;
		}

		/* --- Construct the linear combination for the accelerated
		   solution, and restore the accelerated solution --  ----- */

		i0 = (Ngs->count - 1) % (Norder + 2);
		for (i = 0; i < Norder; i++) {
			ip = (Ngs->count - 2 - i) % (Norder + 2);
			for (k = 0; k < Ngs->N; k++) {
				solution[k] += Ngs->b[i] *
							   (Ngs->previous[ip][k] - Ngs->previous[i0][k]);
			}
		}
		/* --- Keep the step only if it keeps the sign of the
			   unaccelerated solution g (still in previous[i0]): no
			   value below 0 where g >= 0, none at 0 where g > 0, and
			   none non-finite where g is finite. Ng's extrapolation has
			   no safeguard, and far from the solution it can overshoot
			   to negative populations, then negative opacities (N-33).
			   Otherwise restore g: the iteration goes on with the plain
			   step. With floating-point traps on (SetFPEtraps) an
			   overflow traps before this check; without them a
			   non-finite step is rejected here --        -------------- */

		Ngs->rejected = 0;
		for (k = 0; k < Ngs->N; k++) {
			g = Ngs->previous[i0][k];
			if ((isfinite(g) && !isfinite(solution[k])) ||
				(g >= 0.0 && solution[k] < 0.0) || (g > 0.0 && solution[k] == 0.0))
				Ngs->rejected++;
		}
		if (Ngs->rejected) {
			for (k = 0; k < Ngs->N; k++) solution[k] = Ngs->previous[i0][k];
		} else {
			for (k = 0; k < Ngs->N; k++) Ngs->previous[i0][k] = solution[k];
		}
		/* --- Clean up the temporary variables --        --------------- */

		free(weight);
		freeMatrix(rh, (void **) Delta);

		getCPU(rh, 4, TIME_POLL, "Accelerate");
		*accelerated = (Ngs->rejected == 0);
	} else {
		Ngs->rejected = 0;
		*accelerated  = FALSE;
	}
	return RH_OK;
}
/* ------- end ---------------------------- Accelerate.c ------------ */

/* ------- begin -------------------------- NgStatus.c -------------- */

const char *NgStatus(const struct Ng *Ngs, bool_t accel, char *buf, size_t size)
{
	/* --- End of a caller's log line after Accelerate(): whether the
		   step was extrapolated, or the extrapolation rejected (N-33).
		   Returns a string constant, or buf (of size bytes, the
		   caller's) with the number of values rejected -- -------- */

	if (accel)
		return " (accelerated)\n";
	if (Ngs->rejected) {
		snprintf(buf, size,
				 " (Ng step rejected: %ld value%s not positive or not finite;"
				 " plain step kept)\n",
				 Ngs->rejected, (Ngs->rejected == 1) ? "" : "s");
		return buf;
	}
	return "\n";
}
/* ------- end ---------------------------- NgStatus.c -------------- */

/* ------- begin -------------------------- NgFree.c ---------------- */

void NgFree(RHState *rh, struct Ng *Ngs)
{
	freeMatrix(rh, (void **) Ngs->previous);

	if (Ngs->Norder > 0) {
		free(Ngs->b);
		freeMatrix(rh, (void **) Ngs->A);
	}

	free(Ngs);
}
/* ------- end ---------------------------- NgFree.c ---------------- */
