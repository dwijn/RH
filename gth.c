/* ------- file: -------------------------- gth.c -------------------

	   Rate equations solved by GTH elimination (N-31).

	   --------------------------                      ----------RH-- */

/* --- The rate equations of statistical equilibrium,

		 sum_{j != i} n_j R_ji - n_i sum_{j != i} R_ij = 0,
		 sum_i n_i = ntotal,

	   with R_ij >= 0 the rate from level i to level j, describe the
	   stationary state of a continuous-time Markov chain. Gaussian
	   elimination of the equations as a linear system (SolveLinearEq) is
	   normwise stable but not componentwise accurate: each diagonal
	   element -sum_j R_ij is formed before the elimination, and where
	   levels are tightly coupled (fine-structure levels bound by large
	   collisional rates) its rounding erases the small rates that leave
	   the group, so the small populations and the populations of such
	   groups can be wrong in their leading digits. No scaling helps
	   (He_large: componentwise condition number 1e14 with any row or
	   column scaling).

	   The GTH variant of the elimination (Grassmann, Taksar & Heyman
	   1985, Oper. Res. 33, 1107) works on the off-diagonal rates only.
	   Eliminating level l folds its rates into those of the remaining
	   levels,

		 S_l   = sum_{j < l} R_lj,
		 R_il <- R_il / S_l,
		 R_ij <- R_ij + R_il R_lj           (i, j < l, i != j),

	   and the populations follow by back substitution,

		 p_root = 1,  p_l = sum_{i < l} p_i R_il.

	   Every pivot S_l is a sum of non-negative rates, and the update
	   only adds non-negative terms, so there is no subtraction and no
	   cancellation: the populations are accurate to a few rounding
	   errors each, whatever the conditioning (O'Cinneide 1993, Numer.
	   Math. 65, 109). The cost is that of Gaussian elimination.

	   It needs non-negative rates. A negative rate can come from the
	   preconditioning of the radiative rates or from bad atomic data,
	   and a zero pivot from a level that has no rates to the levels
	   left; SolveGTH then returns without changing n and the caller
	   falls back to SolveLinearEq.

	   The level eliminated last (root, p = 1) is the one the others are
	   measured against. The caller passes the most populated level, so
	   that p does not overflow when the populations span a very wide
	   range (several ionization stages); p is rescaled in the back
	   substitution if it grows beyond GTH_BIG regardless.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>

#include "rh.h"

#define GTH_BIG 1.0E+150

/* ------- begin -------------------------- SolveGTH.c -------------- */

enum GTHStatus SolveGTH(int N, double **R, int root, double ntotal, double *n,
						int *i_bad, int *j_bad)
{
	register int a, b, l;

	int           *lev, i, j;
	double         S, f, *p, total;
	enum GTHStatus status = GTH_OK;

	*i_bad = *j_bad = -1;

	/* --- Check the rates: off-diagonal, finite and non-negative -- -- */

	for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (i == j) continue;
			if (!isfinite(R[i][j])) {
				*i_bad = i;
				*j_bad = j;
				return GTH_NONFINITE;
			}
			if (R[i][j] < 0.0) {
				*i_bad = i;
				*j_bad = j;
				return GTH_NEGATIVE;
			}
		}
	}
	if (N == 1) {
		n[0] = ntotal;
		return GTH_OK;
	}
	/* --- Elimination order: the root first (eliminated last), the
		   other levels in their own order --            -------------- */

	lev    = (int *) malloc(N * sizeof(int));
	p      = (double *) malloc(N * sizeof(double));
	lev[0] = root;
	for (i = 0, a = 1; i < N; i++)
		if (i != root) lev[a++] = i;

	for (l = N - 1; l > 0; l--) {
		S = 0.0;
		for (b = 0; b < l; b++) S += R[lev[l]][lev[b]];
		if (!(S > 0.0)) {
			*i_bad = lev[l];
			status = GTH_REDUCIBLE;
			break;
		}
		for (a = 0; a < l; a++) R[lev[a]][lev[l]] /= S;
		for (a = 0; a < l; a++) {
			f = R[lev[a]][lev[l]];
			if (f == 0.0) continue;
			for (b = 0; b < l; b++)
				if (b != a) R[lev[a]][lev[b]] += f * R[lev[l]][lev[b]];
		}
	}

	if (status == GTH_OK) {
		p[0] = 1.0;
		for (l = 1; l < N; l++) {
			p[l] = 0.0;
			for (a = 0; a < l; a++) p[l] += p[a] * R[lev[a]][lev[l]];
			if (p[l] > GTH_BIG)
				for (a = 0; a <= l; a++) p[a] /= GTH_BIG;
		}
		total = 0.0;
		for (l = 0; l < N; l++) total += p[l];
		if (!(isfinite(total) && total > 0.0))
			status = GTH_NONFINITE;
		else
			for (l = 0; l < N; l++) n[lev[l]] = p[l] * (ntotal / total);
	}
	free(lev);
	free(p);

	return status;
}
/* ------- end ---------------------------- SolveGTH.c -------------- */
