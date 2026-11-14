/* ------- file: -------------------------- solvelong.c -------------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu May 31 11:17:51 2018 --

	   --------------------------                      ----------RH-- */

/* --- Piecewise integration of source function for
	   long characteristics. --                        -------------- */

#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "inputs.h"
#include "bezier.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern InputData input;

/* ------- begin -------------------------- SolveLong.c ------------- */

/* --- Solves transfer equation along a long characteristic specified in
	   LongChar structure pointed to by lc. Used when, for periodic
	   boundary conditions, no short characteristic crosses a horizontal
	   grid line in the upwind direction.

	   The characteristic starts on the horizontal grid line
	   (stencil[0]), crosses vertical grid lines at stencil[1 .. Nlc - 1],
	   the last of which is the upwind point of the local short
	   characteristic, and ends at the local point. Section n runs from
	   point n - 1 to point n over stencil[n - 1].ds[UPWIND]; the
	   intensity is returned at the upwind point (B-24).

	   Each section is integrated with S linear in tau or with the cubic
	   Bezier step of the short characteristics, as
	   S_INTERPOLATION_LONGCHAR sets. Both keep I between the extremes of
	   the upwind intensity and the source function; the parabolic
	   integration used before did not (N-39, see rhsc3d/solvelong.c)
	   --                                              -------------- */

double SolveLong(LongChar *lc, int local, double *chi, double *S, double *I)
{
	register int ls;

	double chi_uw, chi_loc, chi_dw, S_uw, S_loc, S_dw, I_uw, ds_uw, ds_dw;

	chi_uw = Interp_2D(chi, &lc->stencil[0], UPWIND);
	S_uw   = Interp_2D(S, &lc->stencil[0], UPWIND);
	I_uw   = Interp_2D(I, &lc->stencil[0], UPWIND);

	chi_loc = Interp_2D(chi, &lc->stencil[1], UPWIND);
	S_loc   = Interp_2D(S, &lc->stencil[1], UPWIND);
	ds_uw   = lc->stencil[0].ds[UPWIND];

	for (ls = 2; ls <= lc->Nlc; ls++) {
		if (ls == lc->Nlc) {
			chi_dw = chi[local];
			S_dw   = S[local];
		} else {
			chi_dw = Interp_2D(chi, &lc->stencil[ls], UPWIND);
			S_dw   = Interp_2D(S, &lc->stencil[ls], UPWIND);
		}
		ds_dw = lc->stencil[ls - 1].ds[UPWIND];

		if (input.S_interpolation_longchar == S_BEZIER3)
			I_uw = LongChar_Bezier3(ds_uw, ds_dw, chi_uw, chi_loc, chi_dw,
									S_uw, S_loc, S_dw, I_uw);
		else
			I_uw = LongChar_Linear(ds_uw, chi_uw, chi_loc, S_uw, S_loc, I_uw);

		chi_uw  = chi_loc;
		S_uw    = S_loc;
		chi_loc = chi_dw;
		S_loc   = S_dw;
		ds_uw   = ds_dw;
	}

	return I_uw;
}
/* ------- end ---------------------------- SolveLong.c ------------- */

/* ------- begin -------------------------- SolveLongStokes.c ------- */

void SolveLongStokes(int nspect, LongChar *lc, int local, double *chi,
					 double **S, double **I, double *I_uw)
{
	register int ls, n, m;

	double c1, c2, dtau_dw, dS_uw[4], dS_dw[4], w[3],
		chi_loc, S_loc[4], dtau_uw, chi_uw, S_uw[4], chi_dw, S_dw[4],
		P[4], Q[4][4], **R, K[4][4], K_uw[4][4];

	R = matrix_double(4, 4);

	/* --- The first point is the intersection with the horizontal
		   grid line --                                  -------------- */

	chi_uw = Interp_2D(chi, &lc->stencil[0], UPWIND);
	StokesK_2D(nspect, &lc->stencil[0], UPWIND, chi_uw, K_uw);

	for (n = 0; n < 4; n++) {
		S_uw[n]  = Interp_2D(S[n], &lc->stencil[0], UPWIND);
		I_uw[n]  = Interp_2D(I[n], &lc->stencil[0], UPWIND);
		S_loc[n] = Interp_2D(S[n], &lc->stencil[1], UPWIND);
	}
	chi_loc = Interp_2D(chi, &lc->stencil[1], UPWIND);
	StokesK_2D(nspect, &lc->stencil[1], UPWIND, chi_loc, K);

	dtau_uw = 0.5 * (chi_uw + chi_loc) * lc->stencil[0].ds[UPWIND];

	/* --- As in SolveLong: the local point is the downwind end of the
		   last section, I_uw is returned at the upwind point -- ------ */

	for (ls = 2; ls <= lc->Nlc; ls++) {
		if (ls == lc->Nlc) {
			chi_dw = chi[local];
			for (n = 0; n < 4; n++) S_dw[n] = S[n][local];
		} else {
			chi_dw = Interp_2D(chi, &lc->stencil[ls], UPWIND);
			for (n = 0; n < 4; n++)
				S_dw[n] = Interp_2D(S[n], &lc->stencil[ls], UPWIND);
		}
		dtau_dw = 0.5 * (chi_loc + chi_dw) * lc->stencil[ls - 1].ds[UPWIND];
		w3(dtau_uw, w);

		for (n = 0; n < 4; n++) {
			dS_uw[n] = (S_uw[n] - S_loc[n]) / dtau_uw;
			dS_dw[n] = (S_loc[n] - S_dw[n]) / dtau_dw;
			c1       = dS_uw[n] * dtau_dw + dS_dw[n] * dtau_uw;
			c2       = dS_uw[n] - dS_dw[n];
			P[n]     = w[0] * S_loc[n] + (w[1] * c1 + w[2] * c2) / (dtau_uw + dtau_dw);
		}
		for (n = 0; n < 4; n++) {
			for (m = 0; m < 4; m++) {
				Q[n][m] = -w[1] / dtau_uw * K_uw[n][m];
				R[n][m] = (w[0] - w[1] / dtau_uw) * K[n][m];
			}
			Q[n][n] = 1.0 - w[0];
			R[n][n] = 1.0;
		}
		for (n = 0; n < 4; n++) {
			for (m = 0; m < 4; m++)
				P[n] += Q[n][m] * I_uw[m];
		}
		/* --- Solve linear equations for I --             -------------- */

		SolveLinearEq(4, R, P, TRUE);

		/* --- Store results for the upwind Stokes vector -- ------------ */

		for (n = 0; n < 4; n++) I_uw[n] = P[n];

		/* --- Reuse upwind quantities --                  -------------- */

		if (ls < lc->Nlc) {
			chi_uw  = chi_loc;
			chi_loc = chi_dw;
			dtau_uw = dtau_dw;

			for (n = 0; n < 4; n++) {
				S_uw[n]  = S_loc[n];
				S_loc[n] = S_dw[n];
				for (m = 0; m < 4; m++) K_uw[n][m] = K[n][m];
			}
			StokesK_2D(nspect, &lc->stencil[ls],
					   UPWIND, chi_dw, K);
		}
	}
	freeMatrix((void **) R);
}
/* ------- end ---------------------------- SolveLongStokes.c ------- */
