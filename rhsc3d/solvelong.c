/* ------- file: -------------------------- solvelong.c -------------

	   Version:       rh2.0, 3-D Cartesian, short characteristics
	   Author:        Han Uitenbroek(huitenbroek@nso.edu)
	   Last modified: Fri Aug 29 15:32:40 2008 --

	   --------------------------                      ----------RH-- */

/* --- Piecewise integration of source function for
	   long characteristic --                          -------------- */

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
	   grid plane in the upwind direction.

	   The characteristic starts on the horizontal plane (stencil[0]),
	   crosses vertical planes at stencil[1 .. Nst - 1], the last of
	   which is the upwind point of the local short characteristic, and
	   ends at the local point. Section n runs from point n - 1 to point
	   n over stencil[n - 1].ds; the intensity is returned at the upwind
	   point, after Nst - 1 sections (B-24).

	   Each section is integrated with S linear in tau or with the cubic
	   Bezier step of the short characteristics, as S_INTERPOLATION_LONGCHAR
	   sets. Both keep I between the extremes of the upwind intensity and
	   the source function. The parabolic integration used before (N-39)
	   did not: the weight of the downwind S is negative for an optically
	   thick upwind section and about -1/dtau_dw for a thin downwind one,
	   and since B-24 the downwind point of the last section is the local
	   point, so that J depended negatively on the local S and the
	   Lambda iteration of background scattering diverged -- ------- */

double SolveLong(Geometry *geometry, Longchar *lc,
				 int k, int l, int m, double *chi, double *S, double *I)
{
	register int ls;

	int    local;
	double chi_uw, chi_loc, chi_dw, S_uw, S_loc, S_dw, I_uw, ds_uw, ds_dw;

	chi_uw = Interpolate_3D(chi, geometry, &lc->stencil[0], l, m);
	S_uw   = Interpolate_3D(S, geometry, &lc->stencil[0], l, m);
	I_uw   = Interpolate_3D(I, geometry, &lc->stencil[0], l, m);

	chi_loc = Interpolate_3D(chi, geometry, &lc->stencil[1], l, m);
	S_loc   = Interpolate_3D(S, geometry, &lc->stencil[1], l, m);
	ds_uw   = lc->stencil[0].ds;

	for (ls = 2; ls <= lc->Nst; ls++) {
		if (ls == lc->Nst) {
			local  = k * geometry->Nplane + m * geometry->Nx + l;
			chi_dw = chi[local];
			S_dw   = S[local];
		} else {
			chi_dw = Interpolate_3D(chi, geometry, &lc->stencil[ls], l, m);
			S_dw   = Interpolate_3D(S, geometry, &lc->stencil[ls], l, m);
		}
		ds_dw = lc->stencil[ls - 1].ds;

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
	/* --- Return the upwind intensity at the nearest plane for the
		   pertinent short characteristic --             -------------- */

	return I_uw;
}
/* ------- end ---------------------------- SolveLong.c ------------- */

/* ------- begin -------------------------- SolveLongStokes.c ------- */

void SolveLongStokes(Geometry *geometry, Longchar *lc,
					 int nspect, int k, int l, int m,
					 double *chi, double **S, double **I,
					 double *I_uw)
{
	register int ls, n, j;

	int    local;
	double c1, c2, dtau_dw, dS_uw[4], dS_dw[4], w[3],
		chi_loc, S_loc[4], dtau_uw, chi_uw, S_uw[4], chi_dw, S_dw[4],
		P[4], Q[4][4], **R, K[4][4], K_uw[4][4];

	R = matrix_double(4, 4);

	/* --- The first point is the intersection with the horizontal
		   grid line --                                  -------------- */

	chi_uw = Interpolate_3D(chi, geometry, &lc->stencil[0], l, m);
	StokesK_3D(nspect, geometry, &lc->stencil[0], l, m, chi_uw, K_uw);
	for (n = 0; n < 4; n++) {
		S_uw[n] = Interpolate_3D(S[n], geometry, &lc->stencil[0], l, m);
		I_uw[n] = Interpolate_3D(I[n], geometry, &lc->stencil[0], l, m);
	}

	chi_loc = Interpolate_3D(chi, geometry, &lc->stencil[1], l, m);
	dtau_uw = 0.5 * (chi_uw + chi_loc) * lc->stencil[0].ds;
	StokesK_3D(nspect, geometry, &lc->stencil[1], l, m, chi_loc, K);
	for (n = 0; n < 4; n++)
		S_loc[n] = Interpolate_3D(S[n], geometry, &lc->stencil[1], l, m);

	/* --- As in SolveLong: the local point is the downwind end of the
		   last section, I_uw is returned at the upwind point -- ------ */

	for (ls = 2; ls <= lc->Nst; ls++) {
		if (ls == lc->Nst) {
			local  = k * geometry->Nplane + m * geometry->Nx + l;
			chi_dw = chi[local];
			for (n = 0; n < 4; n++) S_dw[n] = S[n][local];
		} else {
			chi_dw = Interpolate_3D(chi, geometry, &lc->stencil[ls], l, m);
			for (n = 0; n < 4; n++)
				S_dw[n] = Interpolate_3D(S[n], geometry, &lc->stencil[ls], l, m);
		}
		dtau_dw = 0.5 * (chi_loc + chi_dw) * lc->stencil[ls - 1].ds;
		w3(dtau_uw, w);

		for (n = 0; n < 4; n++) {
			dS_uw[n] = (S_uw[n] - S_loc[n]) / dtau_uw;
			dS_dw[n] = (S_loc[n] - S_dw[n]) / dtau_dw;
			c1       = dS_uw[n] * dtau_dw + dS_dw[n] * dtau_uw;
			c2       = dS_uw[n] - dS_dw[n];
			P[n]     = w[0] * S_loc[n] + (w[1] * c1 + w[2] * c2) / (dtau_uw + dtau_dw);
		}
		for (n = 0; n < 4; n++) {
			for (j = 0; j < 4; j++) {
				Q[n][j] = -w[1] / dtau_uw * K_uw[n][j];
				R[n][j] = (w[0] - w[1] / dtau_uw) * K[n][j];
			}
			Q[n][n] = 1.0 - w[0];
			R[n][n] = 1.0;
		}
		for (n = 0; n < 4; n++) {
			for (j = 0; j < 4; j++)
				P[n] += Q[n][j] * I_uw[j];
		}
		/* --- Solve linear equations for I --             -------------- */

		SolveLinearEq(4, R, P, TRUE);

		/* --- Store results for the upwind Stokes vector -- ------------ */

		for (n = 0; n < 4; n++) I_uw[n] = P[n];

		/* --- Reuse upwind quantities --                  -------------- */

		if (ls < lc->Nst) {
			chi_uw  = chi_loc;
			chi_loc = chi_dw;
			dtau_uw = dtau_dw;

			for (n = 0; n < 4; n++) {
				S_uw[n]  = S_loc[n];
				S_loc[n] = S_dw[n];
				for (j = 0; j < 4; j++) K_uw[n][j] = K[n][j];
			}
			StokesK_3D(nspect, geometry, &lc->stencil[ls], l, m, chi_dw, K);
		}
	}
	freeMatrix((void **) R);
}
/* ------- end ---------------------------- SolveLongStokes.c ------- */
