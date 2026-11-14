/* ------- file: -------------------------- solvelong.c -------------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu May 31 11:17:51 2018 --

	   --------------------------                      ----------RH-- */

/* --- Piecewise integration of source function for
	   long characteristics. --                        -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

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
	/* --- The polarized counterpart of SolveLong: the same points, and
		   each section integrated with DELO-linear (LongChar_Stokes_Linear)
		   or DELO-Bezier (LongChar_Stokes_Bezier3), as
		   S_INTERPOLATION_LONGCHAR sets. The parabolic integration used
		   before had the downwind weight of N-39, here of the local
		   point (N-41). K' at the local point is needed only by the
		   Bezier step. I_uw is returned at the upwind point -- ------ */

	register int ls, n;

	bool_t bezier = (input.S_interpolation_longchar == S_BEZIER3);
	double chi_uw, chi_loc, chi_dw, ds_uw, ds_dw, S_uw[4], S_loc[4], S_dw[4],
		I_loc[4], K_uw[4][4], K_loc[4][4], K_dw[4][4];

	chi_uw = Interp_2D(chi, &lc->stencil[0], UPWIND);
	StokesK_2D(nspect, &lc->stencil[0], UPWIND, chi_uw, K_uw);
	chi_loc = Interp_2D(chi, &lc->stencil[1], UPWIND);
	StokesK_2D(nspect, &lc->stencil[1], UPWIND, chi_loc, K_loc);
	for (n = 0; n < 4; n++) {
		S_uw[n]  = Interp_2D(S[n], &lc->stencil[0], UPWIND);
		I_uw[n]  = Interp_2D(I[n], &lc->stencil[0], UPWIND);
		S_loc[n] = Interp_2D(S[n], &lc->stencil[1], UPWIND);
	}
	ds_uw = lc->stencil[0].ds[UPWIND];

	for (ls = 2; ls <= lc->Nlc; ls++) {
		if (ls == lc->Nlc) {
			chi_dw = chi[local];
			if (bezier) StokesK(nspect, local, chi_dw, K_dw);
			for (n = 0; n < 4; n++) S_dw[n] = S[n][local];
		} else {
			chi_dw = Interp_2D(chi, &lc->stencil[ls], UPWIND);
			StokesK_2D(nspect, &lc->stencil[ls], UPWIND, chi_dw, K_dw);
			for (n = 0; n < 4; n++)
				S_dw[n] = Interp_2D(S[n], &lc->stencil[ls], UPWIND);
		}
		ds_dw = lc->stencil[ls - 1].ds[UPWIND];

		if (bezier)
			LongChar_Stokes_Bezier3(ds_uw, ds_dw, chi_uw, chi_loc, chi_dw,
									K_uw, K_loc, K_dw, S_uw, S_loc, S_dw,
									I_uw, I_loc);
		else
			LongChar_Stokes_Linear(ds_uw, chi_uw, chi_loc, K_uw, K_loc,
								   S_uw, S_loc, I_uw, I_loc);

		memcpy(I_uw, I_loc, sizeof(I_loc));
		if (ls == lc->Nlc) break;

		memcpy(S_uw, S_loc, sizeof(S_loc));
		memcpy(S_loc, S_dw, sizeof(S_dw));
		memcpy(K_uw, K_loc, sizeof(K_loc));
		memcpy(K_loc, K_dw, sizeof(K_dw));
		chi_uw  = chi_loc;
		chi_loc = chi_dw;
		ds_uw   = ds_dw;
	}
}
/* ------- end ---------------------------- SolveLongStokes.c ------- */
