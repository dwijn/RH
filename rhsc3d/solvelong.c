/* ------- file: -------------------------- solvelong.c -------------

	   Version:       rh2.0, 3-D Cartesian, short characteristics
	   Author:        Han Uitenbroek(huitenbroek@nso.edu)
	   Last modified: Fri Aug 29 15:32:40 2008 --

	   --------------------------                      ----------RH-- */

/* --- Piecewise integration of source function for
	   long characteristic --                          -------------- */

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
	/* --- The polarized counterpart of SolveLong: the same points, and
		   each section integrated with DELO-linear (LongChar_Stokes_Linear)
		   or DELO-Bezier (LongChar_Stokes_Bezier3), as
		   S_INTERPOLATION_LONGCHAR sets. The parabolic integration used
		   before had the downwind weight of N-39, here of the local
		   point (N-41). K' at the local point is needed only by the
		   Bezier step. I_uw is returned at the upwind point -- ------ */

	register int ls, n;

	int    local;
	bool_t bezier = (input.S_interpolation_longchar == S_BEZIER3);
	double chi_uw, chi_loc, chi_dw, ds_uw, ds_dw, S_uw[4], S_loc[4], S_dw[4],
		I_loc[4], K_uw[4][4], K_loc[4][4], K_dw[4][4];

	chi_uw = Interpolate_3D(chi, geometry, &lc->stencil[0], l, m);
	StokesK_3D(nspect, geometry, &lc->stencil[0], l, m, chi_uw, K_uw);
	chi_loc = Interpolate_3D(chi, geometry, &lc->stencil[1], l, m);
	StokesK_3D(nspect, geometry, &lc->stencil[1], l, m, chi_loc, K_loc);
	for (n = 0; n < 4; n++) {
		S_uw[n]  = Interpolate_3D(S[n], geometry, &lc->stencil[0], l, m);
		I_uw[n]  = Interpolate_3D(I[n], geometry, &lc->stencil[0], l, m);
		S_loc[n] = Interpolate_3D(S[n], geometry, &lc->stencil[1], l, m);
	}
	ds_uw = lc->stencil[0].ds;

	for (ls = 2; ls <= lc->Nst; ls++) {
		if (ls == lc->Nst) {
			local  = k * geometry->Nplane + m * geometry->Nx + l;
			chi_dw = chi[local];
			if (bezier) StokesK(nspect, local, chi_dw, K_dw);
			for (n = 0; n < 4; n++) S_dw[n] = S[n][local];
		} else {
			chi_dw = Interpolate_3D(chi, geometry, &lc->stencil[ls], l, m);
			StokesK_3D(nspect, geometry, &lc->stencil[ls], l, m, chi_dw, K_dw);
			for (n = 0; n < 4; n++)
				S_dw[n] = Interpolate_3D(S[n], geometry, &lc->stencil[ls], l, m);
		}
		ds_dw = lc->stencil[ls - 1].ds;

		if (bezier)
			LongChar_Stokes_Bezier3(ds_uw, ds_dw, chi_uw, chi_loc, chi_dw,
									K_uw, K_loc, K_dw, S_uw, S_loc, S_dw,
									I_uw, I_loc);
		else
			LongChar_Stokes_Linear(ds_uw, chi_uw, chi_loc, K_uw, K_loc,
								   S_uw, S_loc, I_uw, I_loc);

		memcpy(I_uw, I_loc, sizeof(I_loc));
		if (ls == lc->Nst) break;

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
