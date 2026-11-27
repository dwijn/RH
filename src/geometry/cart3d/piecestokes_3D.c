/* ------- file: -------------------------- piecestokes_3D.c --------

	   Version:       rh2.0, 3-D Cartesian, short characteristics
	   Author:        Han Uitenbroek(huitenbroek@nso.edu)
	   Last modified: Tue Sep 23 12:29:54 2025 --

	   --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "inputs.h"
#include "error.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- ShortChar_Stokes.c ------ */

void ShortChar_Stokes(RHState *rh, Geometry *geometry, int nspect, int mu,
					  bool_t  to_observer,
					  double *chi, double **S, double **I, double *Psi)
{
	const char   routineName[] = "ShortChar_Stokes";
	register int k, l, m, n;

	int Nz = geometry->Nz, Nplane = geometry->Nplane, quadrant,
		dl, lstart, lend, dm, mstart, mend, dk, kstart, kend, local;
	double     I_uw[4], T_uw, Bnu_uw, Bnu, chi_uw, dtau_uw;
	Stencil   *st_uw, *st_dw;
	StokesGrid sg;

	if (to_observer) {

		/* --- Boundary condition at the BOTTOM of the grid -- ---------- */

		switch (geometry->z_boundary_value[BOTTOM]) {
		case IRRADIATED:
			for (l = 0; l < Nplane; l++)
				I[0][(Nz - 1) * Nplane + l] = geometry->Ibottom[nspect][l];
			break;

		case ZERO:
			for (l = (Nz - 1) * Nplane; l < Nz * Nplane; l++)
				I[0][l] = 0.0;
			break;

		case THERMALIZED:
			if (to_observer) {
				local = (geometry->Nz - 1) * geometry->Nplane;
				st_uw = &geometry->stencil_uw[mu][geometry->Nz - 1];

				for (m = 0; m < geometry->Ny; m++) {
					for (l = 0; l < geometry->Nx; l++) {
						T_uw = Interpolate_3D(rh, rh->atmos.T, geometry, st_uw, l, m);
						Planck(1, &T_uw, rh->spectrum.lambda[nspect], &Bnu_uw);
						Planck(1, &rh->atmos.T[local], rh->spectrum.lambda[nspect], &Bnu);

						chi_uw  = Interpolate_3D(rh, chi, geometry, st_uw, l, m);
						dtau_uw = 0.5 * (chi_uw + chi[local]) * st_uw->ds;

						I[0][local] = Bnu - (Bnu_uw - Bnu) / dtau_uw;
						local++;
					}
				}
			} else {
				Error(rh, ERROR_LEVEL_2, routineName,
					  "Boundary condition THERMALIZED not implemented for TOP");
			}
			break;

		case REFLECTIVE:
			Error(rh, ERROR_LEVEL_2, routineName,
				  "Boundary condition REFLECTIVE not implemented");
			break;
		}
		/* --- Assume irradiation is unpolarized --        -------------- */

		for (n = 1; n < 4; n++)
			for (l = (Nz - 1) * Nplane; l < Nz * Nplane; l++) I[n][l] = 0.0;

		if (Psi)
			for (l = (Nz - 1) * Nplane; l < Nz * Nplane; l++) Psi[l] = 0.0;

		dk     = -1;
		kstart = Nz - 2;
		kend   = 0;

		quadrant = geometry->stencil_uw[mu][kstart].quadrant;
	} else {

		/* --- Boundary condition at the TOP of the grid -- ------------- */

		switch (geometry->z_boundary_value[TOP]) {
		case IRRADIATED:
			for (l = 0; l < Nplane; l++) I[0][l] = geometry->Itop[nspect][l];
			break;

		case ZERO:
			for (l = 0; l < Nplane; l++) I[0][l] = 0.0;
			break;

		case THERMALIZED:
			Planck(Nplane, rh->atmos.T, rh->spectrum.lambda[nspect], I[0]);
			break;

		case REFLECTIVE:
			Error(rh, ERROR_LEVEL_2, routineName,
				  "Boundary condition REFLECTIVE not implemented");
			break;
		}
		/* --- Assume irradiation is unpolarized --        -------------- */

		for (n = 1; n < 4; n++)
			for (l = 0; l < Nplane; l++) I[n][l] = 0.0;

		if (Psi)
			for (l = 0; l < Nplane; l++) Psi[l] = 0.0;

		dk     = 1;
		kstart = 1;
		kend   = Nz - 1;

		quadrant = geometry->stencil_dw[mu][kstart].quadrant;
	}
	dl = (quadrant == 1 || quadrant == 4) ? 1 : -1;
	dm = (quadrant == 1 || quadrant == 2) ? 1 : -1;

	/* --- Grid arrays of the cone-preserving interpolation (B-63),
		   with the boundary values of I set -- ------------------------ */

	StokesGrid_init(rh, &sg, nspect, rh->atmos.Nspace, chi, S, I);

	/* --- Figure out in which direction the grids have to be crossed - */

	if (dl == 1) {
		lstart = 0;
		lend   = geometry->Nx - 1;
	} else {
		lstart = geometry->Nx - 1;
		lend   = 0;
	}
	if (dm == 1) {
		mstart = 0;
		mend   = geometry->Ny - 1;
	} else {
		mstart = geometry->Ny - 1;
		mend   = 0;
	}

	for (k = kstart; k != kend + dk; k += dk) {
		if (to_observer) {
			st_uw = &geometry->stencil_dw[mu][k];
			st_dw = &geometry->stencil_uw[mu][k];
		} else {
			st_uw = &geometry->stencil_uw[mu][k];
			st_dw = &geometry->stencil_dw[mu][k];
		}
		if (st_uw->plane == XY || st_uw->plane == XZ) {

			/* --- Inner loop is in the x-direction --       -------------- */

			for (m = mstart; m != mend + dm; m += dm) {
				for (l = lstart; l != lend + dl; l += dl) {
					if (st_uw->plane == XZ && m == mstart)
						SolveLongStokes(rh, geometry, st_uw->longchar, &sg,
										k, l, m, I_uw);
					else
						Interpolate_Stokes_3D(rh, &sg, STOKES_INTENSITY, geometry,
											  st_uw, l, m, I_uw);
					if (rh->input.S_interpolation_stokes == DELO_BEZIER3) {
						Piece_Stokes_Bezier3_3D(rh, geometry, st_uw, st_dw, &sg,
												k, kend, l, m, I_uw, Psi);
					} else {
						Error(rh, ERROR_LEVEL_1, routineName, "Unknown polarization solver: %d",
							  rh->input.S_interpolation_stokes);
					}
				}
			}
		} else {

			/* --- Inner loop is in the y-direction --       -------------- */

			for (l = lstart; l != lend + dl; l += dl) {
				for (m = mstart; m != mend + dm; m += dm) {
					if (l == lstart)
						SolveLongStokes(rh, geometry, st_uw->longchar, &sg,
										k, l, m, I_uw);
					else
						Interpolate_Stokes_3D(rh, &sg, STOKES_INTENSITY, geometry,
											  st_uw, l, m, I_uw);
					if (rh->input.S_interpolation_stokes == DELO_BEZIER3) {
						Piece_Stokes_Bezier3_3D(rh, geometry, st_uw, st_dw, &sg,
												k, kend, l, m, I_uw, Psi);
					} else {
						Error(rh, ERROR_LEVEL_1, routineName, "Unknown polarization solver: %d",
							  rh->input.S_interpolation_stokes);
					}
				}
			}
		}
	}
	StokesGrid_free(&sg);
}
/* ------- end ---------------------------- ShortChar_Stokes.c ------ */

/* ------- begin -------------------------- StokesK_3D_convex.c ----- */

/* --- Where the interpolation weights are convex and the same for
	   every quantity (bilinear: the XZ and YZ planes, and every plane
	   with LINEAR_3D), each component is interpolated on its own as
	   before: eta' is then a convex combination of the eta' of the
	   nodes (weights w_i chi_I,i / sum w_j chi_I,j), and S and I of
	   theirs, inside the cone, and a linear field is reproduced
	   exactly. Only the bicubic interpolation of the XY planes, whose
	   weights have negative lobes and whose clamp acts on each
	   component, goes through the rapidity (stokesgrid.h, B-63) -- */

static inline bool_t convex_weights(RHState *rh, Stencil *st)
{
	return st->plane != XY || rh->input.interpolate_3D == LINEAR_3D;
}

static void StokesK_3D_convex(RHState *rh, StokesGrid *sg, Geometry *geometry, Stencil *st,
							  int l, int m, double chi_I, double K[4][4])
{
	register int i, j;

	int        nspect = sg->nspect;
	ActiveSet *as;

	/* --- Return the elements of the reduced 4x4 Stokes
		   opacity matrix K', which is defined as:

			 =           =         =              =
			 K' = (chi_c*1 + chi_l*Phi) / chi_I - 1,

		 for wavelength# nspect, spatial point k, and ray mu.

	  See: Rees, Murphy, & Durrant, 1989, ApJ 339, 1093-1106.

		   This is the 3-D short characteristics version.
		   --                                            -------------- */
	as = &rh->spectrum.as[nspect];

	for (j = 0; j < 4; j++)
		for (i = 0; i < 4; i++) K[j][i] = 0.0;

	/* --- First the contribution from the active set -- -------------- */

	if (containsPolarized(rh, as)) {
		K[0][1] = Interpolate_3D(rh, as->chi + rh->atmos.Nspace,
								 geometry, st, l, m);
		K[0][2] = Interpolate_3D(rh, as->chi + 2 * rh->atmos.Nspace,
								 geometry, st, l, m);
		K[0][3] = Interpolate_3D(rh, as->chi + 3 * rh->atmos.Nspace,
								 geometry, st, l, m);

		if (rh->input.magneto_optical) {
			K[1][2] = Interpolate_3D(rh, as->chip + 2 * rh->atmos.Nspace,
									 geometry, st, l, m);
			K[1][3] = Interpolate_3D(rh, as->chip + rh->atmos.Nspace,
									 geometry, st, l, m);
			K[2][3] = Interpolate_3D(rh, as->chip,
									 geometry, st, l, m);
		}
	}
	/* --- Add possible contrubution from the background -- ----------- */

	if (rh->atmos.backgrflags[nspect].ispolarized) {
		K[0][1] += Interpolate_3D(rh, as->chi_c + rh->atmos.Nspace,
								  geometry, st, l, m);
		K[0][2] += Interpolate_3D(rh, as->chi_c + 2 * rh->atmos.Nspace,
								  geometry, st, l, m);
		K[0][3] += Interpolate_3D(rh, as->chi_c + 3 * rh->atmos.Nspace,
								  geometry, st, l, m);

		if (rh->input.magneto_optical) {
			K[1][2] += Interpolate_3D(rh, as->chip_c + 2 * rh->atmos.Nspace,
									  geometry, st, l, m);
			K[1][3] += Interpolate_3D(rh, as->chip_c + rh->atmos.Nspace,
									  geometry, st, l, m);
			K[2][3] += Interpolate_3D(rh, as->chip_c,
									  geometry, st, l, m);
		}
	}
	/* --- Divide by Stokes I opacity and fill lower diagonal part -- - */

	for (j = 0; j < 3; j++) {
		for (i = j + 1; i < 4; i++) {
			K[j][i] /= chi_I;
			K[i][j] = K[j][i];
		}
	}
	/* --- Anti-symmetric magneto-optical elements --    -------------- */

	if (rh->input.magneto_optical) {
		K[1][3] *= -1.0;
		K[2][1] *= -1.0;
		K[3][2] *= -1.0;
	}
}
/* ------- end ---------------------------- StokesK_3D_convex.c ----- */

/* ------- begin -------------------------- StokesK_3D.c ------------ */

void StokesK_3D(RHState *rh, StokesGrid *sg, Geometry *geometry, Stencil *st,
				int l, int m, double chi_I, double K[4][4])
{
	/* --- The reduced 4x4 Stokes propagation matrix K',

			 K' = (chi_c*1 + chi_l*Phi) / chi_I - 1,

		   (Rees, Murphy, & Durrant, 1989, ApJ 339, 1093-1106) at the
		   point of stencil st, for the wavelength and ray of the grid
		   sg; chi_I is the interpolated opacity there. With convex
		   weights from the components (StokesK_3D_convex); otherwise
		   eta' = (chi_Q, chi_U, chi_V) / chi_I from the rapidity of
		   the grid, |eta'| < 1 by construction, and rho' with the
		   same weights (stokesgrid.h, B-63). This is the 3-D short
		   characteristics version --                  -------------- */

	double w[3], v[4];

	if (convex_weights(rh, st)) {
		StokesK_3D_convex(rh, sg, geometry, st, l, m, chi_I, K);
		return;
	}
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++) K[j][i] = 0.0;
	if (!sg->polarizedK) return;

	for (int n = 0; n < 3; n++) w[n] = Interpolate_3D_linear(rh, sg->wK[n], geometry, st, l, m);
	Stokes_from_rapidity(1.0, w, v);
	for (int n = 1; n < 4; n++) K[0][n] = K[n][0] = v[n];

	if (sg->magneto_optical) {
		K[1][2] = Interpolate_3D_linear(rh, sg->rK[0], geometry, st, l, m);
		K[3][1] = Interpolate_3D_linear(rh, sg->rK[1], geometry, st, l, m);
		K[2][3] = Interpolate_3D_linear(rh, sg->rK[2], geometry, st, l, m);
		K[2][1] = -K[1][2];
		K[1][3] = -K[3][1];
		K[3][2] = -K[2][3];
	}
}
/* ------- end ---------------------------- StokesK_3D.c ------------ */

/* ------- begin -------------------------- Interpolate_Stokes_3D.c - */

void Interpolate_Stokes_3D(RHState *rh, StokesGrid *sg, enum stokes_vector v,
						   Geometry *geometry, Stencil *st, int l, int m,
						   double out[4])
{
	/* --- The source vector (v = STOKES_SOURCE) or the Stokes vector
		   (STOKES_INTENSITY) of the grid at the point of stencil st:
		   with convex weights each component (Interpolate_3D),
		   otherwise the magnitude with Interpolate_3D and the
		   polarization from the rapidity of the grid with
		   Interpolate_3D_linear: inside the cone by construction
		   either way (stokesgrid.h, B-63) -- ------------------------- */

	double **f  = (v == STOKES_SOURCE) ? sg->S : sg->I,
		   **wg = (v == STOKES_SOURCE) ? sg->wS : sg->wI, w[3];

	if (convex_weights(rh, st)) {
		for (int n = 0; n < 4; n++) out[n] = Interpolate_3D(rh, f[n], geometry, st, l, m);
		return;
	}
	for (int n = 0; n < 3; n++) w[n] = Interpolate_3D_linear(rh, wg[n], geometry, st, l, m);
	Stokes_from_rapidity(Interpolate_3D(rh, f[0], geometry, st, l, m), w, out);
}
/* ------- end ---------------------------- Interpolate_Stokes_3D.c - */
