/* ------- file: -------------------------- bezier_2D.c -------------

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.

   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
			   (Derivatives) Fritsch & Butland (1984),

   1D version Coded by J. de la Cruz Rodriguez (ISP-SU 2017),
   adapted to 2D by Han Uitenbroek

	   Last modified: Tue Sep 23 10:46:46 2025 --

	   --------------------------                      ----------RH-- */

#include <math.h>
#include <string.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "bezier.h"
#include "state.h"

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- Piece_Stokes_Bezier3_2D.c */

void Piece_Stokes_Bezier3_2D(RHState *rh, Geometry *geometry, int nspect, int mu,
							 bool_t to_observer, double *chi, double **S,
							 double **I, double *Psi)
{
	/* --- Cubic DELO-Bezier solver for polarized light
		   Coded by J. de la Cruz Rodriguez (ISP-SU 2017).
		   Adapted to 2-D version by Han Uitenbroek.

		   Reference(s):
		   J. de la Cruz Rodriguez & N. Piskunov (2013)
		   --                                        ------------------ */

	const char   routineName[] = "Piecewise_Bezier3_2D";
	register int l, lp, k, n, m;

	enum boundval  bvalue;
	enum direction upwind, downwind;
	enum sweep     sweep;
	int            Nz = geometry->Nz;
	int            local, dl, lstart, row_left, row_right, dk, kstart, kend;
	double         mux = geometry->mux[mu], I_uw[4], dtau_uw, dtau_dw, c1, c2;
	double         chi_uw, chi_dw, Bnu_uw, Bnu, T_uw;
	double         dchi_uw, dchi_c, dchi_dw, ds03, psi;
	double         K_uw[4][4], K_loc[4][4], K_dw[4][4];
	double         dK_uw[4][4], dK_c[4][4];
	double         S_uw[4], S_loc[4], S_dw[4], dS_uw[4], dS_c[4];
	double         I_loc[4];

	Stencil   *st;
	LongChar  *lc = NULL;
	StokesGrid sg;

	/* --- Boundary conditions in row at the top or bottom -- --------- */

	if (to_observer) {
		row_left = (geometry->Nz - 1) * geometry->Nx;
		bvalue   = geometry->bvalue[BOTTOM];
	} else {
		row_left = 0;
		bvalue   = geometry->bvalue[TOP];
	}
	row_right = row_left + geometry->Nx;

	switch (bvalue) {
	case IRRADIATED:
		if (to_observer) {
			for (l = 0; l < geometry->Nx; l++)
				I[0][row_left + l] = geometry->Ibottom[nspect][l];
		} else {
			for (l = 0; l < geometry->Nx; l++)
				I[0][row_left + l] = geometry->Itop[nspect][l];
		}
		break;

	case ZERO:
		for (l = row_left; l < row_right; l++) I[0][l] = 0.0;
		break;

	case THERMALIZED:
		if (to_observer) {
			local = (geometry->Nz - 1) * geometry->Nx;
			for (l = 0; l < geometry->Nx; l++) {
				st = &geometry->stencil[mu][local];

				T_uw = Interp_2D(rh->atmos.T, st, UPWIND);
				Planck(1, &T_uw, rh->spectrum.lambda[nspect], &Bnu_uw);
				Planck(1, &rh->atmos.T[local], rh->spectrum.lambda[nspect], &Bnu);

				chi_uw  = Interp_2D(chi, st, UPWIND);
				dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[UPWIND];

				I[0][local] = Bnu - (Bnu_uw - Bnu) / dtau_uw;
				local++;
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
	if (Psi)
		for (l = row_left; l < row_right; l++) Psi[l] = 0.0;

	/* --- Assume irradiation is unpolarized --          -------------- */

	for (n = 1; n < 4; n++) {
		for (l = row_left; l < row_right; l++) I[n][l] = 0.0;
	}

	if (to_observer) {
		if (mux >= 0)
			dl = 1;
		else
			dl = -1;
		sweep    = UP;
		dk       = -1;
		kstart   = Nz - 2;
		kend     = 0;
		upwind   = DOWNWIND;
		downwind = UPWIND;
	} else {
		if (mux <= 0)
			dl = 1;
		else
			dl = -1;
		sweep    = DOWN;
		dk       = 1;
		kstart   = 1;
		kend     = Nz - 1;
		upwind   = UPWIND;
		downwind = DOWNWIND;
	}
	/* --- Grid arrays of the cone-preserving interpolation (B-63),
		   with the boundary values of I set -- ------------------------ */

	StokesGrid_init(rh, &sg, nspect, rh->atmos.Nspace, chi, S, I);

	/* --- Finally, go through the grid to solve the non-local transfer
		   equation --                                   -------------- */

	for (k = kstart; k != kend + dk; k += dk) {
		lc     = &geometry->longchar[mu][sweep * geometry->Nz + k];
		lstart = lc->lstart;

		for (lp = 0, l = lstart; lp < geometry->Nx; lp++, l += dl) {
			local = k * geometry->Nx + MODULO(l, geometry->Nx);
			st    = &geometry->stencil[mu][local];

			if (l == lc->lstart && lc->Nlc > 0)
				SolveLongStokes(rh, &sg, lc, local, I_uw);
			else
				Interp_Stokes_2D(&sg, STOKES_INTENSITY, st, upwind, I_uw);

			chi_uw  = Interp_2D(chi, st, upwind);
			dchi_uw = (chi[local] - chi_uw) / st->ds[upwind];
			StokesK_2D(rh, &sg, st, upwind, chi_uw, K_uw);
			Interp_Stokes_2D(&sg, STOKES_SOURCE, st, upwind, S_uw);
			StokesK(rh, nspect, local, chi[local], K_loc);
			for (n = 0; n < 4; n++) S_loc[n] = S[n][local];

			if (k == kend) {

				dtau_uw = 0.5 * (chi[local] + chi_uw) * st->ds[upwind];

				/* --- Piecewise linear at end of ray --       -------------- */

				for (n = 0; n < 4; n++) {
					dS_uw[n] = (S_loc[n] - S_uw[n]) / dtau_uw;
					dS_c[n]  = dS_uw[n];

					for (m = 0; m < 4; m++) {
						dK_uw[n][m] = (K_loc[n][m] - K_uw[n][m]) / dtau_uw;
						dK_c[n][m]  = dK_uw[n][m];
					}
				}
			} else {

				/* --- Bezier piecewise elsewhere --          --------------- */

				chi_dw  = Interp_2D(chi, st, downwind);
				dchi_dw = (chi_dw - chi[local]) / st->ds[downwind];
				StokesK_2D(rh, &sg, st, downwind, chi_dw, K_dw);

				dchi_c = cent_deriv(st->ds[upwind], st->ds[downwind],
									chi_uw, chi[local], chi_dw);

				/* --- Upwind slope of the quadratic Bezier curve through the
				   upwind and local points with the local derivative:
				   2 x secant - local derivative (N-24). The secant alone
				   is exact only for linear functions -- -------------- */

				dchi_uw = 2.0 * dchi_uw - dchi_c;

				/* --- Upwind optical path length --           -------------- */

				ds03    = st->ds[upwind] * 0.3333333333333333333333333;
				c1      = chi_uw + ds03 * dchi_uw;
				c2      = chi[local] - ds03 * dchi_c;
				dtau_uw = st->ds[upwind] * (chi[local] + chi_uw + c1 + c2) * 0.25;

				/* --- Downwind optical path length --          ------------- */

				ds03    = st->ds[downwind] * 0.3333333333333333333333333;
				c1      = chi_dw - ds03 * dchi_dw;
				c2      = chi[local] + ds03 * dchi_c;
				dtau_dw = st->ds[downwind] * (chi[local] + chi_dw + c1 + c2) * 0.25;

				Interp_Stokes_2D(&sg, STOKES_SOURCE, st, downwind, S_dw);
				for (n = 0; n < 4; n++) {
					dS_uw[n] = (S_loc[n] - S_uw[n]) / dtau_uw;

					for (m = 0; m < 4; m++)
						dK_uw[n][m] = (K_loc[n][m] - K_uw[n][m]) / dtau_uw;
				}

				cent_deriv_vec(dS_c, dtau_uw, dtau_dw, S_uw, S_loc, S_dw);
				cent_deriv_mat(dK_c, dtau_uw, dtau_dw, K_uw, K_loc, K_dw);

				/* --- Quadratic Bezier upwind slopes (N-24) -- ------- */

				for (n = 0; n < 4; n++) {
					dS_uw[n] = 2.0 * dS_uw[n] - dS_c[n];
					for (m = 0; m < 4; m++)
						dK_uw[n][m] = 2.0 * dK_uw[n][m] - dK_c[n][m];
				}
			}
			psi = Bezier3_Stokes_step(dtau_uw, K_uw, K_loc, dK_uw, dK_c, S_uw,
									  S_loc, dS_uw, dS_c, I_uw, I_loc);
			if (Psi) Psi[local] = psi;
			for (n = 0; n < 4; n++) I[n][local] = I_loc[n];
			StokesGrid_update_I(&sg, local);
		}
	}
	StokesGrid_free(&sg);
}
/* ------- end ------------------------- Piece_Stokes_Bezier3_2D.c -- */

/* ------- begin -------------------------- Piecewise_Bezier3_2D.c -- */

void Piecewise_Bezier3_2D(RHState *rh, Geometry *geometry, int nspect, int mu,
						  bool_t to_observer, double *chi, double *S,
						  double *I, double *Psi)
{

	/* --- Cubic short-char Bezier solver.
		   Coded by J. de la Cruz Rodriguez (ISP-SU 2017).
		   Adapted to 2-D version by Han Uitenbroek.

		   Reference(s):
		   J. de la Cruz Rodriguez & N. Piskunov (2013)
		   --                                        ------------------ */

	const char   routineName[] = "Piecewise_Bezier3_2D";
	register int l, k, n;

	enum boundval  bvalue;
	enum direction upwind, downwind;
	enum sweep     sweep;
	int            Nx = geometry->Nx, Nz = geometry->Nz, local, dl, lstart;
	int            Nrow, row_left, row_right, dk, kstart, kend;
	double         mux = geometry->mux[mu], I_uw, *Iboundary, dtau_uw, dtau_dw;
	double         c1, c2, dS_uw, S_uw, chi_uw, S_dw, chi_dw, dt03;
	double         eps = 0, alpha = 0, beta = 0, gamma = 0, theta = 0;
	double         dchi_dw, dchi_uw, dchi_c, dS_c, ds03;
	double         Bnu_uw, Bnu, T_uw;

	Stencil  *st;
	LongChar *lc = NULL;

	/* --- In case of FIXED left and right boundary conditions treat
		   the left and right columns --                 -------------- */

	if (geometry->hboundary == FIXED) {
		if ((to_observer && mux > 0.0) || (!to_observer && mux < 0.0)) {

			/* --- Rays that begin at the LEFT boundary and end at the RIGHT
			   boundary --                               -------------- */

			for (k = 0; k < Nz; k++) I[k * Nx] = geometry->Ileft[nspect][k];
			if (Psi)
				for (k = 0; k < Nz; k++) Psi[k * Nx] = 0.0;
		} else if ((to_observer && mux < 0.0) || (!to_observer && mux > 0.0)) {

			/* --- Rays that begin at the RIGHT boundary and end at the
				   LEFT boundary --                          -------------- */

			for (k = 0; k < Nz; k++)
				I[(k + 1) * Nx - 1] = geometry->Iright[nspect][k];
			if (Psi)
				for (k = 0; k < Nz; k++) Psi[(k + 1) * Nx - 1] = 0.0;
		}
	}
	/* --- Boundary conditions in row at the top or bottom -- --------- */

	if (to_observer) {
		row_left = (Nz - 1) * Nx;
		bvalue   = geometry->bvalue[BOTTOM];
	} else {
		row_left = 0;
		bvalue   = geometry->bvalue[TOP];
	}
	row_right = row_left + Nx;

	switch (bvalue) {
	case IRRADIATED:
		if (to_observer)
			Iboundary = geometry->Ibottom[nspect];
		else
			Iboundary = geometry->Itop[nspect];
		for (l = 0; l < Nx; l++) I[row_left + l] = Iboundary[l];
		break;

	case ZERO:
		for (l = row_left; l < row_right; l++) I[l] = 0.0;
		break;

	case THERMALIZED:
		if (to_observer) {
			local = (Nz - 1) * Nx;
			for (l = 0; l < Nx; l++) {
				st = &geometry->stencil[mu][local];

				T_uw = Interp_2D(rh->atmos.T, st, UPWIND);
				Planck(1, &T_uw, rh->spectrum.lambda[nspect], &Bnu_uw);
				Planck(1, &rh->atmos.T[local], rh->spectrum.lambda[nspect], &Bnu);

				chi_uw  = Interp_2D(chi, st, UPWIND);
				dtau_uw = 0.5 * (chi_uw + chi[local]) * st->ds[UPWIND];

				I[local] = Bnu - (Bnu_uw - Bnu) / dtau_uw;
				local++;
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

	if (Psi)
		for (l = row_left; l < row_right; l++) Psi[l] = 0.0;

	/* --- Finally, go through the grid to solve the non-local transfer
		   equation --                                   -------------- */

	if (to_observer) {
		if (mux >= 0)
			dl = 1;
		else
			dl = -1;
		sweep    = UP;
		dk       = -1;
		kstart   = Nz - 2;
		kend     = 0;
		upwind   = DOWNWIND;
		downwind = UPWIND;
	} else {
		if (mux <= 0)
			dl = 1;
		else
			dl = -1;
		sweep    = DOWN;
		dk       = 1;
		kstart   = 1;
		kend     = Nz - 1;
		upwind   = UPWIND;
		downwind = DOWNWIND;
	}

	for (k = kstart; k != kend + dk; k += dk) {

		if (geometry->hboundary == PERIODIC) {
			lc     = &geometry->longchar[mu][sweep * Nz + k];
			lstart = lc->lstart;
			Nrow   = Nx;
		} else {
			if (mux == 0.0) {
				lstart = 0;
				Nrow   = Nx;
			} else {
				if ((to_observer && mux > 0.0) || (!to_observer && mux < 0.0))
					lstart = 1;
				else
					lstart = Nx - 2;

				Nrow = Nx - 1;
			}
		}

		for (n = 0, l = lstart; n < Nrow; n++, l += dl) {
			if (geometry->hboundary == PERIODIC)
				local = k * Nx + MODULO(l, Nx);
			else
				local = k * Nx + l;
			st = &geometry->stencil[mu][local];

			if (geometry->hboundary == PERIODIC &&
				l == lc->lstart && lc->Nlc > 0)
				I_uw = SolveLong(rh, lc, local, chi, S, I);
			else
				I_uw = Interp_2D(I, st, upwind);

			chi_uw  = Interp_2D(chi, st, upwind);
			dchi_uw = (chi[local] - chi_uw) / st->ds[upwind];
			S_uw    = Interp_2D(S, st, upwind);

			if (k == kend ||
				(geometry->hboundary == FIXED && n == Nrow - 1)) {

				dtau_uw = 0.5 * (chi[local] + chi_uw) * st->ds[upwind];

				/* --- Piecewise linear at end of ray --       -------------- */

				dS_c = (S[local] - S_uw) / dtau_uw;
			} else {

				/* --- Bezier integration --                   -------------- */

				chi_dw  = Interp_2D(chi, st, downwind);
				dchi_dw = (chi_dw - chi[local]) / st->ds[downwind];

				dchi_c = cent_deriv(st->ds[upwind], st->ds[downwind],
									chi_uw, chi[local], chi_dw);

				/* --- Upwind slope of the quadratic Bezier curve through the
				   upwind and local points with the local derivative:
				   2 x secant - local derivative (N-24). The secant alone
				   is exact only for linear functions -- -------------- */

				dchi_uw = 2.0 * dchi_uw - dchi_c;

				/* --- Upwind optical path length --          --------------- */

				ds03    = st->ds[upwind] * 0.333333333333333333333333333;
				c1      = chi[local] - ds03 * dchi_c;
				c2      = chi_uw + ds03 * dchi_uw;
				dtau_uw = st->ds[upwind] * (chi[local] + chi_uw + c1 + c2) * 0.25;

				/* --- Downwind optical path length --          ------------- */

				ds03    = st->ds[downwind] * 0.333333333333333333333333333;
				c1      = chi[local] + ds03 * dchi_c;
				c2      = chi_dw - ds03 * dchi_dw;
				dtau_dw = st->ds[downwind] * (chi[local] + chi_dw + c1 + c2) * 0.25;

				/* --- dS/dt at central point --               -------------- */

				S_dw = Interp_2D(S, st, downwind);
				dS_c = cent_deriv(dtau_uw, dtau_dw, S_uw, S[local], S_dw);
			}
			/* --- Source function control points --       ---------------- */

			dt03  = dtau_uw * 0.333333333333333333333333333;
			dS_uw = 2.0 * (S[local] - S_uw) / dtau_uw - dS_c; /* N-24 */
			c1    = S[local] - dt03 * dS_c;
			c2    = S_uw + dt03 * dS_uw;

			/* --- Compute interpolation parameters --     ---------------- */

			Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);

			/* --- Solve integral in this interval --      ---------------- */

			I[local] = I_uw * eps + beta * S[local] + alpha * S_uw +
					   theta * c1 + gamma * c2;

			/* --- Diagonal operator --                    ---------------- */

			if (Psi) Psi[local] = Bezier3_psi(beta, theta, S[local], c1); /* B-66 */
		}
	}
}
/* ------- end ---------------------------- Piecewise_Bezier3_2D.c -- */
