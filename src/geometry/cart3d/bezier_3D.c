/* ------- file: -------------------------- bezier_3D.c -- ----------

	   Version:       rh2.0, 3-D Cartesian, short characteristics

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.

   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
			   (Derivatives) Fritsch & Butland (1984),

   1D version Coded by J. de la Cruz Rodriguez (ISP-SU 2017),
   adapted to 3D by Han Uitenbroek

	   Last modified: Mon Jan  6 14:59:35 2020 --

	   --------------------------                      ----------RH-- */

#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "inputs.h"
#include "error.h"
#include "bezier.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern char messageStr[];

/* ------- begin -------------------------- Piece_Stokes_Bezier3_3D.c */

void Piece_Stokes_Bezier3_3D(Geometry *geometry, Stencil *st_uw,
							 Stencil *st_dw, StokesGrid *sg,
							 int k, int kend, int l, int m,
							 double *I_uw, double *Psi)
{
	/* --- Cubic short-char Bezier Polarization solver.
		   Coded by J. de la Cruz Rodriguez (ISP-SU 2017).
		   Adapted to 3-D version by Han Uitenbroek.

		   Reference(s):
		   J. de la Cruz Rodriguez & N. Piskunov (2013)
		   --                                        ------------------ */

	register int n, p;

	int    local;
	double chi_uw, chi_dw, S_uw[4], S_dw[4], dS_uw[4];
	double dtau_uw, dtau_dw, S_loc[4], dS_c[4], dchi_uw, dchi_dw, dchi_c;
	double K_loc[4][4], K_uw[4][4], K_dw[4][4], c1, c2, ds03, psi;
	double dK_uw[4][4], dK_c[4][4];
	double I_loc[4];

	int     nspect = sg->nspect;
	double *chi = sg->chi, **S = sg->S, **I = sg->I;

	local = k * geometry->Nplane + m * geometry->Nx + l;

	/* --- The upwind quantities --                      -------------- */

	chi_uw  = Interpolate_3D(chi, geometry, st_uw, l, m);
	dchi_uw = (chi[local] - chi_uw) / st_uw->ds;
	StokesK_3D(sg, geometry, st_uw, l, m, chi_uw, K_uw);
	Interpolate_Stokes_3D(sg, STOKES_SOURCE, geometry, st_uw, l, m, S_uw);
	StokesK(nspect, local, chi[local], K_loc);
	for (n = 0; n < 4; n++) S_loc[n] = S[n][local];

	if (k == kend) {

		/* --- Piecewise linear at end of ray --           -------------- */

		dtau_uw = 0.5 * (chi[local] + chi_uw) * st_uw->ds;

		for (n = 0; n < 4; n++) {
			dS_uw[n] = (S_loc[n] - S_uw[n]) / dtau_uw;
			dS_c[n]  = dS_uw[n];

			for (p = 0; p < 4; p++) {
				dK_uw[n][p] = (K_loc[n][p] - K_uw[n][p]) / dtau_uw;
				dK_c[n][p]  = dK_uw[n][p];
			}
		}
	} else {

		/* --- Bezier piecewise elsewhere --               -------------- */

		chi_dw  = Interpolate_3D(chi, geometry, st_dw, l, m);
		dchi_dw = (chi_dw - chi[local]) / st_dw->ds;
		StokesK_3D(sg, geometry, st_dw, l, m, chi_dw, K_dw);

		dchi_c = cent_deriv(st_uw->ds, st_dw->ds,
							chi_uw, chi[local], chi_dw);

		/* --- Upwind slope of the quadratic Bezier curve through the
		   upwind and local points with the local derivative:
		   2 x secant - local derivative (N-24). The secant alone
		   is exact only for linear functions -- -------------- */

		dchi_uw = 2.0 * dchi_uw - dchi_c;

		/* --- Upwind optical path length --               -------------- */

		ds03    = st_uw->ds * 0.3333333333333333333333333;
		c1      = chi_uw + ds03 * dchi_uw;
		c2      = chi[local] - ds03 * dchi_c;
		dtau_uw = st_uw->ds * (chi[local] + chi_uw + c1 + c2) * 0.25;

		/* --- Downwind optical path length --             ------------- */

		ds03    = st_dw->ds * 0.3333333333333333333333333;
		c1      = chi_dw - ds03 * dchi_dw;
		c2      = chi[local] + ds03 * dchi_c;
		dtau_dw = st_dw->ds * (chi[local] + chi_dw + c1 + c2) * 0.25;

		Interpolate_Stokes_3D(sg, STOKES_SOURCE, geometry, st_dw, l, m, S_dw);
		for (n = 0; n < 4; n++) {
			dS_uw[n] = (S_loc[n] - S_uw[n]) / dtau_uw;

			for (p = 0; p < 4; p++)
				dK_uw[n][p] = (K_loc[n][p] - K_uw[n][p]) / dtau_uw;
		}

		cent_deriv_vec(dS_c, dtau_uw, dtau_dw, S_uw, S_loc, S_dw);
		cent_deriv_mat(dK_c, dtau_uw, dtau_dw, K_uw, K_loc, K_dw);

		/* --- Quadratic Bezier upwind slopes (N-24) --    -------------- */

		for (n = 0; n < 4; n++) {
			dS_uw[n] = 2.0 * dS_uw[n] - dS_c[n];
			for (p = 0; p < 4; p++)
				dK_uw[n][p] = 2.0 * dK_uw[n][p] - dK_c[n][p];
		}
	}
	psi = Bezier3_Stokes_step(dtau_uw, K_uw, K_loc, dK_uw, dK_c, S_uw,
							  S_loc, dS_uw, dS_c, I_uw, I_loc);
	if (Psi) Psi[local] = psi;
	for (n = 0; n < 4; n++) I[n][local] = I_loc[n];
	StokesGrid_update_I(sg, local);
}
/* ------- end ---------------------------- Piece_Stokes_Bezier3_3D.c */

/* ------- begin -------------------------- Piecewise_Bezier3_3D.c -- */

void Piecewise_Bezier3_3D(Geometry *geometry, Stencil *st_uw, Stencil *st_dw,
						  int k, int kend, int l, int m, double I_uw,
						  double *chi, double *S, double *I, double *Psi)
{
	/* --- Cubic short-char Bezier solver.
		   Coded by J. de la Cruz Rodriguez (ISP-SU 2017).
		   Adapted to 3-D version by Han Uitenbroek.

		   Reference(s):
		   J. de la Cruz Rodriguez & N. Piskunov (2013)
		   --                                        ------------------ */

	int    local;
	double dtau_uw, dtau_dw;
	double c1, c2, dS_uw, S_uw, chi_uw, S_dw, chi_dw, dt03;
	double eps = 0, alpha = 0, beta = 0, gamma = 0, theta = 0;
	double dchi_dw, dchi_uw, dchi_c, dS_c, ds03;

	local = k * geometry->Nplane + m * geometry->Nx + l;

	/* --- The upwind quantities --                      -------------- */

	chi_uw  = Interpolate_3D(chi, geometry, st_uw, l, m);
	dchi_uw = (chi[local] - chi_uw) / st_uw->ds;
	S_uw    = Interpolate_3D(S, geometry, st_uw, l, m);

	if (k == kend) {

		dtau_uw = 0.5 * st_uw->ds * (chi[local] + chi_uw);

		/* --- Piecewise linear at end of ray --           -------------- */

		dS_c = (S[local] - S_uw) / dtau_uw;
	} else {

		/* --- Bezier interpolation --                     -------------- */

		chi_dw  = Interpolate_3D(chi, geometry, st_dw, l, m);
		dchi_dw = (chi_dw - chi[local]) / st_dw->ds;

		dchi_c = cent_deriv(st_uw->ds, st_dw->ds,
							chi_uw, chi[local], chi_dw);

		/* --- Upwind slope of the quadratic Bezier curve through the
		   upwind and local points with the local derivative:
		   2 x secant - local derivative (N-24). The secant alone
		   is exact only for linear functions -- -------------- */

		dchi_uw = 2.0 * dchi_uw - dchi_c;

		/* --- Upwind optical path length --             ---------------- */

		ds03    = st_uw->ds * 0.333333333333333333333333333;
		c1      = chi[local] - ds03 * dchi_c;
		c2      = chi_uw + ds03 * dchi_uw;
		dtau_uw = st_uw->ds * (chi[local] + chi_uw + c1 + c2) * 0.25;

		/* --- Downwind optical path length --             ------------- */

		ds03    = st_dw->ds * 0.333333333333333333333333333;
		c1      = chi[local] + ds03 * dchi_c;
		c2      = chi_dw - ds03 * dchi_dw;
		dtau_dw = st_dw->ds * (chi[local] + chi_dw + c1 + c2) * 0.25;
		S_dw    = Interpolate_3D(S, geometry, st_dw, l, m);

		/* --- dS/dt at central point --                   -------------- */

		dS_c = cent_deriv(dtau_uw, dtau_dw, S_uw, S[local], S_dw);
	}
	/* --- Source function control points --             -------------- */

	dt03  = dtau_uw * 0.333333333333333333333333333;
	dS_uw = 2.0 * (S[local] - S_uw) / dtau_uw - dS_c; /* N-24 */
	c1    = S[local] - dt03 * dS_c;
	c2    = S_uw + dt03 * dS_uw;

	/* --- Compute interpolation parameters --           -------------- */

	Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);

	/* --- Solve integral in this interval --            -------------- */

	I[local] = I_uw * eps + beta * S[local] + alpha * S_uw +
			   theta * c1 + gamma * c2;

	/* --- Diagonal operator --                          -------------- */

	if (Psi) Psi[local] = Bezier3_psi(beta, theta, S[local], c1); /* B-66 */
}
/* ------- end ---------------------------- Piecewise_Bezier3_3D.c -- */
