/* ------- file: -------------------------- bezier_1D.c -------------

   Cubic DELO-Bezier (polarized) and cubic short-char Bezier solvers.

   References: de la Cruz Rodriguez & Piskunov (2013), Auer (2003)
			   (Derivatives) Fritsch & Butland (1984),

   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

   Modifications:
		   2017-03-12, JdlCR: Created!
   Adapted to RH:
		   Han Uitenbroek

	   Last modified: Mon Jan  6 15:46:35 2020 --

	   --------------------------                      ----------RH-- */

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "bezier.h"
#include "state.h"

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- Piece_Stokes_Bezier3_1D.c */

void Piece_Stokes_Bezier3_1D(int nspect, int mu, bool_t to_obs,
							 double *chi, double **S, double **I,
							 double *Psi)
{
	/* --- Cubic DELO-Bezier solver for polarized light
		   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

		   Reference(s):
		   J. de la Cruz Rodriguez & N. Piskunov (2013)
		   --                                        ------------------ */

	static const char routineName[] = "Piece_Stokes_Bezier3_1D";
	static const int  siz_mat       = 16 * sizeof(double),
					 siz_vec        = 4 * sizeof(double);

	register int k, n, m, i;

	int    Ndep = rh->geometry->Ndep, k_start, k_end, dk, k_last;
	double dtau_uw, dtau_dw = 0.0, c1, c2, dsdn2, dchi_dn = 0.0,
					I_upw[4], Bnu[2];
	double  dchi_up, dchi_c, psi, I_loc[4];
	double  dsup, dsdn;
	double  Ku[4][4], K0[4][4], Kd[4][4], dKu[4][4], dK0[4][4];
	double  Su[4], S0[4], Sd[4], dSu[4], dS0[4];
	double  imu = 1.0 / rh->geometry->muz[mu];
	double *z   = rh->geometry->height;

	if (to_obs) {
		dk      = -1;
		k_start = Ndep - 1;
		k_end   = 0;
	} else {
		dk      = 1;
		k_start = 0;
		k_end   = Ndep - 1;
	}
	dtau_uw = 0.5 * imu * (chi[k_start] + chi[k_start + dk]) *
			  fabs(z[k_start] - z[k_start + dk]);

	/* --- Boundary conditions --                        -------------- */

	if (to_obs) {
		switch (rh->geometry->vboundary[BOTTOM]) {
		case ZERO:
			for (n = 0; n < 4; n++) I_upw[n] = 0.0;
			break;
		case THERMALIZED:
			Planck(2, &rh->atmos.T[Ndep - 2], rh->spectrum.lambda[nspect], Bnu);
			I_upw[0] = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
			for (n = 1; n < 4; n++) I_upw[n] = 0.0;
			break;
		case IRRADIATED:
			I_upw[0] = rh->geometry->Ibottom[nspect][mu];
			for (n = 1; n < 4; n++) I_upw[n] = 0.0;
			break;
		case REFLECTIVE:
			Error(ERROR_LEVEL_2, routineName, "Boundary condition not implemented: %d",
				  rh->geometry->vboundary[BOTTOM]);
		}
	} else {
		switch (rh->geometry->vboundary[TOP]) {
		case ZERO:
			for (n = 0; n < 4; n++) I_upw[n] = 0.0;
			break;
		case IRRADIATED:
			I_upw[0] = rh->geometry->Itop[nspect][mu];
			for (n = 1; n < 4; n++) I_upw[n] = 0.0;
			break;
		default:
			Error(ERROR_LEVEL_2, routineName, "Boundary condition not implemented: %d",
				  rh->geometry->vboundary[TOP]);
		}
	}

	for (n = 0; n < 4; n++) I[n][k_start] = I_upw[n];
	if (Psi) Psi[k_start] = 0.0;

	k       = k_start + dk;
	dsup    = fabs(z[k] - z[k - dk]) * imu;
	dsdn    = fabs(z[k + dk] - z[k]) * imu;
	dchi_up = (chi[k] - chi[k - dk]) / dsup;

	/* ---  dchi/ds at central point--               ------------------ */

	dchi_c = cent_deriv(dsup, dsdn, chi[k - dk], chi[k], chi[k + dk]);

	/* --- Upwind path_length (BEzier3 integration) -- ---------------- */

	c2 = chi[k] - (dsup * 0.3333333333333333333333333) * dchi_c;
	c1 = chi[k - dk] + (dsup * 0.3333333333333333333333333) * dchi_up;

	dtau_uw = 0.25 * dsup * (chi[k] + chi[k - dk] + c1 + c2);

	/* --- Ku, K0 and dKu, dSu -                     ------------------ */

	StokesK(nspect, k_start, chi[k_start], Ku);
	StokesK(nspect, k_start + dk, chi[k_start + dk], K0);

	Svec(k_start, S, Su);
	Svec(k_start + dk, S, S0);

	/* --- Assume side derivative in the first interval -- ------------ */

	for (n = 0; n < 4; n++) {
		dSu[n] = (S0[n] - Su[n]) / dtau_uw;

		for (m = 0; m < 4; m++)
			dKu[n][m] = (K0[n][m] - Ku[n][m]) / dtau_uw;
	}

	/* --- Solve transfer along ray --                   -------------- */

	k_last = k_end + dk;

	for (k = k_start + dk; k != k_last; k += dk) {

		/* --- dchi/ds at downwind point --                -------------- */

		if (k != k_end) {

			dsdn = fabs(z[k + dk] - z[k]) * imu;

			if (abs(k - k_end) > 1) {
				dsdn2   = fabs(z[k + 2 * dk] - z[k + dk]) * imu;
				dchi_dn = cent_deriv(dsdn, dsdn2, chi[k], chi[k + dk],
									 chi[k + 2 * dk]);
			} else
				dchi_dn = (chi[k + dk] - chi[k]) / dsdn;

			/* --- Make sure that c1 and c2 don't do below zero -- -------- */

			c1 = (chi[k + dk] - (dsdn * 0.333333333333333333333333333) * dchi_dn);
			c2 = (chi[k] + (dsdn * 0.333333333333333333333333333) * dchi_c);

			/* --- Bezier3 integrated dtau --            ------------------ */

			dtau_dw = 0.25 * dsdn * (chi[k] + chi[k + dk] + c1 + c2);

			/* ---- get algebra in place --              ------------------ */

			StokesK(nspect, k + dk, chi[k + dk], Kd);
			Svec(k + dk, S, Sd);

			cent_deriv_mat(dK0, dtau_uw, dtau_dw, Ku, K0, Kd);
			cent_deriv_vec(dS0, dtau_uw, dtau_dw, Su, S0, Sd);
		} else {

			/* --- Last interval, assume linear dependence for the derivatives
			   at the central point. In that case all the info related to the
			   downwind point can be ditched --          -------------- */

			for (n = 0; n < 4; n++) {
				dS0[n] = (S0[n] - Su[n]) / dtau_uw;
				for (m = 0; m < 4; m++) dK0[n][m] = (K0[n][m] -
													 Ku[n][m]) /
													dtau_uw;
			}
		}

		for (n = 0; n < 4; n++) I_upw[n] = I[n][k - dk];
		psi = Bezier3_Stokes_step(dtau_uw, Ku, K0, dKu, dK0, Su, S0, dSu, dS0,
								  I_upw, I_loc);
		if (Psi) Psi[k] = psi;
		for (i = 0; i < 4; i++) I[i][k] = I_loc[i];

		/* --- Shift values for next depth --          ------------------ */

		memcpy(Su, S0, siz_vec);
		memcpy(S0, Sd, siz_vec);
		memcpy(dSu, dS0, siz_vec);

		memcpy(Ku[0], K0[0], siz_mat);
		memcpy(K0[0], Kd[0], siz_mat);
		memcpy(dKu[0], dK0[0], siz_mat);

		dtau_uw = dtau_dw;
		dsup    = dsdn;
		dchi_up = dchi_c;
		dchi_c  = dchi_dn;
	}
}
/* ------- end ------------------------- Piece_Stokes_Bezier3_1D.c -- */

/* ------- begin ----------------------- steffenCoef.c -------------- */

static void steffenCoef(double odx, double dx, double yu, double y0,
						double yd, double c[3])
{
	/* --- The coefficients of the Steffen (1990) derivative of
		   cent_deriv on (yu, y0, yd): on each branch of its limiter
		   the derivative is linear in the three values, y' = c[0] yu
		   + c[1] y0 + c[2] yd (N-51, tridiagonal operator). The branch
		   is chosen as cent_deriv chooses it --        -------------- */

	const double S0 = (yd - y0) / dx;
	const double Su = (y0 - yu) / odx;
	const double P0 = fabs((Su * dx + S0 * odx) / (odx + dx)) * 0.5;
	const double m  = fmin(fabs(Su), fmin(fabs(S0), P0));

	c[0] = c[1] = c[2] = 0.0;
	if (signFortran(S0) + signFortran(Su) == 0.0) return;
	if (m == fabs(Su)) { /* y' = 2 Su */
		c[0] = -2.0 / odx;
		c[1] = 2.0 / odx;
	} else if (m == fabs(S0)) { /* y' = 2 S0 */
		c[1] = -2.0 / dx;
		c[2] = 2.0 / dx;
	} else { /* y' = (Su dx + S0 odx) / (odx + dx) */
		c[0] = -(dx / odx) / (odx + dx);
		c[2] = (odx / dx) / (odx + dx);
		c[1] = -c[0] - c[2];
	}
}
/* ------- end ------------------------- steffenCoef.c -------------- */

/* ------- begin ----------------------- Piecewise_Bezier3_1D.c ----- */

void Piecewise_Bezier3_1D(int nspect, int mu, bool_t to_obs,
						  double *chi, double *S, double *I, double *Psi)
{
	Piecewise_Bezier3_1D_tri(nspect, mu, to_obs, chi, S, I, Psi, NULL);
}
/* ------- end ------------------------- Piecewise_Bezier3_1D.c ----- */

/* ------- begin ----------------------- Piecewise_Bezier3_1D_tri.c - */

void Piecewise_Bezier3_1D_tri(int nspect, int mu, bool_t to_obs,
							  double *chi, double *S, double *I, double *Psi,
							  double *Tri)
{

	/* --- Cubic Bezier solver for unpolarized light
		   Coded by J. de la Cruz Rodriguez (ISP-SU 2017)

		   Reference:
		   J. de la Cruz Rodriguez & N. Piskunov (2013), Auer (2003)

		   Comments:
			 JdlCR: We only check that the control points of the opacity
				and source function are always above zero to avoid having
				a negative interpolant.
		   --                                            -------------- */

	register int k;
	const char   routineName[] = "Piecewise_Bezier3_1D";

	int    k_start, k_end, dk, Ndep = rh->geometry->Ndep;
	double dtau_uw, dtau_dw = 0.0, I_upw = 0.0, c1, c2, zmu, Bnu[2];
	double dsup, dsdn, dt03, eps = 0, alpha = 0, beta = 0, gamma = 0, theta = 0;
	double dS_up, dS_c, dchi_up, dchi_c, dchi_dn = 0.0, dsdn2;

	/* --- Tri (N-51): with it, the exact row of the operator of this
		   ray, Tri[3k] = dI_k/dS_k-dk (upwind), Tri[3k + 1] = dI_k/dS_k
		   (not limited as Psi is, B-66), Tri[3k + 2] = dI_k/dS_k+dk.
		   cSc and cSup are the coefficients of dS_c and dS_up on
		   (S_k-dk, S_k, S_k+dk) of their own centre; T0, Tp those of
		   the upwind point --                          -------------- */

	double cSc[3] = {0.0, 0.0, 0.0}, cSup[3], T0 = 0.0, Tp = 0.0;

	zmu = 1.0 / rh->geometry->muz[mu];

	/* --- Distinguish between rays going from BOTTOM to TOP
		   (to_obs == TRUE), and vice versa --           -------------- */

	if (to_obs) {
		dk      = -1;
		k_start = Ndep - 1;
		k_end   = 0;
	} else {
		dk      = 1;
		k_start = 0;
		k_end   = Ndep - 1;
	}

	dtau_uw = 0.5 * zmu * (chi[k_start] + chi[k_start + dk]) *
			  fabs(rh->geometry->height[k_start] - rh->geometry->height[k_start + dk]);

	/* --- Boundary conditions --                        -------------- */

	if (to_obs) {
		switch (rh->geometry->vboundary[BOTTOM]) {
		case ZERO:
			I_upw = 0.0;
			break;
		case THERMALIZED:
			Planck(2, &rh->atmos.T[Ndep - 2], rh->spectrum.lambda[nspect], Bnu);
			I_upw = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
			break;
		case IRRADIATED:
			I_upw = rh->geometry->Ibottom[nspect][mu];
			break;
		case REFLECTIVE:
			Error(ERROR_LEVEL_2, routineName, "Boundary condition not implemented: %d",
				  rh->geometry->vboundary[BOTTOM]);
		}
	} else {
		switch (rh->geometry->vboundary[TOP]) {
		case ZERO:
			I_upw = 0.0;
			break;
		case THERMALIZED:
			Planck(2, &rh->atmos.T[0], rh->spectrum.lambda[nspect], Bnu);
			I_upw = Bnu[0] - (Bnu[1] - Bnu[0]) / dtau_uw;
			break;
		case IRRADIATED:
			I_upw = rh->geometry->Itop[nspect][mu];
			break;
		case REFLECTIVE:
			Error(ERROR_LEVEL_2, routineName, "Boundary condition not implemented: %d",
				  rh->geometry->vboundary[BOTTOM]);
		}
	}

	I[k_start] = I_upw;
	if (Psi) Psi[k_start] = 0.0;
	if (Tri) Tri[3 * k_start] = Tri[3 * k_start + 1] = Tri[3 * k_start + 2] = 0.0;

	/* --- Set variables for first iteration to allow simple
		   shift for all next iterations --              -------------- */

	k       = k_start + dk;
	dsup    = fabs(rh->geometry->height[k] - rh->geometry->height[k - dk]) * zmu;
	dsdn    = fabs(rh->geometry->height[k + dk] - rh->geometry->height[k]) * zmu;
	dchi_up = (chi[k] - chi[k - dk]) / dsup;

	/* --- dchi/ds at central point --                   -------------- */

	dchi_c = cent_deriv(dsup, dsdn, chi[k - dk], chi[k], chi[k + dk]);

	/* --- upwind path_length (Bezier3 integration) --   -------------- */

	c1 = (chi[k] - (dsup * 0.333333333333333333) * dchi_c);
	c2 = (chi[k - dk] + (dsup * 0.333333333333333333) * dchi_up);

	dtau_uw = dsup * (chi[k] + chi[k - dk] + c1 + c2) * 0.25;

	/* --- dS/dtau at upwind point --                    -------------- */

	dS_up   = (S[k] - S[k - dk]) / dtau_uw;
	cSup[0] = 0.0;
	cSup[1] = -1.0 / dtau_uw;
	cSup[2] = 1.0 / dtau_uw;

	/* --- Solve transfer along ray --                   -------------- */

	for (k = k_start + dk; k != k_end + dk; k += dk) {

		if (k != k_end) {

			/* --- Downwind path length --                   -------------- */

			dsdn = fabs(rh->geometry->height[k + dk] - rh->geometry->height[k]) * zmu;

			/* --- dchi/ds at downwind point --             -------------- */

			if (abs(k - k_end) > 1) {
				dsdn2 = fabs(rh->geometry->height[k + 2 * dk] -
							 rh->geometry->height[k + dk]) *
						zmu;
				dchi_dn = cent_deriv(dsdn, dsdn2, chi[k], chi[k + dk], chi[k + 2 * dk]);
			} else {
				dchi_dn = (chi[k + dk] - chi[k]) / dsdn;
			}

			/* --- Make sure that c1 and c2 don't go below zero -- ------- */

			c1 = (chi[k] + (dsdn * 0.3333333333333333333) * dchi_c);
			c2 = (chi[k + dk] - (dsdn * 0.3333333333333333333) * dchi_dn);

			/* --- Downwind optical path length --          -------------- */

			dtau_dw = dsdn * (chi[k] + chi[k + dk] + c1 + c2) * 0.25;

			/* --- dS/dt at central point --                -------------- */

			dS_c = cent_deriv(dtau_uw, dtau_dw, S[k - dk], S[k], S[k + dk]);
			if (Tri) steffenCoef(dtau_uw, dtau_dw, S[k - dk], S[k], S[k + dk], cSc);

		} else {
			dS_c   = (S[k] - S[k - dk]) / dtau_uw;
			cSc[0] = -1.0 / dtau_uw;
			cSc[1] = 1.0 / dtau_uw;
			cSc[2] = 0.0;
		}

		/* --- Compute interpolation parameters --         -------------- */

		dt03 = dtau_uw * 0.3333333333333333333;
		Bezier3_coeffs(dtau_uw, &alpha, &beta, &gamma, &theta, &eps);

		/* --- Source function control points --           -------------- */

		c1 = (S[k] - dt03 * dS_c);
		c2 = (S[k - dk] + dt03 * dS_up);

		/* --- Solve integral in this interval --          -------------- */

		I[k] = I_upw * eps + alpha * S[k - dk] + beta * S[k] + theta * c1 +
			   gamma * c2;

		/* --- Diagonal operator --                        -------------- */

		if (Psi) Psi[k] = Bezier3_psi(beta, theta, S[k], c1); /* B-66 */

		/* --- Row of the operator: I_k depends on S_k-dk, S_k, S_k+dk
			   directly, through c1 (dS_c) and c2 (dS_up, centred on
			   k - dk), and through I_upw on S_k-dk, S_k -- ------ */

		if (Tri) {
			double Tm_k, T0_k, Tp_k;

			Tm_k = eps * T0 + alpha - theta * dt03 * cSc[0] +
				   gamma * (1.0 + dt03 * cSup[1]);
			T0_k = eps * Tp + beta + theta * (1.0 - dt03 * cSc[1]) +
				   gamma * dt03 * cSup[2];
			Tp_k = -theta * dt03 * cSc[2];

			Tri[3 * k]     = Tm_k;
			Tri[3 * k + 1] = T0_k;
			Tri[3 * k + 2] = Tp_k;
			T0             = T0_k;
			Tp             = Tp_k;
			memcpy(cSup, cSc, sizeof(cSup));
		}

		/* --- Re-use downwind quantities for next upwind position -- --- */

		I_upw   = I[k];
		dsup    = dsdn;
		dchi_up = dchi_c;
		dchi_c  = dchi_dn;
		dtau_uw = dtau_dw;
		dS_up   = dS_c;
	}
}
/* ------- end ---------------------------- Piecewise_Bezier3_1D_tri.c */
