/* ------- file: -------------------------- piecewise_1D.c ----------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu May 31 13:26:58 2018 --

	   --------------------------                      ----------RH-- */

/* --- Piecewise integration of transfer equation in one dimension.

	-- For boundary condition THERMALIZED use the relation

		  I ~= B - \mu dB / d\tau,

	   where \tau is taken in the direction of the ray (i.e. NOT in the
	   sense of optical depth, e.g. see Mihalas, 1978, p. 51).
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern char messageStr[];

/* ------- begin -------------------------- Piecewise_Linear_1D.c --- */

void Piecewise_Linear_1D(int nspect, int mu, bool_t to_obs,
						 double *chi, double *S, double *I, double *Psi)
{
	Piecewise_Linear_1D_tri(nspect, mu, to_obs, chi, S, I, Psi, NULL);
}
/* ------- end ---------------------------- Piecewise_Linear_1D.c --- */

/* ------- begin -------------------------- Piecewise_Linear_1D_tri.c */

void Piecewise_Linear_1D_tri(int nspect, int mu, bool_t to_obs,
							 double *chi, double *S, double *I, double *Psi,
							 double *Tri)
{
	/* --- Tri (N-51): the row of the operator of this ray,
		   Tri[3k] = dI_k/dS_k-dk (upwind), Tri[3k + 1] = dI_k/dS_k,
		   Tri[3k + 2] = 0 (no downwind dependence) -- -------------- */

	const char   routineName[] = "Piecewise_Linear_1D";
	register int k;

	int    k_start, k_end, dk, Ndep = rh->geometry->Ndep;
	double dtau_uw, dS_uw, I_uw = 0.0, w[2], zmu, Bnu[2], T0 = 0.0;

	zmu = 0.5 / rh->geometry->muz[mu];

	/* --- Distinguish between rays going from BOTTOM to TOP
		   (to_obs == TRUE), and vice versa --      -------------- */

	if (to_obs) {
		dk      = -1;
		k_start = Ndep - 1;
		k_end   = 0;
	} else {
		dk      = 1;
		k_start = 0;
		k_end   = Ndep - 1;
	}
	dtau_uw = zmu * (chi[k_start] + chi[k_start + dk]) *
			  fabs(rh->geometry->height[k_start] - rh->geometry->height[k_start + dk]);
	dS_uw = (S[k_start] - S[k_start + dk]) / dtau_uw;

	/* --- Boundary conditions --                        -------------- */

	if (to_obs) {
		switch (rh->geometry->vboundary[BOTTOM]) {
		case ZERO:
			I_uw = 0.0;
			break;
		case THERMALIZED:
			Planck(2, &rh->atmos.T[Ndep - 2], rh->spectrum.lambda[nspect], Bnu);
			I_uw = Bnu[1] - (Bnu[0] - Bnu[1]) / dtau_uw;
			break;
		case IRRADIATED:
			I_uw = rh->geometry->Ibottom[nspect][mu];
			break;
		case REFLECTIVE:
			sprintf(messageStr, "Boundary condition not implemented: %d",
					rh->geometry->vboundary[BOTTOM]);
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
	} else {
		switch (rh->geometry->vboundary[TOP]) {
		case ZERO:
			I_uw = 0.0;
			break;
		case THERMALIZED:
			Planck(2, &rh->atmos.T[0], rh->spectrum.lambda[nspect], Bnu);
			I_uw = Bnu[0] - (Bnu[1] - Bnu[0]) / dtau_uw;
			break;
		case IRRADIATED:
			I_uw = rh->geometry->Itop[nspect][mu];
			break;
		case REFLECTIVE:
			sprintf(messageStr, "Boundary condition not implemented: %d",
					rh->geometry->vboundary[BOTTOM]);
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
	}

	I[k_start] = I_uw;
	if (Psi) Psi[k_start] = 0.0;
	if (Tri) Tri[3 * k_start] = Tri[3 * k_start + 1] = Tri[3 * k_start + 2] = 0.0;

	/* --- Solve transfer along ray --                   -------------- */

	for (k = k_start + dk; k != k_end; k += dk) {
		w2(dtau_uw, w);

		I[k] = (1.0 - w[0]) * I_uw + w[0] * S[k] + w[1] * dS_uw;
		if (Psi) Psi[k] = w[0] - w[1] / dtau_uw;
		if (Tri) {
			Tri[3 * k]     = (1.0 - w[0]) * T0 + w[1] / dtau_uw;
			Tri[3 * k + 1] = T0 = w[0] - w[1] / dtau_uw;
			Tri[3 * k + 2]      = 0.0;
		}

		dtau_uw = zmu * (chi[k] + chi[k + dk]) *
				  fabs(rh->geometry->height[k] - rh->geometry->height[k + dk]);
		dS_uw = (S[k] - S[k + dk]) / dtau_uw;

		I_uw = I[k];
	}
	w2(dtau_uw, w);
	I[k_end] = (1.0 - w[0]) * I_uw + w[0] * S[k_end] + w[1] * dS_uw;
	if (Psi) Psi[k_end] = w[0] - w[1] / dtau_uw;
	if (Tri) {
		Tri[3 * k_end]     = (1.0 - w[0]) * T0 + w[1] / dtau_uw;
		Tri[3 * k_end + 1] = w[0] - w[1] / dtau_uw;
		Tri[3 * k_end + 2] = 0.0;
	}
}
/* ------- end ---------------------------- Piecewise_Linear_1D_tri.c */
