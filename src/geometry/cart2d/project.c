/* ------- file: -------------------------- vproject.c --------------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Feb 12 13:32:43 2003 --

	   --------------------------                      ----------RH-- */

/* --- Calculate projected velocity along ray mu at depth k -- ------ */

#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- vproject.c -------------- */

double vproject(RHState *rh, int k, int mu)
{
	return rh->geometry->mux[mu] * rh->geometry->vx[k] +
		   rh->geometry->muz[mu] * rh->geometry->vz[k];
}
/* ------- end ---------------------------- vproject.c -------------- */

/* ------- begin -------------------------- Bproject.c -------------- */

void Bproject(RHState *rh)
{
	/* --- Compute the cosine of gamma, the angle of B with the line of
		   sight, and the sine and cosine of 2*chi, the angle of B
		   with e1 --                                    -------------- */

	register int k, mu;

	double bx, by, bz, b1, b2, b3, csc_theta, sin_gamma;

	rh->atmos.cos_gamma = matrix_double(rh, rh->atmos.Nrays, rh->atmos.Nspace);
	rh->atmos.cos_2chi  = matrix_double(rh, rh->atmos.Nrays, rh->atmos.Nspace);
	rh->atmos.sin_2chi  = matrix_double(rh, rh->atmos.Nrays, rh->atmos.Nspace);

	for (mu = 0; mu < rh->atmos.Nrays; mu++) {
		if (rh->geometry->muz[mu] == 1.0) {
			for (k = 0; k < rh->atmos.Nspace; k++) {
				rh->atmos.cos_gamma[mu][k] = cos(rh->atmos.gamma_B[k]);
				rh->atmos.cos_2chi[mu][k]  = cos(2.0 * rh->atmos.chi_B[k]);
				rh->atmos.sin_2chi[mu][k]  = sin(2.0 * rh->atmos.chi_B[k]);
			}
		} else {
			csc_theta = 1.0 / sqrt(1.0 - SQ(rh->geometry->muz[mu]));
			for (k = 0; k < rh->atmos.Nspace; k++) {
				sin_gamma = sin(rh->atmos.gamma_B[k]);
				bx        = sin_gamma * cos(rh->atmos.chi_B[k]);
				by        = sin_gamma * sin(rh->atmos.chi_B[k]);
				bz        = cos(rh->atmos.gamma_B[k]);

				b3 = rh->geometry->mux[mu] * bx + rh->geometry->muy[mu] * by +
					 rh->geometry->muz[mu] * bz;
				b1 = csc_theta * (bz - rh->geometry->muz[mu] * b3);
				b2 = csc_theta * (rh->geometry->muy[mu] * bx - rh->geometry->mux[mu] * by);

				rh->atmos.cos_gamma[mu][k] = b3;
				/* --- No transverse field (B along the ray): the angle chi is
					   undefined, and the terms it enters carry a factor
					   sin^2(gamma) = 0 (B-28) --         -------------- */

				if (1.0 - SQ(b3) > 0.0) {
					rh->atmos.cos_2chi[mu][k] = (SQ(b1) - SQ(b2)) / (1.0 - SQ(b3));
					rh->atmos.sin_2chi[mu][k] = 2.0 * b1 * b2 / (1.0 - SQ(b3));
				} else {
					rh->atmos.cos_2chi[mu][k] = 1.0;
					rh->atmos.sin_2chi[mu][k] = 0.0;
				}
			}
		}
	}
}
/* ------- end ---------------------------- Bproject.c -------------- */
