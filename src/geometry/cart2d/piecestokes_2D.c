/* ------- file: -------------------------- piecestokes_2D.c --------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Sep 23 12:23:18 2025 --

	   --------------------------                      ----------RH-- */

/* --- Formal 2-D polarized radiative transfer solver using
	   the method of short characteristics. The transfer equation is
	   solved in one direction along the ray only. Use keyword
	   to_observer=TRUE to solve in the upward direction, and
	   to_observer=FALSE to solve in opposite direction.

  See: P. B. Kunasz and  L. H. Auer 1988, JQSRT 39, 67-79

	   For now it implements fixed and periodic boundary
	   conditions in the horizontal direction.

  See: L. H. Auer, P. Fabiani Bendicho and J. Trujillo Bueno 1994,
	   A&A, 292, 599-615

	   At the top and bottom only a fixed boundary is allowed.

  See: L. H. Auer 2003, Formal Solution: EXPLICIT Answers, in
	   Stellar Atmosphere Modeling, ASP Conference Proceedings, Vol. 288,
	   eds. Hubeny, Mihalas & Werner.


 Note: Rays with to_observer == TRUE go from BOTTOM to TOP.

	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "inputs.h"
#include "error.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- StokesK_2D_convex.c ----- */

/* --- Where the interpolation is linear (SC_LINEAR: LINEAR_2D, and
	   the stencils that have no second neighbour), the weights are
	   convex and the same for every quantity: each component is
	   interpolated on its own as before, inside the cone (eta' a
	   convex combination of the eta' of the nodes), and a linear field
	   is reproduced exactly. The Steffen-limited Bezier interpolation,
	   whose weights depend on each component, goes through the
	   rapidity (stokesgrid.h, B-63) --                -------------- */

static void StokesK_2D_convex(RHState *rh, StokesGrid *sg, Stencil *st, enum direction direction,
							  double chi_I, double K[4][4])
{
	register int i, j;

	int        nspect = sg->nspect;
	ActiveSet *as;

	/* --- Return the elements of the reduced 4x4 Stokes
		   opacity matrix K', which is defined as:

			 =           =         =              =
			 K' = (chi_c*1 + chi_l*Phi) / chi_I - 1,

		 for wavelength# nspect, spatial point k, and ray mu.

		   This is the 2-D short characteristics version.

	  See: Rees, Murphy, & Durrant, 1989, ApJ 339, 1093-1106.

		   --                                            -------------- */
	as = &rh->spectrum.as[nspect];

	for (j = 0; j < 4; j++)
		for (i = 0; i < 4; i++) K[j][i] = 0.0;

	/* --- First the contribution from the active set -- -------------- */

	if (containsPolarized(rh, as)) {
		K[0][1] = Interp_2D(as->chi + rh->atmos.Nspace,
							st, direction);
		K[0][2] = Interp_2D(as->chi + 2 * rh->atmos.Nspace,
							st, direction);
		K[0][3] = Interp_2D(as->chi + 3 * rh->atmos.Nspace,
							st, direction);

		if (rh->input.magneto_optical) {
			K[1][2] = Interp_2D(as->chip + 2 * rh->atmos.Nspace,
								st, direction);
			K[1][3] = Interp_2D(as->chip + rh->atmos.Nspace,
								st, direction);
			K[2][3] = Interp_2D(as->chip,
								st, direction);
		}
	}
	/* --- Add possible contribution from the background -- ----------- */

	if (rh->atmos.backgrflags[nspect].ispolarized) {
		K[0][1] += Interp_2D(as->chi_c + rh->atmos.Nspace,
							 st, direction);
		K[0][2] += Interp_2D(as->chi_c + 2 * rh->atmos.Nspace,
							 st, direction);
		K[0][3] += Interp_2D(as->chi_c + 3 * rh->atmos.Nspace,
							 st, direction);

		if (rh->input.magneto_optical) {
			K[1][2] += Interp_2D(as->chip_c + 2 * rh->atmos.Nspace,
								 st, direction);
			K[1][3] += Interp_2D(as->chip_c + rh->atmos.Nspace,
								 st, direction);
			K[2][3] += Interp_2D(as->chip_c,
								 st, direction);
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
/* ------- end ---------------------------- StokesK_2D_convex.c ----- */

/* ------- begin -------------------------- StokesK_2D.c ------------ */

void StokesK_2D(RHState *rh, StokesGrid *sg, Stencil *st, enum direction direction,
				double chi_I, double K[4][4])
{
	/* --- The reduced 4x4 Stokes propagation matrix K',

			 K' = (chi_c*1 + chi_l*Phi) / chi_I - 1,

		   (Rees, Murphy, & Durrant, 1989, ApJ 339, 1093-1106) at the
		   up- or downwind point of stencil st, for the wavelength and
		   ray of the grid sg; chi_I is the interpolated opacity there.
		   With linear interpolation from the components
		   (StokesK_2D_convex); otherwise eta' = (chi_Q, chi_U, chi_V) /
		   chi_I from the rapidity of the grid, |eta'| < 1 by
		   construction, and rho' with the same weights (stokesgrid.h,
		   B-63). This is the 2-D short characteristics version -- --- */

	double w[3], v[4];

	if (st->order[direction] == SC_LINEAR) {
		StokesK_2D_convex(rh, sg, st, direction, chi_I, K);
		return;
	}
	for (int j = 0; j < 4; j++)
		for (int i = 0; i < 4; i++) K[j][i] = 0.0;
	if (!sg->polarizedK) return;

	for (int n = 0; n < 3; n++) w[n] = Interp_2D_linear(sg->wK[n], st, direction);
	Stokes_from_rapidity(1.0, w, v);
	for (int n = 1; n < 4; n++) K[0][n] = K[n][0] = v[n];

	if (sg->magneto_optical) {
		K[1][2] = Interp_2D_linear(sg->rK[0], st, direction);
		K[3][1] = Interp_2D_linear(sg->rK[1], st, direction);
		K[2][3] = Interp_2D_linear(sg->rK[2], st, direction);
		K[2][1] = -K[1][2];
		K[1][3] = -K[3][1];
		K[3][2] = -K[2][3];
	}
}
/* ------- end ---------------------------- StokesK_2D.c ------------ */

/* ------- begin -------------------------- Interp_Stokes_2D.c ------ */

void Interp_Stokes_2D(StokesGrid *sg, enum stokes_vector v, Stencil *st,
					  enum direction direction, double out[4])
{
	/* --- The source vector (v = STOKES_SOURCE) or the Stokes vector
		   (STOKES_INTENSITY) of the grid at the up- or downwind point
		   of stencil st: with linear interpolation each component,
		   otherwise the magnitude with Interp_2D and the polarization
		   from the rapidity of the grid with Interp_2D_linear: inside
		   the cone by construction either way (stokesgrid.h, B-63)
		   --                                              -------------- */

	double **f  = (v == STOKES_SOURCE) ? sg->S : sg->I,
		   **wg = (v == STOKES_SOURCE) ? sg->wS : sg->wI, w[3];

	if (st->order[direction] == SC_LINEAR) {
		for (int n = 0; n < 4; n++) out[n] = Interp_2D(f[n], st, direction);
		return;
	}
	for (int n = 0; n < 3; n++) w[n] = Interp_2D_linear(wg[n], st, direction);
	Stokes_from_rapidity(Interp_2D(f[0], st, direction), w, out);
}
/* ------- end ---------------------------- Interp_Stokes_2D.c ------ */
