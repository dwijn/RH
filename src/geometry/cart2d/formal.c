/* ------- file: -------------------------- formal.c ----------------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Sep 23 10:50:09 2025 --

	   --------------------------                      ----------RH-- */

/* --- The transfer operations of the 2-D Cartesian geometry for the
	   shared driver (Formal(), src/transfer/formal.c, Q-42): the short
	   characteristics, the negative optical depth check and the
	   emergent intensity --                          -------------- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "error.h"
#include "background.h"
#include "inputs.h"
#include "bgscatter.h"
#include "negtau.h"
#include "geometryops.h"
#include "constant.h"
#include "statistics.h"
#include "xdr.h"
#include "state.h"

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- bgSupported.c ----------- */

static const char *bgSupported(void)
{
	/* --- The background scattering iteration (N-51) is not done in
		   2-D and 3-D (N-68) -- --                    -------------- */

	return "it is implemented in 1-D only (N-68)";
}
/* ------- end ---------------------------- bgSupported.c ----------- */

static bool_t bgLinearFallback(void)
{
	return FALSE;
}

/* ------- begin -------------------------- solve_2D.c -------------- */

/* --- Unpolarized short characteristics along ray mu; no tridiagonal
	   operator and no linear fallback (Tri and lin are never set) -- */

static void solve_2D(int nspect, int mu, bool_t to_obs, double *chi,
					 double *S, double *I, double *Psi, double *Tri,
					 bool_t lin)
{
	const char routineName[] = "solve_2D";

	(void) Tri;
	(void) lin;
	if (rh->input.S_interpolation == S_LINEAR) {
		Piecewise_Linear_2D(rh->geometry, nspect, mu, to_obs, chi, S, I, Psi);
	} else if (rh->input.S_interpolation == S_BEZIER3) {
		Piecewise_Bezier3_2D(rh->geometry, nspect, mu, to_obs, chi, S, I, Psi);
	} else {
		Error(rh, ERROR_LEVEL_1, routineName, "Unknown radiation solver: %d",
			  rh->input.S_interpolation);
	}
}
/* ------- end ---------------------------- solve_2D.c -------------- */

/* ------- begin -------------------------- solveStokes_2D.c -------- */

static void solveStokes_2D(int nspect, int mu, bool_t to_obs, double *chi,
						   double **S, double **I, double *Psi)
{
	const char routineName[] = "solveStokes_2D";

	if (rh->input.S_interpolation_stokes == DELO_BEZIER3) {
		Piece_Stokes_Bezier3_2D(rh->geometry, nspect, mu, to_obs, chi, S, I, Psi);
	} else {
		Error(rh, ERROR_LEVEL_1, routineName, "Unknown polarization solver: %d",
			  rh->input.S_interpolation_stokes);
	}
}
/* ------- end ---------------------------- solveStokes_2D.c -------- */

/* ------- begin -------------------------- checkTau_2D.c ----------- */

static bool_t checkTau_2D(int nspect, int mu, int to_obs, const double *chi,
						  enum negtau_quadrature quad)
{
	return negTauColumns(nspect, mu, to_obs, rh->geometry->Nz, rh->geometry->Nx, chi,
						 rh->geometry->z, rh->geometry->muz[mu], quad);
}
/* ------- end ---------------------------- checkTau_2D.c ----------- */

/* ------- begin -------------------------- emergent_2D.c ----------- */

/* --- The emergent intensity of ray mu at the top of each column -- */

static void emergent_2D(int nspect, int mu, const double *I, double **Ipol)
{
	register int l;

	int n = nspect * rh->atmos.Nrays + mu;

	for (l = 0; l < rh->geometry->Nx; l++)
		rh->spectrum.I[n][l] = I[l];
	if (Ipol) {
		for (l = 0; l < rh->geometry->Nx; l++) {
			rh->spectrum.Stokes_Q[n][l] = Ipol[1][l];
			rh->spectrum.Stokes_U[n][l] = Ipol[2][l];
			rh->spectrum.Stokes_V[n][l] = Ipol[3][l];
		}
	}
}
/* ------- end ---------------------------- emergent_2D.c ----------- */

static double muz_2D(int mu)
{
	return rh->geometry->muz[mu];
}
static double wmu_2D(int mu)
{
	return rh->geometry->wmu[mu];
}

/* --- The operations of the 2-D geometry for the driver -- -------- */

const GeometryOps geometryOps = {
	.name               = "2-D Cartesian",
	.muz                = muz_2D,
	.wmu                = wmu_2D,
	.solve              = solve_2D,
	.solve_stokes       = solveStokes_2D,
	.check_tau          = checkTau_2D,
	.negtau_bezier      = NEGTAU_BEZIER_SC,
	.emergent           = emergent_2D,
	.bg_supported       = bgSupported,
	.bg_linear_fallback = bgLinearFallback,
	.bg_max_albedo      = NULL,
	.bg_add_operator    = NULL};

/* ------- begin -------------------------- Hydrostatic.c ----------- */

void Hydrostatic(int NmaxIter, double iterLimit)
{
	const char routineName[] = "Hydrostatic";

	(void) NmaxIter;
	(void) iterLimit;

	if (rh->atmos.hydrostatic) {
		Error(rh, ERROR_LEVEL_2, routineName, "Can only establish hydrostatic equilibrium in 1-D geometry");
	}
}
/* ------- end ---------------------------- Hydrostatic.c ----------- */
