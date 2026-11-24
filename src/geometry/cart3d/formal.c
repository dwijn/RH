/* ------- file: -------------------------- formal.c ----------------

	   Version:       rh2.0, 3-D Cartesian, short characteristics
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Feb 21 16:24:24 2024 --

	   --------------------------                      ----------RH-- */

/* --- The transfer operations of the 3-D Cartesian geometry for the
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

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry   geometry;
extern Spectrum   spectrum;
extern InputData  input;
extern char       messageStr[];

/* ------- begin -------------------------- bgSupported.c ----------- */

static const char *bgSupported(void)
{
	/* --- The background scattering iteration (N-51) is not done in
		   2-D and 3-D (N-68) -- --                    -------------- */

	return "it is implemented in 1-D and spherical symmetry only (N-68)";
}
/* ------- end ---------------------------- bgSupported.c ----------- */

static bool_t bgLinearFallback(void)
{
	return FALSE;
}

/* ------- begin -------------------------- solve_3D.c -------------- */

/* --- Unpolarized short characteristics along ray mu; no tridiagonal
	   operator and no linear fallback (Tri and lin are never set) -- */

static void solve_3D(int nspect, int mu, bool_t to_obs, double *chi,
					 double *S, double *I, double *Psi, double *Tri,
					 bool_t lin)
{
	(void) Tri;
	(void) lin;
	ShortChar(&geometry, nspect, mu, to_obs, chi, S, I, Psi);
}
/* ------- end ---------------------------- solve_3D.c -------------- */

/* ------- begin -------------------------- solveStokes_3D.c -------- */

static void solveStokes_3D(int nspect, int mu, bool_t to_obs, double *chi,
						   double **S, double **I, double *Psi)
{
	ShortChar_Stokes(&geometry, nspect, mu, to_obs, chi, S, I, Psi);
}
/* ------- end ---------------------------- solveStokes_3D.c -------- */

/* ------- begin -------------------------- checkTau_3D.c ----------- */

static bool_t checkTau_3D(int nspect, int mu, int to_obs, const double *chi,
						  enum negtau_quadrature quad)
{
	return negTauColumns(nspect, mu, to_obs, geometry.Nz, geometry.Nplane, chi,
						 geometry.z, geometry.muz[mu], quad);
}
/* ------- end ---------------------------- checkTau_3D.c ----------- */

/* ------- begin -------------------------- emergent_3D.c ----------- */

/* --- The emergent intensity of ray mu at the top of each column -- */

static void emergent_3D(int nspect, int mu, const double *I, double **Ipol)
{
	register int l;

	int n = nspect * atmos.Nrays + mu;

	for (l = 0; l < geometry.Nplane; l++)
		spectrum.I[n][l] = I[l];
	if (Ipol) {
		for (l = 0; l < geometry.Nplane; l++) {
			spectrum.Stokes_Q[n][l] = Ipol[1][l];
			spectrum.Stokes_U[n][l] = Ipol[2][l];
			spectrum.Stokes_V[n][l] = Ipol[3][l];
		}
	}
}
/* ------- end ---------------------------- emergent_3D.c ----------- */

static double muz_3D(int mu)
{
	return geometry.muz[mu];
}
static double wmu_3D(int mu)
{
	return geometry.wmu[mu];
}

/* --- The operations of the 3-D geometry for the driver -- -------- */

const GeometryOps geometryOps = {
	.name               = "3-D Cartesian",
	.muz                = muz_3D,
	.wmu                = wmu_3D,
	.solve              = solve_3D,
	.solve_stokes       = solveStokes_3D,
	.check_tau          = checkTau_3D,
	.negtau_bezier      = NEGTAU_BEZIER_SC,
	.emergent           = emergent_3D,
	.bg_supported       = bgSupported,
	.bg_linear_fallback = bgLinearFallback,
	.bg_max_albedo      = NULL,
	.bg_add_operator    = NULL,
	.formal             = NULL};

/* ------- begin -------------------------- Hydrostatic.c ----------- */

void Hydrostatic(int NmaxIter, double iterLimit)
{
	const char routineName[] = "Hydrostatic";

	(void) NmaxIter;
	(void) iterLimit;

	if (atmos.hydrostatic) {
		sprintf(messageStr,
				"Can only establish hydrostatic equilibrium in 1-D geometry");
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
}
/* ------- end ---------------------------- Hydrostatic.c ----------- */
