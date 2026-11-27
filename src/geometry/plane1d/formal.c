/* ------- file: -------------------------- formal.c ----------------

	   Version:       rh1.0, 1-D plane-parallel
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Feb 20 15:18:13 2024 --

	   --------------------------                      ----------RH-- */

/* --- The transfer operations of the 1-D plane-parallel geometry for
	   the shared driver (Formal(), src/transfer/formal.c, Q-42): the
	   ray solvers, the negative optical depth check, the emergent
	   intensity, and the background scattering iteration -- ------ */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "constant.h"
#include "background.h"
#include "inputs.h"
#include "error.h"
#include "xdr.h"
#include "bezier.h"
#include "negtau.h"
#include "bgscatter.h"
#include "geometryops.h"
#include "state.h"

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- Piecewise_1D.c ---------- */

/* --- Unpolarized formal solution along one ray with the method chosen
	   by S_INTERPOLATION, or the linear one where lin is set (N-52) - */

static void Piecewise_1D(int nspect, int mu, bool_t to_obs, double *chi,
						 double *S, double *I, double *Psi, double *Tri,
						 bool_t lin)
{
	const char routineName[] = "Piecewise_1D";

	/* --- Tri: NULL, or the row of the operator of the ray (N-51) -- */

	if (lin || rh->input.S_interpolation == S_LINEAR) {
		Piecewise_Linear_1D_tri(nspect, mu, to_obs, chi, S, I, Psi, Tri);
	} else if (rh->input.S_interpolation == S_BEZIER3) {
		Piecewise_Bezier3_1D_tri(nspect, mu, to_obs, chi, S, I, Psi, Tri);
	} else {
		Error(ERROR_LEVEL_1, routineName, "Unknown radiation solver: %d",
			  rh->input.S_interpolation);
	}
}
/* ------- end ---------------------------- Piecewise_1D.c ---------- */

/* ------- begin -------------------------- bgSupported.c ----------- */

static const char *bgSupported(void)
{
	/* --- The background scattering iteration (N-51) needs a positive
		   operator: the cubic Bezier solver's is not (N-52), so with
		   S_BEZIER3 the iterated wavelengths are solved with the linear
		   solver, unless BACKGROUND_LINEAR_ALBEDO >= 1 turns that off - */

	if (rh->input.S_interpolation != S_LINEAR && rh->input.bg_linearAlbedo >= 1.0)
		return "with S_INTERPOLATION = S_BEZIER3 it needs"
			   " BACKGROUND_LINEAR_ALBEDO < 1 (the operator of S_BEZIER3 is"
			   " not positive, N-52)";
	return NULL;
}
/* ------- end ---------------------------- bgSupported.c ----------- */

/* ------- begin -------------------------- bgLinearFallback.c ------ */

static bool_t bgLinearFallback(void)
{
	return rh->input.S_interpolation != S_LINEAR;
}
/* ------- end ---------------------------- bgLinearFallback.c ------ */

/* ------- begin -------------------------- maxAlbedo.c ------------- */

static double maxAlbedo(const ActiveSet *as)
{
	register int k;
	double       a = 0.0, tau = 0.0, c, cprev = 0.0;

	/* --- Largest background scattering albedo sca_c / chi below the
		   vertical optical depth BACKGROUND_LINEAR_TAU (N-52): above it
		   the photons escape and the lambda iteration converges -- -- */

	for (k = 0; k < rh->atmos.Nspace; k++) {
		c = as->chi[k] + as->chi_c[k];
		if (k > 0)
			tau += 0.5 * (c + cprev) *
				   fabs(rh->geometry->height[k] - rh->geometry->height[k - 1]);
		cprev = c;
		if (tau >= rh->input.bg_linearTau && c > 0.0 && as->sca_c[k] / c > a)
			a = as->sca_c[k] / c;
	}
	return a;
}
/* ------- end ---------------------------- maxAlbedo.c ------------- */

/* ------- begin -------------------------- addTri.c ---------------- */

static void addTri(int Nspace, bool_t to_obs, double wmu, const double *Tri,
				   const double *sca_c, const double *chi, double *A)
{
	register int k;

	double dm, dp;

	/* --- Add the row of the operator of one ray (Tri: upwind, local,
		   downwind) to A[3k + i] = dJ_k / du_k-1+i, with S_j =
		   ... + sca_c,j u_j / chi_j (N-51). Upwind is k + 1 for rays to
		   the observer, k - 1 for the others -- --       -------------- */

	for (k = 0; k < Nspace; k++) {
		dm = (to_obs) ? Tri[3 * k + 2] : Tri[3 * k];
		dp = (to_obs) ? Tri[3 * k] : Tri[3 * k + 2];
		if (k > 0) A[3 * k] += wmu * dm * sca_c[k - 1] / chi[k - 1];
		A[3 * k + 1] += wmu * Tri[3 * k + 1] * sca_c[k] / chi[k];
		if (k < Nspace - 1) A[3 * k + 2] += wmu * dp * sca_c[k + 1] / chi[k + 1];
	}
}
/* ------- end ---------------------------- addTri.c ---------------- */

/* ------- begin -------------------------- solveStokes_1D.c -------- */

static void solveStokes_1D(int nspect, int mu, bool_t to_obs, double *chi,
						   double **S, double **I, double *Psi)
{
	const char routineName[] = "solveStokes_1D";

	if (rh->input.S_interpolation_stokes == DELO_BEZIER3) {
		Piece_Stokes_Bezier3_1D(nspect, mu, to_obs, chi, S, I, Psi);
	} else {
		Error(ERROR_LEVEL_1, routineName, "Unknown polarization solver: %d",
			  rh->input.S_interpolation_stokes);
	}
}
/* ------- end ---------------------------- solveStokes_1D.c -------- */

/* ------- begin -------------------------- checkTau_1D.c ----------- */

static bool_t checkTau_1D(int nspect, int mu, int to_obs, const double *chi,
						  enum negtau_quadrature quad)
{
	return negTauColumns(nspect, mu, to_obs, rh->atmos.Nspace, 1, chi,
						 rh->geometry->height, rh->geometry->muz[mu], quad);
}
/* ------- end ---------------------------- checkTau_1D.c ----------- */

/* ------- begin -------------------------- emergent_1D.c ----------- */

static void emergent_1D(int nspect, int mu, const double *I, double **Ipol)
{
	rh->spectrum.I[nspect][mu] = I[0];
	if (Ipol) {
		rh->spectrum.Stokes_Q[nspect][mu] = Ipol[1][0];
		rh->spectrum.Stokes_U[nspect][mu] = Ipol[2][0];
		rh->spectrum.Stokes_V[nspect][mu] = Ipol[3][0];
	}
}
/* ------- end ---------------------------- emergent_1D.c ----------- */

static double muz_1D(int mu)
{
	return rh->geometry->muz[mu];
}
static double wmu_1D(int mu)
{
	return rh->geometry->wmu[mu];
}

/* --- The operations of the 1-D geometry for the driver -- -------- */

const GeometryOps geometryOps = {
	.name               = "1-D plane",
	.muz                = muz_1D,
	.wmu                = wmu_1D,
	.solve              = Piecewise_1D,
	.solve_stokes       = solveStokes_1D,
	.check_tau          = checkTau_1D,
	.negtau_bezier      = NEGTAU_BEZIER3_1D,
	.emergent           = emergent_1D,
	.bg_supported       = bgSupported,
	.bg_linear_fallback = bgLinearFallback,
	.bg_max_albedo      = maxAlbedo,
	.bg_add_operator    = addTri};
