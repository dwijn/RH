/* ------- file: -------------------------- formal_1D.c --------------

	   Version:       rh1.0, 1-D plane-parallel
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Feb 20 15:18:13 2024 --

	   --------------------------                      ----------RH-- */

/* --- Formal solution with given source function, and allowing for
	   a PRD emission profile, polarized line radiation, polarized
	   background lines, and scattering background polarization -- -- */

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

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry   geometry;
extern Spectrum   spectrum;
extern InputData  input;
extern char       messageStr[];

/* ------- begin -------------------------- Piecewise_1D.c ---------- */

/* --- Unpolarized formal solution along one ray with the method chosen
	   by S_INTERPOLATION, or the linear one where lin is set (N-52) - */

static void Piecewise_1D(int nspect, int mu, bool_t to_obs, double *chi,
						 double *S, double *I, double *Psi, double *Tri,
						 bool_t lin)
{
	const char routineName[] = "Piecewise_1D";

	/* --- Tri: NULL, or the row of the operator of the ray (N-51) -- */

	if (lin || input.S_interpolation == S_LINEAR) {
		Piecewise_Linear_1D_tri(nspect, mu, to_obs, chi, S, I, Psi, Tri);
	} else if (input.S_interpolation == S_BEZIER3) {
		Piecewise_Bezier3_1D_tri(nspect, mu, to_obs, chi, S, I, Psi, Tri);
	} else {
		sprintf(messageStr, "Unknown radiation solver: %d",
				input.S_interpolation);
		Error(ERROR_LEVEL_1, routineName, messageStr);
	}
}
/* ------- end ---------------------------- Piecewise_1D.c ---------- */

/* ------- begin -------------------------- bgScatterSupported.c ---- */

const char *bgScatterSupported(void)
{
	/* --- The background scattering iteration (N-51) needs a positive
		   operator: the cubic Bezier solver's is not (N-52), so with
		   S_BEZIER3 the iterated wavelengths are solved with the linear
		   solver, unless BACKGROUND_LINEAR_ALBEDO >= 1 turns that off - */

	if (input.S_interpolation != S_LINEAR && input.bg_linearAlbedo >= 1.0)
		return "with S_INTERPOLATION = S_BEZIER3 it needs"
			   " BACKGROUND_LINEAR_ALBEDO < 1 (the operator of S_BEZIER3 is"
			   " not positive, N-52)";
	return NULL;
}
/* ------- end ---------------------------- bgScatterSupported.c ---- */

/* ------- begin -------------------------- maxAlbedo.c ------------- */

static double maxAlbedo(const ActiveSet *as)
{
	register int k;
	double       a = 0.0, tau = 0.0, c, cprev = 0.0;

	/* --- Largest background scattering albedo sca_c / chi below the
		   vertical optical depth BACKGROUND_LINEAR_TAU (N-52): above it
		   the photons escape and the lambda iteration converges -- -- */

	for (k = 0; k < atmos.Nspace; k++) {
		c = as->chi[k] + as->chi_c[k];
		if (k > 0)
			tau += 0.5 * (c + cprev) *
				   fabs(geometry.height[k] - geometry.height[k - 1]);
		cprev = c;
		if (tau >= input.bg_linearTau && c > 0.0 && as->sca_c[k] / c > a)
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

/* ------- begin -------------------------- bgSweep_1D.c ------------ */

/* --- One formal solution of all rays at nspect for the background
	   scattering iteration (N-51), with the background scattering
	   source sca_c u and the populations of this sweep: J, the
	   tridiagonal A of dJ/du and, if not NULL, the angle average of
	   sca_c / chi; no rates, no output. With angle_dep the opacities
	   are those of each ray again (readBackground, Opacity: both only
	   assign); otherwise those already in as and chi -- -------- */

typedef struct {
	int        nspect;
	bool_t     angle_dep, boundbound;
	double    *chi, *S, *I, *Tri;
	ActiveSet *as;
} BgCtx1D;

static void bgSweep_1D(void *ctx, const double *u, double *J, double *A,
					   double *albedo)
{
	register int k, mu;

	BgCtx1D   *c  = (BgCtx1D *) ctx;
	ActiveSet *as = c->as;
	bool_t     to_obs, initialize;
	int        Nspace = atmos.Nspace;
	double     wmu;

	for (k = 0; k < Nspace; k++) J[k] = 0.0;
	for (k = 0; k < 3 * Nspace; k++) A[k] = 0.0;
	if (albedo)
		for (k = 0; k < Nspace; k++) albedo[k] = 0.0;
	for (mu = 0; mu < atmos.Nrays; mu++) {
		wmu = 0.5 * geometry.wmu[mu];
		for (to_obs = 0; to_obs <= 1; to_obs++) {
			if (c->angle_dep) {
				initialize = (mu == 0 && to_obs == 0);
				if (initialize || atmos.backgrflags[c->nspect].hasline)
					readBackground(c->nspect, mu, to_obs);
				if (initialize || c->boundbound)
					Opacity(c->nspect, mu, to_obs, initialize);
				for (k = 0; k < Nspace; k++) c->chi[k] = as->chi[k] + as->chi_c[k];
			}
			if (c->angle_dep || (mu == 0 && to_obs == 0)) {
				for (k = 0; k < Nspace; k++)
					c->S[k] = (as->eta[k] + as->eta_c[k] + as->sca_c[k] * u[k]) /
							  c->chi[k];
			}
			Piecewise_1D(c->nspect, mu, to_obs, c->chi, c->S, c->I, NULL, c->Tri,
						 TRUE);
			for (k = 0; k < Nspace; k++) J[k] += wmu * c->I[k];
			addTri(Nspace, to_obs, wmu, c->Tri, as->sca_c, c->chi, A);
			if (albedo)
				for (k = 0; k < Nspace; k++)
					albedo[k] += wmu * as->sca_c[k] / c->chi[k];
		}
	}
}
/* ------- end ---------------------------- bgSweep_1D.c ------------ */

/* ------- begin -------------------------- linearChoice.c ---------- */

/* --- With S_BEZIER3, whether this wavelength is solved with the
	   linear solver (N-52): where its background scattering albedo
	   (maxAlbedo) exceeds BACKGROUND_LINEAR_ALBEDO, decided in every
	   sweep that iterates the background scattering; other formal
	   solutions keep the last choice (bgScatterLinear). Only those
	   wavelengths iterate the background scattering (*bg, *Tri) -- */

static bool_t linearChoice(int nspect, bool_t *bg, double **Tri,
						   const ActiveSet *as)
{
	bool_t lin;

	if (input.S_interpolation == S_LINEAR) return FALSE;
	lin = bgScatterLinear(nspect, (*bg) ? maxAlbedo(as) : -1.0);
	if (!lin) {
		*bg  = FALSE;
		*Tri = NULL;
	}
	return lin;
}
/* ------- end ---------------------------- linearChoice.c ---------- */

/* ------- begin -------------------------- Formal.c ---------------- */

double Formal(int nspect, bool_t eval_operator, bool_t redistribute)
{
	const char   routineName[] = "Formal";
	register int k, mu, n;

	bool_t initialize, boundbound, polarized_as, polarized_c,
		PRD_angle_dep, to_obs, solveStokes, angle_dep;
	int                    Nspace = atmos.Nspace, Nrays = atmos.Nrays, muslant;
	enum negtau_quadrature quad, quad_ray;
	double                *I, *chi, *S, **Ipol = NULL, **Spol = NULL, *Psi, *Jdag, wmu,
						 dJmax = 0.0, dJ, *J20dag = NULL, musq, threemu1 = 0.0, threemu2 = 0.0,
						 *J, *J20 = NULL, *Abg = NULL, *albedo = NULL, *Tri = NULL,
						 *Jsrc, *Jstar = NULL;
	bool_t     bg, pending = FALSE, lin = FALSE;
	BgCtx1D    bgc;
	ActiveSet *as;

	/* --- Retrieve active set as of transitions at wavelength nspect - */

	as = &spectrum.as[nspect];
	alloc_as(nspect, eval_operator);

	/* --- Check whether current active set includes a bound-bound
		   and/or polarized transition and/or angle-dependent PRD
		   transition, and/or polarization through background scattering.
		   Otherwise, only angle-independent opacity and source functions
		   are needed --                                 -------------- */

	/* --- Check for bound-bound transition in active set -- ---------- */

	boundbound = containsBoundBound(as);

	/* --- Check for line with angle-dependent PRD in set -- ---------- */

	PRD_angle_dep = (containsPRDline(as) &&
					 (input.PRD_angle_dep == PRD_ANGLE_DEP ||
					  input.PRD_angle_dep == PRD_ANGLE_APPROX));

	/* --- Check for polarized bound-bound transition in active set - - */

	polarized_as = containsPolarized(as);

	/* --- Check for polarized bound-bound transition in background - - */

	polarized_c = atmos.backgrflags[nspect].ispolarized;

	/* --- Determine if we solve for I, or for I, Q, U, V -- ---------- */

	solveStokes = (input.StokesMode == FULL_STOKES &&
				   (polarized_as || polarized_c || input.backgr_pol));

	/* --- The negative-opacity detector integrates chi as the formal
		   solver will (B-68): Bezier for S_BEZIER3 and for the
		   polarized solver (DELO_BEZIER3), trapezoid for S_LINEAR -- */

	quad     = (input.S_interpolation == S_BEZIER3) ? NEGTAU_BEZIER3_1D : NEGTAU_LINEAR;
	quad_ray = (solveStokes) ? NEGTAU_BEZIER3_1D : quad;

	/* --- Determine if we have to do angle-dependent opacity and
		   emissivity --                                 -------------- */

	angle_dep = (polarized_as || polarized_c || PRD_angle_dep || input.backgr_pol ||
				 (atmos.moving &&
				  (boundbound || atmos.backgrflags[nspect].hasline)));

	/* --- Iterate the background scattering at this wavelength after
		   the main sweep (N-51): needs the diagonal operator; not with
		   Stokes or background polarization (J20), which keep the
		   lambda iteration -- --                       -------------- */

	bg = (spectrum.updateJ && bgScatterActive() && !solveStokes &&
		  !input.backgr_pol);

	/* --- Allocate temporary space --                   -------------- */

	if (bg) {
		Abg    = (double *) calloc(7 * Nspace, sizeof(double));
		albedo = Abg + 3 * Nspace;
		Tri    = albedo + Nspace;
	}
	if (eval_operator)
		Psi = (double *) malloc(Nspace * sizeof(double));
	else
		Psi = NULL;

	if (solveStokes) {
		Ipol = matrix_double(4, Nspace);
		I    = Ipol[0];
		Spol = matrix_double(4, Nspace);
		S    = Spol[0];
	} else {
		I = (double *) malloc(Nspace * sizeof(double));
		S = (double *) malloc(Nspace * sizeof(double));
	}
	chi = (double *) malloc(Nspace * sizeof(double));

	/* --- Store current mean intensity, initialize new one to zero - - */

	Jdag = (double *) malloc(Nspace * sizeof(double));
	if (input.limit_memory) {
		J = (double *) malloc(atmos.Nspace * sizeof(double));
		readJlambda(nspect, Jdag);
	} else {
		J = spectrum.J[nspect];
		for (k = 0; k < Nspace; k++) Jdag[k] = J[k];
	}
	if (spectrum.updateJ)
		for (k = 0; k < Nspace; k++) J[k] = 0.0;

	/* --- The background scattering source: from Jdag, or, where the
		   previous sweep found that the background scattering at this
		   wavelength still needs it, from J iterated with the
		   populations fixed (N-51) before this sweep, so that J, the
		   rates, Jgas and the emergent intensity all come from the
		   same formal solution --                      -------------- */

	Jsrc = Jdag;
	if (bg) {
		bgc.nspect     = nspect;
		bgc.angle_dep  = angle_dep;
		bgc.boundbound = boundbound;
		bgc.chi        = chi;
		bgc.S          = S;
		bgc.I          = I;
		bgc.Tri        = Tri;
		bgc.as         = as;
		pending        = bgScatterPending(nspect);
	}
	if (bg && angle_dep && input.S_interpolation != S_LINEAR) {

		/* --- With S_BEZIER3, whether this wavelength is solved and
			   iterated linear (N-52), from the opacities of the first
			   ray, before its iteration --              -------------- */

		readBackground(nspect, 0, 0);
		Opacity(nspect, 0, 0, TRUE);
		lin = linearChoice(nspect, &bg, &Tri, as);
	}
	if (bg && pending && angle_dep) {
		Jstar = (double *) malloc(Nspace * sizeof(double));
		bgScatterBefore(nspect, Nspace, Jdag, Jstar, bgSweep_1D, &bgc);
		Jsrc = Jstar;
	}

	/* --- Store current anisotropy, initialize new one to zero ---- -- */

	if (input.backgr_pol) {
		J20dag = (double *) malloc(Nspace * sizeof(double));
		if (input.limit_memory) {
			J20 = (double *) malloc(Nspace * sizeof(double));
			readJ20lambda(nspect, J20dag);
		} else {
			J20 = spectrum.J20[nspect];
			for (k = 0; k < Nspace; k++)
				J20dag[k] = J20[k];
		}
		if (spectrum.updateJ)
			for (k = 0; k < Nspace; k++) J20[k] = 0.0;
	}
	/* --- Case of angle-dependent opacity and source function -- ----- */

	if (angle_dep) {
		for (mu = 0; mu < Nrays; mu++) {
			wmu = 0.5 * geometry.wmu[mu];
			if (input.backgr_pol) {
				musq     = SQ(geometry.muz[mu]);
				threemu1 = TWOSQRTTWO * (3.0 * musq - 1.0);
				threemu2 = (3.0 * TWOSQRTTWO) * (musq - 1.0);
			}
			for (to_obs = 0; to_obs <= 1; to_obs++) {
				initialize = (mu == 0 && to_obs == 0);

				if (initialize || atmos.backgrflags[nspect].hasline)
					readBackground(nspect, mu, to_obs);

				if (initialize || boundbound)
					Opacity(nspect, mu, to_obs, initialize);
				if (initialize) lin = linearChoice(nspect, &bg, &Tri, as);

				if (eval_operator) addtoCoupling(nspect);
				for (k = 0; k < Nspace; k++) {
					chi[k] = as->chi[k] + as->chi_c[k];
					S[k]   = as->eta[k] + as->eta_c[k] + as->sca_c[k] * Jsrc[k];
				}
				if (negTauColumns(nspect, mu, to_obs, Nspace, 1, chi, geometry.height,
								  geometry.muz[mu], quad_ray))
					goto cleanup;
				if (input.backgr_pol) {
					for (k = 0; k < Nspace; k++) {
						S[k] += threemu1 * as->sca_c[k] * J20dag[k];
					}
				}
				if (solveStokes) {
					for (k = Nspace; k < 4 * Nspace; k++) Spol[0][k] = 0.0;

					/* --- Add emissivity due to active set for Q, U, V -- ---- */

					if (polarized_as) {
						for (k = Nspace; k < 4 * Nspace; k++)
							Spol[0][k] += as->eta[k];
					}
					/* --- Add emissivity due to background lines -- ---------- */

					if (polarized_c) {
						for (k = Nspace; k < 4 * Nspace; k++)
							Spol[0][k] += as->eta_c[k];
					}
					/* --- Add emissivity due to background scattering -- ----- */

					if (input.backgr_pol) {
						for (k = 0; k < Nspace; k++) {
							Spol[1][k] += threemu2 * as->sca_c[k] * J20dag[k];
						}
					}
					for (n = 0; n < 4; n++) {
						for (k = 0; k < Nspace; k++)
							Spol[n][k] /= chi[k];
					}
					/* --- Polarized transfer --                 -------------- */

					if (input.S_interpolation_stokes == DELO_BEZIER3) {
						Piece_Stokes_Bezier3_1D(nspect, mu, to_obs,
												chi, Spol, Ipol, Psi);
					} else {
						sprintf(messageStr,
								"Unknown polarization solver: %d",
								input.S_interpolation_stokes);
						Error(ERROR_LEVEL_1, routineName, messageStr);
					}
				} else {

					for (k = 0; k < Nspace; k++)
						S[k] /= chi[k];

					/* --- Intensity only --                     -------------- */

					Piecewise_1D(nspect, mu, to_obs, chi, S, I, Psi, Tri, lin);
				}

				if (Psi) {
					for (k = 0; k < Nspace; k++) Psi[k] /= chi[k];
					if (eval_operator) addtoGamma(nspect, wmu, I, Psi);
				}
				if (bg) {
					addTri(Nspace, to_obs, wmu, Tri, as->sca_c, chi, Abg);
					for (k = 0; k < Nspace; k++)
						albedo[k] += wmu * as->sca_c[k] / chi[k];
				}

				if (spectrum.updateJ) {

					/* --- Accumulate mean intensity and rates -- ----------- */

					for (k = 0; k < Nspace; k++)
						J[k] += wmu * I[k];
					addtoRates(nspect, mu, to_obs, wmu, I, redistribute);

					/* --- Accumulate anisotropy --            -------------- */

					if (input.backgr_pol) {
						if (solveStokes) {
							for (k = 0; k < Nspace; k++) {
								J20[k] += (threemu1 * Ipol[0][k] +
										   threemu2 * Ipol[1][k]) *
										  wmu;
							}
						} else {
							for (k = 0; k < Nspace; k++) {
								J20[k] += threemu1 * I[k] * wmu;
							}
						}
					}

					/* --- Accumulate gas-frame mean intensity --   ----------- */

					if (atmos.NPRDactive > 0 &&
						input.PRD_angle_dep == PRD_ANGLE_APPROX) {

						addJgasRay(nspect, mu, to_obs, wmu, I);
					}

					if (containsPRDline(as) &&
						input.PRD_angle_dep == PRD_ANGLE_DEP)
						writeImu(nspect, mu, to_obs, I);
				}
			}
			/* --- Save emergent intensity --              -------------- */

			spectrum.I[nspect][mu] = I[0];
			if (solveStokes) {
				spectrum.Stokes_Q[nspect][mu] = Ipol[1][0];
				spectrum.Stokes_U[nspect][mu] = Ipol[2][0];
				spectrum.Stokes_V[nspect][mu] = Ipol[3][0];
			}
		}
	} else {

		/* --- The angle-independent case: opacity and source function
			   are the same for every ray, but each direction is still
			   solved with the method of S_INTERPOLATION (N-23) -- ---- */

		readBackground(nspect, 0, 0);
		Opacity(nspect, 0, 0, initialize = TRUE);
		lin = linearChoice(nspect, &bg, &Tri, as);
		if (eval_operator) addtoCoupling(nspect);

		for (k = 0; k < Nspace; k++) chi[k] = as->chi[k] + as->chi_c[k];
		for (mu = 1, muslant = 0; mu < Nrays; mu++) /* the most slanted ray */
			if (geometry.muz[mu] < geometry.muz[muslant]) muslant = mu;
		if (negTauColumns(nspect, muslant, -1, Nspace, 1, chi, geometry.height,
						  geometry.muz[muslant], quad))
			goto cleanup;
		if (bg && pending) {
			Jstar = (double *) malloc(Nspace * sizeof(double));
			bgScatterBefore(nspect, Nspace, Jdag, Jstar, bgSweep_1D, &bgc);
			Jsrc = Jstar;
		}
		for (k = 0; k < Nspace; k++)
			S[k] = (as->eta[k] + as->eta_c[k] + as->sca_c[k] * Jsrc[k]) / chi[k];

		for (mu = 0; mu < Nrays; mu++) {
			wmu = 0.5 * geometry.wmu[mu];
			for (to_obs = 0; to_obs <= 1; to_obs++) {
				Piecewise_1D(nspect, mu, to_obs, chi, S, I, Psi, Tri, lin);

				if (Psi) {
					for (k = 0; k < Nspace; k++) Psi[k] /= chi[k];
					if (eval_operator) addtoGamma(nspect, wmu, I, Psi);
				}
				if (bg) {
					addTri(Nspace, to_obs, wmu, Tri, as->sca_c, chi, Abg);
					for (k = 0; k < Nspace; k++)
						albedo[k] += wmu * as->sca_c[k] / chi[k];
				}

				if (spectrum.updateJ) {
					for (k = 0; k < Nspace; k++) J[k] += wmu * I[k];
					addtoRates(nspect, mu, to_obs, wmu, I, redistribute);

					/* --- Accumulate gas-frame mean intensity, which is the
						   same as J in the angle-independent case -- ----- */

					if (atmos.NPRDactive > 0 &&
						input.PRD_angle_dep == PRD_ANGLE_APPROX) {
						addJgasRow(nspect, wmu, I);
					}
				}
			}
			spectrum.I[nspect][mu] = I[0];
		}
	}
	/* --- Write new J for current position in the spectrum -- -------- */

	dJmax = 0.0;
	if (spectrum.updateJ) {
		for (k = 0; k < Nspace; k++) {
			/* --- J can be zero where every emissivity along the rays has
			   underflowed (far in the EUV): count it as no change -- */

			if (J[k] != 0.0)
				dJ = fabs(1.0 - Jdag[k] / J[k]);
			else
				dJ = (Jdag[k] != 0.0) ? 1.0 : 0.0;
			dJmax = MAX(dJmax, dJ);
		}
		/* --- Whether the next sweep iterates the background
			   scattering at this wavelength first (N-51) -- -------- */

		if (bg && !pending) bgScatterAfter(nspect, Nspace, Jsrc, J, Abg, albedo);
		if (input.limit_memory) {
			writeJlambda(nspect, J);
			if (input.backgr_pol)
				writeJ20lambda(nspect, J20);
		}
	}
	/* --- Clean up; the formal solution of this wavelength stops here
		   when the opacity is negative enough to overflow it (N-40;
		   the caller stops the run or rejects the trial point) -- --- */

cleanup:
	free_as(nspect, eval_operator);
	if (eval_operator) free(Psi);
	free(Abg);

	free(chi);
	if (solveStokes) {
		freeMatrix((void **) Ipol);
		freeMatrix((void **) Spol);
	} else {
		free(I);
		free(S);
	}

	free(Jdag);
	free(Jstar);
	if (input.limit_memory) free(J);
	if (input.backgr_pol) {
		free(J20dag);
		if (input.limit_memory) free(J20);
	}

	return dJmax;
}
/* ------- end ---------------------------- Formal.c ---------------- */
