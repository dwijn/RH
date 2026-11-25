/* ------- file: -------------------------- formal.c ----------------

	   The transfer driver: the formal solution at one wavelength, for
	   every geometry, through the operations of the geometry
	   (geometryOps, geometryops.h; docs/LIBRARY.md §3.4, Q-42). It was
	   the Formal() of rhf1d, rhsc2d and rhsc3d (Han Uitenbroek), which
	   differed only in what the geometry now provides.

	   --                                              -------------- */

/* --- Formal solution with given source function, and allowing for
	   a PRD emission profile, polarized line radiation, polarized
	   background lines, and scattering background polarization -- -- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "background.h"
#include "inputs.h"
#include "error.h"
#include "negtau.h"
#include "bgscatter.h"
#include "geometryops.h"

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum   spectrum;
extern InputData  input;

/* ------- begin -------------------------- bgSweep.c --------------- */

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
} BgCtx;

static void bgSweep(void *ctx, const double *u, double *J, double *A,
					double *albedo)
{
	register int k, mu;

	BgCtx     *c  = (BgCtx *) ctx;
	ActiveSet *as = c->as;
	bool_t     to_obs, initialize;
	int        Nspace = atmos.Nspace;
	double     wmu;

	for (k = 0; k < Nspace; k++) J[k] = 0.0;
	for (k = 0; k < 3 * Nspace; k++) A[k] = 0.0;
	if (albedo)
		for (k = 0; k < Nspace; k++) albedo[k] = 0.0;
	for (mu = 0; mu < atmos.Nrays; mu++) {
		wmu = 0.5 * geometryOps.wmu(mu);
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
			geometryOps.solve(c->nspect, mu, to_obs, c->chi, c->S, c->I, NULL,
							  c->Tri, TRUE);
			for (k = 0; k < Nspace; k++) J[k] += wmu * c->I[k];
			geometryOps.bg_add_operator(Nspace, to_obs, wmu, c->Tri, as->sca_c,
										c->chi, A);
			if (albedo)
				for (k = 0; k < Nspace; k++)
					albedo[k] += wmu * as->sca_c[k] / c->chi[k];
		}
	}
}
/* ------- end ---------------------------- bgSweep.c --------------- */

/* ------- begin -------------------------- linearChoice.c ---------- */

/* --- Where the geometry falls back to the linear solver (rhf1d with
	   S_BEZIER3), whether this wavelength is solved with it (N-52):
	   where its background scattering albedo (bg_max_albedo) exceeds
	   BACKGROUND_LINEAR_ALBEDO, decided in every sweep that iterates
	   the background scattering; other formal solutions keep the last
	   choice (bgScatterLinear). Only those wavelengths iterate the
	   background scattering (*bg, *Tri) -- --          -------------- */

static bool_t linearChoice(int nspect, bool_t *bg, double **Tri,
						   const ActiveSet *as)
{
	bool_t lin;

	if (!geometryOps.bg_linear_fallback()) return FALSE;
	lin = bgScatterLinear(nspect, (*bg) ? geometryOps.bg_max_albedo(as) : -1.0);
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
	BgCtx      bgc;
	ActiveSet *as;

	/* --- A geometry with a Formal() of its own (rhsphere) -- ------ */

	if (geometryOps.formal)
		return geometryOps.formal(nspect, eval_operator, redistribute);

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

	quad     = (input.S_interpolation == S_BEZIER3) ? geometryOps.negtau_bezier
													: NEGTAU_LINEAR;
	quad_ray = (solveStokes) ? geometryOps.negtau_bezier : quad;

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
	if (bg && angle_dep && geometryOps.bg_linear_fallback()) {

		/* --- With S_BEZIER3, whether this wavelength is solved and
			   iterated linear (N-52), from the opacities of the first
			   ray, before its iteration --              -------------- */

		readBackground(nspect, 0, 0);
		Opacity(nspect, 0, 0, TRUE);
		lin = linearChoice(nspect, &bg, &Tri, as);
	}
	if (bg && pending && angle_dep) {
		Jstar = (double *) malloc(Nspace * sizeof(double));
		bgScatterBefore(nspect, Nspace, Jdag, Jstar, bgSweep, &bgc);
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
			wmu = 0.5 * geometryOps.wmu(mu);
			if (input.backgr_pol) {
				musq     = SQ(geometryOps.muz(mu));
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
				if (geometryOps.check_tau(nspect, mu, to_obs, chi, quad_ray))
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

					geometryOps.solve_stokes(nspect, mu, to_obs, chi, Spol, Ipol, Psi);
				} else {

					for (k = 0; k < Nspace; k++)
						S[k] /= chi[k];

					/* --- Intensity only --                     -------------- */

					geometryOps.solve(nspect, mu, to_obs, chi, S, I, Psi, Tri, lin);
				}

				if (Psi) {
					for (k = 0; k < Nspace; k++) Psi[k] /= chi[k];
					if (eval_operator) addtoGamma(nspect, wmu, I, Psi);
				}
				if (bg) {
					geometryOps.bg_add_operator(Nspace, to_obs, wmu, Tri, as->sca_c, chi, Abg);
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

			geometryOps.emergent(nspect, mu, I, (solveStokes) ? Ipol : NULL);
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
			if (geometryOps.muz(mu) < geometryOps.muz(muslant)) muslant = mu;
		if (geometryOps.check_tau(nspect, muslant, -1, chi, quad))
			goto cleanup;
		if (bg && pending) {
			Jstar = (double *) malloc(Nspace * sizeof(double));
			bgScatterBefore(nspect, Nspace, Jdag, Jstar, bgSweep, &bgc);
			Jsrc = Jstar;
		}
		for (k = 0; k < Nspace; k++)
			S[k] = (as->eta[k] + as->eta_c[k] + as->sca_c[k] * Jsrc[k]) / chi[k];

		for (mu = 0; mu < Nrays; mu++) {
			wmu = 0.5 * geometryOps.wmu(mu);
			for (to_obs = 0; to_obs <= 1; to_obs++) {
				geometryOps.solve(nspect, mu, to_obs, chi, S, I, Psi, Tri, lin);

				if (Psi) {
					for (k = 0; k < Nspace; k++) Psi[k] /= chi[k];
					if (eval_operator) addtoGamma(nspect, wmu, I, Psi);
				}
				if (bg) {
					geometryOps.bg_add_operator(Nspace, to_obs, wmu, Tri, as->sca_c, chi, Abg);
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
			geometryOps.emergent(nspect, mu, I, NULL);
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
