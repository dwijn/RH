/* ------- file: -------------------------- formal.c ----------------

	   Version:       rh2.0, 1-D spherically symmetric
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Dec 14 16:40:20 2011 --

	   --------------------------                      ----------RH-- */

/* --- Formal solution with given source function and PRD emission
	   profile. --                                     -------------- */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "error.h"
#include "inputs.h"
#include "bgscatter.h"
#include "negtau.h"
#include "background.h"
#include "constant.h"
#include "statistics.h"
#include "xdr.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Geometry   geometry;
extern Spectrum   spectrum;
extern InputData  input;
extern char       messageStr[];

/* ------- begin -------------------------- bgScatterSupported.c ---- */

const char *bgScatterSupported(void)
{
	/* --- The Feautrier operator is positive (an M-matrix, see
		   feautrier.c), so the background scattering iteration (N-51)
		   needs no fallback here -- --                 -------------- */

	return NULL;
}
/* ------- end ---------------------------- bgScatterSupported.c ---- */

/* ------- begin -------------------------- bgScatterLinearFallback.c */

bool_t bgScatterLinearFallback(void)
{
	return FALSE;
}
/* ------- end ---------------------------- bgScatterLinearFallback.c */

/* ------- begin -------------------------- addTri_sphere.c --------- */

static void addTri_sphere(const Ray *ray, const double *Tri,
						  const double *sca_c, const double *chi, double *A,
						  double *albedo)
{
	register int k;

	double w;

	/* --- Add the rows of the operator of one ray (Tri[3k + i] =
		   dP_k / dS_k-1+i, Feautrier) to A[3k + i] = dJ_k / du_k-1+i,
		   with S_j = ... + sca_c,j u_j / chi_j (N-51), and its weight
		   of sca_c / chi to the angle-averaged albedo (if not NULL). The ray passes
		   the shells 0 .. Ns - 1 -- --                  -------------- */

	for (k = 0; k < ray->Ns; k++) {
		w = ray->wmu[k];
		if (k > 0) A[3 * k] += w * Tri[3 * k] * sca_c[k - 1] / chi[k - 1];
		A[3 * k + 1] += w * Tri[3 * k + 1] * sca_c[k] / chi[k];
		if (k < ray->Ns - 1)
			A[3 * k + 2] += w * Tri[3 * k + 2] * sca_c[k + 1] / chi[k + 1];
		if (albedo) albedo[k] += w * sca_c[k] / chi[k];
	}
}
/* ------- end ---------------------------- addTri_sphere.c --------- */

/* ------- begin -------------------------- bgSweep_sphere.c -------- */

/* --- One formal solution of all rays at nspect for the background
	   scattering iteration (N-51), with the background scattering
	   source sca_c u and the opacities already in as and chi (the
	   spherical case is angle independent): J, the tridiagonal A of
	   dJ/du and, if not NULL, the angle average of sca_c / chi; no
	   rates, no output --                                       -------------- */

typedef struct {
	int        nspect;
	double    *chi, *S, *P, *Tri;
	ActiveSet *as;
} BgCtxSphere;

static void bgSweep_sphere(void *ctx, const double *u, double *J, double *A,
						   double *albedo)
{
	register int k, mu;

	BgCtxSphere *c      = (BgCtxSphere *) ctx;
	ActiveSet   *as     = c->as;
	int          Nspace = atmos.Nspace;
	Ray         *ray;

	for (k = 0; k < Nspace; k++) {
		J[k]    = 0.0;
		c->S[k] = (as->eta[k] + as->eta_c[k] + as->sca_c[k] * u[k]) / c->chi[k];
		if (albedo) albedo[k] = 0.0;
	}
	for (k = 0; k < 3 * Nspace; k++) A[k] = 0.0;
	for (mu = 0; mu < atmos.Nrays; mu++) {
		ray = geometry.rays + mu;
		Feautrier(c->nspect, mu, c->chi, c->S, STANDARD, c->P, NULL, c->Tri);
		for (k = 0; k < ray->Ns; k++) J[k] += c->P[k] * ray->wmu[k];
		addTri_sphere(ray, c->Tri, as->sca_c, c->chi, A, albedo);
	}
}
/* ------- end ---------------------------- bgSweep_sphere.c -------- */

/* ------- begin -------------------------- Formal.c ---------------- */

double Formal(int nspect, bool_t eval_operator, bool_t redistribute)
{
	const char   routineName[] = "Formal";
	register int k, mu;

	bool_t initialize, boundbound, PRD_angle_dep, polarized, bg,
		pending = FALSE;
	enum FeautrierOrder F_order;
	int                 Nspace = atmos.Nspace, Nrays = atmos.Nrays;
	double             *I, *chi, *S, *Psi, *Jdag, dJ, dJmax = 0.0, *J, *Abg = NULL,
										  *albedo = NULL, *Tri = NULL, *Jsrc = NULL, *Jstar = NULL;
	ActiveSet  *as;
	Ray        *ray;
	BgCtxSphere bgc;

	/* --- Retrieve active set as of transitions at wavelength nspect - */

	as = &spectrum.as[nspect];
	alloc_as(nspect, eval_operator);

	/* --- Check whether current active set includes a bound-bound
		   and/or angledependent PRD transition. Otherwise, only
		   angle-independent opacity and source functions are needed -- */

	boundbound    = containsBoundBound(as);
	PRD_angle_dep = (containsPRDline(as) && input.PRD_angle_dep);
	polarized     = containsPolarized(as);

	/* --- Allocate temporary space --                   -------------- */

	chi = (double *) malloc(Nspace * sizeof(double));

	if (eval_operator)
		Psi = (double *) malloc(Nspace * sizeof(double));
	else
		Psi = NULL;

	S = (double *) malloc(Nspace * sizeof(double));
	I = (double *) malloc(Nspace * sizeof(double));

	/* --- Iterate the background scattering at this wavelength (N-51):
		   with the tridiagonal of dJ/du from every ray -- ---------- */

	bg = (spectrum.updateJ && bgScatterActive());
	if (bg) {
		Abg    = (double *) calloc(7 * Nspace, sizeof(double));
		albedo = Abg + 3 * Nspace;
		Tri    = albedo + Nspace;
	}

	/* --- Store current mean intensities, initialize new ones to zero- */

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

	/* --- Case of angle-dependent opacity and source function -- ----- */

	if (polarized || PRD_angle_dep ||
		(atmos.moving && (boundbound || atmos.backgrflags[nspect].hasline))) {
		Error(ERROR_LEVEL_2, routineName,
			  "Angle dependent source functions and opacities\n"
			  " not yet implemented in spherical symmetry");
	} else {

		/* --- The angle-independent case --               -------------- */

		readBackground(nspect, 0, 0);
		Opacity(nspect, 0, 0, initialize = TRUE);
		if (eval_operator) addtoCoupling(nspect);

		for (k = 0; k < Nspace; k++) chi[k] = as->chi[k] + as->chi_c[k];

		/* --- N-40: along each ray, through the radii 0 .. Ns - 1 at
			   path positions s. The Feautrier solution has no
			   exponential to overflow, but its difference equations lose
			   their meaning just the same --          -------------- */

		for (mu = 0; mu < Nrays; mu++)
			if (negTauColumns(nspect, mu, -1, geometry.rays[mu].Ns, 1, chi,
							  geometry.rays[mu].s, 1.0, NEGTAU_LINEAR))
				goto cleanup;
		/* --- The background scattering source: from Jdag, or, where
			   the previous sweep found that the background scattering
			   at this wavelength still needs it, from J iterated with
			   the populations fixed (N-51), so that J, the rates and
			   the emergent intensity all come from the same formal
			   solution --                              -------------- */

		Jsrc = Jdag;
		if (bg && (pending = bgScatterPending(nspect))) {
			bgc.nspect = nspect;
			bgc.chi    = chi;
			bgc.S      = S;
			bgc.P      = I;
			bgc.Tri    = Tri;
			bgc.as     = as;
			Jstar      = (double *) malloc(Nspace * sizeof(double));
			bgScatterBefore(nspect, Nspace, Jdag, Jstar, bgSweep_sphere, &bgc);
			Jsrc = Jstar;
		}
		for (k = 0; k < Nspace; k++)
			S[k] = (as->eta[k] + as->eta_c[k] + as->sca_c[k] * Jsrc[k]) / chi[k];

		for (mu = 0; mu < Nrays; mu++) {
			spectrum.I[nspect][mu] =
				Feautrier(nspect, mu, chi, S, F_order = STANDARD, I, Psi, Tri);

			ray = geometry.rays + mu;
			if (bg) addTri_sphere(ray, Tri, as->sca_c, chi, Abg, albedo);
			if (eval_operator) {
				for (k = 0; k < ray->Ns; k++) Psi[k] /= chi[k];
				addtoGamma_sphere(nspect, ray, I, Psi);
			}
			if (spectrum.updateJ) {
				for (k = 0; k < ray->Ns; k++)
					J[k] += I[k] * ray->wmu[k];
				addtoRates_sphere(nspect, ray, I, redistribute);
			}
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
		if (input.limit_memory) writeJlambda(nspect, J);
	}
	/* --- Clean up; the formal solution of this wavelength stops here
		   when the opacity is negative enough to break it (N-40; the
		   caller stops the run or rejects the trial point) -- ------- */

cleanup:

	if (input.limit_memory) free(J);
	free_as(nspect, eval_operator);

	free(chi);
	free(Jdag);
	free(Jstar);
	free(Abg);
	free(I);
	free(S);
	if (eval_operator) free(Psi);

	return dJmax;
}
/* ------- end ---------------------------- Formal.c ---------------- */

/* ------- begin -------------------------- writeFlux.c ------------- */

bool_t writeFlux(char *fileName)
{
	const char   routineName[] = "writeFlux";
	register int mu, nspect;

	bool_t  result = TRUE;
	int     Nrays  = geometry.Nrays;
	double *flux, *wmuz, scale_to_surface;
	FILE   *fp_flux;
	XDR     xdrs;

	/* --- Write the radiative flux in the z-direction -- ------------ */

	if ((fp_flux = fopen(fileName, "w")) == NULL) {
		sprintf(messageStr, "Unable to open output file %s", fileName);
		Error(ERROR_LEVEL_1, routineName, messageStr);
		return FALSE;
	}
	xdrstdio_create(&xdrs, fp_flux, XDR_ENCODE);

	flux = (double *) calloc(spectrum.Nspect, sizeof(double));
	wmuz = (double *) malloc(Nrays * sizeof(double));

	wmuz[0] = 0.25 * (SQ(geometry.rays[1].xmu[0]) - SQ(geometry.rays[0].xmu[0]));
	for (mu = 1; mu < Nrays - 1; mu++) {
		wmuz[mu] = 0.25 * (SQ(geometry.rays[mu + 1].xmu[0]) -
						   SQ(geometry.rays[mu - 1].xmu[0]));
	}
	wmuz[Nrays - 1] = 0.25 * (SQ(geometry.rays[Nrays - 1].xmu[0]) -
							  SQ(geometry.rays[Nrays - 2].xmu[0]));

	scale_to_surface = SQ((geometry.r[0] + geometry.Radius) / geometry.Radius);
	for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
		for (mu = 0; mu < geometry.Nrays; mu++)
			flux[nspect] += spectrum.I[nspect][mu] * wmuz[mu];
		flux[nspect] *= 2.0 * PI * scale_to_surface;
	}
	result &= xdr_vector(&xdrs, (char *) flux, spectrum.Nspect,
						 sizeof(double), (xdrproc_t) rh_xdr_double);

	xdr_destroy(&xdrs);
	fclose(fp_flux);
	free(wmuz);
	free(flux);

	return result;
}
/* ------- end ---------------------------- writeFlux.c ------------- */

/* ------- begin -------------------------- addtoRates_sphere.c ----- */

void addtoRates_sphere(int nspect, Ray *ray, double *P,
					   bool_t redistribute)
{
	register int n, k, nact;

	int              nt;
	double           twohnu3_c2, twohc, wmula, *Rij, *Rji;
	ActiveSet       *as;
	Atom            *atom;
	AtomicLine      *line;
	AtomicContinuum *continuum;
	pthread_mutex_t *rate_lock;

	twohc = 2.0 * HPLANCK * CLIGHT / CUBE(NM_TO_M);

	as = &spectrum.as[nspect];
	nt = nspect % input.Nthreads;

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];

		for (n = 0; n < as->Nactiveatomrt[nact]; n++) {
			switch (as->art[nact][n].type) {
			case ATOMIC_LINE:
				line = as->art[nact][n].ptype.line;
				if (redistribute && !line->PRD)
					Rij = NULL;
				else {
					Rij        = line->Rij;
					Rji        = line->Rji;
					twohnu3_c2 = line->Aji / line->Bji;

					rate_lock = &line->rate_lock;
				}
				break;

			case ATOMIC_CONTINUUM:
				if (redistribute)
					Rij = NULL;
				else {
					continuum  = as->art[nact][n].ptype.continuum;
					Rij        = continuum->Rij;
					Rji        = continuum->Rji;
					twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);

					rate_lock = &continuum->rate_lock;
				}
				break;

			default:
				Rij = NULL;
			}
			/* --- Convention: Rij is the rate for transition i -> j -- --- */

			if (Rij != NULL) {
				if (input.Nthreads > 1) pthread_mutex_lock(rate_lock);

				for (k = 0; k < ray->Ns; k++) {
					wmula = atom->rhth[nt].Vij[n][k] *
							atom->rhth[nt].wla[n][k] * ray->wmu[k];
					Rij[k] += P[k] * wmula;
					Rji[k] += atom->rhth[nt].gij[n][k] * (twohnu3_c2 + P[k]) * wmula;
				}
				if (input.Nthreads > 1) pthread_mutex_unlock(rate_lock);
			}
		}
	}
}
/* -------- end --------------------------- addtoRates_sphere.c ----- */

/* ------- begin -------------------------- addtoGamma_sphere.c ----- */

void addtoGamma_sphere(int nspect, Ray *ray, double *P, double *Psi)
{
	const char   routineName[] = "addtoGamma_sphere";
	register int n, k, nact;

	int    i, j, ij, ji, nt;
	double twohnu3_c2, twohc, wmula;

	Atom            *atom;
	Molecule        *molecule;
	AtomicLine      *line;
	AtomicContinuum *continuum;
	MolecularLine   *mrt;
	ActiveSet       *as;

	twohc = 2.0 * HPLANCK * CLIGHT / CUBE(NM_TO_M);

	as = &spectrum.as[nspect];
	nt = nspect % input.Nthreads;

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];

		for (n = 0; n < as->Nactiveatomrt[nact]; n++) {
			switch (as->art[nact][n].type) {
			case ATOMIC_LINE:
				line       = as->art[nact][n].ptype.line;
				i          = line->i;
				j          = line->j;
				twohnu3_c2 = line->Aji / line->Bji;
				break;

			case ATOMIC_CONTINUUM:
				continuum  = as->art[nact][n].ptype.continuum;
				i          = continuum->i;
				j          = continuum->j;
				twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);
				break;

			default:
				sprintf(messageStr, "Invalid transition type");
				Error(ERROR_LEVEL_1, routineName, messageStr);
				twohnu3_c2 = 0.0;
			}
			/* --- In case of atomic transitions, with the cross-coupling
				   terms (see addPreconditionedRates) --      -------------- */

			if (input.Nthreads > 1) pthread_mutex_lock(&atom->Gamma_lock);

			addPreconditionedRates(atom->Nlevel, i, j, twohnu3_c2,
								   atom->rhth[nt].Vij[n], atom->rhth[nt].wla[n],
								   atom->rhth[nt].gij[n], atom->n[i], atom->n[j],
								   atom->rhth[nt].Uji_down[i],
								   atom->rhth[nt].Uji_down[j], P, Psi, 1.0,
								   ray->wmu, ray->Ns, atom->Gamma);

			if (input.Nthreads > 1) pthread_mutex_unlock(&atom->Gamma_lock);
		}
	}
	for (nact = 0; nact < atmos.Nactivemol; nact++) {
		molecule = atmos.activemols[nact];

		for (n = 0; n < as->Nactivemolrt[nact]; n++) {
			switch (as->mrt[nact][n].type) {

			case VIBRATION_ROTATION:
				mrt        = as->mrt[nact][n].ptype.vrline;
				i          = mrt->vi;
				j          = mrt->vj;
				twohnu3_c2 = mrt->Aji / mrt->Bji;
				break;

			default:
				sprintf(messageStr, "Invalid transition type");
				Error(ERROR_LEVEL_1, routineName, messageStr);
				twohnu3_c2 = 0.0;
			}

			if (input.Nthreads > 1) pthread_mutex_lock(&molecule->Gamma_lock);

			/* --- In case of molecular vibration-rotation transitions -- - */

			ij = i * molecule->Nv + j;
			ji = j * molecule->Nv + i;

			for (k = 0; k < ray->Ns; k++) {
				if (molecule->n[k]) {
					wmula = molecule->rhth[nt].Vij[n][k] *
							molecule->rhth[nt].wla[n][k] * ray->wmu[k];
					molecule->Gamma[ji][k] += P[k] * wmula;
					molecule->Gamma[ij][k] += molecule->rhth[nt].gij[n][k] *
											  (twohnu3_c2 + P[k]) * wmula;
				}
			}
			if (input.Nthreads > 1) pthread_mutex_unlock(&molecule->Gamma_lock);
		}
	}
}
/* ------- end ---------------------------- addtoGamma_sphere.c ----- */

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
