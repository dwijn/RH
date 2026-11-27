/* ------- file: -------------------------- stokesgrid.c ------------

	   Grid arrays for the cone-preserving interpolation of the
	   polarized 2-D and 3-D formal solvers (B-63); see stokesgrid.h.

	   --------------------------                      ----------RH-- */

#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"
#include "stokesgrid.h"
#include "state.h"

/* --- |p| of a grid value is taken at most P_MAX (atanh = 11.9):
	   only a vector already outside the cone (or on it) reaches it,
	   and its direction is kept -- -------------------------------- */

#define P_MAX (1.0 - 1.0E-10)

/* ------- begin -------------------------- Stokes_to_rapidity.c ---- */

void Stokes_to_rapidity(double I, double Q, double U, double V, double w[3])
{
	/* --- w = atanh(|p|) p / |p| with p = (Q, U, V) / I; w = 0 where
		   I <= 0, where there is no physical fraction (a negative
		   source function of an inverted line) -- -------------------- */

	double p, f;

	if (!(I > 0.0)) {
		w[0] = w[1] = w[2] = 0.0;
		return;
	}
	p = sqrt(Q * Q + U * U + V * V) / I;
	if (p < 1.0E-8)
		f = 1.0;
	else
		f = atanh(MIN(p, P_MAX)) / p;
	f /= I;
	w[0] = f * Q;
	w[1] = f * U;
	w[2] = f * V;
}
/* ------- end ---------------------------- Stokes_to_rapidity.c ---- */

/* ------- begin -------------------------- Stokes_from_rapidity.c -- */

void Stokes_from_rapidity(double I, const double w[3], double v[4])
{
	/* --- v = (I, I tanh(|w|) w / |w|): |(Q, U, V)| < |I| for every w
		   --                                              -------------- */

	double r = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]),
		   f = (r < 1.0E-4) ? 1.0 - r * r / 3.0 : tanh(r) / r;

	v[0] = I;
	v[1] = I * f * w[0];
	v[2] = I * f * w[1];
	v[3] = I * f * w[2];
}
/* ------- end ---------------------------- Stokes_from_rapidity.c -- */

/* ------- begin -------------------------- StokesGrid_init.c ------- */

void StokesGrid_init(RHState *rh, StokesGrid *sg, int nspect, int Nspace, double *chi,
					 double **S, double **I)
{
	/* --- The grid arrays of one formal solution: the rapidities of S
		   and I, and eta' and rho' from StokesK, at every grid point.
		   I must hold the boundary values; the rest is updated as the
		   sweep goes (StokesGrid_update_I) -- -------------------------- */

	ActiveSet *as = &rh->spectrum.as[nspect];
	double     K[4][4], w[3];

	sg->nspect = nspect;
	sg->Nspace = Nspace;
	sg->chi    = chi;
	sg->S      = S;
	sg->I      = I;

	sg->polarizedK      = containsPolarized(rh, as) || rh->atmos.backgrflags[nspect].ispolarized;
	sg->magneto_optical = sg->polarizedK && rh->input.magneto_optical;

	sg->buffer = (double *) malloc(12 * (size_t) Nspace * sizeof(double));
	for (int n = 0; n < 3; n++) {
		sg->wS[n] = sg->buffer + n * Nspace;
		sg->wI[n] = sg->buffer + (3 + n) * Nspace;
		sg->wK[n] = sg->buffer + (6 + n) * Nspace;
		sg->rK[n] = sg->buffer + (9 + n) * Nspace;
	}
	for (int k = 0; k < Nspace; k++) {
		Stokes_to_rapidity(S[0][k], S[1][k], S[2][k], S[3][k], w);
		for (int n = 0; n < 3; n++) sg->wS[n][k] = w[n];
		StokesGrid_update_I(sg, k);

		if (sg->polarizedK) {
			StokesK(rh, nspect, k, chi[k], K);
			Stokes_to_rapidity(1.0, K[0][1], K[0][2], K[0][3], w);
			for (int n = 0; n < 3; n++) sg->wK[n][k] = w[n];

			/* --- rho' as StokesK has it, signs included: K[1][2],
				   K[3][1] and K[2][3] -- ----------------------------- */

			sg->rK[0][k] = K[1][2];
			sg->rK[1][k] = K[3][1];
			sg->rK[2][k] = K[2][3];
		} else {
			for (int n = 0; n < 3; n++) sg->wK[n][k] = sg->rK[n][k] = 0.0;
		}
	}
}
/* ------- end ---------------------------- StokesGrid_init.c ------- */

/* ------- begin -------------------------- StokesGrid_update_I.c --- */

void StokesGrid_update_I(StokesGrid *sg, int k)
{
	double w[3];

	Stokes_to_rapidity(sg->I[0][k], sg->I[1][k], sg->I[2][k], sg->I[3][k], w);
	for (int n = 0; n < 3; n++) sg->wI[n][k] = w[n];
}
/* ------- end ---------------------------- StokesGrid_update_I.c --- */

/* ------- begin -------------------------- StokesGrid_free.c ------- */

void StokesGrid_free(StokesGrid *sg)
{
	free(sg->buffer);
	sg->buffer = NULL;
}
/* ------- end ---------------------------- StokesGrid_free.c ------- */
