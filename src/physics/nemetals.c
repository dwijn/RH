/* ------- file: -------------------------- nemetals.c --------------

	   Version:       rh1.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Oct 13 10:49:53 2009 --

	   --------------------------                      ----------RH-- */

#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "background.h"
#include "state.h"

#define NEFRACTION 0.01

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- neMetals.c -------------- */

void FMetals(RHState *rh, double *F)
{
	/* --- Calculates the effective ionization fraction F_met for all
		   background metals (i.e. excluding hydrogen) under the
		   assumption of fixed temperature and total electron density.
		   i.e. F_met = F_met(T, ne)

		   ne_metal = nHtot * Sum_n (A_n Sum_l (l * f_nl)) = nHtot * F_met

	 Note: Populations of background metals may be either LTE or NonLTE!
		   --                                            -------------- */

	register int k, i, n;

	double *f_n;
	Atom   *metal;

	f_n = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	for (k = 0; k < rh->atmos.Nspace; k++) F[k] = 0.0;

	for (n = 1; n < rh->atmos.Natom; n++) {
		metal = &rh->atmos.atoms[n];
		for (k = 0; k < rh->atmos.Nspace; k++) f_n[k] = 0.0;

		for (i = 0; i < metal->Nlevel; i++) {
			if (metal->stage[i] > 0) {
				for (k = 0; k < rh->atmos.Nspace; k++)
					f_n[k] += metal->stage[i] * metal->n[i][k];
			}
		}

		for (k = 0; k < rh->atmos.Nspace; k++)
			F[k] += metal->abundance * f_n[k] / metal->ntotal[k];
	}

	free(f_n);
}
/* ------- end ---------------------------- neMetals.c -------------- */

/* ------- begin -------------------------- dFMetals.c -------------- */

void dFMetals(RHState *rh, double *dFdne)
{
	register int n, k;

	bool_t  Debeye;
	double *neold, *Fmin, *Fplus;

	/* --- Numerically evaluate the derivative dF/dne -- -------------- */

	neold = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	Fmin  = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	Fplus = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));

	/* --- Store electron densities, evaluate LTE populations for reduced
		   electron density --                           -------------- */

	for (k = 0; k < rh->atmos.Nspace; k++) {
		neold[k]        = rh->atmos.ne[k];
		rh->atmos.ne[k] = (1.0 - NEFRACTION) * neold[k];
	}
	for (n = 1; n < rh->atmos.Natom; n++)
		LTEpops(rh, &rh->atmos.atoms[n], Debeye = TRUE);

	/* --- Get reduced ionization fraction F --          -------------- */

	FMetals(rh, Fmin);

	/* --- evaluate LTE populations for enhanced electron density -- -- */

	for (k = 0; k < rh->atmos.Nspace; k++)
		rh->atmos.ne[k] = (1.0 + NEFRACTION) * neold[k];
	for (n = 1; n < rh->atmos.Natom; n++)
		LTEpops(rh, &rh->atmos.atoms[n], Debeye = TRUE);

	/* --- Get enhanced ionization fraction F --         -------------- */

	FMetals(rh, Fplus);

	/* --- get the numerical derivative, restore electron densities - - */

	for (k = 0; k < rh->atmos.Nspace; k++) {
		dFdne[k]        = (Fplus[k] - Fmin[k]) / (2.0 * NEFRACTION * neold[k]);
		rh->atmos.ne[k] = neold[k];
	}
	/* --- Restore the LTE populations --                 ------------- */

	for (n = 1; n < rh->atmos.Natom; n++)
		LTEpops(rh, &rh->atmos.atoms[n], Debeye = TRUE);

	free(neold);
	free(Fmin);
	free(Fplus);
}
/* ------- end ---------------------------- dFMetals.c -------------- */
