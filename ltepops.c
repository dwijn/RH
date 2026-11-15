/* ------- file: -------------------------- ltepops.c ---------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu May  6 15:43:59 2021 --

	   --------------------------                      ----------RH-- */

/* --- Various routines to calculate LTE populations -- ------------- */

#include <math.h>
#include <stdlib.h>
#include <float.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"
#include "statistics.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern char       messageStr[];

/* ------- begin -------------------------- lognormalize.c --------- */

static void lognormalize(double **x, int N, int k, double total, bool_t clip)
{
	/* --- Turns the logarithms x[i][k], i = 0..N-1, of a set of
		   populations known up to a common factor into the populations
		   that sum to total: exp(x - max x) is at most 1, so that the
		   largest is exact and the others underflow only on their
		   own, whatever the spread of the logs. With clip set, no
		   population is smaller than DBL_MIN: the last guard against
		   exact zeros in ratios of populations (but not where total is
		   zero, then all are zero).
		   --                                            -------------- */

	register int i;

	double xmax = x[0][k], sum = 0.0;

	for (i = 1; i < N; i++) xmax = MAX(xmax, x[i][k]);
	for (i = 0; i < N; i++) {
		x[i][k] = exp(x[i][k] - xmax);
		sum += x[i][k];
	}
	for (i = 0; i < N; i++) {
		x[i][k] = (total > 0.0) ? x[i][k] * (total / sum) : 0.0;
		if (clip && total > 0.0) x[i][k] = MAX(x[i][k], DBL_MIN);
	}
}
/* ------- end ---------------------------- lognormalize.c --------- */

/* ------- begin -------------------------- LTEpops.c --------------- */

void LTEpops(Atom *atom, bool_t Debeye)
{
	register int k, i, m;

	char   labelStr[MAX_LINE_SIZE];
	int    Z, *dZ, *nDebeye, Nspace = atmos.Nspace;
	double lcNe_T, dE_kT, dEion = 0.0, c1, c2 = 0.0, *lgi0, *dE;

	/* --- Computes LTE populations of a given atom.
		   Takes account of Debeye shielding and lowering of the ionization
		   potential when Debeye is set to TRUE.

		   The populations are computed in logarithms, relative to the
		   largest term (lognormalize), so that they do not underflow or
		   overflow in a many-stage atom when the stages differ by many
		   orders of magnitude: the dominant level is then exact and a
		   minor one underflows only on its own account.

			 dE_ion = -Z * (e^2/4PI*EPSILON_0) / D           [J]
				  D = sqrt(EPSILON_0/(2e^2)) * sqrt(kT/ne)   [m]

	  See: Mihalas (78), pp. 293-295

		 --                                              -------------- */

	getCPU(3, TIME_START, NULL);

	/* --- Depth-independent constants --               --------------- */

	c1 = (HPLANCK / (2.0 * PI * M_ELECTRON)) * (HPLANCK / KBOLTZMANN);

	/* --- Determine the total lowering of ionization potential due
		   to Debeye shielding (zero without), and the quantities that
		   only depend on the level --                  --------------- */

	if (Debeye)
		c2 = sqrt(8.0 * PI / KBOLTZMANN) *
			 pow(SQ(Q_ELECTRON) / (4.0 * PI * EPSILON_0), 1.5);

	nDebeye = (int *) calloc(atom->Nlevel, sizeof(int));
	dZ      = (int *) malloc(atom->Nlevel * sizeof(int));
	lgi0    = (double *) malloc(atom->Nlevel * sizeof(double));
	dE      = (double *) malloc(atom->Nlevel * sizeof(double));
	for (i = 0; i < atom->Nlevel; i++) {
		dZ[i]   = atom->stage[i] - atom->stage[0];
		dE[i]   = atom->E[i] - atom->E[0];
		lgi0[i] = log(atom->g[i] / atom->g[0]);
		if (Debeye && i > 0) {
			Z = atom->stage[i];
			for (m = 1; m <= dZ[i]; m++, Z++) nDebeye[i] += Z;
		}
	}
	/* --- Solve Saha-Boltzmann equilibrium equations --  ------------- */

	for (k = 0; k < Nspace; k++) {
		if (atmos.ne[k] <= 0.0) {
			sprintf(messageStr,
					"Electron density must be positive for LTE populations "
					"(atom %s, ne[%d] = %g)\n",
					atom->ID, k, atmos.ne[k]);
			Error(ERROR_LEVEL_2, "LTEpops", messageStr);
		}
		if (Debeye) dEion = c2 * sqrt(atmos.ne[k] / atmos.T[k]);
		lcNe_T = log(0.5 * atmos.ne[k]) + 1.5 * log(c1 / atmos.T[k]);

		/* --- The logarithm of nstar[i] / nstar[0] --  -------------- */

		atom->nstar[0][k] = 0.0;
		for (i = 1; i < atom->Nlevel; i++) {
			dE_kT = (dE[i] - nDebeye[i] * dEion) / (KBOLTZMANN * atmos.T[k]);

			atom->nstar[i][k] = lgi0[i] - dE_kT;
			if (dZ[i] != 0) atom->nstar[i][k] -= dZ[i] * lcNe_T;
		}
		lognormalize(atom->nstar, atom->Nlevel, k, atom->ntotal[k], TRUE);
	}

	free(nDebeye);
	free(dZ);
	free(lgi0);
	free(dE);

	sprintf(labelStr, "LTEpops %2s", atom->ID);
	getCPU(3, TIME_POLL, labelStr);
}
/* ------- end ---------------------------- LTEpops.c --------------- */

/* ------- begin -------------------------- LTEpops_elem.c ---------- */

void LTEpops_elem(Element *element)
{
	register int k, i;

	bool_t  hunt;
	double *Uk, *Ukp1, C1, *lCT_ne, ln2 = log(2.0);

	getCPU(4, TIME_START, NULL);

	C1 = (HPLANCK / (2.0 * PI * M_ELECTRON)) * (HPLANCK / KBOLTZMANN);

	lCT_ne = (double *) malloc(atmos.Nspace * sizeof(double));
	Uk     = (double *) malloc(atmos.Nspace * sizeof(double));
	Ukp1   = (double *) malloc(atmos.Nspace * sizeof(double));

	for (k = 0; k < atmos.Nspace; k++) {
		if (atmos.ne[k] <= 0.0) {
			sprintf(messageStr,
					"Electron density must be positive for LTE populations "
					"(element %s, ne[%d] = %g)\n",
					element->ID, k, atmos.ne[k]);
			Error(ERROR_LEVEL_2, "LTEpops_elem", messageStr);
		}
		/* --- ln(CT_ne) of the Saha equation; n[i] holds ln(n_i/n_0) -- */

		lCT_ne[k]        = ln2 - 1.5 * log(C1 / atmos.T[k]) - log(atmos.ne[k]);
		element->n[0][k] = 0.0;
	}

	Linear(atmos.Npf, atmos.Tpf, element->pf[0],
		   atmos.Nspace, atmos.T, Uk, hunt = TRUE);

	for (i = 1; i < element->Nstage; i++) {
		Linear(atmos.Npf, atmos.Tpf, element->pf[i],
			   atmos.Nspace, atmos.T, Ukp1, hunt = TRUE);

		for (k = 0; k < atmos.Nspace; k++) {
			element->n[i][k] = element->n[i - 1][k] + lCT_ne[k] + Ukp1[k] - Uk[k] -
							   element->ionpot[i - 1] / (KBOLTZMANN * atmos.T[k]);
		}
		SWAPPOINTER(Uk, Ukp1);
	}

	/* --- Normalize relative to the most populated stage, so that
		   neither the products of many Saha factors nor their sum can
		   underflow or overflow: the dominant stages are exact, the
		   minor ones underflow only on their own -- -------------- */

	for (k = 0; k < atmos.Nspace; k++)
		lognormalize(element->n, element->Nstage, k,
					 element->abund * atmos.nHtot[k], FALSE);

	free(lCT_ne);
	free(Uk);
	free(Ukp1);
}
/* ------- end ---------------------------- LTEpops_elem.c ---------- */

/* ------- begin -------------------------- LTEmolecule.c ----------- */

void LTEmolecule(Molecule *molecule)
{
	/* --- Calculate partition functions for each molecular vibrational
		   state v of the molecule. LTE populations are then given by:

		   nv^*[k] = molecule->n * pfv[v][k] / pf[k].

	 Note: The actual LTE populations are calculated (in initSolution)
		   only after chemical equilibrium has been established.
		   --                                            -------------- */

	register int k, v, J, kr;

	char           labelStr[MAX_LINE_SIZE];
	double         gJ, **E;
	MolecularLine *mrt;

	if (!molecule->active) {
		sprintf(messageStr, "Molecule must be active: %s\n", molecule->ID);
		Error(ERROR_LEVEL_2, "LTEmolecule", messageStr);
	}
	/* --- Fill energy matrix --                         -------------- */

	E = matrix_double(molecule->Nv, molecule->NJ);
	for (kr = 0; kr < molecule->Nrt; kr++) {
		mrt = molecule->mrt + kr;

		E[mrt->vi][(int) (mrt->gi - 1) / 2] = mrt->Ei;
		E[mrt->vj][(int) (mrt->gj - 1) / 2] = mrt->Ej;
	}

	for (k = 0; k < atmos.Nspace; k++)
		molecule->pf[k] = 0.0;

	for (v = 0; v < molecule->Nv; v++) {
		for (J = 0; J < molecule->NJ; J++) {
			gJ = 2 * J + 1;
			for (k = 0; k < atmos.Nspace; k++)
				molecule->pfv[v][k] +=
					gJ * exp(-E[v][J] / (KBOLTZMANN * atmos.T[k]));
		}
		/*  --- Also store the total partition function here -- --------- */

		for (k = 0; k < atmos.Nspace; k++)
			molecule->pf[k] += molecule->pfv[v][k];
	}

	freeMatrix((void **) E);

	sprintf(labelStr, "LTEpops %3s", molecule->ID);
	getCPU(4, TIME_POLL, labelStr);
}
/* ------- end ---------------------------- LTEmolecule.c ----------- */

/* ------- begin -------------------------- SetLTEQuantities.c ------ */

void SetLTEQuantities(void)
{
	register int n;

	bool_t Debeye = TRUE;
	Atom  *atom;

	for (n = 0; n < atmos.Natom; n++) {
		atom = &atmos.atoms[n];

		/* --- Get LTE populations for each atom --        -------------- */

		LTEpops(atom, Debeye);

		if (atom->active) {

			/* --- Read the collisional data (in MULTI's GENCOL format).
				   The input file of an active atom stays open, and this is
				   called again whenever the electron density changes
				   (SOLVE_NE = ITERATION), so rewind to the start of the
				   collisional data first --                 -------------- */

			if (fseek(atom->fp_input, atom->offset_coll, SEEK_SET) != 0) {
				sprintf(messageStr, "Cannot rewind to the collisional data of atom %s",
						atom->ID);
				Error(ERROR_LEVEL_2, "SetLTEQuantities", messageStr);
			}
			CollisionRate(atom, atom->fp_input);

			/* --- Compute the fixed rates and store in Cij -- ------------ */

			if (atom->Nfixed > 0) FixedRate(atom);
		}
	}
}
/* ------- end ---------------------------- SetLTEQuantities.c ------ */
