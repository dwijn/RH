/* ------- file: -------------------------- hydrostat.c -------------

	   Version:       rh2.0, 1-D plane-parallel
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Feb 24 15:30:40 2011 --

	   --------------------------                      ----------RH-- */

/* --- Solve for hydrostatic equilibrium --            -------------- */

#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "background.h"
#include "constant.h"
#include "accelerate.h"
#include "error.h"
#include "statistics.h"
#include "state.h"

/* --- Acceleration parameters --                      -------------- */

#define NG_HSE_DELAY  0
#define NG_HSE_ORDER  0
#define NG_HSE_PERIOD 0

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- Hydrostatic.c ----------- */

static void solveHydrostatic(RHState *rh, int NmaxIter, double iterLimit)
{
	const char   routineName[] = "Hydrostatic";
	register int k, n, m, i;

	char    label[MAX_MESSAGE_LENGTH];
	bool_t  Debeye, quiet;
	int     niter, Ngdelay, Ngperiod, Ngorder, Nhse;
	double *n_k, **dfdn, C1, C2, C3, beta, Phi_H, Phi_H2, *F, *dFdne, *f,
		dnmax = 0.0, nHtot_old, *np, *nH2, fH;
	struct Ng *Nghse;
	Atom      *atom;

	C1 = AMU / (2.0 * KBOLTZMANN);
	C2 = rh->atmos.gravity / KBOLTZMANN;
	if (rh->atmos.Stokes) C3 = 1.0 / (2.0 * MU_0 * KBOLTZMANN);

	F     = (double *) malloc(rh->atmos.Nspace * sizeof(double));
	dFdne = (double *) malloc(rh->atmos.Nspace * sizeof(double));

	/* --- Get the electron fraction F and its numerical derivative.
		   ne^{met} = F(T, ne) n_H^{tot} --              -------------- */

	FMetals(rh, F);
	dFMetals(rh, dFdne);

	/* --- Go through all depth points and iterate on the number
		   densities of nH, ne, np, nH2, and nHtot --    -------------- */

	np  = rh->atmos.H->n[rh->atmos.H->Nlevel - 1];
	nH2 = rh->atmos.H2->n;

	for (k = 0; k < rh->atmos.Nspace; k++) {

		/* --- H2 is part of the system at all temperatures: outside its
			   fit range the equilibrium constant is extrapolated, so nH2
			   is small but nonzero (N-53) --              -------------- */

		Nhse = 5;
		f    = (double *) malloc(Nhse * sizeof(double));
		n_k  = (double *) malloc(Nhse * sizeof(double));
		dfdn = matrix_double(rh, Nhse, Nhse);

		/* --- Get the equilibrium constants for hydrogen ionization and
			   H2 association and dissociation. This can be done because we
			   keep temperature and radiation field constant during the
			   hydrostatic equilibrium iterations --       -------------- */

		beta = (rh->atmos.totalAbund - 1.0) +
			   rh->atmos.wght_per_H * C1 * SQ(rh->atmos.vturb[k]) / rh->atmos.T[k];

		/* --- Starting solution:
			   Note: n_k[0] is the amount of neutral atomic hydrogen. -- */

		n_k[0]        = rh->atmos.H->ntotal[k] - np[k];
		n_k[1]        = rh->atmos.ne[k];
		n_k[2]        = np[k];
		n_k[3]        = nH2[k];
		n_k[Nhse - 1] = rh->atmos.nHtot[k];

		Phi_H = n_k[0] / (n_k[1] * n_k[2]);

		/* --- H2 forms from neutral atomic hydrogen, as in the chemical
			   equilibrium: n_H2 = Phi_H2 n_H^2 (N-44) --  -------------- */

		Phi_H2 = n_k[3] / SQ(n_k[0]);

		/* --- Fraction fH of nHtot in atomic hydrogen, protons and H2.
			   Each H2 molecule holds two hydrogen nuclei. The rest (H-
			   and hydrogen in other molecules) is kept at the fraction
			   of nHtot that the chemical equilibrium found, as F keeps
			   the electrons of the metals (N-44) --       -------------- */

		fH = (n_k[0] + n_k[2] + 2.0 * nH2[k]) / n_k[Nhse - 1];

		/* --- Initialize structure for Ng convergence acceleration -- -- */

		Nghse = NgInit(rh, Nhse, Ngdelay = NG_HSE_DELAY,
					   Ngorder = NG_HSE_ORDER, Ngperiod = NG_HSE_PERIOD, n_k);

		niter = 1;
		while (niter <= NmaxIter) {

			for (n = 0; n < Nhse; n++) {
				f[n] = 0.0;
				for (m = 0; m < Nhse; m++) dfdn[n][m] = 0.0;
			}
			/* --- Row for mass conservation --              -------------- */

			for (m = 0; m < Nhse - 1; m++) {
				f[0] += n_k[m];
				dfdn[0][m] = 1.0;
			}
			f[0] += beta * n_k[Nhse - 1];
			dfdn[0][Nhse - 1] = beta;

			/* --- Row for charge conservation --            -------------- */

			f[1]              = n_k[1] - n_k[2] - F[k] * n_k[Nhse - 1];
			dfdn[1][1]        = 1.0 - n_k[Nhse - 1] * dFdne[k];
			dfdn[1][2]        = -1.0;
			dfdn[1][Nhse - 1] = -F[k];

			/* --- Row for hydrogen number conservation --   -------------- */

			f[2]       = n_k[0] + n_k[2] - fH * n_k[Nhse - 1];
			dfdn[2][2] = dfdn[2][0] = 1.0;
			dfdn[2][Nhse - 1]       = -fH;

			f[2] += 2.0 * n_k[3];
			dfdn[2][3] = 2.0;

			/* --- Row for H2 equilibrium --                 -------------- */

			f[3]       = -SQ(n_k[0]) * Phi_H2 + n_k[3];
			dfdn[3][0] = -2.0 * n_k[0] * Phi_H2;
			dfdn[3][3] = 1.0;

			/* --- Row for hydrogen ionization --            -------------- */

			f[Nhse - 1]       = n_k[0] - n_k[1] * n_k[2] * Phi_H;
			dfdn[Nhse - 1][0] = 1.0;
			dfdn[Nhse - 1][1] = -n_k[2] * Phi_H;
			dfdn[Nhse - 1][2] = -n_k[1] * Phi_H;

			/* --- Fill right hand side --                   -------------- */

			f[0] -= C2 * rh->geometry->cmass[k] / rh->atmos.T[k];
			if (rh->atmos.Stokes) f[0] += C3 * SQ(rh->atmos.B[k]) / rh->atmos.T[k];

			/* --- Solve linearized set --                   -------------- */

			SolveLinearEq(rh, Nhse, dfdn, f, TRUE);
			for (n = 0; n < Nhse; n++) n_k[n] -= f[n];

			/* --- Check convergence and accelerate if appropriate -- ----- */

			Accelerate(rh, Nghse, n_k);
			snprintf(label, sizeof(label),
					 "\n%s-- Hydrostatic equilibrium: depth %3.3d, iteration %d",
					 (niter == 1) ? "\n" : "", k, niter);

			if ((dnmax = MaxChange(rh, Nghse, label, quiet = TRUE)) <= iterLimit)
				break;
			niter++;
		}
		if (dnmax > iterLimit) {
			Warning(rh, routineName,
					"Hydrostatic equilibrium iteration not converged:\n"
					" temperature: %6.1f [K], \n"
					" hydrogen density: %9.3E [m^-3],\n dnmax: %9.3E\n",
					rh->atmos.T[k], rh->atmos.nHtot[k], dnmax);
		}
		/* --- Store results --                            -------------- */

		nHtot_old = rh->atmos.H->ntotal[k];

		rh->atmos.H->ntotal[k] = n_k[0] + n_k[2];
		rh->atmos.ne[k]        = n_k[1];
		np[k]                  = n_k[2];
		rh->atmos.nHtot[k]     = n_k[Nhse - 1];

		nH2[k] = n_k[3];

		/* --- Adjust Non-LTE hydrogen level popultions -- -------------- */

		for (i = 0; i < rh->atmos.H->Nlevel; i++)
			rh->atmos.H->n[i][k] *= rh->atmos.H->ntotal[k] / nHtot_old;

		NgFree(rh, Nghse);
		free(f);
		free(n_k);
		freeMatrix(rh, (void **) dfdn);
	}
	/* --- Adjust total populations of background metals and recalculate
		   all the LTE populations, including hydrogen -- ------------- */

	LTEpops(rh, rh->atmos.H, Debeye = TRUE);

	for (n = 1; n < rh->atmos.Natom; n++) {
		atom = &rh->atmos.atoms[n];
		for (k = 0; k < rh->atmos.Nspace; k++)
			atom->ntotal[k] = atom->abundance * rh->atmos.nHtot[k];

		LTEpops(rh, atom, Debeye = TRUE);
	}

	/* --- Clean up --                                   -------------- */

	free(F);
	free(dFdne);
}

/* --- Hydrostatic() with its calls and wall time counted for the work
	   statistics (S-1, SOLVER.md section 4.1) --       -------------- */

void Hydrostatic(RHState *rh, int NmaxIter, double iterLimit)
{
	double t_start = wallClock();

	solveHydrostatic(rh, NmaxIter, iterLimit);
	workAdd(rh, WORK_HYDROSTATIC, wallClock() - t_start);
}
/* ------- end ---------------------------- Hydrostatic.c ----------- */
