/* ------- file: -------------------------- chemequil.c -------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Feb  6 10:33:22 2018 --

	   --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "background.h"
#include "constant.h"
#include "error.h"
#include "statistics.h"
#include "inputs.h"
#include "negtau.h"
#include "state.h"

#define COMMENT_CHAR "#"

/* --- Solver of the equilibrium (N-61, see below): the largest change
	   of a log density in one Newton step, the relative residual below
	   which full Newton steps are taken without a line search, and the
	   regularization of the scaled Hessian --        -------------- */

#define CHEM_STEP_MAX   30.0
#define CHEM_FULL_STEPS 1.0E-3
#define CHEM_REGULARIZE 1.0E-12

/* --- Evaluate chemical equilibrium for set of molecules and the complete
	   and sufficient set of their constituent nuclei.

	   Equations to solve:

	1) Conservation of number density for each constituent:

		n_i + Sum_m {N_i^m * n_m} = A_i * nHtot,

	   n_i are the atomic population number densities (all atoms and ions
	   of one species not bound in molecules), abund_i the element abundance,
	   n_m the molecular number densities of molecules containing element i,
	   and N_i^m is the number of nuclei of element i in molecule m.


	2) Chemical equilibrium for each molecule (Saha):

		n_m = Prod_i {(f0_i * n_i)^N_i^m} * Phi_m(T),

		where f0_i is the fraction of atoms i in the neutral stage,
		and Phi_m(T) is the equilibrium constant for molecule m.
		If molecule n_m has a charge of +1, one of the f0_i is replaced
		by f1_i, the fraction of element i in the first ionization stage.

  Note: Molecules with other charge than 0 or +1 are currently not allowed
		and are rejected in routine readMolecule.c


	-- The Hminus population number is added in the hydrogen conservatiom
	   equation, which is always the first equation.
	   Hmin formation is given by:

		nHmin = ne * nH * PhiHmin,

	   where PhiHmin = 1/4 * (h^2/(2PI m_e kT))^3/2 exp(Ediss/kT)


 Note: The total hydrogen number density (including protons, H-, and
	   nuclei that are part of molecules like H2 and H2+ is stored in
	   atmos->nHtot. The total number not in molecules or H- (i.e. the
	   neutral atoms plus the protons is stored in atmos->H.ntotal.

 Solution (N-61): the unknowns are the logarithms u_i = ln n_i of the
	   densities of the free nuclei only. Every other species s (H-,
	   each molecule) follows from them by the Saha equation, and so does
	   every free nucleus itself:

		n_s(u) = c_s exp(Sum_i nu_si u_i),

	   with nu_si the number of nuclei i in s, c_s = 1 for a free nucleus,
	   c = ne fn0_H PhiHmin for H-, and c_m = Phi_m Prod_i fn0_i^nu_mi /
	   ne^charge for molecule m. The conservation equations are then the
	   gradient of the strictly convex function

		F(u) = Sum_s n_s(u) - Sum_i a_i u_i,       a_i = A_i nHtot,

	   whose Hessian Sum_s n_s nu_s nu_s^T is positive definite (every
	   free nucleus contributes n_i e_i e_i^T). The equilibrium is the
	   unique minimum of F, and Newton's method with a line search on F
	   converges from any starting point; the densities are positive by
	   construction. The start is everything dissociated (u_i = ln a_i);
	   far from the solution the step is capped, shortened until F
	   decreases enough (Armijo) or lengthened while F keeps decreasing
	   (from the dissociated start at low temperature the Newton step on
	   an exponential is short). Once the residual max |g_i| / a_i is
	   below iterLimit, one more full Newton step is taken, so that the
	   result is converged to rounding whatever the number of steps (the
	   map G stays a function of its state, S-2b). The old Newton
	   iteration on the densities themselves (no damping, 10 steps, only
	   a warning) left negative molecular densities below about 2000 K.
	   A point that does not converge in NmaxIter steps invalidates a
	   trial point of the solvers (negTauInvalid) and stops the run
	   otherwise.
	   --                                              -------------- */

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

double chemLastChangeReset(RHState *rh)
{
	double d = rh->chemequil.chemLastChange;

	rh->chemequil.chemLastChange = 0.0;
	return d;
}

/* ------- begin -------------------------- ChemicalEquilibrium.c --- */

static void chemicalEquilibrium(RHState *rh, int NmaxIter, double iterLimit,
								bool_t ofState);

void ChemicalEquilibrium(RHState *rh, int NmaxIter, double iterLimit)
{
	chemicalEquilibrium(rh, NmaxIter, iterLimit, FALSE);
}

/* --- For the map G with ne in the state (S-2b, docs/SOLVER.md): a
	   function of the state only. The neutral fraction of an active
	   atom comes from its populations (normalized by their sum, not by
	   ntotal of an earlier call), the equilibrium is converged to
	   rounding (N-61), and the populations of active atoms are not
	   rescaled in place: only ntotal and the LTE populations follow the
	   new equilibrium --                               -------------- */

void ChemicalEquilibriumOfState(RHState *rh, int NmaxIter, double iterLimit)
{
	chemicalEquilibrium(rh, NmaxIter, iterLimit, TRUE);
}

/* ------- begin -------------------------- chemSolve.c ------------- */

/* --- The equilibrium of one location in log densities of the free
	   nuclei (N-61, see the header). Species s = 0 .. Ns-1 with
	   stoichiometry nu[s][0 .. Nn-1] and ln c_s in lnc[s] (-HUGE_VAL
	   for an absent species: an equilibrium constant or neutral fraction
	   of 0). On return u holds the solution, *res the relative residual
	   before the last step and *step the largest change of a log
	   density in it; returns the number of steps, or -1 if not
	   converged in maxIter --                         -------------- */

static double chemF(ChemSystem *cs, const double *a, const double *u)
{
	int    s, i;
	double F = 0.0, e;

	for (s = 0; s < cs->Ns; s++) {
		if (cs->lnc[s] == -HUGE_VAL) continue;
		e = cs->lnc[s];
		for (i = 0; i < cs->Nn; i++) e += cs->nu[s][i] * u[i];
		F += exp(e);
	}
	for (i = 0; i < cs->Nn; i++) F -= a[i] * u[i];
	return isfinite(F) ? F : HUGE_VAL;
}

int chemSolve(ChemSystem *cs, const double *a, double *u, int maxIter,
			  double tol, double *res, double *step)
{
	int    Nn = cs->Nn, s, i, j, l, it;
	double r, e, dmax, t, tmax, F0, F1, F2, slope, sum;

	for (it = 0; it <= maxIter; it++) {

		/* --- Densities, gradient (the conservation residual) and the
			   Hessian of F --                            -------------- */

		for (i = 0; i < Nn; i++) {
			cs->g[i] = -a[i];
			for (j = 0; j < Nn; j++) cs->H[i][j] = 0.0;
		}
		for (s = 0; s < cs->Ns; s++) {
			if (cs->lnc[s] == -HUGE_VAL) {
				cs->n[s] = 0.0;
				continue;
			}
			e = cs->lnc[s];
			for (i = 0; i < Nn; i++) e += cs->nu[s][i] * u[i];
			cs->n[s] = exp(e);
			for (i = 0; i < Nn; i++) {
				if (cs->nu[s][i] == 0.0) continue;
				cs->g[i] += cs->nu[s][i] * cs->n[s];
				for (j = 0; j < Nn; j++) cs->H[i][j] += cs->n[s] * cs->nu[s][i] * cs->nu[s][j];
			}
		}
		r = 0.0;
		for (i = 0; i < Nn; i++) r = MAX(r, fabs(cs->g[i]) / a[i]);
		if (!isfinite(r)) return -1;
		*res = r;
		if (it == maxIter && r >= tol) return -1;

		/* --- Newton step: Cholesky of the Hessian scaled to unit
			   diagonal, with a small regularization (one species can
			   dominate several nuclei far from the solution) -- ---- */

		for (i = 0; i < Nn; i++) cs->sc[i] = 1.0 / sqrt(MAX(cs->H[i][i], 1.0E-300));
		for (i = 0; i < Nn; i++) {
			for (j = 0; j < Nn; j++) cs->H[i][j] *= cs->sc[i] * cs->sc[j];
			cs->H[i][i] += CHEM_REGULARIZE;
		}
		for (j = 0; j < Nn; j++) {
			sum = cs->H[j][j];
			for (l = 0; l < j; l++) sum -= SQ(cs->H[j][l]);
			if (!(sum > 0.0)) return -1;
			cs->H[j][j] = sqrt(sum);
			for (i = j + 1; i < Nn; i++) {
				sum = cs->H[i][j];
				for (l = 0; l < j; l++) sum -= cs->H[i][l] * cs->H[j][l];
				cs->H[i][j] = sum / cs->H[j][j];
			}
		}
		for (i = 0; i < Nn; i++) {
			sum = -cs->g[i] * cs->sc[i];
			for (l = 0; l < i; l++) sum -= cs->H[i][l] * cs->du[l];
			cs->du[i] = sum / cs->H[i][i];
		}
		for (i = Nn - 1; i >= 0; i--) {
			sum = cs->du[i];
			for (l = i + 1; l < Nn; l++) sum -= cs->H[l][i] * cs->du[l];
			cs->du[i] = sum / cs->H[i][i];
		}
		dmax = 0.0;
		for (i = 0; i < Nn; i++) {
			cs->du[i] *= cs->sc[i];
			dmax = MAX(dmax, fabs(cs->du[i]));
		}

		/* --- Converged: one more full step takes the solution to
			   rounding --                                -------------- */

		if (r < tol || r < CHEM_FULL_STEPS) {
			for (i = 0; i < Nn; i++) u[i] += cs->du[i];
			*step = dmax;
			if (r < tol) {
				for (s = 0; s < cs->Ns; s++) {
					if (cs->lnc[s] == -HUGE_VAL) continue;
					e = cs->lnc[s];
					for (i = 0; i < Nn; i++) e += cs->nu[s][i] * u[i];
					cs->n[s] = exp(e);
				}
				return it + 1;
			}
			continue;
		}

		/* --- Far from the solution: capped step, Armijo backtracking
			   on F, or doubling while F decreases -- ------------ */

		slope = 0.0;
		for (i = 0; i < Nn; i++) slope += cs->g[i] * cs->du[i];
		tmax = CHEM_STEP_MAX / dmax;
		t    = MIN(1.0, tmax);
		F0   = chemF(cs, a, u);
		for (i = 0; i < Nn; i++) cs->ut[i] = u[i] + t * cs->du[i];
		F1 = chemF(cs, a, cs->ut);
		if (F1 <= F0 + 1.0E-4 * t * slope) {
			while (2.0 * t <= tmax) {
				for (i = 0; i < Nn; i++) cs->ut[i] = u[i] + 2.0 * t * cs->du[i];
				F2 = chemF(cs, a, cs->ut);
				if (!(F2 < F1)) break;
				t *= 2.0;
				F1 = F2;
			}
		} else {
			while (F1 > F0 + 1.0E-4 * t * slope && t > 1.0E-12) {
				t *= 0.5;
				for (i = 0; i < Nn; i++) cs->ut[i] = u[i] + t * cs->du[i];
				F1 = chemF(cs, a, cs->ut);
			}
		}
		for (i = 0; i < Nn; i++) u[i] += t * cs->du[i];
		*step = t * dmax;
	}
	return -1;
}
/* ------- end ---------------------------- chemSolve.c ------------- */

static void chemicalEquilibrium(RHState *rh, int NmaxIter, double iterLimit, bool_t ofState)
{
	const char   routineName[] = "ChemicalEquilibrium";
	register int k, i, j, nu, pass;

	char    tmpStr[13], msg[MAX_MESSAGE_LENGTH];
	int     Nequation, **nucl_index, niter, Nnuclei, Nmaxstage, s;
	double *a, *n, *u, *Phi, PhiHmin, fHmin, CI, res = 0.0, step = 0.0, sum, lnc,
												 *fn0, fraction, *fjk = NULL, *dfjk = NULL;
	ChemSystem cs;
	Atom      *atom;
	Molecule  *molecule;
	Element  **nuclei;

	getCPU(rh, 3, TIME_START, NULL);

	/* --- Constant for Saha equation Hminus --          -------------- */

	CI = (HPLANCK / (2.0 * PI * M_ELECTRON)) * (HPLANCK / KBOLTZMANN);

	/* --- Solve the chemical equilibrium equations. First,
		   determine for what atoms and molecules the equilibrium
		   equations have to be solved. --               -------------- */

	/* --- Collect pointers to elements that can be bound in molecules. */

	Nnuclei = 0;
	nuclei  = (Element **) mallocOrAbort(rh->atmos.Nelem * sizeof(Element *));
	for (i = 0; i < rh->atmos.Nelem; i++) {
		if (rh->atmos.elements[i].Nmolecule > 0) {
			nuclei[Nnuclei++] = &rh->atmos.elements[i];
		}
	}
	nuclei = (Element **) reallocOrAbort(nuclei, Nnuclei * sizeof(Element *));

	/* --- Check that first Nucleus is hydrogen --       -------------- */

	if (!strstr(nuclei[0]->ID, "H ")) {
		Error(rh, ERROR_LEVEL_2, routineName, "First nucleus must be H not %s "
											  "(check H2.molecule)",
			  nuclei[0]->ID);
	}

	Nmaxstage = 0;
	for (j = 0; j < Nnuclei; j++) {
		if (nuclei[j]->model != NULL) {
			if (nuclei[j]->model->stage[0] > 0) {
				Error(rh, ERROR_LEVEL_2, routineName, "Model for element %s does not have a neutral stage\n"
													  " needed for molecular formation\n",
					  nuclei[j]->ID);
			}
		} else
			Nmaxstage = MAX(Nmaxstage, nuclei[j]->Nstage);
	}
	if (Nmaxstage) {
		fjk  = (double *) mallocOrAbort(Nmaxstage * sizeof(double));
		dfjk = (double *) mallocOrAbort(Nmaxstage * sizeof(double));
	}

	/* --- Quantity nucl_index[i][j] stores the index (in array nuclei)
		   of the jth element bound in molecule i --     -------------- */

	nucl_index = (int **) mallocOrAbort(rh->atmos.Nmolecule * sizeof(int *));
	for (i = 0; i < rh->atmos.Nmolecule; i++) {
		molecule      = &rh->atmos.molecules[i];
		nucl_index[i] = (int *) mallocOrAbort(molecule->Nelement * sizeof(int));
		for (j = 0; j < molecule->Nelement; j++) {
			for (nu = 0; nu < Nnuclei; nu++) {
				if (nuclei[nu] == &rh->atmos.elements[molecule->pt_index[j]]) {
					nucl_index[i][j] = nu;
					break;
				}
			}
		}
	}
	/* --- Number of equations --                        -------------- */

	Nequation = Nnuclei + rh->atmos.Nmolecule;

	/* --- Allocate temporary storage space.
		   Quantities for the Newton Raphson --          -------------- */

	n = (double *) callocOrAbort(Nequation, sizeof(double));
	a = (double *) callocOrAbort(Nequation, sizeof(double));
	u = (double *) callocOrAbort(Nnuclei, sizeof(double));

	/* --- Species of the log-density solver: the free nuclei, H- and
		   the molecules, with their stoichiometry (N-61) -- -------- */

	cs.Nn  = Nnuclei;
	cs.Ns  = Nnuclei + 1 + rh->atmos.Nmolecule;
	cs.nu  = matrix_double(rh, cs.Ns, Nnuclei);
	cs.H   = matrix_double(rh, Nnuclei, Nnuclei);
	cs.lnc = (double *) mallocOrAbort(cs.Ns * sizeof(double));
	cs.n   = (double *) mallocOrAbort(cs.Ns * sizeof(double));
	cs.g   = (double *) mallocOrAbort(Nnuclei * sizeof(double));
	cs.du  = (double *) mallocOrAbort(Nnuclei * sizeof(double));
	cs.ut  = (double *) mallocOrAbort(Nnuclei * sizeof(double));
	cs.sc  = (double *) mallocOrAbort(Nnuclei * sizeof(double));
	for (s = 0; s < cs.Ns; s++)
		for (i = 0; i < Nnuclei; i++) cs.nu[s][i] = 0.0;
	for (i = 0; i < Nnuclei; i++) cs.nu[i][i] = 1.0;
	cs.nu[Nnuclei][0] = 1.0; /* H- */
	for (i = 0; i < rh->atmos.Nmolecule; i++) {
		molecule = &rh->atmos.molecules[i];
		for (j = 0; j < molecule->Nelement; j++)
			cs.nu[Nnuclei + 1 + i][nucl_index[i][j]] += molecule->pt_count[j];
	}

	/* --- Equilibrium constant for each molecule, and neutral
		   fraction for each nucleus --                  -------------- */

	Phi = (double *) mallocOrAbort(rh->atmos.Nmolecule * sizeof(double));
	fn0 = (double *) mallocOrAbort(Nnuclei * sizeof(double));

	/* --- Go through spatial grid and solve (local) equations -- ----- */

	for (k = 0; k < rh->atmos.Nspace; k++) {

		/* --- Collect for each atom the population fraction of
			   the neutral stage --                        -------------- */

		for (i = 0; i < Nnuclei; i++) {
			if ((atom = nuclei[i]->model) == NULL) {

				/* --- If no atomic model is present for this nucleus -- ---- */

				getfjk(rh, nuclei[i], rh->atmos.ne[k], k, fjk, dfjk);
				fn0[i] = fjk[0];
				a[i]   = nuclei[i]->abund * rh->atmos.nHtot[k];
			} else {

				/* --- Atomic model has been read for this nucleus -- ------- */

				/* --- N-36: the neutral fraction of an atom comes from its
					   populations n, normalized by their sum, as in the map
					   G. An active atom whose populations have not been
					   initialized yet (atom->popsInit, set by initSolution;
					   the setup Background() runs before it, with n still
					   zero or the atmosphere-file populations of H) takes the
					   LTE populations nstar, as a fresh run always did. A sum
					   that is zero or not finite also falls back to nstar, so
					   that fn0 is never NaN. The legacy iteration used to take
					   nstar for an active atom unless OLD_POPULATIONS, so its
					   SOLVE_NE = ITERATION fixed point differed from that of
					   the G solvers (H2, H-, ne by 1e-6) and a restart
					   changed the result -- --------------------------- */

				double **pop = atom->n;
				if (atom->active && !ofState && !atom->popsInit)
					pop = atom->nstar;
				for (pass = 0; pass < 2; pass++) {
					fn0[i] = sum = 0.0;
					for (j = 0; j < atom->Nlevel; j++) {
						if (atom->stage[j] == 0) fn0[i] += pop[j][k];
						sum += pop[j][k];
					}
					if (sum > 0.0 && isfinite(sum)) break;
					pop = atom->nstar;
				}
				fn0[i] /= sum;
				a[i] = atom->abundance * rh->atmos.nHtot[k];
			}
		}
		PhiHmin = 0.25 * pow(CI / rh->atmos.T[k], 1.5) *
				  exp(E_ION_HMIN / (KBOLTZMANN * rh->atmos.T[k]));
		fHmin = rh->atmos.ne[k] * fn0[0] * PhiHmin;

		/* --- Equilibrium constant for each molecule at this location -- */

		for (i = 0; i < rh->atmos.Nmolecule; i++)
			Phi[i] = equilconstant(rh, &rh->atmos.molecules[i], rh->atmos.T[k]);

		/* --- ln c of every species at this location: 0 for a free
			   nucleus, H- from the neutral hydrogen density, a molecule
			   from its equilibrium constant, the neutral fractions of
			   its nuclei and its charge (N-61) --       -------------- */

		for (i = 0; i < Nnuclei; i++) cs.lnc[i] = 0.0;
		cs.lnc[Nnuclei] = (fHmin > 0.0) ? log(fHmin) : -HUGE_VAL;
		for (i = 0; i < rh->atmos.Nmolecule; i++) {
			molecule = &rh->atmos.molecules[i];
			lnc      = (Phi[i] > 0.0) ? log(Phi[i]) : -HUGE_VAL;
			if (molecule->charge != 0 && lnc != -HUGE_VAL)
				lnc = (rh->atmos.ne[k] > 0.0) ? lnc - molecule->charge * log(rh->atmos.ne[k]) : -HUGE_VAL;
			for (j = 0; j < molecule->Nelement && lnc != -HUGE_VAL; j++) {
				nu  = nucl_index[i][j];
				lnc = (fn0[nu] > 0.0) ? lnc + molecule->pt_count[j] * log(fn0[nu])
									  : -HUGE_VAL;
			}
			cs.lnc[Nnuclei + 1 + i] = lnc;
		}

		/* --- Start from everything dissociated, solve --  ------------ */

		for (i = 0; i < Nnuclei; i++) u[i] = log(a[i]);
		niter = chemSolve(&cs, a, u, NmaxIter, iterLimit, &res, &step);
		if (niter < 0) {
			snprintf(msg, sizeof(msg),
					 "\n  Chemical equilibrium did not converge in %d Newton"
					 " steps at location %d: T = %.1f K, nHtot = %.4E m^-3,"
					 " ne = %.4E m^-3, relative residual %.3E",
					 NmaxIter, k, rh->atmos.T[k], rh->atmos.nHtot[k], rh->atmos.ne[k], res);
			if (negTauInTrial(rh)) {
				negTauInvalid(rh, msg + 3);
			} else
				Error(rh, ERROR_LEVEL_2, routineName, "%s", msg);
		}
		for (i = 0; i < Nnuclei; i++) n[i] = cs.n[i];
		for (i = 0; i < rh->atmos.Nmolecule; i++) n[Nnuclei + i] = cs.n[Nnuclei + 1 + i];
		if (ofState) rh->chemequil.chemLastChange = MAX(rh->chemequil.chemLastChange, step);

		/* --- Store population numbers nuclei --          -------------- */

		for (i = 0; i < Nnuclei; i++) {
			if ((atom = nuclei[i]->model) != NULL) {
				fraction = n[i] / atom->ntotal[k];

				for (j = 0; j < atom->Nlevel; j++) {
					atom->nstar[j][k] *= fraction;
					if (atom->n != atom->nstar && !ofState)
						atom->n[j][k] *= fraction;
				}
				atom->ntotal[k] = n[i];
			}
		}
		/* --- Store Hmin density, from neutral hydrogen as in the
			   conservation equation above (B-40) --       -------------- */

		rh->atmos.nHmin[k] = fHmin * n[0];

		/* --- Store molecular densities --                -------------- */

		for (i = 0; i < rh->atmos.Nmolecule; i++)
			rh->atmos.molecules[i].n[k] = n[Nnuclei + i];
	}

	/* --- Check whether active atom, if present, is in list of nuclei.
		   If so print out a warning --                  -------------- */

	for (nu = 0; nu < Nnuclei; nu++) {
		atom = nuclei[nu]->model;
		if (atom && atom->active) {
			snprintf(msg, sizeof(msg), "\nReduced number density of"
									   " active atom %s due to molecule%s\n  ",
					 atom->ID, (nuclei[nu]->Nmolecule > 1) ? "s" : "");
			for (i = 0; i < nuclei[nu]->Nmolecule; i++) {
				sprintf(tmpStr, "%s%s",
						rh->atmos.molecules[nuclei[nu]->mol_index[i]].ID,
						(i == nuclei[nu]->Nmolecule - 1) ? "\n\n" : ", ");
				snprintf(msg + strlen(msg), sizeof(msg) - strlen(msg), "%s", tmpStr);
			}
			Message(rh, "%s", msg);
		}
	}

	/* --- Clean up --                                   -------------- */

	free(n);
	free(a);
	free(u);
	free(Phi);
	free(fn0);
	freeMatrix(rh, (void **) cs.nu);
	freeMatrix(rh, (void **) cs.H);
	free(cs.lnc);
	free(cs.n);
	free(cs.g);
	free(cs.du);
	free(cs.ut);
	free(cs.sc);

	free(nuclei);
	for (i = 0; i < rh->atmos.Nmolecule; i++) free(nucl_index[i]);
	free(nucl_index);

	if (Nmaxstage) {
		free(fjk);
		free(dfjk);
	}
	getCPU(rh, 3, TIME_POLL, "Chemical equilibrium");
}
/* ------- end ---------------------------- ChemicalEquilibrium.c --- */

/* ------- begin -------------------------- partfunction.c ---------- */

double partfunction(RHState *rh, struct Molecule *molecule, double T)
{
	register int i;

	double pf = 0.0, t;

	/* --- Outside the range of the fit hold the partition function at
		   its value at the nearest edge of [Tmin, Tmax] (N-53). The
		   caller reports the number of depth points outside the range
		   once per molecule (checkMoleculeTrange) --  ---------------- */

	if (T < molecule->Tmin)
		T = molecule->Tmin;
	else if (T > molecule->Tmax)
		T = molecule->Tmax;

	/* --- Evaluate polynomial in temperature for partition function -- */

	switch (molecule->fit) {
	case KURUCZ_70:
		pf = molecule->pf_coef[0];
		for (i = 1; i < molecule->Npf; i++)
			pf = pf * T + molecule->pf_coef[i];

		pf = exp(pf);
		break;

	case KURUCZ_85:
		t  = T * 1.0E-4;
		pf = molecule->pf_coef[0];
		for (i = 1; i < molecule->Npf; i++)
			pf = pf * t + molecule->pf_coef[i];

		pf = exp(pf);
		break;

	case SAUVAL_TATUM_84:
		t  = log10(THETA0 / T);
		pf = molecule->pf_coef[0];
		for (i = 1; i < molecule->Npf; i++)
			pf = pf * t + molecule->pf_coef[i];

		pf = POW10(pf);
		break;

	case IRWIN_81:
		t  = log(T);
		pf = molecule->pf_coef[0];
		for (i = 1; i < molecule->Npf; i++)
			pf = pf * t + molecule->pf_coef[i];

		pf = exp(pf);
		break;

	case TSUJI_73:
		break;
	default:
		Error(rh, ERROR_LEVEL_2, "partfunction", "Unknown method for calculation of partition function\n"
												 "for molecule %s: %d",
			  molecule->ID, molecule->fit);
	}
	return pf;
}
/* ------- end ---------------------------- partfunction.c ---------- */

/* ------- begin -------------------------- equilconstant.c --------- */

/* --- ln of the equilibrium constant [SI] from the fit, evaluated
	   without regard to the range [Tmin, Tmax] of the fit. Used only for
	   the extrapolation outside that range, see equilconstant. -- ----- */

static double lnEquilConstant(RHState *rh, struct Molecule *molecule, double T)
{
	register int i;

	int    mk;
	double kT, t, lneqc = 0.0, theta;

	switch (molecule->fit) {
	case KURUCZ_70:
	case KURUCZ_85:
		t  = (molecule->fit == KURUCZ_85) ? T * 1.0E-4 : T;
		kT = KBOLTZMANN * T;
		mk = molecule->Nnuclei - 1 - molecule->charge;

		lneqc = molecule->eqc_coef[0];
		for (i = 1; i < molecule->Neqc; i++)
			lneqc = lneqc * t + molecule->eqc_coef[i];
		lneqc += molecule->Ediss / kT - 1.5 * mk * log(T) +
				 mk * 3.0 * log(CM_TO_M);
		break;

	case IRWIN_81:
	case SAUVAL_TATUM_84:
		theta = THETA0 / T;
		t     = log10(theta);
		kT    = KBOLTZMANN * T;

		lneqc = molecule->eqc_coef[0];
		for (i = 1; i < molecule->Neqc; i++)
			lneqc = lneqc * t + molecule->eqc_coef[i];
		lneqc = log(10.0) * ((molecule->Ediss / EV) * theta - lneqc) + log(kT);
		break;

	case TSUJI_73:
		theta = THETA0 / T;
		kT    = KBOLTZMANN * T;

		lneqc = molecule->eqc_coef[0];
		for (i = 1; i < molecule->Neqc; i++)
			lneqc = lneqc * theta + molecule->eqc_coef[i];
		lneqc = 2.0 * log(kT) - log(10.0) * lneqc +
				(molecule->Nnuclei - 1) * log(CUBE(CM_TO_M) / ERG_TO_JOULE);
		break;

	default:
		Error(rh, ERROR_LEVEL_2, "lnEquilConstant", "Unknown method for calculation of equilibrium constant\n"
													"for molecule %s: %d",
			  molecule->ID, molecule->fit);
	}
	return lneqc;
}
/* ------- end ---------------------------- lnEquilConstant -------- */

/* --- Outside [Tmin, Tmax] ln K is extrapolated linearly in
	   theta = 5040/T from the nearer edge, using ln K and d lnK / d theta
	   at the edge (van 't Hoff: ln K is close to linear in 1/T, with
	   the dissociation energy as slope). The slope is a centred
	   difference over +-0.1 percent of theta_edge: the fits are smooth
	   (polynomials in T, ln T or log theta), so the truncation error is
	   negligible and no per-fit analytic derivative is needed. The
	   result is capped at ln K = 700 to stay below the double range
	   (N-53) --                                     ----------------- */

static double extrapolatedEquilConstant(RHState *rh, struct Molecule *molecule, double T)
{
	double Tedge, theta, theta_edge, dtheta, lnK_edge, slope, lnK;

	Tedge      = (T < molecule->Tmin) ? molecule->Tmin : molecule->Tmax;
	theta      = THETA0 / T;
	theta_edge = THETA0 / Tedge;
	dtheta     = 1.0E-3 * theta_edge;

	lnK_edge = lnEquilConstant(rh, molecule, Tedge);
	slope    = (lnEquilConstant(rh, molecule, THETA0 / (theta_edge + dtheta)) -
             lnEquilConstant(rh, molecule, THETA0 / (theta_edge - dtheta))) /
			(2.0 * dtheta);

	lnK = lnK_edge + slope * (theta - theta_edge);
	return exp(MIN(lnK, 700.0));
}

double equilconstant(RHState *rh, struct Molecule *molecule, double T)
{
	register int i;

	int    mk;
	double kT, t, eqc = 0.0, theta, cgs_to_SI = 1.0;

	if (T < molecule->Tmin || T > molecule->Tmax)
		return extrapolatedEquilConstant(rh, molecule, T);

	/* --- Evaluate polynomial in temperature for equilibrium constant.
		   Constant should have units of [m^3,6,9, etc] -- ------------ */

	switch (molecule->fit) {
	case KURUCZ_70:
		kT = KBOLTZMANN * T;
		mk = molecule->Nnuclei - 1 - molecule->charge;

		eqc = molecule->eqc_coef[0];
		for (i = 1; i < molecule->Neqc; i++)
			eqc = eqc * T + molecule->eqc_coef[i];
		eqc = exp(molecule->Ediss / kT + eqc - 1.5 * mk * log(T));

		cgs_to_SI = pow(CUBE(CM_TO_M), mk);
		break;

	case KURUCZ_85:
		t  = T * 1.0E-4;
		kT = KBOLTZMANN * T;
		mk = molecule->Nnuclei - 1 - molecule->charge;

		eqc = molecule->eqc_coef[0];
		for (i = 1; i < molecule->Neqc; i++)
			eqc = eqc * t + molecule->eqc_coef[i];
		eqc = exp(molecule->Ediss / kT + eqc - 1.5 * mk * log(T));

		cgs_to_SI = pow(CUBE(CM_TO_M), mk);
		break;

	case IRWIN_81:

		/* ---- Use SAUVAL_TATUM_84 for equilibrium constant -- ---------- */
		;
	case SAUVAL_TATUM_84:
		theta = THETA0 / T;
		t     = log10(theta);
		kT    = KBOLTZMANN * T;

		eqc = molecule->eqc_coef[0];
		for (i = 1; i < molecule->Neqc; i++)
			eqc = eqc * t + molecule->eqc_coef[i];
		eqc = POW10((molecule->Ediss / EV) * theta - eqc) * kT;

		/* --- Constants are already given in SI units by Sauval & Tatum - */

		cgs_to_SI = 1.0;
		break;

	case TSUJI_73:
		theta = THETA0 / T;
		kT    = KBOLTZMANN * T;

		eqc = molecule->eqc_coef[0];
		for (i = 1; i < molecule->Neqc; i++)
			eqc = eqc * theta + molecule->eqc_coef[i];
		eqc = SQ(kT) * POW10(-eqc);

		cgs_to_SI = pow(CUBE(CM_TO_M) / ERG_TO_JOULE, molecule->Nnuclei - 1);
		break;

	default:
		Error(rh, ERROR_LEVEL_2, "equilconstant", "Unknown method for calculation of equilibrium constant\n"
												  "for molecule %s: %d",
			  molecule->ID, molecule->fit);
	}

	return eqc * cgs_to_SI;
}
/* ------- end ---------------------------- equilconstant.c --------- */
