/* ------- file: -------------------------- fillgamma.c -------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu May 21 17:09:54 2020 --

	   --------------------------                      ----------RH-- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "error.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "atom.h"
#include "inputs.h"
#include "constant.h"

/* --- Routines for wavelength- and angle-integrated contributions to
	   the crosscoupling, Gamma matrix, and radiative rates.

	   Convention: \Gamma_ij = Gamma[i][j] represents the
				   transition j --> i, so that \Gamma_ij * n_j
		   is the rate per sec out of level j to level i.
	   --                                              -------------- */

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum   spectrum;
extern InputData  input;
extern char       messageStr[];

/* ------- begin -------------------------- initGammaAtom.c --------- */

void initGammaAtom(Atom *atom, int niter)
{
	register int ij, k;

	double CR_factor;

	/* --- Add the fixed rates into Gamma --             -------------- */

	for (ij = 0; ij < SQ(atom->Nlevel); ij++) {
		for (k = 0; k < atmos.Nspace; k++)
			atom->Gamma[ij][k] = atom->C[ij][k];
	}
	/* --- Apply Collisional Relaxation (CR) factor --   -------------- */

	if (input.CR_Nstep > 0 && niter <= input.CR_Nstep) {
		CR_factor = pow(input.CR_factor,
						(1.0 - (double) input.CR_Nstep / niter));

		for (ij = 0; ij < SQ(atom->Nlevel); ij++) {
			for (k = 0; k < atmos.Nspace; k++)
				atom->Gamma[ij][k] *= CR_factor;
		}
	}
}
/* ------- end ---------------------------- initGammaAtom.c --------- */

/* ------- begin -------------------------- initGammaMolecule.c ----- */

void initGammaMolecule(Molecule *molecule)
{
	register int ij, ji, k, vi, vj;

	/* --- Add the fixed rates into Gamma --             -------------- */

	for (vi = 0; vi < molecule->Nv - 1; vi++) {
		vj = vi + 1;
		ij = vi * molecule->Nv + vj;
		ji = vj * molecule->Nv + vi;
		for (k = 0; k < atmos.Nspace; k++) {
			if (molecule->n[k]) {
				molecule->Gamma[ij][k] = molecule->C_ul[k];
				molecule->Gamma[ji][k] = molecule->C_ul[k] *
										 molecule->nvstar[vj][k] / molecule->nvstar[vi][k];
			}
		}
	}
}
/* ------- end ---------------------------- initGammaMolecule.c ----- */

/* ------- begin -------------------------- zeroGammaThreads.c ------ */

/* --- Zero every thread's private Gamma accumulator at the start of
	   an iteration. On the first iteration a thread's rhth[nt].Gamma
	   may not be allocated yet (lazily allocated in alloc_as() on
	   first use for this atom) -- skip it; a freshly calloc'd buffer
	   already starts at zero, so there is nothing to clear. -- ---- */

void zeroGammaThreads(Atom *atom)
{
	register int nt, ij, k;

	if (input.Nthreads <= 1) return; /* one thread uses atom->Gamma */

	for (nt = 0; nt < input.Nthreads; nt++) {
		if (atom->rhth[nt].Gamma == NULL) continue;
		for (ij = 0; ij < SQ(atom->Nlevel); ij++) {
			for (k = 0; k < atmos.Nspace; k++)
				atom->rhth[nt].Gamma[ij][k] = 0.0;
		}
	}
}
/* ------- end ---------------------------- zeroGammaThreads.c ------ */

/* ------- begin -------------------------- zeroGammaThreadsMolecule.c */

void zeroGammaThreadsMolecule(Molecule *molecule)
{
	register int nt, ij, k;

	if (input.Nthreads <= 1) return;

	for (nt = 0; nt < input.Nthreads; nt++) {
		if (molecule->rhth[nt].Gamma == NULL) continue;
		for (ij = 0; ij < SQ(molecule->Nv); ij++) {
			for (k = 0; k < atmos.Nspace; k++)
				molecule->rhth[nt].Gamma[ij][k] = 0.0;
		}
	}
}
/* ------- end ---------------------------- zeroGammaThreadsMolecule.c */

/* ------- begin -------------------------- mergeGammaThreads.c ----- */

/* --- Merge every thread's private Gamma contributions (accumulated
	   by addtoGamma() during the wavelength sweep, with no locking,
	   since each thread only ever touches its own rhth[nt].Gamma)
	   into the shared atom->Gamma, which initGammaAtom() has already
	   seeded with the collisional rates for this iteration. Call
	   after solveSpectrum()'s eval_operator pass completes, before
	   the result is used to update populations. --      -------------- */

void mergeGammaThreads(Atom *atom)
{
	register int nt, ij, k;

	if (input.Nthreads <= 1) return; /* already accumulated in place */

	for (nt = 0; nt < input.Nthreads; nt++) {
		if (atom->rhth[nt].Gamma == NULL) continue;
		for (ij = 0; ij < SQ(atom->Nlevel); ij++) {
			for (k = 0; k < atmos.Nspace; k++)
				atom->Gamma[ij][k] += atom->rhth[nt].Gamma[ij][k];
		}
	}
}
/* ------- end ---------------------------- mergeGammaThreads.c ----- */

/* ------- begin -------------------------- mergeGammaThreadsMolecule.c */

void mergeGammaThreadsMolecule(Molecule *molecule)
{
	register int nt, ij, k;

	if (input.Nthreads <= 1) return;

	for (nt = 0; nt < input.Nthreads; nt++) {
		if (molecule->rhth[nt].Gamma == NULL) continue;
		for (ij = 0; ij < SQ(molecule->Nv); ij++) {
			for (k = 0; k < atmos.Nspace; k++)
				molecule->Gamma[ij][k] += molecule->rhth[nt].Gamma[ij][k];
		}
	}
}
/* ------- end ---------------------------- mergeGammaThreadsMolecule.c */

/* ------- begin -------------------------- addPreconditionedRates.c */

void addPreconditionedRates(int Nlevel, int i, int j, double twohnu3_c2,
							const double *Vij, const double *wla,
							const double *gij, const double *n_i,
							const double *n_j, const double *U_i,
							const double *U_j, const double *I,
							const double *Psi, double wmu, const double *wk, int Nspace,
							double **Gam)
{
	register int k;

	int    ij = i * Nlevel + j, ji = j * Nlevel + i;
	double Ieff, chi, w, wlamu;

	/* --- Radiative rates of transition i -> j (i lower) with the
		   preconditioning of Rybicki & Hummer (1992, A&A 262, 209),
		   for angle weight wmu, times wk[k] if wk is not NULL,
		   restricted to the two levels of the transition (their sect.
		   2.4, extended to overlapping transitions).

		   U_a is the emission coefficient summed over the active
		   transitions with upper level a (Uji_down), n_i, n_j the
		   populations of the current (old) iteration, chi the old
		   opacity of this transition. Writing

			 I = Ieff + Psi (n_j U_j + n_i U_i),
			 Ieff = I - Psi (n_j U_j + n_i U_i)   (old populations),

		   and keeping the populations in the bracket as unknowns adds
		   to the net rate from i to j the term chi Psi (n_j U_j + n_i U_i):
		   -chi Psi U_j to Gamma[i][j] (j -> i) and +chi Psi U_i to
		   Gamma[j][i] (i -> j). Emission from other levels stays in Ieff.

		   At convergence Ieff + Psi (n_j U_j + n_i U_i) = I exactly, so
		   the populations solve statistical equilibrium. Before (B-6)
		   Ieff subtracted the emission of all active transitions while
		   Gamma represented only part of it, which moved the converged
		   populations, by up to 30 % in Mg II. Coupling every pair of
		   levels instead (full preconditioning, their eq. 2.19) is also
		   exact but puts negative rates between levels that no transition
		   connects, and diverges for multi-ion atoms. --  -------------- */

	for (k = 0; k < Nspace; k++) {
		Ieff  = I[k] - Psi[k] * (n_j[k] * U_j[k] + n_i[k] * U_i[k]);
		chi   = Vij[k] * wla[k] * (n_i[k] - gij[k] * n_j[k]);
		w     = (wk != NULL) ? wmu * wk[k] : wmu;
		wlamu = Vij[k] * wla[k] * w;

		Gam[ji][k] += Ieff * wlamu + chi * Psi[k] * U_i[k] * w;
		Gam[ij][k] += (twohnu3_c2 + Ieff) * gij[k] * wlamu -
					  chi * Psi[k] * U_j[k] * w;
	}
}
/* ------- end ---------------------------- addPreconditionedRates.c */

/* ------- begin -------------------------- addtoGamma.c ------------ */

void addtoGamma(int nspect, double wmu, double *I, double *Psi)
{
	const char   routineName[] = "addtoGamma";
	register int nact, n, k;

	int    i, j, ij, ji, nt;
	double twohnu3_c2, twohc, wlamu, **Gam;

	Atom            *atom;
	AtomicLine      *line;
	AtomicContinuum *continuum;
	Molecule        *molecule;
	MolecularLine   *mrt;
	ActiveSet       *as;

	twohc = 2.0 * HPLANCK * CLIGHT / CUBE(NM_TO_M);

	as = &spectrum.as[nspect];
	nt = nspect % input.Nthreads;

	/* --- Contributions from the active transitions in atoms -- ------ */

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];

		/* --- With a single thread, accumulate straight into the shared
			   Gamma, as a serial run always did. With more threads each
			   thread has its own copy (no lock needed), merged into the
			   shared one after the sweep. --              -------------- */

		Gam = (input.Nthreads > 1) ? atom->rhth[nt].Gamma : atom->Gamma;

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

			/* --- Gam is either the shared Gamma (one thread) or this
				   thread's own copy rhth[nt].Gamma, which only thread nt
				   touches and which is merged into atom->Gamma once the
				   sweep is done (see iterate.c). Cross-coupling terms,
				   currently only for Stokes_I --             -------------- */

			addPreconditionedRates(atom->Nlevel, i, j, twohnu3_c2,
								   atom->rhth[nt].Vij[n], atom->rhth[nt].wla[n],
								   atom->rhth[nt].gij[n], atom->n[i], atom->n[j],
								   atom->rhth[nt].Uji_down[i],
								   atom->rhth[nt].Uji_down[j], I, Psi, wmu, NULL,
								   atmos.Nspace, Gam);
		}
	}
	/* --- Add the active molecular contributions --     -------------- */

	for (nact = 0; nact < atmos.Nactivemol; nact++) {
		molecule = atmos.activemols[nact];
		Gam      = (input.Nthreads > 1) ? molecule->rhth[nt].Gamma : molecule->Gamma;

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

			/* --- In case of molecular vibration-rotation transitions -- - */

			ij = i * molecule->Nv + j;
			ji = j * molecule->Nv + i;

			for (k = 0; k < atmos.Nspace; k++) {
				if (molecule->n[k]) {
					wlamu = molecule->rhth[nt].Vij[n][k] *
							molecule->rhth[nt].wla[n][k] * wmu;
					Gam[ji][k] += I[k] * wlamu;
					Gam[ij][k] += molecule->rhth[nt].gij[n][k] *
								  (twohnu3_c2 + I[k]) * wlamu;
				}
			}
		}
	}
}
/* ------- end ---------------------------- addtoGamma.c ------------ */

/* ------- begin -------------------------- addtoCoupling.c --------- */

void addtoCoupling(int nspect)
{
	const char   routineName[] = "addtoCoupling";
	register int nact, n, k;

	int              i, j, nt;
	double           twohnu3_c2, twohc;
	Atom            *atom;
	AtomicLine      *line;
	AtomicContinuum *continuum;
	ActiveSet       *as;

	twohc = 2.0 * HPLANCK * CLIGHT / CUBE(NM_TO_M);

	as = &spectrum.as[nspect];
	nt = nspect % input.Nthreads;

	/* --- Emission coefficient Uji_down[a] summed over the active
		   transitions with upper level a, for the preconditioning in
		   addPreconditionedRates. It is zero for levels that are only
		   lower levels at this wavelength --             -------------- */

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];

		for (n = 0; n < as->Nlower[nact]; n++) {
			i = as->lower_levels[nact][n];
			for (k = 0; k < atmos.Nspace; k++)
				atom->rhth[nt].Uji_down[i][k] = 0.0;
		}
		for (n = 0; n < as->Nupper[nact]; n++) {
			j = as->upper_levels[nact][n];
			for (k = 0; k < atmos.Nspace; k++)
				atom->rhth[nt].Uji_down[j][k] = 0.0;
		}
	}
	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];

		for (n = 0; n < as->Nactiveatomrt[nact]; n++) {
			switch (as->art[nact][n].type) {
			case ATOMIC_LINE:
				line       = as->art[nact][n].ptype.line;
				j          = line->j;
				twohnu3_c2 = line->Aji / line->Bji;
				break;

			case ATOMIC_CONTINUUM:
				continuum  = as->art[nact][n].ptype.continuum;
				j          = continuum->j;
				twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);
				break;

			default:
				sprintf(messageStr, "Invalid transition type");
				Error(ERROR_LEVEL_1, routineName, messageStr);
				j          = 0;
				twohnu3_c2 = 0.0;
			}
			if (twohnu3_c2) {
				for (k = 0; k < atmos.Nspace; k++) {
					atom->rhth[nt].Uji_down[j][k] +=
						twohnu3_c2 * atom->rhth[nt].gij[n][k] * atom->rhth[nt].Vij[n][k];
				}
			}
		}
	}
}
/* ------- end ---------------------------- addtoCoupling.c --------- */

/* ------- begin -------------------------- zeroRates.c ------------- */

void zeroRates(bool_t redistribute)
{
	register int kr, k, n, nt;

	Atom *atom;

	/* --- Zero the shared line/continuum Rij/Rji, as always, and also
		   every thread's private accumulator (written by addtoRates()
		   with no locking, one slot per thread) ahead of a wavelength
		   sweep. The shared arrays have to be zeroed here because
		   geometries whose formal solution writes them directly
		   (rhsphere) never touch the private ones, and mergeRates() only
		   adds to them.

		   When redistribute == TRUE only the rates of PRD lines are
		   zeroed, matching the set of lines addtoRates() actually
		   updates during a PRD sub-iteration.          -------------- */

	for (n = 0; n < atmos.Natom; n++) {
		atom = &atmos.atoms[n];
		if (atom->active) {
			for (kr = 0; kr < atom->Nline; kr++) {
				if (!redistribute || (redistribute && atom->line[kr].PRD)) {
					for (k = 0; k < atmos.Nspace; k++) {
						atom->line[kr].Rij[k] = 0.0;
						atom->line[kr].Rji[k] = 0.0;
					}
					for (nt = 0; input.Nthreads > 1 && nt < input.Nthreads; nt++) {
						for (k = 0; k < atmos.Nspace; k++) {
							atom->line[kr].Rij_thread[nt][k] = 0.0;
							atom->line[kr].Rji_thread[nt][k] = 0.0;
						}
					}
				}
			}
			if (!redistribute) {
				for (kr = 0; kr < atom->Ncont; kr++) {
					for (k = 0; k < atmos.Nspace; k++) {
						atom->continuum[kr].Rij[k] = 0.0;
						atom->continuum[kr].Rji[k] = 0.0;
					}
					for (nt = 0; input.Nthreads > 1 && nt < input.Nthreads; nt++) {
						for (k = 0; k < atmos.Nspace; k++) {
							atom->continuum[kr].Rij_thread[nt][k] = 0.0;
							atom->continuum[kr].Rji_thread[nt][k] = 0.0;
						}
					}
				}
			}
		}
	}
}
/* ------- end ---------------------------- zeroRates.c ------------- */

/* ------- begin -------------------------- mergeRates.c ------------ */

/* --- Add every thread's private Rij/Rji contributions to the shared
	   line->Rij/Rji, which zeroRates() zeroed at the start of the
	   sweep (and which rhsphere fills directly). Call after a
	   wavelength sweep completes, with the same redistribute flag
	   passed to zeroRates() for that sweep, before anything
	   (PRDScatter, RADRATE_OUTPUT, ...) reads the result.
	   --                                               -------------- */

void mergeRates(bool_t redistribute)
{
	register int kr, k, n, nt;

	Atom *atom;

	if (input.Nthreads <= 1) return; /* rates were accumulated in place */

	for (n = 0; n < atmos.Natom; n++) {
		atom = &atmos.atoms[n];
		if (atom->active) {
			for (kr = 0; kr < atom->Nline; kr++) {
				if (!redistribute || (redistribute && atom->line[kr].PRD)) {
					for (nt = 0; nt < input.Nthreads; nt++) {
						for (k = 0; k < atmos.Nspace; k++) {
							atom->line[kr].Rij[k] += atom->line[kr].Rij_thread[nt][k];
							atom->line[kr].Rji[k] += atom->line[kr].Rji_thread[nt][k];
						}
					}
				}
			}
			if (!redistribute) {
				for (kr = 0; kr < atom->Ncont; kr++) {
					for (nt = 0; nt < input.Nthreads; nt++) {
						for (k = 0; k < atmos.Nspace; k++) {
							atom->continuum[kr].Rij[k] +=
								atom->continuum[kr].Rij_thread[nt][k];
							atom->continuum[kr].Rji[k] +=
								atom->continuum[kr].Rji_thread[nt][k];
						}
					}
				}
			}
		}
	}
}
/* ------- end ---------------------------- mergeRates.c ------------ */

/* ------- begin -------------------------- addtoRates.c ------------ */

void addtoRates(int nspect, int mu, bool_t to_obs, double wmu,
				double *I, bool_t redistribute)
{
	register int nact, n, k;

	int    nt;
	double twohnu3_c2, twohc, wlamu, *Rij, *Rji;

	ActiveSet       *as;
	Atom            *atom;
	AtomicLine      *line;
	AtomicContinuum *continuum;

	(void) mu;
	(void) to_obs;

	/* --- Calculate the radiative rates for atomic transitions.

		   When redistribute == TRUE only radiative rates of PRD lines
		   are evaluated. These are needed in the iterative update of the
		   emission profile ratio \rho.
		   --                                            -------------- */

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
					Rij        = (input.Nthreads > 1) ? line->Rij_thread[nt] : line->Rij;
					Rji        = (input.Nthreads > 1) ? line->Rji_thread[nt] : line->Rji;
					twohnu3_c2 = line->Aji / line->Bji;
				}
				break;

			case ATOMIC_CONTINUUM:
				if (redistribute)
					Rij = NULL;
				else {
					continuum  = as->art[nact][n].ptype.continuum;
					Rij        = (input.Nthreads > 1) ? continuum->Rij_thread[nt] : continuum->Rij;
					Rji        = (input.Nthreads > 1) ? continuum->Rji_thread[nt] : continuum->Rji;
					twohnu3_c2 = twohc / CUBE(spectrum.lambda[nspect]);
				}
				break;

			default:
				Rij = NULL;
			}
			/* --- Convention: Rij is the rate for transition i -> j.
				   With one thread this is the shared Rij/Rji. With more, it
				   is this thread's own copy, which only thread nt touches,
				   merged into the shared arrays once per sweep by
				   mergeRates() after all threads have been joined -- ----- */

			if (Rij != NULL) {
				for (k = 0; k < atmos.Nspace; k++) {
					wlamu =
						atom->rhth[nt].Vij[n][k] * atom->rhth[nt].wla[n][k] * wmu;
					Rij[k] += I[k] * wlamu;
					Rji[k] += atom->rhth[nt].gij[n][k] * (twohnu3_c2 + I[k]) * wlamu;
				}
			}
		}
	}
}
/* ------- end ---------------------------- addtoRates.c ------------ */
