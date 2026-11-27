/* ------- file: -------------------------- statequil.c -------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Apr  1 14:09:35 2009 --

	   --------------------------                      ----------RH-- */

/* --- Solves equations of statistical equilibrium with given rate
	   matrix Gamma[Nlevel*Nlevel][Nspace], where Gamma[i*Nlevel + j]
	   is the rate from level j to level i.

	   At each spatial location the equations are solved by GTH
	   elimination (gth.c, N-31), which is accurate to rounding in every
	   population whatever the conditioning of the rate matrix. GTH needs
	   non-negative rates; where a rate is negative (bad atomic data, or
	   the preconditioning of the radiative rates) or a level has no
	   rates to the others, the equations are solved by LU decomposition
	   (SolveLinearEq) as before, and the run log says where and why,
	   once per call. The number of such locations is also counted for
	   the work summary (work.c).

	   The LU populations are not guaranteed non-negative, and a
	   negative or non-finite one (N-8) is not allowed to go on: during
	   a trial evaluation of G the trial is invalidated (negTauInvalid,
	   as for a negative total opacity), otherwise the run stops. A
	   population is negative below -SE_NEG_TOL times the total, which
	   leaves the rounding of the LU solution of a nearly empty level.

	   isum is the reference level: GTH eliminates it last, and LU
	   replaces its equation by particle conservation. If set to -1 then
	   at each spatial location the level with the largest population
	   from the previous iteration is used. -- -------- */

#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "accelerate.h"
#include "statistics.h"
#include "error.h"
#include "negtau.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* --- The first location of a call where GTH could not be used, for
	   the message --                                  -------------- */

#define SE_NEG_TOL 1e-10

/* --- The first negative or non-finite LU population of a call
	   (N-8) is kept too --                            -------------- */

typedef struct {
	int            Nlu, k, i, j;
	enum GTHStatus status;
	double         rate;
	int            Nbad, kbad, ibad;
	double         nbad;
} SEFallback;

#define SE_FALLBACK_INIT {0, -1, -1, -1, GTH_OK, 0.0, 0, -1, -1, 0.0}

/* ------- begin -------------------------- solveRates.c ------------ */

static void solveRates(RHState *rh, int N, double **Gamma, int k, int isum,
					   double ntotal, double *n, double **R,
					   double **A, SEFallback *fb)
{
	register int i, j;

	int            root, i_bad, j_bad;
	double         nmax, GamDiag;
	enum GTHStatus status;

	/* --- Reference level: the most populated one, or isum -- ------ */

	if (isum == -1) {
		root = 0;
		nmax = 0.0;
		for (i = 0; i < N; i++) {
			if (n[i] > nmax) {
				nmax = n[i];
				root = i;
			}
		}
	} else
		root = isum;

	/* --- R[i][j] is the rate from level i to level j -- ------------ */

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) R[i][j] = (i == j) ? 0.0 : Gamma[j * N + i][k];

	status = SolveGTH(N, R, root, ntotal, n, &i_bad, &j_bad);
	if (status == GTH_OK) return;

	if (fb->Nlu++ == 0) {
		fb->k      = k;
		fb->i      = i_bad;
		fb->j      = j_bad;
		fb->status = status;
		fb->rate   = (j_bad >= 0) ? Gamma[j_bad * N + i_bad][k] : 0.0;
	}
	/* --- LU: for each column i sum over rows to get the diagonal, then
		   close the homogeneous set with particle conservation -- ---- */

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) A[i][j] = Gamma[i * N + j][k];
	for (i = 0; i < N; i++) {
		GamDiag = 0.0;
		A[i][i] = 0.0;
		n[i]    = 0.0;

		for (j = 0; j < N; j++) GamDiag += A[j][i];
		A[i][i] = -GamDiag;
	}
	n[root] = ntotal;
	for (j = 0; j < N; j++) A[root][j] = 1.0;

	SolveLinearEq(rh, N, A, n, TRUE);

	for (i = 0; i < N; i++) {
		if (isfinite(n[i]) && n[i] >= -SE_NEG_TOL * ntotal) continue;
		if (fb->Nbad++ == 0) {
			fb->kbad = k;
			fb->ibad = i;
			fb->nbad = n[i];
		}
		break;
	}
}
/* ------- end ---------------------------- solveRates.c ------------ */

/* ------- begin -------------------------- reportFallback.c -------- */

static void reportFallback(RHState *rh, const char *ID, const SEFallback *fb,
						   double **C, int N)
{
	char why[160];

	switch (fb->status) {
	case GTH_NEGATIVE:
		sprintf(why, "negative rate %d -> %d (%.2E s^-1 at k = %d%s)", fb->i,
				fb->j, fb->rate, fb->k,
				(C == NULL)                           ? ""
				: (C[fb->j * N + fb->i][fb->k] < 0.0) ? ", collisional"
													  : ", radiative");
		break;
	case GTH_NONFINITE:
		sprintf(why, "non-finite rate or result at k = %d", fb->k);
		break;
	case GTH_REDUCIBLE:
		sprintf(why, "level %d without rates to the others at k = %d",
				fb->i, fb->k);
		break;
	default:
		why[0] = '\0';
	}
	Message(rh, " %s: rate equations solved with LU at %d of %d"
				" locations (first: %s)\n",
			ID, fb->Nlu, rh->atmos.Nspace, why);
}
/* ------- end ---------------------------- reportFallback.c -------- */

/* ------- begin -------------------------- badPopulations.c ------- */

static void badPopulations(RHState *rh, const char *ID, const char *level, long Nsolved,
						   const SEFallback *fb)
{
	const char routineName[] = "badPopulations";
	char       what[440];

	/* --- Only reached when LU produced a population that is negative
		   or not finite (N-8): invalidates a trial point, else stops
		   the run --                                   -------------- */

	/* --- A negative population comes from a negative rate (the MALI
		   preconditioning in a thick line core, or bad atomic data), a
		   non-finite one from a non-finite rate or an overflow. The
		   initial solution (initial_xdr.c) can also be the source; it
		   is never a trial, so it stops with this message -- -------- */

	snprintf(what, sizeof(what),
			 "%s: the LU solution of the rate equations has a %s population"
			 " (%.4E) of level %d%s at location %d (%ld of %ld locations"
			 " solved have one); %s",
			 ID, isfinite(fb->nbad) ? "negative" : "non-finite", fb->nbad,
			 fb->ibad, level, fb->kbad, (long) fb->Nbad, Nsolved,
			 isfinite(fb->nbad)
				 ? "negative populations come from negative rates (the MALI"
				   " preconditioning in thick line cores, or bad atomic data)"
				 : "a non-finite population comes from a non-finite rate (bad"
				   " atomic data, overflow)");

	if (negTauInTrial(rh)) {
		negTauInvalid(rh, what);
		return;
	}
	Error(rh, ERROR_LEVEL_2, routineName, "\n  %s", what);
}
/* ------- end ---------------------------- badPopulations.c ------- */

/* ------- begin -------------------------- statEquil.c ------------- */

void statEquil(RHState *rh, Atom *atom, int isum)
{
	register int i, k;

	int        Nlevel;
	double    *n_k, **R, **A;
	SEFallback fb = SE_FALLBACK_INIT;

	getCPU(rh, 3, TIME_START, NULL);

	Nlevel = atom->Nlevel;

	/* --- Need temporary storage because Gamma has to be solved spatial
		   point by spatial point while depth is normally the fastest
		   running index --                              -------------- */

	n_k = (double *) malloc(Nlevel * sizeof(double));
	R   = matrix_double(rh, Nlevel, Nlevel);
	A   = matrix_double(rh, Nlevel, Nlevel);

	for (k = 0; k < rh->atmos.Nspace; k++) {
		for (i = 0; i < Nlevel; i++) n_k[i] = atom->n[i][k];

		solveRates(rh, Nlevel, atom->Gamma, k, isum, atom->ntotal[k], n_k, R, A,
				   &fb);
		for (i = 0; i < Nlevel; i++) atom->n[i][k] = n_k[i];
	}
	workRates(rh, rh->atmos.Nspace, fb.Nlu);
	if (fb.Nlu) reportFallback(rh, atom->ID, &fb, atom->C, Nlevel);
	if (fb.Nbad) {
		char level[ATOM_LABEL_WIDTH + 8];

		sprintf(level, " (%s)", atom->label[fb.ibad]);
		badPopulations(rh, atom->ID, level, rh->atmos.Nspace, &fb);
	}

	free(n_k);
	freeMatrix(rh, (void **) R);
	freeMatrix(rh, (void **) A);

	getCPU(rh, 3, TIME_POLL, "Stat Equil");
}
/* ------- end ---------------------------- statEquil.c ------------- */

/* ------- begin -------------------------- statEquilMolecule.c ----- */

void statEquilMolecule(RHState *rh, struct Molecule *molecule, int isum)
{
	register int vi, k;

	int        Nlevel, Nsolved = 0;
	double    *n_k, **R, **A;
	SEFallback fb = SE_FALLBACK_INIT;

	getCPU(rh, 3, TIME_START, NULL);

	Nlevel = molecule->Nv;

	n_k = (double *) malloc(Nlevel * sizeof(double));
	R   = matrix_double(rh, Nlevel, Nlevel);
	A   = matrix_double(rh, Nlevel, Nlevel);

	for (k = 0; k < rh->atmos.Nspace; k++) {
		if (molecule->n[k] > 0.0) {
			for (vi = 0; vi < Nlevel; vi++) n_k[vi] = molecule->nv[vi][k];

			solveRates(rh, Nlevel, molecule->Gamma, k, isum, molecule->n[k], n_k,
					   R, A, &fb);
			for (vi = 0; vi < Nlevel; vi++) molecule->nv[vi][k] = n_k[vi];
			Nsolved++;
		} else
			for (vi = 0; vi < Nlevel; vi++) molecule->nv[vi][k] = 0.0;
	}
	workRates(rh, Nsolved, fb.Nlu);
	if (fb.Nlu) reportFallback(rh, molecule->ID, &fb, NULL, Nlevel);
	if (fb.Nbad) badPopulations(rh, molecule->ID, "", Nsolved, &fb);

	free(n_k);
	freeMatrix(rh, (void **) R);
	freeMatrix(rh, (void **) A);

	getCPU(rh, 3, TIME_POLL, "Stat Equil");
}
/* ------- end ---------------------------- statEquilMolecule.c ----- */

/* ------- begin -------------------------- updatePopulations.c ----- */

double updatePopulations(RHState *rh, int niter)
{
	register int nact;

	char      label[MAX_MESSAGE_LENGTH];
	bool_t    accel, quiet;
	double    dpops, dpopsmax = 0.0;
	Atom     *atom;
	Molecule *molecule;

	(void) niter;

	/* --- Update active atoms --                        -------------- */

	for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
		atom = rh->atmos.activeatoms[nact];

		statEquil(rh, atom, rh->input.isum);
		accel = Accelerate(rh, atom->Ng_n, atom->n[0]);

		snprintf(label, sizeof(label), " %s,", atom->ID);
		dpops = MaxChange(rh, atom->Ng_n, label, quiet = FALSE);
		Message(rh, "%s", NgStatus(atom->Ng_n, accel));

		dpopsmax = MAX(dpops, dpopsmax);
	}
	/* --- Update active molecules --                    -------------- */

	for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
		molecule = rh->atmos.activemols[nact];

		statEquilMolecule(rh, molecule, 0);
		accel = Accelerate(rh, molecule->Ng_nv, molecule->nv[0]);

		snprintf(label, sizeof(label), " %s ,", molecule->ID);
		dpops = MaxChange(rh, molecule->Ng_nv, label, quiet = FALSE);
		Message(rh, "%s", NgStatus(molecule->Ng_nv, accel));

		dpopsmax = MAX(dpops, dpopsmax);
	}

	return dpopsmax;
}
/* ------- end ---------------------------- updatePopulations.c ----- */
