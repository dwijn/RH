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

	   isum is the reference level: GTH eliminates it last, and LU
	   replaces its equation by particle conservation. If set to -1 then
	   at each spatial location the level with the largest population
	   from the previous iteration is used. -- -------- */

#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "accelerate.h"
#include "statistics.h"
#include "error.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern InputData  input;
extern char       messageStr[];

/* --- The first location of a call where GTH could not be used, for
	   the message --                                  -------------- */

typedef struct {
	int            Nlu, k, i, j;
	enum GTHStatus status;
	double         rate;
} SEFallback;

/* ------- begin -------------------------- solveRates.c ------------ */

static void solveRates(int N, double **Gamma, int k, int isum,
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

	SolveLinearEq(N, A, n, TRUE);
}
/* ------- end ---------------------------- solveRates.c ------------ */

/* ------- begin -------------------------- reportFallback.c -------- */

static void reportFallback(const char *ID, const SEFallback *fb,
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
	sprintf(messageStr, " %s: rate equations solved with LU at %d of %d"
						" locations (first: %s)\n",
			ID, fb->Nlu, atmos.Nspace, why);
	Error(MESSAGE, NULL, messageStr);
}
/* ------- end ---------------------------- reportFallback.c -------- */

/* ------- begin -------------------------- statEquil.c ------------- */

void statEquil(Atom *atom, int isum)
{
	register int i, k;

	int        Nlevel;
	double    *n_k, **R, **A;
	SEFallback fb = {0, -1, -1, -1, GTH_OK, 0.0};

	getCPU(3, TIME_START, NULL);

	Nlevel = atom->Nlevel;

	/* --- Need temporary storage because Gamma has to be solved spatial
		   point by spatial point while depth is normally the fastest
		   running index --                              -------------- */

	n_k = (double *) malloc(Nlevel * sizeof(double));
	R   = matrix_double(Nlevel, Nlevel);
	A   = matrix_double(Nlevel, Nlevel);

	for (k = 0; k < atmos.Nspace; k++) {
		for (i = 0; i < Nlevel; i++) n_k[i] = atom->n[i][k];

		solveRates(Nlevel, atom->Gamma, k, isum, atom->ntotal[k], n_k, R, A,
				   &fb);
		for (i = 0; i < Nlevel; i++) atom->n[i][k] = n_k[i];
	}
	workRates(atmos.Nspace, fb.Nlu);
	if (fb.Nlu) reportFallback(atom->ID, &fb, atom->C, Nlevel);

	free(n_k);
	freeMatrix((void **) R);
	freeMatrix((void **) A);

	getCPU(3, TIME_POLL, "Stat Equil");
}
/* ------- end ---------------------------- statEquil.c ------------- */

/* ------- begin -------------------------- statEquilMolecule.c ----- */

void statEquilMolecule(struct Molecule *molecule, int isum)
{
	register int vi, k;

	int        Nlevel, Nsolved = 0;
	double    *n_k, **R, **A;
	SEFallback fb = {0, -1, -1, -1, GTH_OK, 0.0};

	getCPU(3, TIME_START, NULL);

	Nlevel = molecule->Nv;

	n_k = (double *) malloc(Nlevel * sizeof(double));
	R   = matrix_double(Nlevel, Nlevel);
	A   = matrix_double(Nlevel, Nlevel);

	for (k = 0; k < atmos.Nspace; k++) {
		if (molecule->n[k] > 0.0) {
			for (vi = 0; vi < Nlevel; vi++) n_k[vi] = molecule->nv[vi][k];

			solveRates(Nlevel, molecule->Gamma, k, isum, molecule->n[k], n_k,
					   R, A, &fb);
			for (vi = 0; vi < Nlevel; vi++) molecule->nv[vi][k] = n_k[vi];
			Nsolved++;
		} else
			for (vi = 0; vi < Nlevel; vi++) molecule->nv[vi][k] = 0.0;
	}
	workRates(Nsolved, fb.Nlu);
	if (fb.Nlu) reportFallback(molecule->ID, &fb, NULL, Nlevel);

	free(n_k);
	freeMatrix((void **) R);
	freeMatrix((void **) A);

	getCPU(3, TIME_POLL, "Stat Equil");
}
/* ------- end ---------------------------- statEquilMolecule.c ----- */

/* ------- begin -------------------------- updatePopulations.c ----- */

double updatePopulations(int niter)
{
	register int nact;

	bool_t    accel, quiet;
	double    dpops, dpopsmax = 0.0;
	Atom     *atom;
	Molecule *molecule;

	(void) niter;

	/* --- Update active atoms --                        -------------- */

	for (nact = 0; nact < atmos.Nactiveatom; nact++) {
		atom = atmos.activeatoms[nact];

		statEquil(atom, input.isum);
		accel = Accelerate(atom->Ng_n, atom->n[0]);

		sprintf(messageStr, " %s,", atom->ID);
		dpops = MaxChange(atom->Ng_n, messageStr, quiet = FALSE);
		Error(MESSAGE, NULL, (accel) ? " (accelerated)\n" : "\n");

		dpopsmax = MAX(dpops, dpopsmax);
	}
	/* --- Update active molecules --                    -------------- */

	for (nact = 0; nact < atmos.Nactivemol; nact++) {
		molecule = atmos.activemols[nact];

		statEquilMolecule(molecule, 0);
		accel = Accelerate(molecule->Ng_nv, molecule->nv[0]);

		sprintf(messageStr, " %s ,", molecule->ID);
		dpops = MaxChange(molecule->Ng_nv, messageStr, quiet = FALSE);
		Error(MESSAGE, NULL, (accel) ? " (accelerated)\n" : "\n");

		dpopsmax = MAX(dpops, dpopsmax);
	}

	return dpopsmax;
}
/* ------- end ---------------------------- updatePopulations.c ----- */
