/* ------- file: -------------------------- gchecks.c ---------------

	   Acceptance checks of the map G (docs/SOLVER.md, S-2), run by
	   the gcheck programs of each geometry.

	   --------------------------                      ----------RH-- */

/* --- After the usual setup a gcheck program iterates x <- G(x)
	   Nsteps times from the starting state x0 (plain G, no Ng), calls
	   the end point x_a, and the state lag steps earlier x_b. Then:

	   twice     G(x_a) twice: the two results must be bit-identical.
	   history   G(x_b), then G(x_a) again: must equal the first G(x_a)
				 bit for bit (state that G keeps from one call to the
				 next shows up here, also if it is consistent along the
				 plain iteration path).
	   ecnoise   (with --ecnoise H) the noise delta of every component
				 of G at x_a, with ECnoise (More & Wild 2011, SIAM J.
				 Sci. Comput. 33, 1292): G at x_a + i H p, i = -M/2..M/2
				 (M = --necn, default 8), p random with |p_k| <= |x_a,k|
				 (a relative step). Reported per block of x, relative to
				 |G_k|: max and RMS over the components where noise was
				 detected, and how many were not (H too small: all
				 values equal; H too large: no consistent level).
	   fd        (with --fd) finite-difference Jacobian products
				 D(eps) = (G(x_a + eps v) - G(x_a)) / eps for eps =
				 F ... F 1e-11 (F = --fdmax, default 0.1), v random
				 with |v_k| <= |x_a,k|. Independent relative steps of
				 the populations can invert closely spaced levels (H_20
				 Rydberg levels), which gives negative line opacity
				 and stops the formal solver: lower --fdmax then.
	   --perturb NAME restricts the ECnoise and FD directions to the
	   blocks of x whose name starts with NAME (e.g. "ne").
				 Reported per block: |D(eps) - D(10 eps)| / |D(eps)|
				 (2-norms). It falls as O(eps) while truncation
				 dominates and rises as O(delta/eps) where noise does.

	   Output lines start with "gcheck:". Nothing is written to the
	   output files of RH.
	   --                                              -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"
#include "solver.h"

/* --- Global variables --                             -------------- */

extern Atmosphere  atmos;
extern InputData   input;
extern CommandLine commandline;
extern char        messageStr[];

#define GCHECK_SEED_P 12345u
#define GCHECK_SEED_V 54321u
#define NEPS          12

/* ------- begin -------------------------- gcheckParse.c ----------- */

void gcheckParse(int *argc, char *argv[], GCheckOptions *opt)
{
	int n, m;

	/* --- Take the gcheck options out of argv, leave the rest for
		   RH's own setOptions() --                      -------------- */

	opt->Nsteps   = 30;
	opt->lag      = 3;
	opt->Necn     = 8;
	opt->hEcnoise = 0.0;
	opt->fd       = FALSE;
	opt->fdmax    = 1.0E-1;
	opt->perturb  = NULL;

	for (n = 1, m = 1; n < *argc; n++) {
		if (!strcmp(argv[n], "--steps") && n + 1 < *argc)
			opt->Nsteps = atoi(argv[++n]);
		else if (!strcmp(argv[n], "--lag") && n + 1 < *argc)
			opt->lag = atoi(argv[++n]);
		else if (!strcmp(argv[n], "--ecnoise") && n + 1 < *argc)
			opt->hEcnoise = atof(argv[++n]);
		else if (!strcmp(argv[n], "--necn") && n + 1 < *argc)
			opt->Necn = atoi(argv[++n]);
		else if (!strcmp(argv[n], "--fd"))
			opt->fd = TRUE;
		else if (!strcmp(argv[n], "--fdmax") && n + 1 < *argc)
			opt->fdmax = atof(argv[++n]);
		else if (!strcmp(argv[n], "--perturb") && n + 1 < *argc)
			opt->perturb = argv[++n];
		else
			argv[m++] = argv[n];
	}
	*argc = m;
	/* --- RH's log is not open yet: report on stderr -- ----------- */

	if (opt->lag < 1 || opt->lag > opt->Nsteps) {
		fprintf(stderr, "gcheck: need 1 <= --lag <= --steps\n");
		exit(RH_EXIT_FATAL);
	}
	if (opt->Necn < 4 || opt->Necn % 2 || opt->Necn > 64) {
		fprintf(stderr, "gcheck: --necn must be even, from 4 to 64\n");
		exit(RH_EXIT_FATAL);
	}
}
/* ------- end ---------------------------- gcheckParse.c ----------- */

/* ------- begin -------------------------- randomDirection.c ------- */

static void randomDirection(const double *x, double *p, long N,
							unsigned int seed)
{
	register long k;

	/* --- p_k = r_k |x_k| with r_k uniform in [-1, 1], from a fixed
		   seed (a small LCG, so that it does not depend on the C
		   library) --                                   -------------- */

	for (k = 0; k < N; k++) {
		seed = 1664525u * seed + 1013904223u;
		p[k] = (2.0 * (seed / 4294967296.0) - 1.0) * fabs(x[k]);
	}
}
/* ------- end ---------------------------- randomDirection.c ------- */

/* ------- begin -------------------------- restrictDirection.c ----- */

static void restrictDirection(double *p, const char *prefix)
{
	register long k;
	register int  n;

	long off;

	/* --- With --perturb NAME: zero the direction outside the blocks
		   whose name starts with NAME (e.g. "ne", "populations",
		   "populations CA", "rho", "J") --              -------------- */

	if (prefix == NULL) return;
	for (n = 0, off = 0; n < stateNblocks(); off += stateBlockSize(n), n++) {
		if (strncmp(stateBlockName(n), prefix, strlen(prefix)))
			for (k = off; k < off + stateBlockSize(n); k++) p[k] = 0.0;
	}
}
/* ------- end ---------------------------- restrictDirection.c ----- */

/* ------- begin -------------------------- compare.c --------------- */

static void compare(const char *what, const double *a, const double *b,
					long N)
{
	register long k;

	long   Ndiff = 0;
	double dmax  = 0.0;

	for (k = 0; k < N; k++) {
		if (memcmp(a + k, b + k, sizeof(double))) {
			Ndiff++;
			if (a[k] != 0.0) dmax = MAX(dmax, fabs((b[k] - a[k]) / a[k]));
		}
	}
	if (Ndiff == 0)
		sprintf(messageStr, "gcheck: %-8s bit-identical (%ld values)\n",
				what, N);
	else
		sprintf(messageStr, "gcheck: %-8s DIFFERS: %ld of %ld values,"
							" max relative difference %.3E\n",
				what, Ndiff, N, dmax);
	Error(MESSAGE, NULL, messageStr);
}
/* ------- end ---------------------------- compare.c --------------- */

/* ------- begin -------------------------- ecnoise.c --------------- */

static int ecnoise(int nf, double *f, double *fnoise)
{
	register int i, j;

	int    k, dsgn[64];
	double level[64], gamma = 1.0, fmin, fmax, emin, emax, sum;

	/* --- ECnoise of More & Wild (2011) for nf equally spaced values
		   f (overwritten). Returns 1 (noise detected, *fnoise), 2 (no
		   noise: step too small) or 3 (step too large) -- --------- */

	*fnoise = 0.0;
	fmin = fmax = f[0];
	for (i = 1; i < nf; i++) {
		fmin = MIN(fmin, f[i]);
		fmax = MAX(fmax, f[i]);
	}
	if (fmax == fmin) return 2;
	if ((fmax - fmin) / MAX(fabs(fmax), fabs(fmin)) > 0.1) return 3;

	for (j = 1; j < nf; j++) {
		for (i = 0; i < nf - j; i++) f[i] = f[i + 1] - f[i];

		gamma *= 0.5 * ((double) j / (2.0 * j - 1.0));
		emin = emax = f[0];
		sum         = 0.0;
		for (i = 0; i < nf - j; i++) {
			sum += f[i] * f[i];
			emin = MIN(emin, f[i]);
			emax = MAX(emax, f[i]);
		}
		level[j - 1] = sqrt(gamma * sum / (nf - j));
		dsgn[j - 1]  = (emin * emax < 0.0);
	}
	for (k = 0; k < nf - 3; k++) {
		emin = MIN(MIN(level[k], level[k + 1]), level[k + 2]);
		emax = MAX(MAX(level[k], level[k + 1]), level[k + 2]);
		if (emax <= 4.0 * emin && dsgn[k]) {
			*fnoise = level[k];
			return 1;
		}
	}
	return 3;
}
/* ------- end ---------------------------- ecnoise.c --------------- */

/* ------- begin -------------------------- gcheckRun.c ------------- */

void gcheckRun(const GCheckOptions *opt)
{
	register long k;
	register int  n, i;

	int     Nb, nf, inform, Nnoise, Nsmall, Nlarge, ne;
	long    N, off;
	double *x, *xb, *y1, *y2, *p, **f, *col, fn, rel, nmax, nsum,
		t_start, *D[2], eps, num, den, dmax;

	N  = stateInit();
	Nb = stateNblocks();
	x  = (double *) malloc(N * sizeof(double));
	xb = (double *) malloc(N * sizeof(double));
	y1 = (double *) malloc(N * sizeof(double));
	y2 = (double *) malloc(N * sizeof(double));

	sprintf(messageStr, "\ngcheck: state vector of %ld values in %d blocks:",
			N, Nb);
	Error(MESSAGE, NULL, messageStr);
	for (n = 0; n < Nb; n++) {
		sprintf(messageStr, " %s (%ld)%s", stateBlockName(n),
				stateBlockSize(n), (n < Nb - 1) ? "," : "\n");
		Error(MESSAGE, NULL, messageStr);
	}
	sprintf(messageStr, "gcheck: %d steps of plain G from x0, x_b = %d steps"
						" before the end; PRD sub-iterations per G: %d;"
						" threads: %d\n",
			opt->Nsteps, opt->lag,
			(atmos.NPRDactive > 0) ? input.PRD_NmaxIter : 0, input.Nthreads);
	Error(MESSAGE, NULL, messageStr);

	/* --- Iterate plain G to x_a, keep x_b --          -------------- */

	stateGather(x);
	t_start = wallClock();
	for (i = 1; i <= opt->Nsteps; i++) {
		if (i == opt->Nsteps - opt->lag + 1) memcpy(xb, x, N * sizeof(double));
		G_eval(x, x);
	}

	/* --- Twice and history --                         -------------- */

	G_eval(x, y1);
	G_eval(x, y2);
	compare("twice", y1, y2, N);
	G_eval(xb, y2);
	G_eval(x, y2);
	compare("history", y1, y2, N);

	for (n = 0, off = 0; n < Nb; off += stateBlockSize(n), n++) {
		dmax = 0.0;
		for (k = off; k < off + stateBlockSize(n); k++)
			if (y1[k]) dmax = MAX(dmax, fabs((y1[k] - x[k]) / y1[k]));
		sprintf(messageStr, "gcheck: at x_a, max |G(x_a) - x_a| / |G(x_a)|"
							" of %s = %.3E\n",
				stateBlockName(n), dmax);
		Error(MESSAGE, NULL, messageStr);
	}

	/* --- ECnoise --                                   -------------- */

	if (opt->hEcnoise > 0.0) {
		nf = opt->Necn + 1;
		p  = (double *) malloc(N * sizeof(double));
		f  = matrix_double(nf, N);
		randomDirection(x, p, N, GCHECK_SEED_P);
		restrictDirection(p, opt->perturb);
		for (i = 0; i < nf; i++) {
			for (k = 0; k < N; k++)
				y2[k] = x[k] + (i - opt->Necn / 2) * opt->hEcnoise * p[k];
			G_eval(y2, f[i]);
		}
		col = (double *) malloc(nf * sizeof(double));
		sprintf(messageStr, "gcheck: ecnoise, %d points, relative step %.1E:"
							" per block max and RMS of delta / |G|, components"
							" with noise / step too small / too large\n",
				nf, opt->hEcnoise);
		Error(MESSAGE, NULL, messageStr);
		for (n = 0, off = 0; n < Nb; off += stateBlockSize(n), n++) {
			Nnoise = Nsmall = Nlarge = 0;
			nmax = nsum = 0.0;
			for (k = off; k < off + stateBlockSize(n); k++) {
				for (i = 0; i < nf; i++) col[i] = f[i][k];
				inform = ecnoise(nf, col, &fn);
				if (inform == 1 && f[opt->Necn / 2][k] != 0.0) {
					rel  = fn / fabs(f[opt->Necn / 2][k]);
					nmax = MAX(nmax, rel);
					nsum += rel * rel;
					Nnoise++;
				} else if (inform == 2)
					Nsmall++;
				else if (inform == 3)
					Nlarge++;
			}
			sprintf(messageStr, "gcheck: ecnoise %-24s max %.3E  RMS %.3E"
								"  (%ld / %d / %d)\n",
					stateBlockName(n), nmax,
					(Nnoise) ? sqrt(nsum / Nnoise) : 0.0, (long) Nnoise,
					Nsmall, Nlarge);
			Error(MESSAGE, NULL, messageStr);
		}
		free(col);
		free(p);
		freeMatrix((void **) f);
	}

	/* --- Finite-difference Jacobian products --       -------------- */

	if (opt->fd) {
		p    = (double *) malloc(N * sizeof(double));
		D[0] = (double *) malloc(N * sizeof(double));
		D[1] = (double *) malloc(N * sizeof(double));
		randomDirection(x, p, N, GCHECK_SEED_V);
		restrictDirection(p, opt->perturb);
		Error(MESSAGE, NULL, "gcheck: fd, |D(eps) - D(10 eps)| / |D(eps)| per block\n");
		sprintf(messageStr, "gcheck: fd %8s", "eps");
		Error(MESSAGE, NULL, messageStr);
		for (n = 0; n < Nb; n++) {
			sprintf(messageStr, " %16.16s", stateBlockName(n));
			Error(MESSAGE, NULL, messageStr);
		}
		Error(MESSAGE, NULL, "\n");
		for (ne = 0, eps = opt->fdmax; ne < NEPS; ne++, eps *= 0.1) {
			for (k = 0; k < N; k++) y2[k] = x[k] + eps * p[k];
			G_eval(y2, D[ne % 2]);
			for (k = 0; k < N; k++) D[ne % 2][k] = (D[ne % 2][k] - y1[k]) / eps;
			if (ne == 0) continue;
			sprintf(messageStr, "gcheck: fd %8.1E", eps);
			Error(MESSAGE, NULL, messageStr);
			for (n = 0, off = 0; n < Nb; off += stateBlockSize(n), n++) {
				num = den = 0.0;
				for (k = off; k < off + stateBlockSize(n); k++) {
					num += SQ(D[ne % 2][k] - D[(ne + 1) % 2][k]);
					den += SQ(D[ne % 2][k]);
				}
				sprintf(messageStr, " %16.3E", (den > 0.0) ? sqrt(num / den) : 0.0);
				Error(MESSAGE, NULL, messageStr);
			}
			Error(MESSAGE, NULL, "\n");
		}
		free(p);
		free(D[0]);
		free(D[1]);
	}
	if (opt->perturb) {
		sprintf(messageStr, "gcheck: ecnoise and fd directions restricted to"
							" blocks \"%s...\"\n",
				opt->perturb);
		Error(MESSAGE, NULL, messageStr);
	}
	sprintf(messageStr, "gcheck: %ld evaluations of G, %.3f s per evaluation\n",
			G_count(), (wallClock() - t_start) / G_count());
	Error(MESSAGE, NULL, messageStr);

	free(x);
	free(xb);
	free(y1);
	free(y2);
}
/* ------- end ---------------------------- gcheckRun.c ------------- */
