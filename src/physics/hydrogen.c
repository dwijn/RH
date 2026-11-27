/* ------- file: -------------------------- hydrogen.c --------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Jan 26 15:47:01 2012 --

	   --------------------------                      ----------RH-- */

/* --- Computes hydrogen related bound-free and free-free opacities,
	   and continuum opacities.

	   Currently includes:

	--  Neutral Hydrogen bound-free and free-free.
	--  H^- bound-free and free-free.
	--  H2^- free-free.
	--  H2^+ free-free.
	--  Rayleigh scattering by molecular H2.

  Note: For Rayleigh scattering by neutral hydrogen the general
		routine rayleigh.c is used (called from Background).


	   Global variables:
		atmos -- Atmos structure for atmospheric data.
		atom  -- Atom structure for current active atomic model.

	   Input:
		lambda -- Wavelength [nm] for which opacity and emissivity
				  are to be calculated.
			 H -- Pointer to Atom structure H containing atomic data
				  for hydrogen in atmosphere.

	   Output:
		chi[Nspace] -- Array for opacities [m^2].
		eta[Nspace] -- Array for emissivities [J s^-1 Hz^-1 sr^-1].

	   --                                              -------------- */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdatomic.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "background.h"
#include "error.h"
#include "inputs.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- distribute_nH.c --------- */

void distribute_nH(RHState *rh)
{
	const char   routineName[] = "distribute_nH";
	register int k, i;

	char    config[4], *ptr;
	int    *quantumNo, Nspace = rh->atmos.Nspace, iq;
	double *g_total, g_i;

	/* --- Redistribute the hydrogen levels for use in the background
		   opacity package.

		   Use the following conventions:

		   1 -- To save memory:
				- let atmos.np point to atmos.H.n[atmos.H.Nlevel-1]
				- let nH2 point to atmos.H2->n (= molecules[0]);

		   2 -- atmos.nHtot represents the total number of hydrogen
				atoms in neutral atoms, protons, H-, and H2 and H2+
				molecules. So in general atmos.nHtot != atmos.H.ntotal.

		   3 -- In case of LTE (set by atmos.H_LTE, see keyword.input)
				let atmos.H.n point to the LTE equivalents atmos.H.nstar.

	   --                                                -------------- */

	if (rh->atmos.H_LTE) {
		rh->atmos.H->NLTEpops = FALSE;
		Message(rh, "\nUsing LTE hydrogen populations for background opacities\n\n");

		/* --- To save memory space let atmos.H->n point to LTE populations
			   atmos.H->nstar --                           -------------- */

		if (!rh->atmos.H->active) rh->atmos.H->n = rh->atmos.H->nstar;

	} else {
		rh->atmos.H->NLTEpops = TRUE;
		if (!rh->atmos.H->active)
			rh->atmos.H->n = matrix_double(rh, rh->atmos.H->Nlevel, rh->atmos.Nspace);

		/* --- Find principal quantum number from label --   ------------ */

		quantumNo = (int *) malloc((rh->atmos.H->Nlevel - 1) * sizeof(int));
		g_total   = (double *) calloc(rh->atmos.NHydr - 1, sizeof(double));

		for (i = 0; i < rh->atmos.H->Nlevel - 1; i++) {
			sscanf(rh->atmos.H->label[i], "H I %s", config);
			ptr = config;
			while (isdigit(*ptr)) ptr++;
			*ptr = ' ';
			sscanf(config, "%d", &quantumNo[i]);

			if (quantumNo[i] < rh->atmos.NHydr)
				g_total[quantumNo[i] - 1] += rh->atmos.H->g[i];
		}
		/* --- Now redistribute atmosphere's Hydrogen populations over
			   the levels of atom H --                     -------------- */

		for (i = 0; i < rh->atmos.H->Nlevel - 1; i++) {
			if (quantumNo[i] < rh->atmos.NHydr) {
				iq  = quantumNo[i] - 1;
				g_i = rh->atmos.H->g[i] / g_total[iq];
				for (k = 0; k < Nspace; k++)
					rh->atmos.H->n[i][k] = g_i * rh->atmos.nH[iq][k];
			} else {
				Warning(rh, routineName, "Too many hydrogen levels (level n = %d)\n"
										 " Background opacity additional levels set to zero%s",
						quantumNo[i], (i == rh->atmos.H->Nlevel - 2) ? "\n\n" : "");
			}
		}
		free(quantumNo);
		free(g_total);

		/* --- The protons come last --                    -------------- */

		for (k = 0; k < Nspace; k++)
			rh->atmos.H->n[rh->atmos.H->Nlevel - 1][k] = rh->atmos.nH[rh->atmos.NHydr - 1][k];
	}
	/* --- Free memory for atmospheric populations --    -------------- */

	freeMatrix(rh, (void **) rh->atmos.nH);
}
/* ------- end ---------------------------- distribute_nH.c --------- */

/* ------- begin -------------------------- Hydrogen_bf.c ----------- */

bool_t Hydrogen_bf(RHState *rh, double lambda, double *chi, double *eta)
{
	/* --- Hydrogen bound-free opacity

	  See: Mihalas (1978) p. 99 --                       -------------- */

	register int k, kr;

	bool_t opaque;
	int    i;
	double lambdaEdge, sigma, sigma0, g_bf, twohnu3_c2, twohc, gijk,
		hc_k, hc_kla, *npstar, expla, n_eff, *np;
	AtomicContinuum *continuum;

	opaque = FALSE;
	for (k = 0; k < rh->atmos.Nspace; k++) {
		chi[k] = 0.0;
		eta[k] = 0.0;
	}
	if (rh->atmos.H->active) return opaque;

	twohc  = (2.0 * HPLANCK * CLIGHT) / CUBE(NM_TO_M);
	hc_k   = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M);
	sigma0 = 32.0 / (3.0 * sqrt(3.0)) * SQ(Q_ELECTRON) / (4.0 * PI * EPSILON_0) /
			 (M_ELECTRON * CLIGHT) * HPLANCK / (2.0 * E_RYDBERG);

	npstar = rh->atmos.H->nstar[rh->atmos.H->Nlevel - 1];

	for (kr = 0; kr < rh->atmos.H->Ncont; kr++) {
		continuum  = rh->atmos.H->continuum + kr;
		lambdaEdge = continuum->lambda0;
		i          = continuum->i;

		if (lambda <= lambdaEdge && lambda >= continuum->lambda[0]) {
			opaque = TRUE;

			/* --- Find the principal quantum number of level i -- -------- */

			n_eff = sqrt(E_RYDBERG /
						 (rh->atmos.H->E[continuum->j] - rh->atmos.H->E[continuum->i]));

			g_bf       = Gaunt_bf(lambda, n_eff, rh->atmos.H->stage[i] + 1);
			sigma      = sigma0 * n_eff * g_bf * CUBE(lambda / lambdaEdge);
			hc_kla     = hc_k / lambda;
			twohnu3_c2 = twohc / CUBE(lambda);

			np = rh->atmos.H->n[rh->atmos.H->Nlevel - 1];
			for (k = 0; k < rh->atmos.Nspace; k++) {
				expla = exp(-hc_kla / rh->atmos.T[k]);
				gijk  = rh->atmos.H->nstar[i][k] / npstar[k] * expla;
				chi[k] += sigma * (1.0 - expla) * rh->atmos.H->n[i][k];
				eta[k] += twohnu3_c2 * gijk * sigma * np[k];
			}
		}
	}
	return opaque;
}
/* ------- end ---------------------------- Hydrogen_bf.c ----------- */

/* ------- begin -------------------------- Hydrogen_ff.c ----------- */

void Hydrogen_ff(RHState *rh, double lambda, double *chi)
{
	/* --- Hydrogen free-free opacity

	  See: Mihalas (1978) p. 101
		  --                                             -------------- */

	register int k;

	int    Nspace = rh->atmos.Nspace;
	double hc_kla, C0, sigma, g_ff, stim, nu3, *np;

	C0    = SQ(Q_ELECTRON) / (4.0 * PI * EPSILON_0) / sqrt(M_ELECTRON);
	sigma = 4.0 / 3.0 * sqrt(2.0 * PI / (3.0 * KBOLTZMANN)) * CUBE(C0) /
			(HPLANCK * CLIGHT);
	nu3    = CUBE((lambda * NM_TO_M) / CLIGHT);
	hc_kla = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M * lambda);

	np = rh->atmos.H->n[rh->atmos.H->Nlevel - 1];
	for (k = 0; k < Nspace; k++) {
		stim   = 1.0 - exp(-hc_kla / rh->atmos.T[k]);
		g_ff   = Gaunt_ff(rh, lambda, 1, rh->atmos.T[k]);
		chi[k] = sigma / sqrt(rh->atmos.T[k]) * nu3 * rh->atmos.ne[k] *
				 np[k] * stim * g_ff;
	}
}
/* ------- end ---------------------------- Hydrogen_ff.c ----------- */

/* ------- begin -------------------------- Gaunt_bf.c -------------- */

double Gaunt_bf(double lambda, double n_eff, int charge)
{
	/* --- M. J. Seaton (1960), Rep. Prog. Phys. 23, 313 -- ----------- */

	double x, x3, nsqx;

	x    = ((HPLANCK * CLIGHT) / (lambda * NM_TO_M)) / (E_RYDBERG * SQ(charge));
	x3   = pow(x, 0.33333333);
	nsqx = 1.0 / (SQ(n_eff) * x);

	return 1.0 + 0.1728 * x3 * (1.0 - 2.0 * nsqx) -
		   0.0496 * SQ(x3) * (1.0 - (1.0 - nsqx) * 0.66666667 * nsqx);
}
/* ------- end ---------------------------- Gaunt_bf.c -------------- */

/* ------- begin -------------------------- Gaunt_ff.c -------------- */

/* --- Thermally averaged free-free Gaunt factor (N-50).

 From: P. A. M. van Hoof, R. J. R. Williams, K. Volk, M. Chatzikos,
	   G. J. Ferland, M. Lykins, R. L. Porter and Y. Wang (2014),
	   MNRAS 444, 420: the non-relativistic table <g_ff>(gamma^2, u),
	   gamma^2 = Z^2 Ry / kT and u = h nu / kT, on a grid of 0.2 dex in
	   both (log gamma^2 = -6 to 10, log u = -16 to 13; relative
	   accuracy 1.0E-5). The file is read at run time (keyword
	   GAUNT_FF_DATA, default Atoms/gauntff.dat, the file of
	   http://data.nublado.org/gauntff/ unaltered, with its licence).

 ln g is interpolated with four-point Lagrange polynomials in
 log gamma^2 and log u (bicubic; one-sided stencils at the edges of
 the table): positive by construction, and within about 5.0E-5 of the
 table between its nodes. Outside the table, which covers 1.6E-5 K to
 1.6E11 K and every wavelength RH uses, the run stops.

 Note: The table is read once in hydrogenPrepare, so that Gaunt_ff
	   can be called from several threads.

 Replaces the expansion in (h nu / Ry)^1/3 of Seaton (1960), which
 was used far outside its range and floored at 1: in the mm it gave
 g = 1 where the table gives 2.2 to 3.5. --        -------------- */

#define GFF_MAGIC    20140210
#define GFF_MAX_LINE 4096

static int    gff_Ng = 0, gff_Nu = 0;
static double gff_lg0, gff_lu0, gff_step, *gff_lng = NULL; /* ln g [u][g] */

/* ------- begin -------------------------- gff_weights.c ----------- */

static void gff_weights(double t, double w[4])
{
	/* --- Four-point Lagrange weights at t for nodes -1, 0, 1, 2. -- */

	w[0] = -t * (t - 1.0) * (t - 2.0) / 6.0;
	w[1] = (t + 1.0) * (t - 1.0) * (t - 2.0) / 2.0;
	w[2] = -(t + 1.0) * t * (t - 2.0) / 2.0;
	w[3] = (t + 1.0) * t * (t - 1.0) / 6.0;
}
/* ------- end ---------------------------- gff_weights.c ----------- */

/* ------- begin -------------------------- gff_stencil.c ----------- */

static int gff_stencil(double x, int N, double w[4])
{
	/* --- First node of the four-point stencil around grid
		   coordinate x (0 <= x <= N - 1), and its weights. -- ---- */

	int i0 = (int) floor(x) - 1;

	if (i0 < 0) i0 = 0;
	if (i0 > N - 4) i0 = N - 4;
	gff_weights(x - (i0 + 1), w);
	return i0;
}
/* ------- end ---------------------------- gff_stencil.c ----------- */

/* ------- begin -------------------------- Gaunt_ff_readTable.c ---- */

int Gaunt_ff_readTable(const char *fileName, char *errmsg, size_t len)
{
	/* --- Read the thermally averaged table of van Hoof et al. (2014)
		   (gauntff.dat: comment lines start with '#'; the magic number
		   20140210; the numbers of points in log gamma^2 and log u;
		   the first log gamma^2, the first log u, the step in dex;
		   then the Gaunt factors, one row per log u. A second block
		   of the same size, the uncertainties, is not read.) Returns
		   0, or -1 with a message in errmsg. Replaces a table read
		   before. --                                    -------------- */

	char   line[GFF_MAX_LINE], *p, *end;
	int    Ng = 0, Nu = 0, n = 0, Ntot = 0, field = 0, lineNo = 0, status = -1;
	long   magic;
	double lg0 = 0.0, lu0 = 0.0, step = 0.0, v, *lng = NULL;
	FILE  *fp;

	if ((fp = fopen(fileName, "r")) == NULL) {
		snprintf(errmsg, len, "Unable to open free-free Gaunt factor table %s",
				 fileName);
		return -1;
	}
	while (fgets(line, sizeof(line), fp) != NULL) {
		lineNo++;
		if (strchr(line, '\n') == NULL && !feof(fp)) {
			snprintf(errmsg, len, "%s:%d: line too long", fileName, lineNo);
			goto done;
		}
		if (line[0] == '#') continue;
		switch (field) {
		case 0:
			if (sscanf(line, "%ld", &magic) != 1) continue; /* blank */
			if (magic != GFF_MAGIC) {
				snprintf(errmsg, len,
						 "%s:%d: magic number %ld, expected %d (the thermally "
						 "averaged table of van Hoof et al. 2014)",
						 fileName, lineNo, magic, GFF_MAGIC);
				goto done;
			}
			field++;
			continue;
		case 1:
			if (sscanf(line, "%d %d", &Ng, &Nu) != 2 || Ng < 4 || Nu < 4) {
				snprintf(errmsg, len,
						 "%s:%d: expected at least 4 x 4 grid points",
						 fileName, lineNo);
				goto done;
			}
			Ntot = Ng * Nu;
			lng  = (double *) malloc(Ntot * sizeof(double));
			field++;
			continue;
		case 2:
		case 3:
		case 4:
			if (sscanf(line, "%lf", &v) != 1 || !isfinite(v)) {
				snprintf(errmsg, len, "%s:%d: expected a number", fileName,
						 lineNo);
				goto done;
			}
			if (field == 2) lg0 = v;
			if (field == 3) lu0 = v;
			if (field == 4) step = v;
			if (field == 4 && !(step > 0.0)) {
				snprintf(errmsg, len, "%s:%d: step must be positive",
						 fileName, lineNo);
				goto done;
			}
			field++;
			continue;
		}
		for (p = line; n < Ntot; p = end) {
			v = strtod(p, &end);
			if (end == p) break;
			if (!(v > 0.0) || !isfinite(v)) {
				snprintf(errmsg, len,
						 "%s:%d: Gaunt factors must be positive and finite",
						 fileName, lineNo);
				goto done;
			}
			lng[n++] = log(v);
		}
		if (n == Ntot) break;
		while (isspace((unsigned char) *p)) p++;
		if (*p != '\0') {
			snprintf(errmsg, len, "%s:%d: not a number: %.20s", fileName,
					 lineNo, p);
			goto done;
		}
	}
	if (field < 5 || n < Ntot) {
		snprintf(errmsg, len, "%s: incomplete table (%d of %d Gaunt factors)",
				 fileName, n, Ntot);
		goto done;
	}
	Gaunt_ff_freeTable();
	gff_Ng   = Ng;
	gff_Nu   = Nu;
	gff_lg0  = lg0;
	gff_lu0  = lu0;
	gff_step = step;
	gff_lng  = lng;
	lng      = NULL;
	status   = 0;

done:
	fclose(fp);
	free(lng);
	return status;
}
/* ------- end ---------------------------- Gaunt_ff_readTable.c ---- */

/* ------- begin -------------------------- Gaunt_ff_freeTable.c ---- */

void Gaunt_ff_freeTable(void)
{
	free(gff_lng);
	gff_lng = NULL;
	gff_Ng = gff_Nu = 0;
}
/* ------- end ---------------------------- Gaunt_ff_freeTable.c ---- */

double Gaunt_ff(RHState *rh, double lambda, int charge, double T)
{
	/* --- <g_ff> at lambda [nm] and temperature T [K] for an ion of
		   charge Z. Needs the table (Gaunt_ff_readTable). -- ----- */

	const char routineName[] = "Gaunt_ff";
	int        i0, j0, i, j;
	double     x, y, wg[4], wu[4], sum = 0.0;

	if (gff_lng == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "free-free Gaunt factor table not read "
											  "(call hydrogenPrepare first)");
	}
	x = (log10(SQ(charge) * E_RYDBERG / (KBOLTZMANN * T)) - gff_lg0) / gff_step;
	y = (log10(HPLANCK * CLIGHT / (lambda * NM_TO_M * KBOLTZMANN * T)) -
		 gff_lu0) /
		gff_step;
	if (!(x >= 0.0 && x <= gff_Ng - 1 && y >= 0.0 && y <= gff_Nu - 1)) {
		Error(rh, ERROR_LEVEL_2, routineName, "lambda = %g nm, T = %g K, Z = %d outside the free-free "
											  "Gaunt factor table",
			  lambda, T, charge);
	}
	i0 = gff_stencil(x, gff_Ng, wg);
	j0 = gff_stencil(y, gff_Nu, wu);
	for (j = 0; j < 4; j++) {
		const double *row = gff_lng + (j0 + j) * gff_Ng + i0;
		double        r   = 0.0;
		for (i = 0; i < 4; i++) r += wg[i] * row[i];
		sum += wu[j] * r;
	}
	return exp(sum);
}
/* ------- end ---------------------------- Gaunt_ff.c -------------- */

/* ------- begin -------------------------- Hminus_bf.c ------------- */

/* --- H-minus bound-free (photodetachment) cross-section (N-35).

 From: B. M. McLaughlin, P. C. Stancil, H. R. Sadeghpour and
	   R. C. Forrey (2017), J. Phys. B 50, 114001, the merged and
	   scaled single-detachment cross-section of their supplementary
	   data (doi:10.1088/1361-6455/aa6c1f).

 Three parts:
  -- lambda_min <= lambda <= lambda_max: linear interpolation in a
	 table read at run time (keyword HMINUS_BF_DATA, default
	 Atoms/Hminus_bf.ecsv: 80-125 nm, with the n = 2, 3, 4
	 resonances).
  -- lambda < lambda_min: sigma = A E^p (1 + (Ec/E)^q)^-s, E in eV,
	 a fit to the source below 85 nm.
  -- lambda > lambda_max: the form of John (1988, A&A 193, 189),
	 sigma = lambda^3 x^(3/2) sum_n C_n x^(n/2), lambda in micron,
	 x = 1/lambda - 1/lambda0, refitted to the source above 120 nm.
 Both fits are in Mb (1.0E-22 m^2). The table must cover the range
 where neither fit holds (85-120 nm). The table and the fits are made
 by docs/opacity/hminus_bf_table.py.

 Note: The table is read once in hydrogenPrepare, so that Hminus_bf
	   can be called from several threads. -------------- */

#define HMBF_LAMBDA0   1643.907988 /* threshold [nm], E0 = 0.754204 eV */
#define HMBF_SHORT_MAX 85.0        /* short-wavelength fit holds below */
#define HMBF_LONG_MIN  120.0       /* long-wavelength fit holds above */
#define HMBF_JOIN_WARN 0.02        /* warn when a join steps more */
#define MBARN_TO_M2    1.0E-22
#define HMBF_MAX_LINE  4096

static const double Hmbf_short[5] = {71425.67213, -3.295011963, 39.10512916,
									 4.530454567, 0.1547112957};
static const double Hmbf_John[6]  = {149.5104613, 43.6244984, -99.21317025,
									 73.67718562, -26.18373202, 3.742690569};

static int     Hmbf_N      = 0;
static double *Hmbf_lambda = NULL, *Hmbf_sigma = NULL;

/* ------- begin -------------------------- Hmbf_short_fit.c -------- */

static double Hmbf_short_fit(double lambda)
{
	/* --- Fit below the table, in Mb. --              -------------- */

	double E = HPLANCK * CLIGHT / (NM_TO_M * lambda) / EV;

	return Hmbf_short[0] * pow(E, Hmbf_short[1]) *
		   pow(1.0 + pow(Hmbf_short[2] / E, Hmbf_short[3]), -Hmbf_short[4]);
}
/* ------- end ---------------------------- Hmbf_short_fit.c -------- */

/* ------- begin -------------------------- Hmbf_long_fit.c --------- */

static double Hmbf_long_fit(double lambda)
{
	/* --- Fit above the table, in Mb; zero at and beyond threshold.
		   --                                            -------------- */

	int    n;
	double lambda_mu = lambda / 1000.0, x, sqrt_x, sum = 0.0, xn = 1.0;

	if (lambda >= HMBF_LAMBDA0) return 0.0;
	x      = 1.0 / lambda_mu - 1000.0 / HMBF_LAMBDA0;
	sqrt_x = sqrt(x);
	for (n = 0; n < 6; n++) {
		sum += Hmbf_John[n] * xn;
		xn *= sqrt_x;
	}
	return CUBE(lambda_mu) * x * sqrt_x * sum;
}
/* ------- end ---------------------------- Hmbf_long_fit.c --------- */

/* ------- begin -------------------------- Hminus_bf_sigma.c ------- */

double Hminus_bf_sigma(double lambda)
{
	/* --- Cross-section [m^2] at lambda [nm]; zero outside
		   0 < lambda < HMBF_LAMBDA0. Needs the table
		   (Hminus_bf_readTable). --                     -------------- */

	int    lo, hi, mid;
	double w;

	if (lambda <= 0.0 || lambda >= HMBF_LAMBDA0) return 0.0;
	if (lambda < Hmbf_lambda[0]) return Hmbf_short_fit(lambda) * MBARN_TO_M2;
	if (lambda > Hmbf_lambda[Hmbf_N - 1])
		return Hmbf_long_fit(lambda) * MBARN_TO_M2;

	lo = 0;
	hi = Hmbf_N - 1;
	while (hi - lo > 1) {
		mid = (lo + hi) / 2;
		if (lambda >= Hmbf_lambda[mid])
			lo = mid;
		else
			hi = mid;
	}
	w = (lambda - Hmbf_lambda[lo]) / (Hmbf_lambda[hi] - Hmbf_lambda[lo]);
	return ((1.0 - w) * Hmbf_sigma[lo] + w * Hmbf_sigma[hi]) * MBARN_TO_M2;
}
/* ------- end ---------------------------- Hminus_bf_sigma.c ------- */

/* ------- begin -------------------------- Hmbf_column.c ----------- */

static bool_t Hmbf_column(const char *line, const char *name,
						  const char *unit)
{
	/* --- Whether an ECSV datatype line declares column name with
		   unit, e.g. "# - {name: lambda, unit: nm, ...}". -- ---- */

	char        key[64];
	const char *p;

	snprintf(key, sizeof(key), "name: %s,", name);
	if ((p = strstr(line, key)) == NULL) return FALSE;
	snprintf(key, sizeof(key), "unit: %s", unit);
	if ((p = strstr(p, key)) == NULL) return FALSE;
	p += strlen(key);
	return (*p == ',' || *p == '}');
}
/* ------- end ---------------------------- Hmbf_column.c ----------- */

/* ------- begin -------------------------- Hminus_bf_readTable.c --- */

int Hminus_bf_readTable(const char *fileName, char *errmsg, size_t len)
{
	/* --- Read the table of Hminus_bf from an ECSV file (the subset
		   written by docs/opacity/hminus_bf_table.py: space-delimited,
		   columns lambda [nm] and sigma [Mbarn]). Returns 0, or -1
		   with a message in errmsg. Replaces a table read before.
		   --                                            -------------- */

	char    line[HMBF_MAX_LINE], first[64], second[64], extra[2];
	bool_t  header = TRUE, has_lambda = FALSE, has_sigma = FALSE;
	int     N = 0, Nmax = 0, lineNo = 0, status = -1;
	double *lambda = NULL, *sigma = NULL, l, s;
	FILE   *fp;

	if ((fp = fopen(fileName, "r")) == NULL) {
		snprintf(errmsg, len, "Unable to open H-minus bound-free table %s",
				 fileName);
		return -1;
	}
	while (fgets(line, sizeof(line), fp) != NULL) {
		lineNo++;
		if (strchr(line, '\n') == NULL && !feof(fp)) {
			snprintf(errmsg, len, "%s:%d: line too long", fileName, lineNo);
			goto done;
		}
		if (lineNo == 1 && strncmp(line, "# %ECSV ", 8) != 0) {
			snprintf(errmsg, len, "%s: not an ECSV file (no \"# %%ECSV\")",
					 fileName);
			goto done;
		}
		if (line[0] == '#') {
			if (!header) continue;
			if (Hmbf_column(line, "lambda", "nm")) has_lambda = TRUE;
			if (Hmbf_column(line, "sigma", "Mbarn")) has_sigma = TRUE;
			if (strstr(line, "delimiter:") != NULL &&
				strstr(line, "delimiter: ' '") == NULL) {
				snprintf(errmsg, len,
						 "%s:%d: only space-delimited ECSV is supported",
						 fileName, lineNo);
				goto done;
			}
			continue;
		}
		if (sscanf(line, "%1s", extra) != 1) continue; /* blank line */

		if (header) {
			/* --- The column-name line ends the header -- -------------- */

			if (!has_lambda || !has_sigma ||
				sscanf(line, "%63s %63s %1s", first, second, extra) != 2 ||
				strcmp(first, "lambda") != 0 || strcmp(second, "sigma") != 0) {
				snprintf(errmsg, len,
						 "%s: need exactly the columns lambda [nm] and sigma "
						 "[Mbarn], in this order",
						 fileName);
				goto done;
			}
			header = FALSE;
			continue;
		}
		if (sscanf(line, "%lf %lf %1s", &l, &s, extra) != 2) {
			snprintf(errmsg, len, "%s:%d: expected two numbers", fileName,
					 lineNo);
			goto done;
		}
		if (s < 0.0 || (N > 0 && l <= lambda[N - 1])) {
			snprintf(errmsg, len,
					 "%s:%d: need sigma >= 0 and lambda strictly increasing",
					 fileName, lineNo);
			goto done;
		}
		if (N == Nmax) {
			Nmax   = (Nmax == 0) ? 64 : 2 * Nmax;
			lambda = (double *) realloc(lambda, Nmax * sizeof(double));
			sigma  = (double *) realloc(sigma, Nmax * sizeof(double));
		}
		lambda[N] = l;
		sigma[N]  = s;
		N++;
	}
	if (N < 2 || lambda[0] > HMBF_SHORT_MAX || lambda[N - 1] < HMBF_LONG_MIN ||
		lambda[N - 1] >= HMBF_LAMBDA0) {
		snprintf(errmsg, len,
				 "%s: the table must have at least two rows and cover "
				 "%.0f-%.0f nm, below threshold (%.2f nm)",
				 fileName, HMBF_SHORT_MAX, HMBF_LONG_MIN, HMBF_LAMBDA0);
		goto done;
	}
	Hminus_bf_freeTable();
	Hmbf_N      = N;
	Hmbf_lambda = lambda;
	Hmbf_sigma  = sigma;
	lambda = sigma = NULL;
	status         = 0;

done:
	fclose(fp);
	free(lambda);
	free(sigma);
	return status;
}
/* ------- end ---------------------------- Hminus_bf_readTable.c --- */

/* ------- begin -------------------------- Hminus_bf_freeTable.c --- */

void Hminus_bf_freeTable(void)
{
	free(Hmbf_lambda);
	free(Hmbf_sigma);
	Hmbf_lambda = Hmbf_sigma = NULL;
	Hmbf_N                   = 0;
}
/* ------- end ---------------------------- Hminus_bf_freeTable.c --- */

/* ------- begin -------------------------- Hminus_bf_joins.c ------- */

void Hminus_bf_joins(double *short_ratio, double *long_ratio)
{
	/* --- Ratio fit / table at both ends of the table. -- ------- */

	*short_ratio = Hmbf_short_fit(Hmbf_lambda[0]) / Hmbf_sigma[0];
	*long_ratio  = Hmbf_long_fit(Hmbf_lambda[Hmbf_N - 1]) / Hmbf_sigma[Hmbf_N - 1];
}
/* ------- end ---------------------------- Hminus_bf_joins.c ------- */

bool_t Hminus_bf(RHState *rh, double lambda, double *chi, double *eta)
{
	const char   routineName[] = "Hminus_bf";
	register int k;

	int    Nspace = rh->atmos.Nspace;
	double hc_kla, stimEmis, twohnu3_c2, alpha_bf;

	if (lambda <= 0.0 || lambda >= HMBF_LAMBDA0) return FALSE;

	if (Hmbf_lambda == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "H-minus bound-free table not read "
											  "(call hydrogenPrepare first)");
	}
	alpha_bf = Hminus_bf_sigma(lambda);

	hc_kla     = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M * lambda);
	twohnu3_c2 = (2.0 * HPLANCK * CLIGHT) / CUBE(NM_TO_M * lambda);

	for (k = 0; k < Nspace; k++) {
		stimEmis = exp(-hc_kla / rh->atmos.T[k]);
		chi[k]   = rh->atmos.nHmin[k] * (1.0 - stimEmis) * alpha_bf;
		eta[k]   = rh->atmos.nHmin[k] * twohnu3_c2 * stimEmis * alpha_bf;
	}

	return TRUE;
}
/* ------- end ---------------------------- Hminus_bf.c ------------- */

/* ------- begin -------------------------- Hminus_ff.c ------------- */

/* --- H-minus free-free absorption coefficient (N-45).

 From: T. L. John (1988), A&A 193, 189-192, Tables 3a and 3b: fits to
	   the R-matrix calculation of K. L. Bell and K. A. Berrington
	   (1987), J. Phys. B 20, 801-806 (theta = 0.5-3.6, lambda >=
	   182.3 nm), which they reproduce to 1.5 %.

	   k = 1.0E-29 sum_{n=1}^{6} theta^((n+1)/2) (A_n lambda^2 + B_n
		   + C_n/lambda + D_n/lambda^2 + E_n/lambda^3 + F_n/lambda^4)
	   in cm^4/dyn per H atom and unit electron pressure (stimulated
	   emission included), theta = 5040/T, lambda in micron.

  -- lambda >= 364.5 nm: Table 3a; 182.3 <= lambda < 364.5 nm:
	 Table 3b. At 364.5 nm they agree to 0.6 % for theta >= 0.5;
	 below, where both are extrapolated, they part (3a/3b = 0.97 at
	 theta = 0.3, 0.81 at 0.1).
  -- lambda < 182.3 nm: k(182.3 nm) (lambda / 182.3 nm)^2. Neither
	 source goes below 182.3 nm, and Table 3b cannot be extrapolated
	 (it turns negative below 150 nm at theta = 0.5).
  -- theta > 3.6 (T < 1400 K) is taken at 3.6: beyond, the fits stop
	 rising with theta (from theta = 4.4) and turn negative near
	 theta = 7. Below theta = 0.5 they are extrapolated; they stay
	 positive and monotonic in theta down to theta = 0.05.

 Replaces the table of Stilley & Callaway (1970) used below 9113 nm.
 Constant tables only, so Hminus_ff can be called from several
 threads. --                                         -------------- */

#define NJOHN           6
#define HMFF_LAMBDA_3A  364.5 /* [nm] Table 3a from here on */
#define HMFF_LAMBDA_MIN 182.3 /* [nm] lambda^2 tail below */
#define HMFF_THETA_MAX  3.6

/* --- Rows n = 1..6, columns A..F. --                 -------------- */

static const double John_3a[NJOHN][NJOHN] = {
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{2483.346, 285.827, -2054.291, 2827.776, -1341.537, 208.952},
	{-3449.889, -1158.382, 8746.523, -11485.632, 5303.609, -812.939},
	{2200.040, 2427.719, -13651.105, 16755.524, -7510.494, 1132.738},
	{-696.271, -1841.400, 8624.970, -10051.530, 4400.067, -655.020},
	{88.283, 444.517, -1863.864, 2095.288, -901.788, 132.985}};

static const double John_3b[NJOHN][NJOHN] = {
	{518.1021, -734.8666, 1021.1775, -479.0721, 93.1373, -6.4285},
	{473.2636, 1443.4137, -1977.3395, 922.3575, -178.9275, 12.3600},
	{-482.2089, -737.1616, 1096.8827, -521.1341, 101.7963, -7.0571},
	{115.5291, 169.6374, -245.6490, 114.2430, -21.9972, 1.5097},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

/* ------- begin -------------------------- Hmff_coeffs.c ----------- */

static double Hmff_coeffs(double lambda, double C_lambda[NJOHN])
{
	/* --- The wavelength-dependent coefficients C_n(lambda), in
		   units of 1.0E-32 m^5/J (1.0E-29 cm^4/dyn). Returns the
		   factor of the lambda^2 tail (1 from 182.3 nm on).
		   --                                            -------------- */

	const double(*T)[NJOHN] = (lambda >= HMFF_LAMBDA_3A) ? John_3a : John_3b;
	double lambda_mu, lambda_inv, scale = 1.0;
	int    n;

	if (lambda < HMFF_LAMBDA_MIN) {
		scale  = SQ(lambda / HMFF_LAMBDA_MIN);
		lambda = HMFF_LAMBDA_MIN;
	}
	lambda_mu  = lambda / MICRON_TO_NM;
	lambda_inv = 1.0 / lambda_mu;
	for (n = 0; n < NJOHN; n++) {
		C_lambda[n] = SQ(lambda_mu) * T[n][0] + T[n][1] +
					  lambda_inv *
						  (T[n][2] +
						   lambda_inv *
							   (T[n][3] +
								lambda_inv * (T[n][4] + lambda_inv * T[n][5])));
	}
	return scale;
}
/* ------- end ---------------------------- Hmff_coeffs.c ----------- */

/* ------- begin -------------------------- Hmff_sum.c -------------- */

static double Hmff_sum(const double C_lambda[NJOHN], double theta)
{
	/* --- sum_{n=1}^{6} theta^((n+1)/2) C_n, theta at most
		   HMFF_THETA_MAX. --                            -------------- */

	double sqrt_theta, theta_n, sum = 0.0;
	int    n;

	if (theta > HMFF_THETA_MAX) theta = HMFF_THETA_MAX;
	sqrt_theta = sqrt(theta);
	theta_n    = theta;
	for (n = 0; n < NJOHN; n++) {
		sum += theta_n * C_lambda[n];
		theta_n *= sqrt_theta;
	}
	return sum;
}
/* ------- end ---------------------------- Hmff_sum.c -------------- */

/* ------- begin -------------------------- Hminus_ff_kappa.c ------- */

double Hminus_ff_kappa(double lambda, double theta)
{
	/* --- Free-free absorption coefficient [m^5/J] per H atom and
		   unit electron pressure at lambda [nm] and theta = 5040/T;
		   zero for lambda <= 0. --                      -------------- */

	double C_lambda[NJOHN], scale;

	if (lambda <= 0.0) return 0.0;
	scale = Hmff_coeffs(lambda, C_lambda);
	return scale * Hmff_sum(C_lambda, theta) * 1.0E-32;
}
/* ------- end ---------------------------- Hminus_ff_kappa.c ------- */

bool_t Hminus_ff(RHState *rh, double lambda, double *chi)
{
	register int k;

	double C_lambda[NJOHN], scale, pe;

	if (lambda <= 0.0) {
		for (k = 0; k < rh->atmos.Nspace; k++) chi[k] = 0.0;
		return TRUE;
	}
	scale = Hmff_coeffs(lambda, C_lambda) * 1.0E-32;
	for (k = 0; k < rh->atmos.Nspace; k++) {
		pe     = rh->atmos.ne[k] * KBOLTZMANN * rh->atmos.T[k];
		chi[k] = rh->atmos.H->n[0][k] * pe * scale *
				 Hmff_sum(C_lambda, THETA0 / rh->atmos.T[k]);
	}
	return TRUE;
}
/* ------- end ---------------------------- Hminus_ff.c ------------- */

/* ------- begin -------------------------- H2minus_ff.c ------------ */

#define NFF_H2    19
#define NTHETA_H2 8

static const double thetaFF_H2m[NTHETA_H2] = {
	0.5, 0.8, 1.0, 1.2, 1.6, 2.0, 2.8, 3.6};

static void H2minus_ff_prepare(RHState *rh)
{
	register int k;

	int    index = 0;
	double theta;

	atomic_flag_clear(&rh->hydrogen.H2mff_beyond_told);
	rh->hydrogen.H2mff_theta_index = (double *) malloc(rh->atmos.Nspace * sizeof(double));
	for (k = 0; k < rh->atmos.Nspace; k++) {
		theta = THETA0 / rh->atmos.T[k];
		if (theta <= thetaFF_H2m[0])
			rh->hydrogen.H2mff_theta_index[k] = 0;
		else if (theta >= thetaFF_H2m[NTHETA_H2 - 1])
			rh->hydrogen.H2mff_theta_index[k] = NTHETA_H2 - 1;
		else {
			Hunt(NTHETA_H2, thetaFF_H2m, theta, &index);
			rh->hydrogen.H2mff_theta_index[k] = index + (theta - thetaFF_H2m[index]) /
															(thetaFF_H2m[index + 1] - thetaFF_H2m[index]);
		}
	}
}

bool_t H2minus_ff(RHState *rh, double lambda, double *chi)
{
	const char   routineName[] = "H2minus_ff";
	register int k;

	/* --- H2-minus Free-Free absorption coefficients (in units of
		   10E-29 m^5/J). Stimulated emission is included.

	 From: Bell, K. L., (1980) J. Phys. B13, 1859.
	 Also: R. Mathisen (1984), Master's thesis, Inst. Theor.
			Astroph., University of Oslo, p. 18

		   The fractional indices for atmospheric temperatures into
		   the theta table are stored in H2mff_theta_index by
		   hydrogenPrepare, so that this routine only reads them and
		   can be called from several threads.

		   Beyond the last tabulated wavelength kappa is extended as
		   lambda^2 at each theta: the fitted exponent between the last
		   two columns is 1.987 to 2.004 over theta = 0.5 to 3.6
		   (N-53). A message is written once per run. The zero
		   row at lambda = 0 interpolates linearly below 350.5 nm and is
		   left as it is. --                           -------------- */

	static const double lambdaFF[NFF_H2] = {
		0.0, 350.5, 414.2, 506.3, 569.6, 650.9, 759.4, 911.3,
		1139.1, 1518.8, 1822.6, 2278.3, 3037.7, 3645.2, 4556.5, 6075.3,
		9113.0, 11391.3, 15188.3};

	static const double kappaFF[NFF_H2 * NTHETA_H2] = {
		/* --- lambda =     0.0 [nm] --                    -------------- */
		0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00, 0.00e+00,
		0.00e+00, 0.00e+00, 0.00e+00,
		/* --- lambda =   350.5 [nm] --                    -------------- */
		4.17e-02, 6.10e-02, 7.34e-02, 8.59e-02, 1.11e-01,
		1.37e-01, 1.87e-01, 2.40e-01,
		/* --- lambda =   414.2 [nm] --                    -------------- */
		5.84e-02, 8.43e-02, 1.01e-01, 1.17e-01, 1.49e-01,
		1.82e-01, 2.49e-01, 3.16e-01,
		/* --- lambda =   506.3 [nm] --                    -------------- */
		8.70e-02, 1.24e-01, 1.46e-01, 1.67e-01, 2.10e-01,
		2.53e-01, 3.39e-01, 4.27e-01,
		/* --- lambda =   569.6 [nm] --                    -------------- */
		1.10e-01, 1.54e-01, 1.80e-01, 2.06e-01, 2.55e-01,
		3.05e-01, 4.06e-01, 5.07e-01,
		/* --- lambda =   650.9 [nm] --                    -------------- */
		1.43e-01, 1.98e-01, 2.30e-01, 2.59e-01, 3.17e-01,
		3.75e-01, 4.92e-01, 6.09e-01,
		/* --- lambda =   759.4 [nm] --                    -------------- */
		1.92e-01, 2.64e-01, 3.03e-01, 3.39e-01, 4.08e-01,
		4.76e-01, 6.13e-01, 7.51e-01,
		/* --- lambda =   911.3 [nm] --                    -------------- */
		2.73e-01, 3.71e-01, 4.22e-01, 4.67e-01, 5.52e-01,
		6.33e-01, 7.97e-01, 9.63e-01,
		/* --- lambda =  1139.1 [nm] --                    -------------- */
		4.20e-01, 5.64e-01, 6.35e-01, 6.97e-01, 8.06e-01,
		9.09e-01, 1.11e+00, 1.32e+00,
		/* --- lambda =  1518.8 [nm] --                    -------------- */
		7.36e-01, 9.75e-01, 1.09e+00, 1.18e+00, 1.34e+00,
		1.48e+00, 1.74e+00, 2.01e+00,
		/* --- lambda =  1822.6 [nm] --                    -------------- */
		1.05e+00, 1.39e+00, 1.54e+00, 1.66e+00, 1.87e+00,
		2.04e+00, 2.36e+00, 2.68e+00,
		/* --- lambda =  2278.3 [nm] --                    -------------- */
		1.63e+00, 2.14e+00, 2.36e+00, 2.55e+00, 2.84e+00,
		3.07e+00, 3.49e+00, 3.90e+00,
		/* --- lambda =  3037.7 [nm] --                    -------------- */
		2.89e+00, 3.76e+00, 4.14e+00, 4.44e+00, 4.91e+00,
		5.28e+00, 5.90e+00, 6.44e+00,
		/* --- lambda =  3645.2 [nm] --                    -------------- */
		4.15e+00, 5.38e+00, 5.92e+00, 6.35e+00, 6.99e+00,
		7.50e+00, 8.32e+00, 9.02e+00,
		/* --- lambda =  4556.5 [nm] --                    -------------- */
		6.47e+00, 8.37e+00, 9.20e+00, 9.84e+00, 1.08e+01,
		1.16e+01, 1.28e+01, 1.38e+01,
		/* --- lambda =  6075.3 [nm] --                    -------------- */
		1.15e+01, 1.48e+01, 1.63e+01, 1.74e+01, 1.91e+01,
		2.04e+01, 2.24e+01, 2.40e+01,
		/* --- lambda =  9113.0 [nm] --                    -------------- */
		2.58e+01, 3.33e+01, 3.65e+01, 3.90e+01, 4.27e+01,
		4.54e+01, 4.98e+01, 5.33e+01,
		/* --- lambda = 11391.3 [nm] --                    -------------- */
		4.03e+01, 5.20e+01, 5.70e+01, 6.08e+01, 6.65e+01,
		7.08e+01, 7.76e+01, 8.30e+01,
		/* --- lambda = 15188.3 [nm] --                    -------------- */
		7.16e+01, 9.23e+01, 1.01e+02, 1.08e+02, 1.18e+02,
		1.26e+02, 1.38e+02, 1.47e+02};

	int    Nspace = rh->atmos.Nspace, index = 0;
	double pe, lambda_index, kappa, *nH2, scale = 1.0;

	if (rh->hydrogen.H2mff_theta_index == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Index table not set up "
											  "(call hydrogenPrepare first)");
	}

	if (lambda >= lambdaFF[NFF_H2 - 1]) {
		lambda_index = NFF_H2 - 1;
		scale        = SQ(lambda / lambdaFF[NFF_H2 - 1]);
		if (lambda > lambdaFF[NFF_H2 - 1] && !atomic_flag_test_and_set(&rh->hydrogen.H2mff_beyond_told)) {
			Warning(rh, routineName, "H2- free-free table ends at %.1f nm: kappa extended as lambda^2 "
									 "beyond it (first at %.1f nm; message once per run)",
					lambdaFF[NFF_H2 - 1], lambda);
		}
	} else {
		Hunt(NFF_H2, lambdaFF, lambda, &index);
		lambda_index = index + (lambda - lambdaFF[index]) /
								   (lambdaFF[index + 1] - lambdaFF[index]);
	}

	nH2 = rh->atmos.H2->n;
	for (k = 0; k < Nspace; k++) {
		if (nH2[k] > 0.0) {
			pe     = rh->atmos.ne[k] * KBOLTZMANN * rh->atmos.T[k];
			kappa  = bilinear(NTHETA_H2, NFF_H2, kappaFF,
							  rh->hydrogen.H2mff_theta_index[k], lambda_index);
			chi[k] = (nH2[k] * 1.0E-29) * pe * (scale * kappa);
		} else
			chi[k] = 0.0;
	}
	return TRUE;
}
/* ------- end ---------------------------- H2minus_ff.c ------------ */

/* ------- begin -------------------------- H2plus_ff.c ------------- */

#define NFF_H2P   15
#define NTEMP_H2P 10

static const double tempFF_H2p[NTEMP_H2P] = {
	2.5E+03, 3.0E+03, 3.5E+03, 4.0E+03, 5.0E+03,
	6.0E+03, 7.0E+03, 8.0E+03, 1.0E+04, 1.2E+04};

static void H2plus_ff_prepare(RHState *rh)
{
	register int k;

	int    index = 0;
	double T;

	atomic_flag_clear(&rh->hydrogen.H2pff_beyond_told);
	rh->hydrogen.H2pff_temp_index = (double *) malloc(rh->atmos.Nspace * sizeof(double));
	for (k = 0; k < rh->atmos.Nspace; k++) {
		T = rh->atmos.T[k];
		if (T <= tempFF_H2p[0])
			rh->hydrogen.H2pff_temp_index[k] = 0;
		else if (T >= tempFF_H2p[NTEMP_H2P - 1])
			rh->hydrogen.H2pff_temp_index[k] = NTEMP_H2P - 1;
		else {
			Hunt(NTEMP_H2P, tempFF_H2p, T, &index);
			rh->hydrogen.H2pff_temp_index[k] = index + (T - tempFF_H2p[index]) /
														   (tempFF_H2p[index + 1] - tempFF_H2p[index]);
		}
	}
}

bool_t H2plus_ff(RHState *rh, double lambda, double *chi)
{
	const char   routineName[] = "H2plus_ff";
	register int k;

	/* --- H2+ Free-Free scattering coefficients in units of
		   1.0E-49 m^-1 / (H atom/m^3) / (proton/M^3). Stimulated emission
	   is included. This represents the following interaction:

		 H + H^+ + \nu ---> H + H^+

	 From: D. R. Bates (1952), MNRAS 112, 40-44
	 Also: R. Mathisen (1984), Master's thesis, Inst. Theor.
			Astroph., University of Oslo, p. 45

		   The fractional indices for atmospheric temperatures into
		   the temperature table are stored in H2pff_temp_index by
		   hydrogenPrepare, so that this routine only reads them and
		   can be called from several threads.

		   Beyond the last tabulated wavelength the last column is
		   held (kappa is flat or falling there, so no power law
		   applies). A message is written once per run (N-53). The
		   zero row at lambda = 0 interpolates linearly below 384.6 nm
		   and is left as it is. --                    -------------- */

	static const double lambdaFF[NFF_H2P] = {
		0.0, 384.6, 555.6, 833.3, 1111.1, 1428.6, 1666.7,
		2000.0, 2500.0, 2857.1, 3333.3, 4000.0, 5000.0, 6666.7, 10000.0};

	static const double kappaFF[NFF_H2P * NTEMP_H2P] = {
		/* --- lambda =      0.0 [nm] --                   -------------- */
		0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
		/* --- lambda =    384.6 [nm] --                   -------------- */
		0.46, 0.46, 0.42, 0.39, 0.36, 0.33, 0.32, 0.30, 0.27, 0.25,
		/* --- lambda =    555.6 [nm] --                   -------------- */
		0.70, 0.62, 0.59, 0.56, 0.51, 0.43, 0.41, 0.39, 0.35, 0.34,
		/* --- lambda =    833.3 [nm] --                   -------------- */
		0.92, 0.86, 0.80, 0.76, 0.70, 0.64, 0.59, 0.55, 0.48, 0.43,
		/* --- lambda =   1111.1 [nm] --                   -------------- */
		1.11, 1.04, 0.96, 0.91, 0.82, 0.74, 0.68, 0.62, 0.53, 0.46,
		/* --- lambda =   1428.6 [nm] --                   -------------- */
		1.26, 1.19, 1.09, 1.02, 0.90, 0.80, 0.72, 0.66, 0.55, 0.48,
		/* --- lambda =   1666.7 [nm] --                   -------------- */
		1.37, 1.25, 1.15, 1.07, 0.93, 0.83, 0.74, 0.67, 0.56, 0.49,
		/* --- lambda =   2000.0 [nm] --                   -------------- */
		1.44, 1.32, 1.21, 1.12, 0.97, 0.84, 0.75, 0.67, 0.56, 0.48,
		/* --- lambda =   2500.0 [nm] --                   -------------- */
		1.54, 1.39, 1.26, 1.15, 0.98, 0.85, 0.75, 0.67, 0.55, 0.46,
		/* --- lambda =   2857.1 [nm] --                   -------------- */
		1.58, 1.42, 1.27, 1.16, 0.98, 0.84, 0.74, 0.66, 0.54, 0.45,
		/* --- lambda =   3333.3 [nm] --                   -------------- */
		1.62, 1.43, 1.28, 1.15, 0.97, 0.83, 0.72, 0.64, 0.52, 0.44,
		/* --- lambda =   4000.0 [nm] --                   -------------- */
		1.63, 1.43, 1.27, 1.14, 0.95, 0.80, 0.70, 0.62, 0.50, 0.42,
		/* --- lambda =   5000.0 [nm] --                   -------------- */
		1.62, 1.40, 1.23, 1.10, 0.90, 0.77, 0.66, 0.59, 0.48, 0.39,
		/* --- lambda =   6666.7 [nm] --                   -------------- */
		1.55, 1.33, 1.16, 1.03, 0.84, 0.71, 0.60, 0.53, 0.43, 0.36,
		/* --- lambda =  10000.0 [nm] --                   -------------- */
		1.39, 1.18, 1.02, 0.90, 0.73, 0.60, 0.52, 0.46, 0.37, 0.31};

	int    Nspace = rh->atmos.Nspace, index = 0;
	double lambda_index, kappa, *np;

	if (rh->hydrogen.H2pff_temp_index == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Index table not set up "
											  "(call hydrogenPrepare first)");
	}

	if (lambda >= lambdaFF[NFF_H2P - 1]) {
		lambda_index = NFF_H2P - 1;
		if (lambda > lambdaFF[NFF_H2P - 1] && !atomic_flag_test_and_set(&rh->hydrogen.H2pff_beyond_told)) {
			Warning(rh, routineName, "H2+ free-free table ends at %.1f nm: last column held beyond it "
									 "(first at %.1f nm; message once per run)",
					lambdaFF[NFF_H2P - 1], lambda);
		}
	} else {
		Hunt(NFF_H2P, lambdaFF, lambda, &index);
		lambda_index = index + (lambda - lambdaFF[index]) /
								   (lambdaFF[index + 1] - lambdaFF[index]);
	}

	np = rh->atmos.H->n[rh->atmos.H->Nlevel - 1];
	for (k = 0; k < Nspace; k++) {
		kappa  = bilinear(NTEMP_H2P, NFF_H2P, kappaFF,
						  rh->hydrogen.H2pff_temp_index[k], lambda_index);
		chi[k] = (rh->atmos.H->n[0][k] * 1.0E-29) * (np[k] * 1.0E-20) * kappa;
	}

	return TRUE;
}
/* ------- end ---------------------------- H2plus_ff.c ------------- */

/* ------- begin -------------------------- Rayleigh_H2.c ----------- */

#define RAYLEIGH_H2_LIMIT 121.57
#define N_RAYLEIGH_H2     21

bool_t Rayleigh_H2(RHState *rh, double lambda, double *scatt)
{
	register int k;

	static const double a[3] = {8.779E+01, 1.323E+06, 2.245E+10};

	static const double lambdaRH2[N_RAYLEIGH_H2] = {
		121.57, 130.00, 140.00, 150.00, 160.00, 170.00, 185.46,
		186.27, 193.58, 199.05, 230.29, 237.91, 253.56, 275.36,
		296.81, 334.24, 404.77, 407.90, 435.96, 546.23, 632.80};

	static const double sigma[N_RAYLEIGH_H2] = {
		2.35E-06, 1.22E-06, 6.80E-07, 4.24E-07, 2.84E-07, 2.00E-07, 1.25E-07,
		1.22E-07, 1.00E-07, 8.70E-08, 4.29E-08, 3.68E-08, 2.75E-08, 1.89E-08,
		1.36E-08, 8.11E-09, 3.60E-09, 3.48E-09, 2.64E-09, 1.04E-09, 5.69E-10};

	/* --- Rayleigh scattering by H2 molecules. Cross-section is given
		   in in units of Mb, 1.0E-22 m^2.

	  See: G. A. Victor and A. Dalgarno (1969), J. Chem. Phys. 50, 2535
		   (for lambda <= 632.80 nm), and
		   S. P. Tarafdar and M. S. Vardya (1973), MNRAS 163, 261
	 Also: R. Mathisen (1984), Master's thesis, Inst. Theor.
		   Astroph., University of Oslo, p. 49
		   --                                            -------------- */

	/* --- Within the table, log sigma is interpolated linearly in log
		   lambda: sigma falls as lambda^-4 to lambda^-9 over its four
		   decades, and a power law between the nodes is exact where
		   linear interpolation in sigma was up to 14 % too high (N-48)
		   --                                              -------------- */

	bool_t hunt;
	double lambda2, sigma_RH2, *nH2, ln_lambda[N_RAYLEIGH_H2], ln_sigma[N_RAYLEIGH_H2], ln_l;

	nH2 = rh->atmos.H2->n;

	if (lambda >= RAYLEIGH_H2_LIMIT) {
		if (lambda <= lambdaRH2[N_RAYLEIGH_H2 - 1]) {
			for (k = 0; k < N_RAYLEIGH_H2; k++) {
				ln_lambda[k] = log(lambdaRH2[k]);
				ln_sigma[k]  = log(sigma[k]);
			}
			ln_l = log(lambda);
			Linear(N_RAYLEIGH_H2, ln_lambda, ln_sigma, 1, &ln_l, &sigma_RH2, hunt = FALSE);
			sigma_RH2 = exp(sigma_RH2);
		} else {
			lambda2   = 1.0 / SQ(lambda);
			sigma_RH2 = (a[0] + (a[1] + a[2] * lambda2) * lambda2) * SQ(lambda2);
		}
		sigma_RH2 *= MEGABARN_TO_M2;

		for (k = 0; k < rh->atmos.Nspace; k++)
			scatt[k] = sigma_RH2 * nH2[k];

		return TRUE;
	} else
		return FALSE;
}
/* ------- end ---------------------------- Rayleigh_H2.c ----------- */

/* ------- begin -------------------------- bilinear.c -------------- */

double bilinear(int Ncol, int Nrow, const double *f, double x, double y)
{
	int    i, j, i1, j1;
	double fx, fy;

	/* --- Bilinear interpolation of the function f on the fractional
		   indices x and y --                            -------------- */

	i  = (int) x;
	fx = x - i;
	if (i == Ncol - 1)
		i1 = i;
	else
		i1 = i + 1;
	j  = (int) y;
	fy = y - j;
	if (j == Nrow - 1)
		j1 = j;
	else
		j1 = j + 1;

	return (1.0 - fx) * (1.0 - fy) * f[j * Ncol + i] +
		   fx * (1.0 - fy) * f[j * Ncol + i1] +
		   (1.0 - fx) * fy * f[j1 * Ncol + i] +
		   fx * fy * f[j1 * Ncol + i1];
}
/* ------- end ---------------------------- bilinear.c -------------- */

/* ------- begin -------------------------- hydrogenPrepare.c ------- */

void hydrogenPrepare(RHState *rh)
{
	/* --- One-time setup for the opacity routines above. After this
		   they only read their tables, and can be called concurrently
		   for different wavelengths. The temperature indices depend on
		   atmos.T, which should not change until hydrogenCleanup. - - */

	const char routineName[] = "hydrogenPrepare";
	char       msg[MAX_MESSAGE_LENGTH];
	double     short_ratio, long_ratio;

	hydrogenCleanup(rh);

	if (gff_lng == NULL &&
		Gaunt_ff_readTable(rh->input.GauntFfData, msg, MAX_MESSAGE_LENGTH) != 0)
		Error(rh, ERROR_LEVEL_2, routineName, "%s", msg);
	if (Hmbf_lambda == NULL) {
		if (Hminus_bf_readTable(rh->input.HminusBfData, msg,
								MAX_MESSAGE_LENGTH) != 0)
			Error(rh, ERROR_LEVEL_2, routineName, "%s", msg);
		Hminus_bf_joins(&short_ratio, &long_ratio);
		if (fabs(short_ratio - 1.0) > HMBF_JOIN_WARN ||
			fabs(long_ratio - 1.0) > HMBF_JOIN_WARN) {
			Warning(rh, routineName, "H-minus bound-free table %s does not join its fits: "
									 "fit/table = %.4f at %.2f nm, %.4f at %.2f nm",
					rh->input.HminusBfData, short_ratio, Hmbf_lambda[0], long_ratio,
					Hmbf_lambda[Hmbf_N - 1]);
		}
	}
	H2minus_ff_prepare(rh);
	H2plus_ff_prepare(rh);
}
/* ------- end ---------------------------- hydrogenPrepare.c ------- */

/* ------- begin -------------------------- hydrogenCleanup.c ------- */

void hydrogenCleanup(RHState *rh)
{
	/* --- Free the temperature index tables. The H-minus bound-free
		   table, which does not depend on the atmosphere, is kept.
		   --                                            -------------- */

	free(rh->hydrogen.H2mff_theta_index);
	free(rh->hydrogen.H2pff_temp_index);
	rh->hydrogen.H2mff_theta_index = rh->hydrogen.H2pff_temp_index = NULL;
}
/* ------- end ---------------------------- hydrogenCleanup.c ------- */
