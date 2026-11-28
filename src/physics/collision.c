/* ------- file: -------------------------- collision.c -------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Sat Sep 19 15:54:07 2009 --

	   --------------------------                      ----------RH-- */

/* --- Reads collisional data from atomic data file and computes
	   collisional rates.

	   Format is similar to MULTI's routine GENCOL:

		 Input line can be either

		   TEMP  Nitem  T[0]     ...   T[Nitem-1]

		 to read in the temperature grid for collisional coefficients, or

		   KEYWORD  i1  i2   coeff[0]    ...   coeff[Nitem-1]

		 to read the coefficients for transitions between i1 and i2.
		 Multiple entries of the temperature grid are allowed, but
		 at least one entry with the correct number of grid points
		 has to precede an entry with coefficients.

		 Allowed keywords are:

		 Keyword    Transition type
		 ----------------------------------------------------

		 TEMP  -->  Temperature grid

		 OMEGA -->  Collisional de-excitation of ions by electrons
		 CE    -->  Collisional de-excitation of neutrals by electrons
		 CI    -->  Collisional ionization by electrons
		 CP    -->  Collisional de-excitation by protons

		 CH0   -->  Charge exchange of ion with neutral hydrogen
		 CH+   -->  Charge exchange of neutral with protons

		 END   -->  End of input data
		 ----------------------------------------------------

 Note: Unit of number density is m^-3.

	   Convention: C_ij = C[i][j] represents the
				   transition j --> i

	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "statistics.h"
#include "state.h"

#define COMMENT_CHAR "#"
#define min(a, b)    (((a) < (b)) ? (a) : (b))

/* --- Function prototypes --                          -------------- */

double                 ftwo(double x);
void                   atomnm(RHState *rh, int anr, char *cseq);
int                    atomnr(RHState *rh, char *ID);
RH_NODISCARD rh_status ar85cea(RHState *rh, int i, int j, int k, struct Atom *atom,
							   double *rate);
void                   rowcol(RHState *rh, int i, int *row, int *col);

/* --- Global variables --                             -------------- */

/* ------- begin ---------------------------rowcol.c ---------------- */

void rowcol(RHState *rh, int anr, int *row, int *col)
{
	/* --- Row (1 to 6) and column in the periodic table of an
		   iso-electronic sequence with anr electrons, 1 <= anr <= 86.
		   The loop used to cover only the first five of the six rows
		   of istart, which left row unset for 55 to 86 electrons -- -- */

	register int i;

	int istart[7] = {0, 2, 10, 18, 36, 54, 86};

	RH_ASSERT(rh, anr > istart[0] && anr <= istart[6], "rowcol",
			  "Iso-electronic sequence with %d electrons is outside the six rows"
			  " of the periodic table that SUMMERS covers",
			  anr);
	*row = *col = 0;
	for (i = 0; i < 6; i++) {
		if (anr > istart[i] && anr <= istart[i + 1]) {
			*row = i + 1;
			*col = anr - istart[i];
			return;
		}
	}
}
/* ------- end ---------------------------- rowcol.c ---------------- */

/* ------- begin -------------------------- fone.c ------------------ */

rh_status fone(RHState *rh, double x, double *f)
{

	/* --- Function f_1(x) = e^x E1(x)

	  Ref: Arnaud & Rothenflug, 1985, A&ASS, 60, 425

	  Above x = 50 the asymptotic series
	  f_1 = sum_k (-1)^k k! / x^(k+1), summed until the terms drop
	  below rounding (B-64). Its leading term 1/x alone is not enough:
	  the AR85-CDI fit cancels f_1 against terms of order x^2
	  (1 + x - x (2 + x) f_1 = 2/x^2 - 12/x^3 + ...), and 1/x made the
	  fit 300 to 750 times too large. The terms decrease while k < x,
	  so above 50 the sum reaches double precision within about 26
	  terms; below 50, E1 (A&S 5.1.56, 2e-8) leaves that bracket within
	  3e-4, and the two sides meet to 5e-5 at x = 50 --  ----------- */

	if (x <= 50.0) {
		double e1;

		RH_TRY(E1(rh, x, &e1));
		*f = exp(x) * e1;
		return RH_OK;
	}

	double term = 1.0 / x, sum = term;

	for (int k = 1; k < 40 && fabs(term) > 1.0E-17 * sum; k++) {
		term *= -k / x;
		sum += term;
	}
	*f = sum;
	return RH_OK;
}
/* ------- end ---------------------------- fone.c ------------------ */

/* ------- begin -------------------------- ftwo.c ------------------ */

#define BRK 4.0

double ftwo(double x)
{
	/* --- Function f_2

	  Ref: Arnaud & Rothenflug, 1985, A&ASS, 60, 425

		   Improved description when x < 4 from:
		   Hummer, 1983, jqsrt, 30 281
		   --                                              -------------- */

	register int i;

	double y;

	double p[15] = {1.0000e+00, 2.1658e+02, 2.0336e+04, 1.0911e+06, 3.7114e+07,
					8.3963e+08, 1.2889e+10, 1.3449e+11, 9.4002e+11, 4.2571e+12,
					1.1743e+13, 1.7549e+13, 1.0806e+13, 4.9776e+11, 0.0000};
	double q[15] = {1.0000e+00, 2.1958e+02, 2.0984e+04, 1.1517e+06, 4.0349e+07,
					9.4900e+08, 1.5345e+10, 1.7182e+11, 1.3249e+12, 6.9071e+12,
					2.3531e+13, 4.9432e+13, 5.7760e+13, 3.0225e+13, 3.3641e+12};

	double px, xfact, qx, gamma, f0x, count, fact, term;

	if (x > BRK) {

		px    = p[0];
		xfact = 1.0;
		for (i = 1; i < 15; i++) {
			xfact /= x;
			px += p[i] * xfact;
		}

		qx    = q[0];
		xfact = 1.0;
		for (i = 1; i < 15; i++) {
			xfact /= x;
			qx += q[i] * xfact;
		}
		y = px / (qx * SQ(x));
	} else {

		gamma = 0.5772156649;
		f0x   = SQ(PI) / 12.0;
		term  = 1.0;
		count = 0.0;
		fact  = 1.0;
		xfact = 1.0;

		while (fabs(term / f0x) > 1.0E-8) {
			count = count + 1.0;
			fact  = fact * count;
			xfact = xfact * (-x);
			term  = xfact / (SQ(count) * fact);

			f0x = f0x + term;
			if (count > 100.0) {
				/*ERROR CODE HERE */
			}
		}
		y = exp(x) * ((log(x) + gamma) * (log(x) + gamma) * 0.5 + f0x);
	}

	return y;
}
/* ------- end ----------------------------- ftwo.c ----------------- */

/* ------- begin --------------------------- atomnr.c --------------- */

void atomnm(RHState *rh, int anr, char *cseq)
{
	/* --- Returns element ID --                         -------------- */

	if (anr < 28) {
		strcpy(cseq, rh->atmos.elements[anr].ID);
	} else {
		strcpy(cseq, "  ");
	}
}
/* ------- end ----------------------------- atomnm.c --------------- */

/* ------- begin --------------------------- ar85cea.c -------------- */

rh_status ar85cea(RHState *rh, int i, int j, int k, struct Atom *atom, double *rate)
{

	/* --- Routine for computing collisional autoionization rates using
		   formalism from Arnaud and Rothenflug 1985, A&ASS, 60, 425

		   94-02-22  new routine: (Philip Judge)
		   96-03-07  modifications: (Philip Judge)
					 Bug fixed: cup initialized to zero
		   --                                              -------------- */

	char   cseq[ATOM_ID_WIDTH + 1];
	int    iz, ichrge, isoseq;
	double zz, cup, bkt, b, zeff, iea = 0.0, y, f1y, a, g;

	(void) j;

	/* --- Initialize output to zero --                  -------------- */

	y   = 0.0;
	f1y = 0.0;
	cup = 0.0;

	/* --- Find element --                               -------------- */

	iz = atomnr(rh, atom->ID) + 1;
	zz = iz;

	if (iz < 1 || iz > 92) {
		/* ERROR CODE HERE */
	}

	/* --- Find iso-electronic sequence --               -------------- */

	ichrge = atom->stage[i];
	isoseq = iz - ichrge;
	atomnm(rh, isoseq - 1, &cseq[0]);

	/* --- Temperature in eV --                          -------------- */

	bkt = KBOLTZMANN * rh->atmos.T[k] / EV;

	/* --- Lithium sequence --                           -------------- */

	if (!strcmp(cseq, "LI")) {

		iea  = 13.6 * (pow(zz - 0.835, 2) - 0.25 * pow(zz - 1.62, 2));
		b    = 1.0 / (1.0 + 2.0E-4 * pow(zz, 3));
		zeff = zz - 0.43;
		y    = iea / bkt;
		RH_TRY(fone(rh, y, &f1y));
		g = 2.22 * f1y + 0.67 * (1.0 - y * f1y) + 0.49 * y * f1y +
			1.2 * y * (1.0 - y * f1y);

		cup = (1.60E-07 * 1.2 * b) / (pow(zeff, 2) * sqrt(bkt)) * exp(-y) * g;

		/* --- Special cases --                            -------------- */

		if (!strcmp(atom->ID, "C")) {

			/* --- C IV - app a ar85 --                      -------------- */

			cup *= 0.6;
		} else if (!strcmp(atom->ID, "N")) {

			/* --- N V  - app a ar85 --                      -------------- */

			cup *= 0.8;
		} else if (!strcmp(atom->ID, "O")) {

			/* --- O VI - app a ar85 --                      -------------- */

			cup *= 1.25;
		}
	} else if (!strcmp(cseq, "NA")) {

		/* --- Sodium sequence --                          -------------- */

		if (iz <= 16) {

			iea = 26.0 * (zz - 10.);
			a   = 2.8E-17 * pow(zz - 11.0, -0.7);
			y   = iea / bkt;

			RH_TRY(fone(rh, y, &f1y));
			cup = 6.69E+7 * a * iea / sqrt(bkt) * exp(-y) * (1.0 - y * f1y);

		} else if (iz >= 18 && iz <= 28) {

			iea = 11.0 * (zz - 10.0) * sqrt(zz - 10.0);
			a   = 1.3E-14 * pow(zz - 10.0, -3.73);
			y   = iea / bkt;
			RH_TRY(fone(rh, y, &f1y));
			cup = 6.69E+7 * a * iea / sqrt(bkt) * exp(-y) *
				  (1.0 - 0.5 * (y - SQ(y) + SQ(y) * y * f1y));

		} else {

			cup = 0.0;
		}
	}

	/* --- Magnesium-sulfur sequences --                 --------------- */

	if (!strcmp(cseq, "MG") || !strcmp(cseq, "AL") || !strcmp(cseq, "SI") ||
		!strcmp(cseq, "P") || !strcmp(cseq, "S")) {

		if (!strcmp(cseq, "MG")) iea = 10.3 * pow(zz - 10.0, 1.52);
		if (!strcmp(cseq, "AL")) iea = 18.0 * pow(zz - 11.0, 1.33);
		if (!strcmp(cseq, "SI")) iea = 18.4 * pow(zz - 12.0, 1.36);
		if (!strcmp(cseq, "P")) iea = 23.7 * pow(zz - 13.0, 1.29);
		if (!strcmp(cseq, "S")) iea = 40.1 * pow(zz - 14.0, 1.1);

		a = 4.0E-13 / (SQ(zz) * iea);
		y = iea / bkt;
		RH_TRY(fone(rh, y, &f1y));
		cup = 6.69E+7 * a * iea / sqrt(bkt) * exp(-y) *
			  (1.0 - 0.5 * (y - SQ(y) + SQ(y) * y * f1y));
	}

	/* --- Special cases --                              -------------- */

	if (!strcmp(atom->ID, "CA") && ichrge == 0) {
		iea = 25.;
		a   = 6.0e-17;
		b   = 1.12;
		y   = iea / bkt;
		RH_TRY(fone(rh, y, &f1y));
		cup = 6.69E+7 * a * iea / sqrt(bkt) * exp(-y) * (1.0 + b * f1y);
	} else if (!strcmp(atom->ID, "CA") && ichrge == 1) {
		a   = 9.8e-17;
		iea = 25.0;
		b   = 1.12;
		y   = iea / bkt;
		RH_TRY(fone(rh, y, &f1y));
		cup = 6.69E+7 * a * iea / sqrt(bkt) * exp(-y) * (1.0 + b * f1y);
	} else if (!strcmp(atom->ID, "FE") && ichrge == 3) {
		a   = 1.8E-17;
		iea = 60.0;
		b   = 1.0;
		y   = iea / bkt;
		RH_TRY(fone(rh, y, &f1y));
		cup = 6.69e+7 * a * iea / sqrt(bkt) * exp(-y) * (1.0 + b * f1y);
	} else if (!strcmp(atom->ID, "FE") && ichrge == 4) {
		a   = 5.0E-18;
		iea = 73.0;
		b   = 1.0;
		y   = iea / bkt;
		RH_TRY(fone(rh, y, &f1y));
		cup = 6.69E+7 * a * iea / sqrt(bkt) * exp(-y) * (1.0 + b * f1y);
	}

	*rate = cup * CUBE(CM_TO_M);
	return RH_OK;
}
/* --------- end --------------------------- ar85cea.c -------------- */

/* ------- begin --------------------------- summers.c -------------- */

rh_status summers(RHState *rh, int i, int j, double nne, struct Atom *atom, double *factor)
{

	/* --- Density sensitive dielectronic recombination
		   22-Jun-1994 changes begin P.G.Judge

	The term adi may be multiplied by a density-sensitive factor
	if needed- this is crucial for Li and B-like ions colliding with
	impacting electrons.

	This simple formulation was derived from a study of the dependence of
	the dielectronic "bump" in the figures of Summers 1974
	(Appleton Laboratory internal memo), and fitting according to the
	parameter Ne / z^7

	This should be accurate to typically +/- 0.1 in log in regions
	where it matters.  Worse case is e.g. C like Neon where it underestimates
	density factor by maybe 0.25 in log.

	June 24, 2006 changes begin P.G.Judge
	original (pre MAY 2006) code
			 define rho = nne/ z^7 where z is charge on recombining ion
		   rho=10.^(alog10(nne) - 7.* alog10(charge))
		   rho0=2.e3
		   if(isos eq 'LI' or isos eq 'NA' or isos eq 'K') then rho0 = 3.e1
		   ne_factor = 1./(1. + rho/rho0)^0.14
		   print,'ne_factor',ne_factor

	June 24, 2006, more accurate version.
	get the row and column of the recombined ion's isoelectronic
	sequence  in the periodic table
	the following parameters mimic the tables 1-19 of
	H. Summers' Appleton Lab Report 367, 1974

	Mar 9 2012 Jorrit Leenaarts: stole routine from Phil Judge's
	DIPER IDL package
	--                                                   -------------- */

	char   cseq[ATOM_ID_WIDTH + 1];
	int    iz, isoseq, row, col;
	double zz, rho0, rhoq, x, beta;

	/* --- Find atomic number of element --              -------------- */

	iz = atomnr(rh, atom->ID) + 1;
	if (iz < 1 || iz > 92) {
		/*ERROR CODE HERE*/
	}

	/* --- Charge of recombining ion --                  ------------- */

	zz = atom->stage[j];

	/* --- Find iso-electronic sequence of recombined ion -- --------- */

	isoseq = iz - atom->stage[i];
	if (isoseq < 1 || isoseq > 86) {
		return Fail(rh, RH_E_DATA, "summers", "Atom %s, SUMMERS for levels %d and %d: the recombined ion has %d"
											  " electrons; the density factor covers 1 to 86",
					atom->ID, i, j, isoseq);
	}
	atomnm(rh, isoseq - 1, &cseq[0]);

	/* --- Row and column in periodic table --           ------------- */

	rowcol(rh, isoseq, &row, &col);

	rhoq = nne * CUBE(CM_TO_M) / pow(zz, 7);
	x    = (0.5 * zz + (col - 1.0)) * row / 3.0;
	beta = -0.2 / log(x + 2.71828);
	rho0 = 30.0 + 50.0 * x;

	*factor = pow(1.0 + rhoq / rho0, beta);
	return RH_OK;
}
/* ------- end ----------------------------- summers.c -------------- */

/* ------- begin --------------------------- atomnr.c --------------- */

int atomnr(RHState *rh, char *ID)
{
	/* --- Returns atomic number of element with name id -- ----------- */

	int i = 0;

	while (strcmp(rh->atmos.elements[i].ID, ID)) i++;

	return i;
}
/* ------- end ---------------------------- atomnr.c ---------------- */

/* ------- begin -------------------------- RateCheck.c ----------- */

/* --- Rates given by formulas are checked at read time for what the
	   coefficients alone decide (finite, no value that makes the formula
	   invalid), and while they are evaluated for the sign of the result,
	   which depends on T. The result is the contribution to the
	   collision matrix C (s^-1): a negative one is set to zero and
	   reported once per record, atom and run (N-34), one that is not
	   finite stops the run, as for the tabulated rates (N-32). A formula
	   that is used outside the range of its fit (AR85-CHP and AR85-CHH,
	   zero there) is reported the same way. -- ------------------- */

typedef struct {
	long   Nneg, Nout;
	double minVal, Tneg, Tout[2];
} RateCheck;

static void rateCheckInit(RateCheck *rc)
{
	rc->Nneg = rc->Nout = 0;
	rc->minVal          = 0.0;
	rc->Tneg            = 0.0;
	rc->Tout[0] = rc->Tout[1] = 0.0;
}

/* --- The records already reported for an atom, kept in the atom
	   (CollisionRate reads the records again whenever the electron
	   density changes, so the state cannot be in the record).
	   Freed with the atom, or by freeRateReports to start again -- -- */

struct RateReports {
	int    N, Nalloc;
	char **key;
};

void freeRateReports(Atom *atom)
{
	struct RateReports *rr = atom->rateReports;
	int                 n;

	if (rr == NULL) return;
	for (n = 0; n < rr->N; n++) free(rr->key[n]);
	free(rr->key);
	free(rr);
	atom->rateReports = NULL;
}

/* --- TRUE the first time that a record is reported for the atom -- */

static bool_t firstReport(Atom *atom, const char *keyword, int i1, int i2,
						  const char *what)
{
	struct RateReports *rr;
	char                key[96];
	int                 n;

	snprintf(key, sizeof(key), "%.20s %d %d %.20s", keyword, i1, i2, what);
	if (atom->rateReports == NULL)
		atom->rateReports = (struct RateReports *) callocOrAbort(1, sizeof(struct RateReports));
	rr = atom->rateReports;
	for (n = 0; n < rr->N; n++)
		if (!strcmp(rr->key[n], key)) return FALSE;
	if (rr->N == rr->Nalloc) {
		rr->Nalloc = (rr->Nalloc == 0) ? 8 : 2 * rr->Nalloc;
		rr->key    = (char **) reallocOrAbort(rr->key, rr->Nalloc * sizeof(char *));
	}
	rr->key[rr->N++] = strdupOrAbort(key);
	return TRUE;
}

/* --- Stop at the first value that is not finite, in a table
	   (T is the grid of the values) or among the coefficients of a
	   formula (T = NULL) -- --------------------------------------- */

static rh_status checkFiniteValues(RHState *rh, const Atom *atom, const char *keyword, int i1,
								   int i2, const double *v, int n, const double *T)
{
	int  m;
	char msg[MAX_MESSAGE_LENGTH];

	for (m = 0; m < n; m++) {
		if (isfinite(v[m])) continue;
		if (T != NULL)
			snprintf(msg, sizeof(msg),
					 "Atom %s: collision record %.20s %d %d has a non-finite"
					 " value (%g) at T = %g K.\n Fix or remove this record in"
					 " the atom file",
					 atom->ID, keyword, i1, i2, v[m], T[m]);
		else
			snprintf(msg, sizeof(msg),
					 "Atom %s: collision record %.20s %d %d has a non-finite"
					 " coefficient (%g, number %d).\n Fix or remove this record"
					 " in the atom file",
					 atom->ID, keyword, i1, i2, v[m], m + 1);
		return Fail(rh, RH_E_DATA, "CollisionRate", "%s", msg);
	}
	return RH_OK;
}

static rh_status rateCheckCoeffError(RHState *rh, const Atom *atom, const char *keyword, int i1,
									 int i2, const char *what, double value)
{
	char msg[MAX_MESSAGE_LENGTH];

	snprintf(msg, sizeof(msg),
			 "Atom %s: collision record %.20s %d %d: %.200s (%g).\n Fix or remove this"
			 " record in the atom file",
			 atom->ID, keyword, i1, i2, what, value);
	return Fail(rh, RH_E_DATA, "CollisionRate", "%s", msg);
}

/* --- A contribution *x (s^-1) to C at temperature T: fails when it is
	   not finite, sets a negative value to zero and counts it -- --- */

static rh_status rateChecked(RHState *rh, double *x, RateCheck *rc, const Atom *atom,
							 const char *keyword, int i1, int i2, double T)
{
	char msg[MAX_MESSAGE_LENGTH];

	if (!isfinite(*x)) {
		snprintf(msg, sizeof(msg),
				 "Atom %s: collision record %.20s %d %d gives a rate that is not"
				 " finite (%g s^-1) at T = %g K.\n Check its coefficients",
				 atom->ID, keyword, i1, i2, *x, T);
		return Fail(rh, RH_E_NUMERIC, "CollisionRate", "%s", msg);
	}
	if (*x < 0.0) {
		if (rc->Nneg == 0 || *x < rc->minVal) {
			rc->minVal = *x;
			rc->Tneg   = T;
		}
		rc->Nneg++;
		*x = 0.0;
	}
	return RH_OK;
}

static void rateCheckOutside(RateCheck *rc, double T)
{
	if (rc->Nout == 0 || T < rc->Tout[0]) rc->Tout[0] = T;
	if (rc->Nout == 0 || T > rc->Tout[1]) rc->Tout[1] = T;
	rc->Nout++;
}

static void rateCheckReport(RHState *rh, const RateCheck *rc, Atom *atom,
							const char *keyword, int i1, int i2, long Nspace,
							double T1, double T2)
{
	char msg[MAX_MESSAGE_LENGTH];

	if (rc->Nneg > 0 && firstReport(atom, keyword, i1, i2, "negative")) {
		snprintf(msg, sizeof(msg),
				 "Atom %s: collision record %.20s %d %d gives a negative"
				 " contribution to the collision rate C at %ld of %ld locations"
				 " (smallest %g s^-1 at T = %g K); set to zero there.\n Check"
				 " the coefficients of this record",
				 atom->ID, keyword, i1, i2, rc->Nneg, Nspace, rc->minVal,
				 rc->Tneg);
		Warning(rh, "CollisionRate", "%s", msg);
	}
	if (rc->Nout > 0 && firstReport(atom, keyword, i1, i2, "range")) {
		snprintf(msg, sizeof(msg),
				 "Atom %s: collision record %.20s %d %d is only defined for"
				 " T = %g to %g K; its rate is zero at %ld of %ld locations"
				 " (T = %g to %g K)",
				 atom->ID, keyword, i1, i2, T1, T2, rc->Nout, Nspace,
				 rc->Tout[0], rc->Tout[1]);
		Warning(rh, "CollisionRate", "%s", msg);
	}
}
/* --- Tables (OMEGA, CE, CI, CP, CH0, CH+, CH) hold their end value
	   outside their temperature range. How often that happens is
	   counted per table type over the records of an atom, and one line
	   per type is reported once per atom and run (N-53 d) -- ------- */

#define NTABLE_TYPES 7

static const char *const tableKeyword[NTABLE_TYPES] = {
	"OMEGA", "CE", "CI", "CP", "CH0", "CH+", "CH"};

typedef struct {
	long   Nrec, Nout, Nbelow, Nabove, Npoints;
	double Ttab[2], Tatm[2];
} TableRange;

static int tableType(const char *keyword)
{
	int n;

	for (n = 0; n < NTABLE_TYPES; n++)
		if (!strcmp(keyword, tableKeyword[n])) return n;
	return -1;
}

/* --- Count the locations of the atmosphere outside the table grid
	   T[0..N-1] (ascending) for one record -- --------------------- */

static void tableRangeCount(RHState *rh, TableRange *tr, const double *T, int N, long Nspace)
{
	long   k, nlo = 0, nhi = 0;
	double Tlo = HUGE_VAL, Thi = 0.0;

	tr->Nrec++;
	for (k = 0; k < Nspace; k++) {
		if (rh->atmos.T[k] < T[0]) {
			nlo++;
			Tlo = MIN(Tlo, rh->atmos.T[k]);
			Thi = MAX(Thi, rh->atmos.T[k]);
		} else if (rh->atmos.T[k] > T[N - 1]) {
			nhi++;
			Tlo = MIN(Tlo, rh->atmos.T[k]);
			Thi = MAX(Thi, rh->atmos.T[k]);
		}
	}
	if (nlo + nhi == 0) return;

	if (tr->Nout == 0) {
		tr->Ttab[0] = T[0];
		tr->Ttab[1] = T[N - 1];
		tr->Tatm[0] = Tlo;
		tr->Tatm[1] = Thi;
	} else {
		tr->Ttab[0] = MIN(tr->Ttab[0], T[0]);
		tr->Ttab[1] = MAX(tr->Ttab[1], T[N - 1]);
		tr->Tatm[0] = MIN(tr->Tatm[0], Tlo);
		tr->Tatm[1] = MAX(tr->Tatm[1], Thi);
	}
	tr->Nout++;
	if (nlo > 0) tr->Nbelow++;
	if (nhi > 0) tr->Nabove++;
	tr->Npoints += nlo + nhi;
}

static void tableRangeReport(RHState *rh, const TableRange *tr, Atom *atom, int type, long Nspace)
{
	char        msg[MAX_MESSAGE_LENGTH];
	const char *rule;

	if (tr->Nout == 0 || !firstReport(atom, tableKeyword[type], 0, 0, "table range"))
		return;
	if (type == 1)
		rule = "the end value of C is held, so the effective collision"
			   " strength grows as T and the rate as sqrt(T)";
	else if (type == 0)
		rule = "the collision strength is held";
	else
		rule = "the end value is held";
	snprintf(msg, sizeof(msg),
			 "Atom %s: %ld of %ld %s tables are used outside their temperature"
			 " range (%ld below, %ld above; %ld table-locations of %ld; table"
			 " grids span T = %g to %g K, atmosphere outside them at T = %g to"
			 " %g K); %s",
			 atom->ID, tr->Nout, tr->Nrec, tableKeyword[type], tr->Nbelow,
			 tr->Nabove, tr->Npoints, tr->Nrec * Nspace, tr->Ttab[0], tr->Ttab[1],
			 tr->Tatm[0], tr->Tatm[1], rule);
	Warning(rh, "CollisionRate", "%s", msg);
}

/* ------- end ---------------------------- RateCheck.c ----------- */

/* ------- begin -------------------------- CollisionRate.c --------- */

#define MSHELL 5

rh_status CollisionRate(RHState *rh, struct Atom *atom, FILE *fp_atom)
{
	const char   routineName[] = "CollisionRate";
	register int k, n, m, ii;

	char inputLine[MAX_LINE_SIZE], keyword[MAX_LINE_SIZE], *pointer,
		labelStr[MAX_LINE_SIZE];
	bool_t exit_on_EOF;
	int    nitem, i1 = 0, i2 = 0, i = 0, j = 0, ij, ji = 0,
			   Nlevel = atom->Nlevel, Nitem, Ntemp = 0, status;
	long       Nspace = rh->atmos.Nspace;
	double     dE, C0, *T, *coeff, *C, Cdown, Cup, gij, *np, xj, fac, fxj;
	PchipTable ctable;
	RateCheck  rc;
	TableRange tabrange[NTABLE_TYPES];
	rh_status  result = RH_OK;

	int      Ncoef = 0, Nrow = 0;
	double **cdi = NULL, **badi = NULL;
	double   acolsh, tcolsh, aradsh, xradsh, adish, bdish, t0sh, t1sh, summrs, tg, cdn, cup;
	double   ar85t1, ar85t2, ar85a, ar85b, ar85c, ar85d, t4;
	double   de, zz, betab, cbar, dekt, dekti, wlog, wb, sumscl, e1, nefac, f1xj;

	getCPU(rh, 3, TIME_START, NULL);

	memset(tabrange, 0, sizeof(tabrange));
	C0 = ((E_RYDBERG / sqrt(M_ELECTRON)) * PI * SQ(RBOHR)) *
		 sqrt(8.0 / (PI * KBOLTZMANN));

	/* --- atom->C is allocated in readAtom; this routine is called
		   again when the electron density changes, so reuse it -- --- */

	if (atom->C == NULL) atom->C = matrix_double(rh, SQ(Nlevel), Nspace);
	for (ij = 0; ij < SQ(Nlevel); ij++) {
		for (k = 0; k < Nspace; k++) {
			atom->C[ij][k] = 0.0;
		}
	}

	T = coeff = NULL;
	C         = (double *) mallocOrAbort(Nspace * sizeof(double));

	/* --- For safety, initialize to 1, since we don't check whether it
		   gets set later on --                          -------------- */

	sumscl = 1.0;

	while ((status = getLine(rh, fp_atom, COMMENT_CHAR,
							 inputLine, exit_on_EOF = FALSE)) != EOF) {
		strcpy(keyword, strtok(inputLine, " "));

		if (!strcmp(keyword, "TEMP")) {

			/* --- Read temperature grid --                  -------------- */

			Nitem = Ntemp = atoi(strtok(NULL, " "));
			T             = (double *) reallocOrAbort(T, Nitem * sizeof(double));
			for (n = 0, nitem = 0; n < Nitem; n++) {
				if ((pointer = strtok(NULL, " ")) == NULL) break;
				nitem += sscanf(pointer, "%lf", T + n);
			}
		} else if (!strcmp(keyword, "OMEGA") || !strcmp(keyword, "CE") ||
				   !strcmp(keyword, "CI") || !strcmp(keyword, "CP") ||
				   !strcmp(keyword, "CH0") || !strcmp(keyword, "CH+") ||
				   !strcmp(keyword, "CH")) {

			/* --- Read level indices and collision coefficients, one per
				   temperature of the last TEMP grid. (Nitem may have been
				   set by a keyword with a fixed number of items since.) -- */

			i1 = atoi(strtok(NULL, " "));
			i2 = atoi(strtok(NULL, " "));
			if (Ntemp <= 0) {
				result = Fail(rh, RH_E_DATA, routineName, "Atom %s: %s collision data for levels %d %d"
														  " before any TEMP line",
							  atom->ID, keyword, i1, i2);
				goto cleanup;
			}
			Nitem = Ntemp;
			coeff = (double *) reallocOrAbort(coeff, Nitem * sizeof(double));

			for (n = 0, nitem = 0; n < Nitem; n++) {
				if ((pointer = strtok(NULL, " ")) == NULL) break;
				nitem += sscanf(pointer, "%lf", coeff + n);
			}
			/* --- Transitions i -> j are stored at index ji, transitions
			   j -> i are stored under ij. --            -------------- */

			i  = MIN(i1, i2);
			j  = MAX(i1, i2);
			ij = i * Nlevel + j;
			ji = j * Nlevel + i;

		} else if (!strcmp(keyword, "AR85-CHP") || !strcmp(keyword, "AR85-CHH")) {

			i1 = atoi(strtok(NULL, " "));
			i2 = atoi(strtok(NULL, " "));

			Nitem = 6;
			coeff = (double *) reallocOrAbort(coeff, Nitem * sizeof(double));

			for (n = 0, nitem = 0; n < Nitem; n++) {
				if ((pointer = strtok(NULL, " ")) == NULL) break;
				nitem += sscanf(pointer, "%lf", coeff + n);
			}

			i  = MIN(i1, i2);
			j  = MAX(i1, i2);
			ij = i * Nlevel + j;
			ji = j * Nlevel + i;

		} else if (!strcmp(keyword, "AR85-CEA") ||
				   !strcmp(keyword, "BURGESS")) {

			i1 = atoi(strtok(NULL, " "));
			i2 = atoi(strtok(NULL, " "));

			Nitem    = 1;
			coeff    = (double *) reallocOrAbort(coeff, Nitem * sizeof(double));
			coeff[0] = atof(strtok(NULL, " "));
			nitem    = 1;

			i  = MIN(i1, i2);
			j  = MAX(i1, i2);
			ij = i * Nlevel + j;
			ji = j * Nlevel + i;

		} else if (!strcmp(keyword, "SHULL82")) {

			i1 = atoi(strtok(NULL, " "));
			i2 = atoi(strtok(NULL, " "));

			Nitem = 8;
			coeff = (double *) reallocOrAbort(coeff, Nitem * sizeof(double));

			for (n = 0, nitem = 0; n < Nitem; n++) {
				if ((pointer = strtok(NULL, " ")) == NULL) break;
				nitem += sscanf(pointer, "%lf", coeff + n);
			}

			i  = MIN(i1, i2);
			j  = MAX(i1, i2);
			ij = i * Nlevel + j;
			ji = j * Nlevel + i;

		} else if (!strcmp(keyword, "BADNELL")) {

			/* --- BADNELL recipe for dielectronic recombination:
				   Bhavna Rathore: 20 Jan 2014
				   --                                        -------------- */

			i1    = atoi(strtok(NULL, " "));
			i2    = atoi(strtok(NULL, " "));
			Ncoef = atoi(strtok(NULL, " "));

			if (Ncoef < 1) {
				result = Fail(rh, RH_E_DATA, routineName, "Atom %s: collision record BADNELL %d %d: the"
														  " number of coefficients must be positive (%d)",
							  atom->ID, i1, i2, Ncoef);
				goto cleanup;
			}
			Nrow  = 2;
			Nitem = Nrow * Ncoef;
			badi  = matrix_double(rh, Nrow, Ncoef);

			for (m = 0, nitem = 0; m < Nrow; m++) {
				status = getLine(rh, fp_atom, COMMENT_CHAR, inputLine,
								 exit_on_EOF = FALSE);

				badi[m][0] = atof(strtok(inputLine, " "));
				nitem++;
				for (n = 1; n < Ncoef; n++) {
					if ((pointer = strtok(NULL, " ")) == NULL) break;
					nitem += sscanf(pointer, "%lf", badi[m] + n);
				}
			}

			i  = MIN(i1, i2);
			j  = MAX(i1, i2);
			ij = i * Nlevel + j;
			ji = j * Nlevel + i;

		} else if (!strcmp(keyword, "SUMMERS")) {

			/* --- Switch for density dependent DR coefficent

				   Give default multiplication factor of summers density
				   dependence of dielectronic recombination:

			   sumscl = 0.0 means there is no density dependence
			   sumscl = 1.0 means full summers density dependence
				   --                                        -------------- */
			Nitem  = 1;
			sumscl = atof(strtok(NULL, " "));
			nitem  = 1;

			/* --- The factor 0 <= f <= 1 of the density dependence keeps
				   sumscl * summers + 1 - sumscl between 0 and 1 (N-34) -- */

			if (!(sumscl >= 0.0 && sumscl <= 1.0)) {
				result = Fail(rh, RH_E_DATA, routineName, "Atom %s: the scale factor of the density dependence"
														  " of the dielectronic recombination (SUMMERS) must be"
														  " between 0 and 1 (%g)",
							  atom->ID, sumscl);
				goto cleanup;
			}

		} else if (!strcmp(keyword, "AR85-CDI")) {

			i1   = atoi(strtok(NULL, " "));
			i2   = atoi(strtok(NULL, " "));
			Nrow = atoi(strtok(NULL, " "));

			if (Nrow < 1) {
				result = Fail(rh, RH_E_DATA, routineName, "Atom %s: collision record AR85-CDI %d %d: the"
														  " number of shells must be positive (%d)",
							  atom->ID, i1, i2, Nrow);
				goto cleanup;
			}
			if (Nrow > MSHELL) {
				result = Fail(rh, RH_E_DATA, routineName, "Nrow: %i greater than mshell %i",
							  Nrow, MSHELL);
				goto cleanup;
			}

			Nitem = Nrow * MSHELL;
			cdi   = matrix_double(rh, Nrow, MSHELL);

			for (m = 0, nitem = 0; m < Nrow; m++) {
				status = getLine(rh, fp_atom, COMMENT_CHAR, inputLine, exit_on_EOF = FALSE);

				cdi[m][0] = atof(strtok(inputLine, " "));
				nitem++;
				for (n = 1; n < MSHELL; n++) {
					if ((pointer = strtok(NULL, " ")) == NULL) break;
					nitem += sscanf(pointer, "%lf", cdi[m] + n);
				}
			}

			i  = MIN(i1, i2);
			j  = MAX(i1, i2);
			ij = i * Nlevel + j;
			ji = j * Nlevel + i;

		} else if (strstr(keyword, "END")) {
			break;
		} else {
			result = Fail(rh, RH_E_DATA, routineName, "Unknown keyword: !%s!", keyword);
			goto cleanup;
		}

		if (nitem != Nitem) {
			result = Fail(rh, RH_E_DATA, routineName, "\n Read %d, not %d items (keyword = %s)\n",
						  nitem, Nitem, keyword);
			goto cleanup;
		}
		/* --- Records with formulas: the coefficients must be finite and
			   have the signs that the formula needs (N-34). The sign of
			   the result may also depend on T, see rateChecked -- ----- */

		if (!strcmp(keyword, "SHULL82")) {
			RH_TRY_GOTO(result, cleanup, checkFiniteValues(rh, atom, keyword, i1, i2, coeff, Nitem, NULL));

			/* --- Stop only where the formula is invalid: a negative A
				   makes the collisional term negative at every T; the other
				   signs and temperatures matter when their term is there
				   (a negative arad or adi is left to the check of the
				   sum) -- --------------------------------------- */

			if (coeff[0] < 0.0)
				RH_TRY_GOTO(result, cleanup, rateCheckCoeffError(rh, atom, keyword, i1, i2, "the collisional ionization coefficient A must not be negative", coeff[0]));
			if (coeff[4] != 0.0 && coeff[6] < 0.0)
				RH_TRY_GOTO(result, cleanup, rateCheckCoeffError(rh, atom, keyword, i1, i2, "the temperature T0 of the dielectronic term must not be negative", coeff[6]));
			if (coeff[4] != 0.0 && coeff[5] != 0.0 && coeff[7] < 0.0)
				RH_TRY_GOTO(result, cleanup, rateCheckCoeffError(rh, atom, keyword, i1, i2, "the temperature T1 of the dielectronic term must not be negative", coeff[7]));
		} else if (!strcmp(keyword, "BADNELL")) {
			RH_TRY_GOTO(result, cleanup, checkFiniteValues(rh, atom, keyword, i1, i2, badi[0], Ncoef, NULL));
			RH_TRY_GOTO(result, cleanup, checkFiniteValues(rh, atom, keyword, i1, i2, badi[1], Ncoef, NULL));
			for (n = 0; n < Ncoef; n++) {
				if (badi[1][n] != 0.0 && badi[0][n] < 0.0)
					RH_TRY_GOTO(result, cleanup, rateCheckCoeffError(rh, atom, keyword, i1, i2, "an energy E_i (first row) with a coefficient must not be negative", badi[0][n]));
			}
		} else if (!strcmp(keyword, "AR85-CDI")) {
			for (m = 0; m < Nrow; m++) {
				RH_TRY_GOTO(result, cleanup, checkFiniteValues(rh, atom, keyword, i1, i2, cdi[m], MSHELL, NULL));

				/* --- A row of zeros (padding) is skipped in the
					   evaluation -- ----------------------------- */

				if (cdi[m][0] <= 0.0 && (cdi[m][1] != 0.0 || cdi[m][2] != 0.0 ||
										 cdi[m][3] != 0.0 || cdi[m][4] != 0.0))
					RH_TRY_GOTO(result, cleanup, rateCheckCoeffError(rh, atom, keyword, i1, i2, "the ionization energy of a shell (first column) must be positive", cdi[m][0]));
			}
		} else if (!strcmp(keyword, "AR85-CEA") || !strcmp(keyword, "BURGESS") ||
				   !strcmp(keyword, "AR85-CHP") || !strcmp(keyword, "AR85-CHH")) {
			RH_TRY_GOTO(result, cleanup, checkFiniteValues(rh, atom, keyword, i1, i2, coeff, Nitem, NULL));
		}

		/* --- End of the reading section. Now filling the collision matrix

			   Interpolation in temperature T for all spatial locations:
			   monotone cubic Hermite in log T and log of the tabulated
			   value (pchip.c), which keeps the rates positive and adds no
			   extrema; the end values outside the table (N-32). Before,
			   a cubic spline in T rang through tables with one steep or
			   bad value and made rates negative. --      -------------- */

		if (!strcmp(keyword, "OMEGA") || !strcmp(keyword, "CE") ||
			!strcmp(keyword, "CI") || !strcmp(keyword, "CP") ||
			!strcmp(keyword, "CH0") || !strcmp(keyword, "CH+") ||
			!strcmp(keyword, "CH")) {

			/* --- Tabulated rates must be finite: a record with Inf or NaN
				   (e.g. a failed van Regemorter evaluation when the file was
				   made) would otherwise crash the spline with no hint of
				   where the problem is --                   -------------- */

			RH_TRY_GOTO(result, cleanup, checkFiniteValues(rh, atom, keyword, i1, i2, coeff, Nitem, T));
			RH_TRY_GOTO(result, cleanup, pchipTableInit(rh, &ctable, Nitem, T, coeff, PCHIP_LOGLOG));
			pchipTableEval(&ctable, Nspace, rh->atmos.T, C);
			pchipTableFree(&ctable);
			tableRangeCount(rh, &tabrange[tableType(keyword)], T, Nitem, Nspace);
		}

		if (!strcmp(keyword, "OMEGA")) {

			/* --- Collisional excitation of ions --         -------------- */

			for (k = 0; k < Nspace; k++) {
				Cdown = C0 * rh->atmos.ne[k] * C[k] /
						(atom->g[j] * sqrt(rh->atmos.T[k]));
				atom->C[ij][k] += Cdown;
				atom->C[ji][k] += Cdown * atom->nstar[j][k] / atom->nstar[i][k];
			}
		} else if (!strcmp(keyword, "CE")) {

			/* --- Collisional excitation of neutrals --     -------------- */

			gij = atom->g[i] / atom->g[j];
			for (k = 0; k < Nspace; k++) {
				Cdown = C[k] * rh->atmos.ne[k] * gij * sqrt(rh->atmos.T[k]);
				atom->C[ij][k] += Cdown;
				atom->C[ji][k] += Cdown * atom->nstar[j][k] / atom->nstar[i][k];
			}
		} else if (!strcmp(keyword, "CI")) {

			/* --- Collisional ionization --                 -------------- */

			dE = atom->E[j] - atom->E[i];
			for (k = 0; k < Nspace; k++) {
				Cup = C[k] * rh->atmos.ne[k] *
					  exp(-dE / (KBOLTZMANN * rh->atmos.T[k])) * sqrt(rh->atmos.T[k]);
				atom->C[ji][k] += Cup;
				atom->C[ij][k] += Cup * atom->nstar[i][k] / atom->nstar[j][k];
			}
		} else if (!strcmp(keyword, "CP")) {

			/* --- Collisions with protons --                -------------- */

			np = rh->atmos.H->n[rh->atmos.H->Nlevel - 1];
			for (k = 0; k < Nspace; k++) {
				Cdown = np[k] * C[k];
				atom->C[ij][k] += Cdown;
				atom->C[ji][k] += Cdown * atom->nstar[j][k] / atom->nstar[i][k];
			}
		} else if (!strcmp(keyword, "CH")) {

			/* --- Collisions with neutral hydrogen --       -------------- */

			for (k = 0; k < Nspace; k++) {
				Cup = rh->atmos.H->n[0][k] * C[k];
				atom->C[ji][k] += Cup;
				atom->C[ij][k] += Cup * atom->nstar[i][k] / atom->nstar[j][k];
			}
		} else if (!strcmp(keyword, "CH0")) {

			/* --- Charge exchange with neutral hydrogen --  -------------- */

			for (k = 0; k < Nspace; k++)
				atom->C[ij][k] += rh->atmos.H->n[0][k] * C[k];

		} else if (!strcmp(keyword, "CH+")) {

			/* --- Charge exchange with protons --           -------------- */

			np = rh->atmos.H->n[rh->atmos.H->Nlevel - 1];
			for (k = 0; k < Nspace; k++)
				atom->C[ji][k] += np[k] * C[k];

		} else if (!strcmp(keyword, "SHULL82")) {

			acolsh = coeff[0];
			tcolsh = coeff[1];
			aradsh = coeff[2];
			xradsh = coeff[3];
			adish  = coeff[4];
			bdish  = coeff[5];
			t0sh   = coeff[6];
			t1sh   = coeff[7];

			rateCheckInit(&rc);
			for (k = 0; k < Nspace; k++) {

				RH_TRY_GOTO(result, cleanup, summers(rh, i, j, rh->atmos.ne[k], atom, &nefac));
				summrs = sumscl * nefac + (1.0 - sumscl);
				tg     = rh->atmos.T[k];

				cdn = aradsh * pow(tg / 1.E4, -xradsh) +
					  summrs * adish / tg / sqrt(tg) * exp(-t0sh / tg) *
						  (1.0 + bdish * (exp(-t1sh / tg)));

				/* --- No collisional ionization without A or T_col (N-34;
					   T_col = 0 divided by zero) -- ------------------- */

				cup = (acolsh != 0.0 && tcolsh != 0.0)
						  ? acolsh * sqrt(tg) * exp(-tcolsh / tg) /
								(1.0 + 0.1 * tg / tcolsh)
						  : 0.0;

				/* --- Convert coefficient from cm^3 s^-1 to m^3 s^-1 -- ---- */

				cdn *= rh->atmos.ne[k] * CUBE(CM_TO_M);
				cup *= rh->atmos.ne[k] * CUBE(CM_TO_M);

				/* --- The sum can be negative: the dielectronic term is,
					   for B < -1, at high T (N-34) -- ------------------ */

				RH_TRY_GOTO(result, cleanup, rateChecked(rh, &cdn, &rc, atom, keyword, i1, i2, tg));
				RH_TRY_GOTO(result, cleanup, rateChecked(rh, &cup, &rc, atom, keyword, i1, i2, tg));

				/* --- 3-body recombination (high density limit) -- -------- */

				cdn += cup * atom->nstar[i][k] / atom->nstar[j][k];

				atom->C[ij][k] += cdn;
				atom->C[ji][k] += cup;
			}
			rateCheckReport(rh, &rc, atom, keyword, i1, i2, Nspace, 0.0, 0.0);
		} else if (!strcmp(keyword, "BADNELL")) {

			/* --- Fit for dielectronic recombination from Badnell

			   Bhavna Rathore Jan-14

			   First line coefficients are the energies in K (ener in Chianti)
			   Second line coefficients are the coefficients (coef in Chianti)
				   --                                        -------------- */

			rateCheckInit(&rc);
			for (k = 0; k < Nspace; k++) {
				RH_TRY_GOTO(result, cleanup, summers(rh, i, j, rh->atmos.ne[k], atom, &nefac));
				summrs = sumscl * nefac + (1.0 - sumscl);
				tg     = rh->atmos.T[k];

				cdn = 0.0;
				for (ii = 0; ii < Ncoef; ii++) {
					cdn += badi[1][ii] * exp(-badi[0][ii] / tg);
				}
				cdn *= pow(tg, -1.5);

				/* --- Convert coefficient from cm^3 s^-1 to m^3 s^-1 -- ---- */

				cdn *= rh->atmos.ne[k] * summrs * CUBE(CM_TO_M);
				RH_TRY_GOTO(result, cleanup, rateChecked(rh, &cdn, &rc, atom, keyword, i1, i2, tg));
				cup = cdn * atom->nstar[j][k] / atom->nstar[i][k];

				/* --- 3-body recombination (high density limit) -- --------- */

				cdn += cup * atom->nstar[i][k] / atom->nstar[j][k];

				atom->C[ij][k] += cdn;
				atom->C[ji][k] += cup;
			}
			rateCheckReport(rh, &rc, atom, keyword, i1, i2, Nspace, 0.0, 0.0);
			freeMatrix(rh, (void **) badi);
			badi = NULL;

		} else if (!strcmp(keyword, "AR85-CDI")) {

			/* --- Direct collionisional ionization --       -------------- */

			rateCheckInit(&rc);
			for (k = 0; k < Nspace; k++) {
				cup = 0.0;
				tg  = rh->atmos.T[k];

				for (m = 0; m < Nrow; m++) {

					if (cdi[m][0] <= 0.0) continue; /* padding, see the read check */

					xj  = cdi[m][0] * EV / (KBOLTZMANN * tg);
					fac = exp(-xj) * sqrt(xj);

					RH_TRY_GOTO(result, cleanup, fone(rh, xj, &f1xj));
					fxj = cdi[m][1] + cdi[m][2] * (1.0 + xj) +
						  (cdi[m][3] - xj * (cdi[m][1] + cdi[m][2] * (2.0 + xj))) * f1xj +
						  cdi[m][4] * xj * ftwo(xj);

					fxj = fxj * fac;
					fac = 6.69E-7 / pow(cdi[m][0], 1.5);
					cup += fac * fxj * CUBE(CM_TO_M);
				}
				cup *= rh->atmos.ne[k];
				RH_TRY_GOTO(result, cleanup, rateChecked(rh, &cup, &rc, atom, keyword, i1, i2, tg));
				cdn = cup * atom->nstar[i][k] / atom->nstar[j][k];

				atom->C[ij][k] += cdn;
				atom->C[ji][k] += cup;
			}
			rateCheckReport(rh, &rc, atom, keyword, i1, i2, Nspace, 0.0, 0.0);
			freeMatrix(rh, (void **) cdi);
			cdi = NULL;

		} else if (!strcmp(keyword, "AR85-CEA")) {

			/* --- Autoionization --                         -------------- */

			rateCheckInit(&rc);
			for (k = 0; k < Nspace; k++) {
				RH_TRY_GOTO(result, cleanup, ar85cea(rh, i, j, k, atom, &fac));
				cup = coeff[0] * fac * rh->atmos.ne[k];
				RH_TRY_GOTO(result, cleanup, rateChecked(rh, &cup, &rc, atom, keyword, i1, i2, rh->atmos.T[k]));
				atom->C[ji][k] += cup;
			}
			rateCheckReport(rh, &rc, atom, keyword, i1, i2, Nspace, 0.0, 0.0);

		} else if (!strcmp(keyword, "AR85-CHP")) {

			/* --- Charge transfer with ionized hydrogen -- --------------- */

			ar85t1 = coeff[0];
			ar85t2 = coeff[1];
			ar85a  = coeff[2];
			ar85b  = coeff[3];
			ar85c  = coeff[4];
			ar85d  = coeff[5];

			/* --- The proton density is the continuum of the hydrogen
				   model, its last level (as for CP and CH+); it was
				   atmos.H->n[5], a bound level of a hydrogen atom with
				   more than six levels and beyond the array for fewer
				   (N-34) -- ------------------------------------ */

			np = rh->atmos.H->n[rh->atmos.H->Nlevel - 1];
			rateCheckInit(&rc);
			for (k = 0; k < Nspace; k++) {
				if (rh->atmos.T[k] >= ar85t1 && rh->atmos.T[k] <= ar85t2) {

					t4  = rh->atmos.T[k] / 1.0E4;
					cup = ar85a * 1e-9 * pow(t4, ar85b) * exp(-ar85c * t4) *
						  exp(-ar85d * EV / KBOLTZMANN / rh->atmos.T[k]) * np[k] *
						  CUBE(CM_TO_M);
					RH_TRY_GOTO(result, cleanup, rateChecked(rh, &cup, &rc, atom, keyword, i1, i2, rh->atmos.T[k]));
					atom->C[ji][k] += cup;
				} else {
					rateCheckOutside(&rc, rh->atmos.T[k]);
				}
			}
			rateCheckReport(rh, &rc, atom, keyword, i1, i2, Nspace, ar85t1, ar85t2);
		} else if (!strcmp(keyword, "AR85-CHH")) {

			/* --- Charge transfer with neutral hydrogen --  -------------- */

			ar85t1 = coeff[0];
			ar85t2 = coeff[1];
			ar85a  = coeff[2];
			ar85b  = coeff[3];
			ar85c  = coeff[4];
			ar85d  = coeff[5];

			rateCheckInit(&rc);
			for (k = 0; k < Nspace; k++) {
				if (rh->atmos.T[k] >= ar85t1 && rh->atmos.T[k] <= ar85t2) {

					t4  = rh->atmos.T[k] / 1.0E4;
					cdn = ar85a * 1E-9 * pow(t4, ar85b) *
						  (1.0 + ar85c * exp(ar85d * t4)) *
						  rh->atmos.H->n[0][k] * CUBE(CM_TO_M);
					RH_TRY_GOTO(result, cleanup, rateChecked(rh, &cdn, &rc, atom, keyword, i1, i2, rh->atmos.T[k]));

					atom->C[ij][k] += cdn;
				} else {
					rateCheckOutside(&rc, rh->atmos.T[k]);
				}
			}
			rateCheckReport(rh, &rc, atom, keyword, i1, i2, Nspace, ar85t1, ar85t2);
		} else if (!strcmp(keyword, "BURGESS")) {

			/* --- Electron impact ionzation following Burgess & Chidichimo 1982,
			   MNRAS, 203, 1269-1280
				   --                                        -------------- */

			de    = (atom->E[j] - atom->E[i]) / EV;
			zz    = atom->stage[i];
			betab = 0.25 * (sqrt((100.0 * zz + 91.0) / (4.0 * zz + 3.0)) - 5.0);
			cbar  = 2.3;

			rateCheckInit(&rc);
			for (k = 0; k < Nspace; k++) {
				dekt  = de * EV / (KBOLTZMANN * rh->atmos.T[k]);
				dekt  = MIN(500, dekt);
				dekti = 1.0 / dekt;
				wlog  = log(1.0 + dekti);
				wb    = pow(wlog, betab / (1.0 + dekti));
				RH_TRY_GOTO(result, cleanup, E1(rh, dekt, &e1));
				cup = 2.1715E-8 * cbar * pow(13.6 / de, 1.5) * sqrt(dekt) *
					  e1 * wb * rh->atmos.ne[k] * CUBE(CM_TO_M);

				/* --- Add fudge factor --                     -------------- */

				cup *= coeff[0];
				RH_TRY_GOTO(result, cleanup, rateChecked(rh, &cup, &rc, atom, keyword, i1, i2, rh->atmos.T[k]));
				cdn = cup * atom->nstar[i][k] / atom->nstar[j][k];

				atom->C[ji][k] += cup;
				atom->C[ij][k] += cdn;
			}
			rateCheckReport(rh, &rc, atom, keyword, i1, i2, Nspace, 0.0, 0.0);
		}
	}

	if (status == EOF) {
		result = Fail(rh, RH_E_DATA, routineName, "Reached end of datafile before all data was read");
		goto cleanup;
	}
	for (n = 0; n < NTABLE_TYPES; n++)
		tableRangeReport(rh, &tabrange[n], atom, n, Nspace);

	sprintf(labelStr, "Collision Rate %2s", atom->ID);
	getCPU(rh, 3, TIME_POLL, labelStr);

	/* --- Clean up, also after a failure --             -------------- */

cleanup:
	if (badi != NULL) freeMatrix(rh, (void **) badi);
	if (cdi != NULL) freeMatrix(rh, (void **) cdi);
	free(C);
	free(T);
	free(coeff);
	return result;
}
/* ------- end ---------------------------- CollisionRate.c --------- */
