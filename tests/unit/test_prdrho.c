/* ------- file: -------------------------- test_prdrho.c -----------

	   The check of rho of the PRD lines after a sub-iteration
	   (PRDRhoCheck in prdrho.c, N-4): rho <= 0 or not finite
	   makes a negative line emissivity.

	   For rho on a grid of wavelengths (and, with PRD_ANGLE_DEP, rays
	   and directions) and locations, the test checks

	   -- all rho positive and finite: valid, RH_OK, nothing recorded;
	   -- a zero, a negative value, a NaN and an infinity each: outside
		  a trial evaluation of G PRDRhoCheck fails with RH_E_NUMERIC;
		  inside one it returns RH_OK with valid FALSE,
		  negTauFound is TRUE until negTauReset and the message names
		  atom, line, wavelength index (and ray and direction),
		  location and value, and says that a non-positive rho gives a
		  negative line emissivity; rho is left as it was (no
		  clipping);
	   -- with several bad values the first one in memory order is
		  reported, and all of them counted.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "error.h"
#include "check.h"
#include "negtau.h"
#include "state.h"

#define NLAMBDA 4
#define NRAYS   2
#define NSPACE  3

static Atom       atom;
static AtomicLine line;
static double   **rho;
static long       Nrow;

static void reset_rho(void)
{
	for (long n = 0; n < Nrow * NSPACE; n++) rho[0][n] = 1.0 + 0.1 * n;
}

/* --- PRDRhoCheck outside a trial: its status and valid; the error
	   is cleared --                                    -------------- */

static rh_status outside_trial(bool_t *valid)
{
	rh_status status = PRDRhoCheck(rh, &atom, 0, valid);

	CHECK(status == RH_OK || (Failed(rh) && rh->error.code == status), "status %d not recorded",
		  status);
	ClearError(rh);
	return status;
}

/* --- Put value at row, location k (and a second one at the end if
	   second is set) and check the behaviour -- ------------------- */

static void check_case(const char *name, long row, int k, double value, bool_t second,
					   const char *where)
{
	char   msg[800];
	bool_t ok;
	double saved;

	reset_rho();
	rho[row][k] = value;
	if (second) rho[Nrow - 1][NSPACE - 1] = -1.0;
	saved = rho[row][k];

	CHECK(outside_trial(&ok) == RH_E_NUMERIC && !ok, "%s: outside a trial PRDRhoCheck does not fail",
		  name);

	negTauReset(rh);
	negTauTrial(rh, TRUE);
	CHECK(PRDRhoCheck(rh, &atom, 0, &ok) == RH_OK, "%s: in a trial: status", name);
	negTauTrial(rh, FALSE);
	CHECK(!ok && negTauFound(rh), "%s: in a trial: ok %d, invalid %d", name, ok,
		  negTauFound(rh));
	negTauMessage(rh, msg, sizeof(msg));
	CHECK(strstr(msg, "atom XX") && strstr(msg, "PRD line 0") && strstr(msg, where) &&
			  strstr(msg, "negative line emissivity"),
		  "%s: message: %s", name, msg);
	{
		char loc[32];

		sprintf(loc, "location %d ", k);
		CHECK(strstr(msg, loc), "%s: no '%s' in: %s", name, loc, msg);
	}
	if (second)
		CHECK(strstr(msg, "(2 of "), "%s: count of bad values: %s", name, msg);
	CHECK(rho[row][k] == saved || (isnan(saved) && isnan(rho[row][k])),
		  "%s: rho changed", name);
	negTauReset(rh);
	CHECK(!negTauFound(rh), "%s: still invalid after negTauReset", name);
}

static void check_all(const char *tag, bool_t angle)
{
	char   name[80], where[80];
	bool_t ok;

	rh->input.PRD_angle_dep = (angle) ? PRD_ANGLE_DEP : PRD_ANGLE_AVER;
	Nrow                    = (angle) ? 2L * NRAYS * NLAMBDA : NLAMBDA;
	rho                     = matrix_double(rh, Nrow, NSPACE);
	line.rho_prd            = rho;

	CHECK(PRDRhoSize(rh, &line) == Nrow * NSPACE, "%s: PRDRhoSize %ld, rows %ld", tag,
		  PRDRhoSize(rh, &line), Nrow);
	reset_rho();
	negTauReset(rh);
	negTauTrial(rh, TRUE);
	CHECK(PRDRhoCheck(rh, &atom, 0, &ok) == RH_OK && ok && !negTauFound(rh),
		  "%s: positive rho rejected", tag);
	negTauTrial(rh, FALSE);
	CHECK(outside_trial(&ok) == RH_OK && ok, "%s: positive rho fails outside a trial", tag);

	/* --- With angle dependence row = 2 * (NRAYS * la + mu) + to_obs -- */

	{
		long row = (angle) ? 2 * (NRAYS * 2 + 1) + 1 : 2; /* wavelength 2, ray 1, up */

		if (angle)
			sprintf(where, "wavelength index 2 (%.6f nm), ray 1 upward", line.lambda[2]);
		else
			sprintf(where, "wavelength index 2 (%.6f nm)", line.lambda[2]);
		sprintf(name, "%s zero", tag);
		check_case(name, row, 1, 0.0, FALSE, where);
		sprintf(name, "%s negative", tag);
		check_case(name, row, 1, -0.3, FALSE, where);
		sprintf(name, "%s NaN", tag);
		check_case(name, row, 1, NAN, FALSE, where);
		sprintf(name, "%s infinity", tag);
		check_case(name, row, 1, INFINITY, FALSE, where);
		sprintf(name, "%s two bad values", tag);
		check_case(name, row, 1, -0.3, TRUE, where);
	}
	freeMatrix(rh, (void **) rho);
}

int main(void)
{
	rh->commandline.logfile = fopen("/dev/null", "w");
	rh->atmos.Nspace        = NSPACE;
	rh->atmos.Nrays         = NRAYS;
	memset(&atom, 0, sizeof(atom));
	memset(&line, 0, sizeof(line));
	strcpy(atom.ID, "XX");
	atom.Nline   = 1;
	atom.line    = &line;
	line.i       = 0;
	line.j       = 1;
	line.PRD     = TRUE;
	line.Nlambda = NLAMBDA;
	line.lambda  = (double *) malloc(NLAMBDA * sizeof(double));
	for (int la = 0; la < NLAMBDA; la++) line.lambda[la] = 500.0 + 0.01 * la;

	check_all("angle independent", FALSE);
	check_all("angle dependent", TRUE);

	free(line.lambda);
	return check_summary("test_prdrho");
}
