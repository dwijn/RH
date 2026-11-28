/* ------- file: -------------------------- test_blends.c -----------

	   Symmetric lines with other wavelengths in their blue half
	   (asymmetricBlends, getlambda.c, N-55). Active atom with lines at
	   498, 500, 600, 400, 800, 900 and 1000 nm (symmetric, 21 points
	   on the red half), a continuum with its edge at 400.5 nm, an
	   inactive atom with a line at 799 nm, and the extra wavelengths
	   899 and 1000 nm:

	   -- 500 has the red half of 498 in its blue half; then 498 has the
		  new blue half of 500 in its own: 498 comes first, so this
		  needs a second pass;
	   -- 400 has continuum wavelengths in its blue half;
	   -- 900 has the extra wavelength 899 in its blue half;
	   -- 600 (nothing near), 800 (only an inactive atom's line near)
		  and 1000 (extra wavelength at line centre) stay symmetric;
	   -- a line made asymmetric has 2N - 1 points with the old N points
		  of its red half unchanged, bit for bit, and the mirror image
		  on the blue side;
	   -- a second call changes nothing.

	   The error paths of the wavelength grid (Q-55), in the process:
	   getLambda with qcore <= 0, and SortLambda with a wavelength
	   table that cannot be opened, is truncated or has a negative
	   count, return RH_E_DATA and record it.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "error.h"
#include "xdr.h"
#include "check.h"
#include "state.h"

#define NLINE   7
#define NHALF   21
#define NCONT_L 30

static const double lambda0[NLINE]     = {498.0, 500.0, 600.0, 400.0, 800.0, 900.0, 1000.0};
static const bool_t expect_asym[NLINE] = {TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE};

static void initLine(AtomicLine *line, double l0, int j, int i)
{
	memset(line, 0, sizeof(*line));
	line->lambda0    = l0;
	line->qcore      = 1.0;
	line->qwing      = 600.0;
	line->Nlambda    = NHALF;
	line->symmetric  = TRUE;
	line->Ncomponent = 1;
	line->c_shift    = (double *) calloc(1, sizeof(double));
	line->j          = j;
	line->i          = i;
	CHECK(getLambda(rh, line) == RH_OK, "getLambda of line %d -> %d failed", j, i);
}

/* --- The error paths (Q-55) --                     -------------- */

static void check_errors(void)
{
	AtomicLine bad;
	char       path[] = "test_blends_waveXXXXXX";
	int        fd, Nwave = 3;
	double     lambda = 500.0;
	FILE      *fp;
	XDR        xdrs;

	memset(&bad, 0, sizeof(bad));
	bad.qcore      = 0.0;
	bad.qwing      = 600.0;
	bad.Nlambda    = NHALF;
	bad.symmetric  = TRUE;
	bad.Ncomponent = 1;
	CHECK(getLambda(rh, &bad) == RH_E_DATA && Failed(rh) && rh->error.code == RH_E_DATA &&
			  bad.lambda == NULL,
		  "getLambda with qcore = 0: not RH_E_DATA");
	ClearError(rh);

	snprintf(rh->input.wavetable_input, sizeof(rh->input.wavetable_input), "%s",
			 "/nonexistent/wavetable");
	CHECK(SortLambda(rh) == RH_E_DATA && Failed(rh) && rh->error.code == RH_E_DATA,
		  "SortLambda, missing wavelength table: not RH_E_DATA");
	ClearError(rh);

	/* --- Three wavelengths announced, one written --  ------------ */

	fd = mkstemp(path);
	fp = (fd >= 0) ? fdopen(fd, "w") : NULL;
	CHECK(fp != NULL, "cannot write %s", path);
	if (fp == NULL) return;
	xdrstdio_create(&xdrs, fp, XDR_ENCODE);
	CHECK(xdr_int(&xdrs, &Nwave) && xdr_double(&xdrs, &lambda), "cannot write %s", path);
	xdr_destroy(&xdrs);
	fclose(fp);
	snprintf(rh->input.wavetable_input, sizeof(rh->input.wavetable_input), "%s", path);
	CHECK(SortLambda(rh) == RH_E_DATA && Failed(rh) && rh->error.code == RH_E_DATA,
		  "SortLambda, truncated wavelength table: not RH_E_DATA");
	ClearError(rh);

	/* --- A negative count --                         ------------ */

	fp    = fopen(path, "w");
	Nwave = -1;
	CHECK(fp != NULL, "cannot write %s", path);
	if (fp == NULL) {
		remove(path);
		return;
	}
	xdrstdio_create(&xdrs, fp, XDR_ENCODE);
	CHECK(xdr_int(&xdrs, &Nwave), "cannot write %s", path);
	xdr_destroy(&xdrs);
	fclose(fp);
	CHECK(SortLambda(rh) == RH_E_DATA && Failed(rh) && rh->error.code == RH_E_DATA,
		  "SortLambda, negative count in the wavelength table: not RH_E_DATA");
	ClearError(rh);
	remove(path);
}

int main(void)
{
	Atom            atoms[2];
	AtomicLine      lines[NLINE], passive;
	AtomicContinuum cont;
	double          red[NLINE][NHALF], extra[2] = {899.0, 1000.0};
	int             kr, la, Nchanged;

	rh->commandline.logfile = stderr;
	rh->commandline.quiet   = TRUE;
	rh->atmos.moving        = FALSE;
	rh->atmos.Stokes        = FALSE;
	rh->atmos.vmicro_char   = 3.0E3; /* q_to_lambda = 1e-5 lambda0: 3 nm at 500 */

	for (kr = 0; kr < NLINE; kr++) {
		initLine(&lines[kr], lambda0[kr], kr + 1, 0);
		CHECK(lines[kr].Nlambda == NHALF && lines[kr].lambda[0] == lambda0[kr],
			  "line %g: %d points from %g", lambda0[kr], lines[kr].Nlambda,
			  lines[kr].lambda[0]);
		memcpy(red[kr], lines[kr].lambda, NHALF * sizeof(double));
	}

	/* --- Continuum up to its edge at 400.5 nm, one point beyond it
		   (SortLambda puts that one at the edge) --     -------------- */

	memset(&cont, 0, sizeof(cont));
	cont.lambda0 = 400.5;
	cont.Nlambda = NCONT_L;
	cont.lambda  = (double *) malloc(NCONT_L * sizeof(double));
	for (la = 0; la < NCONT_L; la++) cont.lambda[la] = 300.0 + 3.5 * la; /* 398, then 401.5 */

	initLine(&passive, 799.0, 1, 0);

	memset(atoms, 0, sizeof(atoms));
	strcpy(atoms[0].ID, "XX");
	atoms[0].active    = TRUE;
	atoms[0].Nline     = NLINE;
	atoms[0].line      = lines;
	atoms[0].Ncont     = 1;
	atoms[0].continuum = &cont;
	strcpy(atoms[1].ID, "YY");
	atoms[1].active = FALSE;
	atoms[1].Nline  = 1;
	atoms[1].line   = &passive;

	rh->atmos.Natom     = 2;
	rh->atmos.atoms     = atoms;
	rh->atmos.Nmolecule = 0;

	CHECK(asymmetricBlends(rh, extra, 2, &Nchanged) == RH_OK, "asymmetricBlends failed");
	CHECK(Nchanged == 4, "%d lines made asymmetric, want 4", Nchanged);

	for (kr = 0; kr < NLINE; kr++) {
		AtomicLine *line = &lines[kr];
		bool_t      asym = !line->symmetric;

		CHECK(asym == expect_asym[kr], "line %g: asymmetric %d, want %d", lambda0[kr],
			  asym, expect_asym[kr]);
		if (!asym) {
			CHECK(line->Nlambda == NHALF && memcmp(line->lambda, red[kr], NHALF * sizeof(double)) == 0,
				  "line %g: symmetric grid changed", lambda0[kr]);
			continue;
		}
		CHECK(line->Nlambda == 2 * NHALF - 1, "line %g: %d points, want %d", lambda0[kr],
			  line->Nlambda, 2 * NHALF - 1);
		CHECK(memcmp(line->lambda + NHALF - 1, red[kr], NHALF * sizeof(double)) == 0,
			  "line %g: red half changed", lambda0[kr]);
		for (la = 1; la < NHALF; la++)
			CHECK_CLOSE(line->lambda0 - line->lambda[NHALF - 1 - la],
						line->lambda[NHALF - 1 + la] - line->lambda0, 0.0, 1e-12,
						"line %g: blue point %d not the mirror of the red one", lambda0[kr], la);
	}

	CHECK(asymmetricBlends(rh, extra, 2, &Nchanged) == RH_OK && Nchanged == 0,
		  "a second call changes lines");

	for (kr = 0; kr < NLINE; kr++) {
		free(lines[kr].lambda);
		free(lines[kr].c_shift);
	}
	free(passive.lambda);
	free(passive.c_shift);
	free(cont.lambda);

	check_errors();

	return check_summary("test_blends");
}
