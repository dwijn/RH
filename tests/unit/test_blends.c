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
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "check.h"

extern Atmosphere  atmos;
extern CommandLine commandline;

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
	getLambda(line);
}

int main(void)
{
	Atom            atoms[2];
	AtomicLine      lines[NLINE], passive;
	AtomicContinuum cont;
	double          red[NLINE][NHALF], extra[2] = {899.0, 1000.0};
	int             kr, la, Nchanged;

	commandline.logfile = stderr;
	commandline.quiet   = TRUE;
	atmos.moving        = FALSE;
	atmos.Stokes        = FALSE;
	atmos.vmicro_char   = 3.0E3; /* q_to_lambda = 1e-5 lambda0: 3 nm at 500 */

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

	atmos.Natom     = 2;
	atmos.atoms     = atoms;
	atmos.Nmolecule = 0;

	Nchanged = asymmetricBlends(extra, 2);
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

	CHECK(asymmetricBlends(extra, 2) == 0, "a second call changes lines");

	for (kr = 0; kr < NLINE; kr++) {
		free(lines[kr].lambda);
		free(lines[kr].c_shift);
	}
	free(passive.lambda);
	free(passive.c_shift);
	free(cont.lambda);

	return check_summary("test_blends");
}
