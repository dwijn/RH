/* ------- file: -------------------------- gcheck.c ----------------

	   Acceptance checks of the map G (docs/SOLVER.md, S-2) for
	   rhf1d: sets up as rhf1d does, then runs gcheckRun() (gchecks.c).

	   usage: gcheck [--steps N] [--lag L] [--ecnoise H] [--necn M] [--fd] [--fdmax F]
					 [--perturb NAME]
					 [rhf1d options]

	   in a run directory whose keyword.input has SOLVER = FIXED_POINT.
	   Writes no output files.
	   --------------------------                      ----------RH-- */

#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "background.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "solver.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

static Geometry geometry_data;
static RHState  rh_state = RH_STATE_INIT(ONE_D_PLANE, &geometry_data);
RHState *const  rh       = &rh_state;

char messageStr[MAX_MESSAGE_LENGTH];

/* ------- begin -------------------------- gcheck.c ---------------- */

int main(int argc, char *argv[])
{
	bool_t        write_analyze_output, equilibria_only;
	GCheckOptions opt;

	gcheckParse(&argc, argv, &opt);
	setOptions(argc, argv);
	getCPU(0, TIME_START, NULL);
	SetFPEtraps();

	readInput();
	if (rh->input.solver != SOLVER_FIXED_POINT) {
		Error(ERROR_LEVEL_2, "gcheck",
			  "gcheck checks the map G of SOLVER = FIXED_POINT: set"
			  " SOLVER = FIXED_POINT\n in keyword.input, with the"
			  " keyword settings it requires");
	}
	rh->spectrum.updateJ = TRUE;

	getCPU(1, TIME_START, NULL);
	MULTIatmos(&rh->atmos, rh->geometry);
	if (rh->atmos.Stokes) Bproject();

	readAtomicModels();
	readMolecularModels();
	SortLambda();

	getBoundary(rh->geometry);

	Background(write_analyze_output = FALSE, equilibria_only = FALSE);
	convertScales(&rh->atmos, rh->geometry);

	getProfiles();
	initSolution();
	initScatter();

	getCPU(1, TIME_POLL, "Total initialize");

	gcheckRun(&opt);

	return ExitStatus();
}
/* ------- end ---------------------------- gcheck.c ---------------- */
