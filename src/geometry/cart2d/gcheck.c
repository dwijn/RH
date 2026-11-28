/* ------- file: -------------------------- gcheck.c ----------------

	   Acceptance checks of the map G (docs/SOLVER.md, S-2) for
	   rhsc2d: sets up as rhsc2d does, then runs gcheckRun() (gchecks.c).

	   usage: gcheck [--steps N] [--lag L] [--ecnoise H] [--necn M] [--fd] [--fdmax F]
					 [--perturb NAME]
					 [rhsc2d options]

	   in a run directory whose keyword.input has SOLVER = FIXED_POINT.
	   Writes no output files.
	   --------------------------                      ----------RH-- */

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

static Geometry       geometry_data;
static RHState        rh_state = RH_STATE_INIT(TWO_D_PLANE, &geometry_data);
static RHState *const rh       = &rh_state;

/* ------- begin -------------------------- gcheck.c ---------------- */

int main(int argc, char *argv[])
{
	bool_t        analyze_output, equilibria_only;
	GCheckOptions opt;

	gcheckParse(&argc, argv, &opt);
	setOptions(rh, argc, argv);
	getCPU(rh, 0, TIME_START, NULL);
	SetFPEtraps(rh);

	readInput(rh);
	if (rh->input.solver != SOLVER_FIXED_POINT) {
		Error(rh, ERROR_LEVEL_2, "gcheck",
			  "gcheck checks the map G of SOLVER = FIXED_POINT: set"
			  " SOLVER = FIXED_POINT\n in keyword.input, with the"
			  " keyword settings it requires");
	}
	rh->spectrum.updateJ = TRUE;

	getCPU(rh, 1, TIME_START, NULL);
	readAtmos(rh, &rh->atmos, rh->geometry);
	if (rh->atmos.Stokes) Bproject(rh);
	fillMesh(rh, rh->geometry);

	readAtomicModels(rh);
	readMolecularModels(rh);
	ExitOnError(rh, SortLambda(rh));

	getBoundary(rh, &rh->atmos, rh->geometry);

	ExitOnError(rh, Background(rh, analyze_output = FALSE, equilibria_only = FALSE));

	ExitOnError(rh, getProfiles(rh));
	ExitOnError(rh, initSolution(rh));
	ExitOnError(rh, initScatter(rh));

	getCPU(rh, 1, TIME_POLL, "Total initialize");

	ExitOnError(rh, gcheckRun(rh, &opt));

	return ExitStatus(rh);
}
/* ------- end ---------------------------- gcheck.c ---------------- */
