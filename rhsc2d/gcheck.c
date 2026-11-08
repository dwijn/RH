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

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

enum Topology topology = TWO_D_PLANE;

Atmosphere   atmos;
Geometry     geometry;
Spectrum     spectrum;
ProgramStats stats;
InputData    input;
CommandLine  commandline;
char         messageStr[MAX_LINE_SIZE];

/* ------- begin -------------------------- gcheck.c ---------------- */

int main(int argc, char *argv[])
{
	bool_t        analyze_output, equilibria_only;
	GCheckOptions opt;

	gcheckParse(&argc, argv, &opt);
	setOptions(argc, argv);
	getCPU(0, TIME_START, NULL);
	SetFPEtraps();

	readInput();
	if (input.solver != SOLVER_FIXED_POINT) {
		Error(ERROR_LEVEL_2, "gcheck",
			  "gcheck checks the map G of SOLVER = FIXED_POINT: set"
			  " SOLVER = FIXED_POINT\n in keyword.input, with the"
			  " keyword settings it requires");
	}
	spectrum.updateJ = TRUE;

	getCPU(1, TIME_START, NULL);
	readAtmos(&atmos, &geometry);
	if (atmos.Stokes) Bproject();
	fillMesh(&geometry);

	readAtomicModels();
	readMolecularModels();
	SortLambda();

	getBoundary(&atmos, &geometry);

	Background(analyze_output = FALSE, equilibria_only = FALSE);

	getProfiles();
	initSolution();
	initScatter();

	getCPU(1, TIME_POLL, "Total initialize");

	gcheckRun(&opt);

	return ExitStatus();
}
/* ------- end ---------------------------- gcheck.c ---------------- */
