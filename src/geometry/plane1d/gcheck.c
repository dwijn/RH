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

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

enum Topology topology = ONE_D_PLANE;

Atmosphere   atmos;
Geometry     geometry;
Spectrum     spectrum;
ProgramStats stats;
InputData    input;
CommandLine  commandline;
char         messageStr[MAX_MESSAGE_LENGTH];

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
	if (input.solver != SOLVER_FIXED_POINT) {
		Error(ERROR_LEVEL_2, "gcheck",
			  "gcheck checks the map G of SOLVER = FIXED_POINT: set"
			  " SOLVER = FIXED_POINT\n in keyword.input, with the"
			  " keyword settings it requires");
	}
	spectrum.updateJ = TRUE;

	getCPU(1, TIME_START, NULL);
	MULTIatmos(&atmos, &geometry);
	if (atmos.Stokes) Bproject();

	readAtomicModels();
	readMolecularModels();
	SortLambda();

	getBoundary(&geometry);

	Background(write_analyze_output = FALSE, equilibria_only = FALSE);
	convertScales(&atmos, &geometry);

	getProfiles();
	initSolution();
	initScatter();

	getCPU(1, TIME_POLL, "Total initialize");

	gcheckRun(&opt);

	return ExitStatus();
}
/* ------- end ---------------------------- gcheck.c ---------------- */
