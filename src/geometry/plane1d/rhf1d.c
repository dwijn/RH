/* ------- file: -------------------------- rhf1d.c -----------------

	   Version:       rh2.0, 1-D plane-parallel
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Fri May  1 08:25:50 2020 --

	   --------------------------                      ----------RH-- */

/* --- Main routine of 1D plane-parallel radiative transfer program.
	   MALI scheme formulated according to Rybicki & Hummer

  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
	   G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215

	   Formal solution is performed with Feautrier difference scheme
	   in static atmospheres, and with piecewise quadratic integration
	   in moving atmospheres.

	   --                                              -------------- */

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
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

static Geometry geometry_data;
static RHState  rh_state = RH_STATE_INIT(ONE_D_PLANE, &geometry_data);
RHState *const  rh       = &rh_state;

/* ------- begin -------------------------- rhf1d.c ----------------- */

int main(int argc, char *argv[])
{
	bool_t write_analyze_output, equilibria_only;
	int    niter, nact;
	double deltaJ;

	Atom     *atom;
	Molecule *molecule;

	/* --- Read input data and initialize --             -------------- */

	setOptions(argc, argv);
	getCPU(0, TIME_START, NULL);
	SetFPEtraps();

	readInput();
	rh->spectrum.updateJ = TRUE;

	getCPU(1, TIME_START, NULL);
	MULTIatmos(&rh->atmos, rh->geometry);
	if (rh->atmos.Stokes) Bproject();

	readAtomicModels();
	readMolecularModels();
	SortLambda();

	getBoundary(rh->geometry);

	Background(write_analyze_output = TRUE, equilibria_only = FALSE);
	convertScales(&rh->atmos, rh->geometry);

	getProfiles();
	initSolution();
	initScatter();

	getCPU(1, TIME_POLL, "Total Initialize");

	/* --- Solve radiative transfer for active ingredients -- --------- */

	Iterate(rh->input.NmaxIter, rh->input.iterLimit);

	/* --- If appropriate final solution(s) should be polarized -- ---- */

	adjustStokesMode();

	niter = 0;
	while (niter < rh->input.NmaxScatter) {
		deltaJ = solveSpectrum(FALSE, FALSE);
		if (!rh->input.backgr_pol && deltaJ <= rh->input.iterLimit) break;

		niter++;
	}

	/* --- Write output files --                         -------------- */

	if (rh->atmos.hydrostatic) {
		rh->geometry->scale = COLUMN_MASS;
		convertScales(&rh->atmos, rh->geometry);
	}
	getCPU(1, TIME_START, NULL);

	writeInput();
	writeAtmos(&rh->atmos);
	writeGeometry(rh->geometry);
	writeSpectrum(&rh->spectrum);
	writeFlux(FLUX_DOT_OUT);

	for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
		atom = rh->atmos.activeatoms[nact];

		writeAtom(atom);
		writePopulations(atom);
		writeRadRate(atom);
		writeCollisionRate(atom);
		writeDamping(atom);
	}
	for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
		molecule = rh->atmos.activemols[nact];
		writeMolPops(molecule);
	}

	writeOpacity();

	getCPU(1, TIME_POLL, "Write output");
	printWork();
	printTotalCPU();

	return ExitStatus();
}
/* ------- end ---------------------------- rhf1d.c ----------------- */
