/* ------- file: -------------------------- rhsc2d.c ----------------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Sep 23 10:57:11 2025 --

	   --------------------------                      ----------RH-- */

/* --- Main routine of 2-D Cartesian radiative transfer program.
	   MALI scheme formulated according to Rybicki & Hummer

  See: G. B. Rybicki and D. G. Hummer 1991, A&A 245, p. 171-181
	   G. B. Rybicki and D. G. Hummer 1992, A&A 263, p. 209-215

	   Formal solution is performed with the method of
	   short characteristics.

  See: P. B. Kunasz and L. H. Auer 1988, JQSRT 39, p. 67

	   --                                              -------------- */

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

static Geometry       geometry_data;
static RHState        rh_state = RH_STATE_INIT(TWO_D_PLANE, &geometry_data);
static RHState *const rh       = &rh_state;

/* ------- begin -------------------------- rhsc2d.c ---------------- */

int main(int argc, char *argv[])
{
	bool_t analyze_output, equilibria_only;
	int    niter, nact;
	double deltaJ;

	Atom     *atom;
	Molecule *molecule;

	/* --- Read input data and initialize --             -------------- */

	setOptions(rh, argc, argv);
	getCPU(rh, 0, TIME_START, NULL);
	SetFPEtraps();

	readInput(rh);
	rh->spectrum.updateJ = TRUE;

	getCPU(rh, 1, TIME_START, NULL);
	readAtmos(rh, &rh->atmos, rh->geometry);
	if (rh->atmos.Stokes) Bproject(rh);
	fillMesh(rh, rh->geometry);

	readAtomicModels(rh);
	readMolecularModels(rh);
	SortLambda(rh);

	getBoundary(rh, &rh->atmos, rh->geometry);

	Background(rh, analyze_output = TRUE, equilibria_only = FALSE);

	getProfiles(rh);
	initSolution(rh);
	initScatter(rh);

	getCPU(rh, 1, TIME_POLL, "Total initialize");

	/* --- Solve radiative transfer for active ingredients -- --------- */

	Iterate(rh, rh->input.NmaxIter, rh->input.iterLimit);

	/* --- If appropriate final solution(s) should be polarized -- ---- */

	adjustStokesMode(rh);

	niter = 0;
	while (niter < rh->input.NmaxScatter) {
		deltaJ = solveSpectrum(rh, FALSE, FALSE);
		if (!rh->input.backgr_pol && deltaJ <= rh->input.iterLimit) break;

		niter++;
	}

	/* --- Write output files --                     ------------------ */

	getCPU(rh, 1, TIME_START, NULL);

	writeInput(rh);
	writeAtmos(rh, &rh->atmos);
	writeGeometry(rh, rh->geometry);
	writeSpectrum(rh, &rh->spectrum);
	writeFlux(rh, FLUX_DOT_OUT);

	for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
		atom = rh->atmos.activeatoms[nact];

		writeAtom(rh, atom);
		writePopulations(rh, atom);
		writeRadRate(rh, atom);
		writeCollisionRate(rh, atom);
		writeDamping(rh, atom);
	}
	for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
		molecule = rh->atmos.activemols[nact];
		writeMolPops(rh, molecule);
	}

	writeOpacity(rh);

	getCPU(rh, 1, TIME_POLL, "Write output");
	printWork(rh);
	printTotalCPU(rh);

	return ExitStatus(rh);
}
/* ------- end ---------------------------- rhsc2d.c ---------------- */
