/* ------- file: -------------------------- backgrcontr.c -----------

	   Version:       rh2.0, 1-D plane-parallel
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Apr 22 09:45:42 2009 --

	   --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "background.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"
#include "state.h"

#define CONTR_INPUT_FILE "contribute.input"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

static Geometry       geometry_data;
static RHState        rh_state = RH_STATE_INIT(ONE_D_PLANE, &geometry_data);
static RHState *const rh       = &rh_state;

/* ------- begin -------------------------- backgrcontr.c ----------- */

int main(int argc, char *argv[])
{
	register int n;

	int     Nlambda, result;
	double *lambda;
	FILE   *fp;

	/* --- Read input data and initialize --             -------------- */

	setOptions(rh, argc, argv);
	SetFPEtraps();

	readInput(rh);
	MULTIatmos(rh, &rh->atmos, rh->geometry);
	readAtomicModels(rh);
	readMolecularModels(rh);

	if ((fp = fopen(CONTR_INPUT_FILE, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Unable to open inputfile %s", CONTR_INPUT_FILE);
	}

	result = fscanf(fp, "%d", &Nlambda);
	if (result != 1) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Unable to read Nlambda from inputfile %s",
			  CONTR_INPUT_FILE);
	}
	lambda = (double *) malloc(Nlambda * sizeof(double));
	for (n = 0; n < Nlambda; n++) {
		result = fscanf(fp, "%lf", &lambda[n]);
		if (result != 1) {
			Error(rh, ERROR_LEVEL_2, argv[0], "Unable to read wavelength %d from inputfile %s",
				  n, CONTR_INPUT_FILE);
		}
	}
	fclose(fp);

	backgrOpac(rh, Nlambda, lambda);

	free(lambda);
}
/* ------- end ---------------------------- backgrcontr.c ----------- */
