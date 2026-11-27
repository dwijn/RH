/* ------- file: -------------------------- backgrcontr.c -----------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Mon Nov 22 15:17:25 2010 --

	   --------------------------                      ----------RH-- */

#include <stdio.h>
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
static RHState        rh_state = RH_STATE_INIT(TWO_D_PLANE, &geometry_data);
static RHState *const rh       = &rh_state;

/* ------- begin -------------------------- backgrcontr.c ----------- */

int main(int argc, char *argv[])
{
	register int n;

	int     Nlambda, Nread;
	double *lambda;
	FILE   *fp;

	/* --- Read input data and initialize --             -------------- */

	setOptions(rh, argc, argv);
	SetFPEtraps();

	readInput(rh);
	readAtmos(rh, &rh->atmos, rh->geometry);
	readAtomicModels(rh);
	readMolecularModels(rh);

	if ((fp = fopen(CONTR_INPUT_FILE, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Unable to open inputfile %s", CONTR_INPUT_FILE);
	}

	Nread = fscanf(fp, "%d", &Nlambda);
	if (Nread != 1) {
		Error(rh, ERROR_LEVEL_2, argv[0], "Unable to read Nlambda from inputfile %s",
			  CONTR_INPUT_FILE);
	}
	lambda = (double *) malloc(Nlambda * sizeof(double));
	for (n = 0; n < Nlambda; n++) {
		Nread = fscanf(fp, "%lf", &lambda[n]);
		if (Nread != 1) {
			Error(rh, ERROR_LEVEL_2, argv[0], "Unable to read wavelength %d from inputfile %s",
				  n, CONTR_INPUT_FILE);
		}
	}
	fclose(fp);

	backgrOpac(rh, Nlambda, lambda);

	free(lambda);
}
/* ------- end ---------------------------- backgrcontr.c ----------- */
