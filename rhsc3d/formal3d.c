/* ------- file: -------------------------- formal2d.c -------------- */

/* Formal solution in 3-D geometry of just one ray.
 *
 * Han Uitenbroek
 * Last modified: Wed Apr 22 09:50:30 2009 --
 */

#include <stdlib.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"
#include "constant.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

enum Topology topology = THREE_D_PLANE;

Atmosphere   atmos;
Spectrum     spectrum;
Geometry     geometry;
ProgramStats stats;
CommandLine  commandline;
InputData    input;
char         messageStr[MAX_LINE_SIZE];

/* ------- begin --------------------------             ------------- */

int main(int argc, char *argv[])
{
	register int k;

	bool_t to_observer;
	int    mu, Nread, Nwrite, Nspace, Nspect = 1, nspect, Nx, Ny, Nz,
								   Nplane;
	double *chi, *S, *I, *Psi;

	(void) argc;
	(void) argv;

	stats.printCPU      = TRUE;
	commandline.quiet   = FALSE;
	commandline.logfile = stderr;
	input.Eddington     = FALSE;
	input.Nthreads      = 0;

	getCPU(0, TIME_START, NULL);
	SetFPEtraps();

	input.interpolate_3D = BICUBIC_3D;

	geometry.x_boundary               = PERIODIC;
	geometry.y_boundary               = PERIODIC;
	geometry.z_boundary_value[TOP]    = IRRADIATED;
	geometry.z_boundary_value[BOTTOM] = IRRADIATED;

	Nread = fread(&atmos.angleSet.set, sizeof(int), 1, stdin);
	if (Nread != 1) {
		sprintf(messageStr, "Unable to read atmos.angleSet.set from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	Nread = fread(&mu, sizeof(int), 1, stdin);
	if (Nread != 1) {
		sprintf(messageStr, "Unable to read mu from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	if (atmos.angleSet.set == SET_GL) {
		Nread = fread(&atmos.angleSet.Ninclination, sizeof(int), 1, stdin);
		if (Nread != 1) {
			sprintf(messageStr, "Unable to read atmos.angleSet.Ninclination from standard input");
			Error(ERROR_LEVEL_2, argv[0], messageStr);
		}
		Nread = fread(&atmos.angleSet.Nazimuth, sizeof(int), 1, stdin);
		if (Nread != 1) {
			sprintf(messageStr, "Unable to read atmos.angleSet.Nazimuth from standard input");
			Error(ERROR_LEVEL_2, argv[0], messageStr);
		}
	}

	Nread = fread(&geometry.Nx, sizeof(int), 1, stdin);
	if (Nread != 1) {
		sprintf(messageStr, "Unable to read geometry.Nx from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	Nx    = geometry.Nx;
	Nread = fread(&geometry.Ny, sizeof(int), 1, stdin);
	if (Nread != 1) {
		sprintf(messageStr, "Unable to read geometry.Ny from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	Ny    = geometry.Ny;
	Nread = fread(&geometry.Nz, sizeof(int), 1, stdin);
	if (Nread != 1) {
		sprintf(messageStr, "Unable to read geometry.Nz from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	Nz     = geometry.Nz;
	Nspace = atmos.Nspace = geometry.Nx * geometry.Ny * geometry.Nz;
	geometry.Nplane = Nplane = Nx * Ny;

	/* --- Get increments in x, store and check for monotonicity -- --- */

	Nread = fread(&geometry.dx, sizeof(double), 1, stdin);
	if (Nread != 1) {
		sprintf(messageStr, "Unable to read geometry.dx from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	Nread = fread(&geometry.dy, sizeof(double), 1, stdin);
	if (Nread != 1) {
		sprintf(messageStr, "Unable to read geometry.dy from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	geometry.dx *= KM_TO_M;
	geometry.dy *= KM_TO_M;

	/* --- Get vertical grid --                          -------------- */

	geometry.z = (double *) malloc(Nz * sizeof(double));
	Nread      = fread(geometry.z, sizeof(double), Nz, stdin);
	if (Nread != Nz) {
		sprintf(messageStr, "Unable to read geometry.z from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	for (k = 0; k < Nz; k++) geometry.z[k] *= KM_TO_M;

	geometry.Itop    = matrix_double(Nspect, geometry.Nplane);
	geometry.Ibottom = matrix_double(Nspect, geometry.Nplane);
	Nread            = fread(geometry.Itop[0], sizeof(double), geometry.Nplane, stdin);
	if (Nread != geometry.Nplane) {
		sprintf(messageStr, "Unable to read geometry.Itop[0] from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	Nread = fread(geometry.Ibottom[0], sizeof(double), geometry.Nplane, stdin);
	if (Nread != geometry.Nplane) {
		sprintf(messageStr, "Unable to read geometry.Ibottom[0] from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}

	chi   = (double *) malloc(Nspace * sizeof(double));
	S     = (double *) malloc(Nspace * sizeof(double));
	Nread = fread(chi, sizeof(double), Nspace, stdin);
	if (Nread != Nspace) {
		sprintf(messageStr, "Unable to read chi from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	Nread = fread(S, sizeof(double), Nspace, stdin);
	if (Nread != Nspace) {
		sprintf(messageStr, "Unable to read S from standard input");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}

	I   = (double *) malloc(Nspace * sizeof(double));
	Psi = (double *) malloc(Nspace * sizeof(double));

	getAngleQuadr(&geometry);
	fillMesh(&geometry);
	getCPU(1, TIME_START, NULL);
	ShortChar(&geometry, nspect = 0, mu, to_observer = FALSE, chi, S, I, Psi);
	getCPU(1, TIME_POLL, "Formal Solution down");
	Nwrite = fwrite(I, sizeof(double), Nspace, stdout);
	if (Nwrite != Nspace) {
		sprintf(messageStr, "Unable to write I to standard output");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	Nwrite = fwrite(Psi, sizeof(double), Nspace, stdout);
	if (Nwrite != Nspace) {
		sprintf(messageStr, "Unable to write Psi to standard output");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}

	getCPU(1, TIME_START, NULL);
	ShortChar(&geometry, nspect = 0, mu, to_observer = TRUE, chi, S, I, Psi);
	getCPU(1, TIME_POLL, "Formal Solution up");
	Nwrite = fwrite(I, sizeof(double), Nspace, stdout);
	if (Nwrite != Nspace) {
		sprintf(messageStr, "Unable to write I to standard output");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}
	Nwrite = fwrite(Psi, sizeof(double), Nspace, stdout);
	if (Nwrite != Nspace) {
		sprintf(messageStr, "Unable to write Psi to standard output");
		Error(ERROR_LEVEL_2, argv[0], messageStr);
	}

	printTotalCPU();
}
/* ------- end ---------------------------- formal2d.c -------------- */
