/* ------- file: -------------------------- globals.c ---------------

	   The global variables that the RH programs define in their main
	   source file and that librh refers to, so that unit tests can link
	   against librh. geometry.h comes from the geometry directory of
	   the test, RH_TOPOLOGY from its CMake definition.
	   --                                              -------------- */

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "spectrum.h"
#include "statistics.h"
#include "inputs.h"

enum Topology topology = RH_TOPOLOGY;

Atmosphere   atmos;
Geometry     geometry;
Spectrum     spectrum;
ProgramStats stats;
InputData    input;
CommandLine  commandline;
char         messageStr[MAX_MESSAGE_LENGTH];
