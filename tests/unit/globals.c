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
#include "state.h"

static Geometry geometry_data;
static RHState  rh_state = RH_STATE_INIT(RH_TOPOLOGY, &geometry_data);
RHState *const  rh       = &rh_state;
