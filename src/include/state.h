/* ------- file: -------------------------- state.h -----------------

	   The state of an RH run in one struct (Q-45, phase 2 of the
	   library track, docs/LIBRARY.md 6.3).

	   --------------------------                      ----------RH-- */

#ifndef __STATE_H__
#define __STATE_H__

/* --- What the programs used to define as separate globals. Each
	   program (and tests/unit/globals.c) defines the one instance
	   with RH_STATE_INIT and the pointer rh to it; the library reaches
	   the state through rh, as in rh->atmos.T. In phase 3 a parameter
	   RHState *rh takes over the name, and the code that uses it stays
	   as it is. rh is a const pointer, so the compiler loads it once
	   per function and keeps it across calls.

	 - The type of the geometry differs per geometry (geometry.h in
	   src/geometry/<geometry>), so the library holds a pointer to it
	   and never looks inside.

	 - messageStr stays a global until errors return (phase 3).
	   --                                              -------------- */

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"
#include "statistics.h"

struct Geometry;

typedef struct RHState {
	enum Topology    topology;
	Atmosphere       atmos;
	struct Geometry *geometry;
	Spectrum         spectrum;
	ProgramStats     stats;
	InputData        input;
	CommandLine      commandline;
} RHState;

/* --- Static initializer of the instance: the topology of the program
	   and the address of its geometry. --             -------------- */

#define RH_STATE_INIT(topo, geom) {.topology = (topo), .geometry = (geom)}

extern RHState *const rh;

#endif /* !__STATE_H__ */

/* ------- end ---------------------------- state.h ----------------- */
