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

	   And what the modules kept in static variables for the length of
	   a run (Q-46): one part per module, of a type declared in the
	   module's header, as rh->solver.Nstate. Their initial values are
	   in the module's <MODULE>_STATE_INIT; the others start at zero.
	   Still static: the tables loaded once (phase 4), the scratch and
	   return buffers and the warn-once flags (phase 3), and the file
	   handles (phase 5).
	   --                                              -------------- */

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"
#include "statistics.h"
#include "error.h"
#include "background.h"
#include "solver.h"
#include "bgscatter.h"
#include "negtau.h"

struct Geometry;

typedef struct RHState {
	enum Topology    topology;
	Atmosphere       atmos;
	struct Geometry *geometry;
	Spectrum         spectrum;
	ProgramStats     stats;
	InputData        input;
	CommandLine      commandline;

	/* --- Kept by the modules (Q-46) --               -------------- */

	SolverState     solver;     /* nlte/solver.c */
	BgScatterState  bgscatter;  /* transfer/bgscatter.c */
	NegTauState     negtau;     /* transfer/negtau.c */
	OpacityState    opacity;    /* transfer/opacity.c */
	BackgroundState background; /* physics/background.c */
	KuruczState     kurucz;     /* physics/kurucz.c */
	HydrogenState   hydrogen;   /* physics/hydrogen.c */
	ChemEquilState  chemequil;  /* physics/chemequil.c */
	ErrorState      error;      /* runtime/error.c */
	GetCPUState     getcpu;     /* runtime/getcpu.c */
	WorkCounters    work;       /* runtime/work.c */
} RHState;

/* --- Static initializer of the instance: the topology of the program,
	   the address of its geometry (NULL in the tools), and the initial
	   values of the modules. --                       -------------- */

#define RH_STATE_INIT(topo, geom) {.topology  = (topo),               \
								   .geometry  = (geom),               \
								   .solver    = SOLVER_STATE_INIT,    \
								   .bgscatter = BGSCATTER_STATE_INIT, \
								   .negtau    = NEGTAU_STATE_INIT,    \
								   .opacity   = OPACITY_STATE_INIT}

extern RHState *const rh;

#endif /* !__STATE_H__ */

/* ------- end ---------------------------- state.h ----------------- */
