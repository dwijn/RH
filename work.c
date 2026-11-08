/* ------- file: -------------------------- work.c ------------------

	   Work counters for comparing solvers (SOLVER.md section 4.1, S-1).

	   --------------------------                      ----------RH-- */

/* --- Counts the unit of work of the iteration: wavelength solves (the
	   formal solution at one wavelength, all rays), split by the phase
	   of the run (initial background-scattering iteration, main
	   iteration, final solution) and by kind (with approximate
	   operator, plain, PRD sub-iteration over the PRD wavelengths only).
	   Calls of Background() and Hydrostatic() and their wall time are
	   counted separately: they do no formal solution but can be a large
	   part of the run (SOLVE_NE = ITERATION).

	   The counters are only touched from the main thread, between
	   sweeps, so they need no locking.
	   --                                              -------------- */

#include <stdio.h>
#include <time.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"
#include "error.h"
#include "statistics.h"

/* --- Global variables --                             -------------- */

extern Spectrum    spectrum;
extern InputData   input;
extern char        messageStr[];

static WorkCounters work = {WORK_INITIAL};

static const char *phaseName[WORK_NPHASE] = {"initial", "iterate", "final"};

/* ------- begin -------------------------- wallClock.c ------------- */

double wallClock(void)
{
	struct timespec ts;

	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + 1.0E-9 * ts.tv_nsec;
}
/* ------- end ---------------------------- wallClock.c ------------- */

/* ------- begin -------------------------- workPhase.c ------------- */

void workPhase(enum WorkPhase phase) { work.phase = phase; }
/* ------- end ---------------------------- workPhase.c ------------- */

/* ------- begin -------------------------- workSweep.c ------------- */

void workSweep(bool_t eval_operator, bool_t redistribute, int Nsolved)
{
	enum WorkKind kind;

	kind = (redistribute) ? WORK_PRD
						  : ((eval_operator) ? WORK_OPERATOR : WORK_PLAIN);
	work.sweeps[work.phase][kind]++;
	work.lambda[work.phase][kind] += Nsolved;
}
/* ------- end ---------------------------- workSweep.c ------------- */

/* ------- begin -------------------------- workAdd.c --------------- */

void workAdd(enum WorkTimer timer, double seconds)
{
	work.calls[timer]++;
	work.seconds[timer] += seconds;
}
/* ------- end ---------------------------- workAdd.c --------------- */

/* ------- begin -------------------------- printWork.c ------------- */

void printWork(void)
{
	static const char *kindName[WORK_NKIND] = {"operator", "plain", "PRD"};
	int                phase, kind;

	/* --- One line per phase, a fixed format that
		   tests/analyze_run.py parses --                -------------- */

	sprintf(messageStr,
			"\n -- Work (wavelength solves; one formal solution = %d)\n",
			spectrum.Nspect);
	Error(MESSAGE, NULL, messageStr);
	for (phase = 0; phase < WORK_NPHASE; phase++) {
		sprintf(messageStr, "    %-8s", phaseName[phase]);
		Error(MESSAGE, NULL, messageStr);
		for (kind = 0; kind < WORK_NKIND; kind++) {
			sprintf(messageStr, "  %s %ld (%ld sweeps)", kindName[kind],
					work.lambda[phase][kind], work.sweeps[phase][kind]);
			Error(MESSAGE, NULL, messageStr);
		}
		Error(MESSAGE, NULL, "\n");
	}
	sprintf(messageStr,
			"    Background %ld calls %.3f s;  Hydrostatic %ld calls %.3f s;"
			"  Iterate %.3f s;  threads %d\n",
			work.calls[WORK_BACKGROUND], work.seconds[WORK_BACKGROUND],
			work.calls[WORK_HYDROSTATIC], work.seconds[WORK_HYDROSTATIC],
			work.seconds[WORK_ITERATION], input.Nthreads);
	Error(MESSAGE, NULL, messageStr);
}
/* ------- end ---------------------------- printWork.c ------------- */
