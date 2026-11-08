/* ------- file: -------------------------- statistics.h ------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Mon Jan 24 15:57:13 2000 --

	   --------------------------                      ----------RH-- */

#ifndef __STATISTICS_H__
#define __STATISTICS_H__

/* --- Keep track of program statistics like CPU time and memory usage
	   (not yet implemented).
	   --                                              -------------- */

#define N_TIME_LEVELS 5

enum CPUaction { TIME_START,
				 TIME_POLL,
				 TIME_ADD };

typedef struct {
	bool_t printCPU;
	long   memory;
	FILE  *fp_CPU;
} ProgramStats;

/* --- Associated function prototypes --               -------------- */

void getCPU(int level, enum CPUaction action, char *label);
void printTotalCPU(void);

/* --- Work counters (work.c): wavelength solves per phase of the run
	   and kind of sweep, calls and wall time of Background() and
	   Hydrostatic(). Needs rh.h for bool_t --         -------------- */

enum WorkPhase { WORK_INITIAL,
				 WORK_ITERATE,
				 WORK_FINAL,
				 WORK_NPHASE };
enum WorkKind { WORK_OPERATOR,
				WORK_PLAIN,
				WORK_PRD,
				WORK_NKIND };
enum WorkTimer { WORK_BACKGROUND,
				 WORK_HYDROSTATIC,
				 WORK_ITERATION,
				 WORK_NTIMER };

typedef struct {
	enum WorkPhase phase;
	long           sweeps[WORK_NPHASE][WORK_NKIND],
		lambda[WORK_NPHASE][WORK_NKIND], calls[WORK_NTIMER];
	double seconds[WORK_NTIMER];
} WorkCounters;

double wallClock(void);
void   workPhase(enum WorkPhase phase);
void   workSweep(bool_t eval_operator, bool_t redistribute, int Nsolved);
void   workAdd(enum WorkTimer timer, double seconds);
void   printWork(void);

#endif /* !__STATISTICS_H__ */

/* ------- end ---------------------------- statistics.h ------------ */
