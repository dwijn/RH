/* ------- file: -------------------------- gettime.c ---------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed May 19 15:56:36 2021 --

	   --------------------------                      ----------RH-- */

/* --- Prints elapsed time in seconds:    action == TIME_POLL.
	   Clock is started if called with:   action == TIME_START.
	   Borrowed heavily from MULTI's ctime().

	   --                                              -------------- */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>

#include "rh.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "state.h"

#define TAB_POSITION 24
#define DOTS         "...................................."
#define SPACES       "                                    "
#define TIME_DOT_OUT "time.out"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- getTime.c --------------- */

void getCPU(RHState *rh, int level, enum CPUaction action, char *label)
{
	register int n;

	int Nblanck, Nspace, Ndot;

	static const double scale = 1.0 / CLOCKS_PER_SEC;
	clock_t             CPUtime;

	if (level == 0 && action == TIME_START && !rh->commandline.showkeywords) {
		rh->stats.fp_CPU = fopen(TIME_DOT_OUT, "w");
		setvbuf(rh->stats.fp_CPU, NULL, _IOLBF, BUFSIZ);
	}

	if (!rh->stats.printCPU || !rh->stats.fp_CPU)
		return;
	else {
		CPUtime = clock();

		/* --- The timing is optional output: switched off when the
			   processor time is not available -- ---------------- */

		if (CPUtime == (clock_t) -1) {
			Warning(rh, "getCPU", "Cannot poll the processor time: timing (PRINT_CPU) switched off");
			rh->stats.printCPU = FALSE;
			return;
		}

		switch (action) {
		case TIME_START:
			rh->getcpu.CPU[level] = CPUtime;
			return;
		case TIME_POLL:
			CPUtime = CPUtime - rh->getcpu.CPU[level];
			break;
		case TIME_ADD:
			break;
		}
	}
	/* --- Format output --                              -------------- */

	Nspace  = (level < 4) ? level + 1 : 4;
	Ndot    = TAB_POSITION - Nspace - strlen(label);
	Nblanck = (level > 2) ? 0 : 3 - level;

	fprintf(rh->stats.fp_CPU, "%.*s%.*s%s: %7.3f [sec]\n", Nspace, SPACES,
			Ndot, DOTS, label, CPUtime * scale);
	for (n = 0; n < Nblanck; n++) fprintf(rh->stats.fp_CPU, "\n");
}

/* ------- end ---------------------------- getTime.c --------------- */

/* ------- begin -------------------------- printTotalCPU ----------- */

void printTotalCPU(RHState *rh)
{
	double     user, system, total;
	struct tms tms;

	times(&tms);
	user   = ((double) tms.tms_utime) / sysconf(_SC_CLK_TCK);
	system = ((double) tms.tms_stime) / sysconf(_SC_CLK_TCK);
	total  = user + system;

	/* --- times() counts in clock ticks (typically 10 ms): a short run
		   can use none at all, and 0/0 would trap with FP exceptions
		   enabled (B-33) --                              -------------- */

	fprintf(rh->stats.fp_CPU, "\nTotal_time____: %8.3f [sec] ", total);
	if (total > 0.0)
		fprintf(rh->stats.fp_CPU, " (user: %4.1f%%, system: %4.1f%%)\n\n",
				100 * user / total, 100 * system / total);
	else
		fprintf(rh->stats.fp_CPU, " (less than one clock tick)\n\n");
}
/* ------- end ---------------------------- printTotalCPU ----------- */
