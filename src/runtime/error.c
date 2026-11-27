/* ------- file: -------------------------- Error.c -----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Dec 23 14:47:16 2008 --

	   --------------------------                      ----------RH-- */

/* --- Messages and warnings (Log, through the macros Message and
	   Warning of error.h), and errors that end the run (Error) --
													   -------------- */

#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- formatMessage.c --------- */

/* --- The message formatted into buf (of size bufsize), or into
	   allocated memory if it is longer: no length limit (B-14, B-41).
	   Free the result with freeMessage. A buffer per call, not a
	   global, so that threads can report at the same time (T-1).
	   --                                              -------------- */

static char *formatMessage(char *buf, size_t bufsize, const char *format,
						   va_list ap)
{
	va_list ap2;
	char   *text;
	int     length;

	va_copy(ap2, ap);
	length = vsnprintf(buf, bufsize, format, ap2);
	va_end(ap2);
	if (length < 0) return strcpy(buf, " (message could not be formatted)");
	if ((size_t) length < bufsize) return buf;

	if ((text = malloc((size_t) length + 1)) == NULL) return buf;
	vsnprintf(text, (size_t) length + 1, format, ap);
	return text;
}

static void freeMessage(char *text, char *buf)
{
	if (text != buf) free(text);
}
/* ------- end ---------------------------- formatMessage.c --------- */

/* ------- begin -------------------------- Log.c ------------------- */

/* --- A message (MESSAGE) or a warning (WARNING); an error level is
	   reported as an error, without ending the run (fatal errors go
	   through Error). Each message is written with one call, so that
	   the messages of threads do not mix --           -------------- */

void Log(RHState *rh, enum errorlevel level, const char *routineName,
		 const char *format, ...)
{
	char    buf[MAX_MESSAGE_LENGTH], *text;
	va_list ap;

	va_start(ap, format);
	if (level == MESSAGE) {
		if (!rh->commandline.quiet)
			vfprintf(rh->commandline.logfile, format, ap);
	} else {
		text = formatMessage(buf, sizeof(buf), format, ap);
		fprintf(rh->commandline.logfile, "\n-%s in routine %s\n %s\n",
				(level == WARNING) ? "WARNING" : "ERROR", routineName, text);
		freeMessage(text, buf);
	}
	va_end(ap);
}
/* ------- end ---------------------------- Log.c ------------------- */

/* ------- begin -------------------------- Error.c ----------------- */

/* --- A fatal error: the message to the log file and to stderr, and
	   exit with status RH_EXIT_FATAL. ERROR_LEVEL_1 and _2 end the run
	   alike (Q-16) --                                 -------------- */

void Error(RHState *rh, enum errorlevel level, const char *routineName,
		   const char *format, ...)
{
	char    buf[MAX_MESSAGE_LENGTH], *text;
	va_list ap;

	va_start(ap, format);
	text = formatMessage(buf, sizeof(buf), format, ap);
	va_end(ap);

	fprintf(rh->commandline.logfile,
			"\a\n\n-TERMINATING_ERROR in routine %s\n %s \n Exiting.....\n",
			routineName, text);
	if (rh->commandline.logfile != stderr)
		fprintf(stderr,
				"\a\n\n-TERMINATING_ERROR in routine %s\n %s \n Exiting.....\n",
				routineName, text);

	(void) level;
	exit(RH_EXIT_FATAL);
}
/* ------- end ---------------------------- Error.c ----------------- */

/* ------- begin -------------------------- ExitStatus.c ------------ */

/* --- Problems of a completed run, see error.h --     -------------- */

void SetExitFlag(int flag)
{
	/* --- Only bits 1-63: higher statuses mean failure -- ---------- */

	rh->error.exitFlags |= (flag & 63);
}

int ExitStatus(void)
{
	return rh->error.exitFlags;
}
/* ------- end ---------------------------- ExitStatus.c ------------ */
