/* ------- file: -------------------------- Error.c -----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Dec 23 14:47:16 2008 --

	   --------------------------                      ----------RH-- */

/* --- Messages and warnings (Log, through the macros Message and
	   Warning of error.h), errors that end the run (Error), and errors
	   that return (Fail and its companions, Q-51) --  -------------- */

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

/* --- The end of a run on an error: the message to the log file and
	   to stderr, and exit with status RH_EXIT_FATAL -- -------------- */

static _Noreturn void terminate(RHState *rh, const char *routineName,
								const char *text)
{
	fprintf(rh->commandline.logfile,
			"\a\n\n-TERMINATING_ERROR in routine %s\n %s \n Exiting.....\n",
			routineName, text);
	if (rh->commandline.logfile != stderr)
		fprintf(stderr,
				"\a\n\n-TERMINATING_ERROR in routine %s\n %s \n Exiting.....\n",
				routineName, text);
	exit(RH_EXIT_FATAL);
}

/* --- A fatal error where the code does not return errors yet.
	   ERROR_LEVEL_1 and _2 end the run alike (Q-16) -- ------------- */

void Error(RHState *rh, enum errorlevel level, const char *routineName,
		   const char *format, ...)
{
	char    buf[MAX_MESSAGE_LENGTH], *text;
	va_list ap;

	va_start(ap, format);
	text = formatMessage(buf, sizeof(buf), format, ap);
	va_end(ap);

	(void) level;
	terminate(rh, routineName, text);
}
/* ------- end ---------------------------- Error.c ----------------- */

/* ------- begin -------------------------- Fail.c ------------------ */

/* --- Records an error and returns its code (error.h). The first
	   caller to move failState from 0 to 1 writes the record and then
	   sets it to 2; any other caller, also in another thread at the
	   same moment, only returns its code --          -------------- */

rh_status Fail(RHState *rh, rh_status code, const char *routineName,
			   const char *format, ...)
{
	char    buf[MAX_MESSAGE_LENGTH], *text;
	int     none = 0;
	va_list ap;

	if (!atomic_compare_exchange_strong(&rh->error.failState, &none, 1))
		return code;

	va_start(ap, format);
	text = formatMessage(buf, sizeof(buf), format, ap);
	va_end(ap);

	rh->error.code = code;
	snprintf(rh->error.routine, sizeof(rh->error.routine), "%s",
			 routineName ? routineName : "");
	rh->error.message = strdup(text);
	freeMessage(text, buf);

	atomic_store(&rh->error.failState, 2);
	return code;
}

bool_t Failed(RHState *rh)
{
	return atomic_load(&rh->error.failState) != 0;
}

void ClearError(RHState *rh)
{
	free(rh->error.message);
	rh->error.message    = NULL;
	rh->error.code       = RH_OK;
	rh->error.routine[0] = '\0';
	atomic_store(&rh->error.failState, 0);
}
/* ------- end ---------------------------- Fail.c ------------------ */

/* ------- begin -------------------------- ExitOnError.c ----------- */

/* --- For the callers not yet converted: a status other than RH_OK
	   ends the run with the recorded error, printed as Error() prints
	   it --                                           -------------- */

void ExitOnError(RHState *rh, rh_status status)
{
	if (status == RH_OK) return;

	if (atomic_load(&rh->error.failState) == 2)
		terminate(rh, rh->error.routine,
				  rh->error.message ? rh->error.message : "(no memory for the message)");
	else {
		char text[MAX_MESSAGE_LENGTH];

		snprintf(text, sizeof(text),
				 "Status %d without a recorded error (a bug in RH)", (int) status);
		terminate(rh, "ExitOnError", text);
	}
}
/* ------- end ---------------------------- ExitOnError.c ----------- */

/* ------- begin -------------------------- AssertFailed.c ---------- */

/* --- A broken invariant (RH_ASSERT): a bug in RH. The message, the
	   condition and where it is go to the log file and to stderr, and
	   the process aborts --                           -------------- */

void AssertFailed(RHState *rh, const char *condition, const char *file, int line,
				  const char *routineName, const char *format, ...)
{
	char        buf[MAX_MESSAGE_LENGTH], *text;
	const char *src = strstr(file, "src/");
	va_list     ap;

	va_start(ap, format);
	text = formatMessage(buf, sizeof(buf), format, ap);
	va_end(ap);

	if (src) file = src;
	if (rh && rh->commandline.logfile && rh->commandline.logfile != stderr)
		fprintf(rh->commandline.logfile,
				"\a\n\n-INTERNAL_ERROR in routine %s\n %s \n"
				" A bug in RH: %s is false (%s:%d)\n Aborting.....\n",
				routineName, text, condition, file, line);
	fprintf(stderr,
			"\a\n\n-INTERNAL_ERROR in routine %s\n %s \n"
			" A bug in RH: %s is false (%s:%d)\n Aborting.....\n",
			routineName, text, condition, file, line);
	fflush(NULL);
	abort();
}
/* ------- end ---------------------------- AssertFailed.c ---------- */

/* ------- begin -------------------------- ExitStatus.c ------------ */

/* --- Problems of a completed run, see error.h --     -------------- */

void SetExitFlag(RHState *rh, int flag)
{
	/* --- Only bits 1-63: higher statuses mean failure -- ---------- */

	rh->error.exitFlags |= (flag & 63);
}

int ExitStatus(RHState *rh)
{
	return rh->error.exitFlags;
}
/* ------- end ---------------------------- ExitStatus.c ------------ */
