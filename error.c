/* ------- file: -------------------------- Error.c -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Dec 23 14:47:16 2008 --

       --------------------------                      ----------RH-- */

/* --- Messages and warnings (ErrorMessage), and errors that end the
       run (ErrorExit). Code calls them through Error() in error.h --
                                                       -------------- */

#include <stdlib.h>
#include <errno.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern CommandLine commandline;


/* ------- begin -------------------------- ErrorMessage.c ---------- */

void ErrorMessage(enum errorlevel level, const char *routineName,
		  const char *messageStr)
{
  if (level == MESSAGE) {
    if (!commandline.quiet)
      fprintf(commandline.logfile, "%s", (messageStr) ? messageStr : "");
  } else {
    fprintf(commandline.logfile, "\n-WARNING in routine %s\n %s\n",
	    routineName, (messageStr) ? messageStr : " (Undocumented)\n");
  }
}
/* ------- end ---------------------------- ErrorMessage.c ---------- */

/* ------- begin -------------------------- ErrorExit.c ------------- */

void ErrorExit(enum errorlevel level, const char *routineName,
	       const char *messageStr)
{
  char errorStr[MAX_MESSAGE_LENGTH];

  sprintf(errorStr, "\a\n\n-TERMINATING_ERROR in routine %s\n %s \n %s\n",
	  routineName,(messageStr) ? messageStr : " (Undocumented)\n",
	  "Exiting.....");

  fprintf(commandline.logfile, "%s", errorStr);
  if (commandline.logfile != stderr) fprintf(stderr, "%s", errorStr);

  (void) level;
  exit(RH_EXIT_FATAL);
}
/* ------- end ---------------------------- ErrorExit.c ------------- */

/* ------- begin -------------------------- ExitStatus.c ------------ */

/* --- Problems of a completed run, see error.h --     -------------- */

static int exitFlags = 0;

void SetExitFlag(int flag)
{
  exitFlags |= flag;
}

int ExitStatus(void)
{
  return exitFlags;
}
/* ------- end ---------------------------- ExitStatus.c ------------ */
