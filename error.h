/* ---------------------------------------- error.h -----------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Wed Dec  8 14:57:23 1999 --

       --------------------------                      ----------RH-- */

#ifndef __ERROR_H__
#define __ERROR_H__

/* --- Defines error numbers, messages and associated actions -- ---- */


/* ------- begin -------------------------- error.h ----------------- */

enum errorlevel {MESSAGE, WARNING, ERROR_LEVEL_1, ERROR_LEVEL_2};

/* --- Exit status of the RH programs (rhf1d, rhsc2d, rhsc3d,
       rhsphere):

         0        the run completed and converged
         1 - 63   the run completed, with problems: a bitmask (several
                  can occur together; the output files are written)
                    1  populations not converged (N_MAX_ITER reached,
                       or stopped with STOP_RH)
                  (PRD, chemical-equilibrium and hydrostatic convergence
                  are not checked yet)
         70       fatal error (EX_SOFTWARE in sysexits.h): the run
                  stopped, see the TERMINATING_ERROR message
         128 + n  killed by signal n, e.g. 136 for a trapped floating-
                  point exception (SIGFPE)

       Bits are set with SetExitFlag() during the run; main() returns
       ExitStatus(). --                                -------------- */

#define RH_EXIT_NOT_CONVERGED  1
#define RH_EXIT_FATAL         70


/* --- Associated function prototypes --               -------------- */

void ErrorMessage(enum errorlevel level, const char *routineName,
		  const char *messageStr);
_Noreturn void ErrorExit(enum errorlevel level, const char *routineName,
			 const char *messageStr);
void SetExitFlag(int flag);
int  ExitStatus(void);

/* --- Print a message (MESSAGE) or warning (WARNING) and return, or
       print an error (ERROR_LEVEL_1, ERROR_LEVEL_2) and exit with status
       RH_EXIT_FATAL. Inline, so that the compiler knows that the
       error levels do not return --                   -------------- */

static inline void Error(enum errorlevel level, const char *routineName,
			 const char *messageStr)
{
  if (level >= ERROR_LEVEL_1) ErrorExit(level, routineName, messageStr);
  ErrorMessage(level, routineName, messageStr);
}

#endif /* !__ERROR_H__ */

/* ------- end ---------------------------- error.h ----------------- */
