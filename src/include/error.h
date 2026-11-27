/* ---------------------------------------- error.h -----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Dec  8 14:57:23 1999 --

	   --------------------------                      ----------RH-- */

#ifndef __ERROR_H__
#define __ERROR_H__

/* --- Defines error numbers, messages and associated actions -- ---- */

/* ------- begin -------------------------- error.h ----------------- */

enum errorlevel { MESSAGE,
				  WARNING,
				  ERROR_LEVEL_1,
				  ERROR_LEVEL_2 };

/* --- Exit status of the RH programs (rhf1d, rhsc2d, rhsc3d):

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

#define RH_EXIT_NOT_CONVERGED 1
#define RH_EXIT_FATAL         70

/* --- Associated function prototypes --               -------------- */

void ErrorMessage(enum errorlevel level, const char *routineName,
				  const char *format, ...)
	__attribute__((format(printf, 3, 4)));
_Noreturn void ErrorExit(enum errorlevel level, const char *routineName,
						 const char *format, ...)
	__attribute__((format(printf, 3, 4)));
void SetExitFlag(int flag);
int  ExitStatus(void);

/* --- Problems of the run so far (SetExitFlag), its part of RHState
	   (rh->error, Q-46) --                            -------------- */

typedef struct {
	int exitFlags;
} ErrorState;

/* --- Error(level, routineName, format, ...): print a message
	   (MESSAGE) or warning (WARNING) and return, or print an error
	   (ERROR_LEVEL_1, ERROR_LEVEL_2) and exit with status
	   RH_EXIT_FATAL. The message is formatted as by printf, and the
	   compiler checks the format against the arguments (Q-48; before,
	   the caller formatted it into the global messageStr). A macro,
	   so that the compiler sees which branch a constant level takes
	   and knows that the error levels do not return.
	   --                                              -------------- */

#define Error(level, routineName, ...)                                         \
	((level) >= ERROR_LEVEL_1 ? ErrorExit((level), (routineName), __VA_ARGS__) \
							  : ErrorMessage((level), (routineName), __VA_ARGS__))

#endif /* !__ERROR_H__ */

/* ------- end ---------------------------- error.h ----------------- */
