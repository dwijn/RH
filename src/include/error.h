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

/* --- Messages, warnings and fatal errors (Q-48, Q-49):

	   Message(rh, format, ...)                a message, unless the
												   command line asks for quiet
	   Warning(rh, routine, format, ...)       a warning
	   Error(rh, level, routine, format, ...)  a fatal error, level
												   ERROR_LEVEL_1 or _2: prints
												   the message and exits with
												   status RH_EXIT_FATAL (until
												   errors return, Q-51)

	   The format is a printf format, which the compiler checks against
	   the arguments. Message and Warning are macros over Log(), the one
	   function that writes messages and warnings (from phase 6 it hands
	   them to the caller's callback). Macros and not inline functions:
	   C cannot pass `...` on, and gcc and clang never inline a function
	   with a variable argument list. --               -------------- */

typedef struct RHState RHState;

void           Log(RHState *rh, enum errorlevel level, const char *routineName,
				   const char *format, ...) __attribute__((format(printf, 4, 5)));
_Noreturn void Error(RHState *rh, enum errorlevel level,
					 const char *routineName, const char *format, ...)
	__attribute__((format(printf, 4, 5)));

#define Message(rh, ...) Log((rh), MESSAGE, NULL, __VA_ARGS__)
#define Warning(rh, routineName, ...) \
	Log((rh), WARNING, (routineName), __VA_ARGS__)

/* --- The one instance of the state, until every function that uses
	   it is passed it (Q-49) --                       -------------- */

extern RHState *const rh;

void SetExitFlag(RHState *rh, int flag);
int  ExitStatus(RHState *rh);

/* --- Problems of the run so far (SetExitFlag), its part of RHState
	   (rh->error, Q-46) --                            -------------- */

typedef struct {
	int exitFlags;
} ErrorState;

#endif /* !__ERROR_H__ */

/* ------- end ---------------------------- error.h ----------------- */
