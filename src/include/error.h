/* ---------------------------------------- error.h -----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Dec  8 14:57:23 1999 --

	   --------------------------                      ----------RH-- */

#ifndef __ERROR_H__
#define __ERROR_H__

#include <stdatomic.h>

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
												   status RH_EXIT_FATAL (in the
												   code not yet converted to
												   errors that return, below)

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

void SetExitFlag(RHState *rh, int flag);
int  ExitStatus(RHState *rh);

/* --- Errors that return (Q-51, docs/LIBRARY.md 3.6). A function that
	   can fail returns an rh_status; one that computes a value returns
	   the status, and the value through an out-parameter. The codes are
	   those of the API (LIBRARY.md 4), rh_status in rh.h:

		 RH_E_INPUT     the configuration is invalid
		 RH_E_DATA      a data file is missing, unreadable or invalid
		 RH_E_NUMERIC   a numerical procedure failed
		 RH_E_TRIAL     the outer solver's trial point is invalid
		 RH_E_NOMEM     an allocation failed
		 RH_E_STOPPED   the caller asked the solver to stop
		 RH_E_STATE     a call in the wrong state (the API)
		 RH_E_INTERNAL  reserved for the API

	   Fail(rh, code, routine, format, ...)
		   records the error in rh->error and returns code, as in
		   return Fail(rh, RH_E_DATA, routineName, "...", ...);
		   The first error is kept, also when worker threads fail at
		   the same time; later ones only return their code.
	   RH_TRY(call)
		   passes a failure of call on: returns its status from the
		   function that contains it unless it is RH_OK.
	   ExitOnError(rh, status)
		   for a caller not yet converted: when status is not RH_OK,
		   prints the recorded error as Error() does and exits with
		   RH_EXIT_FATAL, so the run behaves as before.
	   RH_ASSERT(rh, condition, routine, format, ...)
		   an invariant: false only through a bug in RH. Prints the
		   message, the condition and where it is, and aborts, also in
		   Release builds. Not for anything that input or data can
		   cause.
	   Failed(rh)      whether an error is recorded (workers poll it to
					   stop early)
	   ClearError(rh)  forgets the recorded error (not while threads run)
	   --                                              -------------- */

rh_status Fail(RHState *rh, rh_status code, const char *routineName,
			   const char *format, ...) __attribute__((format(printf, 4, 5)));
void      ExitOnError(RHState *rh, rh_status status);
bool_t    Failed(RHState *rh);
void      ClearError(RHState *rh);

_Noreturn void AssertFailed(RHState *rh, const char *condition, const char *file,
							int line, const char *routineName, const char *format, ...)
	__attribute__((format(printf, 6, 7)));

#define RH_TRY(call)                       \
	do {                                   \
		rh_status rh_try_status_ = (call); \
		if (rh_try_status_ != RH_OK)       \
			return rh_try_status_;         \
	} while (0)

#define RH_ASSERT(rh, condition, routineName, ...)             \
	do {                                                       \
		if (!(condition))                                      \
			AssertFailed((rh), #condition, __FILE__, __LINE__, \
						 (routineName), __VA_ARGS__);          \
	} while (0)

/* --- Problems of the run so far (SetExitFlag), and the recorded
	   error (Fail): its part of RHState (rh->error, Q-46). failState is
	   0 with no error, 1 while Fail records one, 2 once it is recorded
	   --                                              -------------- */

typedef struct {
	int        exitFlags;
	atomic_int failState;
	rh_status  code;
	char       routine[64];
	char      *message; /* allocated by Fail, freed by ClearError */
} ErrorState;

#endif /* !__ERROR_H__ */

/* ------- end ---------------------------- error.h ----------------- */
