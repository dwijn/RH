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


/* --- Associated function prototypes --               -------------- */

void ErrorMessage(enum errorlevel level, const char *routineName,
		  const char *messageStr);
_Noreturn void ErrorExit(enum errorlevel level, const char *routineName,
			 const char *messageStr);

/* --- Print a message (MESSAGE) or warning (WARNING) and return, or
       print an error (ERROR_LEVEL_1, ERROR_LEVEL_2) and exit with the
       level as exit status. Inline, so that the compiler knows that the
       error levels do not return --                   -------------- */

static inline void Error(enum errorlevel level, const char *routineName,
			 const char *messageStr)
{
  if (level >= ERROR_LEVEL_1) ErrorExit(level, routineName, messageStr);
  ErrorMessage(level, routineName, messageStr);
}

#endif /* !__ERROR_H__ */

/* ------- end ---------------------------- error.h ----------------- */
