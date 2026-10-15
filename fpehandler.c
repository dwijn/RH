/* ------- file: -------------------------- fpehandler.c ------------

       Version:       rh2.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Thu May 20 11:29:23 2021 --

       --------------------------                      ----------RH-- */

/* --- Trap floating point exceptions on various machines:

       --                                              -------------- */

#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "error.h"

extern char messageStr[];


#if defined(Darwin)

void SetFPEtraps(void)
{
  /* --- Explicitly do not set traps --                -------------- */

  Error(MESSAGE, "SetFPEtraps", 
	"\n Darwin: FPE traps have not been set explicitly\n");
}
/* --- end SETNOTRAPS --                               -------------- */

#elif defined(Linux)


#define _GNU_SOURCE 1
#include <execinfo.h>
#include <fenv.h>
#include <signal.h>
#include <string.h>
#include <unistd.h>

int feenableexcept(int excepts);

#define FPE_MAX_FRAMES 32

/* --- Write a string to stderr from a signal handler (write() is
       async-signal-safe, stdio is not) --              -------------- */

static void fpeWrite(const char *str)
{
  ssize_t n = write(STDERR_FILENO, str, strlen(str));
  (void) n;
}

static void fpeWriteHex(unsigned long value)
{
  char buf[2 + 2*sizeof(unsigned long) + 1], *p = buf + sizeof(buf) - 1;

  *p = '\0';
  do {
    *--p = "0123456789abcdef"[value & 0xf];
    value >>= 4;
  } while (value);
  *--p = 'x';  *--p = '0';
  fpeWrite(p);
}

/* ------- begin -------------------------- fpeHandler.c ------------ */

static void fpeHandler(int sig, siginfo_t *info, void *context)
{
  /* --- Say what happened and where before the run dies, instead of
         the bare "Floating point exception" of the default action.
         The signal is then raised again with the default action, so
         the exit status is unchanged --               -------------- */

  void *frames[FPE_MAX_FRAMES];
  int   Nframes;
  const char *what;

  (void) context;
  switch (info->si_code) {
  case FPE_FLTDIV:  what = "floating-point division by zero";  break;
  case FPE_FLTINV:  what = "invalid floating-point operation (NaN)";  break;
  case FPE_FLTOVF:  what = "floating-point overflow";  break;
  case FPE_INTDIV:  what = "integer division by zero";  break;
  default:          what = "arithmetic exception";  break;
  }
  fpeWrite("\n-FATAL ERROR: ");
  fpeWrite(what);
  fpeWrite(" at address ");
  fpeWriteHex((unsigned long) info->si_addr);
  fpeWrite("\n Backtrace (resolve with: addr2line -f -e <program> <address>):\n");

  Nframes = backtrace(frames, FPE_MAX_FRAMES);
  backtrace_symbols_fd(frames, Nframes, STDERR_FILENO);

  signal(sig, SIG_DFL);
  raise(sig);
}
/* ------- end ---------------------------- fpeHandler.c ------------ */

void SetFPEtraps(void)
{
  struct sigaction action;
  void  *frame[1];

  /* --- Enable some exceptions.
         At startup all exceptions are masked. --      -------------- */

  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

  /* --- Report trapped exceptions. backtrace() is called once here
         so that its library is loaded before it is needed inside the
         signal handler --                             -------------- */

  backtrace(frame, 1);
  memset(&action, 0, sizeof(action));
  action.sa_sigaction = fpeHandler;
  action.sa_flags     = SA_SIGINFO;
  sigemptyset(&action.sa_mask);
  sigaction(SIGFPE, &action, NULL);
}

#else

/* --- unknown --                                      -------------- */

void SetFPEtraps(void)
{
  Error(MESSAGE, "SetFPEtraps",
	"\nUnsupported CPU and/or OS: cannot set FPE traps explicitly\n");
}

#endif

/* ------- end ------------------------------------------------------ */
