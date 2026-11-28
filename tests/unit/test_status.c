/* ------- file: -------------------------- test_status.c -----------

	   Errors that return (error.c, Q-51): Fail, RH_TRY, RH_TRY_GOTO,
	   ExitOnError, RH_ASSERT, Failed and ClearError.

	   -- Fail records the code, the routine and the whole message (also
		  one far longer than MAX_MESSAGE_LENGTH) and returns the code;
		  a second Fail returns its own code and leaves the first record;
	   -- RH_TRY returns a failure from the function that contains it
		  and goes on after RH_OK;
	   -- RH_TRY_GOTO goes on after RH_OK; on a failure it stores the
		  status and jumps to the label;
	   -- ClearError forgets the record;
	   -- ExitOnError returns on RH_OK; on a failure it prints the
		  recorded error exactly as Error() prints the same routine and
		  message, to the log file and to stderr, and exits with status
		  RH_EXIT_FATAL (in a child process), also when the failure was
		  not recorded;
	   -- RH_ASSERT does nothing when its condition holds; otherwise it
		  prints the message, the condition and where it is, and aborts
		  (SIGABRT, in a child process);
	   -- threads that fail at the same moment leave one whole record,
		  the text of one of them with its own code.
	   --                                              -------------- */

#include <signal.h>
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/wait.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"
#include "state.h"
#include "check.h"

/* --- Everything that can be read from fd until end of file -- */

static char *readAll(int fd)
{
	size_t  size = 0, cap = 4096;
	char   *s = malloc(cap);
	ssize_t n;

	while ((n = read(fd, s + size, cap - size - 1)) > 0) {
		size += (size_t) n;
		if (cap - size < 2) s = realloc(s, cap *= 2);
	}
	s[size] = '\0';
	return s;
}

/* --- What a child writes to the log file and to stderr, and how it
	   ends: its exit status, or 128 + the signal that killed it -- */

static int inChild(void (*body)(void), char **logOut, char **errOut)
{
	int   plog[2], perr[2], status;
	pid_t pid;

	if (pipe(plog) || pipe(perr)) {
		*logOut = calloc(1, 1);
		*errOut = calloc(1, 1);
		return -1;
	}
	fflush(NULL);
	if ((pid = fork()) == 0) {
		close(plog[0]);
		close(perr[0]);
		dup2(perr[1], STDERR_FILENO);
		rh->commandline.logfile = fdopen(plog[1], "w");
		body();
		fflush(NULL);
		_exit(0);
	}
	close(plog[1]);
	close(perr[1]);
	*logOut = readAll(plog[0]);
	*errOut = readAll(perr[0]);
	close(plog[0]);
	close(perr[0]);
	waitpid(pid, &status, 0);
	if (WIFEXITED(status)) return WEXITSTATUS(status);
	return WIFSIGNALED(status) ? 128 + WTERMSIG(status) : -1;
}

/* --- The bodies of the children --                  -------------- */

static void errorBody(void)
{
	Error(rh, ERROR_LEVEL_2, "testRoutine", "value %d out of range [%g, %g]", 12, 0.5, 2.0);
}

static void exitOnErrorBody(void)
{
	ExitOnError(rh, Fail(rh, RH_E_DATA, "testRoutine", "value %d out of range [%g, %g]", 12, 0.5, 2.0));
}

static void unrecordedBody(void)
{
	ClearError(rh);
	ExitOnError(rh, RH_E_NUMERIC);
}

static void assertBody(void)
{
	int n = 3;

	RH_ASSERT(rh, n > 0, "testRoutine", "never printed");
	RH_ASSERT(rh, n < 0, "testRoutine", "n = %d", n);
}

/* --- RH_TRY: inner fails or not; outer passes a failure on -- */

static int reached;

static rh_status inner(bool_t fail)
{
	return fail ? Fail(rh, RH_E_NUMERIC, "inner", "inner failed") : RH_OK;
}

static rh_status outer(bool_t fail)
{
	reached = 0;
	RH_TRY(inner(fail));
	reached = 1;
	return RH_OK;
}

/* --- RH_TRY_GOTO: reached stays 0 on a failure, cleaned is set on
	   both paths -- ------------------------------------------------ */

static int cleaned;

static rh_status outerCleanup(bool_t fail)
{
	rh_status result = RH_OK;

	reached = cleaned = 0;
	RH_TRY_GOTO(result, cleanup, inner(fail));
	reached = 1;

cleanup:
	cleaned = 1;
	return result;
}

/* --- Threads failing at once --                      -------------- */

#define NTHREADS 8

/* --- The threads wait until all have started, then fail together
	   (a spin on a counter: macOS has no pthread barriers) -- ------ */

static atomic_int started;

static void *failAtOnce(void *arg)
{
	int t = *(int *) arg;

	atomic_fetch_add(&started, 1);
	while (atomic_load(&started) < NTHREADS);
	Fail(rh, (rh_status) (RH_E_INPUT + t % 4), "failAtOnce", "thread %d code %d",
		 t, RH_E_INPUT + t % 4);
	return NULL;
}

int main(void)
{
	char     *logOut, *errOut, *logRef, *errRef, *longText;
	int       status;
	rh_status s;

	/* --- Fail records the first error and returns its code -- -- */

	CHECK(!Failed(rh), "no error before any Fail");
	s = Fail(rh, RH_E_DATA, "testRoutine", "file %s: %d values", "x.atom", 3);
	CHECK(s == RH_E_DATA, "Fail returns its code, %d", s);
	CHECK(Failed(rh), "Failed after Fail");
	CHECK(rh->error.code == RH_E_DATA, "recorded code %d", rh->error.code);
	CHECK(strcmp(rh->error.routine, "testRoutine") == 0, "routine %s", rh->error.routine);
	CHECK(strcmp(rh->error.message, "file x.atom: 3 values") == 0, "message %s",
		  rh->error.message);

	s = Fail(rh, RH_E_NUMERIC, "otherRoutine", "second");
	CHECK(s == RH_E_NUMERIC, "a second Fail returns its own code, %d", s);
	CHECK(rh->error.code == RH_E_DATA && strcmp(rh->error.message, "file x.atom: 3 values") == 0,
		  "the first record stays: %d %s", rh->error.code, rh->error.message);

	ClearError(rh);
	CHECK(!Failed(rh) && rh->error.message == NULL && rh->error.code == RH_OK,
		  "ClearError forgets the record");

	/* --- A message far longer than MAX_MESSAGE_LENGTH, whole -- -- */

	longText = malloc(5 * MAX_MESSAGE_LENGTH + 1);
	memset(longText, 'x', 5 * MAX_MESSAGE_LENGTH);
	longText[5 * MAX_MESSAGE_LENGTH] = '\0';
	Fail(rh, RH_E_INPUT, "testRoutine", "<%s>", longText);
	CHECK(strlen(rh->error.message) == 5 * MAX_MESSAGE_LENGTH + 2, "long message length %zu",
		  strlen(rh->error.message));
	ClearError(rh);

	/* --- A routine name longer than 64 characters, whole: the
		   programs pass argv[0], a path, as the routine of checkNread
		   (Q-56) --                                     -------------- */

	longText[200] = '\0';
	Fail(rh, RH_E_DATA, longText, "short");
	CHECK(strcmp(rh->error.routine, longText) == 0, "long routine name length %zu",
		  strlen(rh->error.routine));
	ClearError(rh);
	CHECK(rh->error.routine == NULL, "ClearError forgets the routine name");
	free(longText);

	/* --- RH_TRY --                                      -------------- */

	CHECK(outer(FALSE) == RH_OK && reached == 1, "RH_TRY goes on after RH_OK");
	CHECK(!Failed(rh), "nothing recorded on success");
	CHECK(outer(TRUE) == RH_E_NUMERIC && reached == 0, "RH_TRY returns the failure");
	CHECK(Failed(rh) && strcmp(rh->error.routine, "inner") == 0, "the inner error is recorded");
	ClearError(rh);

	/* --- RH_TRY_GOTO --                                 -------------- */

	CHECK(outerCleanup(FALSE) == RH_OK && reached == 1 && cleaned == 1,
		  "RH_TRY_GOTO goes on after RH_OK");
	CHECK(outerCleanup(TRUE) == RH_E_NUMERIC && reached == 0 && cleaned == 1,
		  "RH_TRY_GOTO jumps to the label with the failure");
	ClearError(rh);

	/* --- ExitOnError --                                 -------------- */

	ExitOnError(rh, RH_OK);
	CHECK(1, "ExitOnError returns on RH_OK");

	status = inChild(errorBody, &logRef, &errRef);
	CHECK(status == RH_EXIT_FATAL, "Error: exit status %d", status);
	status = inChild(exitOnErrorBody, &logOut, &errOut);
	CHECK(status == RH_EXIT_FATAL, "ExitOnError: exit status %d", status);
	CHECK(strcmp(logOut, logRef) == 0, "ExitOnError logs what Error logs:\n%s\nvs\n%s", logOut, logRef);
	CHECK(strcmp(errOut, errRef) == 0, "ExitOnError writes to stderr what Error writes:\n%s\nvs\n%s",
		  errOut, errRef);
	CHECK(strstr(errOut, "-TERMINATING_ERROR in routine testRoutine\n value 12 out of range [0.5, 2]") != NULL,
		  "stderr %s", errOut);
	free(logOut);
	free(errOut);
	free(logRef);
	free(errRef);

	status = inChild(unrecordedBody, &logOut, &errOut);
	CHECK(status == RH_EXIT_FATAL, "unrecorded failure: exit status %d", status);
	CHECK(strstr(errOut, "Status 3 without a recorded error") != NULL, "stderr %s", errOut);
	free(logOut);
	free(errOut);

	/* --- RH_ASSERT --                                   -------------- */

	status = inChild(assertBody, &logOut, &errOut);
	CHECK(status == 128 + SIGABRT, "RH_ASSERT: status %d, want %d", status, 128 + SIGABRT);
	CHECK(strstr(errOut, "-INTERNAL_ERROR in routine testRoutine\n n = 3 \n") != NULL,
		  "stderr %s", errOut);
	CHECK(strstr(errOut, "A bug in RH: n < 0 is false (") != NULL, "condition in %s", errOut);
	CHECK(strstr(errOut, "never printed") == NULL, "a holding assertion prints nothing");
	CHECK(strcmp(logOut, errOut) == 0, "the same to the log file: %s", logOut);
	free(logOut);
	free(errOut);

	/* --- Threads failing at once: one whole record -- ------------- */

	for (int round = 0; round < 200; round++) {
		pthread_t thread[NTHREADS];
		int       id[NTHREADS], t, code;

		atomic_store(&started, 0);
		for (t = 0; t < NTHREADS; t++) {
			id[t] = t;
			pthread_create(&thread[t], NULL, failAtOnce, &id[t]);
		}
		for (t = 0; t < NTHREADS; t++) pthread_join(thread[t], NULL);

		CHECK(Failed(rh) && rh->error.message &&
				  sscanf(rh->error.message, "thread %d code %d", &t, &code) == 2 &&
				  code == (int) rh->error.code && code == RH_E_INPUT + t % 4,
			  "round %d: one consistent record (%d, %s)", round, rh->error.code,
			  rh->error.message ? rh->error.message : "NULL");
		ClearError(rh);
	}

	return check_summary("test_status");
}
