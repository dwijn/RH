/* ------- file: -------------------------- test_error.c ------------

	   Messages, warnings and terminating errors (error.c, Q-48).
	   Error(level, routine, format, ...) formats its message as printf
	   does, without a global buffer, and writes it with one call:

	   -- a MESSAGE is the formatted text, exactly, and nothing when the
		  command line asks for quiet;
	   -- a WARNING is "\n-WARNING in routine R\n TEXT\n";
	   -- a message far longer than the old 512-byte buffers comes out
		  whole (B-14, B-41);
	   -- warnings from several threads at once each come out whole and
		  with their own text (T-1: they shared messageStr);
	   -- ERROR_LEVEL_1 and ERROR_LEVEL_2 print "\a\n\n-TERMINATING_ERROR
		  in routine R\n TEXT \n Exiting.....\n" to the log file and to
		  stderr and exit with status RH_EXIT_FATAL (in a child process).
	   --                                              -------------- */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"
#include "state.h"
#include "check.h"

/* --- The log goes to a memory stream; the text written so far -- */

static char  *logText;
static size_t logSize;

static void logOpen(void)
{
	rh->commandline.logfile = open_memstream(&logText, &logSize);
}

static char *logClose(void)
{
	fclose(rh->commandline.logfile);
	rh->commandline.logfile = stderr;
	return logText;
}

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

/* --- A terminating error in a child: its exit status and what it
	   wrote to the log file and to stderr --          -------------- */

static int fatalChild(enum errorlevel level, char **logOut, char **errOut)
{
	int   plog[2], perr[2], status;
	pid_t pid;

	*logOut = *errOut = NULL;
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
		Error(level, "testRoutine", "value %d out of range [%g, %g] in %s",
			  12, 0.5, 2.0, "test");
		_exit(0); /* not reached */
	}
	close(plog[1]);
	close(perr[1]);
	*logOut = readAll(plog[0]);
	*errOut = readAll(perr[0]);
	close(plog[0]);
	close(perr[0]);
	waitpid(pid, &status, 0);
	return WIFEXITED(status) ? WEXITSTATUS(status) : -1;
}

/* --- Warnings from several threads at once --       -------------- */

#define NTHREADS 4
#define NWARN    500

static void *warnMany(void *arg)
{
	int t = *(int *) arg;

	for (int n = 0; n < NWARN; n++)
		Error(WARNING, "warnMany", "thread %d warning %d: %s", t, n,
			  "the same long tail on every line, so that a shared buffer would show");
	return NULL;
}

int main(void)
{
	char *text, *expect, *logOut, *errOut, *longText;
	int   status;

	/* --- MESSAGE: the formatted text as it is --      -------------- */

	logOpen();
	Error(MESSAGE, NULL, "%d %s %.3f%%\n", 7, "x", 1.5);
	Error(MESSAGE, "ignored", "no arguments\n");
	text = logClose();
	CHECK(strcmp(text, "7 x 1.500%\nno arguments\n") == 0, "message %s", text);
	free(text);

	rh->commandline.quiet = TRUE;
	logOpen();
	Error(MESSAGE, NULL, "%s\n", "quiet");
	text = logClose();
	CHECK(strlen(text) == 0, "quiet message %s", text);
	free(text);
	rh->commandline.quiet = FALSE;

	/* --- WARNING --                                    -------------- */

	logOpen();
	Error(WARNING, "testRoutine", "%d of %d", 3, 4);
	text = logClose();
	CHECK(strcmp(text, "\n-WARNING in routine testRoutine\n 3 of 4\n") == 0,
		  "warning %s", text);
	free(text);

	/* --- A message of 5000 characters (the buffers were 512) -- -- */

	longText = malloc(5001);
	for (int k = 0; k < 5000; k++) longText[k] = 'a' + k % 26;
	longText[5000] = '\0';
	logOpen();
	Error(WARNING, "testRoutine", "%s|%d", longText, 42);
	text   = logClose();
	expect = malloc(5100);
	sprintf(expect, "\n-WARNING in routine testRoutine\n %s|42\n", longText);
	CHECK(strcmp(text, expect) == 0, "long warning: %zu characters, %zu expected",
		  strlen(text), strlen(expect));
	free(text);
	free(expect);

	/* --- Threads: every warning whole, with its own text -- ------ */

	{
		pthread_t thread[NTHREADS];
		int       id[NTHREADS], seen[NTHREADS][NWARN] = {{0}}, bad = 0;
		char     *line, *save;

		logOpen();
		for (int t = 0; t < NTHREADS; t++) {
			id[t] = t;
			pthread_create(&thread[t], NULL, warnMany, &id[t]);
		}
		for (int t = 0; t < NTHREADS; t++) pthread_join(thread[t], NULL);
		text = logClose();

		/* --- Each warning is three lines: empty, header, text -- */

		for (line = strtok_r(text, "\n", &save); line; line = strtok_r(NULL, "\n", &save)) {
			int t, n;
			if (strcmp(line, "-WARNING in routine warnMany") == 0) continue;
			if (sscanf(line, " thread %d warning %d:", &t, &n) == 2 &&
				t >= 0 && t < NTHREADS && n >= 0 && n < NWARN &&
				strcmp(strchr(line, ':'),
					   ": the same long tail on every line, so that a shared buffer would show") == 0)
				seen[t][n]++;
			else
				bad++;
		}
		for (int t = 0; t < NTHREADS; t++)
			for (int n = 0; n < NWARN; n++)
				if (seen[t][n] != 1) bad++;
		CHECK(bad == 0, "%d lines or warnings wrong of %d warnings from %d threads",
			  bad, NTHREADS * NWARN, NTHREADS);
		free(text);
	}

	/* --- The error levels: message, stderr, exit status -- ------- */

	expect = "\a\n\n-TERMINATING_ERROR in routine testRoutine\n"
			 " value 12 out of range [0.5, 2] in test \n Exiting.....\n";
	for (enum errorlevel level = ERROR_LEVEL_1; level <= ERROR_LEVEL_2; level++) {
		status = fatalChild(level, &logOut, &errOut);
		CHECK(status == RH_EXIT_FATAL, "level %d: exit status %d", level, status);
		CHECK(strcmp(logOut, expect) == 0, "level %d: log %s", level, logOut);
		CHECK(strcmp(errOut, expect) == 0, "level %d: stderr %s", level, errOut);
		free(logOut);
		free(errOut);
	}
	free(longText);

	return check_summary("test_error");
}
