/* ------- file: -------------------------- test_readatom.c ---------

	   readAtom on a model atom file that does not open stops with
	   "Unable to open input file" and the path (B-74; it printed the
	   path alone). In a run, readAtomicModels opens the file first, in
	   getAtomID, which stops with its own message, so readAtom is
	   called here directly, in a child process.
	   --                                              -------------- */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>

#include "rh.h"
#include "atom.h"
#include "error.h"
#include "inputs.h"
#include "state.h"
#include "check.h"

#define MISSING "no_such_directory/no_such.atom"

int main(void)
{
	static Atom atom;
	char        text[4096];
	int         perr[2], status;
	ssize_t     n, size = 0;
	pid_t       pid;

	CHECK(pipe(perr) == 0, "pipe failed");
	fflush(NULL);
	if ((pid = fork()) == 0) {
		close(perr[0]);
		dup2(perr[1], STDERR_FILENO);
		rh->commandline.logfile = fopen("/dev/null", "w");
		readAtom(rh, &atom, MISSING);
		_exit(0); /* not reached */
	}
	close(perr[1]);
	while (size < (ssize_t) sizeof(text) - 1 &&
		   (n = read(perr[0], text + size, sizeof(text) - 1 - size)) > 0)
		size += n;
	text[size] = '\0';
	close(perr[0]);
	waitpid(pid, &status, 0);

	CHECK(WIFEXITED(status) && WEXITSTATUS(status) == RH_EXIT_FATAL,
		  "readAtom on a missing file: status %d, expected %d",
		  WIFEXITED(status) ? WEXITSTATUS(status) : -1, RH_EXIT_FATAL);
	CHECK(strstr(text, "-TERMINATING_ERROR in routine readAtom\n"
					   " Unable to open input file " MISSING " \n") != NULL,
		  "readAtom on a missing file: stderr %s", text);

	return check_summary("test_readatom");
}
