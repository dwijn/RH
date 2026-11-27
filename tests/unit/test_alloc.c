/* ------- file: -------------------------- test_alloc.c ------------

	   Allocations that cannot fail quietly (alloc.c, matrix.c, Q-53):

	   -- mallocOrAbort, callocOrAbort, reallocOrAbort and strdupOrAbort
		  return what malloc, calloc, realloc and strdup return when
		  those succeed, and a zero size is not a failure;
	   -- a failed allocation (a size no machine has) prints
		  "-OUT_OF_MEMORY in routine R" with the size and the
		  expression to stderr and aborts (SIGABRT, in a child), for
		  each helper and for matrix_double.

	   The sanitizer build runs this with allocator_may_return_null=1
	   (tests/unit/CMakeLists.txt), so that ASan lets malloc fail.
	   --                                              -------------- */

#include <signal.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>

#include "rh.h"
#include "check.h"

#define HUGE_SIZE (SIZE_MAX / 4)

/* --- What a child writes to stderr, and how it ends: its exit
	   status, or 128 + the signal that killed it -- -------------- */

static int inChild(void (*body)(void), char **errOut)
{
	int     perr[2], status;
	size_t  size = 0, cap = 4096;
	ssize_t n;
	pid_t   pid;

	*errOut = calloc(cap, 1);
	if (pipe(perr)) return -1;
	fflush(NULL);
	if ((pid = fork()) == 0) {
		close(perr[0]);
		dup2(perr[1], STDERR_FILENO);
		body();
		_exit(0);
	}
	close(perr[1]);
	while ((n = read(perr[0], *errOut + size, cap - size - 1)) > 0) {
		size += (size_t) n;
		if (cap - size < 2) *errOut = realloc(*errOut, cap *= 2);
	}
	(*errOut)[size] = '\0';
	close(perr[0]);
	waitpid(pid, &status, 0);
	if (WIFEXITED(status)) return WEXITSTATUS(status);
	return WIFSIGNALED(status) ? 128 + WTERMSIG(status) : -1;
}

static void hugeMalloc(void)
{
	volatile size_t size = HUGE_SIZE;
	free(mallocOrAbort(size));
}

static void hugeCalloc(void)
{
	volatile size_t n = HUGE_SIZE;
	free(callocOrAbort(n, sizeof(double)));
}

static void hugeRealloc(void)
{
	volatile size_t size = HUGE_SIZE;
	free(reallocOrAbort(mallocOrAbort(8), size));
}

static void hugeMatrix(void)
{
	freeMatrix(rh, (void **) matrix_double(rh, 1 << 30, 1 << 30));
}

int main(void)
{
	char   *err, *s;
	double *x;
	int     status;

	/* --- Success: what the plain calls give --        -------------- */

	x = mallocOrAbort(10 * sizeof(double));
	CHECK(x != NULL, "mallocOrAbort");
	x = reallocOrAbort(x, 20 * sizeof(double));
	CHECK(x != NULL, "reallocOrAbort");
	free(x);
	x = callocOrAbort(10, sizeof(double));
	CHECK(x != NULL && x[0] == 0.0 && x[9] == 0.0, "callocOrAbort zeroes");
	free(x);
	s = strdupOrAbort("RH");
	CHECK(s != NULL && strcmp(s, "RH") == 0, "strdupOrAbort");
	free(s);
	free(mallocOrAbort(0));
	free(callocOrAbort(0, 8));
	CHECK(1, "a zero size does not abort");

	/* --- Failure: the message, and SIGABRT --          -------------- */

	status = inChild(hugeMalloc, &err);
	CHECK(status == 128 + SIGABRT, "mallocOrAbort: status %d", status);
	CHECK(strstr(err, "-OUT_OF_MEMORY in routine hugeMalloc\n Could not allocate ") != NULL &&
			  strstr(err, " bytes (size) \n Aborting.....\n") != NULL,
		  "mallocOrAbort: stderr %s", err);
	free(err);

	status = inChild(hugeCalloc, &err);
	CHECK(status == 128 + SIGABRT && strstr(err, "routine hugeCalloc") != NULL &&
			  strstr(err, "(n x sizeof(double))") != NULL,
		  "callocOrAbort: status %d, stderr %s", status, err);
	free(err);

	status = inChild(hugeRealloc, &err);
	CHECK(status == 128 + SIGABRT && strstr(err, "Could not reallocate to") != NULL,
		  "reallocOrAbort: status %d, stderr %s", status, err);
	free(err);

	status = inChild(hugeMatrix, &err);
	CHECK(status == 128 + SIGABRT &&
			  strstr(err, "-OUT_OF_MEMORY in routine matrix_double\n Could not allocate a "
						  "1073741824 x 1073741824 matrix") != NULL,
		  "matrix_double: status %d, stderr %s", status, err);
	free(err);

	return check_summary("test_alloc");
}
