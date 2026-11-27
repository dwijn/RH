/* ------- file: -------------------------- test_negpop.c -----------

	   The backstop of statEquil for the LU fallback (N-8): populations
	   that come out of LU negative or not finite.

	   A two-level atom with a negative rate 0 -> 1 has the LU
	   population n_1 = -N a / (|a| ... ) < 0 (GTH needs non-negative
	   rates). The test checks

	   -- outside a trial evaluation of G, statEquil stops the run
		  (ERROR_LEVEL_2, in a child process);
	   -- inside one (negTauTrial(TRUE)) it returns, and negTauFound is
		  TRUE with a message that names atom, level and location and
		  says that negative populations come from negative
		  preconditioned rates, until negTauReset;
	   -- the same for a non-finite rate (a NaN, from LU a NaN);
	   -- a negative rate with non-negative LU populations, and plain
		  non-negative rates, give no stop and no invalidation.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "error.h"
#include "check.h"
#include "negtau.h"
#include "state.h"

/* --- LeakSanitizer would report the buffers that Error() leaves in
	   the child that it ends, and change its exit status -- -------- */

static int in_child = 0;

int __lsan_is_turned_off(void);
int __lsan_is_turned_off(void)
{
	return in_child;
}

#define N      2
#define NSPACE 3

static Atom atom;

/* --- Rates (0 -> 1, 1 -> 0) a, b at location 1, 1.0 and 1.0 at the
	   others --                                       -------------- */

static void set_rates(double a, double b)
{
	for (int k = 0; k < NSPACE; k++) {
		atom.Gamma[1 * N + 0][k] = (k == 1) ? a : 1.0; /* 0 -> 1 */
		atom.Gamma[0 * N + 1][k] = (k == 1) ? b : 1.0; /* 1 -> 0 */
		atom.n[0][k]             = 0.5;
		atom.n[1][k]             = 0.5;
		atom.ntotal[k]           = 1.0;
	}
}

static void run_invalid_outside(void)
{
	negTauInvalid(rh, "x");
}

static void run_stop(void)
{
	statEquil(rh, &atom, -1);
}

/* --- Run f in a child with its output discarded; TRUE if it ended
	   with RH's fatal exit status --                   -------------- */

static int is_fatal(void (*f)(void))
{
	pid_t pid;
	int   status;

	fflush(stdout);
	fflush(stderr);
	if ((pid = fork()) == 0) {
		in_child = 1;
		if (!freopen("/dev/null", "w", stdout) || !freopen("/dev/null", "w", stderr)) _exit(1);
		rh->commandline.logfile = stderr;
		f();
		_exit(0);
	}
	if (pid < 0 || waitpid(pid, &status, 0) != pid) return -1; /* fork failed: not "not fatal" */
	return WIFEXITED(status) && WEXITSTATUS(status) == RH_EXIT_FATAL;
}

static void check_case(const char *name, double a, double b, bool_t bad,
					   const char *level, const char *tail)
{
	char msg[512];

	set_rates(a, b);
	{
		int fatal = is_fatal(run_stop);
		CHECK(fatal == bad, "%s: outside a trial, stops: %d (-1: fork failed), expected %d",
			  name, fatal, bad);
	}

	set_rates(a, b);
	negTauReset(rh);
	negTauTrial(rh, TRUE);
	statEquil(rh, &atom, -1);
	negTauTrial(rh, FALSE);
	CHECK(negTauFound(rh) == bad, "%s: in a trial, invalid: %d, expected %d", name,
		  negTauFound(rh), bad);
	if (bad) {
		negTauMessage(rh, msg, sizeof(msg));
		CHECK(strstr(msg, "XX") && strstr(msg, level) &&
				  strstr(msg, "location 1") &&
				  strstr(msg, "of 3 locations solved") && strstr(msg, tail),
			  "%s: message: %s", name, msg);
		negTauReset(rh);
		CHECK(!negTauFound(rh), "%s: still invalid after negTauReset", name);
	}
}

int main(void)
{
	rh->commandline.logfile = fopen("/dev/null", "w");
	rh->input.Nthreads      = 1;
	rh->atmos.Nspace        = NSPACE;
	memset(&atom, 0, sizeof(atom));
	strcpy(atom.ID, "XX");
	atom.Nlevel   = N;
	atom.label    = (char **) malloc(N * sizeof(char *));
	atom.label[0] = strdup("first");
	atom.label[1] = strdup("second");
	atom.Gamma    = matrix_double(rh, N * N, NSPACE);
	atom.C        = matrix_double(rh, N * N, NSPACE);
	atom.n        = matrix_double(rh, N, NSPACE);
	atom.ntotal   = (double *) malloc(NSPACE * sizeof(double));

	/* --- n_1 = a / (a + b) N: -1 for a = -0.5, b = 1 -- ------------ */

	check_case("negative rate, negative population", -0.5, 1.0, TRUE, "level 1 (second)", "negative rates");
	check_case("NaN rate", NAN, 1.0, TRUE, "non-finite population", "non-finite rate");

	/* --- Plain non-negative rates: GTH, nothing to report -- ------ */

	CHECK(is_fatal(run_invalid_outside) == 1,
		  "negTauInvalid outside a trial does not stop (or fork failed)");
	check_case("non-negative rates", 0.5, 1.0, FALSE, "", "");

	free(atom.label[0]);
	free(atom.label[1]);
	free(atom.label);
	freeMatrix(rh, (void **) atom.Gamma);
	freeMatrix(rh, (void **) atom.C);
	freeMatrix(rh, (void **) atom.n);
	free(atom.ntotal);
	return check_summary("test_negpop");
}
