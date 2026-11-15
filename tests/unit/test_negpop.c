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

extern Atmosphere  atmos;
extern InputData   input;
extern CommandLine commandline;

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

static void run_stop(void)
{
	statEquil(&atom, -1);
}

/* --- Run f in a child with its output discarded; TRUE if it ended
	   with RH's fatal exit status --                   -------------- */

static bool_t is_fatal(void (*f)(void))
{
	pid_t pid;
	int   status;

	fflush(stdout);
	fflush(stderr);
	if ((pid = fork()) == 0) {
		in_child = 1;
		if (!freopen("/dev/null", "w", stdout) || !freopen("/dev/null", "w", stderr)) _exit(1);
		commandline.logfile = stderr;
		f();
		_exit(0);
	}
	if (pid < 0 || waitpid(pid, &status, 0) != pid) return FALSE;
	return WIFEXITED(status) && WEXITSTATUS(status) == RH_EXIT_FATAL;
}

static void check_case(const char *name, double a, double b, bool_t bad,
					   const char *level)
{
	char msg[512];

	set_rates(a, b);
	CHECK(is_fatal(run_stop) == bad, "%s: outside a trial, stops: %d, expected %d",
		  name, !bad, bad);

	set_rates(a, b);
	negTauReset();
	negTauTrial(TRUE);
	statEquil(&atom, -1);
	negTauTrial(FALSE);
	CHECK(negTauFound() == bad, "%s: in a trial, invalid: %d, expected %d", name,
		  negTauFound(), bad);
	if (bad) {
		negTauMessage(msg, sizeof(msg));
		CHECK(strstr(msg, "XX") && strstr(msg, level) &&
				  strstr(msg, "location 1") &&
				  strstr(msg, "negative preconditioned rates"),
			  "%s: message: %s", name, msg);
		negTauReset();
		CHECK(!negTauFound(), "%s: still invalid after negTauReset", name);
	}
}

int main(void)
{
	commandline.logfile = fopen("/dev/null", "w");
	input.Nthreads      = 1;
	atmos.Nspace        = NSPACE;
	memset(&atom, 0, sizeof(atom));
	strcpy(atom.ID, "XX");
	atom.Nlevel   = N;
	atom.label    = (char **) malloc(N * sizeof(char *));
	atom.label[0] = strdup("first");
	atom.label[1] = strdup("second");
	atom.Gamma    = matrix_double(N * N, NSPACE);
	atom.C        = matrix_double(N * N, NSPACE);
	atom.n        = matrix_double(N, NSPACE);
	atom.ntotal   = (double *) malloc(NSPACE * sizeof(double));

	/* --- n_1 = a / (a + b) N: -1 for a = -0.5, b = 1 -- ------------ */

	check_case("negative rate, negative population", -0.5, 1.0, TRUE, "level 1 (second)");
	check_case("NaN rate", NAN, 1.0, TRUE, "non-finite population");

	/* --- Plain non-negative rates: GTH, nothing to report -- ------ */

	check_case("non-negative rates", 0.5, 1.0, FALSE, "");

	free(atom.label[0]);
	free(atom.label[1]);
	free(atom.label);
	freeMatrix((void **) atom.Gamma);
	freeMatrix((void **) atom.C);
	freeMatrix((void **) atom.n);
	free(atom.ntotal);
	return check_summary("test_negpop");
}
