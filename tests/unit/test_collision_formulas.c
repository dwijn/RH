/* ------- file: ------------------- test_collision_formulas.c -------

	   Collision rates given by formulas (N-34): CollisionRate with
	   synthetic atom files.

	   -- the formula records of Atoms/C_I+II_9.atom (SHULL82, BADNELL,
		  AR85-CDI, AR85-CEA, SUMMERS) give finite, non-negative rates at
		  61 temperatures from 1e3 to 1e6 K, without a message;
	   -- a formula that is negative at some T (SHULL82 with B < -1,
		  BADNELL with a negative coefficient, AR85-CDI with a negative
		  sum, AR85-CHH with 1 + c exp(d T/1e4) < 0) gives zero there,
		  with one message per record and run, and the rates elsewhere
		  are unchanged;
	   -- coefficients that cannot be right (not finite, negative
		  energy, T1 > T2, ...) stop the run when read;
	   -- AR85-CHP takes the proton density from the last level of the
		  hydrogen model, whatever its number of levels (it was level 5),
		  and AR85-CHP and AR85-CHH report the locations outside their
		  fit range.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"
#include "check.h"

double summers(int i, int j, double nne, struct Atom *atom);

extern Atmosphere  atmos;
extern InputData   input;
extern CommandLine commandline;

static int in_child = 0;

int __lsan_is_turned_off(void);
int __lsan_is_turned_off(void)
{
	return in_child;
}

#define NLEVEL 9
#define NT     61

static Atom atom, hydrogen;
static int  verbose = 0;

static const double g_lev[NLEVEL]     = {9, 5, 1, 2, 4, 12, 6, 4, 1};
static const int    stage_lev[NLEVEL] = {0, 0, 0, 1, 1, 1, 1, 1, 2};

static const char *elementID[28] = {"H", "HE", "LI", "BE", "B", "C", "N", "O", "F", "NE",
									"NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA",
									"SC", "TI", "V", "CR", "MN", "FE", "CO", "NI"};

/* --- The formula records of Atoms/C_I+II_9.atom -- --------------- */

static const char *c_records =
	"SUMMERS 1.0\n"
	"SHULL82   0   3   0.00e+00 1.31e+05 4.70e-13 6.24e-01 0.00e-03 4.42e-02 1.57e+05 3.74e+05\n"
	"SHULL82   0   4   0.00e+00 1.31e+05 4.70e-13 6.24e-01 0.00e-03 4.42e-02 1.57e+05 3.74e+05\n"
	"BADNELL   0   3   7\n"
	"0.000e+00  1.2170e+01  7.3800e+01  1.5230e+04  1.2070e+05  2.1440e+05  0.000e+00\n"
	"0.000e+00  6.3460e-09  9.7930e-09  1.6340e-06  8.3690e-04  3.3550e-04  0.000e+00\n"
	"BADNELL   0   4   7\n"
	"0.000e+00  1.2170e+01  7.3800e+01  1.5230e+04  1.2070e+05  2.1440e+05  0.000e+00\n"
	"0.000e+00  6.3460e-09  9.7930e-09  1.6340e-06  8.3690e-04  3.3550e-04  0.000e+00\n"
	"SHULL82   3   8   0.00e+00 2.83e+05 1.38e-12 6.45e-01 0.00e+00 5.88e-02 1.41e+05 1.41e+05\n"
	"SHULL82   4   8   0.00e+00 2.83e+05 9.20e-13 6.45e-01 0.00e+00 5.88e-02 1.41e+05 1.41e+05\n"
	"BADNELL   3   8   7\n"
	"0.000e+00 2.6600e+03  3.7560e+03  2.5660e+04  1.4000e+05  1.8010e+06  4.3070e+06\n"
	"0.000e+00 2.0934e-06  1.3332e-07  1.1724e-05  2.5272e-03  1.2222e-04  1.7616e-04\n"
	"BADNELL   4   8   8\n"
	"0.000e+00  2.6600e+03  3.7560e+03  2.5660e+04  1.4000e+05  1.8010e+06  4.3070e+06  0.000e+00\n"
	"0.000e+00  1.3956e-06  8.8880e-08  7.8160e-06  1.6848e-03  8.1480e-05  1.1744e-04  0.000e+00\n"
	"AR85-CDI   0   3   2\n"
	"    11.30     3.60    -9.60     7.20    -9.06\n"
	"    16.60    14.58    -4.68     1.50   -14.40\n"
	"AR85-CEA   0   3  6.00e-01\n"
	"AR85-CDI   0   4   2\n"
	"    11.30     2.40    -6.40     4.80    -6.04\n"
	"    16.60     9.72    -3.12     1.00    -9.60\n"
	"AR85-CEA   0   4  4.00e-01\n"
	"AR85-CDI   3   8   2\n"
	"    24.40    17.30    -9.00     2.50   -13.60\n"
	"    30.90    23.70    -7.60     2.50   -21.70\n"
	"AR85-CEA  3   8  1.00e+00\n"
	"AR85-CDI    4   8   2\n"
	"   24.39    17.30    -9.00     2.50   -13.60\n"
	"    30.89    23.70    -7.60     2.50   -21.70\n"
	"AR85-CEA   4   8  1.00e+00\n"
	"END\n";

/* --- Run CollisionRate on a file with the records; the log goes to
	   logPath (empty file when no message) -- ------------------- */

static void run(const char *records)
{
	FILE *fp = tmpfile();

	fputs(records, fp);
	fputs("END\n", fp);
	rewind(fp);
	CollisionRate(&atom, fp);
	fclose(fp);
}

static char logText[8192];

static void run_logged(const char *records)
{
	FILE *log = tmpfile();
	FILE *old = commandline.logfile;
	long  n;

	commandline.logfile = log;
	run(records);
	commandline.logfile = old;
	fflush(log);
	n = ftell(log);
	if (n >= (long) sizeof(logText)) n = sizeof(logText) - 1;
	rewind(log);
	n          = (long) fread(logText, 1, n, log);
	logText[n] = '\0';
	fclose(log);
}

static int count(const char *needle)
{
	int         n = 0;
	const char *p = logText;

	while ((p = strstr(p, needle)) != NULL) {
		n++;
		p += strlen(needle);
	}
	return n;
}

static int is_fatal(const char *records)
{
	pid_t pid;
	int   status;

	fflush(stdout);
	fflush(stderr);
	if ((pid = fork()) == 0) {
		in_child = 1;
		if (!freopen("/dev/null", "w", stdout) || !freopen("/dev/null", "w", stderr)) _exit(1);
		commandline.logfile = stderr;
		run(records);
		_exit(0);
	}
	if (pid < 0 || waitpid(pid, &status, 0) != pid) return -1;
	return WIFEXITED(status) && WEXITSTATUS(status) == RH_EXIT_FATAL;
}

static void set_temperatures(double Tlo, double Thi, int n)
{
	atmos.Nspace = n;
	for (int k = 0; k < n; k++) {
		atmos.T[k]  = (n > 1) ? Tlo * pow(Thi / Tlo, (double) k / (n - 1)) : Tlo;
		atmos.ne[k] = 1e18;
	}
}

static void check_all_finite_nonneg(const char *what)
{
	int bad = 0;

	for (int ij = 0; ij < NLEVEL * NLEVEL; ij++) {
		double mn = 1e300, mx = 0.0;

		for (int k = 0; k < atmos.Nspace; k++) {
			double c = atom.C[ij][k];

			if (!isfinite(c) || c < 0.0) {
				if (bad++ < 5)
					fprintf(stderr, "%s: C[%d][%d] at T = %g: %g\n", what, ij / NLEVEL,
							ij % NLEVEL, atmos.T[k], c);
			}
			if (c < mn) mn = c;
			if (c > mx) mx = c;
		}
		if (verbose && mx > 0.0)
			fprintf(stderr, "%s: C[%d][%d] (%d <- %d) min %.3e max %.3e\n", what, ij / NLEVEL, ij % NLEVEL,
					ij / NLEVEL, ij % NLEVEL, mn, mx);
	}
	CHECK(bad == 0, "%s: %d rates negative or not finite", what, bad);
}

static void setup(void)
{
	commandline.logfile = fopen("/dev/null", "w");
	input.Nthreads      = 1;
	atmos.Nspace        = NT;
	atmos.T             = (double *) calloc(NT, sizeof(double));
	atmos.ne            = (double *) calloc(NT, sizeof(double));
	atmos.elements      = (Element *) calloc(30, sizeof(Element));
	for (int n = 0; n < 28; n++) strcpy(atmos.elements[n].ID, elementID[n]);

	memset(&atom, 0, sizeof(atom));
	strcpy(atom.ID, "C");
	atom.Nlevel = NLEVEL;
	atom.g      = (double *) malloc(NLEVEL * sizeof(double));
	atom.E      = (double *) calloc(NLEVEL, sizeof(double));
	atom.stage  = (int *) malloc(NLEVEL * sizeof(int));
	atom.nstar  = matrix_double(NLEVEL, NT);
	for (int i = 0; i < NLEVEL; i++) {
		atom.g[i]     = g_lev[i];
		atom.stage[i] = stage_lev[i];
		for (int k = 0; k < NT; k++) atom.nstar[i][k] = g_lev[i];
	}

	/* --- Hydrogen of 8 levels, so that level 5 is not the continuum;
		   a bound level is 1e10 and 2e10 times level 0 and the
		   continuum (level 7) is 3e10 -- ------------------------- */

	memset(&hydrogen, 0, sizeof(hydrogen));
	hydrogen.Nlevel = 8;
	hydrogen.n      = matrix_double(8, NT);
	for (int i = 0; i < 8; i++)
		for (int k = 0; k < NT; k++)
			hydrogen.n[i][k] = (i == 7) ? 3e10 : (i == 5) ? 2e10
														  : 1e10 + i;
	atmos.H = &hydrogen;
}

static void change_hydrogen(int Nlevel, double np)
{
	freeMatrix((void **) hydrogen.n);
	hydrogen.Nlevel = Nlevel;
	hydrogen.n      = matrix_double(Nlevel, NT);
	for (int i = 0; i < Nlevel; i++)
		for (int k = 0; k < NT; k++)
			hydrogen.n[i][k] = (i == Nlevel - 1) ? np : 1e10 + i;
}

int main(int argc, char **argv)
{
	char rec[2048];

	verbose = (argc > 1 && !strcmp(argv[1], "-v"));
	setup();

	/* --- The formula records of C_I+II_9, 1e3 to 1e6 K -- --------- */

	set_temperatures(1e3, 1e6, NT);
	run_logged(c_records);
	check_all_finite_nonneg("C_I+II_9");
	CHECK(strlen(logText) == 0, "C_I+II_9: message: %s", logText);

	/* --- SHULL82 with B = -3 (negative from T1 / ln 3 = 1.05e5 K up);
		   the dielectronic term dominates there: A_rad = 0 -- ----- */

	strcpy(atom.ID, "N");
	snprintf(rec, sizeof(rec), "SHULL82 0 3 0 1e5 0 0.5 1e-10 -3 1e4 1e5\n");
	{
		int nzero = 0;

		set_temperatures(1e3, 1e6, NT);
		run_logged(rec);
		CHECK(count("negative rate") == 1, "SHULL82 B = -3: messages: %d\n%s", count("negative rate"), logText);
		CHECK(strstr(logText, "SHULL82 0 3") != NULL, "SHULL82 message names the record");
		for (int k = 0; k < NT; k++) {
			double T      = atmos.T[k];
			double sm     = summers(0, 3, atmos.ne[k], &atom);
			double expect = sm * 1e-10 / T / sqrt(T) * exp(-1e4 / T) * (1.0 - 3.0 * exp(-1e5 / T));
			expect        = (expect < 0.0) ? 0.0 : expect * atmos.ne[k] * 1e-6; /* cm3 -> m3 */
			CHECK(atom.C[0 * NLEVEL + 3][k] >= 0.0, "SHULL82: C[0][3] at T = %g: %g", T, atom.C[0 * NLEVEL + 3][k]);
			CHECK_CLOSE(atom.C[0 * NLEVEL + 3][k], expect, 1e-9, 0.0, "SHULL82 at T = %g", T);
			nzero += (expect == 0.0);
		}
		CHECK(nzero > 5 && nzero < NT, "SHULL82: %d locations clipped", nzero);

		/* --- Called again (electron density changed): no second message -- */

		run_logged(rec);
		CHECK(count("negative rate") == 0, "second call repeats the message: %s", logText);
	}

	/* --- BADNELL with a negative coefficient: sum of c_i exp(-E_i/T)
		   negative at low T where the second term dominates -- ----- */

	strcpy(atom.ID, "O");
	run_logged("BADNELL 0 3 2\n0 1e4\n1e-8 -1e-6\n");
	CHECK(count("negative rate") == 1 && strstr(logText, "BADNELL 0 3"), "BADNELL: %s", logText);
	for (int k = 0; k < NT; k++) {
		double T    = atmos.T[k];
		double s    = (1e-8 - 1e-6 * exp(-1e4 / T)) * pow(T, -1.5);
		double want = (s < 0.0) ? 0.0 : s;

		/* --- C[0][3] gets the recombination coefficient only through
			   ne and the Summers factor (SUMMERS 1: density dependent), so
			   check the sign and the zero only -- ------------------ */

		CHECK(atom.C[0 * NLEVEL + 3][k] >= 0.0, "BADNELL at T = %g", T);
		CHECK((want == 0.0) == (atom.C[0 * NLEVEL + 3][k] == 0.0), "BADNELL zero at T = %g: %g", T,
			  atom.C[0 * NLEVEL + 3][k]);
	}

	/* --- AR85-CDI with a sum that is negative at all T: a single shell
		   with a < 0 -- ----------------------------------------- */

	strcpy(atom.ID, "NE");
	run_logged("AR85-CDI 0 3 1\n 11.3 -1.0 0.0 0.0 0.0\n");
	CHECK(count("negative rate") == 1 && strstr(logText, "AR85-CDI 0 3"), "AR85-CDI: %s", logText);
	for (int k = 0; k < NT; k++)
		CHECK(atom.C[3 * NLEVEL + 0][k] == 0.0 && atom.C[0 * NLEVEL + 3][k] == 0.0, "AR85-CDI clipped at T = %g", atmos.T[k]);

	/* --- AR85-CHH, fit range 2e3 to 2e5, c < -exp(-d t4) for
		   t4 > 0.1 (d = 5 -> c = -3: 1 - 3 exp(5 t4) < 0 always) -- */

	strcpy(atom.ID, "MG");
	run_logged("AR85-CHH 0 3 2e3 2e5 1.0 0.5 -3.0 5.0\n");
	CHECK(count("negative rate") == 1 && strstr(logText, "AR85-CHH 0 3"), "AR85-CHH negative: %s", logText);
	CHECK(count("only defined for") == 1, "AR85-CHH range: %s", logText);
	for (int k = 0; k < NT; k++) CHECK(atom.C[0 * NLEVEL + 3][k] == 0.0, "AR85-CHH at T = %g", atmos.T[k]);

	/* --- AR85-CHH positive: unchanged, zero outside the range, with a
		   message for the range only -- ------------------------- */

	strcpy(atom.ID, "SI");
	run_logged("AR85-CHH 0 3 2e3 2e5 1.0 0.5 0.1 -1.0\n");
	CHECK(count("negative rate") == 0 && count("only defined for") == 1, "AR85-CHH positive: %s", logText);
	for (int k = 0; k < NT; k++) {
		double T = atmos.T[k], t4 = T / 1e4;
		double want = (T >= 2e3 && T <= 2e5)
						  ? 1.0 * 1e-9 * pow(t4, 0.5) * (1.0 + 0.1 * exp(-t4)) * hydrogen.n[0][k] * 1e-6
						  : 0.0;

		CHECK_CLOSE(atom.C[0 * NLEVEL + 3][k], want, 1e-12, 0.0, "AR85-CHH at T = %g", T);
	}

	/* --- AR85-CHP: proton density is the last level of the hydrogen
		   model (8 levels here: level 5 differs), and for a hydrogen of
		   two levels the read is within the array -- ---------------- */

	strcpy(atom.ID, "S");
	for (int nl = 8; nl >= 2; nl -= 6) {
		change_hydrogen(nl, 4e10);
		set_temperatures(1e3, 1e6, NT);
		run_logged("AR85-CHP 0 3 5e3 5e5 2.0 0.3 0.2 0.1\n");
		CHECK(count("negative rate") == 0, "AR85-CHP: %s", logText);
		for (int k = 0; k < NT; k++) {
			double T = atmos.T[k], t4 = T / 1e4;
			double want = (T >= 5e3 && T <= 5e5)
							  ? 2.0 * 1e-9 * pow(t4, 0.3) * exp(-0.2 * t4) * exp(-0.1 * EV / KBOLTZMANN / T) *
									4e10 * 1e-6
							  : 0.0;

			CHECK_CLOSE(atom.C[3 * NLEVEL + 0][k], want, 1e-12, 0.0, "AR85-CHP, %d levels, at T = %g", nl, T);
		}
		strcpy(atom.ID, nl == 8 ? "CA" : "FE");
	}

	/* --- Read-time stops -- ------------------------------------ */

	CHECK(is_fatal("SHULL82 0 3 0 1e5 1e-13 0.5 1e-10 -3 1e4 1e5\n") == 0, "SHULL82 B = -3 is not a read-time error");
	CHECK(is_fatal("SHULL82 0 3 0 1e5 1e-13 0.5 1e-10 nan 1e4 1e5\n") == 1, "SHULL82 NaN");
	CHECK(is_fatal("SHULL82 0 3 0 1e5 -1e-13 0.5 1e-10 0 1e4 1e5\n") == 1, "SHULL82 negative arad");
	CHECK(is_fatal("SHULL82 0 3 1e-9 0 1e-13 0.5 1e-10 0 1e4 1e5\n") == 1, "SHULL82 Tcol = 0 with A > 0");
	CHECK(is_fatal("SHULL82 0 3 0 0 1e-13 0.5 1e-10 0 -1e4 1e5\n") == 1, "SHULL82 negative T0");
	CHECK(is_fatal("BADNELL 0 3 2\n0 -1e4\n1e-8 1e-6\n") == 1, "BADNELL negative energy");
	CHECK(is_fatal("BADNELL 0 3 2\n0 1e4\n1e-8 inf\n") == 1, "BADNELL infinite coefficient");
	CHECK(is_fatal("AR85-CDI 0 3 1\n 0.0 1 1 1 1\n") == 1, "AR85-CDI zero shell energy");
	CHECK(is_fatal("AR85-CDI 0 3 1\n 11.3 1 1 nan 1\n") == 1, "AR85-CDI NaN");
	CHECK(is_fatal("AR85-CEA 0 3 -0.5\n") == 1, "AR85-CEA negative factor");
	CHECK(is_fatal("BURGESS 0 3 nan\n") == 1, "BURGESS NaN");
	CHECK(is_fatal("AR85-CHP 0 3 5e5 5e3 2.0 0.3 0.2 0.1\n") == 1, "AR85-CHP T1 > T2");
	CHECK(is_fatal("AR85-CHH 0 3 5e3 5e5 2.0 0.3 inf 0.1\n") == 1, "AR85-CHH infinite");
	CHECK(is_fatal("AR85-CHH 0 3 5e3 5e5 2.0 0.3 0.2 0.1\n") == 0, "AR85-CHH valid record is not fatal");

	freeMatrix((void **) hydrogen.n);
	freeMatrix((void **) atom.nstar);
	freeMatrix((void **) atom.C);
	free(atom.g);
	free(atom.E);
	free(atom.stage);
	free(atmos.T);
	free(atmos.ne);
	free(atmos.elements);
	return check_summary("test_collision_formulas");
}
