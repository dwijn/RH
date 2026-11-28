/* ------- file: --------------------- test_collision_tables.c ------

	   Tabulated collision rates (OMEGA, CE, CI, CP, CH0, CH+, CH)
	   outside the temperature range of their table (N-53 d):
	   CollisionRate with synthetic atom files.

	   -- inside the range, a CE table that is C0 Omega / (g_i T) of a
		  power-law Omega gives the same rate as the OMEGA table of that
		  Omega (both are straight lines in log-log, which PCHIP
		  reproduces);
	   -- outside it every type holds the value at the end of its table:
		  for CE the rate then grows as sqrt(T), for OMEGA it falls as
		  T^-1/2;
	   -- one summary line per table type, atom and run, with the
		  counts, and none when all tables cover the atmosphere; no
		  repeat on the next call; again after freeRateReports.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"
#include "check.h"
#include "state.h"

#define NLEVEL 3
#define NT     41

static Atom atom, hydrogen;
static char logText[8192];

static const double g_lev[NLEVEL] = {2, 4, 4}, E_lev[NLEVEL] = {0.0, 3.0e-19, 4.0e-19};

static void run_logged(const char *records)
{
	char  *buf  = NULL, text[4096];
	size_t size = 0, n;
	FILE  *log  = open_memstream(&buf, &size), *fp;
	FILE  *old  = rh->commandline.logfile;

	n = (size_t) snprintf(text, sizeof(text), "%sEND\n", records);
	CHECK(n < sizeof(text), "%s", "records too long");
	rh->commandline.logfile = log;
	fp                      = fmemopen(text, n, "r");
	CHECK(CollisionRate(rh, &atom, fp) == RH_OK, "CollisionRate failed: %s", rh->error.message);
	fclose(fp);
	rh->commandline.logfile = old;
	fclose(log);
	snprintf(logText, sizeof(logText), "%s", buf);
	free(buf);
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

int main(void)
{
	static const double Ttab[4] = {5.0e3, 1.0e4, 2.0e4, 4.0e4};
	double              C0, Om[4], ce[4], T;
	char                rec[2048], num[64];
	int                 k, n;

	C0 = ((E_RYDBERG / sqrt(M_ELECTRON)) * PI * SQ(RBOHR)) *
		 sqrt(8.0 / (PI * KBOLTZMANN));

	rh->commandline.logfile = fopen("/dev/null", "w");
	rh->input.Nthreads      = 1;
	rh->atmos.Nspace        = NT;
	rh->atmos.T             = (double *) calloc(NT, sizeof(double));
	rh->atmos.ne            = (double *) calloc(NT, sizeof(double));
	memset(&atom, 0, sizeof(atom));
	strcpy(atom.ID, "X");
	atom.Nlevel = NLEVEL;
	atom.g      = (double *) malloc(NLEVEL * sizeof(double));
	atom.E      = (double *) malloc(NLEVEL * sizeof(double));
	atom.stage  = (int *) calloc(NLEVEL, sizeof(int));
	atom.nstar  = matrix_double(rh, NLEVEL, NT);
	for (n = 0; n < NLEVEL; n++) {
		atom.g[n] = g_lev[n];
		atom.E[n] = E_lev[n];
		for (k = 0; k < NT; k++) atom.nstar[n][k] = g_lev[n] * (1.0 + 0.1 * n);
	}
	memset(&hydrogen, 0, sizeof(hydrogen));
	hydrogen.Nlevel = 2;
	hydrogen.n      = matrix_double(rh, 2, NT);
	for (k = 0; k < NT; k++) {
		hydrogen.n[0][k] = 1.0e20;
		hydrogen.n[1][k] = 1.0e18;
	}
	rh->atmos.H = &hydrogen;

	/* --- 1e3 to 2e5 K: 5e3 to 4e4 K is the table --   ------------ */

	for (k = 0; k < NT; k++) {
		rh->atmos.T[k]  = 1.0e3 * pow(200.0, (double) k / (NT - 1));
		rh->atmos.ne[k] = 1.0e18;
	}

	/* --- Omega of a power law in T; CE = C0 Omega / (g_i T) -- ------ */

	for (n = 0; n < 4; n++) {
		Om[n] = 0.5 * pow(Ttab[n] / 1.0e4, 0.3);
		ce[n] = C0 * Om[n] / (atom.g[0] * Ttab[n]);
	}
	strcpy(rec, "TEMP 4 5e3 1e4 2e4 4e4\n");
	strcat(rec, "OMEGA 0 1");
	for (n = 0; n < 4; n++) {
		snprintf(num, sizeof(num), " %.17g", Om[n]);
		strcat(rec, num);
	}
	strcat(rec, "\nCE 0 2");
	for (n = 0; n < 4; n++) {
		snprintf(num, sizeof(num), " %.17g", ce[n]);
		strcat(rec, num);
	}
	strcat(rec, "\n");

	freeRateReports(&atom);
	run_logged(rec);
	for (k = 0; k < NT; k++) {
		T = rh->atmos.T[k];
		if (T >= 5.0e3 && T <= 4.0e4)
			CHECK_CLOSE(atom.C[0 * NLEVEL + 2][k], atom.C[0 * NLEVEL + 1][k], 1.0e-9, 0.0,
						"CE vs OMEGA of the same Omega at T = %g K", T);
	}

	/* --- Outside: C and Omega held at the end of the table -- ------ */

	for (k = 0; k < NT; k++) {
		int e;

		T = rh->atmos.T[k];
		if (T >= 5.0e3 && T <= 4.0e4) continue;
		e = (T > 4.0e4) ? 3 : 0;
		CHECK_CLOSE(atom.C[0 * NLEVEL + 2][k],
					ce[e] * rh->atmos.ne[k] * (atom.g[0] / atom.g[2]) * sqrt(T), 1.0e-9, 0.0,
					"CE end value held at T = %g K", T);
		CHECK_CLOSE(atom.C[0 * NLEVEL + 1][k], C0 * rh->atmos.ne[k] * Om[e] / (atom.g[1] * sqrt(T)),
					1.0e-9, 0.0, "OMEGA end value held at T = %g K", T);
	}

	/* --- Messages: CE one line, OMEGA none of its own beyond its table
		   (OMEGA is also outside the table: its own line) -- ------ */

	CHECK(count("CE tables") == 1, "CE messages: %d\n%s", count("CE tables"), logText);
	CHECK(count("OMEGA tables") == 1, "OMEGA messages: %d\n%s", count("OMEGA tables"), logText);
	CHECK(strstr(logText, "1 of 1 CE tables") != NULL, "counts: %s", logText);
	CHECK(strstr(logText, "1 below, 1 above") != NULL, "below, above: %s", logText);
	CHECK(strstr(logText, "T = 5000 to 40000 K") != NULL, "table range: %s", logText);
	CHECK(strstr(logText, "T = 1000 to 200000 K") != NULL, "atmosphere range: %s", logText);

	/* --- No repeat on the next call; again after a new run -- ------- */

	run_logged(rec);
	CHECK(strlen(logText) == 0, "second call repeats: %s", logText);
	freeRateReports(&atom);
	run_logged(rec);
	CHECK(count("CE tables") == 1, "new run: %s", logText);

	/* --- Tables that cover the atmosphere: no message -- -------- */

	freeRateReports(&atom);
	run_logged("TEMP 3 500 3e3 3e5\nCE 0 2 1e-14 2e-15 1e-16\n");
	CHECK(strlen(logText) == 0, "covering table: message %s", logText);

	/* --- The other types hold the value at the end of the table -- */

	{
		double      Cval[4] = {1.0e-20, 2.0e-20, 3.0e-20, 4.0e-20}, dE = E_lev[1] - E_lev[0];
		const char *kw[] = {"CI", "CP", "CH", "CH0", "CH+"};
		int         ij, ji;

		for (n = 0; n < 5; n++) {
			freeRateReports(&atom);
			snprintf(rec, sizeof(rec), "TEMP 4 5e3 1e4 2e4 4e4\n%s 0 1 1e-20 2e-20 3e-20 4e-20\n",
					 kw[n]);
			run_logged(rec);
			CHECK(count("tables are used outside") == 1 && count(kw[n]) >= 1,
				  "%s messages: %s", kw[n], logText);
			ij = 0 * NLEVEL + 1;
			ji = 1 * NLEVEL + 0;
			for (k = 0; k < NT; k++) {
				double Tk = rh->atmos.T[k], Ck = Cval[(Tk > 4.0e4) ? 3 : 0], want, got;

				if (Tk >= 5.0e3 && Tk <= 4.0e4) continue;
				if (!strcmp(kw[n], "CI")) {
					want = Ck * rh->atmos.ne[k] * exp(-dE / (KBOLTZMANN * Tk)) * sqrt(Tk);
					got  = atom.C[ji][k];
				} else if (!strcmp(kw[n], "CP")) {
					want = hydrogen.n[1][k] * Ck;
					got  = atom.C[ij][k];
				} else if (!strcmp(kw[n], "CH")) {
					want = hydrogen.n[0][k] * Ck;
					got  = atom.C[ji][k];
				} else if (!strcmp(kw[n], "CH0")) {
					want = hydrogen.n[0][k] * Ck;
					got  = atom.C[ij][k];
				} else {
					want = hydrogen.n[1][k] * Ck;
					got  = atom.C[ji][k];
				}
				CHECK_CLOSE(got, want, 1.0e-9, 0.0, "%s held at T = %g K", kw[n], Tk);
			}
		}
	}
	return check_summary("test_collision_tables");
}
/* ------- end ----------------------- test_collision_tables.c ----- */
