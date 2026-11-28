/* ------- file: -------------------------- test_gaunt_ff.c ---------

	   Thermally averaged free-free Gaunt factor (hydrogen.c: Gaunt_ff,
	   Gaunt_ff_readTable; N-50).

	   -- The table of the distribution (Atoms/gauntff.dat, van Hoof
		  et al. 2014) reads, and Gaunt_ff gives its values at the
		  nodes, also next to its edges.
	   -- Far below the frequency of kT and for gamma^2 >> 1, the
		  classical low-frequency limit (Kramers normalisation):
		  g = sqrt(3)/pi [ln((2kT)^(3/2) / (pi e^2 m^(1/2) nu))
		  - 5 gamma_E / 2] (cgs), between the nodes.
	   -- In the mm, where the expansion RH used before gave 1, the
		  values of the table (2.2-3.5).
	   -- A synthetic table whose ln g is cubic in log gamma^2 and log u
		  is reproduced exactly between the nodes, also with the
		  one-sided stencils at its edges, and the charge enters
		  through gamma^2 = Z^2 Ry / kT.
	   -- Malformed tables are refused with a message, and the table
		  read before stays in use.
	   --                                              -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"
#include "background.h"
#include "error.h"
#include "check.h"
#include "state.h"

/* ------- begin -------------------------- at.c -------------------- */

/* --- Gaunt_ff with its status checked --         -------------- */

static double gaunt(double lambda, int Z, double T)
{
	double    g = 0.0;
	rh_status s = Gaunt_ff(rh, lambda, Z, T, &g);

	CHECK(s == RH_OK, "Gaunt_ff(%g nm, Z = %d, %g K): status %d", lambda, Z, T, s);
	return g;
}

static double at(double lg2, double lu, int Z)
{
	/* --- Gaunt_ff at log gamma^2 and log u. --      -------------- */

	double T      = SQ(Z) * E_RYDBERG / (KBOLTZMANN * pow(10.0, lg2));
	double lambda = HPLANCK * CLIGHT / (pow(10.0, lu) * KBOLTZMANN * T) / NM_TO_M;

	return gaunt(lambda, Z, T);
}
/* ------- end ---------------------------- at.c -------------------- */

/* ------- begin -------------------------- write_file.c ------------ */

static const char *write_file(const char *name, const char *text)
{
	FILE *fp = fopen(name, "w");

	fputs(text, fp);
	fclose(fp);
	return name;
}
/* ------- end ---------------------------- write_file.c ------------ */

/* ------- begin -------------------------- refuse.c ---------------- */

static void refuse(const char *what, const char *text)
{
	char errmsg[512] = "";
	char file[64];

	snprintf(file, sizeof(file), "gff_%s.dat", what);
	CHECK(Gaunt_ff_readTable(write_file(file, text), errmsg,
							 sizeof(errmsg)) == -1,
		  "%s: table accepted", what);
	CHECK(strlen(errmsg) > 0, "%s: no message", what);
	remove(file);
}
/* ------- end ---------------------------- refuse.c ---------------- */

/* --- Synthetic table: 6 x 7 nodes, log gamma^2 = -1 (0.5) 1.5,
	   log u = -2 (0.5) 1, ln g cubic in both --       -------------- */

#define SYN_G0   -1.0
#define SYN_U0   -2.0
#define SYN_STEP 0.5
#define SYN_NG   6
#define SYN_NU   7

static double syn_lng(double a, double b)
{
	return 0.3 + 0.2 * a - 0.1 * b + 0.05 * a * a * a - 0.04 * a * b * b +
		   0.03 * b * b * b + 0.01 * a * a * b * b * b;
}

static void syn_table(char *text, size_t len)
{
	int    i, j;
	size_t n;

	n = snprintf(text, len,
				 "# synthetic\n 20140210 # magic number\n %d %d\n %g\n %g\n %g\n",
				 SYN_NG, SYN_NU, SYN_G0, SYN_U0, SYN_STEP);
	for (j = 0; j < SYN_NU; j++) {
		for (i = 0; i < SYN_NG; i++)
			n += snprintf(text + n, len - n, " %.17e",
						  exp(syn_lng(SYN_G0 + i * SYN_STEP, SYN_U0 + j * SYN_STEP)));
		n += snprintf(text + n, len - n, "\n");
	}
}

int main(void)
{
	char   errmsg[512] = "", text[8192];
	double k = 1.380658E-16, e = 4.80320471E-10, m = 9.1093897E-28,
		   h = 6.6260755E-27, gammaE = 0.5772156649015329;
	double lg2, lu, T, nu, g_cl;
	int    n;

	/* --- The table of the distribution, at its nodes -- --------- */

	CHECK(Gaunt_ff_readTable(RH_SOURCE_DIR "/Atoms/gauntff.dat", errmsg,
							 sizeof(errmsg)) == 0,
		  "reading Atoms/gauntff.dat: %s", errmsg);

	CHECK_CLOSE(at(2.0, -4.0, 1), 3.83181090, 1.0E-9, 0.0, "%s", "node (2, -4)");
	CHECK_CLOSE(at(1.0, 0.0, 1), 1.19866181, 1.0E-9, 0.0, "%s", "node (1, 0)");
	CHECK_CLOSE(at(-1.0, 2.0, 1), 0.228112124, 1.0E-9, 0.0, "%s", "node (-1, 2)");
	CHECK_CLOSE(at(0.4, -2.6, 1), 3.11685474, 1.0E-9, 0.0, "%s",
				"node (0.4, -2.6)");
	CHECK_CLOSE(at(-5.8, -15.8, 1), 20.5036480, 1.0E-9, 0.0, "%s",
				"node next to the low corner");
	CHECK_CLOSE(at(9.8, 12.8, 1), 0.198702432, 1.0E-9, 0.0, "%s",
				"node next to the high corner");

	/* --- Classical low-frequency limit, between the nodes -- ---- */

	for (n = 0; n < 4; n++) {
		lg2  = 1.3 + 0.37 * n;
		lu   = -6.1 - 1.13 * n;
		T    = E_RYDBERG / (KBOLTZMANN * pow(10.0, lg2));
		nu   = pow(10.0, lu) * k * T / h;
		g_cl = sqrt(3.0) / PI *
			   (log(pow(2.0 * k * T, 1.5) / (PI * SQ(e) * sqrt(m) * nu)) -
				2.5 * gammaE);
		CHECK_CLOSE(at(lg2, lu, 1), g_cl, 2.0E-3, 0.0,
					"classical limit at log gamma^2 = %.2f, log u = %.2f", lg2, lu);
	}

	/* --- The mm (the expansion RH used gave 1 there) -- ---------- */

	CHECK_CLOSE(gaunt(1.0E6, 1, 3000.0), 2.177, 2.0E-3, 0.0, "%s",
				"1 mm, 3000 K");
	CHECK_CLOSE(gaunt(1.0E6, 1, 1.0E4), 2.958, 2.0E-3, 0.0, "%s",
				"1 mm, 10^4 K");
	CHECK_CLOSE(gaunt(3.0E6, 1, 1.0E4), 3.497, 2.0E-3, 0.0, "%s",
				"3 mm, 10^4 K");

	/* --- Outside the table (log u = 14, the table ends at 13):
		   RH_E_NUMERIC, recorded -- ---------------------------- */
	{
		double g = 0.0, T = 1.0E4,
			   lambda = HPLANCK * CLIGHT / (1.0E14 * KBOLTZMANN * T) / NM_TO_M;

		CHECK(Gaunt_ff(rh, lambda, 1, T, &g) == RH_E_NUMERIC && Failed(rh) &&
				  rh->error.code == RH_E_NUMERIC,
			  "%s", "log u = 14: not RH_E_NUMERIC");
		ClearError(rh);
	}

	/* --- Synthetic table: bicubic in ln g is exact -- ------------ */

	syn_table(text, sizeof(text));
	CHECK(Gaunt_ff_readTable(write_file("gff_syn.dat", text), errmsg,
							 sizeof(errmsg)) == 0,
		  "synthetic table: %s", errmsg);
	remove("gff_syn.dat");
	for (lg2 = -0.97; lg2 <= 1.47; lg2 += 0.11) {
		for (lu = -1.98; lu <= 0.98; lu += 0.13) {
			CHECK_CLOSE(log(at(lg2, lu, 1)), syn_lng(lg2, lu), 1.0E-12, 1.0E-12,
						"synthetic ln g at (%.2f, %.2f)", lg2, lu);
			CHECK_CLOSE(log(at(lg2, lu, 2)), syn_lng(lg2, lu), 1.0E-12, 1.0E-12,
						"synthetic ln g at (%.2f, %.2f), Z = 2", lg2, lu);
		}
	}

	/* --- Malformed tables, refused; the table read before stays -- */

	refuse("magic", "# x\n 20150101\n 4 4\n 0\n 0\n 0.2\n 1 1 1 1\n");
	refuse("small", "20140210\n 3 4\n 0\n 0\n 0.2\n 1 1 1\n");
	refuse("step", "20140210\n 4 4\n 0\n 0\n 0\n 1 1 1 1\n");
	refuse("zero", "20140210\n 4 4\n 0\n 0\n 0.2\n 1 1 0 1\n");
	refuse("short", "20140210\n 4 4\n 0\n 0\n 0.2\n 1 1 1 1\n 1 1 1 1\n");
	refuse("word", "20140210\n 4 4\n 0\n 0\n 0.2\n 1 1 x 1\n");
	refuse("header", "20140210\n 4 4\n 0\n");
	CHECK(Gaunt_ff_readTable("no_such_gauntff.dat", errmsg, sizeof(errmsg)) == -1,
		  "%s", "missing file accepted");
	CHECK_CLOSE(log(at(0.3, -0.7, 1)), syn_lng(0.3, -0.7), 1.0E-12, 1.0E-12,
				"%s", "synthetic table still in use after the refusals");

	Gaunt_ff_freeTable();
	return check_summary("test_gaunt_ff");
}
