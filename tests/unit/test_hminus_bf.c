/* ------- file: -------------------------- test_hminus_bf.c --------

	   H-minus bound-free cross-section (hydrogen.c: Hminus_bf_sigma,
	   Hminus_bf_readTable; N-35).

	   -- The table of the distribution (Atoms/Hminus_bf.ecsv) reads,
		  and the cross-section is continuous where the fits take over
		  at both ends of the table.
	   -- The fits give the values of docs/opacity/hminus_bf_table.py
		  (computed there with the constants of hydrogen.c).
	   -- Inside the table the cross-section interpolates linearly.
	   -- Zero at and beyond threshold and for lambda <= 0.
	   -- Malformed tables are refused with a message, and the table
		  read before stays in use.
	   --                                              -------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "background.h"
#include "check.h"

#define MB 1.0E-22

/* ------- begin -------------------------- write_file.c ------------ */

static const char *write_file(const char *name, const char *text)
{
	static char path[512];
	FILE       *fp;

	snprintf(path, sizeof(path), "%s", name);
	fp = fopen(path, "w");
	fputs(text, fp);
	fclose(fp);
	return path;
}
/* ------- end ---------------------------- write_file.c ------------ */

static const char header[] =
	"# %ECSV 1.0\n"
	"# ---\n"
	"# datatype:\n"
	"# - {name: lambda, unit: nm, datatype: float64}\n"
	"# - {name: sigma, unit: Mbarn, datatype: float64}\n"
	"# schema: astropy-2.0\n";

/* ------- begin -------------------------- refuse.c ---------------- */

static void refuse(const char *what, const char *text)
{
	char errmsg[512] = "";
	char file[64];

	snprintf(file, sizeof(file), "hmbf_%s.ecsv", what);
	CHECK(Hminus_bf_readTable(write_file(file, text), errmsg,
							  sizeof(errmsg)) == -1,
		  "%s: table accepted", what);
	CHECK(strlen(errmsg) > 0, "%s: no message", what);
	remove(file);
}
/* ------- end ---------------------------- refuse.c ---------------- */

int main(void)
{
	char   errmsg[512] = "", text[2048];
	double short_ratio, long_ratio, lo, hi, eps = 1.0E-9;

	/* --- The table of the distribution --          -------------- */

	CHECK(Hminus_bf_readTable(RH_SOURCE_DIR "/Atoms/Hminus_bf.ecsv", errmsg,
							  sizeof(errmsg)) == 0,
		  "reading Atoms/Hminus_bf.ecsv: %s", errmsg);

	Hminus_bf_joins(&short_ratio, &long_ratio);
	CHECK_CLOSE(short_ratio, 1.0, 1.0E-5, 0.0, "%s", "fit/table at 80 nm");
	CHECK_CLOSE(long_ratio, 1.0, 1.0E-5, 0.0, "%s", "fit/table at 125 nm");
	lo = Hminus_bf_sigma(80.0 - eps);
	hi = Hminus_bf_sigma(80.0 + eps);
	CHECK_CLOSE(lo, hi, 1.0E-5, 0.0, "%s", "continuity at 80 nm");
	lo = Hminus_bf_sigma(125.0 - eps);
	hi = Hminus_bf_sigma(125.0 + eps);
	CHECK_CLOSE(lo, hi, 1.0E-5, 0.0, "%s", "continuity at 125 nm");

	/* --- Fits, against the generator --             -------------- */

	CHECK_CLOSE(Hminus_bf_sigma(20.0), 0.08714242014623798 * MB, 1.0E-9, 0.0,
				"%s", "short fit at 20 nm");
	CHECK_CLOSE(Hminus_bf_sigma(50.0), 1.2959942545909797 * MB, 1.0E-9, 0.0,
				"%s", "short fit at 50 nm");
	CHECK_CLOSE(Hminus_bf_sigma(79.0), 4.314250362202004 * MB, 1.0E-9, 0.0,
				"%s", "short fit at 79 nm");
	CHECK_CLOSE(Hminus_bf_sigma(126.0), 5.560354133740078 * MB, 1.0E-9, 0.0,
				"%s", "long fit at 126 nm");
	CHECK_CLOSE(Hminus_bf_sigma(300.0), 16.27975573184688 * MB, 1.0E-9, 0.0,
				"%s", "long fit at 300 nm");
	CHECK_CLOSE(Hminus_bf_sigma(1000.0), 37.34888167011167 * MB, 1.0E-9, 0.0,
				"%s", "long fit at 1000 nm");
	CHECK_CLOSE(Hminus_bf_sigma(1640.0), 0.036766059185120215 * MB, 1.0E-8,
				0.0, "%s", "long fit at 1640 nm");

	/* --- Table knots of the distribution (shape resonance) -- -- */

	CHECK_CLOSE(Hminus_bf_sigma(112.95), 76.3388 * MB, 1.0E-12, 0.0,
				"%s", "knot at 112.95 nm");

	/* --- Threshold and lambda <= 0 --              -------------- */

	CHECK(Hminus_bf_sigma(1643.907988) == 0.0, "sigma at threshold");
	CHECK(Hminus_bf_sigma(2000.0) == 0.0, "sigma beyond threshold");
	CHECK(Hminus_bf_sigma(0.0) == 0.0, "sigma at lambda = 0");
	CHECK(Hminus_bf_sigma(-1.0) == 0.0, "sigma at lambda < 0");

	/* --- A synthetic table: linear interpolation --  -------------- */

	snprintf(text, sizeof(text), "%slambda sigma\n70.0 2.0\n100.0 8.0\n130.0 5.0\n",
			 header);
	CHECK(Hminus_bf_readTable(write_file("hmbf_ok.ecsv", text), errmsg,
							  sizeof(errmsg)) == 0,
		  "synthetic table: %s", errmsg);
	remove("hmbf_ok.ecsv");
	CHECK_CLOSE(Hminus_bf_sigma(85.0), 5.0 * MB, 1.0E-12, 0.0, "%s", "interpolation at 85 nm");
	CHECK_CLOSE(Hminus_bf_sigma(115.0), 6.5 * MB, 1.0E-12, 0.0, "%s", "interpolation at 115 nm");
	CHECK_CLOSE(Hminus_bf_sigma(100.0), 8.0 * MB, 1.0E-12, 0.0, "%s", "value at a knot");
	CHECK_CLOSE(Hminus_bf_sigma(50.0), 1.2959942545909797 * MB, 1.0E-9, 0.0,
				"%s", "short fit below a synthetic table");
	CHECK_CLOSE(Hminus_bf_sigma(300.0), 16.27975573184688 * MB, 1.0E-9, 0.0,
				"%s", "long fit above a synthetic table");

	/* --- Refused tables; the synthetic table stays in use -- --- */

	CHECK(Hminus_bf_readTable("hmbf_no_such_file.ecsv", errmsg,
							  sizeof(errmsg)) == -1,
		  "missing file accepted");
	refuse("not_ecsv", "lambda sigma\n80.0 1.0\n125.0 2.0\n");
	snprintf(text, sizeof(text),
			 "# %%ECSV 1.0\n# ---\n# datatype:\n"
			 "# - {name: lambda, unit: Angstrom, datatype: float64}\n"
			 "# - {name: sigma, unit: Mbarn, datatype: float64}\n"
			 "lambda sigma\n800.0 1.0\n1250.0 2.0\n");
	refuse("unit", text);
	snprintf(text, sizeof(text),
			 "# %%ECSV 1.0\n# ---\n# delimiter: ','\n# datatype:\n"
			 "# - {name: lambda, unit: nm, datatype: float64}\n"
			 "# - {name: sigma, unit: Mbarn, datatype: float64}\n"
			 "lambda,sigma\n80.0,1.0\n125.0,2.0\n");
	refuse("comma", text);
	snprintf(text, sizeof(text), "%ssigma lambda\n80.0 1.0\n125.0 2.0\n", header);
	refuse("order", text);
	snprintf(text, sizeof(text), "%slambda sigma\n80.0 1.0\n100.0 2.0\n90.0 3.0\n125.0 2.0\n",
			 header);
	refuse("monotonic", text);
	snprintf(text, sizeof(text), "%slambda sigma\n80.0 1.0\n125.0 -2.0\n", header);
	refuse("negative", text);
	snprintf(text, sizeof(text), "%slambda sigma\n90.0 1.0\n125.0 2.0\n", header);
	refuse("short_end", text);
	snprintf(text, sizeof(text), "%slambda sigma\n80.0 1.0\n115.0 2.0\n", header);
	refuse("long_end", text);
	snprintf(text, sizeof(text), "%slambda sigma\n80.0 1.0 3.0\n125.0 2.0\n", header);
	refuse("columns", text);

	CHECK_CLOSE(Hminus_bf_sigma(85.0), 5.0 * MB, 1.0E-12, 0.0,
				"%s", "table after refused reads");

	Hminus_bf_freeTable();
	return check_summary("test_hminus_bf");
}
/* ------- end ---------------------------- test_hminus_bf.c -------- */
