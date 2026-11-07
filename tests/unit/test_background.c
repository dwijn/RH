/* ------- file: -------------------------- test_background.c ------

	   The background opacity file (readj.c: writeBackground,
	   readBackground; brs_xdr.c: writeBRS, readBRS), single precision
	   since M-1/M-7:

	   -- a round trip of unpolarized and polarized (magneto-optical)
		  blocks through the float file: every opacity comes back as its
		  float rounding, every emissivity (stored times 2^BACKGR_ETA_EXP)
		  as the float rounding of the scaled value, also far below the
		  smallest float (1e-60); the file holds 4 bytes per value;
	   -- a file in the old format (double, unscaled) reads back exactly;
	   -- an emissivity beyond the range of the file is a fatal error;
	   -- brs.out round trip: the format ints are written and read back;
		  a brs.out without them (older versions) means double, scale 1;
		  an unknown format is a fatal error.
	   --                                              -------------- */

#include <fcntl.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "inputs.h"
#include "error.h"
#include "check.h"

extern Atmosphere  atmos;
extern Spectrum    spectrum;
extern InputData   input;
extern CommandLine commandline;

void writeBRS(void);
void readBRS(void);
int  backgroundRecords(int nspect);
int  writeBackground(int nspect, int mu, bool_t to_obs, double *chi_c,
					 double *eta_c, double *sca_c, double *chip_c);
void readBackground(int nspect, int mu, bool_t to_obs);

#define NSPACE 7
#define NSPECT 3

/* --- Value i of a record: log-uniform between lo and hi, sign
	   alternating if signed --                         -------------- */

static double value(int i, int rec, double lo, double hi, bool_t sgn)
{
	double f = fmod(0.618033988749895 * (i + 13 * rec + 1), 1.0);
	double v = lo * pow(hi / lo, f);

	return (sgn && (i + rec) % 2) ? -v : v;
}

static double chi[NSPECT][4 * NSPACE], eta[NSPECT][4 * NSPACE],
	sca[NSPECT][NSPACE], chip[NSPECT][3 * NSPACE];

static void setup(void)
{
	int nspect, k, rec = 0;

	atmos.Nspace = NSPACE;
	atmos.Nrays  = 1;
	atmos.moving = atmos.Stokes = FALSE;
	strcpy(atmos.ID, "test_background");
	spectrum.Nspect       = NSPECT;
	input.magneto_optical = TRUE;
	input.StokesMode      = FULL_STOKES;

	atmos.backgrflags                = calloc(NSPECT, sizeof(flags));
	atmos.backgrrecno                = calloc(NSPECT, sizeof(int));
	spectrum.as                      = calloc(NSPECT, sizeof(ActiveSet));
	atmos.backgrflags[1].ispolarized = TRUE;
	for (nspect = 0; nspect < NSPECT; nspect++) {
		atmos.backgrrecno[nspect] = rec;
		rec += backgroundRecords(nspect);
		spectrum.as[nspect].chi_c  = calloc(4 * NSPACE, sizeof(double));
		spectrum.as[nspect].eta_c  = calloc(4 * NSPACE, sizeof(double));
		spectrum.as[nspect].sca_c  = calloc(NSPACE, sizeof(double));
		spectrum.as[nspect].chip_c = calloc(3 * NSPACE, sizeof(double));
		for (k = 0; k < 4 * NSPACE; k++) {
			chi[nspect][k] = value(k, nspect, 1.0e-20, 1.0e5, k >= NSPACE);
			eta[nspect][k] = value(k, nspect + 7, 1.0e-60, 1.0e2, k >= NSPACE);
		}
		for (k = 0; k < NSPACE; k++) sca[nspect][k] = value(k, nspect + 3, 1.0e-25, 1.0, FALSE);
		for (k = 0; k < 3 * NSPACE; k++) chip[nspect][k] = value(k, nspect + 5, 1.0e-20, 1.0e3, TRUE);
	}
}

static int open_file(const char *name)
{
	return open(name, O_RDWR | O_CREAT | O_TRUNC, 0644);
}

static void write_all(void)
{
	for (int nspect = 0; nspect < NSPECT; nspect++)
		writeBackground(nspect, 0, FALSE, chi[nspect], eta[nspect], sca[nspect],
						chip[nspect]);
}

/* --- What the file should give back: the value itself (double file)
	   or the float rounding of value * scale, divided by scale -- --- */

static double expect(double v, int exp2)
{
	if (atmos.backgrbytes == sizeof(double)) return v;
	return (double) (float) ldexp(v, exp2) / ldexp(1.0, exp2);
}

static void check_read_back(const char *what)
{
	int nspect, k, Nrec, e = atmos.backgretaexp;

	for (nspect = 0; nspect < NSPECT; nspect++) {
		ActiveSet *as = &spectrum.as[nspect];

		memset(as->chi_c, 0, 4 * NSPACE * sizeof(double));
		memset(as->eta_c, 0, 4 * NSPACE * sizeof(double));
		memset(as->sca_c, 0, NSPACE * sizeof(double));
		memset(as->chip_c, 0, 3 * NSPACE * sizeof(double));
		readBackground(nspect, 0, FALSE);

		Nrec = atmos.backgrflags[nspect].ispolarized ? 4 : 1;
		for (k = 0; k < Nrec * NSPACE; k++) {
			CHECK(as->chi_c[k] == expect(chi[nspect][k], 0),
				  "%s: chi_c[%d] of wavelength %d: %.9g, want %.9g", what, k, nspect,
				  as->chi_c[k], expect(chi[nspect][k], 0));
			CHECK(as->eta_c[k] == expect(eta[nspect][k], e),
				  "%s: eta_c[%d] of wavelength %d: %.9g, want %.9g", what, k, nspect,
				  as->eta_c[k], expect(eta[nspect][k], e));
			CHECK_CLOSE(as->eta_c[k], eta[nspect][k], 0x1p-24, 0.0,
						"%s: eta_c[%d], wavelength %d", what, k, nspect);
		}
		for (k = 0; k < NSPACE; k++)
			CHECK(as->sca_c[k] == expect(sca[nspect][k], 0),
				  "%s: sca_c[%d] of wavelength %d", what, k, nspect);
		if (Nrec == 4)
			for (k = 0; k < 3 * NSPACE; k++)
				CHECK(as->chip_c[k] == expect(chip[nspect][k], 0),
					  "%s: chip_c[%d] of wavelength %d", what, k, nspect);
	}
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
		if (!freopen("/dev/null", "w", stdout) || !freopen("/dev/null", "w", stderr)) _exit(1);
		commandline.logfile = stderr;
		f();
		_exit(0);
	}
	if (pid < 0 || waitpid(pid, &status, 0) != pid) return FALSE;
	return WIFEXITED(status) && WEXITSTATUS(status) == RH_EXIT_FATAL;
}

static void write_huge_eta(void)
{
	eta[0][2] = 1.0e10; /* times 2^100: beyond FLT_MAX */
	writeBackground(0, 0, FALSE, chi[0], eta[0], sca[0], chip[0]);
}

static void read_bad_brs(void)
{
	readBRS();
}

int main(void)
{
	char        dir[] = "/tmp/rh_test_background.XXXXXX";
	struct stat st;
	int         Nrec = 0, nspect;
	FILE       *fp;
	long        size;

	commandline.logfile = stderr;
	if (!mkdtemp(dir) || chdir(dir)) {
		perror("test_background: temporary directory");
		return 1;
	}
	setup();
	for (nspect = 0; nspect < NSPECT; nspect++) Nrec += backgroundRecords(nspect);

	/* --- Single-precision file, emissivities scaled -- ------------ */

	atmos.backgrbytes   = sizeof(float);
	atmos.backgretaexp  = BACKGR_ETA_EXP;
	atmos.fd_background = open_file("background.dat");
	write_all();
	check_read_back("float file");
	fstat(atmos.fd_background, &st);
	CHECK(st.st_size == (off_t) Nrec * NSPACE * 4, "float file has %lld bytes, want %d",
		  (long long) st.st_size, Nrec * NSPACE * 4);

	/* --- The emissivities below the smallest float come back nonzero,
		   to float precision (checked above) --         -------------- */

	{
		int below = 0, zero = 0;

		for (nspect = 0; nspect < NSPECT; nspect++) {
			int Nread = (atmos.backgrflags[nspect].ispolarized ? 4 : 1) * NSPACE;

			for (int k = 0; k < Nread; k++)
				if (fabs(eta[nspect][k]) < FLT_MIN) {
					below++;
					if (spectrum.as[nspect].eta_c[k] == 0.0) zero++;
				}
		}
		CHECK(below >= 5 && zero == 0, "%d emissivities below FLT_MIN, %d read back as 0",
			  below, zero);
	}

	/* --- Old (double) file --                           -------------- */

	close(atmos.fd_background);
	atmos.backgrbytes   = sizeof(double);
	atmos.backgretaexp  = 0;
	atmos.fd_background = open_file("background_double.dat");
	write_all();
	check_read_back("double file");
	fstat(atmos.fd_background, &st);
	CHECK(st.st_size == (off_t) Nrec * NSPACE * 8, "double file has %lld bytes",
		  (long long) st.st_size);

	/* --- Emissivity beyond the range of the float file -- -------- */

	close(atmos.fd_background);
	atmos.backgrbytes   = sizeof(float);
	atmos.backgretaexp  = BACKGR_ETA_EXP;
	atmos.fd_background = open_file("background.dat");
	CHECK(is_fatal(write_huge_eta), "%s", "an emissivity of 1e10 was written without error");

	/* --- brs.out: the format ints round trip --        -------------- */

	writeBRS();
	atmos.backgrbytes = atmos.backgretaexp = -1;
	readBRS();
	CHECK(atmos.backgrbytes == sizeof(float) && atmos.backgretaexp == BACKGR_ETA_EXP,
		  "brs.out read back as %d bytes, scale 2^%d", atmos.backgrbytes,
		  atmos.backgretaexp);
	CHECK(atmos.backgrrecno[1] == backgroundRecords(0), "record of wavelength 1: %d",
		  atmos.backgrrecno[1]);

	/* --- A brs.out of an older version, without the ints: double -- */

	fp = fopen("brs.out", "r+");
	fseek(fp, 0, SEEK_END);
	size = ftell(fp);
	fclose(fp);
	CHECK(truncate("brs.out", size - 8) == 0, "%s", "truncate brs.out");
	atmos.backgrbytes = atmos.backgretaexp = -1;
	readBRS();
	CHECK(atmos.backgrbytes == sizeof(double) && atmos.backgretaexp == 0,
		  "old brs.out read as %d bytes, scale 2^%d", atmos.backgrbytes,
		  atmos.backgretaexp);

	/* --- An unknown format is fatal --                  -------------- */

	atmos.backgrbytes  = 5;
	atmos.backgretaexp = 0;
	writeBRS();
	CHECK(is_fatal(read_bad_brs), "%s", "brs.out with 5 bytes per value was accepted");

	close(atmos.fd_background);
	unlink("background.dat");
	unlink("background_double.dat");
	unlink("brs.out");
	if (chdir("/") == 0) rmdir(dir);

	return check_summary("test_background");
}
