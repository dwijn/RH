/* ------- file: -------------------------- readj.c -----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Feb 20 17:12:25 2024 --

	   --------------------------                      ----------RH-- */

/* --- Routines to read and write angle-averaged mean intensity and
	   background opacities and emissivities. --       -------------- */

/* --- Note: pread and pwrite are only able to deal with chunks
			 that are less than 2^31 - 1 bytes. The routines
			 pread_rh and pwrite_rh read/write chunks of PWRITE_SLICE
			 sequentially to overcome this. --         -------------- */

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "constant.h"
#include "inputs.h"
#include "error.h"
#include "prdh.h"

#define PWRITE_SLICE 1073741824

/* --- Function prototypes --                          -------------- */

size_t pwrite_rh(int fd, const void *buf, size_t count, off_t offset);
size_t pread_rh(int fd, const void *buf, size_t count, off_t offset);

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum   spectrum;
extern InputData  input;
extern char       messageStr[];

/* ------- begin -------------------------- readJlambda.c ----------- */

void readJlambda(int nspect, double *J)
{
	const char routineName[] = "readJlambda";

	bool_t result = TRUE;
	size_t recordsize;
	off_t  offset;

	recordsize = atmos.Nspace * sizeof(double);
	offset     = recordsize * nspect;

	result &= (pread_rh(spectrum.fd_J, J, recordsize, offset) == recordsize);

	if (!result) {
		sprintf(messageStr,
				"Error reading file: offset = %lld, recordsize = %zu (%s)",
				(long long) offset, recordsize, strerror(errno));
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
}
/* ------- end ---------------------------- readJlambda.c ----------- */

/* ------- begin -------------------------- writeJlambda.c ---------- */

void writeJlambda(int nspect, double *J)
{
	const char routineName[] = "writeJlambda";

	bool_t result = TRUE;
	size_t recordsize;
	off_t  offset;

	recordsize = atmos.Nspace * sizeof(double);
	offset     = recordsize * nspect;

	result &= (pwrite_rh(spectrum.fd_J, J, recordsize, offset) == recordsize);

	if (!result) {
		sprintf(messageStr,
				"Error writing file: offset = %lld, recordsize = %zu (%s)",
				(long long) offset, recordsize, strerror(errno));
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
}
/* ------- end ---------------------------- writeJlambda.c ---------- */

/* ------- begin -------------------------- readJgas.c -------------- */

void readJgas(double **Jgas)
{
	const char routineName[] = "readJgas";

	size_t recordsize;
	FILE  *fp;

	recordsize = (size_t) atmos.Nspace * spectrum.Nspect;
	if ((fp = fopen("Jgas.dat", "r")) == NULL) {
		sprintf(messageStr, "Cannot open Jgas.dat for reading (%s)",
				strerror(errno));
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	if (fread(&(Jgas[0][0]), sizeof(double), recordsize, fp) != recordsize) {
		sprintf(messageStr, "Error reading Jgas.dat (%zu values expected)",
				recordsize);
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}

	fclose(fp);
}
/* ------- end ---------------------------- readJgas.c ------------- */

/* ------- begin -------------------------- writeJgas.c ------------ */

void writeJgas(double **Jgas)
{
	const char routineName[] = "writeJgas";

	size_t recordsize;
	FILE  *fp;

	recordsize = (size_t) atmos.Nspace * spectrum.Nspect;
	if ((fp = fopen("Jgas.dat", "w")) == NULL) {
		sprintf(messageStr, "Cannot open Jgas.dat for writing (%s)",
				strerror(errno));
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	if (fwrite(&(Jgas[0][0]), sizeof(double), recordsize, fp) != recordsize) {
		sprintf(messageStr, "Error writing Jgas.dat (%zu values expected)",
				recordsize);
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}

	fclose(fp);
}
/* ------- end ------------------------------- writeJgas.c ---------- */

/* ------- begin -------------------------- readJ20lambda.c --------- */

void readJ20lambda(int nspect, double *J20)
{
	const char routineName[] = "readJ20lambda";

	bool_t result = TRUE;
	size_t recordsize;
	off_t  offset;

	recordsize = atmos.Nspace * sizeof(double);
	offset     = recordsize * nspect;

	result &= (pread_rh(spectrum.fd_J20, J20,
						recordsize, offset) == recordsize);

	if (!result) {
		sprintf(messageStr,
				"Error reading file: offset = %lld, recordsize = %zu (%s)",
				(long long) offset, recordsize, strerror(errno));
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
}
/* ------- end ---------------------------- readJ20lambda.c --------- */

/* ------- begin -------------------------- writeJ20lambda.c -------- */

void writeJ20lambda(int nspect, double *J20)
{
	const char routineName[] = "writeJ20lambda";

	bool_t result = TRUE;
	size_t recordsize;
	off_t  offset;

	recordsize = atmos.Nspace * sizeof(double);
	offset     = recordsize * nspect;

	result &= (pwrite_rh(spectrum.fd_J20, J20,
						 recordsize, offset) == recordsize);

	if (!result) {
		sprintf(messageStr,
				"Error writing file: offset = %lld, recordsize = %zu (%s)",
				(long long) offset, recordsize, strerror(errno));
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
}
/* ------- end ---------------------------- writeJ20lambda.c -------- */

/* ------- begin -------------------------- readImu.c --------------- */

void readImu(int nspect, int mu, bool_t to_obs, double *I)
{
	const char routineName[] = "readImu";

	bool_t result = TRUE;
	int    index;
	size_t recordsize;
	off_t  offset;

	recordsize = atmos.Nspace * sizeof(double);

	index  = spectrum.PRDindex[nspect];
	offset = 2 * (index * atmos.Nrays + mu) * recordsize;
	if (to_obs) offset += recordsize;

	result &= (pread_rh(spectrum.fd_Imu, I,
						recordsize, offset) == recordsize);

	if (!result) {
		sprintf(messageStr,
				"Error reading file: offset = %lld, recordsize = %zu (%s)",
				(long long) offset, recordsize, strerror(errno));
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
}
/* ------- end ---------------------------- readImu.c --------------- */

/* ------- begin -------------------------- writeImu.c -------------- */

void writeImu(int nspect, int mu, bool_t to_obs, double *I)
{
	const char routineName[] = "writeImu";

	bool_t result = TRUE;
	int    index;
	size_t recordsize;
	off_t  offset;

	recordsize = atmos.Nspace * sizeof(double);

	index  = spectrum.PRDindex[nspect];
	offset = 2 * (index * atmos.Nrays + mu) * recordsize;
	if (to_obs) offset += recordsize;

	result &= (pwrite_rh(spectrum.fd_Imu, I,
						 recordsize, offset) == recordsize);

	if (!result) {
		sprintf(messageStr,
				"Error writing file: offset = %lld, recordsize = %zu (%s)",
				(long long) offset, recordsize, strerror(errno));
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
}
/* ------- end ---------------------------- writeImu.c -------------- */

/* ------- begin -------------------------- mapBackground.c --------- */

/* --- The background file read-only mapped into memory (P-18). The
	   formal solution reads the background of every wavelength, and in
	   a moving atmosphere of every ray and direction, in every sweep:
	   with pread that is a system call and a copy into a buffer
	   (allocated for each read) per record block, 6 % of the time of
	   the formal solution on a 1-D Bifrost column. The mapping takes
	   the pages from the page cache (or /dev/shm), which hold them
	   anyway, so it needs no memory of its own. Reads fall back to
	   pread if the file cannot be mapped --          -------------- */

static const char *bg_map;
static size_t      bg_maplen;

void mapBackground(void)
{
	struct stat st;
	void       *p;

	unmapBackground();
	if (fstat(atmos.fd_background, &st) != 0 || st.st_size <= 0) return;

	p = mmap(NULL, (size_t) st.st_size, PROT_READ, MAP_SHARED, atmos.fd_background, 0);
	if (p == MAP_FAILED) {
		sprintf(messageStr, "Cannot map the background file (%s); reading it with pread",
				strerror(errno));
		Error(WARNING, "mapBackground", messageStr);
		return;
	}
	bg_map    = (const char *) p;
	bg_maplen = (size_t) st.st_size;
}

void unmapBackground(void)
{
	/* --- Before the file is written again or closed -- ------------ */

	if (bg_map != NULL) {
		munmap((void *) bg_map, bg_maplen);
		bg_map    = NULL;
		bg_maplen = 0;
	}
}
/* ------- end ---------------------------- mapBackground.c --------- */

/* ------- begin -------------------------- backgroundIO.c ---------- */

static float *backgroundBuffer(void)
{
	float *buffer = NULL;

	/* --- Conversion buffer for backgroundIO: the largest part of a
		   block is 4 records (the Stokes opacity or emissivity) -- -- */

	if (atmos.backgrbytes == sizeof(float) && bg_map == NULL) {
		buffer = (float *) malloc(4 * (size_t) atmos.Nspace * sizeof(float));
		if (buffer == NULL)
			Error(ERROR_LEVEL_2, "backgroundBuffer", "Out of memory");
	}
	return buffer;
}

static bool_t backgroundIO(double *values, int Nrec, off_t recordno,
						   bool_t write_it, bool_t emissivity, float *buffer)
{
	const char routineName[] = "backgroundIO";

	size_t Nvalues = (size_t) Nrec * atmos.Nspace, n, count;
	off_t  offset;
	double scale = 1.0, v;

	/* --- Transfer Nrec records (of atmos.Nspace values each) starting
		   at record recordno of the background file. Values are double
		   in memory and atmos.backgrbytes wide on disk: float since M-1,
		   so that the file is half the size, double in files written
		   by older versions (read with OLD_BACKGROUND). The rounding to
		   float changes the background opacities by at most 6e-8
		   relative.

		   Emissivities in SI units go far below the smallest float
		   (1.2e-38) at short wavelengths and low temperatures (1e-50 in
		   FAL-C below 100 nm), and a zero emissivity can make J vanish.
		   They are stored multiplied by 2^atmos.backgretaexp (exact, a
		   power of two), which maps 1e-68 to the smallest float and
		   leaves room up to 2.7e8 before overflow. A value beyond the
		   range of float is an error rather than an infinity in the
		   file. buffer holds at least Nrec records (backgroundBuffer).
		   --                                            -------------- */

	if (emissivity && atmos.backgretaexp != 0)
		scale = ldexp(1.0, atmos.backgretaexp);

	count  = Nvalues * atmos.backgrbytes;
	offset = recordno * (off_t) atmos.Nspace * atmos.backgrbytes;

	/* --- Read from the mapped file: the same conversion as below,
		   without the copy into buffer --               -------------- */

	if (!write_it && bg_map != NULL) {
		if (offset < 0 || (size_t) offset + count > bg_maplen) return FALSE;
		if (atmos.backgrbytes == sizeof(double) && scale == 1.0) {
			memcpy(values, bg_map + offset, count);
			return TRUE;
		} else if (atmos.backgrbytes == sizeof(float)) {
			const float *src = (const float *) (bg_map + offset);

			for (n = 0; n < Nvalues; n++) values[n] = src[n] / scale;
			return TRUE;
		}
	}

	if (atmos.backgrbytes == sizeof(double) && scale == 1.0) {
		if (write_it)
			return (pwrite_rh(atmos.fd_background, values, count, offset) == count);
		else
			return (pread_rh(atmos.fd_background, values, count, offset) == count);
	} else if (atmos.backgrbytes != sizeof(float) || buffer == NULL) {
		sprintf(messageStr, "Unknown format of the background file: %d bytes"
							" per value, emissivity scale 2^%d",
				atmos.backgrbytes, atmos.backgretaexp);
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	if (write_it) {
		for (n = 0; n < Nvalues; n++) {
			v = values[n] * scale;
			if (fabs(v) > FLT_MAX && isfinite(v)) {
				sprintf(messageStr, "Background %s %.3e out of the range of"
									" the (single-precision) background file",
						(emissivity) ? "emissivity" : "opacity", values[n]);
				Error(ERROR_LEVEL_2, routineName, messageStr);
			}
			buffer[n] = (float) v;
		}
		return (pwrite_rh(atmos.fd_background, buffer, count, offset) == count);
	} else {
		if (pread_rh(atmos.fd_background, buffer, count, offset) != count)
			return FALSE;
		for (n = 0; n < Nvalues; n++) values[n] = buffer[n] / scale;
		return TRUE;
	}
}
/* ------- end ---------------------------- backgroundIO.c ---------- */

/* ------- begin -------------------------- readBackground.c -------- */

void readBackground(int nspect, int mu, bool_t to_obs)
{
	const char routineName[] = "readBackground";

	int        recordno, NrecStokes, NskipStokes;
	bool_t     result = TRUE;
	float     *buffer = backgroundBuffer();
	ActiveSet *as;

	/* --- Read background opacity, emissivity and scattering opacity
		   into ActiveSet structure as. This cannot be done sequentially
		   with relative offsets because the file may be read by
		   different threads at the same time. Therefore, we use absolute
		   offsets here and in writeBackground. --       -------------- */

	as = &spectrum.as[nspect];

	if (atmos.moving || atmos.Stokes)
		recordno = atmos.backgrrecno[2 * (nspect * atmos.Nrays + mu) + to_obs];
	else
		recordno = atmos.backgrrecno[nspect];

	/* --- Read emissivity and opacity for Q, U, and V only if we are
		   solving explicitly for all four Stokes quantities, but always
		   skip the proper amount of records --           ------------- */

	NrecStokes = 1;
	if (atmos.backgrflags[nspect].ispolarized) {
		NskipStokes = 4;
		if (input.StokesMode == FULL_STOKES) NrecStokes = 4;
	} else
		NskipStokes = 1;

	/* --- Read background opacity --                    -------------- */

	result &= backgroundIO(as->chi_c, NrecStokes, recordno, FALSE, FALSE, buffer);
	recordno += NskipStokes;

	/* --- Read off-diagonal elements propagation matrix K -- --------- */

	if (atmos.backgrflags[nspect].ispolarized && input.magneto_optical) {
		if (input.StokesMode == FULL_STOKES)
			result &= backgroundIO(as->chip_c, 3, recordno, FALSE, FALSE, buffer);
		recordno += 3;
	}
	/* --- Read background emissivity --                 -------------- */

	result &= backgroundIO(as->eta_c, NrecStokes, recordno, FALSE, TRUE, buffer);
	recordno += NskipStokes;

	/* --- Read background scattering opacity --         -------------- */

	result &= backgroundIO(as->sca_c, 1, recordno, FALSE, FALSE, buffer);
	free(buffer);

	/* --- Exit if reading is unsuccessful --            -------------- */

	if (!result) Error(ERROR_LEVEL_2, routineName, "Error reading file");
}
/* ------- end ---------------------------- readBackground.c -------- */

/* ------- begin -------------------------- backgroundRecords.c ----- */

int backgroundRecords(int nspect)
{
	int NrecStokes, Nrec;

	/* --- Number of records (of atmos.Nspace values) that a block
		   of background data for wavelength nspect occupies: opacity,
		   possibly the off-diagonal elements of the propagation matrix,
		   emissivity, and scattering opacity. This is the one place that
		   defines the layout; writeBackground and the reservation of
		   record numbers in Background() both use it. --  ------------ */

	NrecStokes = (atmos.backgrflags[nspect].ispolarized) ? 4 : 1;

	Nrec = NrecStokes;
	if (atmos.backgrflags[nspect].ispolarized && input.magneto_optical)
		Nrec += 3;
	Nrec += NrecStokes;
	Nrec += 1;

	return Nrec;
}
/* ------- end ---------------------------- backgroundRecords.c ----- */

/* ------- begin -------------------------- writeBackground.c ------- */

int writeBackground(int nspect, int mu, bool_t to_obs,
					double *chi_c, double *eta_c, double *sca_c,
					double *chip_c)
{
	const char routineName[] = "writeBackground";

	int    recordno, Nwrite, NrecStokes;
	bool_t result = TRUE;
	float *buffer = backgroundBuffer();

	/* --- Writes background opacity, emissivity, and scattering
		   opacity. Returns Nwrite, the number of records written,
		   each of atmos.Nspace values (see backgroundIO) -- --------- */

	if (atmos.moving || atmos.Stokes)
		recordno = atmos.backgrrecno[2 * (nspect * atmos.Nrays + mu) + to_obs];
	else
		recordno = atmos.backgrrecno[nspect];

	if (atmos.backgrflags[nspect].ispolarized)
		NrecStokes = 4;
	else
		NrecStokes = 1;

	Nwrite = 0;

	/* --- Write background opacity --                   -------------- */

	result &= backgroundIO(chi_c, NrecStokes, recordno + Nwrite, TRUE, FALSE, buffer);
	Nwrite += NrecStokes;

	/* --- Write off-diagonal elements of propagation matrix K -- ----- */

	if (atmos.backgrflags[nspect].ispolarized && input.magneto_optical) {
		result &= backgroundIO(chip_c, 3, recordno + Nwrite, TRUE, FALSE, buffer);
		Nwrite += 3;
	}
	/* --- Write background emissivity --                -------------- */

	result &= backgroundIO(eta_c, NrecStokes, recordno + Nwrite, TRUE, TRUE, buffer);
	Nwrite += NrecStokes;

	/* --- Write background scattering opacity --        -------------- */

	result &= backgroundIO(sca_c, 1, recordno + Nwrite, TRUE, FALSE, buffer);
	Nwrite += 1;
	free(buffer);

	if (!result)
		Error(ERROR_LEVEL_2, routineName, "Error writing file");

	if (Nwrite != backgroundRecords(nspect))
		Error(ERROR_LEVEL_2, routineName,
			  "Number of records written differs from backgroundRecords");

	/* --- Return the number of written records --       -------------- */

	return Nwrite;
}
/* ------- end ---------------------------- writeBackground.c ------- */

/* ------- begin -------------------------- readProfile.c ----------- */

void readProfile(AtomicLine *line, int lamu, float *phi)
{
	const char routineName[] = "readProfile";

	int    Nrecphi, NrecSkip;
	bool_t result = TRUE;
	size_t recordsize;
	off_t  offset;

	if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
		if (input.magneto_optical)
			NrecSkip = 7;
		else
			NrecSkip = 4;
	} else
		NrecSkip = 1;

	if (line->polarizable && input.StokesMode == FULL_STOKES) {
		if (input.magneto_optical)
			Nrecphi = 7;
		else
			Nrecphi = 4;
	} else
		Nrecphi = 1;

	recordsize = Nrecphi * atmos.Nspace * sizeof(float);
	offset     = NrecSkip * atmos.Nspace * sizeof(float) * (off_t) lamu;

	result &= (pread_rh(line->fd_profile, phi, recordsize, offset) ==
			   recordsize);

	if (!result) Error(ERROR_LEVEL_2, routineName, "Error reading file");
}
/* ------- end ---------------------------- readProfile.c ----------- */

/* ------- begin -------------------------- writeProfile.c ---------- */

void writeProfile(AtomicLine *line, int lamu, float *phi)
{
	const char routineName[] = "writeProfile";

	int    Nrecphi;
	bool_t result = TRUE;
	size_t recordsize;
	off_t  offset;

	if (line->polarizable && (input.StokesMode > FIELD_FREE)) {
		if (input.magneto_optical)
			Nrecphi = 7;
		else
			Nrecphi = 4;
	} else
		Nrecphi = 1;

	recordsize = Nrecphi * atmos.Nspace * sizeof(float);
	offset     = recordsize * (off_t) lamu;

	result &= (pwrite_rh(line->fd_profile, phi, recordsize, offset) ==
			   recordsize);

	if (!result) Error(ERROR_LEVEL_2, routineName, "Error writing file");
}
/* ------- end ---------------------------- writeProfile.c ---------- */

/* ------- begin -------------------------- fullIO.c ---------------- */

static size_t fullIO(int fd, char *buf, size_t count, off_t offset,
					 bool_t write_it)
{
	/* --- Transfer count bytes, in pieces of at most PWRITE_SLICE, and
		   retry partial transfers and interrupted calls. Returns count
		   on success and 0 on failure, with errno set --  ------------ */

	size_t  left = count, chunk;
	ssize_t n;

	while (left > 0) {
		chunk = MIN(left, (size_t) PWRITE_SLICE);
		n     = (write_it) ? pwrite(fd, buf, chunk, offset) : pread(fd, buf, chunk, offset);
		if (n < 0 && errno == EINTR) continue;
		if (n <= 0) {
			if (n == 0) errno = EIO; /* unexpected end of file */
			return 0;
		}
		buf += n;
		offset += n;
		left -= n;
	}
	return count;
}
/* ------- end ---------------------------- fullIO.c ---------------- */

/* ------- begin -------------------------- pread_rh.c -------------- */

size_t pread_rh(int fd, const void *buffer, size_t count, off_t offset)
{
	return fullIO(fd, (char *) buffer, count, offset, FALSE);
}
/* ------- end ---------------------------- pread_rh.c -------------- */

/* ------- begin -------------------------- pwrite_rh.c ------------- */

size_t pwrite_rh(int fd, const void *buffer, size_t count, off_t offset)
{
	return fullIO(fd, (char *) buffer, count, offset, TRUE);
}
/* ------- end ---------------------------- pwrite_rh.c ------------- */

/* ------- begin -------------------------- addJgas.c --------------- */

/* --- The gas-frame mean intensity Jgas of the hybrid PRD approximation
	   (PRD_ANGLE_APPROX). The formal solution of one wavelength adds its
	   intensity, shifted to the gas frame, to the rows of neighbouring
	   wavelengths, and with several threads neighbouring wavelengths are
	   solved at the same time. Each row therefore has a lock (T-4).
	   A ray's contributions are first summed in a buffer of the thread
	   (they fall in a narrow band of rows around its own wavelength) and
	   then added row by row: a few lock operations per ray instead of
	   one per point and coefficient. --               -------------- */

#define JGAS_MAXROWS 64 /* widest band that is buffered */

static __thread double      *jgas_buf;
static __thread char        *jgas_used;
static __thread int          jgas_bufrows;
static __thread PRDhBracket *jgas_br; /* brackets of one ray, PRDH_LIMIT_MEM */
static __thread int          jgas_brN;

void initJgasLocks(void)
{
	int nspect;

	spectrum.Jgas_lock =
		(pthread_mutex_t *) malloc(spectrum.Nspect * sizeof(pthread_mutex_t));
	if (spectrum.Jgas_lock == NULL)
		Error(ERROR_LEVEL_2, "initJgasLocks", "Out of memory");
	for (nspect = 0; nspect < spectrum.Nspect; nspect++)
		pthread_mutex_init(&spectrum.Jgas_lock[nspect], NULL);
}

void freeJgasBuffer(void)
{
	/* --- Worker threads end after each formal solution of the
		   spectrum: free their buffers --               -------------- */

	free(jgas_buf);
	free(jgas_used);
	free(jgas_br);
	jgas_buf     = NULL;
	jgas_used    = NULL;
	jgas_br      = NULL;
	jgas_bufrows = 0;
	jgas_brN     = 0;
}

/* --- A PRD sub-iteration re-solves only the wavelengths with a PRD
	   line, but Jgas of those wavelengths also holds the contributions of
	   their neighbours without one (N-27). In a moving atmosphere the
	   contributions of the wavelengths without PRD line to the rows of
	   wavelengths with PRD line therefore go into their own array,
	   Jgas_fixed, which the sub-iterations leave alone: Jgas of a PRD
	   wavelength = Jgas_fixed + the contributions of the wavelengths with
	   PRD line. Jgas_fixed has a row only for the PRD wavelengths that a
	   wavelength without PRD line can reach (spectrum.PRDrow, see
	   initJgasRows). Everything else, and everything without
	   spectrum.Jgas_fixed, goes into Jgas. --            -------------- */

void initJgasRows(const char *isPRD)
{
	int         nspect, mu, k, Nrows = 0, moving = FALSE;
	long        n;
	double      vmax = 0.0, shift_lo, shift_hi;
	char       *reach;
	PRDhBracket lo, hi;

	/* --- PRDrow: -1 without PRD line; a row for a PRD wavelength that a
		   wavelength without PRD line can reach, -2 for one it cannot.
		   The shift of a ray at a location is 1 +- v_los / c; the points
		   a wavelength gives to run from r[0] to r[2] of its bracket, and
		   both ends do not decrease with the shift, so the range at the
		   smallest shift to the one at the largest holds all of them.
		   v_los is fixed in a run. --                   -------------- */

	free(spectrum.PRDrow);
	if (spectrum.Jgas_fixed != NULL) freeMatrix((void **) spectrum.Jgas_fixed);
	spectrum.Jgas_fixed = NULL;
	spectrum.NJgasFixed = 0;
	spectrum.PRDrow     = (int *) malloc(spectrum.Nspect * sizeof(int));
	if (spectrum.PRDrow == NULL) Error(ERROR_LEVEL_2, "initJgasRows", "Out of memory");
	for (mu = 0; mu < atmos.Nrays; mu++)
		for (k = 0; k < atmos.Nspace; k++)
			vmax = MAX(vmax, fabs(spectrum.v_los[mu][k]));
	moving   = (vmax != 0.0);
	shift_lo = 1.0 - vmax / CLIGHT;
	shift_hi = 1.0 + vmax / CLIGHT;

	reach = (char *) calloc(spectrum.Nspect, 1);
	if (reach == NULL) Error(ERROR_LEVEL_2, "initJgasRows", "Out of memory");
	if (moving) {
		for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
			if (isPRD[nspect]) continue;
			prdhBracket(nspect, shift_lo, &lo);
			prdhBracket(nspect, shift_hi, &hi);
			for (n = lo.r[0]; n < hi.r[2]; n++) reach[n] = TRUE;
		}
	}
	for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
		if (!isPRD[nspect])
			spectrum.PRDrow[nspect] = -1;
		else
			spectrum.PRDrow[nspect] = reach[nspect] ? Nrows++ : -2;
	}
	free(reach);
	if (Nrows > 0) spectrum.Jgas_fixed = matrix_double(Nrows, atmos.Nspace);
	spectrum.NJgasFixed = Nrows;
}

static double *jgasRow(int src, int row)
{
	if (spectrum.Jgas_fixed != NULL && spectrum.PRDrow[src] == -1 && spectrum.PRDrow[row] != -1) {
		if (spectrum.PRDrow[row] == -2) {

			/* --- Called by worker threads, never with a lock held; the
				   message goes in a buffer of its own, not the shared
				   messageStr --                         -------------- */

			char msg[256];

			snprintf(msg, sizeof(msg),
					 "Wavelength %d (%g nm) without PRD line gives to wavelength %d (%g nm) "
					 "with PRD line, which the reach computed at the start excludes",
					 src, spectrum.lambda[src], row, spectrum.lambda[row]);
			Error(ERROR_LEVEL_2, "jgasRow", msg);
		}
		return spectrum.Jgas_fixed[spectrum.PRDrow[row]];
	}
	return spectrum.Jgas[row];
}

void beginJgas(bool_t redistribute)
{
	int    nspect;
	size_t n = (size_t) spectrum.Nspect * atmos.Nspace;

	/* --- Before a formal solution of all wavelengths (redistribute ==
		   FALSE) Jgas and Jgas_fixed start at zero; before one of the PRD
		   wavelengths only, the rows of Jgas of the PRD wavelengths that
		   have one in Jgas_fixed start from it --       -------------- */

	memset(spectrum.Jgas[0], 0, n * sizeof(double));
	if (spectrum.Jgas_fixed == NULL) return;
	if (!redistribute) {
		memset(spectrum.Jgas_fixed[0], 0, (size_t) spectrum.NJgasFixed * atmos.Nspace * sizeof(double));
		return;
	}
	for (nspect = 0; nspect < spectrum.Nspect; nspect++)
		if (spectrum.PRDrow[nspect] >= 0)
			memcpy(spectrum.Jgas[nspect], spectrum.Jgas_fixed[spectrum.PRDrow[nspect]],
				   atmos.Nspace * sizeof(double));
}

void endJgas(bool_t redistribute)
{
	int nspect, k;

	/* --- After a full solution add the part of the wavelengths without
		   PRD line to the rows of the PRD wavelengths --  -------------- */

	if (spectrum.Jgas_fixed == NULL || redistribute) return;
	for (nspect = 0; nspect < spectrum.Nspect; nspect++) {
		if (spectrum.PRDrow[nspect] < 0) continue;
		for (k = 0; k < atmos.Nspace; k++)
			spectrum.Jgas[nspect][k] += spectrum.Jgas_fixed[spectrum.PRDrow[nspect]][k];
	}
}

static void addJgas(double *row, int nspect, int k, double value)
{
	pthread_mutex_lock(&spectrum.Jgas_lock[nspect]);
	row[k] += value;
	pthread_mutex_unlock(&spectrum.Jgas_lock[nspect]);
}

void addJgasRow(int nspect, double wmu, double *I)
{
	int     k;
	double *row = jgasRow(nspect, nspect);

	pthread_mutex_lock(&spectrum.Jgas_lock[nspect]);
	for (k = 0; k < atmos.Nspace; k++) row[k] += wmu * I[k];
	pthread_mutex_unlock(&spectrum.Jgas_lock[nspect]);
}
/* ------- end ---------------------------- addJgas.c --------------- */

/* ------- begin -------------------------- prdhIndex.c ------------- */

long prdhIndex(int nspect, int mu, int to_obs, int k)
{
	/* --- Index of (wavelength, ray, direction, location) in the hybrid
		   PRD table spectrum.nc, in long arithmetic: the table has
		   2 Nrays Nspect Nspace entries (B-32) --         -------------- */

	return (((long) nspect * atmos.Nrays + mu) * 2 + to_obs) * atmos.Nspace + k;
}
/* ------- end ---------------------------- prdhIndex.c ------------- */

/* ------- begin -------------------------- addJgasRay.c ------------ */

/* --- Grid points (a range per location) and weights of one ray: from
	   the table spectrum.nc, iprdh and cprdh, or, with PRDH_LIMIT_MEM,
	   from prdhBracket and prdhWeight, which the table holds (P-17) -- */

static inline long jgasFirst(long lamuk0, int k)
{
	if (input.prdh_limit_mem) return jgas_br[k].r[0];
	return (lamuk0 + k == 0) ? 0 : spectrum.nc[lamuk0 + k - 1];
}

static inline long jgasLast(long lamuk0, int k)
{
	if (input.prdh_limit_mem) return jgas_br[k].r[2];
	return spectrum.nc[lamuk0 + k];
}

static inline void jgasPoint(int k, long n, int *idx, double *w)
{
	/* --- n runs from jgasFirst to jgasLast: a grid index with
		   PRDH_LIMIT_MEM, a position in the table otherwise -- ---- */

	if (input.prdh_limit_mem) {
		*idx = (int) n;
		*w   = prdhWeight(&jgas_br[k], n);
	} else {
		*idx = spectrum.iprdh[n];
		*w   = spectrum.cprdh[n];
	}
}

void addJgasRay(int nspect, int mu, int to_obs, double wmu, double *I)
{
	/* --- Add the intensity I of ray (nspect, mu, to_obs) to Jgas.
		   With and without PRDH_LIMIT_MEM the same points get the same
		   weights, added in the same order --           -------------- */

	int      Nspace = atmos.Nspace, k, row, imin, imax, Nrow, idx;
	long int lamuk0 = 0, n, n0, n1;
	double  *buf, *jrow, w;

	if (input.prdh_limit_mem) {
		if (jgas_brN < Nspace) {
			free(jgas_br);
			jgas_br  = (PRDhBracket *) malloc(Nspace * sizeof(PRDhBracket));
			jgas_brN = Nspace;
			if (jgas_br == NULL) Error(ERROR_LEVEL_2, "addJgasRay", "Out of memory");
		}
		for (k = 0; k < Nspace; k++)
			prdhBracket(nspect, prdhShift(mu, to_obs, k), &jgas_br[k]);
	} else
		lamuk0 = prdhIndex(nspect, mu, to_obs, 0);

	imin = INT_MAX;
	imax = -1;
	for (k = 0; k < Nspace; k++) {
		n0 = jgasFirst(lamuk0, k);
		n1 = jgasLast(lamuk0, k);
		for (n = n0; n < n1; n++) {
			jgasPoint(k, n, &idx, &w);
			if (idx < imin) imin = idx;
			if (idx > imax) imax = idx;
		}
	}
	if (imax < 0) return;
	Nrow = imax - imin + 1;

	/* --- A band wider than JGAS_MAXROWS (not seen in practice: the
		   widest in the test cases has 34 rows) is added point by
		   point, to bound the buffer --                 -------------- */

	if (Nrow > JGAS_MAXROWS) {

		/* --- The target row of each index is resolved on first use -- */

		double **rowp = (double **) calloc(Nrow, sizeof(double *));

		if (rowp == NULL) Error(ERROR_LEVEL_2, "addJgasRay", "Out of memory");
		for (k = 0; k < Nspace; k++) {
			n1 = jgasLast(lamuk0, k);
			for (n = jgasFirst(lamuk0, k); n < n1; n++) {
				jgasPoint(k, n, &idx, &w);
				if (rowp[idx - imin] == NULL) rowp[idx - imin] = jgasRow(nspect, idx);
				addJgas(rowp[idx - imin], idx, k, wmu * w * I[k]);
			}
		}
		free(rowp);
		return;
	}

	if (Nrow > jgas_bufrows) {
		free(jgas_buf);
		free(jgas_used);
		jgas_buf     = (double *) malloc((size_t) Nrow * Nspace * sizeof(double));
		jgas_used    = (char *) malloc(Nrow);
		jgas_bufrows = Nrow;
		if (jgas_buf == NULL || jgas_used == NULL)
			Error(ERROR_LEVEL_2, "addJgasRay", "Out of memory");
	}
	memset(jgas_used, 0, Nrow);

	for (k = 0; k < Nspace; k++) {
		n1 = jgasLast(lamuk0, k);
		for (n = jgasFirst(lamuk0, k); n < n1; n++) {
			jgasPoint(k, n, &idx, &w);
			row = idx - imin;
			if (!jgas_used[row]) {

				/* --- Zero a row when the ray first reaches it -- ---- */

				memset(jgas_buf + (size_t) row * Nspace, 0, Nspace * sizeof(double));
				jgas_used[row] = 1;
			}
			jgas_buf[(size_t) row * Nspace + k] += wmu * w * I[k];
		}
	}
	for (row = 0; row < Nrow; row++) {
		if (!jgas_used[row]) continue;
		buf = jgas_buf + (size_t) row * Nspace;

		jrow = jgasRow(nspect, imin + row);
		pthread_mutex_lock(&spectrum.Jgas_lock[imin + row]);
		for (k = 0; k < Nspace; k++) jrow[k] += buf[k];
		pthread_mutex_unlock(&spectrum.Jgas_lock[imin + row]);
	}
}
/* ------- end ---------------------------- addJgasRay.c ------------ */
