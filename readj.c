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
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "inputs.h"
#include "error.h"

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

/* ------- begin -------------------------- backgroundIO.c ---------- */

static float *backgroundBuffer(void)
{
	float *buffer = NULL;

	/* --- Conversion buffer for backgroundIO: the largest part of a
		   block is 4 records (the Stokes opacity or emissivity) -- -- */

	if (atmos.backgrbytes == sizeof(float)) {
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

static __thread double *jgas_buf;
static __thread char   *jgas_used;
static __thread int     jgas_bufrows;

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
	jgas_buf     = NULL;
	jgas_used    = NULL;
	jgas_bufrows = 0;
}

void addJgas(int nspect, int k, double value)
{
	pthread_mutex_lock(&spectrum.Jgas_lock[nspect]);
	spectrum.Jgas[nspect][k] += value;
	pthread_mutex_unlock(&spectrum.Jgas_lock[nspect]);
}

void addJgasRow(int nspect, double wmu, double *I)
{
	int k;

	pthread_mutex_lock(&spectrum.Jgas_lock[nspect]);
	for (k = 0; k < atmos.Nspace; k++) spectrum.Jgas[nspect][k] += wmu * I[k];
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

void addJgasRay(int nspect, int mu, int to_obs, double wmu, double *I)
{
	/* --- Add the intensity I of ray (nspect, mu, to_obs) to Jgas with
		   the precomputed coefficients spectrum.nc, iprdh and cprdh -- */

	int      Nspace = atmos.Nspace, k, row, imin, imax, Nrow;
	long int lamuk0, lamuk, idx, idx0, idx1;
	double  *buf;

	lamuk0 = prdhIndex(nspect, mu, to_obs, 0);
	idx0   = (lamuk0 == 0) ? 0 : spectrum.nc[lamuk0 - 1];
	idx1   = spectrum.nc[lamuk0 + Nspace - 1];
	if (idx1 <= idx0) return;

	imin = imax = spectrum.iprdh[idx0];
	for (idx = idx0; idx < idx1; idx++) {
		if (spectrum.iprdh[idx] < imin) imin = spectrum.iprdh[idx];
		if (spectrum.iprdh[idx] > imax) imax = spectrum.iprdh[idx];
	}
	Nrow = imax - imin + 1;

	/* --- A band wider than JGAS_MAXROWS (not seen in practice: the
		   widest in the test cases has 34 rows) is added point by
		   point, to bound the buffer --                 -------------- */

	if (Nrow > JGAS_MAXROWS) {
		for (k = 0; k < Nspace; k++) {
			lamuk = lamuk0 + k;
			for (idx = (lamuk == 0) ? 0 : spectrum.nc[lamuk - 1];
				 idx < spectrum.nc[lamuk]; idx++)
				addJgas(spectrum.iprdh[idx], k, wmu * spectrum.cprdh[idx] * I[k]);
		}
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
		lamuk = lamuk0 + k;
		for (idx = (lamuk == 0) ? 0 : spectrum.nc[lamuk - 1];
			 idx < spectrum.nc[lamuk]; idx++) {
			row = spectrum.iprdh[idx] - imin;
			if (!jgas_used[row]) {

				/* --- Zero a row when the ray first reaches it -- ---- */

				memset(jgas_buf + (size_t) row * Nspace, 0, Nspace * sizeof(double));
				jgas_used[row] = 1;
			}
			jgas_buf[(size_t) row * Nspace + k] += wmu * spectrum.cprdh[idx] * I[k];
		}
	}
	for (row = 0; row < Nrow; row++) {
		if (!jgas_used[row]) continue;
		buf = jgas_buf + (size_t) row * Nspace;

		pthread_mutex_lock(&spectrum.Jgas_lock[imin + row]);
		for (k = 0; k < Nspace; k++) spectrum.Jgas[imin + row][k] += buf[k];
		pthread_mutex_unlock(&spectrum.Jgas_lock[imin + row]);
	}
}
/* ------- end ---------------------------- addJgasRay.c ------------ */
