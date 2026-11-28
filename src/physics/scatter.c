/* ------- file: -------------------------- scatter.c ---------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Feb 22 17:02:07 2024 --

	   --------------------------                      ----------RH-- */

/* --- Evaluate scattering integral for PRD transition pointed to
	   by transition pointer PRDline.

	   Adapted for Cross-Redistribution by:
		  Eliza Miller-Ricci (Middlebury College), Jun 29 2001


 Note: Scratch files for the redistribution weights are now written
	   in a location determined from PRD_FILE_TEMPLATE. The files are
	   no longer automatically deleted. Careful, they can be quite big.
	   --                                              -------------- */

#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "inputs.h"
#include "constant.h"
#include "error.h"
#include "statistics.h"
#include "parallel.h"
#include "prdfold.h"
#include "state.h"

#define PRD_FILE_TEMPLATE "PRD_%s_%d-%d.dat"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- prdInflow.c -------------- */

static void prdInflow(RHState *rh, Atom *atom, AtomicLine *PRDline, AtomicLine **XRD,
					  int Nsubordinate, double *Yj)
{
	register int k, kr, ip, kxrd;
	bool_t       subordinate;

	AtomicLine      *line;
	AtomicContinuum *continuum;

	/* --- Rate Yj into the upper level j of PRDline by every process
		   except radiative excitation in the subordinate lines XRD,
		   plus the elastic collisions n_j Q_elast. In statistical
		   equilibrium n_j P_j = Yj + sum_x n_i(x) B_ij(x) Jbar(x), so the
		   emission profile ratio

			 rho = (1 + sum_x g_x <J>_R(x)) / (1 + sum_x g_x Jbar(x)),
			 g_x = n_i(x) B_ij(x) / Yj,

		   equals the familiar 1 + sum_x gamma_x (<J>_R(x) - Jbar(x)),
		   gamma_x = n_i(x) B_ij(x) / (n_j P_j). Unlike that form it
		   stays positive when the populations lag the radiation field
		   during the PRD sub-iterations (B-29). --        -------------- */

	for (k = 0; k < rh->atmos.Nspace; k++)
		Yj[k] = atom->n[PRDline->j][k] * PRDline->Qelast[k];
	for (ip = 0; ip < atom->Nlevel; ip++) {
		if (ip == PRDline->j) continue;
		for (k = 0; k < rh->atmos.Nspace; k++)
			Yj[k] += atom->n[ip][k] * atom->C[PRDline->j * atom->Nlevel + ip][k];
	}
	for (kr = 0; kr < atom->Nline; kr++) {
		line = &atom->line[kr];
		if (line->j == PRDline->j) {
			subordinate = FALSE;
			for (kxrd = 0; kxrd < Nsubordinate; kxrd++)
				if (XRD[kxrd] == line) subordinate = TRUE;
			if (!subordinate)
				for (k = 0; k < rh->atmos.Nspace; k++)
					Yj[k] += atom->n[line->i][k] * line->Rij[k];
		}
		if (line->i == PRDline->j)
			for (k = 0; k < rh->atmos.Nspace; k++)
				Yj[k] += atom->n[line->j][k] * line->Rji[k];
	}
	for (kr = 0; kr < atom->Ncont; kr++) {
		continuum = &atom->continuum[kr];
		if (continuum->j == PRDline->j)
			for (k = 0; k < rh->atmos.Nspace; k++)
				Yj[k] += atom->n[continuum->i][k] * continuum->Rij[k];
		if (continuum->i == PRDline->j)
			for (k = 0; k < rh->atmos.Nspace; k++)
				Yj[k] += atom->n[continuum->j][k] * continuum->Rji[k];
	}
}
/* ------- end ---------------------------- prdInflow.c -------------- */

/* ------- begin -------------------------- folded weights ---------- */

/* --- The scattering integrals (P-8, P-14). For every emission
	   wavelength the integral

		 <J> = sum_q' g(q') J(q') / sum g

	   runs over a fine absorption grid q' (step PRD_DQ, see prdGrid),
	   with J linearly interpolated from the line's own wavelength grid
	   (constant beyond its ends, as in Linear). That makes it linear in
	   the J values on the line's grid:

		 <J> = sum_m W(m) J_m,   W(m) = sum_q' g(q') c(q', m) / sum g,

	   c(q', m) the interpolation coefficients. g and c do not change
	   between calls (the absorption grid of a location is fixed), so W
	   is built once, when the weights are first computed, and every
	   later call is one short dot product per emission wavelength
	   instead of an interpolation onto the fine grid. W(m) is non-zero
	   only for m in a range [m0, m0 + len) of the line's grid, which
	   prdFoldRange computes from the ends of the fine grid on every
	   call (two bisections), so only W itself is stored, plus the size
	   of each location's block (line->prd_block).

	   W is kept in a scratch file per line in single precision (the
	   dot product is accumulated in double and divided by the sum of
	   the stored weights of the row, N-56). The file holds one block
	   per location; workers handle contiguous ranges of locations and
	   read or write their blocks in windows of several MB.
	   --                                              -------------- */

typedef float prdw_t; /* stored folded weights */

/* ------- begin -------------------------- prdFileIO.c ------------- */

static bool_t prdFileIO(int fd, prdw_t *buf, size_t Nvalue, off_t offset,
						bool_t write_it)
{
	char   *p    = (char *) buf;
	size_t  left = Nvalue * sizeof(prdw_t);
	ssize_t n;

	while (left > 0) {

		/* --- A single call may not take more than INT_MAX bytes on
			   some systems, and may be interrupted --     -------------- */

		size_t chunk = (left < (size_t) 1 << 30) ? left : (size_t) 1 << 30;

		n = (write_it) ? pwrite(fd, p, chunk, offset) : pread(fd, p, chunk, offset);
		if (n < 0 && errno == EINTR) continue;
		if (n <= 0) return FALSE;
		p += n;
		offset += n;
		left -= n;
	}
	return TRUE;
}
/* ------- end ---------------------------- prdFileIO.c ------------- */

/* --- Windows of consecutive locations, several MB at a time: large
	   requests are what the disk and the OS cache handle well. With
	   one small block per location, runs whose weight files did not
	   stay in the cache became limited by disk latency. -- --------- */

#define PRD_IO_WINDOW (2 * 1024 * 1024)

typedef struct {
	int     fd;
	char    name[MAX_LINE_SIZE]; /* file name, for messages */
	bool_t  initialize;
	long   *blocksize; /* number of weights in each location's block */
	off_t  *kstart;    /* byte offset of each location's block       */
	int     kbeg, kend, kwin_start, kwin_end;
	size_t  cap;
	prdw_t *win;
} PRDWindow;

static void prdFileName(AtomicLine *PRDline, char *filename)
{
	Atom *atom = PRDline->atom;

	sprintf(filename, (atom->ID[1] == ' ') ? "PRD_%.1s_%d-%d.dat" : "PRD_%s_%d-%d.dat",
			atom->ID, PRDline->j, PRDline->i);
}

static void prdWindowOpen(PRDWindow *w, AtomicLine *PRDline, bool_t initialize,
						  long *blocksize, off_t *kstart, int kbeg, int kend)
{
	w->fd = fileno(PRDline->fp_GII);
	prdFileName(PRDline, w->name);
	w->initialize = initialize;
	w->blocksize  = blocksize;
	w->kstart     = kstart;
	w->kbeg       = kbeg;
	w->kend       = kend;
	w->kwin_start = w->kwin_end = kbeg;
	w->cap                      = 0;
	w->win                      = NULL;
}

static rh_status prdWindowFlush(RHState *rh, PRDWindow *w, const char *routineName)
{
	size_t N;

	if (!w->initialize || w->kwin_end <= w->kwin_start) return RH_OK;
	N = (w->kstart[w->kwin_end - 1] - w->kstart[w->kwin_start]) / sizeof(prdw_t) +
		w->blocksize[w->kwin_end - 1];
	if (!prdFileIO(w->fd, w->win, N, w->kstart[w->kwin_start], TRUE)) {
		return Fail(rh, RH_E_DATA, routineName, "Unable to write the redistribution weights of locations %d to %d"
												" to %s",
					w->kwin_start, w->kwin_end - 1, w->name);
	}
	return RH_OK;
}

/* --- The weights of location k: starts a new window (writing the
	   finished one first, or reading the next one) when needed -- --- */

static rh_status prdWindowBlock(RHState *rh, PRDWindow *w, int k, const char *routineName,
								prdw_t **block)
{
	size_t  bytes, N;
	prdw_t *win;

	if (k == w->kwin_end) {
		RH_TRY(prdWindowFlush(rh, w, routineName));
		w->kwin_start = k;
		bytes         = 0;
		w->kwin_end   = k;
		while (w->kwin_end < w->kend &&
			   (w->kwin_end == k ||
				bytes + w->blocksize[w->kwin_end] * sizeof(prdw_t) <= PRD_IO_WINDOW)) {
			bytes += w->blocksize[w->kwin_end] * sizeof(prdw_t);
			w->kwin_end++;
		}
		N = (w->kstart[w->kwin_end - 1] - w->kstart[k]) / sizeof(prdw_t) +
			w->blocksize[w->kwin_end - 1];

		/* --- One buffer per worker, sized for a full window and kept:
			   growing it window by window, with sizes that differ a
			   little each time, made the allocator hold on to several
			   windows' worth of freed blocks. --        -------------- */

		if (N > w->cap) {
			w->cap = MAX(N, PRD_IO_WINDOW / sizeof(prdw_t));
			win    = (prdw_t *) realloc(w->win, w->cap * sizeof(prdw_t));
			if (win == NULL) return Fail(rh, RH_E_NOMEM, routineName, "Out of memory");
			w->win = win;
		}
		if (!w->initialize && !prdFileIO(w->fd, w->win, N, w->kstart[k], FALSE)) {
			return Fail(rh, RH_E_DATA, routineName, "Unable to read the redistribution weights of locations %d to %d"
													" from %s",
						k, w->kwin_end - 1, w->name);
		}
	}
	*block = w->win + (w->kstart[k] - w->kstart[w->kwin_start]) / sizeof(prdw_t);
	return RH_OK;
}

/* --- A location's rows do not add up to its block: the ranges
	   differ from when the weights were computed --    -------------- */

static rh_status prdBlockError(RHState *rh, const char *routineName, int k, long n, long size)
{
	return Fail(rh, RH_E_NUMERIC, routineName, "Inconsistent block of redistribution weights at"
											   " location %d: %ld instead of %ld",
				k, n, size);
}

/* --- Block offsets in the file from the block sizes -- ---------- */

static off_t prdOffsets(RHState *rh, long *blocksize, off_t *kstart, off_t offset)
{
	int k;

	for (k = 0; k < rh->atmos.Nspace; k++) {
		kstart[k] = offset;
		offset += blocksize[k] * sizeof(prdw_t);
	}
	return offset;
}

/* --- Opens the scratch file of the weights of PRDline on the first
	   call; *initialize is TRUE when the weights still have to be
	   computed --                                      -------------- */

static rh_status prdOpenFile(RHState *rh, AtomicLine *PRDline, const char *routineName,
							 bool_t *initialize)
{
	char filename[MAX_LINE_SIZE];

	*initialize = FALSE;
	if (PRDline->fp_GII != NULL) return RH_OK;
	prdFileName(PRDline, filename);
	if ((PRDline->fp_GII = fopen(filename, "w+")) == NULL) {
		return Fail(rh, RH_E_DATA, routineName, "Unable to open temporary file %s", filename);
	}
	*initialize = TRUE;
	return RH_OK;
}
/* ------- end ---------------------------- folded weights ---------- */

/* ------- begin -------------------------- PRDScatter.c ------------ */

/* --- After a failure in the call that opened the scratch file of the
	   weights: the file is closed and the block sizes dropped, so that
	   the next call computes the weights again -- ---------------- */

static void prdDropFile(AtomicLine *PRDline)
{
	fclose(PRDline->fp_GII);
	PRDline->fp_GII = NULL;
	free(PRDline->prd_block);
	PRDline->prd_block = NULL;
}

/* --- Angle-averaged and hybrid (PRD_ANGLE_APPROX) PRD. The loop over
	   locations is threaded: the integral at one location uses only
	   that location's data and writes only rho_prd[.][k]. The file holds
	   the subordinate lines one after the other, and within each the
	   blocks of the locations, each with the rows of all emission
	   wavelengths in order. --                         -------------- */

typedef struct {
	Atom       *atom;
	AtomicLine *PRDline, *XRDline;
	bool_t      initialize;
	double      waveratio, *Yj, *adamp;
	long       *blocksize;
	off_t      *kstart;
} PRDjob;

/* --- Absorption grid of XRDline at location k, in Doppler units -- */

static void prdAbsGrid(Atom *atom, AtomicLine *XRDline, int k, double *q_abs)
{
	int la;

	for (la = 0; la < XRDline->Nlambda; la++)
		q_abs[la] = (XRDline->lambda[la] - XRDline->lambda0) * CLIGHT /
					(XRDline->lambda0 * atom->vbroad[k]);
}

/* ------- begin -------------------------- prdBlockSizes.c --------- */

static void prdBlockSizes(RHState *rh, int kbeg, int kend, int nt, void *arg)
{
	PRDjob     *job     = (PRDjob *) arg;
	AtomicLine *PRDline = job->PRDline, *XRDline = job->XRDline;
	int         la, k, Np, m0;
	long        n;
	double      q_emit, q0, qN, *q_abs;

	(void) rh;
	(void) nt;

	q_abs = (double *) mallocOrAbort(XRDline->Nlambda * sizeof(double));
	for (k = kbeg; k < kend; k++) {
		prdAbsGrid(job->atom, XRDline, k, q_abs);
		n = 0;
		for (la = 0; la < PRDline->Nlambda; la++) {
			q_emit = (PRDline->lambda[la] - PRDline->lambda0) * CLIGHT /
					 (PRDline->lambda0 * job->atom->vbroad[k]);
			Np = prdGrid(q_emit, job->waveratio, &q0, &qN);
			n += prdFoldRange(XRDline->Nlambda, q_abs, q0, Np, XRDline->symmetric, &m0);
		}
		job->blocksize[k] = n;
	}
	free(q_abs);
}
/* ------- end ---------------------------- prdBlockSizes.c --------- */

/* ------- begin -------------------------- prdIntegrate.c ---------- */

static rh_status prdIntegrateRange(RHState *rh, int kbeg, int kend, int nt, void *arg)
{
	const char routineName[] = "PRDScatter";

	PRDjob      *job     = (PRDjob *) arg;
	Atom        *atom    = job->atom;
	AtomicLine  *PRDline = job->PRDline, *XRDline = job->XRDline;
	register int la, k, lap;
	int          Np, m0, len, m;
	long         pos;
	double       q_emit, q0, qN, *q_abs, *qp = NULL, *qu = NULL, *gw = NULL, *W, *J_k,
								   gnorm, g, scatInt, wsum;
	prdw_t   *wbuf;
	PRDWindow win;
	rh_status result = RH_OK;

	(void) nt;

	J_k   = (double *) mallocOrAbort(XRDline->Nlambda * sizeof(double));
	q_abs = (double *) mallocOrAbort(XRDline->Nlambda * sizeof(double));
	W     = (double *) mallocOrAbort(XRDline->Nlambda * sizeof(double));
	prdWindowOpen(&win, job->PRDline, job->initialize, job->blocksize, job->kstart, kbeg, kend);

	for (k = kbeg; k < kend; k++) {
		RH_TRY_GOTO(result, cleanup, prdWindowBlock(rh, &win, k, routineName, &wbuf));
		g = atom->n[XRDline->i][k] * XRDline->Bij / job->Yj[k];

		/* --- Local mean intensity (in the gas frame for hybrid PRD)
			   and absorption grid in Doppler units --     -------------- */

		for (la = 0; la < XRDline->Nlambda; la++)
			J_k[la] = (rh->input.PRD_angle_dep == PRD_ANGLE_APPROX)
						  ? rh->spectrum.Jgas[XRDline->Nblue + la][k]
						  : rh->spectrum.J[XRDline->Nblue + la][k];
		prdAbsGrid(atom, XRDline, k, q_abs);

		/* --- Loop over emission wavelengths --         -------------- */

		pos = 0;
		for (la = 0; la < PRDline->Nlambda; la++) {
			q_emit = (PRDline->lambda[la] - PRDline->lambda0) * CLIGHT /
					 (PRDline->lambda0 * atom->vbroad[k]);
			Np  = prdGrid(q_emit, job->waveratio, &q0, &qN);
			len = prdFoldRange(XRDline->Nlambda, q_abs, q0, Np, XRDline->symmetric, &m0);
			if (pos + len > job->blocksize[k])
				RH_TRY_GOTO(result, cleanup, prdBlockError(rh, routineName, k, pos + len, job->blocksize[k]));

			if (job->initialize) {

				/* --- Redistribution weights on the fine grid, with
					   the integration weights of Press et al., Numerical
					   Recipes, p. 107, eq. 4.1.12, folded with the
					   interpolation coefficients and normalized -- --- */

				gw = (double *) reallocOrAbort(gw, Np * sizeof(double));
				qp = (double *) reallocOrAbort(qp, Np * sizeof(double));
				qu = (double *) reallocOrAbort(qu, Np * sizeof(double));
				prdFineGrid(q0, Np, XRDline->symmetric, qp);
				prdFineGrid(q0, Np, FALSE, qu); /* GII: unfolded grid */
				gnorm = 0.0;
				for (lap = 0; lap < Np; lap++) {
					gw[lap] = GII(job->adamp[k], job->waveratio, q_emit, qu[lap]) *
							  prdQuadWeight(lap, Np);
					gnorm += gw[lap];
				}
				for (m = 0; m < len; m++) W[m] = 0.0;
				prdFold(rh, XRDline->Nlambda, q_abs, Np, qp, gw, m0, len, W);
				for (m = 0; m < len; m++) wbuf[pos + m] = (prdw_t) (W[m] / gnorm);
			}

			/* --- The integral: a dot product with J on the line's
				   grid, accumulated in double, divided by the sum of
				   the stored weights: the row is 1 to double
				   rounding, not to that of prdw_t (N-56) -- -------------- */

			scatInt = 0.0;
			wsum    = 0.0;
			for (m = 0; m < len; m++) {
				scatInt += (double) wbuf[pos + m] * J_k[m0 + m];
				wsum += (double) wbuf[pos + m];
			}
			PRDline->rho_prd[la][k] += g * scatInt / wsum;
			pos += len;
		}
		if (pos != job->blocksize[k])
			RH_TRY_GOTO(result, cleanup, prdBlockError(rh, routineName, k, pos, job->blocksize[k]));
	}
	result = prdWindowFlush(rh, &win, routineName);

cleanup:
	free(win.win);
	free(J_k);
	free(q_abs);
	free(W);
	free(qp);
	free(qu);
	free(gw);
	return result;
}

/* --- The work function of parallelFor: a failure is recorded (Fail)
	   and stops the loop --                            -------------- */

static void prdIntegrate(RHState *rh, int kbeg, int kend, int nt, void *arg)
{
	(void) prdIntegrateRange(rh, kbeg, kend, nt, arg);
}
/* ------- end ---------------------------- prdIntegrate.c ---------- */

rh_status PRDScatter(RHState *rh, AtomicLine *PRDline)
{
	const char   routineName[] = "PRDScatter";
	register int la, k, kxrd;
	char         label[MAX_MESSAGE_LENGTH];
	int          Nsubordinate;
	off_t        offset;
	double       cDop, *adamp, *Yj, *denom;
	Atom        *atom;
	AtomicLine **XRD;
	PRDjob       job;
	rh_status    result = RH_OK;

	/* --- Routine for angle-averaged and hybrid angle-dependent
		   approximation. In the latter Jgas is used instead of J -- -- */

	atom = PRDline->atom;
	RH_ASSERT(rh, PRDline->PRD, routineName, "Line %d -> %d of %2s is not a PRD line",
			  PRDline->j, PRDline->i, atom->ID);
	getCPU(rh, 3, TIME_START, NULL);

	/* --- Set XRD line array --                         -------------- */

	if (rh->input.XRD)
		Nsubordinate = 1 + PRDline->Nxrd;
	else
		Nsubordinate = 1;

	XRD    = (AtomicLine **) mallocOrAbort(Nsubordinate * sizeof(AtomicLine *));
	XRD[0] = PRDline;
	if (rh->input.XRD) {
		for (kxrd = 0; kxrd < PRDline->Nxrd; kxrd++)
			XRD[kxrd + 1] = PRDline->xrd[kxrd];
	}

	/* --- Initialize the emission profile ratio rho --  -------------- */

	for (la = 0; la < PRDline->Nlambda; la++) {
		for (k = 0; k < rh->atmos.Nspace; k++) {
			PRDline->rho_prd[la][k] = 1.0;
		}
	}
	Yj    = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	denom = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	adamp = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	cDop  = (NM_TO_M * PRDline->lambda0) / (4.0 * PI);

	prdInflow(rh, atom, PRDline, XRD, Nsubordinate, Yj);
	for (k = 0; k < rh->atmos.Nspace; k++) {
		adamp[k] = (PRDline->Grad + PRDline->Qelast[k]) * cDop / atom->vbroad[k];

		/* --- Denominator 1 + sum_x g_x Jbar(x) of rho (see prdInflow) - */

		denom[k] = 1.0;
		for (kxrd = 0; kxrd < Nsubordinate; kxrd++)
			denom[k] += atom->n[XRD[kxrd]->i][k] * XRD[kxrd]->Rij[k] / Yj[k];
	}

	job.atom    = atom;
	job.PRDline = PRDline;
	job.Yj      = Yj;
	job.adamp   = adamp;
	job.kstart  = (off_t *) mallocOrAbort(rh->atmos.Nspace * sizeof(off_t));
	RH_TRY_GOTO(result, cleanup, prdOpenFile(rh, PRDline, routineName, &job.initialize));

	/* --- Sizes of the locations' blocks of weights, for each
		   subordinate line: computed with the weights, kept with the
		   line --                                       -------------- */

	if (job.initialize) {
		free(PRDline->prd_block);
		PRDline->prd_block = (long *) mallocOrAbort((size_t) Nsubordinate * rh->atmos.Nspace * sizeof(long));
	}

	/* --- Loop over subordinate lines --                -------------- */

	offset = 0;
	for (kxrd = 0; kxrd < Nsubordinate; kxrd++) {
		job.XRDline   = XRD[kxrd];
		job.waveratio = job.XRDline->lambda0 / PRDline->lambda0;
		job.blocksize = PRDline->prd_block + (size_t) kxrd * rh->atmos.Nspace;

		/* --- Where each location's block of weights sits in the file -- */

		if (job.initialize)
			RH_TRY_GOTO(result, cleanup,
						parallelFor(rh, rh->atmos.Nspace, rh->input.Nthreads, prdBlockSizes, &job));
		offset = prdOffsets(rh, job.blocksize, job.kstart, offset);

		/* --- Loop over all spatial locations --          -------------- */

		RH_TRY_GOTO(result, cleanup,
					parallelFor(rh, rh->atmos.Nspace, rh->input.Nthreads, prdIntegrate, &job));
	}
	for (la = 0; la < PRDline->Nlambda; la++) {
		for (k = 0; k < rh->atmos.Nspace; k++) PRDline->rho_prd[la][k] /= denom[k];
	}

	snprintf(label, sizeof(label), "Scatter Int %5.1f", PRDline->lambda0);
	getCPU(rh, 3, TIME_POLL, label);

cleanup:
	if (result != RH_OK && job.initialize) prdDropFile(PRDline);
	free(job.kstart);
	free(XRD);
	free(Yj);
	free(denom);
	free(adamp);
	return result;
}
/* ------- end ---------------------------- PRDScatter.c ------------ */

/* ------- begin -------------------------- PRDAngleScatter.c ------- */

/* --- Angle-dependent PRD (P-14). For every emission wavelength,
	   direction and location the integral runs over all absorption
	   directions and, on a fine grid, over absorption wavelength, with
	   the intensity of each direction linearly interpolated from the
	   line's grid. It is folded into weights on the line's grid as in
	   PRDScatter, normalized by the sum of the redistribution weights
	   over all absorption directions, and stored per location: rows in
	   the order emission wavelength, emission ray and direction,
	   absorption ray and direction. The intensities of all absorption
	   directions are read once per call (they used to be read again for
	   every emission wavelength and direction), and the loop over
	   locations is threaded. --                        -------------- */

typedef struct {
	Atom       *atom;
	AtomicLine *PRDline;
	bool_t      initialize;
	int         Ndir;
	double     *g, *denom, *adamp, *sv, **v_los, *Iall;
	long       *blocksize;
	off_t      *kstart;
} PRDAngleJob;

/* --- Emission wavelength in Doppler units for (la, mu, to_obs), and
	   absorption grid of direction d = 2 mup + to_obs_p, at location k */

static double prdAngleEmit(PRDAngleJob *job, int k, int la, int mu, bool_t to_obs)
{
	AtomicLine *line = job->PRDline;
	double      q    = (line->lambda[la] - line->lambda0) * CLIGHT /
			   (line->lambda0 * job->atom->vbroad[k]);

	return (to_obs) ? q + job->v_los[mu][k] : q - job->v_los[mu][k];
}

static void prdAngleAbsGrid(PRDAngleJob *job, int k, int d, double *v_abs)
{
	AtomicLine *line = job->PRDline;
	int         lap, mup = d / 2;
	bool_t      to_obs_p = d % 2;

	for (lap = 0; lap < line->Nlambda; lap++) {
		v_abs[lap] = (line->lambda[lap] - line->lambda0) * CLIGHT /
					 (line->lambda0 * job->atom->vbroad[k]);
		v_abs[lap] += (to_obs_p) ? job->v_los[mup][k] : -job->v_los[mup][k];
	}
}

static void prdAngleBlockSizes(RHState *rh, int kbeg, int kend, int nt, void *arg)
{
	PRDAngleJob *job  = (PRDAngleJob *) arg;
	AtomicLine  *line = job->PRDline;
	int          k, la, mu, to_obs, d, Np, m0;
	long         n;
	double       v_emit, v0, vN, *v_abs;

	(void) nt;

	v_abs = (double *) mallocOrAbort((size_t) job->Ndir * line->Nlambda * sizeof(double));
	for (k = kbeg; k < kend; k++) {
		for (d = 0; d < job->Ndir; d++) prdAngleAbsGrid(job, k, d, v_abs + d * line->Nlambda);
		n = 0;
		for (la = 0; la < line->Nlambda; la++) {
			for (mu = 0; mu < rh->atmos.Nrays; mu++) {
				for (to_obs = 0; to_obs <= 1; to_obs++) {
					v_emit = prdAngleEmit(job, k, la, mu, to_obs);
					Np     = prdGrid(v_emit, 1.0, &v0, &vN);
					for (d = 0; d < job->Ndir; d++)
						n += prdFoldRange(line->Nlambda, v_abs + d * line->Nlambda, v0, Np, FALSE,
										  &m0);
				}
			}
		}
		job->blocksize[k] = n;
	}
	free(v_abs);
}

static rh_status prdAngleIntegrateRange(RHState *rh, int kbeg, int kend, int nt, void *arg)
{
	const char routineName[] = "PRDAngleScatter";

	PRDAngleJob *job  = (PRDAngleJob *) arg;
	AtomicLine  *line = job->PRDline;
	int          Nl   = line->Nlambda, k, la, mu, mup, to_obs, d, lap, lamu, Np, m, *m0, *len,
		Nrow;
	long   pos, start;
	double v_emit, v0, vN, *vp = NULL, *v_abs, *rii = NULL, *W, norm, wmup, s, wsum,
						   *I_d, phi;
	prdw_t   *wbuf;
	PRDWindow win;
	rh_status result = RH_OK;

	(void) nt;

	v_abs = (double *) mallocOrAbort((size_t) job->Ndir * Nl * sizeof(double));
	W     = (double *) mallocOrAbort((size_t) job->Ndir * Nl * sizeof(double));
	m0    = (int *) mallocOrAbort(job->Ndir * sizeof(int));
	len   = (int *) mallocOrAbort(job->Ndir * sizeof(int));
	prdWindowOpen(&win, job->PRDline, job->initialize, job->blocksize, job->kstart, kbeg, kend);

	for (k = kbeg; k < kend; k++) {
		RH_TRY_GOTO(result, cleanup, prdWindowBlock(rh, &win, k, routineName, &wbuf));
		for (d = 0; d < job->Ndir; d++) prdAngleAbsGrid(job, k, d, v_abs + d * Nl);

		pos = 0;
		for (la = 0; la < Nl; la++) {
			for (mu = 0; mu < rh->atmos.Nrays; mu++) {
				for (to_obs = 0; to_obs <= 1; to_obs++) {
					lamu   = 2 * (rh->atmos.Nrays * la + mu) + to_obs;
					v_emit = prdAngleEmit(job, k, la, mu, to_obs);
					Np     = prdGrid(v_emit, 1.0, &v0, &vN);

					Nrow = 0;
					for (d = 0; d < job->Ndir; d++) {
						len[d] = prdFoldRange(Nl, v_abs + d * Nl, v0, Np, FALSE, &m0[d]);
						Nrow += len[d];
					}
					if (pos + Nrow > job->blocksize[k])
						RH_TRY_GOTO(result, cleanup,
									prdBlockError(rh, routineName, k, pos + Nrow, job->blocksize[k]));
					if (job->initialize) {

						/* --- RII(x, mu, x', mu') / phi(x, mu) on the fine
							   grid, with the integration weights of Press et
							   al., Numerical Recipes, p. 107, eq. 4.1.12, and
							   the weight of the absorption direction; folded
							   per direction, normalized over all of them -- */

						phi = (rh->atmos.moving || (line->polarizable && rh->input.StokesMode == FULL_STOKES))
								  ? line->phi[lamu][k]
								  : line->phi[la][k];
						rii = (double *) reallocOrAbort(rii, Np * sizeof(double));
						vp  = (double *) reallocOrAbort(vp, Np * sizeof(double));
						prdFineGrid(v0, Np, FALSE, vp);
						norm = 0.0;
						for (m = 0; m < Nrow; m++) W[m] = 0.0;
						for (d = 0, start = 0; d < job->Ndir; start += len[d], d++) {
							mup  = d / 2;
							wmup = 0.5 * rh->atmos.wmu[mup];
							for (lap = 0; lap < Np; lap++) {
								rii[lap] = RII(rh, v_emit, vp[lap], job->adamp[k], mu, mup) *
										   (job->sv[k] / phi) * prdQuadWeight(lap, Np) * wmup;
								norm += rii[lap];
							}
							prdFold(rh, Nl, v_abs + d * Nl, Np, vp, rii, m0[d], len[d], W + start);
						}
						for (m = 0; m < Nrow; m++) wbuf[pos + m] = (prdw_t) (W[m] / norm);
					}

					/* --- Sum over absorption directions of the dot
						   products with the intensity --     -------------- */

					s    = 0.0;
					wsum = 0.0;
					for (d = 0, start = pos; d < job->Ndir; start += len[d], d++) {
						I_d = job->Iall + ((size_t) k * job->Ndir + d) * Nl;
						for (m = 0; m < len[d]; m++) {
							s += (double) wbuf[start + m] * I_d[m0[d] + m];
							wsum += (double) wbuf[start + m];
						}
					}
					s /= wsum; /* the row sums to 1 in double (N-56) */
					line->rho_prd[lamu][k] = (1.0 + job->g[k] * s) / job->denom[k];
					pos += Nrow;
				}
			}
		}
		if (pos != job->blocksize[k])
			RH_TRY_GOTO(result, cleanup, prdBlockError(rh, routineName, k, pos, job->blocksize[k]));
	}
	result = prdWindowFlush(rh, &win, routineName);

cleanup:
	free(win.win);
	free(v_abs);
	free(W);
	free(m0);
	free(len);
	free(vp);
	free(rii);
	return result;
}

/* --- The work function of parallelFor: a failure is recorded (Fail)
	   and stops the loop --                            -------------- */

static void prdAngleIntegrate(RHState *rh, int kbeg, int kend, int nt, void *arg)
{
	(void) prdAngleIntegrateRange(rh, kbeg, kend, nt, arg);
}

rh_status PRDAngleScatter(RHState *rh, AtomicLine *PRDline)
{
	const char   routineName[] = "PRDAngleScatter";
	register int k, mu, lap, d;
	char         label[MAX_MESSAGE_LENGTH];
	double       cDop, Jbar, *Yj, *Irow;
	Atom        *atom;
	PRDAngleJob  job;
	rh_status    result = RH_OK;

	/* --- Calculate the angle-dependent scattering integral when
		   angle-dependent scattering is requested. --   -------------- */

	atom = PRDline->atom;

	RH_ASSERT(rh, PRDline->PRD, routineName, "Line %d -> %d of %2s is not a PRD line",
			  PRDline->j, PRDline->i, atom->ID);
	getCPU(rh, 3, TIME_START, NULL);

	cDop = (NM_TO_M * PRDline->lambda0) / (4.0 * PI);

	job.atom = atom;
	RH_TRY(prdOpenFile(rh, PRDline, routineName, &job.initialize));
	job.PRDline = PRDline;
	job.Ndir    = 2 * rh->atmos.Nrays;
	job.g       = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	job.denom   = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	job.adamp   = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	job.sv      = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	job.v_los   = matrix_double(rh, rh->atmos.Nrays, rh->atmos.Nspace);
	job.kstart  = (off_t *) mallocOrAbort(rh->atmos.Nspace * sizeof(off_t));
	if (job.initialize) {
		free(PRDline->prd_block);
		PRDline->prd_block = (long *) mallocOrAbort(rh->atmos.Nspace * sizeof(long));
	}
	job.blocksize = PRDline->prd_block;
	Yj            = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));

	/* --- Inflow Yj into the upper level and the weight g of the
		   scattering integral (see prdInflow) --         -------------- */

	prdInflow(rh, atom, PRDline, &PRDline, 1, Yj);
	for (k = 0; k < rh->atmos.Nspace; k++) {
		job.g[k]     = atom->n[PRDline->i][k] * PRDline->Bij / Yj[k];
		Jbar         = PRDline->Rij[k] / PRDline->Bij;
		job.denom[k] = 1.0 + job.g[k] * Jbar;
		job.sv[k]    = 1.0 / (SQRTPI * atom->vbroad[k]);
		job.adamp[k] = (PRDline->Grad + PRDline->Qelast[k]) * cDop / atom->vbroad[k];
	}

	/* --- Line-of-sight velocity in Doppler units --    -------------- */

	for (mu = 0; mu < rh->atmos.Nrays; mu++) {
		for (k = 0; k < rh->atmos.Nspace; k++)
			job.v_los[mu][k] = vproject(rh, k, mu) / atom->vbroad[k];
	}

	/* --- Specific intensity of every absorption direction d = 2 mup
		   + to_obs_p and wavelength of the line, read once, and stored
		   per location (Iall[k][d][lap]), so that each location's
		   intensities are contiguous. That is 2 Nrays times what the
		   old code held (one direction at a time), the same order as
		   the line's per-ray profiles --                 -------------- */

	job.Iall = (double *) malloc((size_t) job.Ndir * PRDline->Nlambda * rh->atmos.Nspace *
								 sizeof(double));
	Irow     = (double *) malloc(rh->atmos.Nspace * sizeof(double));
	if (job.Iall == NULL || Irow == NULL) {
		result = Fail(rh, RH_E_NOMEM, routineName, "Out of memory");
		free(Irow);
		goto cleanup;
	}
	for (d = 0; d < job.Ndir; d++)
		for (lap = 0; lap < PRDline->Nlambda; lap++) {
			readImu(rh, PRDline->Nblue + lap, d / 2, d % 2, Irow);
			for (k = 0; k < rh->atmos.Nspace; k++)
				job.Iall[((size_t) k * job.Ndir + d) * PRDline->Nlambda + lap] = Irow[k];
		}
	free(Irow);

	/* --- RII is called from several threads; a geometry that does not
		   implement it stops here, in this thread --     -------------- */

	if (job.initialize) (void) RII(rh, 0.0, 0.0, job.adamp[0], 0, 0);

	if (job.initialize)
		RH_TRY_GOTO(result, cleanup,
					parallelFor(rh, rh->atmos.Nspace, rh->input.Nthreads, prdAngleBlockSizes, &job));
	prdOffsets(rh, job.blocksize, job.kstart, 0);
	RH_TRY_GOTO(result, cleanup,
				parallelFor(rh, rh->atmos.Nspace, rh->input.Nthreads, prdAngleIntegrate, &job));

	snprintf(label, sizeof(label), "Scatter Int %5.1f", PRDline->lambda0);
	getCPU(rh, 3, TIME_POLL, label);

cleanup:
	if (result != RH_OK && job.initialize) prdDropFile(PRDline);
	free(job.Iall);
	free(job.g);
	free(job.denom);
	free(job.adamp);
	free(job.sv);
	freeMatrix(rh, (void **) job.v_los);
	free(job.kstart);
	free(Yj);
	return result;
}
/* ------- end ---------------------------- PRDAngleScatter.c ------- */
