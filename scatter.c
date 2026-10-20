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

#define TENSION           8.0
#define PRD_FILE_TEMPLATE "PRD_%s_%d-%d.dat"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum   spectrum;
extern InputData  input;
extern char       messageStr[];

/* ------- begin -------------------------- prdInflow.c -------------- */

static void prdInflow(Atom *atom, AtomicLine *PRDline, AtomicLine **XRD,
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

	for (k = 0; k < atmos.Nspace; k++)
		Yj[k] = atom->n[PRDline->j][k] * PRDline->Qelast[k];
	for (ip = 0; ip < atom->Nlevel; ip++) {
		if (ip == PRDline->j) continue;
		for (k = 0; k < atmos.Nspace; k++)
			Yj[k] += atom->n[ip][k] * atom->C[PRDline->j * atom->Nlevel + ip][k];
	}
	for (kr = 0; kr < atom->Nline; kr++) {
		line = &atom->line[kr];
		if (line->j == PRDline->j) {
			subordinate = FALSE;
			for (kxrd = 0; kxrd < Nsubordinate; kxrd++)
				if (XRD[kxrd] == line) subordinate = TRUE;
			if (!subordinate)
				for (k = 0; k < atmos.Nspace; k++)
					Yj[k] += atom->n[line->i][k] * line->Rij[k];
		}
		if (line->i == PRDline->j)
			for (k = 0; k < atmos.Nspace; k++)
				Yj[k] += atom->n[line->j][k] * line->Rji[k];
	}
	for (kr = 0; kr < atom->Ncont; kr++) {
		continuum = &atom->continuum[kr];
		if (continuum->j == PRDline->j)
			for (k = 0; k < atmos.Nspace; k++)
				Yj[k] += atom->n[continuum->i][k] * continuum->Rij[k];
		if (continuum->i == PRDline->j)
			for (k = 0; k < atmos.Nspace; k++)
				Yj[k] += atom->n[continuum->j][k] * continuum->Rji[k];
	}
}
/* ------- end ---------------------------- prdInflow.c -------------- */

/* ------- begin -------------------------- PRDScatter.c ------------ */

/* --- The loop over spatial locations is parallelized: the scattering
	   integral at one location uses only that location's data and
	   writes only rho_prd[.][k] for that k. The redistribution weights
	   live in a scratch file with one block per location (all emission
	   wavelengths of the line, in order), so each worker reads or
	   writes the block of its own locations at a known offset with
	   pread/pwrite instead of sharing a file position. The file layout
	   is the same as before: subordinate line, then location, then
	   emission wavelength.  --                        -------------- */

typedef struct {
	Atom              *atom;
	AtomicLine        *PRDline, *XRDline;
	enum Interpolation representation;
	bool_t             initialize;
	int                fd;
	double             waveratio, *Yj, *adamp;
	long              *blocksize; /* number of doubles in each block      */
	off_t             *kstart;    /* byte offset of each location's block */
} PRDjob;

/* ------- begin -------------------------- prdGrid.c --------------- */

static int prdGrid(double q_emit, double waveratio, double *q0, double *qN)
{
	/* --- Establish integration limits over absorption wavelength,
		   using only regions where the redistribution function is
		   non-zero, and the number of points of the fine grid.
		   (See also function GII.) --                   -------------- */

	if (fabs(q_emit) < PRD_QCORE) {
		*q0 = -PRD_QWING;
		*qN = PRD_QWING;
	} else {
		if (fabs(q_emit) < PRD_QWING) {
			if (q_emit > 0.0) {
				*q0 = -PRD_QWING;
				*qN = waveratio * (q_emit + PRD_QSPREAD);
			} else {
				*q0 = waveratio * (q_emit - PRD_QSPREAD);
				*qN = PRD_QWING;
			}
		} else {
			*q0 = waveratio * (q_emit - PRD_QSPREAD);
			*qN = waveratio * (q_emit + PRD_QSPREAD);
		}
	}
	return (int) ((*qN - *q0) / PRD_DQ) + 1;
}
/* ------- end ---------------------------- prdGrid.c --------------- */

/* ------- begin -------------------------- prdFileIO.c ------------- */

static bool_t prdFileIO(int fd, double *buf, size_t Ndouble, off_t offset,
						bool_t write_it)
{
	char   *p    = (char *) buf;
	size_t  left = Ndouble * sizeof(double);
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

/* --- The weights are read and written in windows of consecutive
	   locations, several MB at a time, instead of one small block per
	   location. The blocks follow each other in the file, and large
	   requests are what the disk and the OS cache handle well: with
	   small ones, runs whose weight files do not stay in the cache
	   became limited by disk latency. --              -------------- */

#define PRD_IO_WINDOW (2 * 1024 * 1024)

static int prdWindowEnd(PRDjob *job, int k, int kend)
{
	size_t bytes = 0;
	int    ke    = k;

	while (ke < kend &&
		   (ke == k ||
			bytes + job->blocksize[ke] * sizeof(double) <= PRD_IO_WINDOW)) {
		bytes += job->blocksize[ke] * sizeof(double);
		ke++;
	}
	return ke;
}

static size_t prdWindowSize(PRDjob *job, int k, int ke)
{
	return (job->kstart[ke - 1] - job->kstart[k]) / sizeof(double) +
		   job->blocksize[ke - 1];
}

/* ------- begin -------------------------- prdBlockSizes.c --------- */

static void prdBlockSizes(int kbeg, int kend, int nt, void *arg)
{
	PRDjob     *job     = (PRDjob *) arg;
	AtomicLine *PRDline = job->PRDline;
	int         la, k;
	long        n;
	double      q_emit, q0, qN;

	(void) nt;

	for (k = kbeg; k < kend; k++) {
		n = 0;
		for (la = 0; la < PRDline->Nlambda; la++) {
			q_emit = (PRDline->lambda[la] - PRDline->lambda0) * CLIGHT /
					 (PRDline->lambda0 * job->atom->vbroad[k]);
			n += prdGrid(q_emit, job->waveratio, &q0, &qN);
		}
		job->blocksize[k] = n;
	}
}
/* ------- end ---------------------------- prdBlockSizes.c --------- */

/* ------- begin -------------------------- prdIntegrate.c ---------- */

static void prdIntegrate(int kbeg, int kend, int nt, void *arg)
{
	const char routineName[] = "PRDScatter";

	PRDjob      *job     = (PRDjob *) arg;
	Atom        *atom    = job->atom;
	AtomicLine  *PRDline = job->PRDline, *XRDline = job->XRDline;
	register int la, k, lap;
	bool_t       hunt;
	int          Np;
	long         pos;
	int          kwin_start = kbeg, kwin_end = kbeg;
	size_t       Nwin = 0, wincap = 0;
	double       q_emit, q0, qN, *q_abs = NULL, *qp = NULL, *wq = NULL,
						   *qpp = NULL, *gii, *win = NULL, *wbuf, *J = NULL, *J_k = NULL,
						   scatInt, gnorm, g;

	(void) nt;

	/* --- Every worker has its own scratch arrays --    -------------- */

	J_k   = (double *) malloc(XRDline->Nlambda * sizeof(double));
	q_abs = (double *) malloc(XRDline->Nlambda * sizeof(double));

	for (k = kbeg; k < kend; k++) {

		/* --- Start a new window of weights when the previous one is
			   used up (writing the finished one out first when the weights
			   are being computed for the first time) --   -------------- */

		if (k == kwin_end) {
			if (job->initialize && k > kbeg) {
				if (!prdFileIO(job->fd, win, Nwin, job->kstart[kwin_start], TRUE)) {
					sprintf(messageStr,
							"Unable to write proper number of redistribution weights\n"
							" Line %d -> %d, locations %d to %d",
							XRDline->j, XRDline->i, kwin_start, kwin_end - 1);
					Error(ERROR_LEVEL_2, routineName, messageStr);
				}
			}
			kwin_start = k;
			kwin_end   = prdWindowEnd(job, k, kend);
			Nwin       = prdWindowSize(job, kwin_start, kwin_end);

			/* --- One buffer per worker, sized for a full window and kept:
				   growing it window by window, with sizes that differ a
				   little each time, made the allocator hold on to several
				   windows' worth of freed blocks. --        -------------- */

			if (Nwin > wincap) {
				wincap = (Nwin > PRD_IO_WINDOW / sizeof(double)) ? Nwin : PRD_IO_WINDOW / sizeof(double);
				win    = (double *) realloc(win, wincap * sizeof(double));
			}
			if (!job->initialize) {
				if (!prdFileIO(job->fd, win, Nwin, job->kstart[kwin_start], FALSE)) {
					sprintf(messageStr,
							"Unable to read proper number of redistribution weights\n"
							" Line %d -> %d, locations %d to %d",
							XRDline->j, XRDline->i, kwin_start, kwin_end - 1);
					Error(ERROR_LEVEL_2, routineName, messageStr);
				}
			}
		}
		wbuf = win + (job->kstart[k] - job->kstart[kwin_start]) / sizeof(double);

		g = atom->n[XRDline->i][k] * XRDline->Bij / job->Yj[k];

		/* --- Get local mean intensity and wavelength in Doppler units - */

		if (input.PRD_angle_dep == PRD_ANGLE_APPROX) {
			for (la = 0; la < XRDline->Nlambda; la++) {
				J_k[la]   = spectrum.Jgas[XRDline->Nblue + la][k];
				q_abs[la] = (XRDline->lambda[la] - XRDline->lambda0) * CLIGHT /
							(XRDline->lambda0 * atom->vbroad[k]);
			}
		} else if (input.PRD_angle_dep == PRD_ANGLE_AVER) {
			for (la = 0; la < XRDline->Nlambda; la++) {
				J_k[la]   = spectrum.J[XRDline->Nblue + la][k];
				q_abs[la] = (XRDline->lambda[la] - XRDline->lambda0) * CLIGHT /
							(XRDline->lambda0 * atom->vbroad[k]);
			}
		}
		switch (job->representation) {
		case LINEAR:
			break;
		case SPLINE:
			splineCoef(XRDline->Nlambda, q_abs, J_k);
			break;
		case EXP_SPLINE:
			exp_splineCoef(XRDline->Nlambda, q_abs, J_k, TENSION);
			break;
		}

		/* --- Outer wavelength loop over emission wavelengths -- ----- */

		pos = 0;
		for (la = 0; la < PRDline->Nlambda; la++) {
			q_emit = (PRDline->lambda[la] - PRDline->lambda0) * CLIGHT /
					 (PRDline->lambda0 * atom->vbroad[k]);
			Np = prdGrid(q_emit, job->waveratio, &q0, &qN);

			qp = (double *) realloc(qp, Np * sizeof(double));
			for (lap = 1, qp[0] = q0; lap < Np; lap++)
				qp[lap] = qp[lap - 1] + PRD_DQ;

			/* --- Fold interpolation of J for symmetric lines. -- ------ */

			if (XRDline->symmetric) {
				qpp = (double *) realloc(qpp, Np * sizeof(double));
				for (lap = 0; lap < Np; lap++) qpp[lap] = fabs(qp[lap]);
			}

			/* --- Interpolate mean intensity onto fine grid. Choose linear,
			   spline or exponential spline interpolation -- -------- */

			J = (double *) realloc(J, Np * sizeof(double));
			switch (job->representation) {
			case LINEAR:
				Linear(XRDline->Nlambda, q_abs, J_k, Np,
					   (XRDline->symmetric) ? qpp : qp, J, hunt = TRUE);
				break;
			case SPLINE:
				splineEval(Np, (XRDline->symmetric) ? qpp : qp, J, hunt = TRUE);
				break;
			case EXP_SPLINE:
				exp_splineEval(Np, (XRDline->symmetric) ? qpp : qp, J, hunt = TRUE);
				break;
			}

			/* --- The redistribution weights: computed and kept in this
				   location's block on the first call, read from it on later
				   calls --                                  -------------- */

			gii = wbuf + pos;
			if (job->initialize) {

				/* --- Integration weights (See: Press et al. Numerical Recipes,
					   p. 107, eq. 4.1.12) --                  -------------- */

				wq    = (double *) realloc(wq, Np * sizeof(double));
				wq[0] = 5.0 / 12.0 * PRD_DQ;
				wq[1] = 13.0 / 12.0 * PRD_DQ;
				for (lap = 2; lap < Np - 2; lap++) wq[lap] = PRD_DQ;
				wq[Np - 1] = 5.0 / 12.0 * PRD_DQ;
				wq[Np - 2] = 13.0 / 12.0 * PRD_DQ;

				for (lap = 0; lap < Np; lap++)
					gii[lap] = GII(job->adamp[k], job->waveratio, q_emit, qp[lap]) *
							   wq[lap];
			}

			/* --- Inner wavelength loop doing actual wavelength integration
				   over absorption wavelengths --            -------------- */

			gnorm   = 0.0;
			scatInt = 0.0;
			for (lap = 0; lap < Np; lap++) {
				gnorm += gii[lap];
				scatInt += J[lap] * gii[lap];
			}
			PRDline->rho_prd[la][k] += g * scatInt / gnorm;
			pos += Np;
		}
	}
	if (job->initialize && kend > kbeg) {
		if (!prdFileIO(job->fd, win, Nwin, job->kstart[kwin_start], TRUE)) {
			sprintf(messageStr,
					"Unable to write proper number of redistribution weights\n"
					" Line %d -> %d, locations %d to %d",
					XRDline->j, XRDline->i, kwin_start, kend - 1);
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
	}
	free(J_k);
	free(q_abs);
	free(qp);
	free(qpp);
	free(wq);
	free(J);
	free(win);
}
/* ------- end ---------------------------- prdIntegrate.c ---------- */

void PRDScatter(AtomicLine *PRDline, enum Interpolation representation)
{
	const char   routineName[] = "PRDScatter";
	register int la, k, kxrd;
	char         filename[MAX_LINE_SIZE];
	bool_t       initialize;
	int          Nsubordinate, nthreads;
	off_t        offset;
	double       cDop, *adamp, *Yj, *denom;
	Atom        *atom;
	AtomicLine **XRD;
	PRDjob       job;

	/* --- Routine for angle-averaged and hybrid angle-dependent
		   approximation. In the latter Jgas is used instead of J -- -- */

	atom = PRDline->atom;
	if (!PRDline->PRD) {
		sprintf(messageStr, "Line %d -> %d of %2s is not a PRD line",
				PRDline->j, PRDline->i, atom->ID);
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	getCPU(3, TIME_START, NULL);

	/* --- Open temporary file for storage of redistribution weights when
		   called for the first time  --                 -------------- */

	initialize = FALSE;
	if (PRDline->fp_GII == NULL) {
		sprintf(filename,
				(atom->ID[1] == ' ') ? "PRD_%.1s_%d-%d.dat" : "PRD_%s_%d-%d.dat",
				atom->ID, PRDline->j, PRDline->i);
		initialize = TRUE;
		if ((PRDline->fp_GII = fopen(filename, "w+")) == NULL) {
			sprintf(messageStr, "Unable to open temporary file %s", filename);
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
	}

	/* --- Set XRD line array --                         -------------- */

	if (input.XRD)
		Nsubordinate = 1 + PRDline->Nxrd;
	else
		Nsubordinate = 1;

	XRD    = (AtomicLine **) malloc(Nsubordinate * sizeof(AtomicLine *));
	XRD[0] = PRDline;
	if (input.XRD) {
		for (kxrd = 0; kxrd < PRDline->Nxrd; kxrd++)
			XRD[kxrd + 1] = PRDline->xrd[kxrd];
	}

	/* --- Initialize the emission profile ratio rho --  -------------- */

	for (la = 0; la < PRDline->Nlambda; la++) {
		for (k = 0; k < atmos.Nspace; k++) {
			PRDline->rho_prd[la][k] = 1.0;
		}
	}
	Yj    = (double *) malloc(atmos.Nspace * sizeof(double));
	denom = (double *) malloc(atmos.Nspace * sizeof(double));
	adamp = (double *) malloc(atmos.Nspace * sizeof(double));
	cDop  = (NM_TO_M * PRDline->lambda0) / (4.0 * PI);

	prdInflow(atom, PRDline, XRD, Nsubordinate, Yj);
	for (k = 0; k < atmos.Nspace; k++) {
		adamp[k] = (PRDline->Grad + PRDline->Qelast[k]) * cDop / atom->vbroad[k];

		/* --- Denominator 1 + sum_x g_x Jbar(x) of rho (see prdInflow) - */

		denom[k] = 1.0;
		for (kxrd = 0; kxrd < Nsubordinate; kxrd++)
			denom[k] += atom->n[XRD[kxrd]->i][k] * XRD[kxrd]->Rij[k] / Yj[k];
	}

	/* --- The spline interpolators keep their coefficients in static
		   storage, so only the linear representation runs in parallel - */

	nthreads = (representation == LINEAR) ? input.Nthreads : 1;

	job.atom           = atom;
	job.PRDline        = PRDline;
	job.representation = representation;
	job.initialize     = initialize;
	job.fd             = fileno(PRDline->fp_GII);
	job.Yj             = Yj;
	job.adamp          = adamp;
	job.blocksize      = (long *) malloc(atmos.Nspace * sizeof(long));
	job.kstart         = (off_t *) malloc(atmos.Nspace * sizeof(off_t));

	/* --- Loop over subordinate lines --                -------------- */

	offset = 0;
	for (kxrd = 0; kxrd < Nsubordinate; kxrd++) {
		job.XRDline   = XRD[kxrd];
		job.waveratio = job.XRDline->lambda0 / PRDline->lambda0;

		/* --- Where each location's block of weights sits in the file -- */

		parallelFor(atmos.Nspace, nthreads, prdBlockSizes, &job);
		for (k = 0; k < atmos.Nspace; k++) {
			job.kstart[k] = offset;
			offset += job.blocksize[k] * sizeof(double);
		}

		/* --- Loop over all spatial locations --          -------------- */

		parallelFor(atmos.Nspace, nthreads, prdIntegrate, &job);
	}
	for (la = 0; la < PRDline->Nlambda; la++) {
		for (k = 0; k < atmos.Nspace; k++) PRDline->rho_prd[la][k] /= denom[k];
	}

	free(job.blocksize);
	free(job.kstart);
	free(XRD);
	free(Yj);
	free(denom);
	free(adamp);

	sprintf(messageStr, "Scatter Int %5.1f", PRDline->lambda0);
	getCPU(3, TIME_POLL, messageStr);
}
/* ------- end ---------------------------- PRDScatter.c ------------ */

/* ------- begin -------------------------- PRDAngleScatter.c ------- */

void PRDAngleScatter(AtomicLine        *PRDline,
					 enum Interpolation representation)
{
	const char   routineName[] = "PRDAngleScatter";
	register int la, k, lap, mu, mup;

	char    filename[MAX_LINE_SIZE];
	bool_t  hunt, initialize, to_obs, to_obs_p;
	int     Np, Nread, Nwrite, lamu;
	double *v_emit, v0, vN, *v_abs = NULL, *vp = NULL, *wv = NULL,
							*rii = NULL, *adamp, *Jbar, cDop, *RIInorm, *I = NULL,
							**Imup, *Ik, *Yj, *g, *denom, **v_los, *phi_emit, wmup, *sv;
	Atom *atom;

	/* --- Calculate the angle-dependent scattering integral when
		   angle-dependent scattering is requested. --   -------------- */

	atom = PRDline->atom;

	if (!PRDline->PRD) {
		sprintf(messageStr, "Line %d -> %d of %2s is not a PRD line",
				PRDline->j, PRDline->i, atom->ID);
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
	getCPU(3, TIME_START, NULL);

	cDop = (NM_TO_M * PRDline->lambda0) / (4.0 * PI);

	initialize = FALSE;
	if (PRDline->fp_GII == NULL) {
		sprintf(filename,
				(atom->ID[1] == ' ') ? "PRD_%.1s_%d-%d.dat" : "PRD_%s_%d-%d.dat",
				atom->ID, PRDline->j, PRDline->i);

		initialize = TRUE;
		if ((PRDline->fp_GII = fopen(filename, "w+")) == NULL) {
			sprintf(messageStr, "Unable to open temporary file %s", filename);
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
	}
	if (!initialize) rewind(PRDline->fp_GII);

	/* --- Temporary storage space --                    -------------- */

	Imup    = matrix_double(PRDline->Nlambda, atmos.Nspace);
	Ik      = (double *) malloc(PRDline->Nlambda * sizeof(double));
	v_emit  = (double *) malloc(atmos.Nspace * sizeof(double));
	v_abs   = (double *) malloc(PRDline->Nlambda * sizeof(double));
	adamp   = (double *) malloc(atmos.Nspace * sizeof(double));
	v_los   = matrix_double(atmos.Nrays, atmos.Nspace);
	g       = (double *) malloc(atmos.Nspace * sizeof(double));
	Yj      = (double *) malloc(atmos.Nspace * sizeof(double));
	denom   = (double *) malloc(atmos.Nspace * sizeof(double));
	Jbar    = (double *) malloc(atmos.Nspace * sizeof(double));
	RIInorm = (double *) malloc(atmos.Nspace * sizeof(double));
	sv      = (double *) malloc(atmos.Nspace * sizeof(double));

	/* --- Inflow Yj into the upper level and the weight g of the
		   scattering integral (see prdInflow) --         -------------- */

	prdInflow(atom, PRDline, &PRDline, 1, Yj);
	for (k = 0; k < atmos.Nspace; k++)
		g[k] = atom->n[PRDline->i][k] * PRDline->Bij / Yj[k];

	/* --- Store line-of-sight velocity in Doppler units to avoid
		   having to recompute it for every wavelength -- ------------- */

	for (mu = 0; mu < atmos.Nrays; mu++) {
		for (k = 0; k < atmos.Nspace; k++)
			v_los[mu][k] = vproject(k, mu) / atom->vbroad[k];
	}
	/* --- Depth-dependent damping parameter --          -------------- */

	for (k = 0; k < atmos.Nspace; k++) {
		Jbar[k]  = PRDline->Rij[k] / PRDline->Bij;
		denom[k] = 1.0 + g[k] * Jbar[k];
		sv[k]    = 1.0 / (SQRTPI * atom->vbroad[k]);
		adamp[k] =
			(PRDline->Grad + PRDline->Qelast[k]) * cDop / atom->vbroad[k];
	}
	/* --- Outer loop over emission wavelength --        -------------- */

	for (la = 0; la < PRDline->Nlambda; la++) {

		/* --- Loop over emission angle (down and up) --   -------------- */

		for (mu = 0; mu < atmos.Nrays; mu++) {
			for (to_obs = 0; to_obs <= 1; to_obs++) {
				lamu = 2 * (atmos.Nrays * la + mu) + to_obs;

				if (atmos.moving ||
					(PRDline->polarizable && input.StokesMode == FULL_STOKES))
					phi_emit = PRDline->phi[lamu];
				else
					phi_emit = PRDline->phi[la];

				for (k = 0; k < atmos.Nspace; k++) {
					PRDline->rho_prd[lamu][k] = 0.0;
					RIInorm[k]                = 0.0;
				}
				/* --- Depth-dependent emission wavelengths for this
					   direction (in Doppler units). --        -------------- */

				if (to_obs) {
					for (k = 0; k < atmos.Nspace; k++)
						v_emit[k] = (PRDline->lambda[la] - PRDline->lambda0) * CLIGHT /
										(atom->vbroad[k] * PRDline->lambda0) +
									v_los[mu][k];
				} else {
					for (k = 0; k < atmos.Nspace; k++)
						v_emit[k] = (PRDline->lambda[la] - PRDline->lambda0) * CLIGHT /
										(atom->vbroad[k] * PRDline->lambda0) -
									v_los[mu][k];
				}

				/* --- Loop over absorption directions --      -------------- */

				for (mup = 0; mup < atmos.Nrays; mup++) {
					wmup = 0.5 * atmos.wmu[mup];
					for (to_obs_p = 0; to_obs_p <= 1; to_obs_p++) {

						/* --- Read specific intensity in this direction for all
								   wavelengths --                      -------------- */

						for (lap = 0; lap < PRDline->Nlambda; lap++)
							readImu(PRDline->Nblue + lap, mup, to_obs_p, Imup[lap]);

						/* --- Loop over space --                  -------------- */

						for (k = 0; k < atmos.Nspace; k++) {
							for (lap = 0; lap < PRDline->Nlambda; lap++) {
								Ik[lap] = Imup[lap][k];

								/* --- Array of absorption wavelengths in this
									   direction --                    -------------- */

								if (to_obs_p)
									v_abs[lap] =
										(PRDline->lambda[lap] - PRDline->lambda0) * CLIGHT /
											(PRDline->lambda0 * atom->vbroad[k]) +
										v_los[mup][k];
								else
									v_abs[lap] =
										(PRDline->lambda[lap] - PRDline->lambda0) * CLIGHT /
											(PRDline->lambda0 * atom->vbroad[k]) -
										v_los[mup][k];
							}
							/* --- Setup spline interpolation coefficients for the
						   integration over absorption intensity -- ------- */

							switch (representation) {
							case SPLINE:
								splineCoef(PRDline->Nlambda, v_abs, Ik);
								break;
							case EXP_SPLINE:
								exp_splineCoef(PRDline->Nlambda, v_abs, Ik, TENSION);
								break;
							case LINEAR: break;
							}
							/* --- Establish integration limits over absorption
									   wavelength, using only regions where the
									   redistribution function is non-zero. -- -------- */

							if (fabs(v_emit[k]) < PRD_QCORE) {
								v0 = -PRD_QWING;
								vN = PRD_QWING;
							} else {
								if (fabs(v_emit[k]) < PRD_QWING) {
									if (v_emit[k] > 0.0) {
										v0 = -PRD_QWING;
										vN = v_emit[k] + PRD_QSPREAD;
									} else {
										v0 = v_emit[k] - PRD_QSPREAD;
										vN = PRD_QWING;
									}
								} else {
									v0 = v_emit[k] - PRD_QSPREAD;
									vN = v_emit[k] + PRD_QSPREAD;
								}
							}
							Np = (int) ((vN - v0) / PRD_DQ) + 1;
							vp = (double *) realloc(vp, Np * sizeof(double));
							for (lap = 1, vp[0] = v0; lap < Np; lap++)
								vp[lap] = vp[lap - 1] + PRD_DQ;

							/* --- Interpolate specific intensity onto fine grid.
									   Choose linear, spline or exponential spline
									   interpolation --                  -------------- */

							I = (double *) realloc(I, Np * sizeof(double));
							switch (representation) {
							case LINEAR:
								Linear(PRDline->Nlambda, v_abs, Ik, Np, vp, I, hunt = TRUE);
								break;
							case SPLINE:
								splineEval(Np, vp, I, hunt = TRUE);
								break;
							case EXP_SPLINE:
								exp_splineEval(Np, vp, I, hunt = TRUE);
								break;
							}
							/* --- Compute the redistribution weights -- ---------- */

							rii = (double *) realloc(rii, Np * sizeof(double));
							if (initialize) {

								/* --- Integration weights (See: Press et al.
									   Numerical Recipes, p. 107, eq. 4.1.12) -- ---- */

								wv    = (double *) realloc(wv, Np * sizeof(double));
								wv[0] = 5.0 / 12.0 * PRD_DQ;
								wv[1] = 13.0 / 12.0 * PRD_DQ;
								for (lap = 2; lap < Np - 2; lap++) wv[lap] = PRD_DQ;
								wv[Np - 1] = 5.0 / 12.0 * PRD_DQ;
								wv[Np - 2] = 13.0 / 12.0 * PRD_DQ;

								/* --- RII(x, mu, x', mu')/phi(x, mu) -- ------------ */

								for (lap = 0; lap < Np; lap++) {
									rii[lap] = RII(v_emit[k], vp[lap], adamp[k], mu, mup) *
											   (sv[k] / phi_emit[k]) * wv[lap] * wmup;
								}
								if ((Nwrite = fwrite(rii, sizeof(double), Np,
													 PRDline->fp_GII)) != Np) {
									sprintf(messageStr,
											"Unable to write proper number of redistribution weights\n"
											" Wrote %d instead of %d.\n Line %d -> %d, la = %d, k = %d",
											Nwrite, Np, PRDline->j, PRDline->i, la, k);
									Error(ERROR_LEVEL_2, routineName, messageStr);
								}
							} else {
								if ((Nread = fread(rii, sizeof(double), Np,
												   PRDline->fp_GII)) != Np) {
									sprintf(messageStr,
											"Unable to read proper number of redistribution weights\n"
											" Read %d instead of %d.\n Line %d -> %d, la = %d, k = %d",
											Nread, Np, PRDline->j, PRDline->i, la, k);
									Error(ERROR_LEVEL_2, routineName, messageStr);
								}
							}
							/* --- Inner wavelength loop doing actual wavelength
									   integration over absorption wavelengths -- ----- */

							for (lap = 0; lap < Np; lap++) {
								RIInorm[k] += rii[lap];
								PRDline->rho_prd[lamu][k] += I[lap] * rii[lap];
							}
						}
					}
				}

				for (k = 0; k < atmos.Nspace; k++) {
					PRDline->rho_prd[lamu][k] =
						(1.0 + g[k] * PRDline->rho_prd[lamu][k] / RIInorm[k]) / denom[k];
				}
			}
		}
	}
	/* --- Clean temporary variable space --             -------------- */

	freeMatrix((void **) Imup);
	freeMatrix((void **) v_los);

	free(Ik);
	free(v_abs);
	free(v_emit);
	free(adamp);
	free(g);
	free(Yj);
	free(denom);
	free(wv);
	free(I);
	free(rii);
	free(vp);
	free(Jbar);
	free(RIInorm);
	free(sv);

	sprintf(messageStr, "Scatter Int %5.1f", PRDline->lambda0);
	getCPU(3, TIME_POLL, messageStr);
}
/* ------- end ---------------------------- PRDAngleScatter.c ------- */
