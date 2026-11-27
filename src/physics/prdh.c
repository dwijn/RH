/* ------- file: -------------------------- prdh.c ------------------

	   Gas-frame interpolation of the hybrid PRD approximation
	   (PRD_ANGLE_APPROX), see prdh.h: the weights with which the
	   intensity of a wavelength is added to the gas-frame mean
	   intensity Jgas of its neighbours, the weights with which rho is
	   interpolated to a Doppler-shifted wavelength, and their tables.

	   The tables used to be built in initSolution, and the
	   PRDH_LIMIT_MEM branches of the formal solutions and of Opacity
	   computed the same weights with code of their own (P-17).
	   --                                              -------------- */

#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "prdh.h"
#include "state.h"

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- prdhShift.c ------------- */

double prdhShift(RHState *rh, int mu, int to_obs, int k)
{
	/* --- Factor that takes an observer's-frame wavelength of ray
		   (mu, to_obs) to the gas frame at location k --  ---------- */

	double sign = (to_obs) ? 1.0 : -1.0;

	return 1.0 + rh->spectrum.v_los[mu][k] * sign / CLIGHT;
}
/* ------- end ---------------------------- prdhShift.c ------------- */

/* ------- begin -------------------------- firstIndex.c ------------ */

static long firstIndex(const double *lambda, long N, double x, bool_t above)
{
	/* --- First index of the ascending array lambda[0..N-1] with
		   lambda[idx] > x (above = TRUE) or lambda[idx] >= x (above =
		   FALSE), N if there is none --                  -------------- */

	long lo = 0, hi = N, mid;

	while (lo < hi) {
		mid = lo + (hi - lo) / 2;
		if (above ? (lambda[mid] > x) : (lambda[mid] >= x))
			hi = mid;
		else
			lo = mid + 1;
	}
	return lo;
}
/* ------- end ---------------------------- firstIndex.c ------------ */

/* ------- begin -------------------------- prdhBracket.c ----------- */

void prdhBracket(RHState *rh, int nspect, double fac, PRDhBracket *b)
{
	/* --- The grid (spectrum.lambda) is sorted and has no duplicates
		   (SortLambda), so the points between two wavelengths are a
		   contiguous range of indices (P-12).

		   At the ends of the grid a neighbour is the point itself
		   (prv == gas for nspect = 0, gas == nxt for the last point)
		   and the points beyond are taken with weight 1: those at or
		   below gas at the lower end, those above gas at the upper end.
		   The point equal to gas, which exists where the shift is 0,
		   then gets weight 1 once. Before P-17 the table took the points
		   at or above gas at the upper end, so that the last point
		   added its own intensity twice where v_los = 0, and the
		   PRDH_LIMIT_MEM code took the points below gas at the lower
		   end, so that the first point did not add it at all.

		   The ends are told by the index, not by prv == gas (N-42): two
		   grid points within about 1 ulp of each other can shift to the
		   same gas-frame wavelength, and the comparison then took every
		   point below (or above) with weight 1. Such an interval is
		   empty, so its weight is never evaluated --     -------------- */

	const double *lambda = rh->spectrum.lambda;
	long          N      = rh->spectrum.Nspect;

	b->prv = lambda[MAX(nspect - 1, 0)] * fac;
	b->gas = lambda[nspect] * fac;
	b->nxt = lambda[MIN(nspect + 1, N - 1)] * fac;

	b->lo_end = (nspect == 0);
	b->hi_end = (nspect == N - 1);

	b->r[0] = (b->lo_end) ? 0 : firstIndex(lambda, N, b->prv, TRUE);
	b->r[1] = firstIndex(lambda, N, b->gas, TRUE);
	b->r[2] = (b->hi_end) ? N : firstIndex(lambda, N, b->nxt, FALSE);

	if (b->r[1] < b->r[0]) b->r[1] = b->r[0];
	if (b->r[2] < b->r[1]) b->r[2] = b->r[1];
}
/* ------- end ---------------------------- prdhBracket.c ----------- */

/* ------- begin -------------------------- prdhWeight.c ------------ */

double prdhWeight(RHState *rh, const PRDhBracket *b, long idx)
{
	/* --- Weight of grid point idx, b->r[0] <= idx < b->r[2]: linear
		   interpolation between the shifted neighbours, 1 beyond the
		   ends of the grid --                            -------------- */

	const double *lambda = rh->spectrum.lambda;

	if (idx < b->r[1])
		return (b->lo_end) ? 1.0 : (lambda[idx] - b->prv) / (b->gas - b->prv);
	else
		return (b->hi_end) ? 1.0 : 1.0 - (lambda[idx] - b->gas) / (b->nxt - b->gas);
}
/* ------- end ---------------------------- prdhWeight.c ------------ */

/* ------- begin -------------------------- prdhRhoWeights.c -------- */

void prdhRhoWeights(const AtomicLine *line, double lag, int *id0, int *id1, double *frac)
{
	/* --- Linear interpolation of rho, on the line's grid, to the
		   gas-frame wavelength lag; constant beyond the grid -- ---- */

	int ilow;

	if (lag <= line->lambda[0]) {
		*frac = 0.0;
		*id0  = 0;
		*id1  = 1;
	} else if (lag >= line->lambda[line->Nlambda - 1]) {
		*frac = 1.0;
		*id0  = line->Nlambda - 2;
		*id1  = line->Nlambda - 1;
	} else {
		Locate(line->Nlambda, line->lambda, lag, &ilow);
		*frac = (lag - line->lambda[ilow]) / (line->lambda[ilow + 1] - line->lambda[ilow]);
		*id0  = ilow;
		*id1  = ilow + 1;
	}
}
/* ------- end ---------------------------- prdhRhoWeights.c -------- */

/* ------- begin -------------------------- initPRDhTables.c -------- */

void initPRDhTables(RHState *rh)
{
	/* --- Tabulate the weights of prdhRhoWeights for every PRD line,
		   wavelength, ray, direction and location, and those of
		   prdhBracket and prdhWeight for every wavelength, ray,
		   direction and location (spectrum.nc, iprdh and cprdh).
		   With PRDH_LIMIT_MEM nothing is tabulated and the weights are
		   computed when needed --                        -------------- */

	const char  routineName[] = "initPRDhTables";
	int         nact, kr, la, mu, to_obs, k, Nlamu, lamu;
	long        lamuk, Nlamuk, lc, idx, Ncoef;
	double      fac;
	Atom       *atom;
	AtomicLine *line;
	PRDhBracket b;

	if (rh->input.prdh_limit_mem) return;

	/* --- Interpolation of rho --                       -------------- */

	for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
		atom = rh->atmos.activeatoms[nact];
		for (kr = 0; kr < atom->Nline; kr++) {
			line = &atom->line[kr];
			if (!line->PRD) continue;

			Nlamu      = 2 * rh->atmos.Nrays * line->Nlambda;
			line->frac = matrix_double(rh, Nlamu, rh->atmos.Nspace);
			line->id0  = matrix_int(rh, Nlamu, rh->atmos.Nspace);
			line->id1  = matrix_int(rh, Nlamu, rh->atmos.Nspace);

			for (la = 0; la < line->Nlambda; la++) {
				for (mu = 0; mu < rh->atmos.Nrays; mu++) {
					for (to_obs = 0; to_obs <= 1; to_obs++) {
						lamu = 2 * (rh->atmos.Nrays * la + mu) + to_obs;
						for (k = 0; k < rh->atmos.Nspace; k++)
							prdhRhoWeights(line, line->lambda[la] * prdhShift(rh, mu, to_obs, k),
										   &line->id0[lamu][k], &line->id1[lamu][k],
										   &line->frac[lamu][k]);
					}
				}
			}
		}
	}

	/* --- Interpolation of Jgas: first the number of points each
		   (wavelength, ray, direction, location) contributes to, as a
		   running sum in spectrum.nc, then the points and weights -- */

	Nlamuk          = 2L * rh->atmos.Nrays * rh->spectrum.Nspect * rh->atmos.Nspace;
	rh->spectrum.nc = (long *) malloc(Nlamuk * sizeof(long));
	if (rh->spectrum.nc == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to allocate the hybrid PRD interpolation table (%ld points)", Nlamuk);
	}
	Ncoef = 0;
	for (la = 0; la < rh->spectrum.Nspect; la++) {
		for (mu = 0; mu < rh->atmos.Nrays; mu++) {
			for (to_obs = 0; to_obs <= 1; to_obs++) {
				for (k = 0; k < rh->atmos.Nspace; k++) {
					prdhBracket(rh, la, prdhShift(rh, mu, to_obs, k), &b);
					Ncoef += b.r[2] - b.r[0];
					rh->spectrum.nc[prdhIndex(rh, la, mu, to_obs, k)] = Ncoef;
				}
			}
		}
	}
	rh->spectrum.iprdh = (int *) malloc(MAX(Ncoef, 1) * sizeof(int));
	rh->spectrum.cprdh = (double *) malloc(MAX(Ncoef, 1) * sizeof(double));
	if (rh->spectrum.iprdh == NULL || rh->spectrum.cprdh == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to allocate the hybrid PRD interpolation table"
											  " (%ld points, %ld coefficients)",
			  Nlamuk, Ncoef);
	}
	for (la = 0; la < rh->spectrum.Nspect; la++) {
		for (mu = 0; mu < rh->atmos.Nrays; mu++) {
			for (to_obs = 0; to_obs <= 1; to_obs++) {
				for (k = 0; k < rh->atmos.Nspace; k++) {
					fac   = prdhShift(rh, mu, to_obs, k);
					lamuk = prdhIndex(rh, la, mu, to_obs, k);
					lc    = (lamuk == 0) ? 0 : rh->spectrum.nc[lamuk - 1];
					prdhBracket(rh, la, fac, &b);
					for (idx = b.r[0]; idx < b.r[2]; idx++, lc++) {
						rh->spectrum.iprdh[lc] = (int) idx;
						rh->spectrum.cprdh[lc] = prdhWeight(rh, &b, idx);
					}
				}
			}
		}
	}
}
/* ------- end ---------------------------- initPRDhTables.c -------- */
