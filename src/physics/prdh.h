/* ------- file: -------------------------- prdh.h ------------------

	   Gas-frame interpolation of the hybrid PRD approximation
	   (PRD_ANGLE_APPROX): the formal solution of each wavelength adds
	   its intensity, shifted to the gas frame, to the gas-frame mean
	   intensity Jgas of neighbouring wavelengths, and the emission
	   profile ratio rho, known on the line's grid in the gas frame, is
	   interpolated to the Doppler-shifted wavelength of each ray. The
	   weights are either tabulated once (initPRDhTables) or, with
	   PRDH_LIMIT_MEM, computed when needed, with the same arithmetic,
	   so that both give the same results to the bit.
	   --                                              -------------- */

#ifndef PRDH_H
#define PRDH_H

#include "atom.h"

/* --- Grid points that wavelength nspect, shifted to the gas frame by
	   the factor fac, contributes to: the gas-frame wavelengths of its
	   neighbours and of itself are prv, gas and nxt; the points in
	   (prv, gas] are r[0] <= idx < r[1] and those in (gas, nxt) are
	   r[1] <= idx < r[2]. At the ends of the grid (lo_end: nspect is
	   the first point, hi_end: the last), where a neighbour is the
	   point itself, all points below (or above) are taken with weight 1
	   (constant extrapolation) --                     -------------- */

typedef struct {
	long   r[3];
	double prv, gas, nxt;
	bool_t lo_end, hi_end;
} PRDhBracket;

double prdhShift(int mu, int to_obs, int k);
void   prdhBracket(int nspect, double fac, PRDhBracket *b);
double prdhWeight(const PRDhBracket *b, long idx);
void   prdhRhoWeights(const AtomicLine *line, double lag, int *id0, int *id1, double *frac);
void   initPRDhTables(void);

#endif /* !PRDH_H */
