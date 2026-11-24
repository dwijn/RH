/* ------- file: -------------------------- mpe.h -------------------

	   Minimal polynomial extrapolation of a vector sequence, for the
	   error estimate of ERROR_LIMIT (N-3), see mpe.c.

	   --------------------------                      ----------RH-- */

#ifndef __MPE_H__
#define __MPE_H__

/* --- Needs rh.h (bool_t) --                      -------------- */

bool_t mpeExtrapolate(int K, long N, double **D);

#endif /* !__MPE_H__ */

/* ------- end ---------------------------- mpe.h ------------------- */
