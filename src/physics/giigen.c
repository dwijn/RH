/* ------- file: -------------------------- giigen.c ----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Jul  1 05:11:13 2010 --

	   --------------------------                      ----------RH-- */

/* --- Gouttebroze's fast approximation GII of the angle-averaged general
	   redistribution function PII for ordinary and cross redistribution
	   of lines with sharp lower and broadened upper levels. (The exact
	   PII by Gauss quadrature, which nothing called, was removed: Q-47.)

  See: D. G. Hummer 1962, MNRAS, 125, 21-37.
	   I. Huben\'y, 1982, JQSRT 27, 593.

	   Parameters:

		 adamp     --  Damping parameter line.
		 waveRatio --  Ratio of wavelengths of lines in which absorption and
					   emission occurs respectively (set to 1.0 for ordinary
					   redistribution).
		 q_abs     --  Frequency absorbed photon (Doppler units)
		 q_emit    --  Frequency emitted photon (Doppler units) -- --- */

#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "atom.h"
#include "constant.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- GII --------------------- */
/*
 * Gouttebroze's fast approximation for
 *  GII(q_abs, q_emit) = PII(q_abs, q_emit) / phi(q_emit)

 * See: P. Gouttebroze, 1986, A&A 160, 195
 *      H. Uitenbroek,  1989, A&A, 216, 310-314 (cross redistribution)
 */

#define GZERO(x) (1.0 / (fabs(x) + sqrt((x) * (x) + 1.273239545)))

double GII(double adamp, double waveratio, double q_emit, double q_abs)
{
	double gii, pcore, aq_emit, umin, epsilon, giiwing, u1, phicore,
		phiwing;

	/* --- Symmetrize with respect to emission frequency --   --------- */

	if (q_emit < 0.0) {
		q_emit = -q_emit;
		q_abs  = -q_abs;
	}
	pcore = 0.0;
	gii   = 0.0;

	/* --- Core region --                                     --------- */

	if (q_emit < PRD_QWING) {
		if ((q_abs < -PRD_QWING) ||
			(q_abs > q_emit + waveratio * PRD_QSPREAD)) return gii;
		if (fabs(q_abs) <= q_emit)
			gii = GZERO(q_emit);
		else
			gii = exp(SQ(q_emit) - SQ(q_abs)) * GZERO(q_abs);

		if (q_emit >= PRD_QCORE) {
			phicore = exp(-SQ(q_emit));
			phiwing = adamp / (SQRTPI * (SQ(adamp) + SQ(q_emit)));
			pcore   = phicore / (phicore + phiwing);
		}
	}
	/* --- Wing region --                                     --------- */

	if (q_emit >= PRD_QCORE) {
		aq_emit = waveratio * q_emit;
		if (q_emit >= PRD_QWING) {
			if (fabs(q_abs - aq_emit) > waveratio * PRD_QSPREAD) return gii;
			pcore = 0.0;
		}
		umin    = fabs((q_abs - aq_emit) / (1.0 + waveratio));
		giiwing = (1.0 + waveratio) * (1.0 - 2.0 * umin * GZERO(umin)) *
				  exp(-SQ(umin));

		if (waveratio == 1.0) {
			epsilon = q_abs / aq_emit;
			giiwing *= (2.75 - (2.5 - 0.75 * epsilon) * epsilon);
		} else {
			u1 = fabs((q_abs - aq_emit) / (waveratio - 1.0));
			giiwing -= fabs(1.0 - waveratio) * (1.0 - 2.0 * u1 * GZERO(u1)) *
					   exp(-SQ(u1));
		}
		/* --- Linear combination of core- and wing contributions ------- */

		giiwing = giiwing / (2.0 * waveratio * SQRTPI);
		gii     = pcore * gii + (1.0 - pcore) * giiwing;
	}
	return gii;
}
/* ------- end ---------------------------- GII --------------------- */
