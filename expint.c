/* ------- file: -------------------------- expint.c ----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Dec  8 16:02:51 1999 --

	   --------------------------                      ----------RH-- */

/* --- Evaluates the exponential integral functions E1(x) and E2(x).

  Ref: Handbook of Mathematical Functions,
	   Applied mathematics series 55 (1964) (ed. Abramovitz and Stegun).
	   Num. Recipes, 2nd ed., Sect. 6.3 (continued fraction for E2).
	   --                                              -------------- */

#include <float.h>
#include <math.h>
#include <stdio.h>

#include "rh.h"
#include "error.h"

/* --- Above this E1 and E2 (~ exp(-x)/x) underflow to 0 -- --------- */

#define X_UNDERFLOW 745.0

/* --- Function prototypes --                          -------------- */

double E1(double x);

/* --- Global variables --                             -------------- */

extern char messageStr[];

/* ------- begin -------------------------- E1.c -------------------- */

#define N_53 6
#define N_56 4

double E1(double x)
{
	/* --- Use Eq 5.1.53 and 5.1.56 --                   ------------- */

	static double a53[N_53] = {-0.57721566, 0.99999193, -0.24991055,
							   0.05519968, -0.00976004, 0.00107857};
	static double a56[N_56] = {8.5733287401, 18.0590169730,
							   8.6347608925, 0.2677737343};
	static double b56[N_56] = {9.5733223454, 25.6329561486,
							   21.0996530827, 3.9584969228};
	double        E1        = 0.0;

	if (x <= 0.0) {
		sprintf(messageStr, "Exponential integral E1 of x = %e", x);
		Error(ERROR_LEVEL_2, "E1", messageStr);
	} else if (x > 0.0 && x <= 1.0) {
		E1 = -log(x) + a53[0] + x * (a53[1] + x * (a53[2] + x * (a53[3] + x * (a53[4] + x * a53[5]))));
	} else if (x <= X_UNDERFLOW) {
		/* --- 5.1.56 holds for all x >= 1 -- --------------------------- */

		E1 = a56[3] / x + a56[2] + x * (a56[1] + x * (a56[0] + x));
		E1 /= b56[3] + x * (b56[2] + x * (b56[1] + x * (b56[0] + x)));
		E1 *= exp(-x);
	}

	return E1;
}
/* ------- end ---------------------------- E1.c -------------------- */

/* ------- begin -------------------------- E2.c -------------------- */

#define E2_MAX_ITER 100

double E2(double x)
{
	register int i;

	double E2 = 0.0, an, b, c, d, delta;

	if (x <= 0.0) {
		sprintf(messageStr, "Exponential integral E2 of x = %e", x);
		Error(ERROR_LEVEL_2, "E2", messageStr);
	} else if (x < 1.0) {

		/* --- Use the recurrence relation En+1 = 1/n (exp(-x) - xEn(x))
			   (Eq. 5.1.14) --                           -------------- */

		E2 = exp(-x) - x * E1(x);
	} else if (x <= X_UNDERFLOW) {

		/* --- For x >= 1 the recurrence cancels (x E1 -> exp(-x)): use
			   the continued fraction of E_n (5.1.22) with n = 2, by the
			   modified Lentz method -- ----------------------------- */

		b  = x + 2.0;
		c  = 1.0 / DBL_MIN;
		d  = 1.0 / b;
		E2 = d;
		for (i = 1; i <= E2_MAX_ITER; i++) {
			an = -i * (i + 1.0);
			b += 2.0;
			d     = 1.0 / (an * d + b);
			c     = b + an / c;
			delta = c * d;
			E2 *= delta;
			if (fabs(delta - 1.0) < DBL_EPSILON) break;
		}
		E2 *= exp(-x);
	}

	return E2;
}
/* ------- end ---------------------------- E2.c -------------------- */
