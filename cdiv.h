/* ------- file: -------------------------- cdiv.h ------------------

       Complex division for the Voigt routines.

       --------------------------                      ----------RH-- */

#ifndef RH_CDIV_H
#define RH_CDIV_H

#include <complex.h>

/* --- Textbook complex division. The / operator of C guards against
       overflow (each compiler in its own way, e.g. Smith's method),
       which is slower and rounds differently from one compiler to the
       next. The arguments in the Voigt routines stay far from overflow,
       so the plain formula is safe, and it gives the same results
       everywhere. --                                  -------------- */

static inline double complex cdiv(double complex a, double complex b)
{
  double d = creal(b)*creal(b) + cimag(b)*cimag(b);

  return CMPLX((creal(a)*creal(b) + cimag(a)*cimag(b)) / d,
	       (cimag(a)*creal(b) - creal(a)*cimag(b)) / d);
}

#endif /* !RH_CDIV_H */

/* ------- end ---------------------------- cdiv.h ------------------ */
