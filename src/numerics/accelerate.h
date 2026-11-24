/* ------- begin -------------------------- Accelerate.h ------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Dec 28 14:02:56 1999 --

	   --------------------------                      ----------RH-- */

#ifndef __ACCELERATE_H__
#define __ACCELERATE_H__

/* --- Defines structure and prototypes for Ng acceleration -- ------ */

struct Ng {
	int      N, Ndelay, Norder, Nperiod, count;
	double **previous, **A, *b, *theStorage;
	/* --- Values the last extrapolation made non-positive or
		   non-finite; it was rejected if > 0 (N-33) --  ------------ */
	long rejected;
};

/* --- Associated function prototypes --               -------------- */

bool_t      Accelerate(struct Ng *Ngs, double *solution);
void        NgFree(struct Ng *Ngs);
struct Ng  *NgInit(int N, int Ndelay, int Norder, int Nperiod,
				   double *solution);
double      MaxChange(struct Ng *Ngs, char *text, bool_t quiet);
const char *NgStatus(struct Ng *Ngs, bool_t accel);

#endif /* !__ACCELERATE_H__ */

/* ------- end ---------------------------- Accelerate.h ------------ */
