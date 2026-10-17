/* ------- file: -------------------------- humlicek.c --------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Dec 13 17:17:04 2001 --

	   --------------------------                      ----------RH-- */

#include "rh.h"
#include <complex.h>

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* -- Voigt function subroutines in different regions.
	  Relative accuracy 1.0E-04. Also calculates Faraday-Voigt
	  function needed in Stokes radiative transfer.

	  C version.

 See: Humlicek 1982, JQSRT 27, p. 437
	  --                                               -------------- */

/* ------- begin -------------------------- Humlicek1.c ------------- */

double complex Humlicek1(double complex z)
{
	/* --- Approximation in region I --                  -------------- */

	double complex z1, z2;

	z1 = 0.5641896 * z;
	z2 = z * z + 0.5;

	return z1 / z2;
}
/* ------- end ---------------------------- Humlicek1.c ------------- */

/* ------- begin -------------------------- Humlicek2.c ------------- */

double complex Humlicek2(double complex z)
{
	/* --- Approximation in region II --                 -------------- */

	double complex z1, z2, u = z * z;

	z1 = z * (0.5641896 * u + 1.410474);
	z2 = u * (u + 3.0) + 0.75;

	return z1 / z2;
}
/* ------- end ---------------------------- Humlicek2.c ------------- */

#define N_HUMLICEK_3 5

/* ------- begin -------------------------- Humlicek3.c ------------- */

double complex Humlicek3(double complex z)
{
	register int n;

	static double a[N_HUMLICEK_3] =
		{0.5642236, 3.778987, 11.96482, 20.20933, 16.4955};
	static double b[N_HUMLICEK_3] =
		{6.699398, 21.69274, 39.27121, 38.82363, 16.4955};

	double complex z1, z2;

	/* --- Approximation in region III --                -------------- */

	z1 = a[0];
	z2 = z + b[0];

	for (n = 1; n < N_HUMLICEK_3; n++) {
		z1 = z1 * z + a[n];
		z2 = z2 * z + b[n];
	}

	return z1 / z2;
}
/* ------- end ---------------------------- Humlicek3.c ------------- */

#define N_HUMLICEK_4 7

/* ------- begin -------------------------- Humlicek4.c ------------- */

double complex Humlicek4(double complex z)
{
	register int n;

	static double a[N_HUMLICEK_4] =
		{0.56419, 1.320522, 35.7668, 219.031, 1540.787, 3321.99, 36183.31};
	static double b[N_HUMLICEK_4] =
		{1.841439, 61.57037, 364.2191, 2186.181, 9022.228, 24322.84, 32066.6};

	double complex z1, z2, u, expu;

	u = -z * z;

	/* --- Approximation in region IV --                 -------------- */

	z1 = a[0];
	z2 = u + b[0];

	for (n = 1; n < N_HUMLICEK_4; n++) {
		z1 = u * z1 + a[n];
		z2 = u * z2 + b[n];
	}
	expu = cexp(-u);

	return expu - (z * z1) / z2;
}
/* ------- end ---------------------------- Humlicek4.c ------------- */
