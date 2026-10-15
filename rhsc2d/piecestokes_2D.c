/* ------- file: -------------------------- piecestokes_2D.c --------

       Version:       rh2.0, 2-D Cartesian
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Sep 23 12:23:18 2025 --

       --------------------------                      ----------RH-- */

/* --- Formal 2-D polarized radiative transfer solver using
       the method of short characteristics. The transfer equation is
       solved in one direction along the ray only. Use keyword
       to_observer=TRUE to solve in the upward direction, and
       to_observer=FALSE to solve in opposite direction.

  See: P. B. Kunasz and  L. H. Auer 1988, JQSRT 39, 67-79

       For now it implements fixed and periodic boundary
       conditions in the horizontal direction.

  See: L. H. Auer, P. Fabiani Bendicho and J. Trujillo Bueno 1994,
       A&A, 292, 599-615

       At the top and bottom only a fixed boundary is allowed.

  See: L. H. Auer 2003, Formal Solution: EXPLICIT Answers, in 
       Stellar Atmosphere Modeling, ASP Conference Proceedings, Vol. 288,
       eds. Hubeny, Mihalas & Werner.


 Note: Rays with to_observer == TRUE go from BOTTOM to TOP.

       --                                              -------------- */

 
#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "geometry.h"
#include "inputs.h"
#include "error.h"


/* --- Function prototypes --                          -------------- */


/* --- Global variables --                             -------------- */

extern Atmosphere atmos;
extern Spectrum spectrum;
extern InputData input;


/* ------- begin -------------------------- StokesK_2D.c ------------ */

void StokesK_2D(int nspect, Stencil *st, enum direction direction,
		bool_t monotonic, double chi_I, double K[4][4])
{
  register int i, j;

  ActiveSet *as;

  /* --- Return the elements of the reduced 4x4 Stokes
         opacity matrix K', which is defined as:

           =           =         =              =
           K' = (chi_c*1 + chi_l*Phi) / chi_I - 1,

	   for wavelength# nspect, spatial point k, and ray mu.

         This is the 2-D short characteristics version.

    See: Rees, Murphy, & Durrant, 1989, ApJ 339, 1093-1106.

         --                                            -------------- */
  as = &spectrum.as[nspect];

  for (j = 0;  j < 4;  j++)
    for (i = 0;  i < 4;  i++) K[j][i] = 0.0;

  /* --- First the contribution from the active set -- -------------- */

  if (containsPolarized(as)) {
    K[0][1] = Quadr(as->chi + atmos.Nspace,
                    st, direction, monotonic);
    K[0][2] = Quadr(as->chi + 2*atmos.Nspace,
                    st, direction, monotonic);
    K[0][3] = Quadr(as->chi + 3*atmos.Nspace,
                    st, direction, monotonic);

    if (input.magneto_optical) {
      K[1][2] = Quadr(as->chip + 2*atmos.Nspace,
                      st, direction, monotonic);
      K[1][3] = Quadr(as->chip + atmos.Nspace,
                      st, direction, monotonic);
      K[2][3] = Quadr(as->chip,
                      st, direction, monotonic);
    }
  }
  /* --- Add possible contribution from the background -- ----------- */

  if (atmos.backgrflags[nspect].ispolarized) {
    K[0][1] += Quadr(as->chi_c + atmos.Nspace,
                     st, direction, monotonic);
    K[0][2] += Quadr(as->chi_c + 2*atmos.Nspace,
                     st, direction, monotonic);
    K[0][3] += Quadr(as->chi_c + 3*atmos.Nspace,
                     st, direction, monotonic);

    if (input.magneto_optical) {
      K[1][2] += Quadr(as->chip_c + 2*atmos.Nspace,
                      st, direction, monotonic);
      K[1][3] += Quadr(as->chip_c + atmos.Nspace,
                      st, direction, monotonic);
      K[2][3] += Quadr(as->chip_c,
                      st, direction, monotonic);
    }
  }
  /* --- Divide by Stokes I opacity and fill lower diagonal part -- - */

  for (j = 0;  j < 3;  j++) {
    for (i = j+1;  i < 4;  i++) {
      K[j][i] /= chi_I;
      K[i][j]  = K[j][i];
    }
  }
  /* --- Anti-symmetric magneto-optical elements --    -------------- */

  if (input.magneto_optical) {
    K[1][3] *= -1.0;
    K[2][1] *= -1.0;
    K[3][2] *= -1.0;
  }
}
/* ------- end ---------------------------- StokesK_2D.c ------------ */


