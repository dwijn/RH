/* ------- file: -------------------------- stokesgrid.h ------------

	   The polarized quantities of one 2-D or 3-D formal solution (one
	   wavelength and ray) on the grid, which the short and long
	   characteristics interpolate at their up- and downwind points:
	   the total opacity chi_I, the source vector S and the Stokes
	   vector I (B-63). The polarized solvers pass this one structure
	   down instead of the three arrays, and interpolate S and I only
	   through Interpolate_Stokes_3D and Interp_Stokes_2D, and K' only
	   through StokesK_3D and StokesK_2D.

	   --------------------------                      ----------RH-- */

#ifndef __STOKESGRID_H__
#define __STOKESGRID_H__

/* --- Which Stokes vector on the grid to interpolate -- ----------- */

enum stokes_vector { STOKES_SOURCE,
					 STOKES_INTENSITY };

typedef struct {
	int      nspect;
	double  *chi; /* chi_I, Nspace */
	double **S;   /* source vector, 4 x Nspace */
	double **I;   /* Stokes vector, 4 x Nspace */
} StokesGrid;

#endif /* !__STOKESGRID_H__ */

/* ---------------------------------------- stokesgrid.h ------------ */
