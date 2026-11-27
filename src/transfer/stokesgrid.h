/* ------- file: -------------------------- stokesgrid.h ------------

	   The polarized quantities of one 2-D or 3-D formal solution (one
	   wavelength and ray) on the grid, which the short and long
	   characteristics interpolate at their up- and downwind points:
	   the total opacity chi_I, the source vector S and the Stokes
	   vector I (B-63). The polarized solvers pass this one structure
	   down instead of the three arrays, and interpolate S and I only
	   through Interpolate_Stokes_3D and Interp_Stokes_2D, and K' only
	   through StokesK_3D and StokesK_2D.

	   Physical Stokes vectors lie in a cone: I >= |(Q, U, V)| for the
	   Stokes vector and the source vector, and chi_I >= |(chi_Q,
	   chi_U, chi_V)|, |eta'| <= 1, for the absorption part of K'. An
	   interpolation of each component on its own leaves the cone where
	   the weights depend on the data (Steffen-limited Bezier in 2-D,
	   the clamped bicubic in 3-D): |eta'| of 100 to 200 and I < 0 at
	   the Ca II H and K cores of a 3-D atmosphere (B-63). So each of
	   these vectors is interpolated as a magnitude (I, S_I: the scalar
	   interpolation as before, which keeps it positive) and the
	   rapidity of its polarization fraction p = (Q, U, V) / I,

			w = atanh(|p|) p / |p|,     p = tanh(|w|) w / |w|,

	   which maps the open unit ball onto all of R^3: any value of an
	   interpolated w is a physical p. w is interpolated with weights
	   that do not depend on the data (the unclamped bicubic in 3-D,
	   Bezier with unlimited derivatives in 2-D), so the result turns
	   with a rotation of the (Q, U) frame. eta' is such a fraction
	   already; the magneto-optical rho' = rho / chi_I, which has no
	   bound, is interpolated with the same weights. The grid arrays
	   of w and rho' are made once per formal solution, and those of I
	   updated at each point as the sweep computes it.

	   --------------------------                      ----------RH-- */

#ifndef __STOKESGRID_H__
#define __STOKESGRID_H__

/* --- Which Stokes vector on the grid to interpolate -- ----------- */

enum stokes_vector { STOKES_SOURCE,
					 STOKES_INTENSITY };

typedef struct {
	int      nspect, Nspace;
	double  *chi; /* chi_I, Nspace */
	double **S;   /* source vector, 4 x Nspace */
	double **I;   /* Stokes vector, 4 x Nspace */

	/* --- Rapidities of the polarization fractions of S and I and of
		   eta', and rho', each 3 x Nspace, on the grid -- ---------- */

	double *wS[3], *wI[3], *wK[3], *rK[3];
	int     polarizedK, magneto_optical;
	double *buffer;
} StokesGrid;

void StokesGrid_init(RHState *rh, StokesGrid *sg, int nspect, int Nspace, double *chi,
					 double **S, double **I);
void StokesGrid_update_I(StokesGrid *sg, int k);
void StokesGrid_free(StokesGrid *sg);
void Stokes_to_rapidity(double I, double Q, double U, double V, double w[3]);
void Stokes_from_rapidity(double I, const double w[3], double v[4]);

#endif /* !__STOKESGRID_H__ */

/* ---------------------------------------- stokesgrid.h ------------ */
