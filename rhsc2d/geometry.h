/* ------- file: -------------------------- geometry.h --------------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Dec 17 1996

	   --------------------------                      ----------RH-- */

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <math.h>

/* --- Defines structures to contain the angle quadrature and
	   interpolation mesh for 2-D Cartesian version. -- ------------- */

enum boundary { FIXED,
				PERIODIC };
enum boundval { IRRADIATED,
				ZERO,
				THERMALIZED,
				REFLECTIVE };
enum horizontal { LEFT = 0,
				  RIGHT };
enum vertical { TOP = 0,
				BOTTOM };
enum direction { UPWIND = 0,
				 DOWNWIND };
enum intersect { VERTICAL = 0,
				 HORIZONTAL };
enum SC_order { SC_LINEAR,
				SC_BEZIER };
enum sweep { DOWN = 0,
			 UP };

/* --- The Stencil structure holds, in the UPWIND and DOWNWIND
	   direction for each ray at each spatial location in the grid, the
	   interval [index[1], index[2]] of grid points that the
	   characteristic crosses, the fraction of the way from index[1],
	   and for the monotonic cubic Bezier interpolation (SC_BEZIER) the
	   neighbours index[0] and index[3] (-1 beyond a non-periodic
	   boundary), with weights from the spacings h[0..2] of the three
	   intervals (setStencilWeights): ih = 1 / h (0 where a neighbour
	   is missing), hc = h[1] / 3, and a[0] = h[1] / (h[0] + h[1]),
	   a[1] = h[2] / (h[1] + h[2]) -- -------------- */

typedef struct {
	enum intersect intersect[2];
	enum SC_order  order[2];
	int            index[2][4];
	double         ih[2][3], a[2][2], hc[2], ds[2], fraction[2];
} Stencil;

typedef struct {
	int      Nlc, lstart;
	Stencil *stencil;
} LongChar;

typedef struct {
	enum boundary hboundary;
	enum boundval bvalue[2];
	int           Nx, Nz, Nrays;
	double       *x, *dx, *z, *dz, *mux, *muy, *muz, *wmu, *vx, *vz, **Itop,
		**Ibottom, **Ileft, **Iright;
	Stencil  **stencil;
	LongChar **longchar;
} Geometry;

/* --- Associated function prototypes --               -------------- */

void   getAngleQuadr(Geometry *geometry);
void   fillMesh(Geometry *geometry);
bool_t writeFlux(char *fileName);
void   writeGeometry(Geometry *geometry);
void   readAtmos(Atmosphere *atmos, Geometry *geometry);
void   getBoundary(Atmosphere *atmos, Geometry *geometry);
void   setStencilWeights(Stencil *stencil, enum direction direction,
						 const double h[3]);
void   StokesK_2D(int nspect, Stencil *st, enum direction direction,
				  double chi_I, double K[4][4]);

/* --- Interpolation at the up- and downwind points, inline: the
	   formal solvers call it five times per point with one stencil
	   --                                              -------------- */

/* ------- begin -------------------------- steffen.c --------------- */

static inline double steffen(double s_up, double s_dn, double p)
{
	/* --- Derivative at a point between intervals with secants s_up
		   and s_dn, of which p is the weighted mean (the derivative of
		   the parabola through the three points), limited as in Steffen
		   (1990): 0 at an extremum, else at most twice the smaller
		   secant and at most |p|. The same as cent_deriv -- --------- */

	/* --- MIN, not fmin: without -ffast-math fmin is a call to libm,
		   which made this interpolation 1.6 times slower -- --------- */

	const double m = MIN(MIN(fabs(s_up), fabs(s_dn)), 0.5 * fabs(p));

	return (s_up * s_dn > 0.0) ? copysign(2.0 * m, s_dn) : 0.0;
}
/* ------- end ---------------------------- steffen.c --------------- */

/* ------- begin -------------------------- end_deriv.c ------------- */

static inline double end_deriv(double s_near, double s_far, double e)
{
	/* --- Derivative at the end of an interval with secant s_near, from
		   the parabola through it and the next interval (secant s_far,
		   e = h_near / (h_near + h_far)), limited as at the ends of a
		   Steffen (1990) spline: between 0 and twice s_near -- -------- */

	const double p = s_near + (s_near - s_far) * e;

	if (p * s_near <= 0.0) return 0.0;
	if (fabs(p) > 2.0 * fabs(s_near)) return 2.0 * s_near;
	return p;
}
/* ------- end ---------------------------- end_deriv.c ------------- */

/* ------- begin -------------------------- Interp_2D.c ------------- */

static inline double Interp_2D(const double *values, const Stencil *stencil,
							   enum direction direction)
{
	const int   *index = stencil->index[direction];
	const double t = stencil->fraction[direction], y0 = values[index[1]],
				 y1 = values[index[2]];

	if (stencil->order[direction] == SC_LINEAR) return (1.0 - t) * y0 + t * y1;

	/* --- Monotonic cubic Bezier interpolation between grid points
		   index[1] and index[2], with the derivatives of Steffen (1990)
		   at both ends, as the Bezier formal solvers use along the ray.
		   The Steffen derivative lies between 0 and twice the secant,
		   so the control points stay between y0 and y1: the result
		   never leaves the range of the two end values, and it is a
		   continuous function of the values. The quadratic
		   interpolation this replaces fell back to linear where it
		   overshot, a jump that left the iteration cycling (N-37).
		   Where the neighbour beyond one end is missing (index -1), the
		   derivative there is the one-sided one of Steffen (1990), from
		   the parabola through the three other points
		   --                                              -------------- */

	const double *ih = stencil->ih[direction], *a = stencil->a[direction];
	const double  s    = (y1 - y0) * ih[1];
	const double  s_up = (index[0] >= 0) ? (y0 - values[index[0]]) * ih[0] : 0.0,
				 s_dn  = (index[3] >= 0) ? (values[index[3]] - y1) * ih[2] : 0.0;
	double d0, d1;

	if (index[0] >= 0)
		d0 = steffen(s_up, s, a[0] * s_up + (1.0 - a[0]) * s);
	else if (index[3] >= 0)
		d0 = end_deriv(s, s_dn, 1.0 - a[1]);
	else
		d0 = s;
	if (index[3] >= 0)
		d1 = steffen(s, s_dn, a[1] * s + (1.0 - a[1]) * s_dn);
	else if (index[0] >= 0)
		d1 = end_deriv(s, s_up, a[0]);
	else
		d1 = s;

	const double hc = stencil->hc[direction], c0 = y0 + hc * d0, c1 = y1 - hc * d1,
				 u = 1.0 - t;

	return u * u * (u * y0 + 3.0 * t * c0) + t * t * (3.0 * u * c1 + t * y1);
}
/* ------- end ---------------------------- Interp_2D.c ------------- */

/* --- Formal solution related --                      -------------- */

double SolveLong(LongChar *lc, int local, double *chi, double *S,
				 double *I);

double Feautrier(int nspect, int lx, double *chi, double *S,
				 enum FeautrierOrder order, double *P, double *Psi);

void Piecewise_Linear_2D(Geometry *geometry, int nspect, int mu,
						 bool_t to_observer, double *chi, double *S,
						 double *I, double *Psi);
void Piecewise_Bezier_2D(Geometry *geometry, int nspect, int mu,
						 bool_t to_observer, double *chi, double *S,
						 double *I, double *Psi);
void Piecewise_Hermite_2D(Geometry *geometry, int nspect, int mu,
						  bool_t to_observer, double *chi, double *S,
						  double *I, double *Psi);

/* --- Prototype for parabolic DELO solver --          -------------- */

void SolveLongStokes(int nspect, LongChar *lc, int local, double *chi,
					 double **S, double **I, double *I_uw);

/* --- Prototypes for cubic Bezier solvers --          -------------- */

void Piecewise_Bezier3_2D(Geometry *geometry, int nspect, int mu,
						  bool_t to_observer, double *chi, double *S,
						  double *I, double *Psi);

void Piece_Stokes_Bezier3_2D(Geometry *geometry, int nspect, int mu,
							 bool_t to_observer, double *chi, double **S,
							 double **I, double *Psi);

#endif /* !__GEOMETRY_H__ */

/* ---------------------------------------- geometry.h -------------- */
