/* ------- file: -------------------------- geometry.h --------------

	   Version:       rh2.0, 3-D Cartesian, short characteristics
	   Authors:       Han Uitenbroek   (huitenbroek@nso.edu)
					  Jorrit Leenaarts (j.leenaarts@astro.uu.nl)

	   Last modified: Fri Jun 15 10:49:52 2018 --

	   --------------------------                     -----------RH-- */

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include "stokesgrid.h"

/* --- Defines structures to contain the angle quadrature and
	   interpolation mesh for 3-D Cartesian version. --  ------------ */

#define NCC 4

enum boundary { FIXED,
				PERIODIC };
enum boundval { IRRADIATED,
				ZERO,
				THERMALIZED,
				REFLECTIVE };
enum vertical { TOP = 0,
				BOTTOM };
enum direction { UPWIND = 0,
				 DOWNWIND };
enum intersect { XY,
				 XZ,
				 YZ };

typedef struct Stencil  Stencil;
typedef struct Longchar Longchar;

struct Stencil {
	enum intersect plane;
	int            quadrant, xbase[NCC], ybase[NCC], zbase[2];
	double         ds, xkernel[NCC], ykernel[NCC], zkernel[2];
	Longchar      *longchar;
};

struct Longchar {
	int      Nst;
	Stencil *stencil;
};

typedef struct Geometry {
	enum boundary x_boundary, y_boundary;
	enum boundval z_boundary_value[2];
	int           Nx, Ny, Nz, Nplane, Nrays, Nlongchar;
	double        dx, dy, *z, *mux, *muy, *muz, *wmu, *vx, *vy, *vz,
		**Itop, **Ibottom;
	Stencil **stencil_uw, **stencil_dw;
} Geometry;

typedef struct {
	enum intersect plane;
	int            index; /* the plane crossed, relative to the local point */
	double         s, x, y, z;
} Intersect;

/* --- Associated function prototypes --               -------------- */

void readAtmos(RHState *rh, Atmosphere *atmos, Geometry *geometry);

void   getAngleQuadr(RHState *rh, Geometry *geometry);
void   fillMesh(RHState *rh, Geometry *geometry);
void   getBoundary(RHState *rh, Atmosphere *atmos, Geometry *geometry);
bool_t writeFlux(RHState *rh, char *fileName);
void   writeGeometry(RHState *rh, Geometry *geometry);
double Interpolate_3D(RHState *rh, double *f, Geometry *geometry, Stencil *st,
					  int l, int m);
double Interpolate_3D_linear(RHState *rh, double *f, Geometry *geometry, Stencil *st,
							 int l, int m);

void ShortChar(RHState *rh, Geometry *geometry, int nspect, int mu,
			   bool_t to_observer, double *chi, double *S, double *I,
			   double *Psi);
void ShortChar_Stokes(RHState *rh, Geometry *geometry, int nspect, int mu,
					  bool_t to_observer, double *chi,
					  double **S, double **I, double *Psi);
void Piecewise_Linear_3D(RHState *rh, Geometry *geometry, Stencil *st_uw, Stencil *st_dw,
						 int k, int kend, int l, int m, double I_uw,
						 double *chi, double *S, double *I, double *Psi);
void Piecewise_Bezier3_3D(RHState *rh, Geometry *geometry, Stencil *st_uw, Stencil *st_dw,
						  int k, int kend, int l, int m, double I_uw,
						  double *chi, double *S, double *I, double *Psi);

void Piece_Stokes_Bezier3_3D(RHState *rh, Geometry *geometry, Stencil *st_uw,
							 Stencil *st_dw, StokesGrid *sg,
							 int k, int kend, int l, int m,
							 double *I_uw, double *Psi);

double SolveLong(RHState *rh, Geometry *geometry, Longchar *lc, int k, int l, int m,
				 double *chi, double *S, double *I);
void   SolveLongStokes(RHState *rh, Geometry *geometry, Longchar *lc, StokesGrid *sg,
					   int k, int l, int m, double *I_uw);

void StokesK_3D(RHState *rh, StokesGrid *sg, Geometry *geometry, Stencil *st,
				int l, int m, double chi_I, double K[4][4]);
void Interpolate_Stokes_3D(RHState *rh, StokesGrid *sg, enum stokes_vector v,
						   Geometry *geometry, Stencil *st, int l, int m,
						   double out[4]);

#endif /* !__GEOMETRY_H__ */

/* ---------------------------------------- geometry.h -------------- */
