/* ------- file: -------------------------- geometry.c --------------

	   Version:       rh2.0, 2-D Cartesian
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Oct 16 11:30:43 2001 --

	   --------------------------                      ----------RH-- */

/* --- Store the up- and downwind weights and indices for interpolation
	   in Geometry structure (defined in geometry.h).
	   Version for rectangular 2-D grid.

	   The Geometry:


											   || || Itop
		   ^ z                                 \/ \/
		   |
	   TOP ========================================================= k=0
		   =   |       |       |    ^ up |              |          =
		   =   |       |       |   /     |              |          =
		   =----------------------.--------------------------------=
		   =   |       |       | / k-1   |              |          =
		   =   |       |     k |/        |              |          =
		   =-------------------.-----------------------------------=
		   =   |       |l-1   /|l        | l+1          |          =
	  LEFT =   |       |     / |         |              |          = RIGHT
		   =   |       |    /  |         |              |dz        =
		   =   |       |   /   |         |              |    dx    =
		   =--------------.----------------------------------------=
		   =   |       | / k+1 |         |              |          =
		   =   |       |/down  |         |              |          =
	BOTTOM =========================================================-> x
									/\ /\
									|| ||  Ibottom


	   Conventions:

		 -- mu[xyz] is the cosine of the angle of a ray with the
		positive [xyz]-axis. (Note: mux^2 + muy^2 + muz^2 = 1.0)

	   --                                              -------------- */

#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "geometry.h"
#include "inputs.h"
#include "statistics.h"
#include "error.h"

/* --- Function prototypes --                          -------------- */

void setStencil(Geometry *geometry, Stencil *stencil,
				enum direction direction, enum intersect intersect,
				double ds, double fraction, int l, int k,
				enum SC_order order, enum sweep sweep);

void setVertStencil(Geometry *geometry, Stencil *stencil,
					enum direction direction, double ds, int l, int k);

void fillLongChar(Geometry *geometry, int mu, int k, enum SC_order order,
				  enum sweep sweep);

/* --- Global variables --                             -------------- */

extern char      messageStr[];
extern InputData input;

/* ------- begin -------------------------- fillMesh.c -------------- */

void fillMesh(Geometry *geometry)
{
	register int mu, k, l;

	enum SC_order order;
	enum boundary hboundary = geometry->hboundary;
	int           Nx = geometry->Nx, Nz = geometry->Nz, Nrays = geometry->Nrays,
		count, Nlongchar = 0, Nstart, Nend;
	double    mux, muz, dx, dz, fraction;
	Stencil  *st;
	LongChar *lc;

	if (input.Eddington) return;

	getCPU(2, TIME_START, NULL);

	/* --- Allocate space for interpolation stencils --  -------------- */

	geometry->stencil = (Stencil **) malloc(Nrays * sizeof(Stencil *));
	for (mu = 0; mu < Nrays; mu++) {
		geometry->stencil[mu] = (Stencil *) malloc(Nx * Nz * sizeof(Stencil));
	}

	/* --- If the horizontal boundary conditions are PERIODIC allocate
		   space for possible long characteristics --    -------------- */

	if (hboundary == PERIODIC) {
		geometry->longchar = (LongChar **) malloc(Nrays * sizeof(LongChar *));
		for (mu = 0; mu < Nrays; mu++)
			geometry->longchar[mu] = (LongChar *) malloc(2 * Nz * sizeof(LongChar));
	} else
		geometry->longchar = NULL;

	order = (input.interpolate_2D == LINEAR_2D) ? SC_LINEAR : SC_BEZIER;

	/* --- For each direction ray go through the x-z grid and store the
		   stencil of interpolation indices and coefficients in the
		   Geometry structure --                         -------------- */

	for (mu = 0; mu < Nrays; mu++) {
		mux = geometry->mux[mu];
		muz = geometry->muz[mu];

		/* --- Distinguish between mux > 0.0, mu = 0.0, and mux < 0.0.
			   Rays with mux > 0.0 run from right to left in the downsweep,
			   ie. with decreasing l. --                   -------------- */

		if (mux > 0.0) {

			/* --- First evaluate the UPWIND interpolation weigths and
				   indices for the DOWN sweep. These will be the DOWNWIND
				   ones in the reverse (UP) sweep of the formal solution
				   --                                        -------------- */

			for (k = 1; k < Nz; k++) {
				dz = geometry->dz[k - 1];

				count = 0;
				if (hboundary == PERIODIC)
					Nend = Nx;
				else
					Nend = Nx - 1;
				for (l = 0; l < Nend; l++) {
					st = &geometry->stencil[mu][k * Nx + l];

					dx       = geometry->dx[l];
					fraction = mux * dz / (muz * dx);

					/* --- Characteristic intersects horizontal grid line? -- - */

					if (fraction <= 1.0) {
						count++;
						setStencil(geometry, st, UPWIND, HORIZONTAL, dz / muz,
								   fraction, l, k - 1, order, DOWN);
					} else {
						setStencil(geometry, st, UPWIND, VERTICAL, dx / mux,
								   1.0 - 1.0 / fraction, l + 1, k - 1, order, DOWN);
					}
				}
				/* --- If no horizontal intersection occurred use long
					   characteristic, otherwise find proper starting point - */

				if (geometry->longchar != NULL) {
					lc = &geometry->longchar[mu][DOWN * Nz + k];
					if (count == 0) {
						lc->lstart = Nx - 1;
						fillLongChar(geometry, mu, k, order, DOWN);
						Nlongchar++;
					} else {
						lc->Nlc    = 0;
						lc->lstart = 0;
						st         = &geometry->stencil[mu][k * Nx + lc->lstart];
						while (st->intersect[UPWIND] == VERTICAL) {
							lc->lstart++;
							st++;
						}
					}
				}
			}
			/* --- Evaluate the DOWNWIND interpolation weights and indices
				   for the DOWN sweep. These will be the UPWIND ones in
				   the reverse (UP) sweep. --                  ------------ */

			for (k = 0; k < Nz - 1; k++) {
				dz = geometry->dz[k];

				count = 0;
				if (hboundary == PERIODIC)
					Nstart = 0;
				else
					Nstart = 1;
				for (l = Nstart; l < Nx; l++) {
					st = &geometry->stencil[mu][k * Nx + l];

					dx       = geometry->dx[MODULO(l - 1, Nx)];
					fraction = mux * dz / (muz * dx);

					/* --- Characteristic intersects horizontal grid line? -- - */

					if (fraction <= 1.0) {
						count++;
						setStencil(geometry, st, DOWNWIND, HORIZONTAL, dz / muz,
								   1.0 - fraction, l - 1, k + 1, order, UP);
					} else {
						setStencil(geometry, st, DOWNWIND, VERTICAL, dx / mux,
								   1.0 / fraction, l - 1, k, order, UP);
					}
				}
				/* --- If no horizontal intersection occurred use long
					   characteristic, otherwise find proper starting point - */

				if (geometry->longchar != NULL) {
					lc = &geometry->longchar[mu][UP * Nz + k];
					if (count == 0) {
						lc->lstart = 0;
						fillLongChar(geometry, mu, k, order, UP);
						Nlongchar++;
					} else {
						lc->Nlc    = 0;
						lc->lstart = Nx - 1;
						st         = &geometry->stencil[mu][k * Nx + lc->lstart];
						while (st->intersect[DOWNWIND] == VERTICAL) {
							lc->lstart--;
							st--;
						}
					}
				}
			}
		} else if (mux == 0.0) {

			/* --- Case for purely vertical rays --          -------------- */

			for (k = 1; k < Nz; k++) {
				dz = geometry->dz[k - 1];
				for (l = 0; l < Nx; l++) {
					st = &geometry->stencil[mu][k * Nx + l];
					setVertStencil(geometry, st, UPWIND, dz, l, k - 1);
				}
				if (geometry->hboundary == PERIODIC) {
					lc         = &geometry->longchar[mu][DOWN * Nz + k];
					lc->Nlc    = 0;
					lc->lstart = 0;
				}
			}
			for (k = 0; k < Nz - 1; k++) {
				dz = geometry->dz[k];
				for (l = 0; l < Nx; l++) {
					st = &geometry->stencil[mu][k * Nx + l];
					setVertStencil(geometry, st, DOWNWIND, dz, l, k + 1);
				}
				if (geometry->hboundary == PERIODIC) {
					lc         = &geometry->longchar[mu][UP * Nz + k];
					lc->Nlc    = 0;
					lc->lstart = 0;
				}
			}
		} else {

			/* --- Rays with mux < 0.0 run from left to right in the
				   downsweep, ie. with increasing l.
				   First evaluate the UPWIND interpolation weigths and
				   indices for the DOWN sweep. These will be the DOWNWIND
				   ones in the reverse (UP) sweep of the formal solution
				   --                                        -------------- */

			for (k = 1; k < Nz; k++) {
				dz = geometry->dz[k - 1];

				count = 0;
				if (hboundary == PERIODIC)
					Nstart = 0;
				else
					Nstart = 1;
				for (l = Nstart; l < Nx; l++) {
					st = &geometry->stencil[mu][k * Nx + l];

					dx       = geometry->dx[MODULO(l - 1, Nx)];
					fraction = -mux * dz / (muz * dx);

					/* --- Characteristic intersects horizontal grid line? -- - */

					if (fraction <= 1.0) {
						count++;
						setStencil(geometry, st, UPWIND, HORIZONTAL, dz / muz,
								   1.0 - fraction, l - 1, k - 1, order, DOWN);
					} else {
						setStencil(geometry, st, UPWIND, VERTICAL, -dx / mux,
								   1.0 - 1.0 / fraction, l - 1, k - 1, order, DOWN);
					}
				}
				/* --- If no horizontal intersection occurred use long
					   characteristic, otherwise find proper starting point - */

				if (geometry->longchar != NULL) {
					lc = &geometry->longchar[mu][DOWN * Nz + k];
					if (count == 0) {
						lc->lstart = 0;
						fillLongChar(geometry, mu, k, order, DOWN);
						Nlongchar++;
					} else {
						lc->Nlc    = 0;
						lc->lstart = Nx - 1;
						st         = &geometry->stencil[mu][k * Nx + lc->lstart];
						while (st->intersect[UPWIND] == VERTICAL) {
							lc->lstart--;
							st--;
						}
					}
				}
			}
			/* --- Evaluate the DOWNWIND interpolation weights and indices
				   for the DOWN sweep. These will be the UPWIND ones in
				   the reverse (UP) sweep. --                  ------------ */

			for (k = 0; k < Nz - 1; k++) {
				dz = geometry->dz[k];

				count = 0;
				if (hboundary == PERIODIC)
					Nend = Nx;
				else
					Nend = Nx - 1;
				for (l = 0; l < Nend; l++) {
					st = &geometry->stencil[mu][k * Nx + l];

					dx       = geometry->dx[l];
					fraction = -mux * dz / (muz * dx);

					/* --- Characteristic intersects horizontal grid line? -- - */

					if (fraction <= 1.0) {
						count++;
						setStencil(geometry, st, DOWNWIND, HORIZONTAL, dz / muz,
								   fraction, l, k + 1, order, UP);
					} else {
						setStencil(geometry, st, DOWNWIND, VERTICAL, -dx / mux,
								   1.0 / fraction, l + 1, k, order, UP);
					}
				}
				/* --- If no horizontal intersection occurred use long
					   characteristic, otherwise find proper starting point - */

				if (geometry->longchar != NULL) {
					lc = &geometry->longchar[mu][UP * Nz + k];
					if (count == 0) {
						lc->lstart = Nx - 1;
						fillLongChar(geometry, mu, k, order, UP);
						Nlongchar++;
					} else {
						lc->Nlc    = 0;
						lc->lstart = 0;
						st         = &geometry->stencil[mu][k * Nx + lc->lstart];
						while (st->intersect[DOWNWIND] == VERTICAL) {
							lc->lstart++;
							st++;
						}
					}
				}
			}
		}
	}

	sprintf(messageStr, "\n-Used %d Long Characteristics\n\n", Nlongchar);
	Error(MESSAGE, NULL, messageStr);
	getCPU(2, TIME_POLL, "Geometry");
}
/* ------- end ---------------------------- fillMesh.c -------------- */

/* ------- begin -------------------------- setVertStencil.c -------- */

void setVertStencil(Geometry *geometry, Stencil *stencil,
					enum direction direction, double ds, int l, int k)
{
	int Nx = geometry->Nx;

	/* --- Purely vertical ray: the up- or downwind point is grid
		   point (l, k) itself --                         -------------- */

	stencil->intersect[direction] = HORIZONTAL;
	stencil->ds[direction]        = ds;
	stencil->fraction[direction]  = 0.0;
	stencil->order[direction]     = SC_LINEAR;

	for (int n = 0; n < 4; n++) stencil->index[direction][n] = k * Nx + l;
	for (int n = 0; n < 3; n++) stencil->ih[direction][n] = 0.0;
	stencil->a[direction][0] = stencil->a[direction][1] = 0.0;
	stencil->hc[direction]                              = 0.0;
}
/* ------- end ---------------------------- setVertStencil.c -------- */

/* ------- begin -------------------------- setStencil.c ------------ */

void setStencil(Geometry *geometry, Stencil *stencil,
				enum direction direction, enum intersect intersect,
				double ds, double fraction, int l, int k,
				enum SC_order order, enum sweep sweep)
{
	const char routineName[] = "setStencil";

	const bool_t periodic = (geometry->hboundary == PERIODIC);
	const int    Nx = geometry->Nx, Nz = geometry->Nz;
	int          p[4];
	double       h[3];

	/* --- The characteristic crosses the grid line between points l
		   and l + 1 of layer k (HORIZONTAL), or between layers k and
		   k + 1 on column line l (VERTICAL), at fraction of the way
		   from l or k. Points p[1] and p[2] bound that interval; the
		   cubic Bezier interpolation also uses the neighbours p[0] and
		   p[3], which are -1 beyond a non-periodic boundary (the top
		   and bottom always, the sides with FIXED boundaries).

		   The stencil interpolates I in the given sweep, which has
		   computed layer k (or k + 1) completely before it reaches the
		   points that use the stencil, but not the layers further down
		   (DOWN) or up (UP). A vertical crossing therefore leaves out
		   the neighbour on that side, as the quadratic interpolation
		   did with its choice of triplet -- -------------- */

	if (fraction < 0.0) {
		Error(ERROR_LEVEL_2, routineName,
			  "Fractional index for interpolation should be > 0.0\n");
	}
	if (direction != UPWIND && direction != DOWNWIND)
		Error(ERROR_LEVEL_2, routineName, "Invalid direction index for Stencil");

	stencil->intersect[direction] = intersect;
	stencil->ds[direction]        = ds;
	stencil->fraction[direction]  = fraction;
	stencil->order[direction]     = order;

	if (intersect == HORIZONTAL) {
		for (int n = 0; n < 4; n++) {
			p[n] = l + n - 1;
			if (periodic)
				p[n] = MODULO(p[n], Nx);
			else if (p[n] < 0 || p[n] >= Nx)
				p[n] = -1;
			stencil->index[direction][n] = (p[n] < 0) ? -1 : k * Nx + p[n];
		}
		for (int n = 0; n < 3; n++)
			h[n] = (p[n] < 0 || p[n + 1] < 0) ? 0.0 : geometry->dx[p[n]];
	} else {
		if (periodic) l = MODULO(l, Nx);
		for (int n = 0; n < 4; n++) {
			p[n]                         = (k + n - 1 < 0 || k + n - 1 >= Nz) ? -1 : k + n - 1;
			stencil->index[direction][n] = (p[n] < 0) ? -1 : p[n] * Nx + l;
		}
		if (sweep == DOWN)
			p[3] = stencil->index[direction][3] = -1;
		else
			p[0] = stencil->index[direction][0] = -1;
		for (int n = 0; n < 3; n++)
			h[n] = (p[n] < 0 || p[n + 1] < 0) ? 0.0 : geometry->dz[p[n]];
	}
	if (p[1] < 0 || p[2] < 0 || h[1] <= 0.0)
		Error(ERROR_LEVEL_2, routineName, "Interpolation interval outside the grid");
	setStencilWeights(stencil, direction, h);
}
/* ------- end ---------------------------- setStencil.c ------------ */

/* ------- begin -------------------------- setStencilWeights.c ----- */

void setStencilWeights(Stencil *stencil, enum direction direction,
					   const double h[3])
{
	/* --- Spacings h[0..2] of the intervals index[0..3] (0 where a
		   point is missing) to the weights Interp_2D uses, so that it
		   needs no division --                           -------------- */

	/* --- Where a point is missing (h = 0) Interp_2D does not read its
		   1/h, which is then 1/DBL_MIN, finite. A guarded division,
		   h[n] > 0.0 ? 1.0 / h[n] : 0.0, traps on 1/0 with the
		   floating-point exceptions of SetFPEtraps: compilers (icx,
		   clang) turn the guard into a select and compute the division
		   for every n anyway, and they do the same with a guarded
		   divisor. fmax is not a select -- ----------------------- */

	for (int n = 0; n < 3; n++)
		stencil->ih[direction][n] = 1.0 / fmax(h[n], DBL_MIN);
	stencil->a[direction][0] = (h[0] > 0.0) ? h[1] / (h[0] + h[1]) : 0.0;
	stencil->a[direction][1] = (h[2] > 0.0) ? h[2] / (h[1] + h[2]) : 0.0;
	stencil->hc[direction]   = h[1] / 3.0;
}
/* ------- end ---------------------------- setStencilWeights.c ----- */

/* ------- begin -------------------------- fillLongChar.c ---------- */

void fillLongChar(Geometry *geometry, int mu, int k, enum SC_order order,
				  enum sweep sweep)
{
	register int l, n;

	int Nx = geometry->Nx, Nz = geometry->Nz, dl, k_upwind, count,
		l_intersect;
	double mux, muz, dx_total, dx_intersect, dz, dz_total, ds_total,
		ds, fraction;
	LongChar *lc;

	/* --- Fills the appropriate LongChar structure of the interpolation
		   mesh with interpolation stencils to solve along a
		   long characteristic when, with periodic or even boundary
		   conditions, no shortcharacteristic intersects a horizontal grid
		   line for a given k. --                        -------------- */

	lc  = geometry->longchar[mu] + (sweep * Nz + k);
	mux = geometry->mux[mu];
	muz = geometry->muz[mu];

	l = lc->lstart;
	if (sweep == DOWN) {
		k_upwind = k - 1;
		dz       = geometry->dz[k_upwind];
		if (mux >= 0.0) {
			dl = 1;
		} else {
			dl = -1;
		}
	} else {
		k_upwind = k + 1;
		dz       = geometry->dz[k];
		if (mux >= 0.0) {
			dl = -1;
		} else {
			dl = 1;
		}
	}
	/* --- Find first intersection with horizontal grid line -- ------- */

	count        = 0;
	dx_total     = 0.0;
	dx_intersect = dz / muz * fabs(mux);

	if (dl > 0) {
		while (dx_total < dx_intersect) {
			dx_total += geometry->dx[l];
			l = MODULO(l + 1, Nx);
			count++;
		}
		l = MODULO(l - 1, Nx);
	} else {
		while (dx_total < dx_intersect) {
			l = MODULO(l - 1, Nx);
			dx_total += geometry->dx[l];
			count++;
		}
	}
	/* --- Allocate space for long characteristic indices, interpolation
		   coefficients, and pathlength.

	 Note: Interpolation data is stored in an array of Stencil structures
		   so that setStencil can determine the indices and weights, and
		   Interp_2D can do the actual interpolation. Only the UPWIND
		   half of each stencil is used because the long characteristic
		   is solved only in one direction, towards the root. -- ------ */

	lc->Nlc     = count;
	lc->stencil = (Stencil *) malloc(count * sizeof(Stencil));

	/* --- Store the interpolation in the horizontal grid line -- ----- */

	if (dl > 0) {
		fraction = 1.0 - (dx_total - dx_intersect) / geometry->dx[l];
		ds       = fraction * geometry->dx[l] / fabs(mux);
	} else {
		fraction = (dx_total - dx_intersect) / geometry->dx[l];
		ds       = (1.0 - fraction) * geometry->dx[l] / fabs(mux);
	}
	setStencil(geometry, &lc->stencil[0], UPWIND, HORIZONTAL, ds,
			   fraction, l, k_upwind, order, sweep);

	/* --- Then segment for segment the vertical crossings -- --------- */

	ds_total = 0.0;
	for (n = 1; n < count; n++) {
		ds_total += lc->stencil[n - 1].ds[UPWIND];
		dz_total = ds_total * muz;
		if (dl > 0)
			l_intersect = l;
		else
			l_intersect = MODULO(l + 1, Nx);

		l  = MODULO(l - dl, Nx);
		ds = geometry->dx[l] / fabs(mux);
		if (sweep == DOWN) {
			fraction = dz_total / dz;
			setStencil(geometry, &lc->stencil[n], UPWIND, VERTICAL, ds,
					   fraction, l_intersect, k_upwind, order, sweep);
		} else {
			fraction = 1.0 - dz_total / dz;
			setStencil(geometry, &lc->stencil[n], UPWIND, VERTICAL, ds,
					   fraction, l_intersect, k, order, sweep);
		}
	}
}
/* ------- end ---------------------------- fillLongChar.c ---------- */
