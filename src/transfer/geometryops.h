/* ------- file: -------------------------- geometryops.h -----------

	   The operations a geometry provides to the transfer driver
	   (Formal(), src/transfer/formal.c): docs/LIBRARY.md §3.4, Q-42.

	   Include after rh.h and spectrum.h (bool_t, ActiveSet).

	   --                                              -------------- */

#ifndef __GEOMETRYOPS_H__
#define __GEOMETRYOPS_H__

#include "negtau.h"

/* --- Each geometry (src/geometry/<name>/formal.c) defines one table,
	   geometryOps, through which the driver does everything that
	   depends on the geometry. The entries marked optional may be
	   NULL --                                         -------------- */

typedef struct GeometryOps {
	const char *name;

	/* --- Ray mu (0 .. atmos.Nrays - 1): its direction cosine to the
		   vertical and its quadrature weight --       -------------- */

	double (*muz)(int mu);
	double (*wmu)(int mu);

	/* --- One ray in one direction (to_obs): the intensity I and, if
		   Psi is not NULL, the diagonal of the operator, from chi and
		   the source function S. If Tri is not NULL, the row of the
		   tridiagonal operator goes there (background scattering,
		   N-51); lin asks for the linear solver (N-52). A geometry
		   without these never gets Tri or lin -- ----------------- */

	void (*solve)(int nspect, int mu, bool_t to_obs, double *chi, double *S,
				  double *I, double *Psi, double *Tri, bool_t lin);

	/* --- The same for I, Q, U and V: S and I are [4][Nspace], Psi the
		   diagonal for I --                           -------------- */

	void (*solve_stokes)(int nspect, int mu, bool_t to_obs, double *chi,
						 double **S, double **I, double *Psi);

	/* --- The negative optical depth check (N-40) of ray mu in the
		   direction to_obs (-1: both, for an angle-independent chi) on
		   the columns of the geometry; and the quadrature with which
		   its Bezier solvers integrate chi (B-68) -- ------------ */

	bool_t (*check_tau)(int nspect, int mu, int to_obs, const double *chi,
						enum negtau_quadrature quad);
	enum negtau_quadrature negtau_bezier;

	/* --- Store the emergent intensity of ray mu in spectrum: I, and
		   Q, U and V from Ipol[1 .. 3] if Ipol is not NULL -- ---- */

	void (*emergent)(int nspect, int mu, const double *I, double **Ipol);

	/* --- Background scattering iteration (N-51): NULL if the geometry
		   has it, else the reason it has not; whether wavelengths it
		   iterates switch to the linear solver (N-52). Optional where it
		   is not supported: the largest albedo below
		   BACKGROUND_LINEAR_TAU (N-52), and adding the row Tri of the
		   operator of one ray to A[3k + i] = dJ_k / du_k-1+i -- -- */

	const char *(*bg_supported)(void);
	bool_t (*bg_linear_fallback)(void);
	double (*bg_max_albedo)(const ActiveSet *as);
	void (*bg_add_operator)(int Nspace, bool_t to_obs, double wmu,
							const double *Tri, const double *sca_c,
							const double *chi, double *A);

	/* --- Optional: a Formal() of its own that replaces the driver
		   (rhsphere, until Q-43) --                   -------------- */

	double (*formal)(int nspect, bool_t eval_operator, bool_t redistribute);
} GeometryOps;

extern const GeometryOps geometryOps;

#endif /* !__GEOMETRYOPS_H__ */

/* ------- end ---------------------------- geometryops.h ----------- */
