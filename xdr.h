/* ------- file: -------------------------- xdr.h -------------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Sep 21 13:36:03 2004 --

	   --------------------------                      ----------RH-- */

#ifndef __XDR_H__
#define __XDR_H__

/* --- To be included in routines that use the XDR (see ``man xdr'').
	   (eXternal Data Representation) format to read/write data.
	   --                                              -------------- */

#include <rpc/xdr.h>

#include "atom.h"
#include "atmos.h"

/* --- Element routines for xdr_vector() and xdr_array(). These have
	   exactly the xdrproc_t signature of Apple's <rpc/xdr.h>,
	   (XDR *, void *, unsigned int), so passing them involves no cast
	   between incompatible function types (calling xdr_double and the
	   like through xdrproc_t would be undefined behaviour). With the
	   variadic xdrproc_t of glibc/libtirpc the (xdrproc_t) cast at the
	   call sites is still needed and is harmless.
	   --                                              -------------- */

static inline bool_t rh_xdr_double(XDR *xdrs, void *p, unsigned int n)
{
	(void) n;
	return xdr_double(xdrs, (double *) p);
}

static inline bool_t rh_xdr_int(XDR *xdrs, void *p, unsigned int n)
{
	(void) n;
	return xdr_int(xdrs, (int *) p);
}

static inline bool_t rh_xdr_bool(XDR *xdrs, void *p, unsigned int n)
{
	(void) n;
	return xdr_bool(xdrs, (bool_t *) p);
}

/* --- Associated function prototypes --               -------------- */

void   convertJ(FILE *fp_J, int Nspect, enum xdr_op xdrOperation);
bool_t readB(Atmosphere *atmos);
bool_t writeFlux(char *fileName);
bool_t writeMetals(char *fileName);
bool_t xdr_atom(XDR *xdrs, Atom *atom);
bool_t xdr_counted_string(XDR *xdrs, char **p);
bool_t xdr_populations(XDR *xdrs, char *atmosID, int Nlevel, int Nspace,
					   double *n, double *nstar);
bool_t xdr_BRS(XDR *xdrs);

#endif /* !__XDR_H__ */

/* ------- end ---------------------------- xdr.h ------------------- */
