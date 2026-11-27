/* ------- file: -------------------------- writen_xdr.c ------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Jan 26 14:50:56 2012 --

	   --------------------------                      ----------RH-- */

/* --- Routines for reading and writing populations from and to file.
	   XDR (external data representation) version. --  -------------- */

#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "error.h"
#include "xdr.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- xdr_populations.c ------- */

bool_t xdr_populations(XDR *xdrs, char *atmosID, int Nlevel, int Nspace,
					   double *n, double *nstar)
{
	const char routineName[] = "xdr_populations";

	char  *ID;
	bool_t result = TRUE;
	int    Npop   = Nlevel * Nspace, Nl, Ns;

	/* --- The actual reading/writing routine. Upon input the values
		   for atmosID, Nlevel and Nspace in the file are checked against
		   the values supplied in the call: populations of another
		   atmosphere (ID) or of another model stop the run (B-61).
		   --                                            -------------- */

	if (xdrs->x_op == XDR_ENCODE) {
		result &= xdr_counted_string(xdrs, &atmosID);
		result &= xdr_int(xdrs, &Nlevel);
		result &= xdr_int(xdrs, &Nspace);
	} else {
		ID = NULL;
		if (!xdr_counted_string(xdrs, &ID)) {
			free(ID);
			return FALSE;
		}
		if (strcmp(ID, atmosID) != 0) {
			Error(rh, ERROR_LEVEL_2, routineName,
				  "Populations were computed with atmosphere \"%s\",\n"
				  " not with the current one, \"%s\"",
				  ID, atmosID);
			free(ID);
		}
		free(ID);

		result &= xdr_int(xdrs, &Nl);
		result &= xdr_int(xdrs, &Ns);
		if (result && ((Nl != Nlevel) || (Ns != Nspace))) {
			Error(rh, ERROR_LEVEL_2, routineName,
				  "Populations have %d levels at %d depths, the model and\n"
				  " atmosphere %d levels at %d depths",
				  Nl, Ns, Nlevel, Nspace);
		}
	}
	/* --- Exit if true populations do not exist --      -------------- */

	if (n != NULL) {
		result &= xdr_vector(xdrs, (char *) n, Npop,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
	} else
		return FALSE;

	/* --- We can live without the LTE values since they can always be
		   constructed with routine LTEpops --           -------------- */

	if (nstar != NULL)
		result &= xdr_vector(xdrs, (char *) nstar, Npop,
							 sizeof(double), (xdrproc_t) rh_xdr_double);

	return result;
}
/* ------- end ---------------------------- xdr_populations.c ------- */

/* ------- begin -------------------------- writePopulations.c ------ */

void writePopulations(Atom *atom)
{
	const char routineName[] = "writePopulations";

	FILE *fp_out;
	XDR   xdrs;

	if ((fp_out = fopen(atom->popsoutFile, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open output file %s",
			  atom->popsoutFile);
		return;
	}
	xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

	if (!xdr_populations(&xdrs, rh->atmos.ID, atom->Nlevel, rh->atmos.Nspace,
						 atom->n[0], atom->nstar[0])) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to write to output file %s",
			  atom->popsoutFile);
	}

	xdr_destroy(&xdrs);
	fclose(fp_out);
}
/* ------- end ---------------------------- writePopulations.c ------ */

/* ------- begin -------------------------- readPopulations.c ------- */

void readPopulations(Atom *atom)
{
	const char routineName[] = "readPopulations";

	FILE *fp_in;
	XDR   xdrs;

	/* --- Read populations from file.

	 Note: readPopulations only reads the true populations and not
		   the LTE populations, which are those of the current
		   atmosphere and electron density (SetLTEQuantities). To this
		   effect it passes a NULL pointer to xdr_populations as its
		   last argument (it passed nstar until B-60).
		   --                                            -------------- */

	if ((fp_in = fopen(atom->popsinFile, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to open input file %s",
			  atom->popsinFile);
	}
	xdrstdio_create(&xdrs, fp_in, XDR_DECODE);

	if (!xdr_populations(&xdrs, rh->atmos.ID, atom->Nlevel, rh->atmos.Nspace,
						 atom->n[0], NULL)) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read from input file %s",
			  atom->popsinFile);
	}
	xdr_destroy(&xdrs);
	fclose(fp_in);
}
/* ------- end ---------------------------- readPopulations.c ------- */

/* ------- begin -------------------------- writeMolPops.c ---------- */

void writeMolPops(struct Molecule *molecule)
{
	const char routineName[] = "writeMolPops";

	FILE *fp_out;
	XDR   xdrs;

	/* --- Write molecular (vibration) populations to file. -- -------- */

	if ((fp_out = fopen(molecule->popsFile, "w")) == NULL) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to open output file %s",
			  molecule->popsFile);
		return;
	}
	xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

	if (!xdr_populations(&xdrs, rh->atmos.ID, molecule->Nv, rh->atmos.Nspace,
						 molecule->nv[0], molecule->nvstar[0])) {
		Error(rh, ERROR_LEVEL_1, routineName, "Unable to write to output file %s",
			  molecule->popsFile);
	}
	xdr_destroy(&xdrs);
	fclose(fp_out);
}
/* ------- end ---------------------------- writeMolPops.c ---------- */

/* ------- begin -------------------------- readMolPops.c ----------- */

void readMolPops(struct Molecule *molecule)
{
	const char routineName[] = "readMolPops";

	FILE *fp_in;
	XDR   xdrs;

	/* --- Read molecular (vibration) populations from file.

	 Note: readMolPops only reads the true populations and not
		   the LTE populations. To this effect it passes a NULL pointer
		   to xdr_populations as its last argument.
		   --                                            -------------- */

	if ((fp_in = fopen(molecule->popsFile, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to open input file %s",
			  molecule->popsFile);
	}
	xdrstdio_create(&xdrs, fp_in, XDR_DECODE);

	if (!xdr_populations(&xdrs, rh->atmos.ID, molecule->Nv, rh->atmos.Nspace,
						 molecule->nv[0], NULL)) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read from input file %s",
			  molecule->popsFile);
	}
	xdr_destroy(&xdrs);
	fclose(fp_in);
}
/* ------- end ---------------------------- readPopulations.c ------- */
