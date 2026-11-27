/* ------- file: -------------------------- backgropac_xdr.c --------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue May 12 17:52:18 2020 --

	   --------------------------                      ----------RH-- */

/* --- Write separate cotributions to background opacity to file -- - */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "background.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "state.h"

#define COMMENT_CHAR "#"

enum opac_type { ABSORPTION,
				 SCATTERING };

/* --- Function prototypes --                          -------------- */

bool_t write_contrib_xdr(XDR *xdrs, char *label, enum opac_type type,
						 int Nspace, double *scatt);

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- backgrOpac.c ------------ */

void backgrOpac(RHState *rh, int Nlambda, double *lambda)
{
	const char   routineName[] = "backgrOpac";
	register int nspect, n, k;

	char metalID[12], format[21], contribution_out[MAX_LINE_SIZE],
		inputLine[MAX_LINE_SIZE];
	bool_t         exit_on_EOF, do_fudge;
	int            Nfudge;
	enum opac_type type;
	double        *chi, *eta, *scatt, *thomson,
		Hmin_fudge, scatt_fudge, metal_fudge, *lambda_fudge = NULL,
											  **fudge = NULL;
	Atom *He;
	FILE *fp_out, *fp_fudge;
	XDR   xdrs;

	/* --- Only evaluate for static case --              -------------- */

	rh->atmos.moving = FALSE;
	rh->atmos.Stokes = FALSE;

	if (rh->input.solve_ne) Solve_ne(rh, rh->atmos.ne, TRUE);

	SetLTEQuantities(rh);

	if (rh->input.NonICE)
		readMolecules(rh, MOLECULAR_CONCENTRATION_FILE);
	else
		ChemicalEquilibrium(rh, N_MAX_CHEM_ITER, CHEM_ITER_LIMIT);

	if (strcmp(rh->input.fudgeData, "none")) {
		do_fudge = TRUE;

		/* --- Read wavelength-dependent fudge factors to compensate for
			   missing UV backround line haze --           -------------- */

		if ((fp_fudge = fopen(rh->input.fudgeData, "r")) == NULL) {
			Error(rh, ERROR_LEVEL_2, routineName, "Unable to open input file %s", rh->input.fudgeData);
		}
		Message(rh, "\n-Fudging background opacities with file\n  %s\n\n",
				rh->input.fudgeData);

		getLine(rh, fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
		sscanf(inputLine, "%d", &Nfudge);
		lambda_fudge = (double *) mallocOrAbort(Nfudge * sizeof(double));
		fudge        = matrix_double(rh, 3, Nfudge);
		for (n = 0; n < Nfudge; n++) {
			getLine(rh, fp_fudge, COMMENT_CHAR, inputLine, exit_on_EOF = TRUE);
			sscanf(inputLine, "%lf %lf %lf %lf", &lambda_fudge[n],
				   &fudge[0][n], &fudge[1][n], &fudge[2][n]);
		}
		for (n = 0; n < 3 * Nfudge; n++) fudge[0][n] += 1.0;
		fclose(fp_fudge);
	} else
		do_fudge = FALSE;

	/* --- Temporary storage for this routine --         -------------- */

	chi     = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	eta     = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	scatt   = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
	thomson = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));

	/* --- Check whether He is present among the metals --  ----------- */

	He = (rh->atmos.elements[1].model) ? rh->atmos.elements[1].model : NULL;

	backgroundPrepare(rh);

	/* --- Go through wavelengths one-by-one --          -------------- */

	for (nspect = 0; nspect < Nlambda; nspect++) {

		/* --- Open output file for fractions of background opacity,
			   emissivity, and scattering --               -------------- */

		sprintf(format, "bopac%%%1.1d.3f",
				(int) log10(lambda[nspect]) + 5);
		sprintf(contribution_out, format, lambda[nspect]);
		if ((fp_out = fopen(contribution_out, "w")) == NULL) {
			Error(rh, ERROR_LEVEL_1, routineName, "Unable to open output file %s",
				  contribution_out);
			return;
		}
		xdrstdio_create(&xdrs, fp_out, XDR_ENCODE);

		printf("%s Processing lambda = %9.3f [nm]\n",
			   (nspect == 0) ? "\n\n " : " ", lambda[nspect]);
		xdr_double(&xdrs, &lambda[nspect]);

		/* --- Thomson scattering by free electrons --     -------------- */

		Thomson(rh, thomson);
		if (do_fudge) {
			Linear(Nfudge, lambda_fudge, fudge[1],
				   1, &lambda[nspect], &scatt_fudge, FALSE);
			for (k = 0; k < rh->atmos.Nspace; k++) thomson[k] *= scatt_fudge;
		}
		write_contrib_xdr(&xdrs, "THOMSON", type = SCATTERING,
						  rh->atmos.Nspace, thomson);

		/* --- Negative hydrogen ion, bound-free and free-free -- ------- */

		if (do_fudge) {
			Linear(Nfudge, lambda_fudge, fudge[0],
				   1, &lambda[nspect], &Hmin_fudge, FALSE);
		}
		if (Hminus_bf(rh, lambda[nspect], chi, eta)) {
			if (do_fudge) {
				for (k = 0; k < rh->atmos.Nspace; k++) chi[k] *= Hmin_fudge;
			}
			write_contrib_xdr(&xdrs, "HMINUS_BF", type = ABSORPTION,
							  rh->atmos.Nspace, chi);
		}
		if (Hminus_ff(rh, lambda[nspect], chi)) {
			if (do_fudge) {
				for (k = 0; k < rh->atmos.Nspace; k++) chi[k] *= Hmin_fudge;
			}
			write_contrib_xdr(&xdrs, "HMINUS_FF", type = ABSORPTION,
							  rh->atmos.Nspace, chi);
		}

		/* --- Bound-free opacities from OH and CH molecules -- --------- */

		if (OH_bf_opac(rh, lambda[nspect], chi, eta)) {
			write_contrib_xdr(&xdrs, "OH_BF", type = ABSORPTION,
							  rh->atmos.Nspace, chi);
		}
		if (CH_bf_opac(rh, lambda[nspect], chi, eta)) {
			write_contrib_xdr(&xdrs, "CH_BF", type = ABSORPTION,
							  rh->atmos.Nspace, chi);
		}

		/* --- Neutral Hydrogen Bound-Free and Free-Free --  ------------ */

		if (Hydrogen_bf(rh, lambda[nspect], chi, eta)) {
			write_contrib_xdr(&xdrs, "HYDROGEN_BF", type = ABSORPTION,
							  rh->atmos.Nspace, chi);
		}
		Hydrogen_ff(rh, lambda[nspect], chi);
		write_contrib_xdr(&xdrs, "HYDROGEN_FF", type = ABSORPTION,
						  rh->atmos.Nspace, chi);

		/* --- Rayleigh scattering by neutral hydrogen --  -------------- */

		if (Rayleigh(rh, lambda[nspect], rh->atmos.H, scatt)) {
			if (do_fudge) {
				for (k = 0; k < rh->atmos.Nspace; k++) scatt[k] *= scatt_fudge;
			}
			write_contrib_xdr(&xdrs, "RAYLEIGH_H", type = SCATTERING,
							  rh->atmos.Nspace, scatt);
		}
		/* --- Rayleigh scattering by neutral helium --    -------------- */

		if (He && Rayleigh(rh, lambda[nspect], He, scatt)) {
			if (do_fudge) {
				for (k = 0; k < rh->atmos.Nspace; k++) scatt[k] *= scatt_fudge;
			}
			write_contrib_xdr(&xdrs, "RAYLEIGH_HE", type = SCATTERING,
							  rh->atmos.Nspace, scatt);
		}
		/* --- Rayleigh scattering and free-free absorption by
			   molecular hydrogen --                       -------------- */

		if (Rayleigh_H2(rh, lambda[nspect], scatt)) {
			if (do_fudge) {
				for (k = 0; k < rh->atmos.Nspace; k++) scatt[k] *= scatt_fudge;
			}
			write_contrib_xdr(&xdrs, "RAYLEIGH_H2", type = SCATTERING,
							  rh->atmos.Nspace, scatt);
		}
		/* --- Absorption by H + H^+ (referred to as H2plus free-free) -- */

		if (H2plus_ff(rh, lambda[nspect], chi))
			write_contrib_xdr(&xdrs, "H2PLUS_FF", type = ABSORPTION,
							  rh->atmos.Nspace, chi);

		if (H2minus_ff(rh, lambda[nspect], chi)) {
			write_contrib_xdr(&xdrs, "H2MINUS_FF", type = ABSORPTION,
							  rh->atmos.Nspace, chi);
		}
		/* --- Bound-Free opacities due to ``metals'' --   -------------- */

		if (do_fudge) {
			Linear(Nfudge, lambda_fudge, fudge[2],
				   1, &lambda[nspect], &metal_fudge, FALSE);
		}
		for (n = 1; n < rh->atmos.Natom; n++) {
			Metal_bf(rh, lambda[nspect], 1, &rh->atmos.atoms[n], chi, eta);
			sprintf(metalID, "METAL_BF_%2s", rh->atmos.atoms[n].ID);
			if (do_fudge) {
				for (k = 0; k < rh->atmos.Nspace; k++) chi[k] *= metal_fudge;
			}
			write_contrib_xdr(&xdrs, metalID, type = ABSORPTION,
							  rh->atmos.Nspace, chi);
		}
		xdr_destroy(&xdrs);
		fclose(fp_out);
	}

	/* --- Clean up --                                     ------------ */

	for (n = 0; n < rh->atmos.Natom; n++) freeAtom(rh, rh->atmos.atoms + n);
	free(rh->atmos.atoms);

	for (n = 0; n < rh->atmos.Nmolecule; n++)
		freeMolecule(rh, rh->atmos.molecules + n);
	if (rh->atmos.Nmolecule > 0) free(rh->atmos.molecules);

	/* --- Free the temporary space set up by backgroundPrepare -- --- */

	backgroundCleanup(rh);

	free(chi);
	free(eta);
	free(scatt);
	free(thomson);

	if (do_fudge) {
		free(lambda_fudge);
		freeMatrix(rh, (void **) fudge);
	}
}
/* ------- end ---------------------------- backgrOpac.c ------------ */

/* ------- begin -------------------------- write_contrib_xdr.c ----- */

bool_t write_contrib_xdr(XDR *xdrs, char *label, enum opac_type type,
						 int Nspace, double *scatt)
{
	bool_t result = TRUE;

	result &= xdr_counted_string(xdrs, &label);
	result &= xdr_enum(xdrs, (enum_t *) &type);
	result &= xdr_vector(xdrs, (char *) scatt, Nspace,
						 sizeof(double), (xdrproc_t) rh_xdr_double);

	return result;
}
/* ------- end ---------------------------- write_contrib_xdr.c ----- */
