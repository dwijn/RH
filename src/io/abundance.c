/* ------- file: -------------------------- abundance.c -------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu May 20 13:52:41 2021 --

	   --------------------------                      ----------RH-- */

/* --- Defines the content of the elements structure of atmosphere atmos.
	   The periodic table index pt_index and atomic weight weight are
	   hard wired for all the elements. The abundance values are read
	   from the file defined by the keyword ABUND_FILE in keyword.input.
	   (see atom.h for a structure definition of struct Element.)

	   struct Element {
		 char     ID[ATOM_ID_WIDTH+1];
		 int     *mol_index, Nstage, Nmolecule;
		 double   weight, abund, *ionpot, **pf, double **n;
		 struct   Atom *model;
	   };

	   ID        -- Two-character atom ID
	   mol_index -- Index of molecules in which element can be bound
	   Nstage    -- Number of stages for which partition function is given
	   Nmolecule -- Size of array mol_index
	   weight    -- Atomic weight in atomic units
	   abund     -- Abundance ratio to hydrogen
	   ionpot    -- Array for ionization potential for successive stages
	   pf        -- Matrix with depth-dependent partition function
					 pf[Nstage][Nspace]
	   n         -- Matrix with populations for each stage n[Nstage][Nspace]
	   model     -- Pointer to atomic model if present in list of metals

	   If the Hydrogen abundance equals 12.0 the routine assumes that
	   abundances are given on a DEX scale. They are then converted to
	   base 10 numbers (i.e. A = 10^A - 12.0).

	   A metallicity factor is applied to elements other than hydrogen.
	   metallicity = POW10(input.metallicity).

	   Also evaluates the mean molecular weight avgWeight and stores it
	   in the Atmos structure.

 Note: List of elements and their atomic weights are taken from the
	   Atlas 9 code of R. L. Kurucz.
	   --                                              -------------- */

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "atomweights.h"
#include "state.h"

#define COMMENT_CHAR          "#"
#define MINIMUM_LOG_ABUNDANCE -8.0

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- readAbundance.c --------- */

void readAbundance(RHState *rh, Atmosphere *atmos)
{
	const char   routineName[] = "readAbundance";
	register int n, k, i;

	char     ID[ATOM_ID_WIDTH + 1], line[MAX_LINE_SIZE], *match;
	bool_t   result = TRUE, DEX = FALSE, exit_on_EOF;
	int      Nread, pti;
	double   abund, totalAbund, avgWeight, metallicity;
	Element *element;
	FILE    *fp_abund, *fp_pf;
	XDR      xdrs;

	atmos->Nelem    = sizeof(atomweight) / sizeof(struct AtomWeight);
	atmos->elements = (Element *) mallocOrAbort(atmos->Nelem * sizeof(Element));
	for (n = 0; n < atmos->Nelem; n++) {
		element = &atmos->elements[n];
		strcpy(element->ID, atomweight[n].ID);
		element->abundance_set = FALSE;
		element->abund         = 0.0;
		element->mol_index     = NULL;
		element->Nstage = element->Nmolecule = 0;
		element->weight                      = atomweight[n].weight;
		element->ionpot                      = NULL;
		element->pf                          = NULL;
		element->n                           = NULL;
		element->model                       = NULL;
	}

	if ((fp_abund = fopen(rh->input.abund_input, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to open input file %s", rh->input.abund_input);
	}
	/* --- Read abundances from file --                  -------------- */

	while (getLine(rh, fp_abund, COMMENT_CHAR, line, exit_on_EOF = FALSE) != EOF) {
		if ((Nread = sscanf(line, "%s %lf", ID, &abund)) != 2) {
			Error(rh, ERROR_LEVEL_2, routineName, "Unable to read input file %s",
				  rh->input.abund_input);
		}
		UpperCase(ID);
		if (strlen(ID) == 1) strcat(ID, " ");

		match = NULL;
		for (n = 0; n < atmos->Nelem; n++) {
			if ((match = strstr(atmos->elements[n].ID, ID))) {
				atmos->elements[n].abund = abund;
				if (strstr(ID, "H ") && (abund == 12.0)) DEX = TRUE;
				atmos->elements[n].abundance_set = TRUE;
				break;
			}
		}
		if (!match) {
			Warning(rh, routineName, "Abundance for element %s not used", ID);
		}
	}
	fclose(fp_abund);

	/* --- Multiply with metallicity factor if different from 1.0 -- -- */

	metallicity = POW10(rh->input.metallicity);
	if (rh->input.metallicity != 0.0) {
		Message(rh, "\nMultiplying metal abundances by metallicity of %5.3f\n\n",
				metallicity);
	}
	/* --- Open the data file with partition functions and first read the
		   temperature interpolation grid --             -------------- */

	if ((fp_pf = fopen(rh->input.pfData, "r")) == NULL) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to open input file %s for partition function data",
			  rh->input.pfData);
	}
	xdrstdio_create(&xdrs, fp_pf, XDR_DECODE);

	result &= xdr_int(&xdrs, &atmos->Npf);
	atmos->Tpf = (double *) mallocOrAbort(atmos->Npf * sizeof(double));
	result &= xdr_vector(&xdrs, (char *) atmos->Tpf, atmos->Npf,
						 sizeof(double), (xdrproc_t) rh_xdr_double);

	totalAbund = avgWeight = 0.0;
	for (n = 0; n < atmos->Nelem; n++) {
		element = atmos->elements + n;

		if (!element->abundance_set) {
			Warning(rh, routineName, "Found no abundance for element %s",
					element->ID);
		} else {

			/* --- Convert if abundances were given on logarithmic scale -- */

			if (DEX) element->abund = POW10(element->abund - 12.0);

			/* --- Apply metallicity factor to elements other than hydrogen */

			if (metallicity != 1.0 && !strstr(element->ID, "H "))
				element->abund *= metallicity;

			totalAbund += element->abund;
			avgWeight += element->abund * element->weight;
		}
		/* --- Read partition function information for all elements
			   to prevent complications in the PF input file processing - */

		result &= xdr_int(&xdrs, &pti);
		result &= xdr_int(&xdrs, &element->Nstage);
		element->pf = matrix_double(rh, element->Nstage, atmos->Npf);
		element->ionpot =
			(double *) mallocOrAbort(element->Nstage * sizeof(double));

		result &= xdr_vector(&xdrs, (char *) element->pf[0],
							 element->Nstage * atmos->Npf,
							 sizeof(double), (xdrproc_t) rh_xdr_double);
		result &= xdr_vector(&xdrs, (char *) element->ionpot,
							 element->Nstage,
							 sizeof(double), (xdrproc_t) rh_xdr_double);

		/* --- Store the logarithmic values of the partition functions
		   to facilitate logarithmic interpolation in temperature
		   for the calculation of population numbers. Do this only
			   if the abundance of this element is actually set -- ---- */

		if (element->abundance_set) {
			for (i = 0; i < element->Nstage; i++) {
				element->ionpot[i] *= (HPLANCK * CLIGHT) / CM_TO_M;
				for (k = 0; k < atmos->Npf; k++)
					element->pf[i][k] = log(element->pf[i][k]);
			}
		}
	}
	atmos->totalAbund = totalAbund;
	atmos->wght_per_H = avgWeight;
	atmos->avgMolWght = avgWeight / totalAbund;

	if (!result) {
		Error(rh, ERROR_LEVEL_2, routineName, "Unable to read partition function data from %s",
			  rh->input.pfData);
	}
	xdr_destroy(&xdrs);
	fclose(fp_pf);
}
/* ------- end ---------------------------- readAbundance.c --------- */
