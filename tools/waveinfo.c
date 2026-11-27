/* ------- file: -------------------------- waveinfo.c --------------

	   Version:       rh1.0
	   Author:        Han Uitenbroek (huitenbroek@sunspot.noao.edu)
	   Last modified: Mon Jul  8 14:28:37 2002 --

	   --------------------------                      ----------RH-- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "spectrum.h"
#include "xdr.h"

#define CHECK_OPTION(string, option, lmin) \
	(strncmp((string), (option), MAX((size_t) (lmin), strlen(string))) == 0)

struct Option {
	char *name;
	int   deflt;
	int   minlen;
	char *cmnt;
};

/* ------- begin -------------------------- ShowUsage.c ------------- */

void ShowUsage(char *programName, struct Option *theOpt)
{
	fprintf(stderr, " Usage:  %s [Options] Filename\n", programName);
	fprintf(stderr, "  where Options can be any of the following:\n\n");

	while (theOpt->name != NULL) {
		fprintf(stderr, "   %s \t<%s> (default: %d)\n",
				theOpt->name, theOpt->cmnt, theOpt->deflt);
		theOpt++;
	}
	fprintf(stderr,
			"\n  Each option can be abreviated to its unambiguous length.\n\n");

	exit(-1);
}
/* ------- end ---------------------------- ShowUsage.c ------------- */

/* ------- begin -------------------------- waveinfo.c -------------- */

int main(int argc, char *argv[])
{
	register int i;

	bool_t  showfull = FALSE, to_air = FALSE;
	int     Nlambda;
	char   *fileName;
	double *lambda;
	FILE   *fp_wave;
	XDR     xdrs;

	/* --- Reads wavelength input file and diplays contents -- -------- */

	struct Option theOpt[] = {
		{"-help", FALSE, 2, "Print this message          "},
		{"-full", FALSE, 2, "Print all wavelengths       "},
		{"-air", FALSE, 2, "Convert to wavelength in air"},
		{NULL, FALSE, 0, NULL}};
	if (argc <= 1) ShowUsage(argv[0], theOpt);

	for (i = 1; i < argc && argv[i][0] == '-'; i++) {

		if (CHECK_OPTION(argv[i], theOpt[0].name, theOpt[0].minlen)) {
			ShowUsage(argv[0], theOpt);
		} else if (CHECK_OPTION(argv[i], theOpt[1].name, theOpt[1].minlen)) {
			showfull = TRUE;
		} else if (CHECK_OPTION(argv[i], theOpt[2].name, theOpt[2].minlen)) {
			to_air = TRUE;
		} else {
			fprintf(stderr, " Invalid or ambiguous option: %s\n", argv[i]);
			ShowUsage(argv[0], theOpt);
		}
	}

	fileName = argv[argc - 1];
	if ((fp_wave = fopen(fileName, "r")) == NULL) {
		fprintf(stderr, " Error opening file: %s\n", fileName);
		exit(EXIT_FAILURE);
	}
	xdrstdio_create(&xdrs, fp_wave, XDR_DECODE);

	/* --- The number of wavelengths, then the wavelengths: a file that
		   ends before either, or gives a negative number, is not a
		   wavelength file (B-75) --                       -------------- */

	if (!xdr_int(&xdrs, &Nlambda)) {
		fprintf(stderr, " Error reading file: %s: no number of wavelengths\n",
				fileName);
		exit(EXIT_FAILURE);
	}
	if (Nlambda < 0) {
		fprintf(stderr, " Error reading file: %s: a negative number of "
						"wavelengths (%d)\n",
				fileName, Nlambda);
		exit(EXIT_FAILURE);
	}
	lambda = (double *) malloc(MAX(Nlambda, 1) * sizeof(double));
	if (lambda == NULL) {
		fprintf(stderr, " Error reading file: %s: no memory for %d wavelengths\n",
				fileName, Nlambda);
		exit(EXIT_FAILURE);
	}
	if (!xdr_vector(&xdrs, (char *) lambda, Nlambda, sizeof(double),
					(xdrproc_t) rh_xdr_double)) {
		fprintf(stderr, " Error reading file: %s: it ends before its %d "
						"wavelengths\n",
				fileName, Nlambda);
		exit(EXIT_FAILURE);
	}
	xdr_destroy(&xdrs);
	fclose(fp_wave);

	if (to_air) vacuum_to_air(Nlambda, lambda, lambda);

	if (Nlambda == 0) {
		fprintf(stderr, "\n Found no wavelengths\n\n");
	} else {
		fprintf(stderr,
				"\n Found %d wavelengths, ranging from %9.3f to %9.3f [nm]%s\n\n",
				Nlambda, lambda[0], lambda[Nlambda - 1],
				(to_air) ? " (in air)" : " (in vacuum)");
	}
	/* --- Five wavelengths to a row, each once (B-76: the last row was
		   always the last five) --                       -------------- */

	if (showfull && Nlambda > 0) {
		for (i = 0; i < Nlambda; i++) {
			printf("%9.3f", lambda[i]);
			if (i == Nlambda - 1)
				printf("\n\n");
			else
				printf("%s", (i % 5 == 4) ? "\n" : "   ");
		}
	}
	free(lambda);
	exit(EXIT_SUCCESS);
}
/* ------- end ---------------------------- waveinfo.c -------------- */
