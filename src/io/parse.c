/* ------- file: -------------------------- parse.c -----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Apr 22 09:00:26 2009 --

	   --------------------------                      ----------RH-- */

#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */

void ShowUsage(char *programName, int Noption, Option *theOptions);

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- parse.c ----------------- */

rh_status parse(RHState *rh, int argc, char *argv[], int Noption, Option *theOptions)
{
	const char   routineName[] = "parse";
	register int n, nopt;

	char    **options, **values;
	bool_t    recognized;
	int       Narg   = argc - 1, Nset;
	rh_status result = RH_OK;

	/* --- Break down the command line in option/value pairs -- ------- */

	options = (char **) mallocOrAbort(Narg * sizeof(char *));
	values  = (char **) mallocOrAbort(Narg * sizeof(char *));

	Nset = 0;
	for (n = 1; n <= Narg; n++) {
		if (argv[n][0] == '-') {
			options[Nset] = argv[n] + 1;
			if (n < Narg && argv[n + 1][0] != '-') {
				n++;
				values[Nset] = argv[n];
			} else
				values[Nset] = NULL;
			Nset++;
		}
	}
	/* --- If called with -help show options and exit -- -------------- */

	for (n = 0; n < Nset; n++) {
		if (CHECK_OPTION(options[n], "help",
						 theOptions[0].minlength)) {
			ShowUsage(argv[0], Noption, theOptions);
		}
	}
	/* --- Match the command line options with internal options -- ---- */

	for (n = 0; n < Nset; n++) {
		recognized = FALSE;
		for (nopt = 1; nopt < Noption; nopt++) {
			if (CHECK_OPTION(options[n], theOptions[nopt].name,
							 theOptions[nopt].minlength)) {
				recognized = TRUE;

				/* --- If a value was specified it is copied to the appropriate
					   tag of the Option structure --          -------------- */

				if (values[n] != NULL)
					strcpy(theOptions[nopt].value, values[n]);
				else {
					if (theOptions[nopt].value_required) {
						/* --- An error results if value is not given when required */

						result = Fail(rh, RH_E_INPUT, routineName, "Option %s requires value",
									  theOptions[nopt].name);
						goto cleanup;
					} else
						/* --- Options that require no value are assumed to be
								   boolean and are set to TRUE --      -------------- */

						strcpy(theOptions[nopt].value, "TRUE");
				}
			}
		}

		if (!recognized) {
			Warning(rh, routineName, "Unknown command line option %s", options[n]);
		}
	}
	/* --- Finally, set the internal options with either the
		   command line specified values or their defaults -- --------- */

	for (nopt = 1; nopt < Noption; nopt++) {
		if (theOptions[nopt].pointer != NULL)
			RH_TRY_GOTO(result, cleanup,
						(*theOptions[nopt].setValue)(rh, theOptions[nopt].value,
													 theOptions[nopt].pointer));
	}

cleanup:
	free(options);
	free(values);
	return result;
}
/* ------- end ---------------------------- parse.c ----------------- */

/* ---------------------------------------- ShowUsage.c ------------- */

void ShowUsage(char *programName, int Noption, Option *theOptions)
{
	register int nopt;

	fprintf(stderr, " Usage:  %s [Options]\n", programName);
	fprintf(stderr, "  where options can be any of the following:\n\n");

	/* --- Print out command line options, their default settings and
		   meaning --                                    -------------- */

	fprintf(stderr, "   -%s\n     [%s]\n\n", theOptions[0].name,
			theOptions[0].message);

	for (nopt = 1; nopt < Noption; nopt++)
		fprintf(stderr, "   -%s %s (default=%s)\n     [%s]\n\n",
				theOptions[nopt].name,
				(theOptions[nopt].value_required) ? "value" : "",
				theOptions[nopt].value,
				theOptions[nopt].message);

	fprintf(stderr,
			"\n  Each option can be abreviated to its unambiguous length.\n\n");

	exit(EXIT_SUCCESS);
}
/* ------- end ---------------------------- ShowUsage.c ------------- */
