/* ------- file: -------------------------- getline.c ---------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Dec  9 11:43:18 1999 --

	   --------------------------                      ----------RH-- */

/* --- Routines for reading formatted input files. --  -------------- */

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "rh.h"
#include "error.h"
#include "inputs.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- getLineOrEOF.c ---------- */

int getLineOrEOF(FILE *inputFile, const char *commentChar, char *line)
{
	char *linePtr, dummy[2];

	/* --- Reads (into char array line) inputFile till first line
		   that does not start with commentChar and returns length of
		   line, or EOF at end of file. Empty lines are also ignored.
		   --                                            -------------- */

	while ((linePtr = fgets(line, MAX_LINE_SIZE, inputFile)) != NULL) {
		if ((sscanf(line, "%1s", dummy) > 0) && (*line != *commentChar))
			break;
	}
	if (linePtr == NULL) return EOF;
	return strlen(line);
}
/* ------- end ---------------------------- getLineOrEOF.c ---------- */

/* ------- begin -------------------------- getLine.c --------------- */

rh_status getLine(RHState *rh, FILE *inputFile, const char *commentChar, char *line)
{
	/* --- As getLineOrEOF, for a line that must be there: the end of
		   the file is an error (Q-56) --                -------------- */

	if (getLineOrEOF(inputFile, commentChar, line) == EOF)
		return Fail(rh, RH_E_DATA, "getLine", "Reached end of input file before all data was read");
	return RH_OK;
}
/* ------- end ---------------------------- getLine.c --------------- */

/* ------- begin -------------------------- checkNread.c ------------ */

rh_status checkNread(RHState *rh, int Nread, int Nrequired, const char *routineName,
					 int checkPoint)
{
	/* --- Check whether Nread igreater or equal to Nrequired and issue
		   error otherwise --                           --------------- */

	if (Nread < Nrequired) {
		return Fail(rh, RH_E_DATA, routineName, "Unable to read input file\n"
												" At checkpoint %d. Needed at least %d, read %d item%s",
					checkPoint, Nrequired, Nread, (Nread > 1) ? "s" : "");
	}
	return RH_OK;
}
/* ------- end ---------------------------- checkNread.c ------------ */

/* ------- end ---------------------------- UpperCase.c ------------- */

void UpperCase(char *string)
{
	register int n;

	/* --- Change string content to upper case --        -------------- */

	for (n = 0; n < (int) strlen(string); n++)
		string[n] = toupper(string[n]);
}
/* ------- end ---------------------------- UpperCase.c ------------- */
