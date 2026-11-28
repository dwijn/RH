/* ------- file: -------------------------- test_config.c -----------

	   The error paths of the configuration code (Q-56), in the process,
	   so that LeakSanitizer also checks that a failure closes and frees
	   what it holds:

	   -- getLine fails with RH_E_DATA at the end of the file (after
		  skipping comments and blank lines, as getLineOrEOF does, which
		  returns EOF there); checkNread fails when too few items were
		  read;
	   -- a keyword setter with a value outside its set returns
		  RH_E_INPUT (SOLVER, ANGLE_SET, S_INTERPOLATION), and a valid
		  one RH_OK with the value stored;
	   -- readValues: a line without a value, a required keyword not
		  set, and a value its setter rejects return RH_E_INPUT;
	   -- parse: an option without its required value returns
		  RH_E_INPUT; with the value, it is stored;
	   -- readInput: a keyword file that cannot be opened returns
		  RH_E_DATA; one with the required keywords and an invalid
		  value returns the setter's RH_E_INPUT, one with a value that
		  a check of readInput rejects (NG_ORDER < 0) that of
		  readInput (and the file is closed).
	   --                                              -------------- */

#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "error.h"
#include "check.h"
#include "state.h"

/* --- A temporary file with the given text; its name in path -- ---- */

static FILE *tempFile(char *path, const char *text)
{
	int   fd = mkstemp(path);
	FILE *fp = (fd >= 0) ? fdopen(fd, "w+") : NULL;

	CHECK(fp != NULL, "cannot write %s", path);
	if (fp == NULL) return NULL;
	fputs(text, fp);
	rewind(fp);
	return fp;
}

/* --- The recorded error is code from routine; clears it -- --------- */

static void expectError(rh_status got, rh_status code, const char *routine,
						const char *what)
{
	CHECK(got == code && Failed(rh) && rh->error.code == code,
		  "%s: status %d, recorded %d, want %d", what, (int) got,
		  (int) rh->error.code, (int) code);
	CHECK(!strcmp(rh->error.routine, routine), "%s: routine %s, want %s", what,
		  rh->error.routine, routine);
	ClearError(rh);
}

static void testGetLine(void)
{
	char  path[] = "test_config_XXXXXX", line[MAX_LINE_SIZE];
	FILE *fp     = tempFile(path, "# comment\n\n   \nfirst line\n#\n");

	if (fp == NULL) return;
	CHECK(getLineOrEOF(fp, "#", line) == (int) strlen("first line\n") &&
			  !strcmp(line, "first line\n"),
		  "getLineOrEOF: read \"%s\"", line);
	CHECK(getLineOrEOF(fp, "#", line) == EOF, "getLineOrEOF: no EOF at the end");
	expectError(getLine(rh, fp, "#", line), RH_E_DATA, "getLine", "getLine at the end");
	rewind(fp);
	CHECK_OK(getLine(rh, fp, "#", line));
	CHECK(!strcmp(line, "first line\n"), "getLine: read \"%s\"", line);
	fclose(fp);
	remove(path);

	expectError(checkNread(rh, 3, 4, "readSomething", 7), RH_E_DATA, "readSomething",
				"checkNread, 3 of 4");
	CHECK_OK(checkNread(rh, 4, 4, "readSomething", 7));
}

static void testSetters(void)
{
	enum solver_type solver = SOLVER_NG;
	enum S_interpol  interpolation;
	AngleSet         angleSet;

	expectError(setsolver(rh, "NEWTON", &solver), RH_E_INPUT, "setsolver", "SOLVER = NEWTON");
	CHECK(solver == SOLVER_NG, "SOLVER changed by a rejected value");
	CHECK_OK(setsolver(rh, "NEWTON_KRYLOV", &solver));
	CHECK(solver == SOLVER_NEWTON_KRYLOV, "SOLVER = NEWTON_KRYLOV not stored");

	expectError(setAngleSet(rh, "SET_GL_0X4", &angleSet), RH_E_INPUT, "setAngleSet",
				"ANGLE_SET = SET_GL_0X4");
	expectError(setAngleSet(rh, "SET_C4", &angleSet), RH_E_INPUT, "setAngleSet",
				"ANGLE_SET = SET_C4");
	CHECK_OK(setAngleSet(rh, "SET_GL_2X3", &angleSet));
	CHECK(angleSet.set == SET_GL && angleSet.Ninclination == 2 && angleSet.Nazimuth == 3,
		  "ANGLE_SET = SET_GL_2X3 not stored");

	expectError(set_S_Interpolation(rh, "S_PARABOLIC", &interpolation), RH_E_INPUT,
				"set_S_Interpolation", "S_INTERPOLATION = S_PARABOLIC");
	CHECK_OK(set_S_Interpolation(rh, "S_BEZIER3", &interpolation));
	CHECK(interpolation == S_BEZIER3, "S_INTERPOLATION = S_BEZIER3 not stored");
}

static void testReadValues(void)
{
	char             path[32];
	FILE            *fp;
	enum solver_type solver        = SOLVER_NG;
	int              Nrays         = 0;
	Keyword          theKeywords[] = {
        {"SOLVER", "NG", FALSE, KEYWORD_DEFAULT, &solver, setsolver},
        {"NRAYS", "", FALSE, KEYWORD_REQUIRED, &Nrays, setintValue}};

	strcpy(path, "test_config_XXXXXX");
	if ((fp = tempFile(path, "NRAYS\n")) == NULL) return;
	expectError(readValues(rh, fp, 2, theKeywords), RH_E_INPUT, "readValues",
				"readValues, a line without a value");
	fclose(fp);
	remove(path);

	strcpy(path, "test_config_XXXXXX");
	if ((fp = tempFile(path, "SOLVER = ANDERSON\n")) == NULL) return;
	expectError(readValues(rh, fp, 2, theKeywords), RH_E_INPUT, "readValues",
				"readValues, a required keyword not set");
	fclose(fp);
	remove(path);

	theKeywords[0].set = FALSE;
	strcpy(path, "test_config_XXXXXX");
	if ((fp = tempFile(path, "NRAYS = 3\nSOLVER = FOO\n")) == NULL) return;
	expectError(readValues(rh, fp, 2, theKeywords), RH_E_INPUT, "setsolver",
				"readValues, SOLVER = FOO");
	fclose(fp);
	remove(path);

	theKeywords[0].set = theKeywords[1].set = FALSE;
	strcpy(theKeywords[0].value, "NG");
	strcpy(path, "test_config_XXXXXX");
	if ((fp = tempFile(path, "NRAYS = 3\n# SOLVER = FOO\n")) == NULL) return;
	CHECK_OK(readValues(rh, fp, 2, theKeywords));
	CHECK(Nrays == 3 && solver == SOLVER_NG, "readValues: NRAYS = %d, SOLVER = %d", Nrays,
		  (int) solver);
	fclose(fp);
	remove(path);
}

static void testParse(void)
{
	char   input[MAX_VALUE_LENGTH] = "";
	Option theOptions[]            = {
        {"help", 1, FALSE, "", NULL, NULL, "Prints this message"},
        {"input", 1, TRUE, "keyword.input", input, setcharValue, "File name for input keywords"}};
	char *without[] = {"prog", "-input"}, *with[] = {"prog", "-input", "my.input"};

	expectError(parse(rh, 2, without, 2, theOptions), RH_E_INPUT, "parse",
				"parse, -input without a value");
	CHECK_OK(parse(rh, 3, with, 2, theOptions));
	CHECK(!strcmp(input, "my.input"), "parse: -input my.input gave \"%s\"", input);
}

/* --- The keywords that readInput requires, and NRAYS in 1-D -- ---- */

#define REQUIRED                                                   \
	"ATMOS_FILE = a\nABUND_FILE = b\nATOMS_FILE = c\n"             \
	"MOLECULES_FILE = d\nI_SUM = 0\nN_MAX_ITER = 10\nJ_FILE = J\n" \
	"BACKGROUND_FILE = bg\nSTARTING_J = NEW_J\n"                   \
	"KURUCZ_PF_DATA = pf\nVMICRO_CHAR = 2.0\nNRAYS = 3\n"

static void testReadInput(void)
{
	char  path[] = "test_config_XXXXXX", path2[] = "test_config_XXXXXX";
	FILE *fp;

	snprintf(rh->commandline.keyword_input, sizeof(rh->commandline.keyword_input), "%s",
			 "/nonexistent/keyword.input");
	expectError(readInput(rh), RH_E_DATA, "readInput", "readInput, missing keyword file");

	if ((fp = tempFile(path, REQUIRED "S_INTERPOLATION = S_PARABOLIC\n")) == NULL) return;
	fclose(fp);
	snprintf(rh->commandline.keyword_input, sizeof(rh->commandline.keyword_input), "%s", path);
	expectError(readInput(rh), RH_E_INPUT, "set_S_Interpolation",
				"readInput, S_INTERPOLATION = S_PARABOLIC");
	remove(path);

	if ((fp = tempFile(path2, REQUIRED "NG_ORDER = -1\n")) == NULL) return;
	fclose(fp);
	snprintf(rh->commandline.keyword_input, sizeof(rh->commandline.keyword_input), "%s", path2);
	CHECK(readInput(rh) == RH_E_INPUT && rh->error.message &&
			  strstr(rh->error.message, "NG_ORDER") != NULL,
		  "readInput, NG_ORDER = -1: %s", rh->error.message ? rh->error.message : "no error");
	expectError(RH_E_INPUT, RH_E_INPUT, "readInput", "readInput, NG_ORDER = -1");
	remove(path2);
}

int main(void)
{
	rh->commandline.logfile = stderr;
	rh->commandline.quiet   = TRUE;

	testGetLine();
	testSetters();
	testReadValues();
	testParse();
	testReadInput();

	return check_summary("test_config");
}
