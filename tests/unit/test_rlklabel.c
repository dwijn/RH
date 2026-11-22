/* ------- file: -------------------------- test_rlklabel.c ---------

	   Kurucz level labels in LS coupling (RLKdet_level, kurucz.c,
	   RLK_EXPLICIT = FALSE):

	   -- "s6D)4d f7D" (the example in kurucz.c): term 7D, so S = 3 and
		  L = 2, outer electron d (l = 2);
	   -- labels that cannot be parsed ("4d f7d": no upper-case term
		  letter; "": no words) return FALSE;
	   -- getWords on labels with the most words for their length
		  ("a b c", "a"): all words found, no write past the array
		  (the sanitizers check that).

	   Under the sanitizers the run also checks that the word array of
	   getWords is freed on every path: before B-72 the FALSE returns
	   leaked it, which LeakSanitizer reports at exit.
	   --                                              -------------- */

#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "inputs.h"
#include "check.h"

extern InputData   input;
extern CommandLine commandline;

bool_t RLKdet_level(char *label, RLK_level *level);

int main(void)
{
	RLK_level level;
	char      ok[] = "s6D)4d f7D", lower[] = "4d f7d", empty[] = "";

	commandline.logfile = stderr;
	commandline.quiet   = TRUE;
	input.RLK_explicit  = FALSE;

	memset(&level, 0, sizeof(level));
	CHECK(RLKdet_level(ok, &level), "\"s6D)4d f7D\" not parsed");
	CHECK(level.S == 3.0 && level.L == 2 && level.l == 2 && level.cpl == LS_COUPLING,
		  "\"s6D)4d f7D\": S = %g, L = %d, l = %d, want 3, 2, 2", level.S, level.L, level.l);

	CHECK(!RLKdet_level(lower, &level), "\"4d f7d\" parsed");
	CHECK(!RLKdet_level(empty, &level), "empty label parsed");

	{
		char abc[] = "a b c", a[] = "a", **words;
		int  count;

		words = getWords(abc, " ", &count);
		CHECK(count == 3 && strcmp(words[2], "c") == 0, "\"a b c\": %d words", count);
		free(words);
		words = getWords(a, " ", &count);
		CHECK(count == 1 && strcmp(words[0], "a") == 0, "\"a\": %d words", count);
		free(words);
	}

	return check_summary("test_rlklabel");
}
