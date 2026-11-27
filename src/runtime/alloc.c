/* ------- file: -------------------------- alloc.c -----------------

	   Allocations that cannot fail quietly (Q-53). The library makes
	   about 600 allocation calls; before Q-53 a dozen tested the
	   result, and a failed allocation ended the run in a segmentation
	   fault. These helpers, used through the macros of alloc.h
	   (mallocOrAbort, callocOrAbort, reallocOrAbort, strdupOrAbort),
	   print the routine, the size and the expression of a failed
	   allocation to stderr and abort. Phase 4 of the library track
	   replaces them by allocations that return RH_E_NOMEM, together
	   with the owners of the data and their free paths. The
	   declarations are in alloc.h.
	   --                                              -------------- */

#include <stdarg.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <stdio.h>

#include "alloc.h"

/* ------- begin -------------------------- outOfMemory.c ----------- */

void outOfMemory(const char *routineName, const char *format, ...)
{
	va_list ap;

	fflush(NULL);
	fprintf(stderr, "\a\n\n-OUT_OF_MEMORY in routine %s\n ", routineName);
	va_start(ap, format);
	vfprintf(stderr, format, ap);
	va_end(ap);
	fprintf(stderr, " \n Aborting.....\n");
	fflush(stderr);
	abort();
}
/* ------- end ---------------------------- outOfMemory.c ----------- */

/* ------- begin -------------------------- checkedMalloc.c --------- */

/* --- A zero size may give NULL: that is not a failure -- -------- */

void *checkedMalloc(size_t size, const char *routineName, const char *expression)
{
	void *p = malloc(size);

	if (p == NULL && size > 0)
		outOfMemory(routineName, "Could not allocate %zu bytes (%s)", size, expression);
	return p;
}

void *checkedCalloc(size_t n, size_t size, const char *routineName, const char *expression)
{
	void *p = calloc(n, size);

	if (p == NULL && n > 0 && size > 0)
		outOfMemory(routineName, "Could not allocate %zu x %zu bytes (%s)", n, size, expression);
	return p;
}

void *checkedRealloc(void *ptr, size_t size, const char *routineName, const char *expression)
{
	void *p = realloc(ptr, size);

	if (p == NULL && size > 0)
		outOfMemory(routineName, "Could not reallocate to %zu bytes (%s)", size, expression);
	return p;
}

char *checkedStrdup(const char *s, const char *routineName)
{
	char *p = strdup(s);

	if (p == NULL)
		outOfMemory(routineName, "Could not copy a string of %zu bytes", strlen(s) + 1);
	return p;
}
/* ------- end ---------------------------- checkedMalloc.c --------- */
