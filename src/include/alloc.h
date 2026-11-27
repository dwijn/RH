/* ------- file: -------------------------- alloc.h -----------------

	   Allocations that cannot fail quietly (alloc.c, Q-53): on a
	   failure they print the routine, the size and the expression, and
	   abort. Phase 4 of the library track replaces them by allocations
	   that return RH_E_NOMEM. A zero size that gives NULL is not a
	   failure. No RH types, so that the generic numerics (anderson.c,
	   newton.c) can use them; rh.h includes it for everything else.
	   --                                              -------------- */

#ifndef __ALLOC_H__
#define __ALLOC_H__

#include <stddef.h>

_Noreturn void outOfMemory(const char *routineName, const char *format, ...)
	__attribute__((format(printf, 2, 3)));
void *checkedMalloc(size_t size, const char *routineName, const char *expression);
void *checkedCalloc(size_t n, size_t size, const char *routineName, const char *expression);
void *checkedRealloc(void *ptr, size_t size, const char *routineName, const char *expression);
char *checkedStrdup(const char *s, const char *routineName);

#define mallocOrAbort(size)       checkedMalloc((size), __func__, #size)
#define callocOrAbort(n, size)    checkedCalloc((n), (size), __func__, #n " x " #size)
#define reallocOrAbort(ptr, size) checkedRealloc((ptr), (size), __func__, #size)
#define strdupOrAbort(s)          checkedStrdup((s), __func__)

#endif /* !__ALLOC_H__ */

/* ------- end ---------------------------- alloc.h ----------------- */
