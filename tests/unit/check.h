/* ------- file: -------------------------- check.h -----------------

	   Minimal unit-test support: checks that count failures and print
	   what went wrong, without stopping the test.

	   int main(void)
	   {
		   CHECK(n == 3, "n = %d", n);
		   CHECK_CLOSE(x, 1.0, 1e-12, 0.0, "x at k = %d", k);
		   return check_summary("my_test");
	   }

	   CHECK_CLOSE passes when |got - want| <= rtol * |want| + atol.
	   check_summary prints the number of checks and failures and
	   returns the exit status for ctest (0 = all passed).
	   --                                              -------------- */

#ifndef CHECK_H
#define CHECK_H

#include <math.h>
#include <stdarg.h>
#include <stdio.h>

static int check_count, check_failed;

/* --- Print at most this many failures; the count includes all -- */

#define CHECK_MAX_REPORTED 20

static inline void check_fail(const char *file, int line, const char *what,
							  const char *fmt, ...)
{
	va_list ap;

	if (++check_failed > CHECK_MAX_REPORTED) return;
	fprintf(stderr, "%s:%d: FAILED %s: ", file, line, what);
	va_start(ap, fmt);
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	fputc('\n', stderr);
}

#define CHECK(cond, ...)                                                 \
	do {                                                                 \
		check_count++;                                                   \
		if (!(cond)) check_fail(__FILE__, __LINE__, #cond, __VA_ARGS__); \
	} while (0)

static inline int check_close(double got, double want, double rtol, double atol)
{
	return fabs(got - want) <= rtol * fabs(want) + atol;
}

#define CHECK_CLOSE(got, want, rtol, atol, fmt, ...)                         \
	do {                                                                     \
		double check_g_ = (got), check_w_ = (want);                          \
		check_count++;                                                       \
		if (!check_close(check_g_, check_w_, (rtol), (atol)))                \
			check_fail(__FILE__, __LINE__, #got,                             \
					   "got %.17g, want %.17g (rel. diff %.3g); " fmt,       \
					   check_g_, check_w_,                                   \
					   (check_g_ - check_w_) / fabs(check_w_), __VA_ARGS__); \
	} while (0)

static inline int check_summary(const char *name)
{
	if (check_failed > CHECK_MAX_REPORTED)
		fprintf(stderr, "... %d more failures not shown\n",
				check_failed - CHECK_MAX_REPORTED);
	printf("%s: %d checks, %d failed\n", name, check_count, check_failed);
	return check_failed ? 1 : 0;
}

#endif /* !CHECK_H */
