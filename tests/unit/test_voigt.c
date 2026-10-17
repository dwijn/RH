/* ------- file: -------------------------- test_voigt.c ------------

	   Voigt and Faraday-Voigt functions of every algorithm against the
	   Faddeeva function w(v + ia) = H(a, v) + i L(a, v), computed with
	   mpmath (voigt_reference.h, from make_voigt_reference.py).

	   Two error measures per algorithm, over the whole reference grid:

	   -- rel:  |H - H_ref| / H_ref, the relative error of the line
				opacity at that frequency, where H_ref >= 1e-12 H_ref(a, 0)
				(a strong line's wing can matter against the continuum
				down to about there);
	   -- peak: |H - H_ref| / H_ref(a, 0), the error relative to line
				center, which is what matters when the wing is weak
				compared with the continuum.

	   The Faraday-Voigt function F returned by HUI_ETAL and HUMLICEK is
	   checked against L the same way, relative to max |L| for the
	   given a.

	   Run with --table to print the errors per damping parameter, with
	   --bench to time the algorithms.
	   --                                              -------------- */

#include <math.h>
#include <string.h>
#include <time.h>

#include "rh.h"
#include "check.h"
#include "voigt_reference.h"

#define REL_RANGE 1.0e-12

/* --- Accuracy of each algorithm over the reference grid -- -------- */

typedef struct {
	const char         *name;
	enum VoigtAlgorithm algorithm;
	bool_t              faraday;
	double              max_rel, max_peak, max_faraday; /* bounds */
} Algorithm;

/* --- Bounds: about twice the maximum errors measured with gcc on
	   x86-64. HUI_ETAL's relative error is not bounded: the rational
	   approximation has no Gaussian wings, so for small a it is wrong by
	   orders of magnitude where exp(-v^2) dominates (but small compared
	   with the peak) --                               -------------- */

static const Algorithm algorithms[] = {
	{"ARMSTRONG", ARMSTRONG, FALSE, 2.5e-4, 3.0e-6, 0.0},
	{"RYBICKI", RYBICKI, FALSE, 4.0e-8, 2.0e-10, 0.0},
	{"HUI_ETAL", HUI_ETAL, TRUE, HUGE_VAL, 8.0e-6, 1.3e-5},
	{"HUMLICEK", HUMLICEK, TRUE, 2.0e-4, 6.0e-5, 1.0e-4}};

#define N_ALGORITHMS (int) (sizeof(algorithms) / sizeof(algorithms[0]))

static void check_algorithm(const Algorithm *alg, bool_t table)
{
	double max_rel = 0.0, max_peak = 0.0, max_faraday = 0.0;

	if (table)
		printf("\n%s\n%10s %10s %10s %10s %8s\n", alg->name, "a",
			   "rel", "peak", "faraday", "worst v");

	for (int ia = 0; ia < VOIGT_REF_NA; ia++) {
		double a = voigt_ref_a[ia], H0 = voigt_ref_w[ia][0][0], Lmax = 0.0;
		double row_rel = 0.0, row_peak = 0.0, row_faraday = 0.0, worst_v = 0.0;

		for (int iv = 0; iv < VOIGT_REF_NV; iv++)
			Lmax = MAX(Lmax, fabs(voigt_ref_w[ia][iv][1]));

		for (int iv = 0; iv < VOIGT_REF_NV; iv++) {
			double v = voigt_ref_v[iv], F = 0.0,
				   H    = Voigt(a, v, alg->faraday ? &F : NULL, alg->algorithm),
				   Href = voigt_ref_w[ia][iv][0],
				   Lref = voigt_ref_w[ia][iv][1];
			double rel  = (Href >= REL_RANGE * H0) ? fabs(H - Href) / Href : 0.0,
				   peak = fabs(H - Href) / H0;

			CHECK(isfinite(H) && isfinite(F), "%s: H(%g, %g) = %g, F = %g",
				  alg->name, a, v, H, F);

			if (rel > row_rel) {
				row_rel = rel;
				worst_v = v;
			}
			row_peak = MAX(row_peak, peak);
			if (alg->faraday && Lmax > 0.0)
				row_faraday = MAX(row_faraday, fabs(F - Lref) / Lmax);

			/* --- H is even in v, F is odd (RYBICKI only to about
				   1e-10) --                                -------------- */

			double Fm = 0.0, Hm = Voigt(a, -v, alg->faraday ? &Fm : NULL,
										alg->algorithm);
			CHECK_CLOSE(Hm, H, 1.0e-9, 0.0, "%s: H(%g, -%g)", alg->name, a, v);
			if (alg->faraday)
				CHECK_CLOSE(-Fm, F, 1.0e-9, 0.0, "%s: F(%g, -%g)",
							alg->name, a, v);
		}
		if (table)
			printf("%10.3g %10.2e %10.2e %10.2e %8g\n", a, row_rel,
				   row_peak, row_faraday, worst_v);

		max_rel     = MAX(max_rel, row_rel);
		max_peak    = MAX(max_peak, row_peak);
		max_faraday = MAX(max_faraday, row_faraday);
	}
	printf("%-10s max rel %.2e  peak %.2e  faraday %.2e\n", alg->name,
		   max_rel, max_peak, max_faraday);

	CHECK(max_rel <= alg->max_rel, "%s: relative error %.3g > %.3g",
		  alg->name, max_rel, alg->max_rel);
	CHECK(max_peak <= alg->max_peak, "%s: error rel. to peak %.3g > %.3g",
		  alg->name, max_peak, alg->max_peak);
	CHECK(max_faraday <= alg->max_faraday, "%s: Faraday error %.3g > %.3g",
		  alg->name, max_faraday, alg->max_faraday);
}

/* --- Time per call on a line-profile-like sample: a from 1e-4 to
	   1, v uniform in [-20, 20] --                    -------------- */

static void bench(void)
{
	enum { NA     = 16,
		   NV     = 4001,
		   REPEAT = 20 };
	double          sum = 0.0;
	struct timespec t0, t1;

	for (int n = 0; n < N_ALGORITHMS; n++) {
		const Algorithm *alg = &algorithms[n];
		double           F;

		clock_gettime(CLOCK_MONOTONIC, &t0);
		for (int r = 0; r < REPEAT; r++)
			for (int ia = 0; ia < NA; ia++) {
				double a = pow(10.0, -4.0 + 4.0 * ia / (NA - 1));
				for (int iv = 0; iv < NV; iv++)
					sum += Voigt(a, -20.0 + 40.0 * iv / (NV - 1),
								 alg->faraday ? &F : NULL, alg->algorithm);
			}
		clock_gettime(CLOCK_MONOTONIC, &t1);
		printf("%-10s %6.1f ns per call\n", alg->name,
			   1e9 * ((t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec)) /
				   ((double) REPEAT * NA * NV));
	}
	if (sum == 0.0) printf("\n"); /* keep the loop */
}

int main(int argc, char *argv[])
{
	bool_t table = (argc > 1 && strcmp(argv[1], "--table") == 0);

	if (argc > 1 && strcmp(argv[1], "--bench") == 0) {
		bench();
		return 0;
	}

	for (int n = 0; n < N_ALGORITHMS; n++)
		check_algorithm(&algorithms[n], table);

	return check_summary("test_voigt");
}
