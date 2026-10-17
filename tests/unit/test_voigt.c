/* ------- file: -------------------------- test_voigt.c ------------

	   Voigt and Faraday-Voigt functions against the Faddeeva function
	   w(v + ia) = H(a, v) + i L(a, v), computed with mpmath
	   (voigt_reference.h, from make_voigt_reference.py).

	   Error measures over the whole reference grid (a from 0 to 100,
	   v from 0 to 1e4):

	   -- rel:     |H - H_ref| / H_ref where H_ref >= 1e-12 H_ref(a, 0)
				   (a strong line's wing can matter against the continuum
				   down to about there);
	   -- peak:    |H - H_ref| / H_ref(a, 0), relative to line center;
	   -- faraday: |F - L_ref| / max |L_ref| for the given a.

	   Also: H even and F odd in v (exactly), finite values, and
	   continuity across the switch between the two approximations in
	   voigt.c at |v| + a = 8. Run with --table to print the errors per
	   damping parameter, with --bench to time Voigt().
	   --                                              -------------- */

#include <math.h>
#include <string.h>
#include <time.h>

#include "rh.h"
#include "check.h"
#include "voigt_reference.h"

#define REL_RANGE 1.0e-12

/* --- Bounds: about twice the maximum errors measured with gcc on
	   x86-64: rel 1.6e-6 (a ~ 1e-6, just outside the Doppler core),
	   peak 2.7e-10, faraday 9.5e-10 --                -------------- */

#define MAX_REL     3.5e-6
#define MAX_PEAK    6.0e-10
#define MAX_FARADAY 2.0e-9

static void check_reference(bool_t table)
{
	double max_rel = 0.0, max_peak = 0.0, max_faraday = 0.0;

	if (table)
		printf("%10s %10s %10s %10s %8s\n", "a", "rel", "peak", "faraday",
			   "worst v");

	for (int ia = 0; ia < VOIGT_REF_NA; ia++) {
		double a = voigt_ref_a[ia], H0 = voigt_ref_w[ia][0][0], Lmax = 0.0;
		double row_rel = 0.0, row_peak = 0.0, row_faraday = 0.0, worst_v = 0.0;

		for (int iv = 0; iv < VOIGT_REF_NV; iv++)
			Lmax = MAX(Lmax, fabs(voigt_ref_w[ia][iv][1]));

		for (int iv = 0; iv < VOIGT_REF_NV; iv++) {
			double v = voigt_ref_v[iv], F = 0.0, Fm = 0.0,
				   H    = Voigt(a, v, &F),
				   Hm   = Voigt(a, -v, &Fm),
				   Href = voigt_ref_w[ia][iv][0],
				   Lref = voigt_ref_w[ia][iv][1];
			double rel  = (Href >= REL_RANGE * H0) ? fabs(H - Href) / Href : 0.0,
				   peak = fabs(H - Href) / H0;

			CHECK(isfinite(H) && isfinite(F), "H(%g, %g) = %g, F = %g", a, v,
				  H, F);
			CHECK(Hm == H && Fm == -F, "H, F(%g, %g) = %.17g, %.17g; at -v: "
									   "%.17g, %.17g",
				  a, v, H, F, Hm, Fm);
			if (rel > row_rel) {
				row_rel = rel;
				worst_v = v;
			}
			row_peak = MAX(row_peak, peak);
			if (Lmax > 0.0)
				row_faraday = MAX(row_faraday, fabs(F - Lref) / Lmax);
		}
		if (table)
			printf("%10.3g %10.2e %10.2e %10.2e %8g\n", a, row_rel, row_peak,
				   row_faraday, worst_v);

		max_rel     = MAX(max_rel, row_rel);
		max_peak    = MAX(max_peak, row_peak);
		max_faraday = MAX(max_faraday, row_faraday);
	}
	printf("max rel %.2e  peak %.2e  faraday %.2e\n", max_rel, max_peak,
		   max_faraday);

	CHECK(max_rel <= MAX_REL, "relative error %.3g > %.3g", max_rel, MAX_REL);
	CHECK(max_peak <= MAX_PEAK, "error rel. to peak %.3g > %.3g", max_peak,
		  MAX_PEAK);
	CHECK(max_faraday <= MAX_FARADAY, "Faraday error %.3g > %.3g", max_faraday,
		  MAX_FARADAY);
}

/* --- The two approximations meet at |v| + a = 8: the jump there
	   must be within their accuracy (relative to H(a, 0) ~ 1 for
	   small a; the values there are ~0.07) --        -------------- */

static void check_switch(void)
{
	const double as[] = {0.0, 1.0e-4, 1.0e-2, 0.5, 2.0, 7.9};

	for (int n = 0; n < (int) (sizeof(as) / sizeof(as[0])); n++) {
		double a = as[n], v = 8.0 - a, Fl, Fr,
			   Hl = Voigt(a, nextafter(v, 0.0), &Fl),
			   Hr = Voigt(a, v, &Fr);

		CHECK(fabs(Hl - Hr) <= MAX_PEAK && fabs(Fl - Fr) <= MAX_FARADAY,
			  "jump at a = %g, v = %g: H %.3g, F %.3g", a, v, Hl - Hr, Fl - Fr);
	}
}

/* --- Time per call on a line-profile-like sample: a from 1e-4 to
	   1, v uniform in [-20, 20] --                    -------------- */

static void bench(void)
{
	enum { NA     = 16,
		   NV     = 4001,
		   REPEAT = 20 };
	double          sum = 0.0, F;
	struct timespec t0, t1;

	clock_gettime(CLOCK_MONOTONIC, &t0);
	for (int r = 0; r < REPEAT; r++)
		for (int ia = 0; ia < NA; ia++) {
			double a = pow(10.0, -4.0 + 4.0 * ia / (NA - 1));
			for (int iv = 0; iv < NV; iv++)
				sum += Voigt(a, -20.0 + 40.0 * iv / (NV - 1), &F) + F;
		}
	clock_gettime(CLOCK_MONOTONIC, &t1);
	printf("Voigt: %.1f ns per call\n",
		   1e9 * ((t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec)) /
			   ((double) REPEAT * NA * NV));
	if (sum == 0.0) printf("\n"); /* keep the loop */
}

int main(int argc, char *argv[])
{
	bool_t table = (argc > 1 && strcmp(argv[1], "--table") == 0);

	if (argc > 1 && strcmp(argv[1], "--bench") == 0) {
		bench();
		return 0;
	}
	check_reference(table);
	check_switch();

	return check_summary("test_voigt");
}
