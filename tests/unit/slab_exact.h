/* ------- file: -------------------------- slab_exact.h ------------

	   Exact intensities in a slab whose source function is a polynomial
	   in optical depth, S(tau) = sum_n s_n tau^n (n < NPOLY), along a ray
	   with direction cosine mu:

		 outgoing (to_obs):  I+(tau) = sum_n s_n sum_j n!/(n-j)! mu^j tau^(n-j)

		 incoming, I-(tau0) = 0 at the top tau0:
		   I-(tau) = P-(tau) - exp(-(tau - tau0)/mu) P-(tau0),
		   P-(tau) = sum_n s_n sum_j n!/(n-j)! (-mu)^j tau^(n-j),

	   for a semi-infinite atmosphere, and also for a finite slab when
	   the lower boundary is irradiated with I+(tau_max). Used by the
	   formal-solver unit tests.
	   --                                              -------------- */

#ifndef SLAB_EXACT_H
#define SLAB_EXACT_H

#include <math.h>

#define NPOLY 4

/* --- P(tau) - P(tau0) for P(tau) = sum_n s[n] sum_j n!/(n-j)! m^j
	   tau^(n-j), without cancellation: tau^p - tau0^p is summed as
	   (tau - tau0) sum_i tau^(p-1-i) tau0^i. With tau0 = 0 this is
	   P(tau) - P(0) --                                -------------- */

static inline double P_diff(const double s[NPOLY], double tau, double tau0, double m)
{
	double P = 0.0;

	for (int n = 0; n < NPOLY; n++) {
		double fac = 1.0, sum = 0.0;

		for (int j = 0; j < n; j++) {
			int    p = n - j;
			double d = 0.0;

			for (int i = 0; i < p; i++)
				d += pow(tau, p - 1 - i) * pow(tau0, i);
			sum += fac * d;
			fac *= (n - j) * m;
		}
		P += s[n] * sum;
	}
	return (tau - tau0) * P;
}

/* --- P(0) = sum_n s[n] n! m^n --                    -------------- */

static inline double P_zero(const double s[NPOLY], double m)
{
	double P = 0.0, fac = 1.0;

	for (int n = 0; n < NPOLY; n++) {
		P += s[n] * fac;
		fac *= (n + 1) * m;
	}
	return P;
}

/* --- Exact intensity for S = sum_n s[n] tau^n, with I- = 0 at the
	   top of the slab, tau_top. I- is written as
	   P-(tau) - P-(tau0) - expm1(-(tau - tau0)/mu) P-(tau0), which has
	   no cancellation just below the top --           -------------- */

static inline double I_exact(const double s[NPOLY], double tau, double tau_top,
							 double mu, int to_obs)
{
	if (to_obs)
		return P_zero(s, mu) + P_diff(s, tau, 0.0, mu);
	else {
		double Ptop = P_zero(s, -mu) + P_diff(s, tau_top, 0.0, -mu);
		return P_diff(s, tau, tau_top, -mu) - expm1(-(tau - tau_top) / mu) * Ptop;
	}
}

#endif /* !SLAB_EXACT_H */
