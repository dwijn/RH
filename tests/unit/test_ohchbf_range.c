/* ------- file: -------------------------- test_ohchbf_range.c -----

	   OH and CH bound-free opacity outside the temperature range of
	   the tables of Kurucz, van Dishoeck & Tarafdar (ohchbf.c, N-53).
	   The tables run from 2000 to 9000 K and the opacity used to be
	   set to 0 outside, while OH and CH exist from 1000 K. Now
	   log(sigma) is held at its value at the nearest edge:

	   -- the cross section kappa = chi / (n (1 - exp(-h nu / kT))) is
		  the same at every T below 2000 K as at 2000 K, and the same
		  above 9000 K as at 9000 K (the Planck factors follow the
		  actual T, so chi and eta themselves keep varying);
	   -- it is positive and finite there, and continuous across the
		  edges (2000 K against 2000 K - 1e-6 K and 9000 K + 1e-6 K);
	   -- the emissivity is chi times the Planck function
		  (eta = chi * B_nu(T)), also outside the range;
	   -- outside the energy range of the table (OH: 2.1 - 15 eV,
		  CH: 0.1 - 10.5 eV) the routine returns FALSE and the
		  opacity stays 0, as before.
	   --                                              -------------- */

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "constant.h"
#include "background.h"
#include "check.h"
#include "state.h"

#define NT 8

int main(void)
{
	static const double T[NT] = {500.0, 1000.0, 1500.0, 2000.0 - 1.0E-6, 2000.0,
								 9000.0, 9000.0 + 1.0E-6, 20000.0};
	Molecule            OH, CH;
	double              chi[NT], eta[NT], n[NT], kap[NT], lambda;
	int                 k, m;

	memset(&OH, 0, sizeof(OH));
	memset(&CH, 0, sizeof(CH));
	for (k = 0; k < NT; k++) n[k] = 1.0E15;
	OH.n             = n;
	CH.n             = n;
	rh->atmos.OH     = &OH;
	rh->atmos.CH     = &CH;
	rh->atmos.Nspace = NT;
	rh->atmos.T      = (double *) T;

	for (m = 0; m < 2; m++) {
		bool_t (*bf)(RHState *, double, double *, double *) = m ? CH_bf_opac : OH_bf_opac;
		const char *name                                    = m ? "CH" : "OH";

		for (lambda = 120.0; lambda < 560.0; lambda += 37.0) {
			double hnu = HPLANCK * CLIGHT / (lambda * NM_TO_M);

			if (!bf(rh, lambda, chi, eta)) continue;
			for (k = 0; k < NT; k++) {
				double stim = exp(-hnu / (KBOLTZMANN * T[k]));
				kap[k]      = chi[k] / (n[k] * (1.0 - stim));
				CHECK(kap[k] > 0.0 && isfinite(kap[k]),
					  "%s kappa positive and finite at T = %g K, lambda = %g nm", name,
					  T[k], lambda);
				/* eta = chi (2 h nu^3 / c^2) stim / (1 - stim) = chi B_nu */

				CHECK_CLOSE(eta[k],
							chi[k] * 2.0 * HPLANCK * CLIGHT / CUBE(lambda * NM_TO_M) *
								stim / (1.0 - stim),
							1.0E-12, 0.0, "%s eta = chi B_nu at T = %g K", name, T[k]);
			}
			/* --- Held at the edges, continuous across them --      */

			CHECK_CLOSE(kap[0], kap[4], 1.0E-12, 0.0, "%s kappa(500 K) = kappa(2000 K)", name);
			CHECK_CLOSE(kap[1], kap[4], 1.0E-12, 0.0, "%s kappa(1000 K) = kappa(2000 K)", name);
			CHECK_CLOSE(kap[2], kap[4], 1.0E-12, 0.0, "%s kappa(1500 K) = kappa(2000 K)", name);
			CHECK_CLOSE(kap[3], kap[4], 1.0E-6, 0.0, "%s kappa continuous at 2000 K", name);
			CHECK_CLOSE(kap[6], kap[5], 1.0E-6, 0.0, "%s kappa continuous at 9000 K", name);
			CHECK_CLOSE(kap[7], kap[5], 1.0E-12, 0.0, "%s kappa(20000 K) = kappa(9000 K)", name);
		}
		/* --- Outside the energy range: no opacity, returns FALSE -- */

		CHECK(!bf(rh, m ? 20000.0 : 700.0, chi, eta), "%s: lambda below/above the energy range", name);
		CHECK(!bf(rh, 60.0, chi, eta), "%s: lambda 60 nm is above 15 eV", name);
	}
	return check_summary("test_ohchbf_range");
}
