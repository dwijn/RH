/* ------- file: -------------------------- test_h2ff.c -------------

	   H2-minus and H2-plus free-free opacities beyond the last
	   tabulated wavelength (hydrogen.c: H2minus_ff, H2plus_ff; N-53).

	   -- H2-minus: kappa goes as lambda^2 beyond 15188.3 nm and joins
		  the table at that wavelength;
	   -- H2-plus: the last column (10000 nm) is held beyond it;
	   -- the table values are returned at the nodes, theta = 0.5
		  (T = 10080 K) for H2-minus and T = 2500 K for H2-plus;
	   -- the opacities are positive and finite at 20, 50, 100 micron.
	   --                                              -------------- */

#include <stdlib.h>
#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "background.h"
#include "constant.h"
#include "inputs.h"
#include "check.h"
#include "state.h"

#ifndef RH_SOURCE_DIR
#error "RH_SOURCE_DIR must be defined"
#endif

int main(void)
{
	enum { NSPACE = 2 };
	static double T[NSPACE]  = {THETA0 / 0.5, 2500.0},
				  ne[NSPACE] = {1.0E19, 1.0E19}, nH[NSPACE] = {1.0E22, 1.0E22},
				  np[NSPACE] = {1.0E18, 1.0E18}, nH2[NSPACE] = {1.0E21, 1.0E21},
				  chi[NSPACE], chi0[NSPACE], chip0[NSPACE];
	static double  *n[2] = {nH, np}, *nm[1] = {nH2};
	static Atom     H;
	static Molecule H2;
	double          lam, want, lambdas[] = {20000.0, 50000.0, 100000.0};
	int             i, k, nm2 = 0, nmp = 0;
	char            line[512];
	FILE           *log = tmpfile();

	H.n              = n;
	H.Nlevel         = 2;
	H2.n             = nm[0];
	rh->atmos.H      = &H;
	rh->atmos.H2     = &H2;
	rh->atmos.Nspace = NSPACE;
	rh->atmos.T      = T;
	rh->atmos.ne     = ne;
	strcpy(rh->input.GauntFfData, RH_SOURCE_DIR "/Atoms/gauntff.dat");
	strcpy(rh->input.HminusBfData, RH_SOURCE_DIR "/Atoms/Hminus_bf.ecsv");
	CHECK(log != NULL, "%s", "tmpfile");
	rh->commandline.logfile = log;
	hydrogenPrepare(rh);

	/* --- Table nodes (theta = 0.5 for depth 0, T = 2500 K for depth 1) --- */

	CHECK(H2minus_ff(rh, 15188.3, chi0), "%s", "H2minus_ff at the last node");
	want = nH2[0] * 1.0E-29 * ne[0] * KBOLTZMANN * T[0] * 71.6;
	CHECK_CLOSE(chi0[0], want, 1.0E-12, 0.0, "%s", "H2- node, theta = 0.5");
	CHECK(H2plus_ff(rh, 10000.0, chi), "%s", "H2plus_ff at the last node");
	want = nH[1] * 1.0E-29 * np[1] * 1.0E-20 * 1.39;
	CHECK_CLOSE(chi[1], want, 1.0E-12, 0.0, "%s", "H2+ node, T = 2500 K");

	/* --- Beyond the tables --                       -------------- */

	for (i = 0; i < 3; i++) {
		lam = lambdas[i];
		CHECK(H2minus_ff(rh, lam, chi), "H2minus_ff returns FALSE at %.0f nm", lam);
		for (k = 0; k < NSPACE; k++) {
			CHECK(chi[k] > 0.0 && isfinite(chi[k]), "H2- %g at %.0f nm", chi[k],
				  lam);
			CHECK_CLOSE(chi[k], chi0[k] * SQ(lam / 15188.3), 1.0E-12, 0.0,
						"H2- lambda^2 at %.0f nm, depth %d", lam, k);
		}
		CHECK(H2plus_ff(rh, lam, chi), "H2plus_ff returns FALSE at %.0f nm", lam);
		CHECK(H2plus_ff(rh, 10000.0, chip0), "%s", "H2plus_ff at the last node");
		for (k = 0; k < NSPACE; k++) {
			CHECK(chi[k] > 0.0 && isfinite(chi[k]), "H2+ %g at %.0f nm", chi[k],
				  lam);
			CHECK_CLOSE(chi[k], chip0[k], 1.0E-12, 0.0,
						"H2+ held at %.0f nm, depth %d", lam, k);
		}
	}

	/* --- Continuous at the last node --              -------------- */

	H2minus_ff(rh, 15188.3 * (1.0 - 1.0E-9), chi);
	H2minus_ff(rh, 15188.3, chi0);
	for (k = 0; k < NSPACE; k++)
		CHECK_CLOSE(chi[k], chi0[k], 1.0E-7, 0.0, "H2- join, depth %d", k);
	H2plus_ff(rh, 10000.0 * (1.0 - 1.0E-9), chi);
	H2plus_ff(rh, 10000.0, chi0);
	for (k = 0; k < NSPACE; k++)
		CHECK_CLOSE(chi[k], chi0[k], 1.0E-7, 0.0, "H2+ join, depth %d", k);

	/* --- One message per routine per run --          -------------- */

	rewind(log);
	while (fgets(line, sizeof(line), log)) {
		if (strstr(line, "WARNING in routine H2minus_ff")) nm2++;
		if (strstr(line, "WARNING in routine H2plus_ff")) nmp++;
	}
	CHECK(nm2 == 1, "H2minus_ff messages: %d", nm2);
	CHECK(nmp == 1, "H2plus_ff messages: %d", nmp);
	fclose(log);

	hydrogenCleanup(rh);
	return check_summary("test_h2ff");
}
/* ------- end ---------------------------- test_h2ff.c ------------- */
