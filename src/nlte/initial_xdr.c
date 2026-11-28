/* ------- file: -------------------------- initial_xdr.c -----------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Mon Sep 22 15:54:37 2025 --

	   --------------------------                      ----------RH-- */

/* --- Reads and/or computes the initial solution (populations and/or
	   mean intensity J).

	   XDR (external data representation) version.

	   Possible options:

		 LTE_POPULATIONS    -- Assume LTE populations initially.
		 ZERO_RADIATION     -- Solve statistical equilibrium with
							   zero radiation field
		 OLD_POPULATIONS    -- Read old populations from file
		 OLD_POPS_AND_J     -- Read both old populations and J from file
		 ESCAPE_PROBABILITY -- Not yet implemented
		 OLD_J              -- Use mean intensities from previous solution
							   (Only implemented for wavelength_table).

	   --                                              -------------- */

#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "accelerate.h"
#include "constant.h"
#include "statistics.h"
#include "error.h"
#include "inputs.h"
#include "xdr.h"
#include "prdh.h"
#include "state.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- initSolution.c ---------- */

rh_status initSolution(RHState *rh)
{
	const char   routineName[] = "initSolution";
	register int k, i, ij, nspect, n, nact, mu;

	char             permission[3];
	bool_t           openJfile;
	int              la, j, Nsr, Nplane, index, oflag;
	double           gijk, wla, twohnu3_c2, hc_k, twoc, fourPI, *J, *J20;
	Atom            *atom;
	Molecule        *molecule;
	ActiveSet       *as;
	AtomicLine      *line;
	AtomicContinuum *continuum;
	FILE            *fp;
	bool_t           openedImu = FALSE;
	rh_status        status    = RH_OK;

	getCPU(rh, 2, TIME_START, NULL);

	/* --- Allocate space for angle-averaged mean intensity -- -------- */

	if (!rh->input.limit_memory) {
		RH_TRY(allocMatrixDouble(rh, rh->spectrum.Nspect, rh->atmos.Nspace, &rh->spectrum.J));

		/* --- If we do background polarization we need space for the
			   anisotropy --                               -------------- */

		if (rh->input.backgr_pol)
			RH_TRY(allocMatrixDouble(rh, rh->spectrum.Nspect, rh->atmos.Nspace, &rh->spectrum.J20));
	}

	/* --- For the PRD angle approximation we need to store J in
		   the gas frame --                              -------------- */

	if (rh->input.PRD_angle_dep == PRD_ANGLE_APPROX &&
		rh->atmos.NPRDactive > 0) {

		RH_TRY(allocMatrixDouble(rh, rh->spectrum.Nspect, rh->atmos.Nspace, &rh->spectrum.Jgas));
		initJgasLocks(rh);
		RH_TRY(allocMatrixDouble(rh, rh->atmos.Nrays, rh->atmos.Nspace, &rh->spectrum.v_los));

		/* --- Calculate line of sight velocity --         -------------- */

		for (mu = 0; mu < rh->atmos.Nrays; mu++) {
			for (k = 0; k < rh->atmos.Nspace; k++) {
				rh->spectrum.v_los[mu][k] = vproject(rh, k, mu); // / vbroad[k];
			}
		}

		/* --- Interpolation weights of rho and Jgas, unless they are
			   computed when needed (PRDH_LIMIT_MEM) --  -------------- */

		RH_TRY(initPRDhTables(rh));
	}

	/* --- Allocate space for the emergent intensity --  -------------- */

	switch (rh->topology) {
	case ONE_D_PLANE:
		RH_TRY(allocMatrixDouble(rh, rh->spectrum.Nspect, rh->atmos.Nrays, &rh->spectrum.I));
		if (rh->atmos.Stokes || rh->input.backgr_pol) {
			RH_TRY(allocMatrixDouble(rh, rh->spectrum.Nspect, rh->atmos.Nrays, &rh->spectrum.Stokes_Q));
			RH_TRY(allocMatrixDouble(rh, rh->spectrum.Nspect, rh->atmos.Nrays, &rh->spectrum.Stokes_U));
			RH_TRY(allocMatrixDouble(rh, rh->spectrum.Nspect, rh->atmos.Nrays, &rh->spectrum.Stokes_V));
		}
		break;
	case TWO_D_PLANE:
		Nsr = rh->spectrum.Nspect * rh->atmos.Nrays;
		RH_TRY(allocMatrixDouble(rh, Nsr, rh->atmos.N[0], &rh->spectrum.I));
		if (rh->atmos.Stokes || rh->input.backgr_pol) {
			RH_TRY(allocMatrixDouble(rh, Nsr, rh->atmos.N[0], &rh->spectrum.Stokes_Q));
			RH_TRY(allocMatrixDouble(rh, Nsr, rh->atmos.N[0], &rh->spectrum.Stokes_U));
			RH_TRY(allocMatrixDouble(rh, Nsr, rh->atmos.N[0], &rh->spectrum.Stokes_V));
		}
		break;
	case THREE_D_PLANE:
		Nsr    = rh->spectrum.Nspect * rh->atmos.Nrays;
		Nplane = rh->atmos.N[0] * rh->atmos.N[1];
		RH_TRY(allocMatrixDouble(rh, Nsr, Nplane, &rh->spectrum.I));
		if (rh->atmos.Stokes || rh->input.backgr_pol) {
			RH_TRY(allocMatrixDouble(rh, Nsr, Nplane, &rh->spectrum.Stokes_Q));
			RH_TRY(allocMatrixDouble(rh, Nsr, Nplane, &rh->spectrum.Stokes_U));
			RH_TRY(allocMatrixDouble(rh, Nsr, Nplane, &rh->spectrum.Stokes_V));
		}
		break;
	default:
		RH_ASSERT(rh, FALSE, routineName, "Unknown topology (%d)", rh->topology);
	}
	/* --- Read angle-averaged intensity from previous run if necessary,
		   and open file for J in case option for limited memory is set.
		   From here on a failure closes the files opened here (failed,
		   below) --                                     -------------- */

	rh->spectrum.fd_J   = -1;
	rh->spectrum.fd_J20 = -1;
	oflag               = 0;
	openJfile           = FALSE;

	if (rh->input.startJ == OLD_J) {
		if (rh->spectrum.updateJ) {
			strcpy(permission, "r+");
			oflag |= O_RDWR;
		} else {
			strcpy(permission, "r");
			oflag |= O_RDONLY;
		}
		openJfile = TRUE;
	} else {
		if (rh->input.limit_memory) {
			strcpy(permission, "w+");
			oflag |= (O_RDWR | O_CREAT);
			openJfile = TRUE;
		}
	}
	if (openJfile) {
		if ((rh->spectrum.fd_J = open(rh->input.JFile, oflag, PERMISSIONS)) == -1) {
			status = Fail(rh, RH_E_DATA, routineName, "Unable to open input file %s with permission %s",
						  rh->input.JFile, permission);
			goto failed;
		}
		if (rh->input.backgr_pol) {
			if ((rh->spectrum.fd_J20 = open(J20_DOT_OUT, oflag,
											PERMISSIONS)) == -1) {
				status = Fail(rh, RH_E_DATA, routineName, "Unable to open input file %s with permission %s",
							  J20_DOT_OUT, permission);
				goto failed;
			}
		}
	}
	if (rh->input.limit_memory) {
		if (oflag & O_CREAT) {
			J = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));

			/* --- Initialize J file with zeroes --          -------------- */

			for (k = 0; k < rh->atmos.Nspace; k++) J[k] = 0.0;
			for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++)
				writeJlambda(rh, nspect, J);

			free(J);

			if (rh->input.backgr_pol) {
				J20 = (double *) mallocOrAbort(rh->atmos.Nspace * sizeof(double));
				for (k = 0; k < rh->atmos.Nspace; k++) J20[k] = 0.0;
				for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++)
					writeJ20lambda(rh, nspect, J20);

				free(J20);
			}
		}
	} else {
		if (rh->input.startJ == OLD_J) {

			/* --- Fill matrix J with old values from previous run ----- -- */

			for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++)
				readJlambda(rh, nspect, rh->spectrum.J[nspect]);

			close(rh->spectrum.fd_J);
			rh->spectrum.fd_J = -1;

			if (rh->input.backgr_pol) {
				for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++)
					readJ20lambda(rh, nspect, rh->spectrum.J20[nspect]);

				close(rh->spectrum.fd_J20);
				rh->spectrum.fd_J20 = -1;
			}
		}
	}

	/* --- Look for Jgas and read, otherwise use spectrum.J --- ------- */

	if (rh->atmos.NPRDactive > 0 &&
		rh->input.PRD_angle_dep == PRD_ANGLE_APPROX) {

		fp = fopen("Jgas.dat", "r");
		if (fp) {

			/* --- File exists --                            -------------- */

			fclose(fp);
			readJgas(rh, rh->spectrum.Jgas);
			Message(rh, "Read spectrum.Jgas from file.");
		} else {

			/* --- File does not exist --                    -------------- */

			Warning(rh, routineName, "Jgas.dat does not exist,setting spectrum.Jgas spectrum.J.");

			for (k = 0; k < rh->atmos.Nspace; k++) {
				for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++) {
					rh->spectrum.Jgas[nspect][k] = rh->spectrum.J[nspect][k];
				}
			}
		}
	}

	/* --- Need storage for angle-dependent specific intensities for
		   angle-dependent PRD --                        -------------- */

	if (rh->atmos.NPRDactive > 0 &&
		rh->input.PRD_angle_dep == PRD_ANGLE_DEP) {

		oflag = 0;
		if (rh->input.startJ == OLD_J) {
			if (rh->spectrum.updateJ) {
				strcpy(permission, "r+");
				oflag |= O_RDWR;
			} else {
				strcpy(permission, "r");
				oflag |= O_RDONLY;
			}
		} else {
			strcpy(permission, "w+");
			oflag |= (O_RDWR | O_CREAT);
		}
		if ((rh->spectrum.fd_Imu = open(IMU_FILENAME, oflag, PERMISSIONS)) == -1) {
			status = Fail(rh, RH_E_DATA, routineName, "Unable to open %s file %s with permission %s",
						  (rh->spectrum.updateJ) ? "update" : "input",
						  IMU_FILENAME, permission);
			goto failed;
		}
		openedImu = TRUE;
		/* --- Fill the index list that keeps track of the location
			   of intensity Imu in file spectrum.fd_Imu at wavelength
			   corresponding to nspect. --                 -------------- */

		rh->spectrum.PRDindex = (int *) mallocOrAbort(rh->spectrum.Nspect * sizeof(int));
		index                 = 0;
		for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++) {
			if (containsPRDline(rh, &rh->spectrum.as[nspect])) {
				rh->spectrum.PRDindex[nspect] = index;
				index++;
			}
		}
	}
	if (rh->input.PRD_angle_dep == PRD_ANGLE_APPROX && rh->atmos.NPRDactive > 0) {

		/* --- Which wavelengths the PRD sub-iterations re-solve, and the
			   rows of Jgas_fixed (N-27) --                  -------------- */

		char *isPRD = (char *) mallocOrAbort(rh->spectrum.Nspect);

		for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++)
			isPRD[nspect] = containsPRDline(rh, &rh->spectrum.as[nspect]);
		initJgasRows(rh, isPRD);
		free(isPRD);
	}
	for (nact = 0; nact < rh->atmos.Nactiveatom; nact++) {
		atom = rh->atmos.activeatoms[nact];

		/* --- Allocate memory for the rate equation matrix -- ---------- */

		RH_TRY_GOTO(status, failed,
					allocMatrixDouble(rh, SQ(atom->Nlevel), rh->atmos.Nspace, &atom->Gamma));

		switch (atom->initial_solution) {
		case LTE_POPULATIONS:
			for (i = 0; i < atom->Nlevel; i++) {
				for (k = 0; k < rh->atmos.Nspace; k++)
					atom->n[i][k] = atom->nstar[i][k];
			}
			break;

		case ZERO_RADIATION:
			hc_k   = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M);
			twoc   = 2.0 * CLIGHT / CUBE(NM_TO_M);
			fourPI = 4.0 * PI;

			initGammaAtom(rh, atom, 1);

			/* --- Then add radiative contributions of active transitions --  */

			for (nspect = 0; nspect < rh->spectrum.Nspect; nspect++) {
				as = rh->spectrum.as + nspect;

				for (n = 0; n < as->Nactiveatomrt[nact]; n++) {
					switch (as->art[nact][n].type) {
					case ATOMIC_LINE:
						line = as->art[nact][n].ptype.line;
						la   = nspect - line->Nblue;
						i    = line->i;
						j    = line->j;
						ij   = i * atom->Nlevel + j;

						if (la == 0) {
							for (k = 0; k < rh->atmos.Nspace; k++)
								atom->Gamma[ij][k] += line->Aji;
						}
						break;

					case ATOMIC_CONTINUUM:
						continuum = as->art[nact][n].ptype.continuum;
						la        = nspect - continuum->Nblue;
						i         = continuum->i;
						j         = continuum->j;
						ij        = i * atom->Nlevel + j;

						wla = fourPI * getwlambda_cont(continuum, la) /
							  continuum->lambda[la];
						twohnu3_c2 = twoc / CUBE(continuum->lambda[la]);
						for (k = 0; k < rh->atmos.Nspace; k++) {
							gijk = atom->nstar[i][k] / atom->nstar[j][k] *
								   exp(-hc_k / (continuum->lambda[la] * rh->atmos.T[k]));
							atom->Gamma[ij][k] += gijk * twohnu3_c2 *
												  continuum->alpha[la] * wla;
						}
						break;
					default:
						break;
					}
				}
			}
			/* --- Solve statistical equilibrium equations --  ------------ */

			RH_TRY_GOTO(status, failed,
						statEquil(rh, atom, (rh->input.isum == -1) ? 0 : rh->input.isum));
			break;

		case OLD_POPULATIONS:
			readPopulations(rh, atom);
			break;

		default:;
			break;
		}
		/* --- From here on the populations are the starting solution;
			   ChemicalEquilibrium() takes its neutral fraction from them
			   (N-36) --                                   -------------- */

		atom->popsInit = TRUE;
	}
	/* --- Now the molecules that are active --          -------------- */

	for (nact = 0; nact < rh->atmos.Nactivemol; nact++) {
		molecule = rh->atmos.activemols[nact];

		/* --- Calculate the LTE vibration level populations here. They
			   cannot be calculated yet in readMolecule since chemical
			   equilibrium has to be established first --  -------------- */

		for (i = 0; i < molecule->Nv; i++) {
			for (k = 0; k < rh->atmos.Nspace; k++)
				molecule->nvstar[i][k] = molecule->n[k] *
										 molecule->pfv[i][k] / molecule->pf[k];
		}
		/* --- Allocate memory for the rate equation matrix -- ---------- */

		RH_TRY_GOTO(status, failed,
					allocMatrixDouble(rh, SQ(molecule->Nv), rh->atmos.Nspace, &molecule->Gamma));

		switch (molecule->initial_solution) {

		case LTE_POPULATIONS:
			for (i = 0; i < molecule->Nv; i++) {
				for (k = 0; k < rh->atmos.Nspace; k++)
					molecule->nv[i][k] = molecule->nvstar[i][k];
			}
			break;

		case OLD_POPULATIONS:
			readMolPops(rh, molecule);
			break;

		default:;
		}

		/* --- Calculate collisions for molecule (must be done here because
			   rotation-vibration transitions are dominated by hydrogen and
			   H2 collisions for which chemical equilibrium needs to be
			   established first --                        -------------- */

		COcollisions(rh, molecule); /* only CO can be active (readMolecule, B-69) */
	}
	return RH_OK;

failed:
	/* --- The files opened here are closed; the arrays made here
		   belong to the run (phase 4) --                -------------- */

	if (rh->spectrum.fd_J >= 0) close(rh->spectrum.fd_J);
	if (rh->spectrum.fd_J20 >= 0) close(rh->spectrum.fd_J20);
	rh->spectrum.fd_J   = -1;
	rh->spectrum.fd_J20 = -1;
	if (openedImu) {
		close(rh->spectrum.fd_Imu);
		rh->spectrum.fd_Imu = -1;
	}
	return status;
}
/* ------- end ---------------------------- initSolution.c ---------- */
