/* ------- file: -------------------------- readinput.c -------------

	   Version:       rh2.0, 1-D plane-parallel
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Mon Sep 22 16:09:49 2025 --

	   --------------------------                      ----------RH-- */

/* --- Reads input data for and defines keywords for 1-D
	   plane-parallel version --                       -------------- */

#include <string.h>

#include "rh.h"
#include "atom.h"
#include "atmos.h"
#include "spectrum.h"
#include "constant.h"
#include "statistics.h"
#include "inputs.h"
#include "error.h"

/* --- Global variables --                             -------------- */

extern enum Topology topology;

extern Atmosphere   atmos;
extern Spectrum     spectrum;
extern InputData    input;
extern CommandLine  commandline;
extern ProgramStats stats;
extern char         messageStr[];

/* --- Function prototypes --                          -------------- */

static bool_t keywordSet(const Keyword *theKeywords, int Nkeyword,
						 const char *name);
static void   setKeywordValue(Keyword *theKeywords, int Nkeyword,
							  const char *name, const char *value);

/* ------- begin -------------------------- readInput.c ------------- */

void readInput()
{
	const char routineName[] = "readInput";
	char       errorStr[MAX_MESSAGE_LENGTH];

	int              Nkeyword, Nconflict, n, p;
	char             sname[MAX_KEYWORD_LENGTH + 16];
	const char      *prefix;
	enum solver_type owner;
	FILE            *fp_keyword;

	Keyword theKeywords[] = {
		{"ATMOS_FILE", "", FALSE, KEYWORD_REQUIRED, input.atmos_input,
		 setcharValue},
		{"ABUND_FILE", "", FALSE, KEYWORD_REQUIRED, input.abund_input,
		 setcharValue},

		{"NRAYS", "0", FALSE, KEYWORD_OPTIONAL, &atmos.Nrays, setintValue},
		{"ANGLE_SET", "NO_SET", FALSE, KEYWORD_OPTIONAL, &atmos.angleSet,
		 setAngleSet},

		{"ATMOS_ITOP", "none", FALSE, KEYWORD_OPTIONAL, input.Itop, setcharValue},

		{"WAVETABLE", "none", FALSE, KEYWORD_OPTIONAL, input.wavetable_input,
		 setcharValue},
		{"ATOMS_FILE", "", FALSE, KEYWORD_REQUIRED, input.atoms_input,
		 setcharValue},
		{"MOLECULES_FILE", "", FALSE, KEYWORD_REQUIRED, input.molecules_input,
		 setcharValue},
		{"NON_ICE", "FALSE", FALSE, KEYWORD_OPTIONAL, &input.NonICE,
		 setboolValue},

		{"N_MAX_SCATTER", "5", FALSE, KEYWORD_OPTIONAL, &input.NmaxScatter,
		 setintValue},

		{"I_SUM", "0", FALSE, KEYWORD_REQUIRED, &input.isum, setintValue},
		{"N_MAX_ITER", "", FALSE, KEYWORD_REQUIRED, &input.NmaxIter,
		 setintValue},
		{"ITER_LIMIT", "", FALSE, KEYWORD_REQUIRED, &input.iterLimit,
		 setdoubleValue},
		{"NG_DELAY", "0", FALSE, KEYWORD_OPTIONAL, &input.Ngdelay, setintValue},
		{"NG_ORDER", "0", FALSE, KEYWORD_OPTIONAL, &input.Ngorder, setintValue},
		{"NG_PERIOD", "1", FALSE, KEYWORD_OPTIONAL, &input.Ngperiod,
		 setintValue},
		{"NG_MOLECULES", "FALSE", FALSE, KEYWORD_DEFAULT, &input.accelerate_mols,
		 setboolValue},
		{"COLL_RELAX_NSTEP", "0", FALSE, KEYWORD_OPTIONAL, &input.CR_Nstep,
		 setintValue},
		{"COLL_RELAX_FACTOR", "1.0", FALSE, KEYWORD_OPTIONAL, &input.CR_factor,
		 setdoubleValue},

		{"PRD_N_MAX_ITER", "3", FALSE, KEYWORD_OPTIONAL, &input.PRD_NmaxIter,
		 setintValue},
		{"PRD_ITER_LIMIT", "1.0E-2", FALSE, KEYWORD_OPTIONAL, &input.PRDiterLimit,
		 setdoubleValue},
		{"PRD_NG_DELAY", "0", FALSE, KEYWORD_OPTIONAL, &input.PRD_Ngdelay,
		 setintValue},
		{"PRD_NG_ORDER", "0", FALSE, KEYWORD_OPTIONAL, &input.PRD_Ngorder,
		 setintValue},
		{"PRD_NG_PERIOD", "1", FALSE, KEYWORD_OPTIONAL, &input.PRD_Ngperiod,
		 setintValue},
		{"PRD_ANGLE_DEP", "FALSE", FALSE, KEYWORD_DEFAULT,
		 &input.PRD_angle_dep, setPRDangle},
		{"PRDH_LIMIT_MEM", "FALSE", FALSE, KEYWORD_OPTIONAL, &input.prdh_limit_mem,
		 setboolValue},
		{"XRD", "FALSE", FALSE, KEYWORD_DEFAULT, &input.XRD, setboolValue},

		{"J_FILE", "", FALSE, KEYWORD_REQUIRED, input.JFile, setcharValue},
		{"BACKGROUND_FILE", "", FALSE, KEYWORD_REQUIRED, input.background_File,
		 setcharValue},
		{"BACKGROUND_RAY_FILE", "background.ray", FALSE,
		 KEYWORD_OPTIONAL, input.background_ray_File,
		 setcharValue},
		{"OLD_BACKGROUND", "FALSE", FALSE, KEYWORD_OPTIONAL,
		 &input.old_background, setboolValue},
		{"BACKGROUND_DOUBLE", "FALSE", FALSE, KEYWORD_DEFAULT,
		 &input.background_double, setboolValue},
		{"STARTING_J", "", FALSE, KEYWORD_REQUIRED, &input.startJ,
		 setstartValue},
		{"HYDROGEN_LTE", "FALSE", FALSE, KEYWORD_DEFAULT, &atmos.H_LTE,
		 setboolValue},
		{"HYDROSTATIC", "FALSE", FALSE, KEYWORD_DEFAULT, &atmos.hydrostatic,
		 setboolValue},
		{"KURUCZ_DATA", "none", FALSE, KEYWORD_OPTIONAL, &input.KuruczData,
		 setcharValue},
		{"RLK_SCATTER", "FALSE", FALSE, KEYWORD_DEFAULT, &input.rlkscatter,
		 setboolValue},
		{"KURUCZ_PF_DATA", "../../Atoms/pf_Kurucz.input", FALSE,
		 KEYWORD_REQUIRED, &input.pfData, setcharValue},
		{"SOLVE_NE", "NONE", FALSE, KEYWORD_DEFAULT, &input.solve_ne,
		 setnesolution},
		/* --- NG is the default in 2-D and 3-D; ANDERSON in 1-D and
			   spherical geometry, set below (S-10) -- -- */
		{"SOLVER", "NG", FALSE, KEYWORD_DEFAULT, &input.solver, setsolver},
		{"ANDERSON_DEPTH", "10", FALSE, KEYWORD_DEFAULT, &input.Anderson_depth,
		 setintValue},
		{"ANDERSON_DELAY", "2", FALSE, KEYWORD_DEFAULT, &input.Anderson_delay,
		 setintValue},
		{"ANDERSON_COND_MAX", "1.0E+08", FALSE, KEYWORD_DEFAULT,
		 &input.Anderson_cond_max, setdoubleValue},
		{"ANDERSON_SAFEGUARD", "1.0", FALSE, KEYWORD_DEFAULT,
		 &input.Anderson_safeguard, setdoubleValue},
		{"ANDERSON_VARIABLES", "LOG", FALSE, KEYWORD_DEFAULT,
		 &input.Anderson_variables, setSolverVariables},
		{"ANDERSON_COEFFICIENTS", "GLOBAL", FALSE, KEYWORD_DEFAULT,
		 &input.Anderson_coefficients, setAndersonCoefficients},
		{"ANDERSON_EXTRAPOLATE", "ALL", FALSE, KEYWORD_DEFAULT,
		 &input.Anderson_extrapolate, setAndersonExtrapolate},
		{"NEWTON_NOISE", "1.0E-10", FALSE, KEYWORD_DEFAULT, &input.Newton_noise,
		 setdoubleValue},
		{"NEWTON_KRYLOV_MAX", "20", FALSE, KEYWORD_DEFAULT,
		 &input.Newton_krylov_max, setintValue},
		{"NEWTON_ETA_MAX", "0.9", FALSE, KEYWORD_DEFAULT, &input.Newton_eta_max,
		 setdoubleValue},
		{"NEWTON_VARIABLES", "LOG", FALSE, KEYWORD_DEFAULT,
		 &input.Newton_variables, setSolverVariables},
		{"OPACITY_FUDGE", "none", FALSE, KEYWORD_OPTIONAL, &input.fudgeData,
		 setcharValue},
		{"METALLICITY", "0.0", FALSE, KEYWORD_DEFAULT, &input.metallicity,
		 setdoubleValue},

		{"ATMOS_OUTPUT", "atmos.out", FALSE, KEYWORD_DEFAULT, input.atmos_output,
		 setcharValue},
		{"GEOMETRY_OUTPUT", "geometry.out", FALSE, KEYWORD_OPTIONAL,
		 input.geometry_output, setcharValue},
		{"SPECTRUM_OUTPUT", "spectrum.out", FALSE, KEYWORD_DEFAULT,
		 input.spectrum_output, setcharValue},
		{"OPACITY_OUTPUT", "none", FALSE, KEYWORD_OPTIONAL, input.opac_output,
		 setcharValue},
		{"RADRATE_OUTPUT", "none", FALSE, KEYWORD_OPTIONAL, input.radrateFile,
		 setcharValue},
		{"COLLRATE_OUTPUT", "none", FALSE, KEYWORD_OPTIONAL, input.collrateFile,
		 setcharValue},
		{"DAMPING_OUTPUT", "none", FALSE, KEYWORD_OPTIONAL, input.dampingFile,
		 setcharValue},
		{"COOLING_OUTPUT", "none", FALSE, KEYWORD_OPTIONAL, input.coolingFile,
		 setcharValue},

		{"VMICRO_CHAR", "", FALSE, KEYWORD_REQUIRED, &atmos.vmicro_char,
		 setdoubleValue},
		{"VMACRO_TRESH", "0.0", FALSE, KEYWORD_OPTIONAL, &atmos.vmacro_tresh,
		 setdoubleValue},
		{"LAMBDA_REF", "500.0", FALSE, KEYWORD_DEFAULT, &atmos.lambda_ref,
		 setdoubleValue},
		{"VACUUM_TO_AIR", "0", FALSE, KEYWORD_OPTIONAL, &spectrum.vacuum_to_air,
		 setboolValue},

		/* --- Magnetic field related inputs go here --     ------------- */

		{"STOKES_INPUT", "none", FALSE, KEYWORD_OPTIONAL, input.Stokes_input,
		 setcharValue},
		{"B_STRENGTH_CHAR", "0.0", FALSE, KEYWORD_DEFAULT, &atmos.B_char,
		 setdoubleValue},
		{"STOKES_MODE", "NO_STOKES", FALSE, KEYWORD_OPTIONAL,
		 &input.StokesMode, setStokesMode},
		{"MAGNETO_OPTICAL", "TRUE", FALSE, KEYWORD_DEFAULT,
		 &input.magneto_optical, setboolValue},
		{"BACKGROUND_POLARIZATION", "FALSE", FALSE, KEYWORD_DEFAULT,
		 &input.backgr_pol, setboolValue},
		{"RLK_EXPLICIT", "FALSE", FALSE, KEYWORD_DEFAULT,
		 &input.RLK_explicit, setboolValue},

		{"XDR_ENDIAN", "TRUE", FALSE, KEYWORD_OPTIONAL,
		 &input.xdr_endian, setboolValue},

		{"S_INTERPOLATION", "S_BEZIER3", FALSE, KEYWORD_DEFAULT,
		 &input.S_interpolation, set_S_Interpolation},
		{"S_INTERPOLATION_STOKES", "DELO_BEZIER3", FALSE, KEYWORD_DEFAULT,
		 &input.S_interpolation_stokes, set_S_interpolation_stokes},

		{"INTERPOLATE_2D", "BEZIER_2D", FALSE, KEYWORD_DEFAULT,
		 &input.interpolate_2D, setInterpolate_2D},
		{"INTERPOLATE_3D", "BICUBIC_3D", FALSE, KEYWORD_DEFAULT,
		 &input.interpolate_3D, setInterpolate_3D},

		{"PRINT_CPU", "0", FALSE, KEYWORD_OPTIONAL, &stats.printCPU,
		 setboolValue},
		{"N_THREADS", "0", FALSE, KEYWORD_OPTIONAL, &input.Nthreads,
		 setThreadValue},

		{"LIMIT_MEMORY", "FALSE", FALSE, KEYWORD_DEFAULT, &input.limit_memory,
		 setboolValue},
		{"ALLOW_PASSIVE_BB", "TRUE", FALSE, KEYWORD_DEFAULT,
		 &input.allow_passive_bb, setboolValue}};
	Nkeyword = sizeof(theKeywords) / sizeof(Keyword);

	/* --- The default solver depends on the geometry (S-10): ANDERSON
		   in 1-D and spherical geometry, NG in 2-D and 3-D, where
		   Anderson's history does not fit in memory yet (S-11) -- -- */

	if (topology == ONE_D_PLANE || topology == SPHERICAL_SYMMETRIC)
		setKeywordValue(theKeywords, Nkeyword, "SOLVER", "ANDERSON");

	/* --- Open the input data file --                    ------------- */

	if ((fp_keyword = fopen(commandline.keyword_input, "r")) == NULL) {
		sprintf(messageStr, "Unable to open inputfile %s",
				commandline.keyword_input);
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}

	readValues(fp_keyword, Nkeyword, theKeywords);
	fclose(fp_keyword);

	/* --- Perform some sanity checks --                 -------------- */

	switch (topology) {
	case ONE_D_PLANE:
		if (atmos.Nrays == 0)
			Error(ERROR_LEVEL_2, routineName,
				  "Must set keyword NRAYS in 1-D plane parallel geometry");

		if (atmos.angleSet.set == SET_VERTICAL ||
			atmos.angleSet.set == SET_GAUSS_LOBATTO ||
			atmos.angleSet.set == SET_EDDINGTON) {
			switch (atmos.angleSet.set) {
			case SET_VERTICAL:

				if (atmos.Nrays != 1) {
					Error(WARNING, routineName,
						  "Ignoring value of NRAYS in 1-D plane geometry when"
						  " SET_VERTICAL is specified\n. Setting NRAYS = 1");
				}
				atmos.Nrays = 1;
				break;

			case SET_GAUSS_LOBATTO:

				if (atmos.Nrays != NRO_GLOB_4 && atmos.Nrays != NRO_GLOB_6) {
					sprintf(errorStr, "NRAYS must be %d or %d for SET_GAUSS_LOBATTO",
							NRO_GLOB_4, NRO_GLOB_6);
					Error(ERROR_LEVEL_2, routineName, errorStr);
				}
				break;

			case SET_EDDINGTON:

				if (atmos.Nrays != 1) {
					Error(WARNING, routineName,
						  "Ignoring value of NRAYS in 1-D plane geometry when"
						  " SET_EDDINGTON is specified\n. Setting NRAYS = 1,"
						  " and muz = 1/sqrt(3)");
				}
				atmos.Nrays = 1;
				break;

			default:
				break;
			}
		} else {
			Error(WARNING, routineName,
				  "Ignoring value of keyword ANGLE_SET in 1-D plane geometry");
		}
		break;

	case SPHERICAL_SYMMETRIC:
		if (atmos.Nrays > 0 || atmos.angleSet.set != NO_SET) {
			Error(WARNING, routineName,
				  "Ignoring value of keywords ANGLE_SET and NRAYS in "
				  "spherical geometry");
		}
		break;

	case TWO_D_PLANE:
		if (input.Eddington && atmos.angleSet.set != NO_SET) {
			Error(WARNING, routineName,
				  "Ignoring value of keywords ANGLE_SET > NO_VERTICAL when\n "
				  " EDDINGTON is set to TRUE\n"
				  " Using SET_VERTICAL with muz = 1/sqrt(3)");
			atmos.angleSet.set = SET_EDDINGTON;
		}
		break;

	case THREE_D_PLANE:
		if (atmos.angleSet.set == NO_SET)
			Error(ERROR_LEVEL_2, routineName,
				  "Must set keyword ANGLE_SET in multi-D plane geometry");
		if (atmos.Nrays > 0)
			Error(WARNING, routineName,
				  "Ignoring value of keyword NRAYS in multi-D plane geometry");
		break;
	}

	/* --- Ng orders and delays are counts (B-45: with order + 2 <= 0
		   the history allocation is corrupted or, at -2, Accelerate
		   divides by zero; a negative delay means nothing), and Ng
		   acceleration is applied every N_PERIOD iterations, so the
		   period must be at least 1 --                  -------------- */

	if (input.Ngorder < 0 || input.Ngdelay < 0)
		Error(ERROR_LEVEL_2, routineName,
			  "Values of NG_ORDER and NG_DELAY should be at least 0");
	if (input.PRD_Ngorder < 0 || input.PRD_Ngdelay < 0)
		Error(ERROR_LEVEL_2, routineName,
			  "Values of PRD_NG_ORDER and PRD_NG_DELAY should be at least 0");
	if (input.Ngorder > 0 && input.Ngperiod < 1)
		Error(ERROR_LEVEL_2, routineName,
			  "Value of NG_PERIOD should be at least 1 when NG_ORDER > 0");
	if (input.PRD_Ngorder > 0 && input.PRD_Ngperiod < 1)
		Error(ERROR_LEVEL_2, routineName,
			  "Value of PRD_NG_PERIOD should be at least 1 when PRD_NG_ORDER > 0");

	/* --- Lambda_ref should be positive since it is used for the conversion
		   of depth scales in routine convertScales --    ------------- */

	switch (topology) {
	case ONE_D_PLANE:
	case SPHERICAL_SYMMETRIC:
		if (atmos.lambda_ref <= 0.0)
			Error(ERROR_LEVEL_2, routineName,
				  "Value of LAMBDA_REF should be larger than 0.0");
		break;
	default:
		if (atmos.lambda_ref < 0.0)
			Error(ERROR_LEVEL_2, routineName,
				  "Value of LAMBDA_REF should be larger than or equal 0.0");
		break;
	}

	if (strcmp(input.Stokes_input, "none")) {

		/* --- Magnetic field is specified --              -------------- */

		switch (topology) {
		case ONE_D_PLANE:
		case TWO_D_PLANE:
		case THREE_D_PLANE:
			if (atmos.B_char == 0.0) {
				Error(WARNING, routineName,
					  "Parameter atmos.B_char not set or set to zero\n"
					  " Wavelength grids of line profiles do not take account "
					  " of Zeeman splitting");
			}
			if (input.StokesMode == NO_STOKES) {
				sprintf(messageStr, "%s",
						"Keyword STOKES_MODE == NO_STOKES.\n"
						" Set to FIELD_FREE or FULL_STOKES\n"
						" when doing polarization calculations");
				Error(ERROR_LEVEL_1, routineName, messageStr);
			}
			break;
		case SPHERICAL_SYMMETRIC:
			Error(ERROR_LEVEL_2, routineName,
				  "Cannot accomodate magnetic fields in this topology");
		}
	} else {

		/* --- No magnetic field input specified --        -------------- */

		if (atmos.B_char != 0.0) {
			Error(WARNING, routineName,
				  "Ignoring value of keyword B_STRENGTH_CHAR when no "
				  "magnetic field is read");
		}
		if ((input.StokesMode == FIELD_FREE ||
			 input.StokesMode == FULL_STOKES) &&
			input.backgr_pol == FALSE) {
			Error(ERROR_LEVEL_2, routineName,
				  "Should not run STOKES_MODE == FIELD_FREE or FULL_STOKES\n "
				  "when no magnetic field is read and "
				  "BACKGROUND_POLARIZATION == FALSE");
		}
	}
	/* --- Hydrostatic equilibrium only in 1-D plane parallel -- ------ */

	switch (topology) {
	case ONE_D_PLANE:
		break;
	default:
		if (atmos.hydrostatic)
			Error(ERROR_LEVEL_2, routineName,
				  "Can only perform hydrostatic equilibrium calculation"
				  " in 1-D Cartesian geometry");
		break;
	}

	/* --- ANDERSON_* keywords are read by SOLVER = ANDERSON only,
		   NEWTON_* by SOLVER = NEWTON_KRYLOV only; set with another
		   solver they would be ignored silently -- ---------------- */

	/* --- The messages say when the solver is the default, so that a
		   keyword file written for the old default (NG) explains
		   itself (S-10) -- -------------------------------------- */

	input.solver_default = !keywordSet(theKeywords, Nkeyword, "SOLVER");
	sprintf(sname, "%s%s", solverName(input.solver),
			input.solver_default ? " (the default)" : "");
	for (p = 0; p < 2; p++) {
		prefix = (p == 0) ? "ANDERSON_" : "NEWTON_";
		owner  = (p == 0) ? SOLVER_ANDERSON : SOLVER_NEWTON_KRYLOV;
		if (input.solver == owner) continue;
		Nconflict = 0;
		for (n = 0; n < Nkeyword; n++) {
			if (!strncmp(theKeywords[n].keyword, prefix, strlen(prefix)) &&
				theKeywords[n].set) {
				sprintf(messageStr,
						"%s is set, but SOLVER = %s does not use it. Set"
						" SOLVER = %s\n or remove the keyword.",
						theKeywords[n].keyword, sname, solverName(owner));
				Error(WARNING, routineName, messageStr);
				Nconflict++;
			}
		}
		if (Nconflict) {
			sprintf(messageStr,
					"%s* keywords without SOLVER = %s; see the"
					" warnings above.",
					prefix, solverName(owner));
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
	}

	/* --- SOLVER = FIXED_POINT, ANDERSON and NEWTON_KRYLOV iterate the
		   map G of docs/SOLVER.md (S-2, S-3, S-4): one MALI step as a
		   function of an explicit state, with a fixed number of PRD
		   sub-iterations and no acceleration inside. Every keyword that
		   is set is checked against the solver, whatever the atoms
		   (SOLVER.md, Anderson as the default solver, item 2): a setting
		   the solver would not honour stops the run here, before any
		   work is done. A default is never a conflict. -- ------------- */

	if (input.solver != SOLVER_NG) {
		/* --- Keywords without an "off" value: an error when set.
			   FIXED_POINT keeps Ng for the "Ng on G" control. (An
			   "off" value states what the solver does and is
			   accepted: NG_ORDER, PRD_NG_ORDER, PRD_ITER_LIMIT and
			   COLL_RELAX_NSTEP 0, COLL_RELAX_FACTOR 1.) -- -- */
		static const struct {
			const char *name;
			bool_t      ng;
		} noOff[] = {{"NG_DELAY", TRUE},
					 {"NG_PERIOD", TRUE},
					 {"PRD_NG_DELAY", FALSE},
					 {"PRD_NG_PERIOD", FALSE}};

		Nconflict = 0;
		for (n = 0; n < (int) (sizeof(noOff) / sizeof(noOff[0])); n++) {
			if (noOff[n].ng && input.solver == SOLVER_FIXED_POINT) continue;
			if (keywordSet(theKeywords, Nkeyword, noOff[n].name)) {
				sprintf(messageStr,
						"%s is set, but SOLVER = %s does not use it. Set"
						" SOLVER = NG\n or remove the keyword.",
						noOff[n].name, sname);
				Error(WARNING, routineName, messageStr);
				Nconflict++;
			}
		}
		if ((input.solver == SOLVER_ANDERSON ||
			 input.solver == SOLVER_NEWTON_KRYLOV) &&
			input.Ngorder != 0) {
			sprintf(messageStr,
					"SOLVER = %s: NG_ORDER = %d would extrapolate the"
					" populations a second\n time after the solver's step."
					" Set SOLVER = NG, or remove the keyword or set it\n to 0.",
					sname, input.Ngorder);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		if (input.solver == SOLVER_ANDERSON &&
			(input.Anderson_depth < 0 || input.Anderson_depth > 50 ||
			 input.Anderson_delay < 0 || !(input.Anderson_cond_max > 1.0) ||
			 !(input.Anderson_safeguard >= 0.0))) {
			sprintf(messageStr,
					"SOLVER = ANDERSON: ANDERSON_DEPTH = %d (0 .. 50),"
					" ANDERSON_DELAY = %d (>= 0),\n ANDERSON_COND_MAX = %g"
					" (> 1), ANDERSON_SAFEGUARD = %g (>= 0, 0 = off):\n"
					" a value is out of range.",
					input.Anderson_depth, input.Anderson_delay,
					input.Anderson_cond_max, input.Anderson_safeguard);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		if (input.solver == SOLVER_NEWTON_KRYLOV &&
			(!(input.Newton_noise > 0.0 && input.Newton_noise < 1.0) ||
			 input.Newton_krylov_max < 1 || input.Newton_krylov_max > 200 ||
			 !(input.Newton_eta_max > 0.0 && input.Newton_eta_max < 1.0))) {
			sprintf(messageStr,
					"SOLVER = NEWTON_KRYLOV: NEWTON_NOISE = %g (0 .. 1),"
					" NEWTON_KRYLOV_MAX = %d (1 .. 200),\n NEWTON_ETA_MAX"
					" = %g (0 .. 1): a value is out of range.",
					input.Newton_noise, input.Newton_krylov_max,
					input.Newton_eta_max);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		/* --- The default 1.0E-2 is the NG solver's; unset, the map
			   does PRD_N_MAX_ITER sub-iterations -- -------------- */
		if (!keywordSet(theKeywords, Nkeyword, "PRD_ITER_LIMIT")) {
			input.PRDiterLimit = 0.0;
			setKeywordValue(theKeywords, Nkeyword, "PRD_ITER_LIMIT", "0.0");
		} else if (input.PRDiterLimit != 0.0) {
			sprintf(messageStr,
					"SOLVER = %s: PRD_ITER_LIMIT = %g lets the PRD"
					" sub-iterations stop early,\n which makes the map"
					" discontinuous. Set SOLVER = NG, or remove the keyword\n or"
					" set it to 0 (G then does PRD_N_MAX_ITER sub-iterations).",
					sname, input.PRDiterLimit);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		if (input.PRD_Ngorder != 0) {
			sprintf(messageStr,
					"SOLVER = %s: PRD_NG_ORDER = %d, but there is no"
					" acceleration inside\n the map. Set SOLVER = NG, or remove"
					" the keyword or set it to 0.",
					sname, input.PRD_Ngorder);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		if (keywordSet(theKeywords, Nkeyword, "COLL_RELAX_FACTOR") &&
			input.CR_factor != 1.0) {
			sprintf(messageStr,
					"SOLVER = %s: COLL_RELAX_FACTOR = %g, but the map has no"
					" collisional\n relaxation. Set SOLVER = NG, or remove the"
					" keyword or set it to 1.",
					sname, input.CR_factor);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		if (input.CR_Nstep != 0) {
			sprintf(messageStr,
					"SOLVER = %s: COLL_RELAX_NSTEP = %d makes the map"
					" depend on the\n iteration number. Set SOLVER = NG, or"
					" remove the keyword or set it to 0.",
					sname, input.CR_Nstep);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		if (atmos.hydrostatic) {
			sprintf(messageStr,
					"SOLVER = %s: HYDROSTATIC = TRUE changes the"
					" atmosphere between\n iterations. Set SOLVER = NG or"
					" HYDROSTATIC = FALSE.",
					sname);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		if (input.limit_memory) {
			sprintf(messageStr,
					"SOLVER = %s: the mean intensity J is part of the"
					" state and must be\n kept in memory. Set SOLVER = NG or"
					" LIMIT_MEMORY = FALSE.",
					sname);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		if (input.solve_ne == ITERATION && input.old_background) {
			sprintf(messageStr,
					"SOLVER = %s with SOLVE_NE = ITERATION recomputes the"
					" background\n from the state. Set SOLVER = NG or"
					" OLD_BACKGROUND = FALSE.",
					sname);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		if (input.solve_ne == ITERATION && input.NonICE) {
			sprintf(messageStr,
					"SOLVER = %s with SOLVE_NE = ITERATION computes the"
					" chemical\n equilibrium from the state. Set SOLVER = NG or"
					" NON_ICE = FALSE.",
					sname);
			Error(WARNING, routineName, messageStr);
			Nconflict++;
		}
		/* --- With ne in the state the background is stored in double
			   (docs/SOLVER.md §9 Q11, decided 2026-11-08): the float
			   file makes G a step function of ne, which limits
			   finite-difference Jacobian products to ~3e-3. -- -- */
		if (input.solve_ne == ITERATION) {
			if (keywordSet(theKeywords, Nkeyword, "BACKGROUND_DOUBLE") &&
				!input.background_double) {
				sprintf(messageStr,
						"SOLVER = %s with SOLVE_NE = ITERATION has ne in"
						" the state, and\n the background must then be stored"
						" in double. Set BACKGROUND_DOUBLE = TRUE\n or leave it"
						" out.",
						sname);
				Error(WARNING, routineName, messageStr);
				Nconflict++;
			}
			input.background_double = TRUE;
		}
		if (Nconflict) {
			sprintf(messageStr,
					"SOLVER = %s: %d keyword setting%s conflict with"
					" it; see the\n warning%s above. docs/SOLVER.md describes"
					" the solvers and their keywords.",
					sname, Nconflict, (Nconflict > 1) ? "s" : "",
					(Nconflict > 1) ? "s" : "");
			Error(ERROR_LEVEL_2, routineName, messageStr);
		}
	}

	/* --- If called with -showkeywords commandline option -- --------- */

	if (commandline.showkeywords) showValues(Nkeyword, theKeywords);

	/* --- Convert to MKSA units where necessary --       ------------- */

	atmos.vmicro_char *= KM_TO_M;
	atmos.vmacro_tresh *= KM_TO_M;
}
/* ------- end ---------------------------- readInput.c ------------- */

/* ------- begin -------------------------- keywordSet -------------- */

/* --- Whether the keyword was given in the keyword input file (its
	   default is never a conflict) -- --------------------------- */

static bool_t keywordSet(const Keyword *theKeywords, int Nkeyword,
						 const char *name)
{
	int n;

	for (n = 0; n < Nkeyword; n++) {
		if (!strcmp(theKeywords[n].keyword, name)) return theKeywords[n].set;
	}
	return FALSE;
}
/* ------- end ---------------------------- keywordSet -------------- */

/* ------- begin -------------------------- setKeywordValue --------- */

/* --- Change the value string of a keyword before or after reading:
	   a default that depends on the program, or the value in use for
	   -showkeywords -- ------------------------------------------- */

static void setKeywordValue(Keyword *theKeywords, int Nkeyword,
							const char *name, const char *value)
{
	int n;

	for (n = 0; n < Nkeyword; n++) {
		if (!strcmp(theKeywords[n].keyword, name)) {
			strcpy(theKeywords[n].value, value);
			return;
		}
	}
}
/* ------- end ---------------------------- setKeywordValue --------- */
