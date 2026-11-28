/* ------- file: -------------------------- inputs.h ----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Feb 20 09:34:35 2024 --

	   --------------------------                      ----------RH-- */

#ifndef __INPUTS_H__
#define __INPUTS_H__

/* --- Include file for command line options and keyword input
	   parameters. --                                  -------------- */

#define MAX_KEYWORD_LENGTH 32
#define MAX_VALUE_LENGTH   160

#define CHECK_OPTION(option, name, length) \
	(strncmp(option, name, MAX((size_t) (length), strlen(option))) == 0)

enum keywordtype { KEYWORD_REQUIRED,
				   KEYWORD_DEFAULT,
				   KEYWORD_OPTIONAL };

enum S_interpol { S_LINEAR,
				  S_BEZIER3 };
enum S_interpol_stokes { DELO_BEZIER3 };

enum order_2D { LINEAR_2D,
				BEZIER_2D };
enum order_3D { LINEAR_3D,
				BICUBIC_3D };
enum ne_solution { NONE,
				   ONCE,
				   ITERATION };
enum solver_type { SOLVER_NG,
				   SOLVER_FIXED_POINT,
				   SOLVER_ANDERSON,
				   SOLVER_NEWTON_KRYLOV };
enum solver_variables { VARIABLES_LOG,
						VARIABLES_RELATIVE };
enum bg_operator { BG_OPERATOR_DIAGONAL,
				   BG_OPERATOR_TRIDIAGONAL }; /* N-51 */
enum anderson_coefficients { ANDERSON_GLOBAL,
							 ANDERSON_PER_BLOCK };
enum anderson_extrapolate { EXTRAPOLATE_ALL,
							EXTRAPOLATE_POPULATIONS,
							EXTRAPOLATE_POPULATIONS_RHO };

typedef struct {
	char   name[MAX_KEYWORD_LENGTH];
	int    minlength;
	bool_t value_required;
	char   value[MAX_VALUE_LENGTH];
	void  *pointer;
	rh_status (*setValue)(RHState *rh, char *value, void *pointer);
	char message[MAX_VALUE_LENGTH];
} Option;

typedef struct {
	char   keyword_input[MAX_VALUE_LENGTH];
	bool_t quiet, showkeywords;
	FILE  *logfile;
} CommandLine;

typedef struct {
	char             keyword[MAX_KEYWORD_LENGTH], value[MAX_VALUE_LENGTH];
	bool_t           set;
	enum keywordtype type;
	void            *pointer;
	rh_status (*setValue)(RHState *rh, char *value, void *pointer);
} Keyword;

typedef struct {
	char atmos_input[MAX_VALUE_LENGTH],
		abund_input[MAX_VALUE_LENGTH],
		wavetable_input[MAX_VALUE_LENGTH],
		atoms_input[MAX_VALUE_LENGTH],
		molecules_input[MAX_VALUE_LENGTH],
		Stokes_input[MAX_VALUE_LENGTH],
		KuruczData[MAX_VALUE_LENGTH],
		pfData[MAX_VALUE_LENGTH],
		fudgeData[MAX_VALUE_LENGTH],
		atmos_output[MAX_VALUE_LENGTH],
		spectrum_output[MAX_VALUE_LENGTH],
		geometry_output[MAX_VALUE_LENGTH],
		opac_output[MAX_VALUE_LENGTH],
		JFile[MAX_VALUE_LENGTH],
		background_File[MAX_VALUE_LENGTH],
		background_ray_File[MAX_VALUE_LENGTH],
		H_atom[MAX_VALUE_LENGTH],
		H2_molecule[MAX_VALUE_LENGTH],
		radrateFile[MAX_VALUE_LENGTH],
		collrateFile[MAX_VALUE_LENGTH],
		dampingFile[MAX_VALUE_LENGTH],
		coolingFile[MAX_VALUE_LENGTH],
		Itop[MAX_VALUE_LENGTH],
		HminusBfData[MAX_VALUE_LENGTH],
		GauntFfData[MAX_VALUE_LENGTH];
	bool_t magneto_optical, XRD,
		backgr_pol, limit_memory, allow_passive_bb, NonICE,
		rlkscatter, xdr_endian, old_background, accelerate_mols,
		background_double,
		RLK_explicit, prdh_limit_mem;
	enum solution startJ;

	enum StokesMode        StokesMode;
	enum S_interpol        S_interpolation;
	enum S_interpol        S_interpolation_longchar; /* N-39 */
	enum S_interpol_stokes S_interpolation_stokes;
	enum PRDangle          PRD_angle_dep;
	enum order_2D          interpolate_2D;
	enum order_3D          interpolate_3D;
	enum ne_solution       solve_ne;
	enum solver_type       solver;
	bool_t                 solver_default;               /* SOLVER not set in the file */
	bool_t                 Anderson_extrapolate_default; /* ANDERSON_EXTRAPOLATE not set */

	enum solver_variables      Anderson_variables;
	enum anderson_coefficients Anderson_coefficients;
	enum anderson_extrapolate  Anderson_extrapolate;
	int                        Anderson_depth, Anderson_delay;
	int                        Anderson_safeguard_window;
	double                     Anderson_cond_max, Anderson_safeguard;

	enum solver_variables Newton_variables;
	int                   Newton_krylov_max;
	double                Newton_noise, Newton_eta_max;

	int isum, Ngdelay, Ngorder, Ngperiod, NmaxIter,
		PRD_NmaxIter, PRD_Ngdelay, PRD_Ngorder, PRD_Ngperiod,
		NmaxScatter, Nthreads, CR_Nstep;
	double iterLimit, PRDiterLimit, PRD_damping, metallicity, CR_factor;

	/* --- Background scattering per wavelength (N-51, bgscatter.c) -- */
	int              bg_NmaxIter, bg_depth;
	double           bg_iterLimit;
	enum bg_operator bg_operator;
	bool_t           bg_NmaxIter_set;               /* BACKGROUND_N_MAX_ITER in the file */
	double           bg_linearAlbedo, bg_linearTau; /* N-52 */

	/* --- ERROR_LIMIT (N-3): the estimated error at which ANDERSON and
		   NEWTON_KRYLOV stop, 0 = not set. Without ITER_LIMIT,
		   iterLimit takes its value -- ------------------------------ */

	double errorLimit;
	bool_t iterLimit_set; /* ITER_LIMIT set in the file */

	pthread_attr_t thread_attr;
} InputData;

/* --- Associated function prototypes --               -------------- */

int                    getLineOrEOF(FILE *inputFile, const char *commentChar, char *line);
RH_NODISCARD rh_status getLine(RHState *rh, FILE *inputFile, const char *commentChar,
							   char *line);
RH_NODISCARD rh_status parse(RHState *rh, int argc, char *argv[], int Noption,
							 Option *theOptions);
RH_NODISCARD rh_status readInput(RHState *rh);
RH_NODISCARD rh_status readValues(RHState *rh, FILE *fp_keyword, int Nkeyword,
								  Keyword *theKeywords);

RH_NODISCARD rh_status setAngleSet(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setcharValue(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setboolValue(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setintValue(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setdoubleValue(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setstartValue(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setnesolution(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setsolver(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setStokesMode(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setPRDangle(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setThreadValue(RHState *rh, char *value, void *pointer);

RH_NODISCARD rh_status setSolverVariables(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setAndersonCoefficients(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setBgOperator(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setAndersonExtrapolate(RHState *rh, char *value, void *pointer);
const char            *solverName(enum solver_type solver);

RH_NODISCARD rh_status setInterpolate_2D(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status setInterpolate_3D(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status set_S_Interpolation(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status set_S_interpolation_longchar(RHState *rh, char *value, void *pointer);
RH_NODISCARD rh_status set_S_interpolation_stokes(RHState *rh, char *value, void *pointer);

void showValues(int Nkeyword, Keyword *theKeywords);
void UpperCase(char *string);

void writeInput(RHState *rh);

#endif /* !__INPUTS_H__ */

/* ------- end ---------------------------- inputs.h ---------------- */
