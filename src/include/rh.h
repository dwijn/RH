/* ------- file: -------------------------- rh.h --------------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Feb 22 17:06:40 2024 --

	   --------------------------                      ----------RH-- */

/* --- General include file for Rybicki & Hummer family of radiative
	   transfer codes. --                              -------------- */

/* --- Include stdio.h here because it defines NULL
	 - Include rpc/types.h here because it defines bool_t
	 - Include pthread.h for POSIX threads
	   --                                              -------------- */

#ifndef __RH_H__
#define __RH_H__

#include <stdio.h>
#include <rpc/types.h>
#include <pthread.h>

#include "alloc.h"

/* --- The state of a run (state.h); a parameter of most functions
	   (Q-49) --                                       -------------- */

typedef struct RHState RHState;

/* --- The status of a function that can fail (Q-51; the codes and how
	   to use them: error.h) --                        -------------- */

typedef enum {
	RH_OK = 0,
	RH_E_INPUT,
	RH_E_DATA,
	RH_E_NUMERIC,
	RH_E_TRIAL,
	RH_E_NOMEM,
	RH_E_STOPPED,
	RH_E_STATE,
	RH_E_INTERNAL
} rh_status;

/* --- A function that returns an rh_status is declared RH_NODISCARD:
	   the compiler warns when a caller drops the status (and -Werror
	   makes that an error) --                         -------------- */

#define RH_NODISCARD __attribute__((warn_unused_result))

/* --- SPHERICAL_SYMMETRIC is no longer a geometry (rhsphere was
	   retired in Q-43); it keeps its value because geometry.out
	   stores the topology as a number (python/rhanalyze, idl) -- */

enum Topology { ONE_D_PLANE,
				TWO_D_PLANE,
				SPHERICAL_SYMMETRIC,
				THREE_D_PLANE };
enum FeautrierOrder { STANDARD,
					  FEAUTRIER_HERMITE };
enum Interpolation { LINEAR,
					 SPLINE,
					 EXP_SPLINE };
enum solution { UNKNOWN = -1,
				LTE_POPULATIONS,
				ZERO_RADIATION,
				OLD_POPULATIONS,
				NEW_J,
				OLD_J };
enum StokesMode { NO_STOKES,
				  FIELD_FREE,
				  FULL_STOKES };
enum PRDangle { PRD_ANGLE_AVER = 0,
				PRD_ANGLE_APPROX,
				PRD_ANGLE_DEP };

#define MAX_LINE_SIZE      512
#define MAX_MESSAGE_LENGTH 512
#define MAX_KEYWORD_SIZE   32

#define LG10 2.30258509299404568402

#ifndef MAX
#define MAX(x, y) (((x) > (y) ? (x) : (y)))
#endif
#ifndef MIN
#define MIN(x, y) (((x) < (y) ? (x) : (y)))
#endif

#define SWAPPOINTER(a, b) \
	{                     \
		void *tmp;        \
		tmp = a;          \
		a   = b;          \
		b   = tmp;        \
	}
#define SWAPDOUBLE(a, b) \
	{                    \
		double tmp;      \
		tmp = a;         \
		a   = b;         \
		b   = tmp;       \
	}

#define SQ(x)    ((x) * (x))
#define CUBE(x)  ((x) * (x) * (x))
#define POW10(x) exp(LG10 *(x))

#define MODULO(n, N) (((n) >= 0) ? ((n) % (N)) : ((((n) % (N)) + (N)) % (N)))

#define FLUX_DOT_OUT "flux.out"
#define J20_DOT_OUT  "J20.out"

#define PERMISSIONS 0666

/* --- Function prototypes --                          -------------- */

void                   SetFPEtraps(RHState *rh);
RH_NODISCARD rh_status checkNread(RHState *rh, int Nread, int Nrequired,
								  const char *routineName, int checkPoint);
RH_NODISCARD rh_status setOptions(RHState *rh, int argc, char *argv[]);
double                 vproject(RHState *rh, int k, int mu);
void                   Bproject(RHState *rh);
bool_t                 StopRequested(RHState *rh);
char                 **getWords(char *label, char *separator, int *count);

/* --- Matrix manipulation --                          -------------- */

char                 **matrix_char(RHState *rh, int Nrow, int Ncol);
int                  **matrix_int(RHState *rh, int Nrow, int Ncol);
double               **matrix_double(RHState *rh, int Nrow, int Ncol);
float                **matrix_float(RHState *rh, int Nrow, int Ncol);
RH_NODISCARD rh_status allocMatrixDouble(RHState *rh, int Nrow, int Ncol, double ***matrix);
RH_NODISCARD rh_status allocMatrixInt(RHState *rh, int Nrow, int Ncol, int ***matrix);
void                   freeMatrix(RHState *rh, void **Matrix);
RH_NODISCARD rh_status SolveLinearEq(RHState *rh, int N, double **A, double *b,
									 bool_t improve);
RH_NODISCARD rh_status LUdecomp(RHState *rh, int N, double **A, int *index, double *d);
void                   LUbacksubst(int N, double **A, int *index, double *b);

/* --- Rate equations by GTH elimination (gth.c, N-31) -- ---------- */

enum GTHStatus { GTH_OK,
				 GTH_NEGATIVE,
				 GTH_NONFINITE,
				 GTH_REDUCIBLE };

enum GTHStatus SolveGTH(int N, double **R, int root, double ntotal, double *n,
						int *i_bad, int *j_bad);

/* --- Interpolation routines --                       -------------- */

int qsascend(const void *v1, const void *v2);
int qsdescend(const void *v1, const void *v2);

void   Hunt(int n, const double *array, double value, int *iLower);
void   Locate(int n, const double *array, double value, int *theIndex);
double bilinear(int Ncol, int Nrow, const double *f, double x, double y);

void Linear(int Ntable, double *xtable, double *ytable,
			int N, double *x, double *y, bool_t hunt);

enum tau_unity_status { TAU_UNITY_OK,
						TAU_NOT_INCREASING,
						TAU_NOT_BRACKETED };

int                    heightAtUnityTau(int N, const double *tau_ref, const double *z, double *z1);
RH_NODISCARD rh_status zeroPointAtUnityTau(RHState *rh, const char *routineName,
										   const char *atmosID, int N, const double *tau_ref,
										   const double *z, double *z0);

/* --- Monotone piecewise cubic Hermite interpolation (pchip.c, N-32),
	   optionally in log x and log y. The table is copied; evaluation
	   keeps no state. -- --------------------------------------------- */

enum PchipScale { PCHIP_LINEAR = 0,
				  PCHIP_LOGX   = 1,
				  PCHIP_LOGY   = 2,
				  PCHIP_LOGLOG = 3 };

typedef struct {
	int     N;
	bool_t  logx, logy;
	double *x, *y, *d;
	double *ylin, *dlin; /* logy with a node <= 0 (N-49), else NULL */
} PchipTable;

RH_NODISCARD rh_status pchipTableInit(RHState *rh, PchipTable *pt, int N, const double *x,
									  const double *y, int scale);
void                   pchipTableEval(const PchipTable *pt, int N, const double *x, double *y);
void                   pchipTableFree(PchipTable *pt);

void   cc_kernel(double s, double *u);
double cubeconvol(int Nx, int Ny, double *f, double x, double y);

/* --- Special functions --                            -------------- */

RH_NODISCARD rh_status E1(RHState *rh, double x, double *e1);
RH_NODISCARD rh_status E2(RHState *rh, double x, double *e2);

void   GaussLeg(double x1, double x2, double *x, double *w, int n);
double Voigt(double a, double v, double *F);
double gammln(double xx);

void w2(double dtau, double *w);
void w3(double dtau, double *w);

void U3(double dtau, double *U);

#endif /* !__RH_H__ */

/* ------- end ---------------------------- rh.h -------------------- */
