/* ------- file: -------------------------- background.h ------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Wed Jan 16 14:10:02 2013 --

	   --------------------------                      ----------RH-- */

#ifndef __BACKGROUND_H__
#define __BACKGROUND_H__

/* --- Include file for background opacities and emissivities. -- --- */

/* --- Chemical equilibrium (N-61): the largest number of Newton steps,
	   and the relative residual of the conservation equations below
	   which one more full step ends the iteration (the result is then
	   converged to rounding). A location that does not converge stops
	   the run or invalidates a trial point --          -------------- */

#define N_MAX_CHEM_ITER 100
#define CHEM_ITER_LIMIT 1.0E-12

/* --- Fixed number of Newton steps for the electron density in the
	   map G with ne in the state (S-2b), from the ne of the state -- */

#define N_NE_ITER_OF_STATE 10

/* --- Power of two by which emissivities are multiplied in the
	   (single-precision) background file, see backgroundIO in readj.c */

#define BACKGR_ETA_EXP 100

/* --- Maximum number of iterations and minimum accuracy for hydrostatic
	   equilibrium calculations --                     -------------- */

#define N_MAX_HSE_ITER 20
#define HSE_ITER_LIMIT 1.0E-2

/* --- Structure definitions --                        -------------- */

struct Linelist {
	bool_t      used;
	double     *adamp;
	AtomicLine *line;
};

/* --- Associated function prototypes --               -------------- */

void ChemicalEquilibrium(int NmaxIter, double iterLimit);
void ChemicalEquilibriumOfState(int NmaxIter, double iterLimit);

/* --- The chemical equilibrium of one location in log densities of the
	   free nuclei (N-61, chemequil.c): Nn nuclei, Ns species with
	   stoichiometry nu[s][0 .. Nn-1] and ln c_s in lnc[s] (-HUGE_VAL:
	   absent); n, g, du, ut, sc and H are work space (Ns, Nn, Nn, Nn,
	   Nn, Nn x Nn). chemSolve returns the number of Newton steps, with
	   n the densities of the solution u, or -1 if not converged in
	   maxIter --                                       -------------- */

typedef struct {
	int      Nn, Ns;
	double **nu, *lnc, *n, *g, *du, *ut, *sc, **H;
} ChemSystem;

int    chemSolve(ChemSystem *cs, const double *a, double *u, int maxIter, double tol,
				 double *res, double *step);
void   BackgroundOfState(bool_t write_analyze_output);
double chemLastChangeReset(void);

/* --- The state of the background opacities between calls, their
	   parts of RHState (Q-46) --                      -------------- */

/* --- rh->background (background.c): whether the Kurucz lines have
	   been read (they are kept when Background() runs more than once),
	   and the ne solves so far (with SOLVE_NE = ITERATION the first
	   starts from scratch) --                         -------------- */

typedef struct {
	bool_t rlk_read;
	int    ne_iter;
} BackgroundState;

/* --- rh->kurucz (kurucz.c): per-point tables that depend only on
	   the atmosphere, the element or the ionization stage, and not on
	   wavelength. They are filled by rlk_prepare, are read-only
	   afterwards, and are released by rlk_cleanup -- -------------- */

typedef struct {
	double ***lnU;    /* [element][stage][k], only for (element, stage) used */
	double  **vbroad; /* [element][k], only for elements used                */
	double   *T03;    /* [k], pow(T, 0.3), only if Unsold lines are present  */
} RLK_Tables;

typedef struct {
	RLK_Tables rlk_tab;
} KuruczState;

/* --- rh->hydrogen (hydrogen.c): per depth point, the fractional
	   index into the theta table of H2- free-free and into the
	   temperature table of H2+ free-free, made by hydrogenPrepare from
	   atmos.T --                                      -------------- */

typedef struct {
	double *H2mff_theta_index;
	double *H2pff_temp_index;
} HydrogenState;

/* --- rh->chemequil (chemequil.c): the largest change of a log
	   density in the last Newton step of the chemical equilibrium of
	   state, since chemLastChangeReset --             -------------- */

typedef struct {
	double chemLastChange;
} ChemEquilState;

double Gaunt_bf(double lambda, double n_eff, int charge);
double Gaunt_ff(double lambda, int charge, double T);
int    Gaunt_ff_readTable(const char *fileName, char *errmsg, size_t len);
void   Gaunt_ff_freeTable(void);
void   FMetals(double *Fmetal);
void   dFMetals(double *dFdne);
void   Hydrostatic(int NmaxIter, double iterLimit);

void   Background(bool_t analyze_output, bool_t equilibria_only);
void   backgrOpac(int Nlambda, double *lambda);
void   backgroundPrepare(void);
void   backgroundCleanup(void);
bool_t duplicateLevel(Atom *atom, char *label);
bool_t duplicateLine(Atom *atom, char *labeli, char *labelj);
void   readBackground(int nspect, int mu, bool_t top_to_bottom);
void   mapBackground(void);
void   unmapBackground(void);
int    backgroundRecords(int nspect);
int    writeBackground(int nspect, int mu, bool_t to_obs,
					   double *chi_c, double *eta_c, double *sca_c,
					   double *chip_c);
void   writeBRS(void);
void   readBRS(void);

void readKuruczLines(char *fileName);
int  rlk_ascend(const void *v1, const void *v2);
void rlk_locate(int N, RLK_Line *lines, double lambda, int *low);
void rlk_prepare(void);
void rlk_cleanup(void);

bool_t Hminus_bf(double lambda, double *chi, double *eta);
double Hminus_bf_sigma(double lambda);
int    Hminus_bf_readTable(const char *fileName, char *errmsg, size_t len);
void   Hminus_bf_freeTable(void);
void   Hminus_bf_joins(double *short_ratio, double *long_ratio);
bool_t Hminus_ff(double lambda, double *chi);
double Hminus_ff_kappa(double lambda, double theta);
bool_t Hydrogen_bf(double lambda, double *chi, double *eta);
void   Hydrogen_ff(double lambda, double *chi);
bool_t H2minus_ff(double lambda, double *chi);
bool_t H2plus_ff(double lambda, double *chi);
void   hydrogenPrepare(void);
void   hydrogenCleanup(void);
bool_t Rayleigh(double lambda, Atom *atom, double *chi);
bool_t Rayleigh_H2(double lambda, double *chi);
void   Thomson(double *chi);
bool_t Metal_bf(double lambda, int Nmetal, Atom *metals,
				double *chi, double *eta);
void   metalPrepare(void);
bool_t OH_bf_opac(RHState *rh, double lambda, double *chi, double *eta);
bool_t CH_bf_opac(RHState *rh, double lambda, double *chi, double *eta);

#endif /* !__BACKGROUND_H__ */

/* ------- end ---------------------------- background.h ------------ */
