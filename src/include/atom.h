/* ------- file: -------------------------- atom.h ------------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Sep 23 10:58:12 2025 --

	   --------------------------                      ----------RH-- */

#ifndef __ATOM_H__
#define __ATOM_H__

/* --- Defines atomic data structures. --              -------------- */

#define ATOM_LABEL_WIDTH  20
#define ATOM_ID_WIDTH     2
#define MOLECULE_ID_WIDTH 10

#define PRD_QCORE   2.0
#define PRD_QWING   4.0
#define PRD_QSPREAD 5.0
#define PRD_DQ      0.25

enum type { ATOMIC_LINE,
			ATOMIC_CONTINUUM,
			VIBRATION_ROTATION,
			MOLECULAR_ELECTRONIC };
enum ftype { FIXED_LINE,
			 FIXED_CONTINUUM };
enum Trad_option { TRAD_ATMOSPHERIC,
				   TRAD_PHOTOSPHERIC,
				   TRAD_CHROMOSPHERIC };
enum vdWaals { UNSOLD,
			   RIDDER_RENSBERGEN,
			   BARKLEM,
			   KURUCZ };
enum fit_type { KURUCZ_70,
				KURUCZ_85,
				SAUVAL_TATUM_84,
				IRWIN_81,
				TSUJI_73 };
enum Hund { CASE_A,
			CASE_B };
enum Barklemtype { SP,
				   PD,
				   DF };
enum orbit_am { S_ORBIT = 0,
				P_ORBIT,
				D_ORBIT,
				F_ORBIT };
enum zeeman_cpl { LS_COUPLING = 0,
				  JK_COUPLING,
				  JJ_COUPLING };

/* --- Structure prototypes --                         -------------- */

typedef struct Atom             Atom;
typedef struct Molecule         Molecule;
typedef struct AtomicTransition AtomicTransition;
typedef struct MolTransition    MolTransition;
typedef struct AtomicLine       AtomicLine;
typedef struct ZeemanMultiplet  ZeemanMultiplet;
typedef struct rhthread         rhthread;
typedef struct Paschenstruct    Paschenstruct;

/* --- Structure defines radiative transition --       -------------- */

struct AtomicLine {
	bool_t       symmetric, polarizable, Voigt, PRD;
	enum vdWaals vdWaals;
	int          i, j, Nlambda, Nblue, Ncomponent, Nxrd, fd_profile,
		**id0, **id1;
	double lambda0, *lambda, isotope_frac, g_Lande_eff,
		Aji, Bji, Bij, *Rij, *Rji, *wphi, *Qelast, Grad, cvdWaals[4],
		cStark, qcore, qwing, **rho_prd, *c_shift, *c_fraction, **frac;

	/* --- Wavelength integration weights (nm) on the merged grid where
		   they differ from the trapezoid; NULL otherwise (N-62,
		   mirrorQuadrature) --                         -------------- */

	double *wquad;

	/* --- Absorption profiles (and the Stokes and magneto-optical
		   profiles), stored in single precision: they are by far the
		   largest arrays in angle-dependent runs (M-1). They are
		   computed in double, rounded once on storage, and wphi
		   normalizes the stored values. --             -------------- */

	float **phi, **phi_Q, **phi_U, **phi_V, **psi_Q, **psi_U, **psi_V;
	FILE   *fp_GII;
	long   *prd_block; /* sizes of the per-location blocks of the
						  folded PRD weights in fp_GII (scatter.c) */
	struct Ng       *Ng_prd;
	Atom            *atom;
	AtomicLine     **xrd;
	ZeemanMultiplet *zm;
	/* --- Per-thread Rij/Rji accumulators (one row per thread, see
		   addtoRates() in fillgamma.c), the same pattern as
		   rhth[nt].Gamma for addtoGamma(), so no lock is needed.
		   Only allocated when N_THREADS > 1; with one thread the shared
		   Rij/Rji are accumulated into directly. -- */
	double **Rij_thread, **Rji_thread;
};

typedef struct {
	bool_t     hydrogenic;
	int        i, j, Nlambda, Nblue;
	double     lambda0, *lambda, isotope_frac, alpha0, *alpha, *Rij, *Rji;
	double    *wquad; /* as in AtomicLine (N-62) */
	Atom      *atom;
	double   **Rij_thread, **Rji_thread;
	PchipTable bftable; /* Cross-sections for passive Metal_bf */
} AtomicContinuum;

typedef struct {
	enum type type;
	enum Hund Hundi, Hundj;
	bool_t    symmetric, Voigt, polarizable;
	char      configi[3], configj[3], parityi[2], parityj[2];
	int       vi, vj, Nlambda, Nblue, subi, subj, Lambdai, Lambdaj, ecnoi, ecnoj;
	double    lambda0, *lambda, isotope_frac, Ei, Ej, gi, gj, Si, Sj,
		Omegai, Omegaj, *wphi, g_Lande_eff,
		Grad, qcore, qwing, Aji, Bji, Bij;
	double          *wquad; /* as in AtomicLine (N-62) */
	float          **phi;   /* Absorption profile, single precision (M-1) */
	Molecule        *molecule;
	ZeemanMultiplet *zm;
} MolecularLine;

typedef struct {
	enum ftype       type;
	enum Trad_option option;
	int              i, j;
	double           lambda0, strength, Trad;
	Atom            *atom;
} FixedTransition;

struct AtomicTransition {
	enum type type;
	union {
		AtomicLine      *line;
		AtomicContinuum *continuum;
	} ptype;
	Atom *atom;
};

struct MolTransition {
	enum type type;
	union {
		MolecularLine *vrline;
	} ptype;
	Molecule *molecule;
};

struct rhthread {
	double **gij, **Vij, **wla, **Uji_down, *eta,
		**Gamma;
	int Nrt_alloc;
};

struct Atom {
	char          ID[ATOM_ID_WIDTH + 1], **label, *popsinFile, *popsoutFile;
	bool_t        active, NLTEpops, popsInit;
	enum solution initial_solution;
	int           Nlevel, Nline, Ncont, Nfixed, Nprd, *stage, periodic_table,
		activeindex;
	long   offset_coll;
	double abundance, weight, *g, *E, **C, *vbroad, **n, **nstar,
		*ntotal, **Gamma;
	AtomicLine         *line;
	AtomicContinuum    *continuum;
	FixedTransition    *ft;
	struct Ng          *Ng_n;
	rhthread           *rhth;
	FILE               *fp_input;
	struct RateReports *rateReports;
};

typedef struct {
	char   ID[ATOM_ID_WIDTH + 1];
	bool_t abundance_set;
	int   *mol_index, Nstage, Nmolecule;
	double weight, abund, *ionpot, **pf, **n;
	Atom  *model;
} Element;

struct Molecule {
	char          ID[MOLECULE_ID_WIDTH + 1], *popsFile, *configs;
	bool_t        active;
	enum fit_type fit;
	enum solution initial_solution;
	int          *pt_index, *pt_count, Nelement, Nnuclei, Npf, Neqc, Nrt,
		charge, Nconfig, Nv, NJ, activeindex;
	double Ediss, Tmin, Tmax, weight, *vbroad, *pf_coef, *eqc_coef,
		*pf, **pfv, *n, **nv, **nvstar, *C_ul, **Gamma;
	MolecularLine *mrt;
	struct Ng     *Ng_nv;
	rhthread      *rhth;
};

typedef struct {
	int             L, L1, l1, l2, l;
	double          g, E, S, J, S1, J1, j1, j2, K, gL, hfs;
	enum zeeman_cpl cpl;
	bool_t          zm_explicit;
} RLK_level;

typedef struct {
	bool_t       polarizable;
	enum vdWaals vdwaals;
	int          pt_index, stage, isotope;
	double       lambda0, Bji, Aji, Bij,
		Grad, GStark, GvdWaals, hyperfine_frac,
		isotope_frac, iso_dl, cross, alpha;
	RLK_level        level_i, level_j;
	ZeemanMultiplet *zm;
} RLK_Line;

struct ZeemanMultiplet {
	int     Ncomponent, *q;
	double *shift, *strength, g_eff;
};

typedef struct {
	int     N1, N2;
	double *neff1, *neff2, **cross, **alpha;
} Barklemstruct;

struct Paschenstruct {
	int      Nj;
	double  *eigenval;
	double **C;
};

/* --- Associated function prototypes --               -------------- */

void initSolution(RHState *rh);
void Iterate(RHState *rh, int NmaxIter, double iterLimit);

void   readAtomicModels(RHState *rh);
void   readMolecularModels(RHState *rh);
void   OHCH_bf_check(RHState *rh);
void   statEquil(RHState *rh, Atom *atom, int isum);
double updatePopulations(RHState *rh, int niter);

void   CollisionRate(RHState *rh, Atom *atom, FILE *atomFile);
void   Damping(RHState *rh, AtomicLine *line, double *adamp);
void   FixedRate(RHState *rh, Atom *atom);
void   freeAtom(RHState *rh, Atom *atom);
void   freeRateReports(Atom *atom);
void   freeAtomicLine(RHState *rh, AtomicLine *line);
void   freeAtomicContinuum(RHState *rh, AtomicContinuum *continuum);
void   getfjk(RHState *rh, Element *element, double ne, int k, double *fjk, double *dfjk);
void   getLambda(RHState *rh, AtomicLine *line);
int    asymmetricBlends(RHState *rh, const double *extra, int Nextra);
int    mirrorQuadrature(RHState *rh);
void   getLambdaCont(AtomicContinuum *continuum, double lambdamin);
double getwlambda_line(AtomicLine *line, int la);
double getwlambda_cont(AtomicContinuum *continuum, int la);
void   initAtom(Atom *atom);
void   initAtomicLine(AtomicLine *line);
void   initAtomicContinuum(AtomicContinuum *continuum);

void   initZeeman(ZeemanMultiplet *zm);
void   freeZeeman(ZeemanMultiplet *zm);
double zm_gamma(double J, double S, double L);

void initGammaAtom(RHState *rh, Atom *atom, int iter);
void initGammaMolecule(RHState *rh, Molecule *molecule);
void zeroGammaThreads(RHState *rh, Atom *atom);
void zeroGammaThreadsMolecule(RHState *rh, Molecule *molecule);
void mergeGammaThreads(RHState *rh, Atom *atom);
void mergeGammaThreadsMolecule(RHState *rh, Molecule *molecule);

void LTEpops(RHState *rh, Atom *atom, bool_t Debeye);
void LTEpops_elem(RHState *rh, Element *element);
void SetLTEQuantities(RHState *rh);
void getProfiles(RHState *rh);
void Profile(RHState *rh, AtomicLine *line);
void readProfile(RHState *rh, AtomicLine *line, int lamu, float *phi);
void writeProfile(RHState *rh, AtomicLine *line, int lamu, float *phi);
void readAtom(RHState *rh, Atom *atom, char *atomFileName);
void readPopulations(RHState *rh, Atom *atom);
void SortLambda(RHState *rh);
void Stark(RHState *rh, AtomicLine *line, double *GStark);
void StarkLinear(RHState *rh, AtomicLine *line, double *GStark);
void VanderWaals(RHState *rh, AtomicLine *line, double *GvdW);
void writeAtom(RHState *rh, Atom *atom);
void writePopulations(RHState *rh, Atom *atom);
void zeroRates(RHState *rh, bool_t redistribute);
void mergeRates(RHState *rh, bool_t redistribute);

bool_t readRadRate(RHState *rh, Atom *atom);
bool_t writeRadRate(RHState *rh, Atom *atom);
bool_t readCollisionRate(RHState *rh, Atom *atom);
bool_t writeCollisionRate(RHState *rh, Atom *atom);
bool_t readDamping(RHState *rh, Atom *atom);
bool_t writeDamping(RHState *rh, Atom *atom);

bool_t readBarklemTable(RHState *rh, enum Barklemtype type, Barklemstruct *bs);
bool_t getBarklemcross(RHState *rh, Barklemstruct *bs, RLK_Line *rlk);
bool_t getBarklemactivecross(RHState *rh, AtomicLine *line);
bool_t getBarklemExplicit(RHState *rh, AtomicLine *line);
bool_t determinate_abo(RHState *rh, char *label, int *l);

/* --- Associated function prototypes --               -------------- */

void statEquilMolecule(RHState *rh, Molecule *molecule, int isum);

void   COcollisions(RHState *rh, Molecule *molecule);
void   MolecularDamping(RHState *rh, MolecularLine *mrt, double *adamp);
double equilconstant(RHState *rh, Molecule *molecule, double T);
void   freeMolecule(RHState *rh, Molecule *molecule);
void   freeMolLine(RHState *rh, MolecularLine *mrt);
double getwlambda_mrt(MolecularLine *mrt, int la);
void   initMolecule(Molecule *molecule);
void   initMolLine(MolecularLine *mrt, enum type line_type);
void   MolecularProfile(RHState *rh, MolecularLine *mrt);
void   mrt_locate(int N, MolecularLine *lines, double lambda, int *low);
void   LTEmolecule(RHState *rh, Molecule *molecule);
double partfunction(RHState *rh, Molecule *molecule, double T);
void   readMolecule(RHState *rh, Molecule *molecule, char *fileName, bool_t active);
void   readMolecularLines(RHState *rh, Molecule *molecule, char *line_data);
void   readMolPops(RHState *rh, Molecule *molecule);
void   writeMolPops(RHState *rh, Molecule *molecule);

/* --- Redistribution function --                      -------------- */

/* --- rh->redistribute (redistribute.c): whether the PRD_DAMPING note
	   is in the log, once per run --                  -------------- */

typedef struct {
	bool_t dampingNoted;
} RedistributeState;

double GII(double adamp, double waveRatio, double q_emit, double q_abs);
double RII(RHState *rh, double v_emit, double v_abs, double adamp, int mu1, int mu2);
void   Redistribute(RHState *rh, int NmaxIter, double iterLimit);
double PRDDamp(RHState *rh, AtomicLine *line, const double *old, long N);
long   PRDRhoSize(RHState *rh, const AtomicLine *line);
bool_t PRDRhoCheck(RHState *rh, const Atom *atom, int kr);
void   PRDScatter(RHState *rh, AtomicLine *line);
void   PRDAngleScatter(RHState *rh, AtomicLine *PRDline);

/* --- Polarization related --                         -------------- */

void   adjustStokesMode(RHState *rh);
bool_t determinate(RHState *rh, char *label, double g, int *n, double *S, int *L,
				   double *J);
double effectiveLande(RHState *rh, AtomicLine *line);
bool_t zeemanPolarizable(RHState *rh, Atom *atom, int i, int j, double g_Lande_eff);
double Lande(double S, int L, double J);
void   StokesProfile(AtomicLine *line);
void   Zeeman(RHState *rh, AtomicLine *line);
void   MolZeeman(RHState *rh, MolecularLine *mrt);
double MolLande_eff(RHState *rh, MolecularLine *mrt);
int    getOrbital(RHState *rh, char orbit);
double ZeemanStrength(RHState *rh, double Ju, double Mu, double Jl, double Ml);

double w3js(RHState *rh, double j1, double j2, double j3,
			double m1, double m2, double m3);
double w6js(RHState *rh, double j1, double j2, double j3,
			double J1, double J2, double J3);

#endif /* !__ATOM_H__ */

/* ------- end ---------------------------- atom.h ------------------ */
