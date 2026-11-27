/* ------- file: -------------------------- spectrum.h --------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Tue Feb 20 09:47:52 2024 --

	   --------------------------                      ----------RH-- */

#ifndef __SPECTRUM_H__
#define __SPECTRUM_H__

/* --- Defines structure containing the overall set of wavelengths,
	   and the transitions that are active at each of the wavelengths.
	   --                                              -------------- */

#define N_MAX_OVERLAP 500

#define VACUUM_TO_AIR_LIMIT 200.0000
#define AIR_TO_VACUUM_LIMIT 199.9352

/* --- File for storage of intensities in case of angle-dependent PRD */

#define IMU_FILENAME "Imu.dat"

/* --- Stores set of active transitions at each wavelength and pointers
	   to temporary storage for opacity and emissivity components --  */

typedef struct {
	int *Nactiveatomrt, *Nactivemolrt,
		*Nlower, *Nupper, **lower_levels, **upper_levels;
	double            *chi, *eta, *chip, *chi_c, *eta_c, *sca_c, *chip_c;
	AtomicTransition **art;
	MolTransition    **mrt;
} ActiveSet;

/* --- Stores emergent intensities and array of active sets -- ------ */

typedef struct {
	bool_t  vacuum_to_air, updateJ;
	int     Nspect, *PRDindex, fd_J, fd_J20, fd_Imu;
	double *lambda, **J, **I, **Stokes_Q, **Stokes_U, **Stokes_V, **J20,
		**Jgas, **Jgas_fixed, **Ilast, **v_los;
	ActiveSet *as;
	int        NJgasFixed; /* number of rows of Jgas_fixed (N-27) */
	int       *PRDrow;     /* N-27: -1 wavelength without PRD line, else the row of Jgas_fixed
							  of a PRD wavelength a wavelength without PRD line can reach,
							  -2 for one it cannot */
	long            *nc;   /* cumulative coefficient counts (B-32) */
	int             *iprdh;
	double          *cprdh;
	pthread_mutex_t *Jgas_lock; /* one per wavelength row of Jgas (T-4) */
} Spectrum;

/* --- Associated function prototypes --               -------------- */

double Formal(RHState *rh, int nspect, bool_t eval_operator, bool_t redistribute);
double solveSpectrum(RHState *rh, bool_t eval_operator, bool_t redistribute);

void addtoGamma(RHState *rh, int nspect, double wmu, double *P, double *Psi);
void addtoRates(RHState *rh, int nspect, int mu, bool_t to_obs, double wmu,
				double *I, bool_t redistribute);
void initScatter(RHState *rh);

void StokesK(RHState *rh, int nspect, int k, double chi_I, double K[4][4]);
void addtoCoupling(RHState *rh, int nspect);
void addPreconditionedRates(int Nlevel, int i, int j, double twohnu3_c2,
							const double *Vij, const double *wla,
							const double *gij, const double *n_i,
							const double *n_j, const double *U_i,
							const double *U_j, const double *I,
							const double *Psi, double wmu, int Nspace,
							double **Gam);

/* --- What type of lines are present in active set as? -- ---------- */

bool_t containsPolarized(RHState *rh, ActiveSet *as);
bool_t containsBoundBound(RHState *rh, ActiveSet *as);
bool_t containsPRDline(RHState *rh, ActiveSet *as);
bool_t containsActive(RHState *rh, ActiveSet *as);

void init_as(ActiveSet *as);
void alloc_as(RHState *rh, int nspect, bool_t crosscoupling);
void free_as(RHState *rh, int nspect, bool_t crosscoupling);

/* --- The opacity buffers of alloc_as, pooled per thread slot (see
	   opacity.c), its part of RHState (rh->opacity, Q-46) -- ------ */

typedef struct {
	double *chi_c, *eta_c, *sca_c, *chip_c, *chi, *eta, *chip;
	int     Nrec_c, Nrec_as;
} AsPool;

typedef struct {
	AsPool         *as_pool; /* input.Nthreads slots, made at first use */
	pthread_mutex_t as_pool_lock;
} OpacityState;

#define OPACITY_STATE_INIT {.as_pool_lock = PTHREAD_MUTEX_INITIALIZER}
void initSpectrum(void);
void Opacity(RHState *rh, int nspect, int mu, bool_t top_to_bottom, bool_t activate);
void Planck(int Nspace, double *T, double lambda0, double *Bnu);

void readJlambda(RHState *rh, int nspect, double *J);
void writeJlambda(RHState *rh, int nspect, double *J);
void readJ20lambda(RHState *rh, int nspect, double *J20);
void writeJ20lambda(RHState *rh, int nspect, double *J20);
void writeJgas(RHState *rh, double **Jgas);
void initJgasRows(RHState *rh, const char *isPRD);
void beginJgas(RHState *rh, bool_t redistribute);
void endJgas(RHState *rh, bool_t redistribute);
void addJgasRow(RHState *rh, int nspect, double wmu, double *I);
void addJgasRay(RHState *rh, int nspect, int mu, int to_obs, double wmu, double *I);
long prdhIndex(RHState *rh, int nspect, int mu, int to_obs, int k);
void initJgasLocks(RHState *rh);
void freeJgasBuffer(void);
void readJgas(RHState *rh, double **Jgas);

void readImu(RHState *rh, int nspect, int mu, bool_t to_obs, double *I);
void writeImu(RHState *rh, int nspect, int mu, bool_t to_obs, double *I);
void writeSpectrum(RHState *rh, Spectrum *spectrum);
void freeSpectrum(RHState *rh, Spectrum *spectrum);
void writeOpacity(RHState *rh);

/* --- Wavelength conversion --                         ------------- */

void vacuum_to_air(int Nlambda, double *lambda_vac, double *lambda_air);
void air_to_vacuum(int Nlambda, double *lambda_air, double *lambda_vac);

#endif /* !__SPECTRUM_H__ */

/* ------- end ---------------------------- spectrum.h -------------- */
