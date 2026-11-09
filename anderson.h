/* ------- file: -------------------------- anderson.h --------------

	   Anderson acceleration (type II) on a history of iterate pairs,
	   see anderson.c and S-3 in docs/SOLVER.md.

	   --------------------------                      ----------RH-- */

#ifndef __ANDERSON_H__
#define __ANDERSON_H__

typedef struct {
	long            N;        /* length of the vectors */
	int             m;        /* depth: at most m difference columns */
	int             Npairs;   /* pairs held, at most m + 1 */
	int             newest;   /* ring index of the newest pair */
	double        **u, **g;   /* [m + 1][N]: iterates and their images */
	unsigned char **valid;    /* [m + 1][N]: entry usable in that pair */
	unsigned char  *mask;     /* [N]: usable in every pair held */
	double        **A;        /* [m][N]: weighted columns, then Q */
	double         *b;        /* [N]: weighted residual */
	double         *R, *Rinv; /* [m * m] */
	double         *r;        /* [m]: scratch */

	/* --- Safeguard state and counters (andersonStep) -- -------- */

	int    trial;     /* the last point returned was extrapolated */
	double fnorm_acc; /* ||W f|| of the newest accepted pair */
	long   Nextrapolated, Nrejected, Nplain, Ndropped;
} Anderson;

/* --- Options of andersonStep: the blocks of the vector (start[0] = 0
	   ... start[Nblock] = N), one coefficient set for all of them or
	   one per block, the condition limit and the safeguard factor
	   c (0: off) --                                  -------------- */

typedef struct {
	int         Nblock;
	const long *start;
	int         per_block;
	double      cond_max, safeguard;
} AndersonOptions;

enum AndersonAction { ANDERSON_PLAIN,
					  ANDERSON_EXTRAPOLATED,
					  ANDERSON_REJECTED };

/* --- Result of one least-squares solve --         -------------- */

typedef struct {
	int    Ncols; /* columns used */
	int    Ndrop; /* oldest columns left out for the condition number */
	double cond;  /* condition number of R (1-norm) with Ncols columns */
} AndersonSolve;

Anderson            *andersonInit(long N, int m);
void                 andersonFree(Anderson *a);
void                 andersonReset(Anderson *a, int keep_newest);
void                 andersonPush(Anderson *a, const double *u, const double *g,
								  const unsigned char *valid);
int                  andersonNcols(const Anderson *a);
AndersonSolve        andersonSolve(Anderson *a, long lo, long hi, const double *w,
								   double cond_max, double *gamma);
void                 andersonApply(const Anderson *a, long lo, long hi, const double *gamma,
								   int Ncols, double theta, double *u_out);
const unsigned char *andersonMask(const Anderson *a);
double               andersonNorm(long N, const double *u, const double *g,
								  const unsigned char *valid, const double *w);
enum AndersonAction  andersonStep(Anderson *a, const AndersonOptions *opt,
								  const double *u, const double *g,
								  const unsigned char *valid, const double *w,
								  int extrapolate, double *u_next,
								  double *gamma, AndersonSolve *info);
void                 andersonDamp(const Anderson *a, const AndersonOptions *opt,
								  const double *gamma, const AndersonSolve *info,
								  double theta, double *u_next);
void                 andersonAbandon(Anderson *a, double *u_next);

#endif /* !__ANDERSON_H__ */

/* ------- end ---------------------------- anderson.h -------------- */
