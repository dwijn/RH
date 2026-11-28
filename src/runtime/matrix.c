/* ------- file: -------------------------- matrix.c ----------------

	   Version:       rh2.0
	   Author:        Han Uitenbroek (huitenbroek@nso.edu)
	   Last modified: Thu Oct  8 16:57:27 2020 --

	   --------------------------                      ----------RH-- */

/* --- Routines to create 2-dimensional arrays that are contiguous in
	   memory. Therefore the whole array pointed to by pointer **matrix_??
	   can be addressed as either 1-dimensional array with pointer
	   matrix_??[0], or as 2-dimensional array as matrix_??[i][j].

 Note: space is initialized to zeros by using calloc rather than malloc.

	   --                                              -------------- */

#include <stdlib.h>

#include "rh.h"
#include "error.h"

/* --- Function prototypes --                          -------------- */

/* --- Global variables --                             -------------- */

/* ------- begin -------------------------- checkAlloc.c ------------ */

/* --- A failed allocation aborts with a message (Q-53) -- -------- */

static void checkAlloc(const char *routineName, void *theMatrix, void *Matrix, int Nrow,
					   int Ncol)
{
	if ((theMatrix == NULL && (size_t) Nrow * Ncol > 0) ||
		(Matrix == NULL && Nrow > 0)) {
		outOfMemory(routineName, "Could not allocate a %d x %d matrix", Nrow, Ncol);
	}
}
/* ------- end ---------------------------- checkAlloc.c ------------ */

/* ------- begin -------------------------- matrix_char.c ----------- */

char **matrix_char(RHState *rh, int Nrow, int Ncol)
{
	register int i;

	char *theMatrix, **Matrix;
	int   typeSize = sizeof(char), pointerSize = sizeof(char *);

	(void) rh; /* for RH_E_NOMEM in phase 4 */
	theMatrix = (char *) calloc((size_t) Nrow * Ncol, typeSize);
	Matrix    = (char **) malloc((size_t) Nrow * pointerSize);
	checkAlloc("matrix_char", theMatrix, Matrix, Nrow, Ncol);
	for (i = 0; i < Nrow; i++, theMatrix += Ncol)
		Matrix[i] = theMatrix;

	return Matrix;
}
/* ------- end ---------------------------- matrix_char.c ----------- */

/* ------- begin -------------------------- matrix_int.c ------------ */

int **matrix_int(RHState *rh, int Nrow, int Ncol)
{
	register int i;

	int *theMatrix, **Matrix, typeSize = sizeof(int),
							  pointerSize = sizeof(int *);

	(void) rh; /* for RH_E_NOMEM in phase 4 */
	theMatrix = (int *) calloc((size_t) Nrow * Ncol, typeSize);
	Matrix    = (int **) malloc((size_t) Nrow * pointerSize);
	checkAlloc("matrix_int", theMatrix, Matrix, Nrow, Ncol);
	for (i = 0; i < Nrow; i++, theMatrix += Ncol)
		Matrix[i] = theMatrix;

	return Matrix;
}
/* ------- end ---------------------------- matrix_int.c ------------ */

/* ------- begin -------------------------- matrix_double.c --------- */

double **matrix_double(RHState *rh, int Nrow, int Ncol)
{
	register int i;

	int     typeSize = sizeof(double), pointerSize = sizeof(double *);
	size_t  memSize;
	double *theMatrix, **Matrix;

	memSize = (size_t) Nrow * Ncol;
	(void) rh; /* for RH_E_NOMEM in phase 4 */
	theMatrix = (double *) calloc(memSize, typeSize);
	Matrix    = (double **) malloc((size_t) Nrow * pointerSize);
	checkAlloc("matrix_double", theMatrix, Matrix, Nrow, Ncol);
	for (i = 0; i < Nrow; i++, theMatrix += Ncol)
		Matrix[i] = theMatrix;

	return Matrix;
}
/* ------- end ---------------------------- matrix_double.c --------- */

/* ------- begin -------------------------- matrix_float.c ---------- */

float **matrix_float(RHState *rh, int Nrow, int Ncol)
{
	register int i;

	float *theMatrix, **Matrix;

	(void) rh; /* for RH_E_NOMEM in phase 4 */
	theMatrix = (float *) calloc((size_t) Nrow * Ncol, sizeof(float));
	Matrix    = (float **) malloc((size_t) Nrow * sizeof(float *));
	checkAlloc("matrix_float", theMatrix, Matrix, Nrow, Ncol);
	for (i = 0; i < Nrow; i++, theMatrix += Ncol)
		Matrix[i] = theMatrix;

	return Matrix;
}
/* ------- end ---------------------------- matrix_float.c ---------- */

/* ------- begin -------------------------- allocMatrix.c ----------- */

/* --- As matrix_double and matrix_int, but a failed allocation
	   returns RH_E_NOMEM (recorded with Fail) and leaves *matrix NULL:
	   for the large arrays (Q-54; all allocations in phase 4) -- --- */

rh_status allocMatrixDouble(RHState *rh, int Nrow, int Ncol, double ***matrix)
{
	register int i;

	double *theMatrix, **Matrix;

	*matrix   = NULL;
	theMatrix = (double *) calloc((size_t) Nrow * Ncol, sizeof(double));
	Matrix    = (double **) malloc((size_t) Nrow * sizeof(double *));
	if ((theMatrix == NULL && (size_t) Nrow * Ncol > 0) || (Matrix == NULL && Nrow > 0)) {
		free(theMatrix);
		free(Matrix);
		return Fail(rh, RH_E_NOMEM, "allocMatrixDouble", "Could not allocate a %d x %d matrix",
					Nrow, Ncol);
	}
	for (i = 0; i < Nrow; i++, theMatrix += Ncol)
		Matrix[i] = theMatrix;
	*matrix = Matrix;
	return RH_OK;
}

rh_status allocMatrixInt(RHState *rh, int Nrow, int Ncol, int ***matrix)
{
	register int i;

	int *theMatrix, **Matrix;

	*matrix   = NULL;
	theMatrix = (int *) calloc((size_t) Nrow * Ncol, sizeof(int));
	Matrix    = (int **) malloc((size_t) Nrow * sizeof(int *));
	if ((theMatrix == NULL && (size_t) Nrow * Ncol > 0) || (Matrix == NULL && Nrow > 0)) {
		free(theMatrix);
		free(Matrix);
		return Fail(rh, RH_E_NOMEM, "allocMatrixInt", "Could not allocate a %d x %d matrix",
					Nrow, Ncol);
	}
	for (i = 0; i < Nrow; i++, theMatrix += Ncol)
		Matrix[i] = theMatrix;
	*matrix = Matrix;
	return RH_OK;
}
/* ------- end ---------------------------- allocMatrix.c ----------- */

/* ------- begin -------------------------- freeMatrix.c ------------ */

void freeMatrix(RHState *rh, void **matrix)
{
	const char routineName[] = "freeMatrix";

	RH_ASSERT(rh, matrix != NULL && matrix[0] != NULL, routineName,
			  "Trying to free NULL pointer");
	free(matrix[0]);
	free(matrix);
}
/* ------- end ---------------------------- freeMatrix.c ------------ */
