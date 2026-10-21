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

extern char messageStr[];

/* ------- begin -------------------------- checkAlloc.c ------------ */

static void checkAlloc(const char *routineName, void *theMatrix, void *Matrix,
					   int Nrow, int Ncol)
{
	if ((theMatrix == NULL && (size_t) Nrow * Ncol > 0) ||
		(Matrix == NULL && Nrow > 0)) {
		sprintf(messageStr, "Unable to allocate a %d x %d matrix", Nrow, Ncol);
		Error(ERROR_LEVEL_2, routineName, messageStr);
	}
}
/* ------- end ---------------------------- checkAlloc.c ------------ */

/* ------- begin -------------------------- matrix_char.c ----------- */

char **matrix_char(int Nrow, int Ncol)
{
	register int i;

	char *theMatrix, **Matrix;
	int   typeSize = sizeof(char), pointerSize = sizeof(char *);

	theMatrix = (char *) calloc((size_t) Nrow * Ncol, typeSize);
	Matrix    = (char **) malloc((size_t) Nrow * pointerSize);
	checkAlloc("matrix_char", theMatrix, Matrix, Nrow, Ncol);
	for (i = 0; i < Nrow; i++, theMatrix += Ncol)
		Matrix[i] = theMatrix;

	return Matrix;
}
/* ------- end ---------------------------- matrix_char.c ----------- */

/* ------- begin -------------------------- matrix_int.c ------------ */

int **matrix_int(int Nrow, int Ncol)
{
	register int i;

	int *theMatrix, **Matrix, typeSize = sizeof(int),
							  pointerSize = sizeof(int *);

	theMatrix = (int *) calloc((size_t) Nrow * Ncol, typeSize);
	Matrix    = (int **) malloc((size_t) Nrow * pointerSize);
	checkAlloc("matrix_int", theMatrix, Matrix, Nrow, Ncol);
	for (i = 0; i < Nrow; i++, theMatrix += Ncol)
		Matrix[i] = theMatrix;

	return Matrix;
}
/* ------- end ---------------------------- matrix_int.c ------------ */

/* ------- begin -------------------------- matrix_double.c --------- */

double **matrix_double(int Nrow, int Ncol)
{
	register int i;

	int     typeSize = sizeof(double), pointerSize = sizeof(double *);
	size_t  memSize;
	double *theMatrix, **Matrix;

	memSize   = (size_t) Nrow * Ncol;
	theMatrix = (double *) calloc(memSize, typeSize);
	Matrix    = (double **) malloc((size_t) Nrow * pointerSize);
	checkAlloc("matrix_double", theMatrix, Matrix, Nrow, Ncol);
	for (i = 0; i < Nrow; i++, theMatrix += Ncol)
		Matrix[i] = theMatrix;

	return Matrix;
}
/* ------- end ---------------------------- matrix_double.c --------- */

/* ------- begin -------------------------- freeMatrix.c ------------ */

void freeMatrix(void **matrix)
{
	const char routineName[] = "freeMatrix";

	if (matrix == NULL || matrix[0] == NULL) {
		sprintf(messageStr, "Trying to free NULL pointer");
		Error(ERROR_LEVEL_2, routineName, messageStr);
	} else {
		/* --- Free the memory allocated for matrix matrix -- ----------- */

		free(matrix[0]);
		free(matrix);
	}
}
/* ------- end ---------------------------- freeMatrix.c ------------ */
